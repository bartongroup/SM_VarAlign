# seems to be required at the top, otherwise get display error...
import argparse
import logging
import os

import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import tqdm
import vcf
from Bio import AlignIO
from matplotlib.backends.backend_pdf import PdfPages
from numpy import vectorize

from varalign import aacon
from varalign import alignments
from varalign import analysis_toolkit
from varalign import ensembl
from varalign import gnomad
from varalign import jalview
from varalign import occ_gmm
from varalign.config import defaults
from varalign.utils import make_dir_if_needed


logging.basicConfig(filename='align_variants.log', format='%(asctime)s %(name)s [%(levelname)-8s] - %(message)s')
log = logging.getLogger(__name__)
log.setLevel('INFO')


def _chunk_alignment(aln, n):
    """
    Return a generator that provides chunks of an alignment (sequence-wise).

    :param aln: Multiple sequence alignment.
    :param n: Chunk size.
    :return: Generator of n-sized alignment chunks.
    """
    return (aln[i:i + n] for i in range(0, len(aln), n))


def _chunk_table(table, n):
    """
    Return a generator that provides chunks of a DataFrame (row-wise).

    :param table: DataFrame.
    :param n: Chunk size.
    :return: Generator of n-sized alignment chunks.
    """
    return (table.iloc[i:i + n] for i in range(0, len(table), n))


def _dump_table_and_log(method, path, what):
    method(path)
    log.info('{} saved to {}'.format(what, path))


def _build_vep_filter(canonical=eval(defaults.canonical), consequences=defaults.consequences,
                      additional=defaults.additional):
    """
    Build a VEP filter string suitable for DataFrame.query().

    :param canonical:
    :param consequences:
    :param additional:
    :return:
    """
    query = []
    if consequences != ['']:
        query.append('Consequence in {}'.format(consequences))
    if canonical:
        query.append('CANONICAL == "YES"')
    if additional != '':
        query.append(additional)
    return ' & '.join(query)


def _default_variant_filter(variants_table):
    """
    Apply standard filters to a alignment derived variant table.

    :param variants_table:
    :return:
    """
    # See version 1 of notebook for other ideas (e.g. Protin_position in UniProt range...)
    # Reduce transcript duplication
    is_canonical = variants_table[('VEP', 'CANONICAL')] == 'YES'
    is_ccds = variants_table[('VEP', 'CCDS')] != ''
    # Only want those that can map to a residue
    is_protein_coding = variants_table[('VEP', 'BIOTYPE')] == 'protein_coding'
    at_protein_position = variants_table[('VEP', 'Protein_position')] != ''
    # Filter least useful effects
    is_not_modifier = variants_table[('VEP', 'IMPACT')] != 'MODIFIER'
    # Source protein filter
    swissprot_matches_source = (variants_table['External', 'SOURCE_ACCESSION'] == variants_table['VEP', 'SWISSPROT'])
    vcontains = vectorize(lambda x, y: x in y)
    trembl_matches_source = vcontains(variants_table[('External', 'SOURCE_ACCESSION')],
                                      variants_table[('VEP', 'TREMBL')])
    trembl_matches_source[:] = False  # OVERRIDE TREMBL TO KEEP ONLY SWISSPROT
    # Apply filter
    filtered_variants = variants_table.loc[is_canonical & is_protein_coding & is_not_modifier & is_ccds &
                                           (swissprot_matches_source | trembl_matches_source) &
                                           at_protein_position].copy()
    return filtered_variants


def _map_uniprot_to_genome(uniprot, species='homo_sapiens', collapse=True):
    """
    Map a UniProt entry to the genome.

    :param uniprot:
    :param species:
    :param collapse:
    :return:
    """
    # Map to Gene IDs
    ensembl_genes = ensembl.get_xrefs(uniprot, species=species,
                                      features='gene')  # TODO: This could be transcrpts or translations...
    if len(ensembl_genes) == 0:
        return None  # no mapping
    ensembl_ranges = [ensembl.get_genomic_range(x) for x in ensembl_genes]
    for i in range(len(ensembl_genes)):
        log.info('Mapped {} to {} on chr: {}, {}-{}'.format(uniprot, ensembl_genes[i], *ensembl_ranges[i]))
    # Identify and remove non-standard sequence regions
    non_standard_ranges = [i for i, x in enumerate(ensembl_ranges) if x[0] not in ensembl.standard_regions]
    if len(non_standard_ranges) > 0:
        non_standard_ranges.sort(reverse=True)
        [ensembl_genes.pop(i) for i in non_standard_ranges]
        [ensembl_ranges.pop(i) for i in non_standard_ranges]
        log.info('Removed %s non-standard sequence regions from %s.', len(non_standard_ranges), uniprot)
    # Check for no mapping
    if len(ensembl_ranges) == 0:
        raise ValueError('Could not map {} to the genome.'.format(uniprot))  # TODO: handle this...
    # Collapse ranges if desired
    if collapse:
        ensembl_ranges = ensembl.merge_ranges(ensembl_ranges, min_gap=1000)
    return ensembl_ranges


def _mapping_table(alignment_info):
    """
    Construct a alignment column to sequence residue mapping table.

    :param alignment_info:
    :return:
    """
    mapping_table = pd.DataFrame(alignment_info['mapping'].tolist(),
                                 index=[alignment_info.index, alignment_info['seq_id']])  # From top of notebook
    mapping_table.reset_index(inplace=True)
    mapping_table = pd.melt(mapping_table, id_vars=['level_0', 'seq_id'])
    mapping_table.dropna(subset=['value'], inplace=True)
    indexed_map_table = pd.DataFrame(mapping_table['value'].tolist(),
                                     columns=['Column', 'Protein_position'],
                                     index=[mapping_table['seq_id']]).reset_index()
    indexed_map_table = indexed_map_table.set_index(['seq_id', 'Protein_position']).sort_index()
    indexed_map_table.index.rename(['SOURCE_ID', 'Protein_position'], inplace=True)
    indexed_map_table.columns = pd.MultiIndex.from_tuples([('Alignment', x) for x in indexed_map_table.columns],
                                                          names=['Type', 'Field'])
    return indexed_map_table


def _occupancy_from_mapping_table(indexed_mapping_table):
    """
    Calculate column occupancy from a mapping table.

    :param indexed_mapping_table:
    :return:
    """
    column_occupancy = indexed_mapping_table[('Alignment', 'Column')].value_counts().sort_index()
    column_occupancy.name = 'occupancy'
    return column_occupancy


def _interpret_regression_results(regression_table, p_threshold=0.05, action=None):
    """
    Provide human interpretation of regression results.
    
    :param regression_table: Output from `analysis_toolkit._comparative_regression` (DataFrame)
    :param p_threshold: Significance threshold for regression p-value (float)
    :param action: Method to call on  interpretations, e.g. `print` or `logging.info`
    :return: Nothing if action is unspecified otherwise a list of strings.
    """
    # Synonymous variant counts should NOT be positively correlated with conservation
    negative_control_p = regression_table.loc['filtered_synonymous', 'pvalue'] > p_threshold
    negative_control_m = regression_table.loc['filtered_synonymous', 'slope'] < 0  # TODO: Why sometimes signif -ve?
    pass_negative = negative_control_p or negative_control_m

    # Missense variant counts should be positively correlated with conservation
    positive_control_p = regression_table.loc['filtered_missense', 'pvalue'] < p_threshold
    positive_control_m = regression_table.loc['filtered_missense', 'slope'] > 0
    pass_positive = positive_control_p and positive_control_m

    # Human readable interpretation
    results = ['Filtered synonymous vs. Shenkin (negative control)... {}'.format('PASS' if pass_negative else 'FAIL'),
               'Filtered missense vs Shenkin (positive control)... {}'.format('PASS' if pass_positive else 'FAIL')]

    # Either return strings...
    if not action:
        return results

    # ...or call action
    for r in results:
        action(r)


def _write_variants_as_features(alignment_variant_table, feature_file_name):
    """
    Write a Jalview feature file marking up the variants in the alignment.

    :param alignment_variant_table:
    :param feature_file_name:
    :return:
    """
    jalview.create_jalview_feature_file({'missense_variant': 'red', 'synonymous_variant': 'blue'}, feature_file_name)
    for (seq_id, consequence), variant_table in alignment_variant_table['VEP'].groupby(['SOURCE_ID', 'Consequence']):
        if consequence in ('missense_variant', 'synonymous_variant'):
            residue_indexes = list(variant_table.index.get_level_values(1))
            #     residue_indexes = [x - 1 + int(seq_id.split('/')[1].split('-')[0]) for x in residue_indexes]
            variant_ids = list(variant_table['Existing_variation'])
            jalview.append_jalview_variant_features(seq_id.split('/')[0], residue_indexes, variant_ids, consequence,
                                                    feature_file_name)
    log.info('Wrote alignment variants to Jalview feature file %s', feature_file_name)


def get_genome_mappings(aln_info_table, species):
    """

    :param aln_info_table:
    :param species:
    :return:
    """
    # TODO: If get transcript ID can use to filter variant table (duplicate)
    genomic_ranges = [
        (row.seq_id, _map_uniprot_to_genome(row.uniprot_id, species=species))
        for row in tqdm.tqdm(aln_info_table.itertuples(), total=len(aln_info_table), desc='Mapping sequences...')
    ]
    if len(genomic_ranges) == 0:
        log.error('Failed to map any sequences to the genome... Are you sure there are human sequences?')
        raise ValueError
    log.info("Mapped {} sequences to genome.".format(len(genomic_ranges)))
    # Format to table
    genomic_mapping_table = pd.DataFrame(genomic_ranges, columns=['seq_id', 'genomic_ranges'])
    return genomic_mapping_table


def get_gnomad_variants(aln_info_table):
    """

    :param aln_info_table:
    :return:
    """
    # NB. gnomad fetcher is packed into a generator, which is extracted in the following list comp.
    sequence_variant_lists = [(row.seq_id, (x for r in row.genomic_ranges for x in gnomad.gnomad.fetch(*r)))
                              for row in aln_info_table.dropna(subset=['genomic_ranges']).itertuples()]
    all_variants = [(variant, seq_id)
                    for seq_id, range_reader in tqdm.tqdm(sequence_variant_lists, desc='Loading variants...')
                    for variant in range_reader]

    # pass VCF records and source_ids
    n = 1000  # chunking seems to interact with redundant rows... Fix by adding chunk ID with `keys`
    variants_table = pd.concat([gnomad.vcf_row_to_table(*list(zip(*all_variants[i:i + n])))
                                for i in tqdm.tqdm(range(0, len(all_variants), n), desc='Parsing variants...')],
                               keys=list(range(0, len(all_variants), n)))

    # Write alignment variants to a VCF
    with open('alignment_variants.vcf', 'w') as vcf_out:  # TODO: add alignment to file name? (needs refactoring...)
        vcf_writer = vcf.Writer(vcf_out, gnomad.gnomad)
        for v, _ in all_variants:
            vcf_writer.write_record(v)

    return variants_table


def map_variants_to_alignment(variants_df, residue_column_map):
    """
    Add alignment column numbers to a variant table.

    :param variants_df: Unaligned variant table (DataFrame)
    :param residue_column_map:
    :return:
    """
    # Coerce Protein_position to correct type
    variants_df.loc[:, ('VEP', 'Protein_position')] = pd.to_numeric(variants_df.loc[:, ('VEP', 'Protein_position')],
                                                                    errors='coerce')
    # Set index for merge
    variants_df.reset_index(['SITE', 'ALLELE_NUM', 'Feature'], inplace=True)
    variants_df.set_index(('VEP', 'Protein_position'), append=True, inplace=True)
    variants_df.index.set_names(['SOURCE_ID', 'Protein_position'], inplace=True)
    variants_df.sort_index(inplace=True)
    # Merge to map
    aligned_variants = residue_column_map.join(variants_df)  # Drops variants that map outside alignment
    aligned_variants.sort_index(inplace=True)
    return aligned_variants


def align_variants(aln_info_table, species='HUMAN'):
    """

    :param species:
    :return:
    """

    # ----- Map sequences to genome -----
    # TODO: If get transcript ID can use to filter variant table
    genomic_mapping_table = get_genome_mappings(aln_info_table, species)
    aln_info_table = aln_info_table.merge(genomic_mapping_table, on=['seq_id'], how='left')

    # ----- Fetch variants for the mapped genomic ranges -----
    variants_table = get_gnomad_variants(aln_info_table)

    # ----- Add source UniProt identifiers to the table -----
    # Create UniProt ID series that shares an index with the variant table
    source_uniprot_ids = aln_info_table.set_index('seq_id')['uniprot_id']
    source_uniprot_ids.name = ('External', 'SOURCE_ACCESSION')
    source_uniprot_ids.index.name = 'SOURCE_ID'
    # Add IDs to variant tables
    variants_table = variants_table.join(source_uniprot_ids)

    # ----- Filter variant table -----
    filtered_variants = _default_variant_filter(variants_table)
    log.info('Redundant rows:\t{}'.format(sum(filtered_variants.reset_index('Feature').index.duplicated())))
    filtered_variants.reset_index(level=0, drop=True, inplace=True)  # Remove chunk ID
    log.info('Total rows:\t{}'.format(len(filtered_variants)))

    # ----- Map variants to columns -----
    # Generate alignment column / sequence residue mapping table
    indexed_map_table = _mapping_table(aln_info_table)
    aligned_variants = map_variants_to_alignment(filtered_variants, indexed_map_table)

    return aligned_variants


def cli(argv=None, logger=log):
    # CLI
    parser = argparse.ArgumentParser(description='Align variants to a Pfam alignment.')
    parser.add_argument('alignment', type=str, help='Path to the alignment.')
    parser.add_argument('--max_gaussians', type=int, default=5,
                        help='Maximum number of Gaussians for occupancy fitting.')
    parser.add_argument('--n_groups', type=int, default=1, help='Top Gaussians to select after occupancy fitting.')
    parser.add_argument('--override', help='Override any previously generated files.', action='store_true')
    parser.add_argument('--species', type=str, help='Species (used for alignment filtering)', default='HUMAN')
    args = parser.parse_args(argv)

    if logger:
        # Log arguments
        for arg, value in sorted(vars(args).items()):
            logger.info("Command line argument %s: %r", arg, value)
        # TODO: or just `log.info(args)`

    return args


def main(args):
    # Results and data will be written in these folders
    results_path = 'results'
    make_dir_if_needed(results_path)
    results_prefix = os.path.join(results_path, args.alignment)
    data_path = os.path.join('.varalign', 'aligned_variants_data')
    make_dir_if_needed(data_path)
    data_prefix = os.path.join(data_path, args.alignment)
    alignment = AlignIO.read(args.alignment, format='stockholm')

    # Check if data is available from previous run
    is_data_available = all([os.path.isfile(data_prefix + '_variants.p.gz'),
                             os.path.isfile(data_prefix + '_info.p.gz'),
                             os.path.isfile(data_prefix + '_mappings.p.gz')])


    # Run align variants pipeline in chunks
    # Parse alignment info
    log.info('Generating alignment info table...')
    alignment_info = alignments.alignment_info_table(alignment, args.species)
    log.info('Alignment info table head:\n%s', alignment_info.head().to_string())
    if args.override or not is_data_available:
        # TODO: Chunk size should be optimised? Also, its effectiveness depends on human sequences in each chunk...
        chunk_size = 500
        vartable_chunks = []
        chunked_info = _chunk_table(alignment_info, chunk_size)
        n_chunks = len(list(range(0, len(alignment_info), chunk_size)))
        for chunk in tqdm.tqdm(chunked_info, desc='Alignment chunks...', total=n_chunks):
            _alignment_variant_table = align_variants(chunk)
            vartable_chunks.append(_alignment_variant_table)
        alignment_variant_table = pd.concat(vartable_chunks)

        indexed_mapping_table = _mapping_table(alignment_info)  # TODO: Should be passed or returned by align_variants?
        # Write data
        _dump_table_and_log(alignment_info.to_pickle, data_prefix + '_info.p.gz', 'Alignment info table pickle')
        _dump_table_and_log(alignment_variant_table.to_pickle, data_prefix + '_variants.p.gz',
                            'Alignment variant table pickle')
        _dump_table_and_log(indexed_mapping_table.to_pickle, data_prefix + '_mappings.p.gz',
                            'Alignment mapping table pickle')
    else:
        log.info('Loading data for {}...'.format(args.alignment))
        alignment_info = pd.read_pickle(data_prefix + '_info.p.gz')
        alignment_variant_table = pd.read_pickle(data_prefix + '_variants.p.gz')
        indexed_mapping_table = pd.read_pickle(data_prefix + '_mappings.p.gz')

    # AACon
    # Run AACon and save results
    alignment_conservation = aacon.get_aacon(alignment)
    _dump_table_and_log(alignment_conservation.to_csv, results_prefix + '_aacon_scores.csv',
                        'Formatted AACons results')

    # The remainder is pretty much all analysis, plotting and formatting (e.g., to Jalview output)

    # Calculate column variant aggregations and save results
    # Count variants over columns
    column_variant_counts = analysis_toolkit.count_column_variant_consequences(alignment_variant_table)
    _dump_table_and_log(column_variant_counts.to_csv, results_prefix + '.col_var_counts.csv',
                        'Column variant counts')
    # Count *rare* variants over columns
    rare_maf_threshold = 0.001
    is_rare = alignment_variant_table[('Allele_INFO', 'AF_POPMAX')] < rare_maf_threshold
    column_rare_counts = analysis_toolkit.count_column_variant_consequences(alignment_variant_table[is_rare])
    _dump_table_and_log(column_rare_counts.to_csv, results_prefix + '.col_rare_counts.csv',
                        'Column rare variant counts')
    # Count ClinVar annotations for *missense* variants over columns
    is_missense = alignment_variant_table[('VEP', 'Consequence')] == 'missense_variant'
    column_missense_clinvar = analysis_toolkit.count_column_clinvar(alignment_variant_table[is_missense])
    _dump_table_and_log(column_missense_clinvar.to_csv, results_prefix + '.col_mis_clinvar.csv',
                        'Column missense variant ClinVar annotation frequencies')
    # Count ClinVar annotations for *synonymous* variants over columns
    is_synonymous = alignment_variant_table[('VEP', 'Consequence')] == 'synonymous_variant'
    column_synonymous_clinvar = analysis_toolkit.count_column_clinvar(alignment_variant_table[is_synonymous])
    _dump_table_and_log(column_synonymous_clinvar.to_csv, results_prefix + '.col_syn_clinvar.csv',
                        'Column synonymous variant ClinVar annotation frequencies')
    # Use mapping table to calculate human residue occupancy
    # TODO: Adjust for unmapped seqs
    column_occupancy = _occupancy_from_mapping_table(indexed_mapping_table)
    # Merge required data for further standard analyses; this is saved after missense scores are added
    column_summary = column_variant_counts.join([column_missense_clinvar, column_occupancy, alignment_conservation])

    # Occupancy GMM
    gmms = occ_gmm._fit_mixture_models(column_summary['occupancy'], args.max_gaussians)
    M_best = occ_gmm._pick_best(gmms['models'], gmms['data'])
    # M_best.means_
    subset_mask_gmm = occ_gmm._core_column_mask(M_best, gmms['data'], args.n_groups)
    column_summary = column_summary.assign(column_gmm_pass=subset_mask_gmm)

    # Regression statistics
    # This checks whether missense and synonymous variant counts are correlated with column occupancy before and
    # after column filtering
    variants_vs_occ = analysis_toolkit._comparative_regression(column_summary, 'occupancy', filter_mask=subset_mask_gmm)
    _dump_table_and_log(variants_vs_occ.to_csv, results_prefix + '.variant_occ_regression.csv',
                        'Variant vs. occupancy regression parameters')
    # TODO: Test variants_vs_occ.loc['filtered_missense', 'pvalue'] > 0.05
    # Conservation plane with Shenkin score
    shenkin_regressions = analysis_toolkit._comparative_regression(column_summary, 'shenkin',
                                                                   filter_mask=subset_mask_gmm)
    _dump_table_and_log(shenkin_regressions.to_csv, results_prefix + '.variant_shenkin_regression.csv',
                        'Variant vs. Shenkin regression parameters')
    _interpret_regression_results(shenkin_regressions, action=log.info)

    # Column variant scores, for block columns only
    missense_scores = analysis_toolkit._column_variant_scores(column_summary[subset_mask_gmm],
                                                              variant_class='missense_variant',
                                                              occupancy='occupancy')
    _dump_table_and_log(missense_scores.to_csv, results_prefix + '.col_missense_scores.csv', 'Column missense scores')
    column_summary = column_summary.join(missense_scores)
    # Add shenkin percentile rank
    column_summary = column_summary.join(column_summary.loc[subset_mask_gmm, 'shenkin'].rank(pct=True),
                                         rsuffix='_percentile')
    _dump_table_and_log(column_summary.to_csv, results_prefix + '.col_summary.csv', 'Column summary data')

    # Plot output
    pdf = PdfPages(results_prefix + '.figures.pdf')
    # PDF metadata
    d = pdf.infodict()
    d['Title'] = 'Aligned Variant Diagnostics Plots for {}'.format(args.alignment)
    d['Author'] = 'align_variants.py'

    # Plot GMM diagnostics
    occ_gmm._gmm_plot(M_best, gmms['models'], gmms['data'])
    pdf.attach_note('Residue Occupancy GMM Diagnostics')
    pdf.savefig()
    plt.close()

    # Plot 1
    fig, axs = plt.subplots(1, 2, figsize=(10, 5), sharex=True, sharey=True)
    column_summary.plot.scatter('occupancy', 'missense_variant', ax=axs[0])
    column_summary.plot.scatter('occupancy', 'synonymous_variant', ax=axs[1])
    axs[0].axvline(column_summary[subset_mask_gmm]['occupancy'].min())
    axs[1].axvline(column_summary[subset_mask_gmm]['occupancy'].min())
    fig.suptitle('N Variants vs. Occupancy')
    pdf.attach_note('N Variants vs. Occupancy')
    pdf.savefig()
    plt.close()

    # Conservation plane plot: Variant counts vs. Shenkin
    fig, axs = plt.subplots(1, 2, figsize=(15, 5), sharex=True, sharey=True)
    sns.regplot(x='shenkin', y='missense_variant', data=column_summary[subset_mask_gmm], ax=axs[0])
    pd.plotting.table(axs[0], shenkin_regressions.loc[['missense', 'filtered_missense']].round(2),
                      loc='upper right', colWidths=[0.12] * 5, zorder=100)
    sns.regplot(x='shenkin', y='synonymous_variant', data=column_summary[subset_mask_gmm], ax=axs[1])
    pd.plotting.table(axs[1], shenkin_regressions.loc[['synonymous', 'filtered_synonymous']].round(2),
                      loc='upper right', colWidths=[0.12] * 5, zorder=100)
    plt.title('N Variants vs. Shenkin')
    pdf.attach_note('N Variants vs. Shenkin')
    pdf.savefig()
    plt.close()

    # Conservation plane plot: Missense Scores vs. Shenkin
    plot_data = column_summary[subset_mask_gmm]
    plot_data = plot_data.assign(pass_alpha=plot_data['pvalue'] < 0.1)
    # TODO: plot.scatter throws AttributeError with pandas 0.22.0 or matplotlib 2.1.2
    ax = plot_data.plot.scatter('shenkin', 'oddsratio', c='pass_alpha',  # Valdar is well correlated...
                                colorbar=False,
                                logy=True, figsize=(10, 10))
    _ = plt.setp(ax.get_xticklabels(), visible=True)
    plt.title('Missense Score vs. Shenkin')
    pdf.attach_note('Missense Score vs. Shenkin')
    pdf.savefig()
    plt.close()

    # Other aggregations, some of these just produce the plot
    # Variants per sequence histogram
    protein_consequences = analysis_toolkit._aggregate_annotation(alignment_variant_table, ('VEP', 'Consequence'),
                                                                  aggregate_by=['SOURCE_ID'])
    protein_consequences.hist(facecolor='black', edgecolor='black')
    plt.title('Variants per Sequence')
    pdf.attach_note('Distribution of variants over alignment sequences')
    pdf.savefig()
    plt.close()
    # Variants per residue and column histograms
    fig, axes = plt.subplots(1, 2)
    residue_counts = alignment_variant_table.pipe(analysis_toolkit._aggregate_annotation,
                                                  ('VEP', 'Consequence'),
                                                  aggregate_by=['SOURCE_ID', 'Protein_position'])
    residue_counts = residue_counts.reindex(indexed_mapping_table.index).fillna(0)  # Fill in residues with no variants
    residue_counts['missense_variant'].astype(int).value_counts().plot.bar(ax=axes[0], width=1, facecolor='black',
                                                                           edgecolor='black')
    axes[0].set_title('Missense Variants per Residue')
    # column_variant_counts['missense_variant'].hist(ax=axes[0])
    column_summary.loc[subset_mask_gmm, 'missense_variant'].hist(ax=axes[1], facecolor='black', edgecolor='black')
    axes[1].set_title('Missense Variants per Column')
    pdf.attach_note('Distribution of variants over residues and alignment columns')
    pdf.savefig()
    plt.close()
    pdf.close()

    # Pick extreme columns and identify residues (useful for follow-up)
    umd_mask = column_summary.eval('shenkin_percentile > 0.75 & oddsratio < 1 & pvalue < 0.1')
    ume_mask = column_summary.eval('shenkin_percentile > 0.75 & oddsratio > 1 & pvalue < 0.1')
    cmd_mask = column_summary.eval('shenkin_percentile < 0.25 & oddsratio < 1 & pvalue < 0.1')
    cme_mask = column_summary.eval('shenkin_percentile < 0.25 & oddsratio > 1 & pvalue < 0.1')
    # Get columns from masks
    umd = umd_mask[umd_mask].index
    ume = ume_mask[ume_mask].index
    cmd = cmd_mask[cmd_mask].index
    cme = cme_mask[cme_mask].index
    # Save residues in selection
    indexed_mapping_table.reset_index().set_index(('Alignment', 'Column')).loc[umd].to_csv(
        results_prefix + '.umdres.csv')
    indexed_mapping_table.reset_index().set_index(('Alignment', 'Column')).loc[ume].to_csv(
        results_prefix + '.umeres.csv')
    indexed_mapping_table.reset_index().set_index(('Alignment', 'Column')).loc[cmd].to_csv(
        results_prefix + '.cmdres.csv')
    indexed_mapping_table.reset_index().set_index(('Alignment', 'Column')).loc[cme].to_csv(
        results_prefix + '.cmeres.csv')
    # Write marking jalview tracks
    alignment_column_index = list(range(1, alignment.get_alignment_length() + 1))
    jalview.marked_columns_track(umd_mask.reindex(alignment_column_index, fill_value=False), 'UMD',
                                 'UMD columns at Shenkin PCR > 0.75 and missense OR < 1, p < 0.1',
                                 results_prefix + '.corners.ann')
    jalview.marked_columns_track(ume_mask.reindex(alignment_column_index, fill_value=False), 'UME',
                                 'UME columns at Shenkin PCR > 0.75 and missense OR > 1, p < 0.1',
                                 results_prefix + '.corners.ann', append=True)
    jalview.marked_columns_track(cmd_mask.reindex(alignment_column_index, fill_value=False), 'CMD',
                                 'CMD columns at Shenkin PCR < 0.25 and missense OR < 1, p < 0.1',
                                 results_prefix + '.corners.ann', append=True)
    jalview.marked_columns_track(cme_mask.reindex(alignment_column_index, fill_value=False), 'CME',
                                 'CME columns at Shenkin PCR < 0.25 and missense OR > 1, p < 0.1',
                                 results_prefix + '.corners.ann', append=True)

    # Write variant jalview feature file
    # Label all variants with sequence features
    _write_variants_as_features(alignment_variant_table, results_prefix + '_variant_features.feat')

    # Log completion
    log.info('DONE.')


if __name__ == '__main__':
    parameters = cli()
    main(parameters)
