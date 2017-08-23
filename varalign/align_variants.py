# seems to be required at the top, otherwise get display error...
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import seaborn as sns

import argparse
import logging
from operator import itemgetter

import pandas as pd
from Bio import AlignIO
from numpy import vectorize
import tqdm

import aacon
import alignments
import analysis_toolkit
import ensembl
import gnomad
import jalview
from config import defaults
from gnomad import tabulate_variant_effects
from varalign.analysis_toolkit import _aggregate_annotation
from varalign import occ_gmm

log = logging.getLogger(__name__)
log.setLevel('INFO')


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
        raise ValueError('Could not map {} to the genome.'.format(uniprot))
    # Collapse ranges if desired
    if collapse:
        ensembl_ranges = ensembl.merge_ranges(ensembl_ranges, min_gap=1000)
    return ensembl_ranges


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


def align_variants(alignment, species='HUMAN'):
    """

    :param alignment_info:
    :param species:
    :return:
    """
    # ----- Parse alignment info -----
    log.info('Generating alignment info table...')
    alignment_info = alignments.alignment_info_table(alignment, species)
    log.info('Alignment info table head:\n%s', alignment_info.head().to_string())

    # ----- Map sequences to genome -----
    # TODO: If get transcript ID can use to filter variant table
    genomic_ranges = [
        (row.seq_id, _map_uniprot_to_genome(row.uniprot_id, species=species))
        for row in tqdm.tqdm(alignment_info.itertuples(), total=len(alignment_info))
    ]
    log.info("Mapped {} sequences to genome.".format(len(genomic_ranges)))

    # Add ranges to alignment info
    genomic_mapping_table = pd.DataFrame(genomic_ranges, columns=['seq_id', 'genomic_ranges'])
    alignment_info = alignment_info.merge(genomic_mapping_table, on=['seq_id'], how='left')

    # ----- Fetch variants for the mapped genomic ranges -----
    sequence_variant_lists = [(row.seq_id, (x for r in row.genomic_ranges for x in gnomad.gnomad.fetch(*r)))
                              for row in alignment_info.dropna(subset=['genomic_ranges']).itertuples()]
    all_variants = ((variant, seq_id) for seq_id, range_reader in tqdm.tqdm(sequence_variant_lists)
                    for variant in range_reader)
    variants_table = gnomad.vcf_row_to_table(*zip(*all_variants))

    # ----- Add source UniProt identifiers to the table -----
    # Create UniProt ID series that shares an index with the variant table
    source_uniprot_ids = alignment_info.set_index('seq_id')['uniprot_id']
    source_uniprot_ids.name = ('External', 'SOURCE_ACCESSION')
    source_uniprot_ids.index.name = 'SOURCE_ID'
    # Add IDs to variant tables
    variants_table = variants_table.join(source_uniprot_ids)

    # ----- Filter variant table -----
    filtered_variants = _default_variant_filter(variants_table)
    log.info('Redundant rows:\t{}'.format(sum(filtered_variants.reset_index('Feature').index.duplicated())))
    log.info('Total rows:\t{}'.format(len(filtered_variants)))

    # ----- Map variants to columns -----
    # Generate alignment column / sequence residue mapping table
    indexed_map_table = _mapping_table(alignment_info)
    # Coerce Protein_position to correct type
    filtered_variants.loc[:, ('VEP', 'Protein_position')] = pd.to_numeric(
        filtered_variants.loc[:, ('VEP', 'Protein_position')],
        errors='coerce')
    # Set index for merge
    filtered_variants.reset_index(['SITE', 'ALLELE_NUM', 'Feature'], inplace=True)
    filtered_variants.set_index(('VEP', 'Protein_position'), append=True, inplace=True)
    filtered_variants.index.set_names(['SOURCE_ID', 'Protein_position'], inplace=True)
    filtered_variants.sort_index(inplace=True)
    # Merge to map
    alignment_variant_table = indexed_map_table.join(filtered_variants)  # Drops variants that map outside alignment
    alignment_variant_table.sort_index(inplace=True)

    return alignment_info, alignment_variant_table


if __name__ == '__main__':
    # CLI
    parser = argparse.ArgumentParser(description='Align variants to a Pfam alignment.')
    parser.add_argument('alignment', type=str, help='Path to the alignment.')
    parser.add_argument('--max_gaussians', type=int, default=5, help='Maximum number of Gaussians for occupancy fitting.')
    parser.add_argument('--n_groups', type=int, default=1, help='Top Gaussians to select after occupancy fitting.')
    args = parser.parse_args()

    # Run align variants pipeline
    alignment = AlignIO.read(args.alignment, format='stockholm')
    alignment_info, alignment_variant_table = align_variants(alignment=alignment)
    indexed_mapping_table = _mapping_table(alignment_info)  # TODO: Should this be passed or returned by align_variants?
    # Write data
    alignment_info.to_pickle(args.alignment+'_info.p.gz')
    alignment_variant_table.to_pickle(args.alignment+'_variants.p.gz')

    # Run aacon
    # Format for AACon and run
    aacon_alignment, orig_col_nums = aacon._reformat_alignment_for_aacon(alignment)
    alignment_conservation = aacon._run_aacon(aacon_alignment, orig_col_nums)

    # Save result
    cons_scores_file = 'aacon_scores.csv'
    alignment_conservation.to_csv(cons_scores_file)
    log.info('Formatted AACons results saved to {}'.format(cons_scores_file))

    # Column aggregations
    # Count variants over columns
    column_variant_counts = _aggregate_annotation(alignment_variant_table, ('VEP', 'Consequence'))
    column_variant_counts.to_csv(args.alignment + '.col_var_counts.csv')

    # Count *rare* variants over columns
    rare_maf_threshold = 0.001
    is_rare = alignment_variant_table[('Allele_INFO', 'AF_POPMAX')] < rare_maf_threshold
    column_rare_counts = _aggregate_annotation(alignment_variant_table[is_rare], ('VEP', 'Consequence'))
    column_rare_counts.to_csv(args.alignment + '.col_rare_counts.csv')

    # Count ClinVar annotations for *missense* variants over columns
    is_missense = alignment_variant_table[('VEP', 'Consequence')] == 'missense_variant'
    column_missense_clinvar = _aggregate_annotation(alignment_variant_table[is_missense], ('VEP', 'CLIN_SIG'))
    column_missense_clinvar.to_csv(args.alignment + '.col_mis_clinvar.csv')

    # Count ClinVar annotations for *synonymous* variants over columns
    is_synonymous = alignment_variant_table[('VEP', 'Consequence')] == 'synonymous_variant'
    column_synonymous_clinvar = _aggregate_annotation(alignment_variant_table[is_synonymous], ('VEP', 'CLIN_SIG'))
    column_synonymous_clinvar.to_csv(args.alignment + '.col_syn_clinvar.csv')

    # Use mapping table to calculate human residue occupancy
    # TODO: Adjust for unmapped seqs
    column_occupancy = indexed_mapping_table[('Alignment', 'Column')].value_counts().sort_index()
    column_occupancy.name = 'occupancy'

    # Merge required data for further standard analyses; this is saved after missense scores are added
    column_summary = column_variant_counts.join([column_missense_clinvar, column_occupancy, alignment_conservation])

    # Occupancy GMM
    gmms = occ_gmm._fit_mixture_models(column_summary['occupancy'], args.max_gaussians)
    M_best = occ_gmm._pick_best(gmms['models'], gmms['data'])
    # M_best.means_
    subset_mask_gmm = occ_gmm._core_column_mask(M_best, gmms['data'], args.n_groups)

    # Regression statistics

    # This checks whether missense and synonymous variant counts are correlated with column occupancy before and
    # after column filtering
    variants_vs_occ = analysis_toolkit._comparative_regression(column_summary, subset_mask_gmm, 'occupancy')
    variants_vs_occ.to_csv(args.alignment + '.variant_occ_regression.csv')
    # TODO: Test variants_vs_occ.loc['filtered_missense', 'pvalue'] > 0.05

    # Conservation plane with Shenkin score
    shenkin_regressions = analysis_toolkit._comparative_regression(column_summary, subset_mask_gmm, 'shenkin')
    shenkin_regressions.to_csv(args.alignment + '.variant_shenkin_regression.csv')

    negative_control_p = shenkin_regressions.loc['filtered_synonymous', 'pvalue'] > 0.05
    positive_control_p = shenkin_regressions.loc['filtered_missense', 'pvalue'] < 0.05
    positive_control_m = shenkin_regressions.loc['filtered_missense', 'slope'] > 0
    print 'Filtered synonymous vs. Shenkin (negative control)... {}'.format('PASS' if negative_control_p else 'FAIL')
    print 'Filtered missense vs Shenkin (positive control)... {}'.format(
        'PASS' if positive_control_p and positive_control_m else 'FAIL')

    # Column variant scores, for block columns only
    missense_scores = analysis_toolkit._column_variant_scores(column_summary[subset_mask_gmm],
                                                              variant_class='missense_variant',
                                                              occupancy='occupancy')
    missense_scores.to_csv(args.alignment + '.col_missense_scores.csv')
    column_summary = column_summary.join(missense_scores)
    column_summary.to_csv(args.alignment + '.col_summary.csv')


    # Plot output
    pdf = PdfPages(args.alignment + '.figures.pdf')
    # PDF metadata
    d = pdf.infodict()
    d['Title'] = 'Aligned Variant Diagnostics Plots for {}'.format(args.alignment)
    d['Author'] = 'align_variants.py'

    # Plot GMM diagnostics
    occ_gmm._gmm_plot(gmms['models'], gmms['data'])
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
                      loc='upper right', colWidths=[0.12]*5, zorder=100)
    sns.regplot(x='shenkin', y='synonymous_variant', data=column_summary[subset_mask_gmm], ax=axs[1])
    pd.plotting.table(axs[1], shenkin_regressions.loc[['synonymous', 'filtered_synonymous']].round(2),
                      loc='upper right', colWidths=[0.12]*5, zorder=100)
    plt.title('N Variants vs. Shenkin')
    pdf.attach_note('N Variants vs. Shenkin')
    pdf.savefig()
    plt.close()

    # Conservation plane plot: Missense Scores vs. Shenkin
    plot_data = column_summary[subset_mask_gmm]
    plot_data = plot_data.assign(pass_alpha=plot_data['pvalue'] < 0.1)
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
    protein_consequences = _aggregate_annotation(alignment_variant_table, ('VEP', 'Consequence'),
                                                 aggregate_by=['SOURCE_ID'])
    protein_consequences.hist(facecolor='black', edgecolor='black')
    plt.title('Variants per Sequence')
    pdf.attach_note('Distribution of variants over alignment sequences')
    pdf.savefig()
    plt.close()

    # Variants per residue and column histograms
    fig, axes = plt.subplots(1, 2)
    residue_counts = alignment_variant_table.pipe(_aggregate_annotation,
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
    column_summary = column_summary.join(column_summary.loc[subset_mask_gmm, 'shenkin'].rank(pct=True),
                                         rsuffix='_percentile')
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
    indexed_mapping_table.reset_index().set_index(('Alignment', 'Column')).loc[umd].to_csv(args.alignment+'.umdres.csv')
    indexed_mapping_table.reset_index().set_index(('Alignment', 'Column')).loc[ume].to_csv(args.alignment+'.umeres.csv')
    indexed_mapping_table.reset_index().set_index(('Alignment', 'Column')).loc[cmd].to_csv(args.alignment+'.cmdres.csv')
    indexed_mapping_table.reset_index().set_index(('Alignment', 'Column')).loc[cme].to_csv(args.alignment+'.cmeres.csv')

    # Write marking jalview tracks
    alignment_column_index = range(1, alignment.get_alignment_length() + 1)
    jalview.marked_columns_track(umd_mask.reindex(alignment_column_index, fill_value=False), 'UMD',
                                 'UMD columns at Shenkin PCR > 0.75 and missense OR < 1, p < 0.1',
                                 args.alignment+'.umd.ann')
    jalview.marked_columns_track(ume_mask.reindex(alignment_column_index, fill_value=False), 'UME',
                                 'UME columns at Shenkin PCR > 0.75 and missense OR > 1, p < 0.1',
                                 args.alignment+'.ume.ann')
    jalview.marked_columns_track(cmd_mask.reindex(alignment_column_index, fill_value=False), 'CMD',
                                 'CMD columns at Shenkin PCR < 0.25 and missense OR < 1, p < 0.1',
                                 args.alignment+'.cmd.ann')
    jalview.marked_columns_track(cme_mask.reindex(alignment_column_index, fill_value=False), 'CME',
                                 'CME columns at Shenkin PCR < 0.25 and missense OR > 1, p < 0.1',
                                 args.alignment+'.cme.ann')

