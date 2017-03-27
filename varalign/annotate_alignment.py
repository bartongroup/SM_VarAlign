import argparse
import code
import logging
import math
import os.path

import numpy as np
import pandas as pd
from Bio import AlignIO, SeqIO

from config import defaults
from fetchers import _fetch_variants, select_uniprot_sequence
from jalview import write_jalview_annotation, append_jalview_variant_features, create_jalview_feature_file
from mapping import get_sequence_column_numbers, map_columns_to_res_nums, map_seq_resnums_or_try_isoforms
from stats import run_fisher_tests, calculate_rvis, fill_variant_count
import uniprot
from utils import filter_alignment, is_missense_variant, is_from_to_variant, is_worse_than_type, \
    is_common_variant, is_non_synonomous

log = logging.getLogger(__name__)


def main(alignment, alignment_name, use_local_alignment, local_uniprot_index, downloads):
    """
    Fetch variants for identified protein sequences in an MSA, map to residues and columns and write Jalview feature
    files with key statistics.

    :return:
    """
    # Some parameters
    UniProt_sequences_downloads = os.path.join(downloads, 'UniProt_sequences')
    annotation_suffix = '_jalview_annotations.csv'
    variant_table_suffix = '_alignment_variant_table.csv'

    # Create downloads directory if required
    if not os.path.exists(downloads):
        os.makedirs(downloads)
    jalview_out_file = alignment_name + annotation_suffix

    # Alignment length if used repeatedly to store here
    alignment_length = alignment.get_alignment_length()

    # Map alignment columns to sequence UniProt residue numbers
    uniprot_ids = []
    mapped_records = []
    for seq in alignment:
        # reset holders
        columns = None
        residues = None

        try:
            uniprot_id = seq.annotations['accession']
            uniprot_seq = uniprot.get_uniprot_fasta(uniprot_id)
            seq_name = uniprot_seq.name
        except AttributeError:
            log.warning('Error when retrieving sequence {}; it may be obsolete. Skipping.'.format(seq.id))
            continue

        # Skip unknown sequences
        if uniprot_seq is None:
            log.error('Did not find UniProt sequence corresponding to {}. Skipping.'.format(seq_name))
            continue

        # Map alignment sequence to UniProt sequence
        residues, uniprot_seq = map_seq_resnums_or_try_isoforms(seq, uniprot_seq, use_local_alignment, downloads)
        # Map non-gap column numbers
        columns = get_sequence_column_numbers(seq)

        # Map columns to residues
        if columns and residues:
            col_num_index, sequence_col_nums = columns
            mapped_records.append(map_columns_to_res_nums(sequence_col_nums, residues))
            uniprot_ids.append(uniprot_id)

    # If we skipped all sequences, log and exit
    if len(uniprot_ids) == 0:
        log.info('All sequences filtered or none mapped. Analysis of {} is not applicable. Exiting.'.format(alignment_name))
        return 1

    # Create and concat mapping tables
    log.debug('Tabulating data...')
    sequence_mappings = pd.concat(mapped_records, ignore_index=True)

    # Fetch variants
    protein_identifiers = [uniprot._strip_version(x) for x in uniprot_ids]  # Ensure prots contains UniProt IDs (could be protein names)
    germline_table = _fetch_variants(protein_identifiers, downloads, alignment_name + variant_table_suffix)

    # Merge variant table and key table
    mapped_variants = pd.merge(sequence_mappings, germline_table,
                            left_on=['UniProt_ID', 'uniprot_res_num'],
                            right_on=['UniProt_dbAccessionId', 'start'])

    # TODO: does any of this need deduped?

    # Useful variant masks
    is_missense = is_missense_variant(mapped_variants)
    is_ED = is_from_to_variant('E', 'D', mapped_variants)
    is_DE = is_from_to_variant('D', 'E', mapped_variants)

    # Column variant counts
    total_variant_counts = mapped_variants['alignment_col_num'].value_counts(sort=False)
    missense_variant_counts = mapped_variants.loc[is_missense, 'alignment_col_num'].value_counts(sort=False)
    missense_exc_DE_counts = mapped_variants.loc[is_missense & ~(is_ED | is_DE), 'alignment_col_num'].value_counts(
        sort=False)

    # Column ordered variant counts
    total_variants_per_column = fill_variant_count(total_variant_counts, alignment_length)
    missense_variants_per_column = fill_variant_count(missense_variant_counts, alignment_length)
    missense_exc_DE_per_column = fill_variant_count(missense_exc_DE_counts, alignment_length)

    # Basic Jalview annotation tracks
    variant_counts = [zip(*total_variants_per_column)[1],
                      zip(*missense_variants_per_column)[1],
                      zip(*missense_exc_DE_per_column)[1]]
    titles = ['Total_Variants',
              'Missense_Variants',
              'Missense_Variants (exc. DE)']
    descriptions = ['Total number of variants in summed over all proteins.',
                    'Total number of missense variants in summed over all proteins.',
                    'Number of missense variants excluding E-D and D-E summed over all proteins.']
    # write_jalview_annotation(variant_counts, jalview_out_file, titles, descriptions)
    write_jalview_annotation(variant_counts[1], jalview_out_file, titles[1], descriptions[1])

    # TODO: This is a good scheme for classifying variants, use elsewhere?

    # RVIS calculation: need to count variants and format data.
    is_common = is_common_variant(mapped_variants, maf=None)
    is_bad_type = is_worse_than_type('missense_variant', mapped_variants)
    is_mutant = is_non_synonomous(mapped_variants)
    is_functional = is_bad_type & is_mutant
    common_functional = mapped_variants.loc[is_common & is_functional, 'alignment_col_num'].value_counts(sort=False)
    common_functional_per_column = fill_variant_count(common_functional, alignment_length)

    y = np.array(zip(*common_functional_per_column)[1])
    x = np.array(zip(*total_variants_per_column)[1])
    pred_error, rvis = calculate_rvis(x, y)

    write_jalview_annotation(tuple(pred_error), jalview_out_file, 'Unstandardised RVIS', '', append=True)
    write_jalview_annotation(rvis, jalview_out_file, 'RVIS', '', append=True)

    # Count pathogenic variants if we have at least one unambiguous pathogenic variant.
    if 'pathogenic' in list(mapped_variants['clinical_significance']):
        try:
            clinical_significance_counts = \
                pd.crosstab(mapped_variants.loc[is_missense, 'alignment_col_num'],
                            mapped_variants.loc[is_missense, 'clinical_significance'])
            pathogenic_variant_counts = clinical_significance_counts['pathogenic']
            pathogenic_column_counts = fill_variant_count(pathogenic_variant_counts, alignment_length)
            tooltips = mapped_variants[mapped_variants['clinical_significance'] == 'pathogenic'].groupby('alignment_col_num')
            tooltips = tooltips['seq_id'].aggregate(lambda x: ';'.join(x.unique()))
            tooltips = zip(*fill_variant_count(tooltips, alignment_length))[1]
            write_jalview_annotation(zip(*pathogenic_column_counts)[1], jalview_out_file,
                                     'Pathogenic_missense_variants',
                                     'Number of missense variants annotated pathogenic by ClinVar.', append=True,
                                     tooltips=tooltips)
        except:
            log.warning('Could not count pathogenic variants (possibly there are none).')

        try:
            #TODO: These big try except blocks are a bad idea, especially when actively developing
            # Label pathogenic variants with sequence features
            pathogenic_tables = mapped_variants[(mapped_variants.clinical_significance == 'pathogenic') & \
                is_missense].groupby('seq_id')
            pathogenic_features_file = alignment_name + '_pathogenic_features.txt'
            create_jalview_feature_file({'pathogenic_variant': 'red'}, pathogenic_features_file)
            for seq_id, sub_table in pathogenic_tables:
                residue_indexes = list(sub_table['sequence_index'])
                residue_indexes = [x - 1 + int(seq_id.split('/')[1].split('-')[0]) for x in residue_indexes]
                variant_ids = list(sub_table['variant_id'])
                append_jalview_variant_features(seq_id.split('/')[0], residue_indexes, variant_ids, 'pathogenic_variant',
                                                pathogenic_features_file)
        except:
            log.warning('Could not write pathogenic variants as features (possibly there are none).')

    # Label all variants with sequence features
    missense_tables = mapped_variants[is_missense].groupby('seq_id')
    missense_features_file = alignment_name + '_missense_features.txt'
    create_jalview_feature_file({'missense_variant': 'blue'}, missense_features_file)
    for seq_id, sub_table in missense_tables:
        residue_indexes = list(sub_table['sequence_index'])
        residue_indexes = [x - 1 + int(seq_id.split('/')[1].split('-')[0]) for x in residue_indexes]
        variant_ids = list(sub_table['variant_id'])
        append_jalview_variant_features(seq_id.split('/')[0], residue_indexes, variant_ids,
                                        'missense_variant', missense_features_file)

    # TODO: %gap threshold? Where ignored columns given the worst value for jalview visualisation...

    # Calculate and write fisher test results to Jalview annotation file.
    fisher_test_results = run_fisher_tests(alignment, is_missense, mapped_variants)
    missense_significance = tuple(1 - x for x in zip(*fisher_test_results)[1])
    phred_significance = tuple(-10 * math.log10(x) for x in zip(*fisher_test_results)[1])
    # missense_ratio = tuple(1./x for x in zip(*fisher_test_results)[0])
    write_jalview_annotation(zip(*fisher_test_results)[1], jalview_out_file,
                             'Missense p-value', '', append=True)
    write_jalview_annotation(missense_significance, jalview_out_file,
                             'Missense "significance" (1 - p)', '', append=True)
    write_jalview_annotation(phred_significance, jalview_out_file,
                             'Phred Missense "significance" (1 - p)', '', append=True)
    # write_jalview_annotation(missense_ratio, jalview_out_file,
    #                          'Missense OR', '', append=True)

    return mapped_variants, fisher_test_results, rvis


if __name__ == '__main__':
    # CLI Arguments
    parser = argparse.ArgumentParser(description='Create variant counts column wise for an alignment.')
    parser.add_argument('fasta_file', type=str, help='Path to fasta file containing MSA.')
    parser.add_argument('--use_local_alignment', action='store_true',
                        help='Align sub-sequences to UniProt rather than enforcing exact match.')
    parser.add_argument('--format', type=str, default='fasta',
                        help='Alignment format.')
    parser.add_argument('--seq_id_filter', type=str,
                        help='An inclusive filter to process only a subset of sequences.')
    parser.add_argument('--downloads', type=str, default=defaults.db_root,
                        help='A directory to store downloaded files.')
    parser.add_argument('--interpreter', action='store_true',
                        help='Drop into interactive python session once analysis is complete.')
    parser.add_argument('--local_uniprot_file', type=str,
                        help='Local uniprot flatfile for sequence lookup.')
    parser.add_argument('--write_filtered_alignment', action='store_true',
                        help='Write the alignment containing only the sequences that pass `seq_id_filter`')
    args = parser.parse_args()

    # Initialise variables and parameters
    msa = AlignIO.read(args.fasta_file, args.format)
    msa_name = args.fasta_file
    use_local = args.use_local_alignment
    downloads = args.downloads

    # Initialise local UniProt if provided
    local_uniprot_file = args.local_uniprot_file
    if local_uniprot_file:
        if args.format != 'stockholm':
            log.error('Can only use local UniProt with Stockholm alignments that have AC annotations.')
            raise TypeError
        local_uniprot_index = SeqIO.index(local_uniprot_file, 'swiss')
    else:
        local_uniprot_index = None

    # Filter unwanted sequences if required
    id_filter = args.seq_id_filter
    write_filtered = args.write_filtered_alignment
    if id_filter:
        msa = filter_alignment(msa, id_filter)
    if write_filtered:
        AlignIO.write(msa, msa_name + '_filtered.sto', 'stockholm')

    # Run analysis
    merged_table, fisher_results, rvis_scores = main(msa, msa_name, use_local, local_uniprot_index, downloads)

    # Write results
    merged_table.to_csv(msa_name + '_merged_table.csv')
    scores = pd.DataFrame(fisher_results)
    scores.columns = ['or', 'p']
    scores['rvis'] = pd.Series(rvis_scores)
    scores.to_csv(msa_name + '_scores.csv')

    if args.interpreter:
        code.interact(local=dict(globals(), **locals()))
