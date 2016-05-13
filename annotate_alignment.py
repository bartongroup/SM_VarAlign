import argparse
import code
import logging
import math
import os.path

import numpy as np
import pandas as pd
from Bio import AlignIO, SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.pairwise2 import format_alignment

from config import defaults
from fetchers import fetch_uniprot_sequences, _fetch_variants
from jalview_writers import write_jalview_annotation, append_jalview_variant_features, create_jalview_feature_file
from stats import run_fisher_tests, calculate_rvis, fill_variant_count
from utils import worse_than, parse_seq_name, filter_alignment

log = logging.getLogger(__name__)


def get_row_residue_numbers(subseq, uniprot_seq, use_local_alignment):
    """
    Map each sequence in an alignment to a longer sequence and return the residue numbers.

    :param align: MultipleSeqAlignment for processing.
    :param uniprot_seq: List of full-length SeqRecords.
    :param use_local_alignment: If True, use a local alignment rather than exact pattern matching.
    :return: A list of tuples with the sequence IDs and residues numbers for each sequence in the alignment.
    """
    # Prepare sequences for alignment to UniProt by removing gaps
    subseq = SeqRecord(Seq(str(subseq.seq).replace('-', '').upper(), subseq.seq.alphabet), id=subseq.id)
    sequence_name = parse_seq_name(subseq.id)  # TODO: No longer needed

    if use_local_alignment:
        # Align input alignment sequences to UniProt Sequences
        log.debug('Aligning sub-sequence: {} to UniProt sequence: {}'.format(subseq.id, uniprot_seq[1].id))
        local_alignment = pairwise2.align.localxs(subseq.seq, uniprot_seq[1].seq, -.5, -.1)
        for pairwise_alignment in local_alignment:
            log.debug('{}'.format('\n' + format_alignment(*pairwise_alignment)))
        alignment = subseq.id, uniprot_seq[1].id, local_alignment

        # Build list of UniProt residue numbers for each non-gap for each sequence
        sub_seq_id, uniprot_seq_id, pairwise = alignment
        seq = str(pairwise[0][0])
        res_nums = [i + 1 for i, s in enumerate(seq) if s != '-']  # TODO: wrong if other seq has gaps too
        seq_indexes = range(1, len(res_nums) + 1)
        align_res_nums = sub_seq_id, uniprot_seq_id, res_nums, seq_indexes
    else:
        str_seq = str(subseq.seq)
        uniprot_str_seq = str(uniprot_seq[1].seq)
        if str_seq in uniprot_str_seq:
            log.debug('Matching sub-sequence: {} to UniProt sequence: {}'.format(subseq.id, uniprot_seq[1].id))
            start = uniprot_str_seq.find(str_seq) + 1
            end = start + len(str_seq)
            res_nums = range(start, end)
            seq_indexes = range(1, len(res_nums) + 1)
            align_res_nums = subseq.id, uniprot_seq[1].id, res_nums, seq_indexes
        else:
            log.error('Sub-sequence: {} does not match UniProt sequence: {}'.format(subseq.id, uniprot_seq[1].id))
            log.error('Sub-sequence: {} does not match UniProt sequence: {}'.format(str_seq, uniprot_str_seq))
            raise TypeError

    return align_res_nums


def get_sequence_column_numbers(sequence):
    """
    Build list of column numbers for each non-gap for each sequence.

    :param sequence: Aligned sequence.
    :return:
    """
    seq_id = sequence.id
    seq = str(sequence.seq)
    col_nums = [i + 1 for i, s in enumerate(seq) if s != '-']
    align_col_nums = seq_id, col_nums

    return align_col_nums


def map_columns_to_residues(alignment_column_numbers, alignment_residue_numbers):
    """
    Map alignment columns to UniProt residue numbers.

    :param alignment_column_numbers:
    :param alignment_residue_numbers:
    :return:
    """
    col_num_index, sequence_col_nums = zip(*alignment_column_numbers)
    mapped = []
    for seq_id, uniprot_seq_id, res_nums, seq_index in alignment_residue_numbers:
        log.debug('Mapping {}...'.format(seq_id))
        col_nums = sequence_col_nums[col_num_index.index(seq_id)]  # Lookup sequence column numbers
        record = {'seq_id': seq_id, 'uniprot_seq_id': uniprot_seq_id, 'uniprot_res_num': res_nums,
                  'alignment_col_num': col_nums, 'sequence_index': seq_index}
        # TODO: This is messy due to dependency on uniprot_seq_id format in `alignment_residue_numbers`
        if '|' in record['uniprot_seq_id']:
            prot_name = record['uniprot_seq_id'].split('|')[1]  # UniProt ID
        else:
            prot_name = record['uniprot_seq_id']
        record.update({'UniProt_ID': prot_name})

        mapped.append(pd.DataFrame(record))

    # Create and concat mapping tables
    log.debug('Tabulating data...')
    mapped_df = pd.concat(mapped, ignore_index=True)

    return mapped_df


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
    uniprot_sequences = []
    alignment_residue_numbers = []
    alignment_column_numbers = []
    for seq in alignment:
        # Identify sequence and retrieve full UniProt
        seq_name = parse_seq_name(seq.id)
        if not local_uniprot_index:
            uniprot_seq = fetch_uniprot_sequences(seq_name, UniProt_sequences_downloads)
        else:
            # TODO: Currently local lookup only working with Stockholm format that has AC annotations
            accession_code = seq.annotations['accession'].split('.')[0]  # Dropping sequence version
            if accession_code in local_uniprot_index:
                uniprot_seq = accession_code, local_uniprot_index[accession_code]
            else:
                uniprot_seq = None

        # Skip unknown sequences
        if uniprot_seq is None:
            log.error('Did not find UniProt sequence corresponding to {}. Skipping.'.format(seq_name))
            continue
        uniprot_sequences.append(uniprot_seq)  # Keep for later too

        # Map alignment sequence to UniProt sequence
        try:
            alignment_residue_numbers.append(get_row_residue_numbers(seq, uniprot_seq, use_local_alignment))
        except TypeError:
            # Maybe it's a different isoform
            canonical_uniprot = uniprot_seq[1].id.split('|')[1]
            for suffix in ('-2', '-3'):
                try:
                    isoform = canonical_uniprot + suffix
                    print isoform
                    uniprot_seq = fetch_uniprot_sequences(isoform, downloads)
                    alignment_residue_numbers.append(
                        get_row_residue_numbers(seq, uniprot_seq, use_local_alignment))
                    break
                except TypeError:
                    continue

        # Map non-gap column numbers
        alignment_column_numbers.append(get_sequence_column_numbers(seq))

    # If we skipped all sequences, log and exit
    if len(uniprot_sequences) == 0:
        log.info('All sequences filtered or none mapped. Analysis of {} is not applicable. Exiting.'.format(alignment_name))
        return 1

    # Map columns to residues
    mapped = map_columns_to_residues(alignment_column_numbers, alignment_residue_numbers)

    # Fetch variants
    protein_identifiers = zip(*uniprot_sequences)[0]  # Ensure prots contains UniProt IDs (could be protein names)
    germline_table = _fetch_variants(protein_identifiers, downloads, alignment_name + variant_table_suffix)

    # Merge variant table and key table
    merged_table = pd.merge(mapped, germline_table,
                            left_on=['UniProt_ID', 'uniprot_res_num'],
                            right_on=['UniProt_dbAccessionId', 'start'])

    # TODO: does any of this need deduped?

    # Useful variant masks
    is_missense = (merged_table['type'] == 'missense_variant') & \
                  (merged_table['from_aa'] != merged_table['to_aa_expanded'])
    is_ED = (merged_table['from_aa'] == 'E') & (merged_table['to_aa_expanded'] == 'D')
    is_DE = (merged_table['from_aa'] == 'D') & (merged_table['to_aa_expanded'] == 'E')

    # Column variant counts
    total_variant_counts = merged_table['alignment_col_num'].value_counts(sort=False)
    missense_variant_counts = merged_table.loc[is_missense, 'alignment_col_num'].value_counts(sort=False)
    missense_exc_DE_counts = merged_table.loc[is_missense & ~(is_ED | is_DE), 'alignment_col_num'].value_counts(
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
    is_common = merged_table['minor_allele_frequency'].notnull()
    is_bad_type = merged_table.type.apply(lambda x: x in worse_than('missense_variant'))
    is_mutant = merged_table['from_aa'] != merged_table['to_aa_expanded']
    is_functional = is_bad_type & is_mutant
    common_functional = merged_table.loc[is_common & is_functional, 'alignment_col_num'].value_counts(sort=False)
    common_functional_per_column = fill_variant_count(common_functional, alignment_length)

    y = np.array(zip(*common_functional_per_column)[1])
    x = np.array(zip(*total_variants_per_column)[1])
    pred_error, rvis = calculate_rvis(x, y)

    write_jalview_annotation(tuple(pred_error), jalview_out_file, 'Unstandardised RVIS', '', append=True)
    write_jalview_annotation(rvis, jalview_out_file, 'RVIS', '', append=True)

    # Count pathogenic variants if we have at least one unambiguous pathogenic variant.
    if 'pathogenic' in list(merged_table['clinical_significance']):
        try:
            clinical_significance_counts = \
                pd.crosstab(merged_table.loc[is_missense, 'alignment_col_num'],
                            merged_table.loc[is_missense, 'clinical_significance'])
            pathogenic_variant_counts = clinical_significance_counts['pathogenic']
            pathogenic_column_counts = fill_variant_count(pathogenic_variant_counts, alignment_length)
            tooltips = merged_table[merged_table['clinical_significance'] == 'pathogenic'].groupby('alignment_col_num')
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
            pathogenic_tables = merged_table[(merged_table.clinical_significance == 'pathogenic') & \
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

    # TODO: %gap threshold? Where ignored columns given the worst value for jalview visualisation...

    # Calculate and write fisher test results to Jalview annotation file.
    fisher_test_results = run_fisher_tests(alignment, is_missense, merged_table)
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

    return merged_table, fisher_test_results, rvis


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

    # Read alignment, initialise variables and run analysis
    msa = AlignIO.read(args.fasta_file, args.format)
    msa_name = args.fasta_file
    id_filter = args.seq_id_filter
    use_local = args.use_local_alignment
    local_uniprot_file = args.local_uniprot_file
    if local_uniprot_file:
        if args.format != 'stockholm':
            log.error('Can only use local UniProt with Stockholm alignments that have AC annotations.')
            raise TypeError
        local_uniprot_index = SeqIO.index(local_uniprot_file, 'swiss')
    else:
        local_uniprot_index = None
    write_filtered = args.write_filtered_alignment
    downloads = args.downloads

    # Filter unwanted sequences
    if id_filter:
        msa = filter_alignment(msa, id_filter)
    if write_filtered:
        AlignIO.write(msa, msa_name + '_filtered.sto', 'stockholm')

    merged_table, fisher_results, rvis_scores = main(msa, msa_name, use_local, local_uniprot_index, downloads)

    merged_table.to_csv(msa_name + '_merged_table.csv')
    scores = pd.DataFrame(fisher_results)
    scores.columns = ['or', 'p']
    scores['rvis'] = pd.Series(rvis_scores)
    scores.to_csv(msa_name + '_scores.csv')

    if args.interpreter:
        code.interact(local=dict(globals(), **locals()))
