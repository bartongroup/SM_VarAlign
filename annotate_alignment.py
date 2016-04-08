import argparse
import code
import os.path
import re
import urllib2

import numpy as np
import pandas as pd
from Bio import AlignIO, SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.pairwise2 import format_alignment
from scipy.stats import fisher_exact, linregress

from utils import urlopen_with_retry, query_uniprot, worse_than
# Use my developement branch of ProteoFAV
import sys

sys.path.extend(['/Users/smacgowan/PycharmProjects/ProteoFAV'])
from proteofav.variants import select_uniprot_variants
from proteofav.analysis.utils import expand_dataframe

import logging

log = logging.getLogger(__name__)


def parse_seq_name(seq_name):
    return re.search('\w*', seq_name).group().strip()


def _non_variant_sequences(alignment, variant_table):
    """
    Identify sequences in an aligment that have no corresponding entry in a variant table.

    :param alignment: Bio.Align.MultipleSeqAlignment
    :param variant_table: Pandas.DataFrame variant table
    :return: List of sequence IDs not found in the variant table
    """
    # Identify missing sequences because entirely non-variant
    all_sequence_ids = [a.id for a in alignment]
    sequences_with_variants = list(variant_table['seq_id'].unique())
    non_variant_sequences = [a for a in all_sequence_ids if a not in sequences_with_variants]
    return non_variant_sequences


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
        align_res_nums = sub_seq_id, uniprot_seq_id, res_nums
    else:
        str_seq = str(subseq.seq)
        uniprot_str_seq = str(uniprot_seq[1].seq)
        if str_seq in uniprot_str_seq:
            log.debug('Matching sub-sequence: {} to UniProt sequence: {}'.format(subseq.id, uniprot_seq[1].id))
            start = uniprot_str_seq.find(str_seq) + 1
            end = start + len(str_seq)
            res_nums = range(start, end)
            align_res_nums = subseq.id, uniprot_seq[1].id, res_nums
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


def fetch_uniprot_sequences(seq_name, downloads=None):
    """
    Retrieve UniProt sequences.

    :param protein_identifiers: List of protein identifiers (UniProt IDs or protein name)
    :return: List of protein sequences.
    """
    url = 'http://www.uniprot.org/uniprot/'
    p = seq_name.strip()
    fasta_file_name = os.path.join(downloads, p + '.fasta')
    remote_fasta = url + p + '.fasta'
    if not os.path.isfile(fasta_file_name):
        print remote_fasta
        try:
            handle = urlopen_with_retry(remote_fasta)
        except urllib2.HTTPError:
            # Will need to query instead
            p = parse_seq_name(p)  # First word only
            # TODO: This should be configurable
            p = query_uniprot(('gene:' + p, 'reviewed:yes', 'organism:human'), first=True)
            remote_fasta = url + p + '.fasta'
            handle = urlopen_with_retry(remote_fasta)
        seq_record = SeqIO.read(handle, "fasta")
        if downloads is not None:
            if not os.path.exists(downloads):
                os.makedirs(downloads)
            SeqIO.write(seq_record, fasta_file_name, "fasta")
    else:
        handle = open(fasta_file_name, 'r')
        seq_record = SeqIO.read(handle, "fasta")
    p = seq_record.id.split('|')[1]  # Extract UniProt ID
    uniprot_sequence = p, seq_record
    return uniprot_sequence


def map_columns_to_residues(alignment_column_numbers, alignment_residue_numbers):
    """
    Map alignment columns to UniProt residue numbers.

    :param alignment_column_numbers:
    :param alignment_residue_numbers:
    :return:
    """
    mapped = []
    for seq_id, uniprot_seq_id, res_nums in alignment_residue_numbers:
        ind = zip(*alignment_column_numbers)[0].index(seq_id)
        col_nums = zip(*alignment_column_numbers)[1][ind]
        mapped.append({'seq_id': seq_id, 'uniprot_seq_id': uniprot_seq_id, 'uniprot_res_num': res_nums,
                       'alignment_col_num': col_nums})

    for i in mapped:
        prot_name = i['uniprot_seq_id'].split('|')[1]  # UniProt ID
        i.update({'UniProt_ID': prot_name})

    # Create and concat mapping tables
    mapped_df = pd.DataFrame()
    for i in mapped:
        mapped_df = mapped_df.append(pd.DataFrame(i), ignore_index=True)

    return mapped_df


def _fetch_variants(prots, downloads=None, save_name=None):
    """
    
    :param prots:
    :param downloads:
    :param save_name:
    :return:
    """
    # Get variant data
    # Get the data with EnsEMBL variants
    table_file_name = os.path.join(downloads, save_name)
    if not os.path.isfile(table_file_name):
        tables = []
        for p in list(set(prots)):
            try:
                variant_table = select_uniprot_variants(p, reduced_annotations=False)  # TODO: Use new variant fetcher?
                variant_table['UniProt_dbAccessionId'] = p
                tables.append(variant_table)
            except (ValueError, KeyError):
                log.error('Could not retrieve variants for {}'.format(p))

        # Concatenate and process all those variant tables
        concat_table = pd.concat(tables, ignore_index=True)
        # Need to expand on 'to_aa' before dedupping
        concat_table['orig_index'] = concat_table.index
        concat_table = expand_dataframe(df=concat_table, expand_column='to_aa', id_column='orig_index')
        concat_table = concat_table.drop('orig_index', 1)
        # Fix or remove list columns
        concat_table = concat_table.drop('to_aa', 1)
        concat_table['clinical_significance'] = concat_table['clinical_significance'].apply(lambda x: ';'.join(x))
        concat_table['clinical_significance'].fillna('')
        # And dedup, bearing in mind the same variant can pop up in different transcripts
        # (so dedupping is only done on certain columns)
        concat_table = concat_table.drop_duplicates(['UniProt_dbAccessionId',
                                                     'start',
                                                     'end',
                                                     'variant_id',
                                                     'to_aa_expanded']).reset_index(drop=True)

        # Write table to file
        concat_table.to_csv(table_file_name)
    else:
        concat_table = pd.read_csv(table_file_name)

    # is_somatic = concat_table['variant_id'].apply(lambda x: x.startswith('COS'))  #TODO: include this?
    is_germline = concat_table['variant_id'].apply(lambda x: x.startswith('rs'))
    # somatic_table = concat_table[is_somatic]
    germline_table = concat_table[is_germline]

    return germline_table


def fill_variant_count(value_counts, length):
    """
    Order an alignment column number value counts Series by alignment column and insert 0s for any unobserved columns.

    :param value_counts: `alignment_col_num`.value_counts() Series
    :param length: The length of the alignment.
    :return:
    """
    variants_per_pos = []
    for i in xrange(length):
        col_pos = i + 1
        try:
            variants_per_pos.append((col_pos, value_counts[col_pos]))
        except KeyError:
            variants_per_pos.append((col_pos, 0))
    return variants_per_pos


def write_jalview_annotation(ordered_values, file_name, title, description, append=False):
    """
    Write data to a Jalview annotation file.

    :param ordered_values: Tuple or list of tuples containing annotation scores for each column in alignment.
    :param file_name: Filename for Jalview annotations
    :param title: Tuple or list of tuples of titles for each annotation track.
    :param description: As above containing descriptions.
    :param append: Wheter to append to an existing Jalview annotation file.
    :return:
    """
    file_mode = 'w'
    if append:
        file_mode = 'a'

    with open(file_name, file_mode) as results_file:
        if not append:
            results_file.write('JALVIEW_ANNOTATION\n')
        if isinstance(ordered_values, tuple) and all(map(lambda x: isinstance(x, str), [title, description])):
            results_file.write('BAR_GRAPH\t{}\t{}\t'.format(title, description) +
                               '|'.join('{},,{}'.format(str(x), str(x)) for x in ordered_values))
            results_file.write('\n')
        elif all(map(lambda x: isinstance(x, list), [ordered_values, title, description])):
            arg_lengths = map(len, [ordered_values, title, description])
            if len(set(arg_lengths)) == 1:
                for v, t, d in zip(ordered_values, title, description):
                    results_file.write('BAR_GRAPH\t{}\t{}\t'.format(t, d) +
                                       '|'.join('{},,{}'.format(str(x), str(x)) for x in v))
                    results_file.write('\n')
            else:
                log.error('List arguments must be of same length')
                raise TypeError
        else:
            log.error('Must provide same number of titles/descriptions as lists of values.')
            raise TypeError

    return 0


def run_fisher_tests(alignment, table_mask, merged_table):
    """
    Test each column in the alignment for being statistically significantly depleted of variants on a per residue basis.

    :param alignment: Original alignment to calculate column statistics.
    :param table_mask: A mask to pre-filter the variant table.
    :param merged_table: The merged alignment residue mapping and variant tables table.
    :return:
    """
    log.info('---Running Tests for Conservation---')

    # Collect alignment level counts
    t = merged_table[table_mask]
    cross_table = pd.crosstab(t['seq_id'], t['alignment_col_num']) #TODO: Only has non-variant residues if protein has variant elsewhere

    # # Drop columns that have a lot of gaps
    # max_gaps = 5
    # for i in range(alignment.get_alignment_length()):
    #     column_string = alignment[:, i]
    #     number_of_gaps = column_string.count('-')
    #     if number_of_gaps > max_gaps and i in cross_table.columns:
    #         cross_table = cross_table.drop(i, axis=1)
    # Count sequences that have no variants anywhere

    non_variant_sequences = _non_variant_sequences(alignment, t)
    n_non_variant_sequences = len(non_variant_sequences)

    # Calculate how many positions are in non-variant columns
    alignment_length = alignment.get_alignment_length()
    non_variant_columns = [i for i in range(1, alignment_length + 1) if i not in cross_table.columns]

    # Count gaps per column
    gaps_per_column = [alignment[:, i].count('-') for i in range(alignment_length)]  # slow

    # Run fisher tests for all columns
    fisher_test_results = []
    for col_num in range(alignment_length):
        # TODO: this doesn't account for number of gaps in a column
        # TODO: also doesn't account for sequences with no variants
        # TODO: Might be double counting residues that are in both non_variant_columns and non_variant_sequences

        col_num += 1  # aligment_col_num is 1-indexed

        # Count gaps
        n_gaps = gaps_per_column[col_num - 1]
        other_columns = range(alignment_length)
        other_columns.remove(col_num - 1)
        n_gaps_other = sum([gaps_per_column[i] for i in other_columns])

        # Calculate positions in other non_variant columns and sequences
        if col_num not in non_variant_columns:
            n_positions_in_non_variant_columns = len(non_variant_columns) * len(alignment)
        else:
            n_positions_in_non_variant_columns = (len(non_variant_columns) - 1) * len(alignment)
        n_positions_in_non_variant_seqs = n_non_variant_sequences * (alignment_length - (1 + len(non_variant_columns)))

        # # Count non-variant sequence residues in and not in column
        # non_variant_sub_alignment = [str(a.seq) for a in sub_alignment if a.id not in sequences_with_variants]
        # n_other_residues_non_variant_seq = sum([len(a) for a in non_variant_sub_alignment])

        # Count variants
        if col_num in cross_table.columns:
            variants_in_column = sum(cross_table.loc[:, col_num])
            non_variant_in_column = sum(cross_table.loc[:, col_num] == 0) + n_non_variant_sequences - n_gaps
            variants_in_other = sum(cross_table.drop(col_num, axis=1).sum())
            non_variant_other = sum((cross_table.drop(col_num, axis=1) == 0).sum()) \
                                + n_positions_in_non_variant_columns \
                                + n_positions_in_non_variant_seqs \
                                - n_gaps_other
        else:
            variants_in_column = 0
            non_variant_in_column = len(alignment) - n_gaps
            variants_in_other = sum(cross_table.sum())
            non_variant_other = sum((cross_table == 0).sum()) \
                                + n_positions_in_non_variant_columns \
                                + n_positions_in_non_variant_seqs \
                                - n_gaps_other

        # Calculate OR and p-value
        # print (variants_in_column, variants_in_other)
        # print (non_variant_in_column, non_variant_other)
        odds_ratio, pvalue = fisher_exact([[variants_in_column, variants_in_other],
                                           [non_variant_in_column, non_variant_other]],
                                          alternative='less')
        fisher_test_results.append((odds_ratio, pvalue))
        log.info('Alignment column: {}, OR = {}, p = {}'.format(col_num, odds_ratio, pvalue))

    return fisher_test_results


def calculate_rvis(x, y):
    """

    :param x:
    :param y:
    :return:
    """
    # RVIS approx.
    slope, intercept, r_value, p_value, slope_std_error = linregress(x, y)
    predict_y = intercept + slope * x
    pred_error = y - predict_y
    rvis = tuple(pred_error / np.std(pred_error))

    # # Proper RVIS
    # x = x.reshape((len(x),1))
    # y = y.reshape((len(y),1))
    # regr = linear_model.LinearRegression()
    # regr.fit(x, y)
    # rvis_int_stud = residuals(regr, x, y, 'standardized')  # Different format...
    # rvis_int_stud = tuple(rvis_int_stud.reshape((len(rvis_int_stud), )))
    # rvis_ext_stud = tuple(residuals(regr, x, y, 'studentized'))
    # write_jalview_annotation(rvis_int_stud, jalview_out_file, 'Int. Stud. RVIS', '', append=True)
    # write_jalview_annotation(rvis_ext_stud, jalview_out_file, 'Ext. Stud. RVIS', '', append=True)

    return pred_error, rvis


def main(args):
    """
    Fetch variants for identified protein sequences in an MSA, map to residues and columns and write Jalview feature
    files with key statistics.

    :return:
    """
    # Some parameters
    downloads = '.VarAlign'
    UniProt_sequences_downloads = os.path.join(downloads, 'UniProt_sequences')
    if not os.path.exists(downloads):
        os.makedirs(downloads)
    jalview_out_file = args.fasta_file + '_jalview_annotations.csv'

    # Read alignment
    alignment = AlignIO.read(args.fasta_file, "fasta")
    alignment_length = alignment.get_alignment_length()

    # Map alignment columns to sequence UniProt residue numbers
    uniprot_sequences = []
    alignment_residue_numbers = []
    alignment_column_numbers = []
    for seq in alignment:
        if args.seq_id_filter is not None and args.seq_id_filter not in seq.id:
            log.info('Filtering sequence {}.'.format(seq.id))
            continue
        seq_name = parse_seq_name(seq.id)
        try:
            uniprot_seq = fetch_uniprot_sequences(seq_name, UniProt_sequences_downloads)
        except ValueError:
            log.error('Could not retrieve sequence for {}, skipping...'.format(seq_name))
            continue
        uniprot_sequences.append(uniprot_seq)  # Keep for later too
        try:
            alignment_residue_numbers.append(get_row_residue_numbers(seq, uniprot_seq, args.use_local_alignment))
        except TypeError:
            # Maybe it's a different isoform
            canonical_uniprot = uniprot_seq[1].id.split('|')[1]
            for suffix in ('-2', '-3'):
                try:
                    isoform = canonical_uniprot + suffix
                    print isoform
                    uniprot_seq = fetch_uniprot_sequences(isoform, downloads)
                    alignment_residue_numbers.append(
                        get_row_residue_numbers(seq, uniprot_seq, args.use_local_alignment))
                    break
                except TypeError:
                    continue
        alignment_column_numbers.append(get_sequence_column_numbers(seq))
    mapped = map_columns_to_residues(alignment_column_numbers, alignment_residue_numbers)  # Map columns to residues

    # Fetch variants
    protein_identifiers = zip(*uniprot_sequences)[0]  # Ensure prots contains UniProt IDs (could be protein names)
    germline_table = _fetch_variants(protein_identifiers, downloads, args.fasta_file + '_alignment_variant_table.csv')

    # Merge variant table and key table
    merged_table = pd.merge(mapped, germline_table,
                            left_on=['UniProt_ID', 'uniprot_res_num'],
                            right_on=['UniProt_dbAccessionId', 'start'])
    merged_table.to_csv(args.fasta_file + '_merged_table.csv')

    # Counting variants and writing Jalview annotations
    # TODO: does any of this need deduped?
    is_missense = (merged_table['type'] == 'missense_variant') & \
                  (merged_table['from_aa'] != merged_table['to_aa_expanded'])
    is_ED = (merged_table['from_aa'] == 'E') & (merged_table['to_aa_expanded'] == 'D')
    is_DE = (merged_table['from_aa'] == 'D') & (merged_table['to_aa_expanded'] == 'E')
    total_variant_counts = merged_table['alignment_col_num'].value_counts(sort=False)
    total_variants_per_column = fill_variant_count(total_variant_counts, alignment_length)
    missense_variant_counts = merged_table.loc[is_missense, 'alignment_col_num'].value_counts(sort=False)
    missense_variants_per_column = fill_variant_count(missense_variant_counts, alignment_length)
    missense_exc_DE_counts = merged_table.loc[is_missense & ~(is_ED | is_DE), 'alignment_col_num'].value_counts(
        sort=False)
    missense_exc_DE_per_column = fill_variant_count(missense_exc_DE_counts, alignment_length)

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

    # RVIS calculation: need to count variants and format data.
    # TODO: This is a good scheme for classifying variants, use elsewhere?
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
            write_jalview_annotation(zip(*pathogenic_column_counts)[1], jalview_out_file,
                                     'Pathogenic_missense_variants',
                                     'Number of missense variants annotated pathogenic by ClinVar.', append=True)
        except:
            log.warning('Count not count pathogenic variants.')

    # Calculate and write fisher test results to Jalview annotation file.
    # TODO: This and other calcs could be run with %gap threshold, ignored columns given the worst value for jalview visualisation
    fisher_test_results = run_fisher_tests(alignment, is_missense, merged_table)
    missense_significance = tuple(1 - x for x in zip(*fisher_test_results)[1])
    # missense_ratio = tuple(1./x for x in zip(*fisher_test_results)[0])
    write_jalview_annotation(zip(*fisher_test_results)[1], jalview_out_file,
                             'Missense p-value', '', append=True)
    write_jalview_annotation(missense_significance, jalview_out_file,
                             'Missense "sginificance" (1 - p)', '', append=True)
    # write_jalview_annotation(missense_ratio, jalview_out_file,
    #                          'Missense OR', '', append=True)

    return merged_table, fisher_test_results


if __name__ == '__main__':
    # CLI Arguments
    parser = argparse.ArgumentParser(description='Create variant counts column wise for an alignment.')
    parser.add_argument('fasta_file', type=str, help='Path to fasta file containing MSA.')
    parser.add_argument('--use_local_alignment', action='store_true',
                        help='Align sub-sequences to UniProt rather than enforcing exact match.')
    parser.add_argument('--seq_id_filter', type=str,
                        help='An inclusive filter to process only a subset of sequences.')
    parser.add_argument('--interpreter', action='store_true',
                        help='Drop into interactive python session once analysis is complete.')
    args = parser.parse_args()

    merged_table, fisher_test_results = main(args)

    if args.interpreter:
        code.interact(local=dict(globals(), **locals()))
