"""
This module contains functions that calculate column-wise statistics on aligned variant tables.
"""
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, linregress

import logging

log = logging.getLogger(__name__)


def run_fisher_tests(alignment, table_mask, merged_table):
    """
    Test each column in the alignment for being significantly depleted of variants on a per residue basis.

    Counts are obtained via the cross-table of the 'alignment_col_num' and 'seq_id' columns of `merged_table`.
    This lets us count the number of variants in a column trivially. Counting the number of non-variant residues is
    complicated by gaps in the alignment and non-variant sequences, the latter of which are absent from `merged_table`.

    :param alignment: Alignment to calculate column statistics.
    :param table_mask: A mask to pre-filter the variant table.
    :param merged_table: The merged alignment residue mapping and variant tables table.
    :return:
    """
    log.info('---Running Tests for Conservation---')

    # Collect alignment level counts
    t = merged_table[table_mask]

    # Count sequences that have no variants anywhere
    non_variant_sequences = _non_variant_sequences(alignment, t)
    n_non_variant_sequences = len(non_variant_sequences)

    # Calculate how many positions are in non-variant columns
    alignment_length = alignment.get_alignment_length()
    non_variant_columns = [i for i in range(1, alignment_length + 1) if i not in t['alignment_col_num'].unique()]

    # Count gaps per column
    gaps_per_column = [alignment[:, i].count('-') for i in range(alignment_length)]  # slow

    # Count variants at each column
    cross_table = pd.crosstab(t['seq_id'], t['alignment_col_num'])  # TODO: Only has non-variant residues if protein has variant elsewhere

    # # Drop columns that have a lot of gaps
    # max_gaps = 5
    # for i in range(alignment.get_alignment_length()):
    #     column_string = alignment[:, i]
    #     number_of_gaps = column_string.count('-')
    #     if number_of_gaps > max_gaps and i in cross_table.columns:
    #         cross_table = cross_table.drop(i, axis=1)

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
            log.debug('Counting variants for column {} using cross-table...'.format(col_num))
            variants_in_column = sum(cross_table.loc[:, col_num])
            non_variant_in_column = sum(cross_table.loc[:, col_num] == 0) + n_non_variant_sequences - n_gaps
            variants_in_other = sum(cross_table.drop(col_num, axis=1).sum())
            non_variant_other = sum((cross_table.drop(col_num, axis=1) == 0).sum()) \
                                + n_positions_in_non_variant_columns \
                                + n_positions_in_non_variant_seqs \
                                - n_gaps_other
        else:
            log.debug('Column {} not in cross-table, assuming 0 variants...'.format(col_num))
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


def calculate_rvis(x, y):
    """
    A simple implementation of Petrovski et al's RVIS calculation extended to provide column level scores for an
    alignment.

    :param x: A column ordered numpy.array containing total number of variants for each column.
    :param y: A column ordered numpy.array containing number of 'common functional variants' for each column.
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


def fill_variant_count(value_counts, length):
    """
    Reformat Series.value_counts().

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