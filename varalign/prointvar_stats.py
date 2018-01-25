import numpy as np
import pandas as pd
from scipy import stats


def _structure_column_counts(aligned_prointvar_table, query, unique_sequences_name, total_interactions_name):
    """

    :param aligned_prointvar_table:
    :param query:
    :param unique_proteins_name:
    :param total_interactions_name:
    :return:
    """
    # Subset structure table with query if specified
    if query:
        sub_table = aligned_prointvar_table.query(query)
    else:
        sub_table = aligned_prointvar_table

    # TODO: always do this for both ATOM_A and _B? only when contacts is lower diagonal form as I use? make optional?
    # Group interactions at each column and summarise
    gr_a = sub_table.groupby(['SOURCE_ID_A', 'Alignment_column_A']).size().groupby('Alignment_column_A')
    n_unique_sequences = gr_a.size()  # How many unique sequences have interactions?
    n_interactions = gr_a.sum()  # How many interactions are there in total?

    # Repeat for ATOM_B
    # Group interactions at each column and summarise
    gr_b = sub_table.groupby(['SOURCE_ID_B', 'Alignment_column_B']).size().groupby('Alignment_column_B')
    n_unique_sequences = n_unique_sequences.add(gr_b.size(), fill_value=0)
    n_interactions = n_interactions.add(gr_b.sum(), fill_value=0)

    # Format results
    n_unique_sequences.name = unique_sequences_name
    n_interactions.name = total_interactions_name

    # Return as DataFrame
    return pd.concat([n_unique_sequences, n_interactions], axis=1).fillna(0)


def _column_ligand_contacts(aligned_prointvar_table):
    """
    Calculate frequency of ligand interactions over alignment columns. Gives both total counts and number of sequences
    with evidence of a ligand interaction.

    :param aligned_prointvar_table:
    :return:
    """
    return _structure_column_counts(aligned_prointvar_table, query='interaction_type == "Protein-Ligand"',
                                    unique_sequences_name='protein_ligand_interactions',
                                    total_interactions_name='total_ligand_interactions')


def _column_total_contacts(aligned_prointvar_table):
    """
    Calculate total number of interactions observed at each alignment column. Gives both total counts and number of
    sequences with at least one interaction.

    :param aligned_prointvar_table:
    :return:
    """
    return _structure_column_counts(aligned_prointvar_table, query=None,
                                    unique_sequences_name='sequences_with_contacts',
                                    total_interactions_name='total_contacts')


def _column_protein_contacts(aligned_prointvar_table):
    """
    Calculate frequency of protein-protein interactions over alignment columns. Gives both total counts and number of
    sequences with evidence of a protein-protein interaction.

    :param aligned_prointvar_table:
    :return:
    """
    return _structure_column_counts(aligned_prointvar_table,
                                    query='interaction_type == "Protein-Protein" & polymer_topology == "Interpolymer"',
                                    unique_sequences_name='protein_protein_interactions',
                                    total_interactions_name='total_protein_interactions')


def _interpolate_index(table):
    """
    Re-index an interger indexed table to include missing values within the range of the original index.

    :param table:
    :return:
    """
    idx = pd.Index(range(int(table.index.min()), int(table.index.max())))
    return table.reindex(idx)


def _count_nan(table):
    """
    Count the number of NaNs in all columns of a DataFrame.

    :param table:
    :return:
    """
    return {col: table.eval('{0} != {0}'.format(col)).sum() for col in table}


def _column_set_feature_enrichment(column_annotation_table, selection_query, feature_column, feature_background):
    """
    Test a subset of a columns for enrichment of a feature compared to the remaining columns.

    :param column_annotation_table: DataFrame
    :param selection_query: Query string to identify subset from `column_annotation_table`.
    :param feature_column: Column name containing feature to test for enrichment.
    :param feature_background: Column name containing background for the test.
    :return:
    """
    # generate 2x2 contingency table
    groups = column_annotation_table.eval(selection_query)
    test_df = column_annotation_table.groupby(groups).sum()[[feature_column, feature_background]]
    test_df.eval('background = {} - {}'.format(feature_background, feature_column), inplace=True)  # TODO: make option?
    test_df = test_df.reindex([0, 1]).fillna(0)  # Deal with zero counts if present
    test_df = test_df.loc[[1, 0], [feature_column, 'background']]  # Enforce layout
    test_df.index = ['selection', 'other']

    # Fishers test and 95% CI e^(ln OR +/- 1.96 * sqrt(1/a+1/b+1/c+1/d))
    oddsratio, pvalue = stats.fisher_exact(test_df)
    ci = 1.96 * np.sqrt((1 / test_df).values.sum())
    lower_ci, upper_ci = np.exp(np.log(oddsratio) - ci), np.exp(np.log(oddsratio) + ci)

    return test_df, (oddsratio, pvalue, lower_ci, upper_ci)


def collect_column_structure_stats(aligned_prointvar_table):
    """
    Produce standard column-structure feature frequency table.

    :param aligned_prointvar_table:
    :return:
    """
    # Generate standard tables
    lig_stats = _column_ligand_contacts(aligned_prointvar_table)
    cov_stats = _column_total_contacts(aligned_prointvar_table)
    pips = _column_protein_contacts(aligned_prointvar_table)

    # Merge
    structure_stats = cov_stats.join(lig_stats).join(pips)

    # Interpolate index (useful for plotting)
    structure_stats = _interpolate_index(structure_stats)  # TODO: index to *full* alignment?

    return structure_stats


def column_set_ppi_enrichment(column_annotation_table, selection_query):
    """
    Test a subset of columns for enrichment of protein-protein interactions.

    :param column_annotation_table:
    :param selection_query:
    :return:
    """
    return _column_set_feature_enrichment(column_annotation_table, selection_query,
                                          feature_column='protein_protein_interactions',
                                          feature_background='sequences_with_contacts')


def column_set_ligand_enrichment(column_annotation_table, selection_query):
    """
    Test a subset of columns for enrichment of protein-ligand interactions.

    :param column_annotation_table:
    :param selection_query:
    :return:
    """
    return _column_set_feature_enrichment(column_annotation_table, selection_query,
                                          feature_column='protein_ligand_interactions',
                                          feature_background='sequences_with_contacts')
