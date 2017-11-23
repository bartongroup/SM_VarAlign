import numpy as np
import pandas as pd
from scipy import stats


def _column_ligand_contacts(aligned_prointvar_table):
    """
    Calculate frequency of ligand interactions over alignment columns. Gives both total counts and number of sequences
    with evidence of a ligand interaction.

    :param aligned_prointvar_table:
    :return:
    """
    # Select ligand interactions from structure table
    ligand_contacts = aligned_prointvar_table.query('interaction_type == "Protein-Ligand"')
    # Count interactions at each column
    ligand_evidence = \
        ligand_contacts.groupby(['SOURCE_ID_A', 'Alignment_column_A']).size().groupby('Alignment_column_A')
    protein_ligand_interactions = ligand_evidence.size()
    pdb_ligand_interactions = ligand_evidence.sum()
    # Format results and return DataFrame
    protein_ligand_interactions.name = 'protein_ligand_interactions'
    pdb_ligand_interactions.name = 'total_ligand_interactions'
    return pd.concat([protein_ligand_interactions, pdb_ligand_interactions], axis=1).fillna(0)


def _column_total_contacts(aligned_prointvar_table):
    """
    Calculate total number of interactions observed at each alignment column. Gives both total counts and number of
    sequences with at least one interaction.

    :param aligned_prointvar_table:
    :return:
    """
    # Count coverage at each column
    column_contacts = \
        aligned_prointvar_table.groupby(['SOURCE_ID_A', 'Alignment_column_A']).size().groupby('Alignment_column_A')
    column_sequence_contacts = column_contacts.size()
    column_total_contacts = column_contacts.sum()
    # Format results and return DataFrame
    column_sequence_contacts.name = 'sequences_with_contacts'
    column_total_contacts.name = 'total_contacts'
    return pd.concat([column_sequence_contacts, column_total_contacts], axis=1).fillna(0)


def _column_protein_contacts(aligned_prointvar_table):
    """
    Calculate frequency of protein-protein interactions over alignment columns. Gives both total counts and number of
    sequences with evidence of a protein-protein interaction.

    :param aligned_prointvar_table:
    :return:
    """
    # Select ligand interactions from structure table
    protein_contacts = aligned_prointvar_table.query('interaction_type == "Protein-Protein" & polymer_topology == "Interpolymer"')
    # Count interactions at each column
    protein_evidence = \
        protein_contacts.groupby(['SOURCE_ID_A', 'Alignment_column_A']).size().groupby('Alignment_column_A')
    protein_protein_interactions = protein_evidence.size()
    pdb_protein_interactions = protein_evidence.sum()
    # Format results and return DataFrame
    protein_protein_interactions.name = 'protein_protein_interactions'
    pdb_protein_interactions.name = 'total_protein_interactions'
    return pd.concat([protein_protein_interactions, pdb_protein_interactions], axis=1).fillna(0)


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
