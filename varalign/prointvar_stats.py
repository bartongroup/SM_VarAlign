import pandas as pd


def _column_ligand_contacts(aligned_prointvar_table):
    """
    Calculate frequency of ligand interactions over alignment columns. Gives both total counts and number of sequences
    with evidence of a ligand interaction.

    :param aligned_prointvar_table:
    :return:
    """
    # Select ligand interactions from structure table
    ligand_contacts = aligned_prointvar_table.query('group_PDB_B == "HETATM" & ATOM_B != "HOH"')
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


def _interpolate_index(table):
    """
    Re-index an interger indexed table to include missing values within the range of the original index.

    :param table:
    :return:
    """
    idx = pd.Index(range(int(table.index.min()), int(table.index.max())))
    return table.reindex(idx)
