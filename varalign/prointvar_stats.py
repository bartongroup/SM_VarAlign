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


def _pdb_protein_coverage(aligned_prointvar_table):
    # Count coverage at each column
    structure_cov = \
        aligned_prointvar_table.groupby(['SOURCE_ID_A', 'Alignment_column_A']).size().groupby('Alignment_column_A')
    protein_structure_cov = structure_cov.size()
    pdb_structure_cov = structure_cov.sum()
    # Format results and return DataFrame
    protein_structure_cov.name = 'protein_structure_coverage'
    pdb_structure_cov.name = 'total_structure_coverage'
    return pd.concat([protein_structure_cov, pdb_structure_cov], axis=1).fillna(0)
