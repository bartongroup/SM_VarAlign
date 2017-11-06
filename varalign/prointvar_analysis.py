"""
This is the script that runs the structural analyses using `ProIntVar`.

It requires Python 3, for `ProIntVar` compatibility, and should be run after `align_variants.py`.
"""

import argparse
import logging
import matplotlib.pyplot as plt
import multiprocessing
import os
import pandas as pd
from prointvar import merger
import subprocess
import tqdm


log = logging.getLogger(__name__)
log.setLevel('INFO')


def _format_structure_data(pdb):
    # Read structure data
    try:
        pdbx, dssp, sifts, contacts = merger.table_generator(pdb_id=pdb, bio=False, contacts=True,
                                                             override=False, residue_agg=True, dssp=False)
    except FileNotFoundError as e:
        print('{} failed with FileNotFoundError:'.format(pdb))
        return None
    # Merge tables
    table = merger.TableMerger(pdbx_table=pdbx, sifts_table=sifts,
                               contacts_table=contacts, dssp_table=dssp).merge()
    return table


def _filter_by_sequence_range(g, alignment_info):
    uniprot_id = g.name
    start_ends = alignment_info.query('uniprot_id == @uniprot_id')['start_end']  # Could be > 1 range
    if start_ends.empty:
        return None
    else:
        in_seq_range = []
        for test_range in start_ends:
            in_seq_range.append(pd.to_numeric(g['UniProt_dbResNum_A']).between(*test_range))
        result = pd.concat(in_seq_range, axis=1).any(axis=1)
        return g[result]


def _swap_ending(x, swap=['_A', '_B']):
    a, b = swap
    if x.endswith(a):
        # return x.replace(*swap)
        return b.join(x.rsplit(a, 1))  # Prevent 'PDB_Annotation_A': 'PDB_Bnnotation_B'
    elif x.endswith(b):
        # return x.replace(*swap[::-1])
        return a.join(x.rsplit(b, 1))
    else:
        return x


if __name__ == '__main__':
    # CLI
    parser = argparse.ArgumentParser(description='Structural properties of alignment columns.')
    parser.add_argument('alignment', type=str, help='Path to the alignment.')
    parser.add_argument('--n_proc', type=int, help='Number of processors.', default=1)
    args = parser.parse_args()

    # Read data produced by `align_variants.py`
    aln_info = pd.read_pickle(args.alignment+'_info.p.gz')
    #alignment_variant_table = pd.read_pickle(args.alignment+'_variants.p.gz')
    indexed_mapping_table = pd.read_pickle(args.alignment+'_mappings.p.gz')
    column_stats = pd.read_csv(args.alignment + '.col_summary.csv')

    # Get SIFTS best and download for all proteins in alignment
    status = []
    with open('prointvar_download.log', 'wb') as splog:
        for pid in tqdm.tqdm(aln_info['uniprot_id'].unique()):
            status.append((pid, subprocess.call(['ProIntVar', 'download',
                                                 '--mmcif', '--pdb', '--sifts',
                                                 '--best_structures', pid],
                                                stdout=splog, stderr=subprocess.STDOUT)))

    # Process all downloaded structural data with ProIntVar
    to_load = [x.split('.')[0] for x in os.listdir('sifts/')]  # NB. ALL pdbs in db_sifts
    # tabs = []
    # errors = []
    # for pdb in tqdm.tqdm_notebook(to_load[:20]):
    #     tabs.append(_format_structure_data(pdb))
    p = multiprocessing.Pool(args.n_proc)
    tabs = list(tqdm.tqdm(p.imap(_format_structure_data, to_load), total=len(to_load)))
    structure_table = pd.concat(tabs)

    # Pre-filter and format the structure table
    log.info('{} atom-atom records created.'.format(len(structure_table)))  # Shouldn't this be even?

    # Filter non-alignment residues from the table
    structure_table = structure_table.groupby('UniProt_dbAccessionId_A').apply(_filter_by_sequence_range,
                                                                               alignment_info=aln_info)
    log.info('{} atom-atom records remain with >=1 residue in alignment sequence.'.format(len(structure_table)))

    # Remove each contact's duplicate row
    atom_id_fields = ['PDB_dbAccessionId', 'PDB_dbChainId', 'PDB_dbResNum', 'PDB_entityId']
    atom_id_fields = ['{}_{}'.format(field, ab) for ab in ['A', 'B'] for field in atom_id_fields]
    is_duplicated = structure_table[atom_id_fields].apply(frozenset, axis=1, raw=True).duplicated()
    log.info('Removing {} duplicates...'.format(sum(is_duplicated)))
    structure_table = structure_table.loc[~is_duplicated]
    log.info('{} atom-atom records after removing A/B, B/A duplicates.'.format(len(structure_table)))

    # Add alignment columns to table
    # Format mapping table for join
    aln_uniprot_ids = aln_info.set_index('seq_id').loc[:, 'uniprot_id']
    aln_uniprot_ids.index.name = 'SOURCE_ID'
    aln_mappings = indexed_mapping_table.join(aln_uniprot_ids)
    aln_mappings.reset_index(inplace=True)
    aln_mappings.loc[:, 'Protein_position'] = aln_mappings.loc[:, 'Protein_position'].astype(str)
    aln_mappings.set_index(['uniprot_id', 'Protein_position'], inplace=True)
    aln_mappings.index.names = ['UniProt_dbAccessionId_A', 'UniProt_dbResNum_A']
    aln_mappings.columns = ['SOURCE_ID', 'Alignment_column']
    # Format structure table for join
    structure_table = structure_table.dropna(subset=['UniProt_dbAccessionId_A', 'UniProt_dbResNum_A'])
    structure_table = structure_table.set_index(['UniProt_dbAccessionId_A', 'UniProt_dbResNum_A'], drop=False)
    # Map ATOM_A contacts to alignment
    structure_table = aln_mappings.join(structure_table, how='inner')  # Left join indicates unmapped residues...
    log.info('{} atom-atom records after adding ATOM_A alignment columns.'.format(len(structure_table)))
    # Map ATOM_B contacts to alignment
    aln_mappings.index.rename(['UniProt_dbAccessionId_B', 'UniProt_dbResNum_B'], inplace=True)
    structure_table = pd.merge(aln_mappings, structure_table,
                               right_on=['UniProt_dbAccessionId_B', 'UniProt_dbResNum_B'],
                               left_index=True, how='right', suffixes=['_B', '_A'])  # Note right join here
    log.info('{} atom-atom records after adding ATOM_B alignment columns.'.format(len(structure_table)))

    # Sort A/B contacts
    rename_dict = {c: _swap_ending(c) for c in structure_table.columns}
    structure_table = pd.concat(
        [structure_table.query('Alignment_column_A <= Alignment_column_B'),
         structure_table.query('Alignment_column_A > Alignment_column_B').rename(columns=rename_dict),
         structure_table.query('Alignment_column_B != Alignment_column_B')
         ]
    )

    # Write structure table
    structure_table.to_pickle(args.alignment+'_prointvar_structure_table.p.gz')
