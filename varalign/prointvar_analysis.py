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


def _download_structure_data(alignment_info_table, logfile='prointvar_download'):
    download_status = []
    for uniprot_id in tqdm.tqdm(alignment_info_table['uniprot_id'].unique()):
        with open(logfile + '.out.' + uniprot_id, 'w') as process_out:
            with open(logfile + '.log.' + uniprot_id, 'w') as process_err:
                exit_code = subprocess.call(
                    ['ProIntVar', 'download', '--mmcif', '--pdb', '--sifts', '--best_structures', uniprot_id],
                    stdout=process_out, stderr=process_err)
                sifts_best_order = pd.read_table(process_out.name, header=None,
                                                 names=['sifts_index', 'pdb_id', 'chain_id'])
                sifts_best_order['query'] = uniprot_id
                download_status.append([uniprot_id, exit_code, process_out.name, process_err.name, sifts_best_order])
    # Format sifts best order
    sifts_result_table = pd.concat([x.pop() for x in download_status])
    return download_status, sifts_result_table


def _format_mapping_table(alignment_info_table, alignment_mapping_table):
    # Format mapping table for join
    aln_uniprot_ids = alignment_info_table.set_index('seq_id').loc[:, 'uniprot_id']
    aln_uniprot_ids.index.name = 'SOURCE_ID'
    mapping_table = alignment_mapping_table.join(aln_uniprot_ids)
    mapping_table.reset_index(inplace=True)
    mapping_table.loc[:, 'Protein_position'] = mapping_table.loc[:, 'Protein_position'].astype(str)
    mapping_table.set_index(['uniprot_id', 'Protein_position'], inplace=True)
    mapping_table.index.names = ['UniProt_dbAccessionId_A', 'UniProt_dbResNum_A']
    mapping_table.columns = ['SOURCE_ID', 'Alignment_column']
    return mapping_table


def _merge_alignment_columns_to_contacts(alignment_mappings, contacts_table):
    # Format structure table for join
    contacts_table = contacts_table.dropna(subset=['UniProt_dbAccessionId_A', 'UniProt_dbResNum_A'])
    contacts_table = contacts_table.set_index(['UniProt_dbAccessionId_A', 'UniProt_dbResNum_A'], drop=False)
    # Map ATOM_A contacts to alignment
    contacts_table = alignment_mappings.join(contacts_table, how='inner')  # Left join indicates unmapped residues...
    log.info('{} atom-atom records after adding ATOM_A alignment columns.'.format(len(contacts_table)))
    # Map ATOM_B contacts to alignment
    alignment_mappings.index.rename(['UniProt_dbAccessionId_B', 'UniProt_dbResNum_B'], inplace=True)
    contacts_table = pd.merge(alignment_mappings, contacts_table,
                              right_on=['UniProt_dbAccessionId_B', 'UniProt_dbResNum_B'],
                              left_index=True, how='right', suffixes=['_B', '_A'])  # Note right join here
    log.info('{} atom-atom records after adding ATOM_B alignment columns.'.format(len(contacts_table)))
    return contacts_table


def _sort_ab_contacts(aligned_contacts_table):
    """
    Sort the entries in an aligned contacts table so that alignment column A <= B.

    :param aligned_contacts_table:
    :return: AB sorted contacts table.
    """
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

    rename_dict = {c: _swap_ending(c) for c in aligned_contacts_table.columns}
    aligned_contacts_table = pd.concat(
        [aligned_contacts_table.query('Alignment_column_A <= Alignment_column_B'),
         aligned_contacts_table.query('Alignment_column_A > Alignment_column_B').rename(columns=rename_dict),
         aligned_contacts_table.query('Alignment_column_B != Alignment_column_B')  # NaN fields, value != value
         ]
    )
    return aligned_contacts_table


def _dedupe_ab_contacts(contacts_table):
    """
    Remove duplicate contacts from a contacts_table or derivative.

    :param contacts_table:
    :return:
    """
    atom_id_fields = ['PDB_dbAccessionId', 'PDB_dbChainId', 'PDB_dbResNum', 'PDB_entityId']
    atom_id_fields = ['{}_{}'.format(field, ab) for ab in ['A', 'B'] for field in atom_id_fields]
    is_duplicated = contacts_table[atom_id_fields].apply(frozenset, axis=1, raw=True).duplicated()
    log.info('Removing {} duplicates...'.format(sum(is_duplicated)))
    contacts_table = contacts_table.loc[~is_duplicated]
    log.info('{} atom-atom records after removing A/B, B/A duplicates.'.format(len(contacts_table)))
    return contacts_table


def _filter_extra_domain_contacts(prointvar_table, alignment_info):
    """
    Filter contacts that do not involve at least one residue in the alignment.

    :param prointvar_table:
    :param alignment_info:
    :return:
    """
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

    prointvar_table = prointvar_table.groupby('UniProt_dbAccessionId_A').apply(_filter_by_sequence_range,
                                                                               alignment_info=alignment_info)
    log.info('{} atom-atom records remain with >=1 residue in alignment sequence.'.format(len(prointvar_table)))
    return prointvar_table


if __name__ == '__main__':
    # CLI
    parser = argparse.ArgumentParser(description='Structural properties of alignment columns.')
    parser.add_argument('alignment', type=str, help='Path to the alignment.')
    parser.add_argument('--n_proc', type=int, help='Number of processors.', default=1)
    args = parser.parse_args()

    # Read data produced by `align_variants.py`
    aln_info = pd.read_pickle(args.alignment+'_info.p.gz')
    indexed_mapping_table = pd.read_pickle(args.alignment+'_mappings.p.gz')
    column_stats = pd.read_csv(args.alignment + '.col_summary.csv')

    # Get SIFTS best and download for all proteins in alignment
    download_logfile = args.alignment + '_prointvar_download.log'
    status, downloaded = _download_structure_data(aln_info, download_logfile)

    # Process all downloaded structural data with ProIntVar
    to_load = downloaded['pdb_id'].dropna().unique()
    p = multiprocessing.Pool(args.n_proc)
    tabs = list(tqdm.tqdm(p.imap(_format_structure_data, to_load), total=len(to_load)))
    structure_table = pd.concat(tabs)
    log.info('{} atom-atom records created.'.format(len(structure_table)))  # Shouldn't this be even?

    # Filter non-alignment residues from the table
    structure_table = _filter_extra_domain_contacts(structure_table, aln_info)

    # Remove each contact's duplicate row
    structure_table = _dedupe_ab_contacts(structure_table)

    # Add alignment columns to table
    aln_mappings = _format_mapping_table(aln_info, indexed_mapping_table)
    structure_table = _merge_alignment_columns_to_contacts(aln_mappings, structure_table)

    # Sort A/B contacts
    structure_table = _sort_ab_contacts(structure_table)

    # Write structure table
    log.info('Writing {} atom-atom records to file...'.format(len(structure_table)))
    structure_table.to_pickle(args.alignment+'_prointvar_structure_table.p.gz')
    log.info('DONE.')
