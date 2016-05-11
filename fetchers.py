"""
This module contains functions for retrieving external data required for the analysis. In particular,
the full-length UniProt sequences for the aligned sequences and the variant tables.
"""
import os.path
import urllib2
import pandas as pd
from Bio import SeqIO

from utils import urlopen_with_retry, query_uniprot, parse_seq_name

# Use my developement branch of ProteoFAV
import sys
sys.path.extend(['/Users/smacgowan/PycharmProjects/ProteoFAV'])
from proteofav.variants import select_uniprot_variants

import logging

log = logging.getLogger(__name__)


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
        try:
            seq_record = SeqIO.read(handle, "fasta")
        except ValueError:
            log.error('Could not retrieve sequence for {}'.format(seq_name))
            return None
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


def _fetch_variants(prots, downloads=None, save_name=None):
    """
    Fetch variants for a set of proteins from EnsEMBL.

    This function is built specifically to produce a multi-protein variant table that is compatible with the
    downstream analysis. It handles:
    - Concatenation of variant tables for each protein.
    - Expanding the 'to_aa' column
    - Parsing 'clinical_significance' column
    - De-duping
    - Separating germline from somatic
    (- Caching the resultant table)

    :param prots: List of UniProt ACs.
    :param downloads: Folder to save whole variant set to.
    :param save_name: Filename for whole variant set.
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
        log.debug('---Concatenating variant tables---')
        concat_table = pd.concat(tables, ignore_index=True)
        # Need to expand on 'to_aa' before dedupping
        log.debug('---Expanding `to_aa` column---')
        to_aa_columns = pd.DataFrame(concat_table.to_aa.tolist(), )
        split = pd.concat([concat_table, to_aa_columns], axis=1)
        concat_table = pd.melt(split, id_vars=list(concat_table.columns), value_name='to_aa_expanded')
        concat_table = concat_table[concat_table.to_aa_expanded.notnull()]
        concat_table = concat_table.drop('variable', 1)  # Remove melt variable
        # Fix or remove list columns
        concat_table = concat_table.drop('to_aa', 1)

        log.debug('---Parsing `clinical_significance`---')
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