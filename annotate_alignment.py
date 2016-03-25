import argparse
from Bio import AlignIO, SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.pairwise2 import format_alignment
import code
import pandas as pd
import re
import urllib2

# Use my developement branch of ProteoFAV
import sys
sys.path.extend(['/Users/smacgowan/PycharmProjects/ProteoFAV'])
from proteofav.variants import select_uniprot_variants
from proteofav.analysis.utils import expand_dataframe

import logging

log = logging.getLogger(__name__)

def get_row_residue_numbers(align, uniprot_seqs, use_local_alignment):
    """
    Map each sequence in an alignment to a longer sequence and return the residue numbers.

    :param align: MultipleSeqAlignment for processing.
    :param uniprot_seqs: List of full-length SeqRecords.
    :param use_local_alignment: If True, use a local alignment rather than exact pattern matching.
    :return: A list of tuples with the sequence IDs and residues numbers for each sequence in the alignment.
    """
    # Prepare sequences for alignment to UniProt by removing gaps
    unaligned_seqs = []
    for a in align:
        unaligned_seqs.append(SeqRecord(Seq(str(a.seq).replace('-', '').upper(), a.seq.alphabet), id=a.id))

    if use_local_alignment:
        # Align input alignment sequences to UniProt Sequences
        alignments = []
        for s, uniprot_seq in zip(unaligned_seqs, uniprot_seqs):
            log.debug('Aligning sub-sequence: {} to UniProt sequence: {}'.format(s.id, uniprot_seq[1].id))
            local_alignment = pairwise2.align.localxs(s.seq, uniprot_seq[1].seq, -.5, -.1)
            for pairwise_alignment in local_alignment:
                log.debug('{}'.format('\n' + format_alignment(*pairwise_alignment)))
            alignments.append((s.id, uniprot_seq[1].id, local_alignment))

        # Build list of UniProt residue numbers for each non-gap for each sequence
        align_res_nums = []
        for sub_seq_id, uniprot_seq_id, pairwise in alignments:
            seq = str(pairwise[0][0])
            res_nums = [i for i, s in enumerate(seq) if s != '-']  # TODO: wrong if other seq has gaps too
            align_res_nums.append((sub_seq_id, uniprot_seq_id, res_nums))
    else:
        align_res_nums = []
        for s, uniprot_seq in zip(unaligned_seqs, uniprot_seqs):
            str_seq = str(s.seq)
            uniprot_str_seq = str(uniprot_seq[1].seq)
            if str_seq in uniprot_str_seq:
                log.debug('Matching sub-sequence: {} to UniProt sequence: {}'.format(s.id, uniprot_seq[1].id))
                start = uniprot_str_seq.find(str_seq)
                end = start + len(str_seq)
                res_nums = range(start, end)
                align_res_nums.append((s.id, uniprot_seq[1].id, res_nums))
            else:
                log.error('Sub-sequence: {} does not match UniProt sequence: {}'.format(s.id, uniprot_seq[1].id))
                log.error('Sub-sequence: {} does not match UniProt sequence: {}'.format(str_seq, uniprot_str_seq))
                raise TypeError

    return align_res_nums


def get_sequence_column_numbers(align):
    """
    Build list of column numbers for each non-gap for each sequence.

    :param align:
    :return:
    """
    align_col_nums = []
    for a in align:
        seq_id = a.id
        seq = str(a.seq)
        col_nums = [i for i, s in enumerate(seq) if s != '-']
        align_col_nums.append((seq_id, col_nums))

    return align_col_nums


def fetch_uniprot_sequences(protein_identifiers):
    """
    Retrieve UniProt sequences.

    :param protein_identifiers: List of protein identifiers (UniProt IDs or protein name)
    :return: List of protein sequences.
    """
    url = 'http://www.uniprot.org/uniprot/'
    uniprot_sequences = []
    for p in protein_identifiers:
        p = p.strip()
        remote_fasta = url + p + '.fasta'
        handle = urllib2.urlopen(remote_fasta)
        for seq_record in SeqIO.parse(handle, "fasta"):
            p = seq_record.id.split('|')[1]  # Extract UniProt ID
            uniprot_sequences.append((p, seq_record))
    return uniprot_sequences


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


def _fetch_variants(prots):
    # Get variant data
    # Get the data with EnsEMBL variants
    tables = []
    for p in list(set(prots)):
        try:
            variant_table = select_uniprot_variants(p, reduced_annotations=False)
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
    # And dedup, bearing in mind the same variant can pop up in different transcripts
    # (so dedupping is only done on certain columns)
    concat_table = concat_table.drop_duplicates(['UniProt_dbAccessionId',
                                                 'start',
                                                 'end',
                                                 'variant_id',
                                                 'to_aa_expanded']).reset_index(drop=True)
    # is_somatic = concat_table['variant_id'].apply(lambda x: x.startswith('COS'))
    is_germline = concat_table['variant_id'].apply(lambda x: x.startswith('rs'))
    # somatic_table = concat_table[is_somatic]
    germline_table = concat_table[is_germline]
    return germline_table


def fill_variant_count(value_counts, length=160):
        variants_per_pos = []
        for i in xrange(length):
            col_pos = i + 1
            try:
                variants_per_pos.append((col_pos, value_counts[col_pos]))
            except KeyError:
                variants_per_pos.append((col_pos, 0))
        return variants_per_pos


def write_jalview_annotation(ordered_values, file_name, title, description):
    with open(file_name, 'w') as results_file:
        results_file.write('JALVIEW_ANNOTATION\n')
        if isinstance(ordered_values, tuple) and all(map(lambda x: isinstance(x, str), [title, description])):
            results_file.write('BAR_GRAPH\t{}\t{}\t'.format(title, description) +
                               '|'.join('{},,{}'.format(str(x), str(x)) for x in ordered_values))
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


def main(args):
    """
    Fetch variants for identified protein sequences in an MSA, map to residues and columns and write Jalview feature
    files with key statistics.

    :return:
    """
    # Read alignment
    alignment = AlignIO.read(args.fasta_file, "fasta")

    # Get UniProt sequences
    protein_identifiers = [re.search('\w*', sequence.id).group().strip() for sequence in alignment]
    uniprot_sequences = fetch_uniprot_sequences(protein_identifiers)
    protein_identifiers = zip(*uniprot_sequences)[0]  # Ensure prots contains UniProt IDs (could be protein names)

    # Map columns to residues
    alignment_column_numbers = get_sequence_column_numbers(alignment)
    alignment_residue_numbers = get_row_residue_numbers(alignment, uniprot_sequences, args.use_local_alignment)
    mapped = map_columns_to_residues(alignment_column_numbers, alignment_residue_numbers)

    # Fetch variants
    germline_table = _fetch_variants(protein_identifiers)

    # Merge the data
    # Merge variant table and key table
    merged_table = pd.merge(mapped, germline_table,
                            left_on=['UniProt_ID', 'uniprot_res_num'],
                            right_on=['UniProt_dbAccessionId', 'start'])

    # Counting variants and writing Jalview annotations
    total_variant_counts = merged_table['alignment_col_num'].value_counts(sort=False)
    total_variants_per_column = fill_variant_count(total_variant_counts)
    write_jalview_annotation(zip(*total_variants_per_column)[1], 'variants_per_column.csv', 'Total_Variants',
                             'Total number of variants in summed over all proteins.')

    # Make some other counts
    is_missense = (merged_table['type'] == 'missense_variant') & \
                  (merged_table['from_aa'] != merged_table['to_aa_expanded'])
    is_ED = (merged_table['from_aa'] == 'E') & (merged_table['to_aa_expanded'] == 'D')
    is_DE = (merged_table['from_aa'] == 'D') & (merged_table['to_aa_expanded'] == 'E')

    missense_variant_counts = merged_table.loc[is_missense, 'alignment_col_num'].value_counts(sort=False)
    missense_variants_per_column = fill_variant_count(missense_variant_counts)

    missense_exc_DE_counts = merged_table.loc[is_missense & ~(is_ED | is_DE), 'alignment_col_num'].value_counts(sort=False)
    missense_exc_DE_per_column = fill_variant_count(missense_exc_DE_counts)

    variant_counts = [zip(*missense_variants_per_column)[1],
                      zip(*missense_exc_DE_per_column)[1]]
    titles = ['Missense_Variants',
              'Missense_Variants (exc. DE)']
    descriptions = ['Total number of missense variants in summed over all proteins.',
                    'Number of missense variants excluding E-D and D-E summed over all proteins.']

    write_jalview_annotation(variant_counts, 'jalview_annotations.csv', titles, descriptions)

    # TODO: Need to check key exists as won't if no pathogenic...
    #raw_pathogenic = pd.crosstab(merged_table['alignment_col_num'], merged_table['clinical_significance'])['pathogenic']
    #patho_per_pos = fill_variant_count(raw_pathogenic)
    #write_jalview_annotation(zip(*patho_per_pos)[1], 'patho_per_column.csv', 'Pathogenic_variants',
    #                         'Number of variants annotated pathogenic by ClinVar.')

    return merged_table


if __name__ == '__main__':
    # CLI Arguments
    parser = argparse.ArgumentParser(description='Create variant counts column wise for an alignment.')
    parser.add_argument('fasta_file', type=str, help='Path to fasta file containing MSA.')
    parser.add_argument('--use_local_alignment', action='store_true',
                        help='Align sub-sequences to UniProt rather than enforcing exact match.')
    parser.add_argument('--interpreter', action='store_true',
                        help='Drop into interactive python session once analysis is complete.')
    args = parser.parse_args()

    merged_table = main(args)

    if args.interpreter:
        code.interact(local=dict(globals(), **locals()))