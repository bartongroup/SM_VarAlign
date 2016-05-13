import pandas as pd
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.pairwise2 import format_alignment

import logging
from utils import parse_seq_name

log = logging.getLogger(__name__)


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
        seq_indexes = range(1, len(res_nums) + 1)
        align_res_nums = sub_seq_id, uniprot_seq_id, res_nums, seq_indexes
    else:
        str_seq = str(subseq.seq)
        uniprot_str_seq = str(uniprot_seq[1].seq)
        if str_seq in uniprot_str_seq:
            log.debug('Matching sub-sequence: {} to UniProt sequence: {}'.format(subseq.id, uniprot_seq[1].id))
            start = uniprot_str_seq.find(str_seq) + 1
            end = start + len(str_seq)
            res_nums = range(start, end)
            seq_indexes = range(1, len(res_nums) + 1)
            align_res_nums = subseq.id, uniprot_seq[1].id, res_nums, seq_indexes
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


def map_columns_to_res_nums(col_nums, aligned_res_nums):
    """
    Map alignment columns to UniProt residue numbers.

    :param col_nums:
    :param res_nums:
    :param seq_id:
    :param seq_index:
    :param uniprot_seq_id:
    :return:
    """
    seq_id, uniprot_seq_id, res_nums, seq_index = aligned_res_nums
    log.debug('Mapping {}...'.format(seq_id))
    record = {'seq_id': seq_id, 'uniprot_seq_id': uniprot_seq_id, 'uniprot_res_num': res_nums,
              'alignment_col_num': col_nums, 'sequence_index': seq_index}
    # TODO: This is messy due to dependency on uniprot_seq_id format in `alignment_residue_numbers`
    if '|' in record['uniprot_seq_id']:
        prot_name = record['uniprot_seq_id'].split('|')[1]  # UniProt ID
    else:
        prot_name = record['uniprot_seq_id']
    record.update({'UniProt_ID': prot_name})

    return pd.DataFrame(record)