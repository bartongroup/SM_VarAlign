def get_accession(sequence):
    """
    Return sequence accession from annotations.

    :param sequence:
    :return:
    """
    return sequence.annotations['accession']


def get_start_end(sequence):
    """
    Return sequence start and end from annotations

    :param sequence:
    :return:
    """
    return sequence.annotations['start'], sequence.annotations['end']


def index_seq_to_alignment(sequence, gap_chars=('-', '.'),
                           zero_index=False):
    """
    Index the residue numbers of a sequence to its alignment.

    :param sequence:
    :param gap_chars:
    :param zero_index:
    :return:
    """
    adj_i = 1
    if zero_index:
        adj_i = 0
    alignment_index = [i+adj_i for i, x in enumerate(sequence.seq) if x not in gap_chars]
    start, end = get_start_end(sequence)
    sequence_index = range(start, end+1)
    assert len(alignment_index) == len(sequence_index)
    return zip(alignment_index, sequence_index)