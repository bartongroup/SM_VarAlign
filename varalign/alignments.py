from uniprot import _strip_version


def get_accession(sequence, strip_version=False):
    """
    Return sequence accession from annotations.

    :param sequence:
    :param strip_version:
    :return:
    """
    accession = sequence.annotations['accession']
    if strip_version:
        accession = _strip_version(accession)
    return accession


def get_start_end(sequence):
    """
    Return sequence start and end from annotations

    :param sequence:
    :return:
    """
    return sequence.annotations['start'], sequence.annotations['end']


def index_seq_to_alignment(sequence, gap_chars=('-', '.'),
                           zero_index=False, reverse_mapping = False):
    """
    Index the residue numbers of a sequence to its alignment.

    :param sequence:
    :param gap_chars:
    :param zero_index:
    :param reverse_mapping:
    :return:
    """
    adj_i = 1
    if zero_index:
        adj_i = 0
    alignment_index = [i+adj_i for i, x in enumerate(sequence.seq) if x not in gap_chars]
    start, end = get_start_end(sequence)
    sequence_index = range(start, end+1)
    assert len(alignment_index) == len(sequence_index)
    if reverse_mapping:
        return zip(sequence_index, alignment_index)
    else:
        return zip(alignment_index, sequence_index)
