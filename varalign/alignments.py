import pandas as pd

from varalign.uniprot import _strip_version


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


def index_seq_to_alignment(sequence, gap_chars=set(['-', '.']),
                           zero_index=False, reverse_mapping=False):
    """
    Index the residue numbers of a sequence to its alignment.

    :param sequence:
    :param gap_chars:
    :param zero_index:
    :param reverse_mapping:
    :return:
    """
    index_start = 1
    if zero_index:
        index_start = 0
    seq_characters = str(sequence.seq)
    alignment_index = [i for i, x in enumerate(seq_characters, start=index_start) if x not in gap_chars]
    start, end = get_start_end(sequence)
    sequence_index = list(range(start, end+1))
    assert len(alignment_index) == len(sequence_index)
    if reverse_mapping:
        return list(zip(sequence_index, alignment_index))
    else:
        return list(zip(alignment_index, sequence_index))


def alignment_info_table(alignment, id_filter=''):
    """
    Build a table with key alignment info.
    
    Scan an alignment and return a Pandas DataFrame with sequence identifiers (seq.id, names and UniProts),
    alignment index to sequence mappings and source species.
    
    :param alignment:
    :param id_filter:
    :return: 
    """
    alignment_info = []
    for sequence in alignment:
        if id_filter in sequence.id:
            alignment_info.append((sequence.id,
                                   sequence.name,
                                   get_accession(sequence),
                                   get_start_end(sequence),
                                   index_seq_to_alignment(sequence)))
    alignment_info = pd.DataFrame(alignment_info, columns=['seq_id', 'name', 'uniprot', 'start_end', 'mapping'])
    alignment_info = alignment_info.assign(species=alignment_info['name'].str.split('_').str[1].values)
    alignment_info = alignment_info.assign(uniprot_id=alignment_info['uniprot'].str.split('.').str[0].values)

    # Add sequence lengths
    seq_lengths = alignment_info['start_end'].apply(lambda x: len(list(range(*x))) + 1)
    seq_lengths.name = 'length'
    alignment_info = alignment_info.join(seq_lengths)

    return alignment_info

