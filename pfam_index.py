from Bio import AlignIO
import os
import StringIO

def index_pfam(pfam_path):
    """
    Write an index containing the start positions of the alignments in a Pfam file.

    :param pfam_path: Path to Pfam file.
    :return: None
    """
    # Create the index
    with open(pfam_path, 'r') as pfam_file:
        position = 0
        file_size = os.stat(pfam_path).st_size
        alignment_starts = []
        family_names = []
        while position < file_size:
            line = pfam_file.readline()
            if line == '# STOCKHOLM 1.0\n':
                alignment_starts.append(position)
            if line.startswith('#=GF AC'):
                family_names.append(line.split()[-1].strip())
            position = pfam_file.tell()

    # Write the index
    with open(pfam_path + '.idx', 'w') as index_file:
        for e in zip(alignment_starts, family_names):
            index_file.write(','.join(map(str, e)) + '\n')


def lookup_index(index_path, family=None):
    """
    Lookup the byte offset start of a Pfam family.

    :param index_path: Path to index file.
    :param family: Family accession code.
    :return: byte offset
    """
    with open(index_path, 'r') as index_file:
        if not (isinstance(family, str) and len(family) >= 7):
            msg = 'Invalid family: {}. (Should resemble PF12345[.12])'.format(family)
            raise ValueError(msg)
        for line in index_file:
            offset, ac = line.strip().split(',')
            if ac.startswith(family):
                return int(offset)
        print '{} not found.'.format(family)
    return None


def read_family(pfam_path, start):
    """
    Read an alignment from the Pfam file.

    :param pfam_path: Path to Pfam file.
    :param start: Byte offset of alignment start.
    :return: Biopython multiple alignment.
    """
    with open(pfam_path, 'r') as pfam:
        pfam.seek(start)
        alignment_handle = StringIO.StringIO()
        for line in pfam:
            alignment_handle.write(line)
            if line == '//\n':
                break
    alignment_handle.seek(0)
    alignment = AlignIO.read(alignment_handle, 'stockholm')
    alignment_handle.close()
    return alignment