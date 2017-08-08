import argparse
from Bio import AlignIO
import gzip
import io
import StringIO

def index_pfam(pfam_path):
    """
    Write an index containing the start positions of the alignments in a Pfam file.

    :param pfam_path: Path to Pfam file.
    :return: None
    """
    # Open file
    if pfam_path.endswith('.gz'):
        pfam_file = gzip.open(pfam_path, 'r')
    else:
        pfam_file = open(pfam_path, 'r')

    # Create the index
    position = 0
    alignment_starts = []
    family_names = []
    while True:
        line = pfam_file.readline()
        if not line:
            break
        if line == '# STOCKHOLM 1.0\n':
            alignment_starts.append(position)
        if line.startswith('#=GF AC'):
            family_names.append(line.split()[-1].strip())
        position = pfam_file.tell()

    pfam_file.close()

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
    # Open file
    if pfam_path.endswith('.gz'):
        pfam = gzip.open(pfam_path, 'r')
    else:
        pfam = open(pfam_path, 'r')

    pfam.seek(start)
    alignment_handle = StringIO.StringIO()
    for line in pfam:
        alignment_handle.write(line)
        if line == '//\n':
            break
    pfam.close()

    alignment_handle.seek(0)
    alignment = AlignIO.read(alignment_handle, 'stockholm')
    alignment_handle.close()
    return alignment


if __name__ == '__main__':
    # CLI
    parser = argparse.ArgumentParser(description='Index Pfam alignments.')
    parser.add_argument('pfam_file', type=str, help='Path to the Pfam file.')
    args = parser.parse_args()

    # Index Pfam file
    index_pfam(args.pfam_file)
