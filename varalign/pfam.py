import argparse
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import gzip
import io
from itertools import count, groupby
import pandas as pd
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


def filter_non_swissprot(aln, swissprot_id_file='/homes/smacgowan/NOBACK/resources/swissprot_ids.tab'):
    """Return new Pfam alignment with only SwissProt sequences.

    :param aln:
    :param swissprot_id_file:
    :return:
    """
    # Read list of Swissprot IDS
    swissprot_ids = pd.read_csv(swissprot_id_file, sep='\t', header=0, comment='#')
    swissprot_seq_names = set(swissprot_ids['Entry name'].tolist())

    # Filter alignment
    # Test each sequence for SwissProt membership
    passing_seqs = [seq for seq in aln if seq.name in swissprot_seq_names]

    # Build new alignment for swissprot sequences
    filtered_alignment = MultipleSeqAlignment(passing_seqs)
    filtered_alignment.annotations = aln.annotations

    # Identify occupied columns
    occupied_cols = []
    for i in xrange(aln.get_alignment_length()):
        if not all([x == '-' for x in filtered_alignment[:, i]]):
            occupied_cols.append(i)

    degapped_alignment = _mask_alignment(filtered_alignment, occupied_cols)

    return degapped_alignment


def _mask_alignment(aln, column_indexes):
    """Return specified columns of an alignment.
    
    :param aln: 
    :param column_indexes: 
    :return: 
    """
    # Convert to ranges
    # Inspired by:
    # https://stackoverflow.com/questions/3429510/pythonic-way-to-convert-a-list-of-integers-into-a-string-of-comma-separated-range/3430231#3430231
    G = (list(x) for _, x in groupby(column_indexes, lambda x, c=count(): next(c) - x))
    mask_ranges = [(g[0], g[-1])[:len(g)] for g in G]
    mask_ranges = [x if len(x) == 2 else x * 2 for x in mask_ranges]
    # Build list of continuous occupied sub-alignments
    MSA_columns = []
    for start, end in mask_ranges:
        MSA_columns.append(aln[:, start:end + 1])

    # Concatenate sub-alignments
    masked_alignment = aln[:, 0:0]
    for x in MSA_columns:
        masked_alignment = masked_alignment + x
    # NB. `reduce` could be faster...
    # Add sequence annotations to new alignment
    for new, old in zip(masked_alignment, aln):
        new.annotations = old.annotations

    return masked_alignment


if __name__ == '__main__':
    # CLI
    parser = argparse.ArgumentParser(description='Index Pfam alignments.')
    parser.add_argument('pfam_file', type=str, help='Path to the Pfam file.')
    args = parser.parse_args()

    # Index Pfam file
    index_pfam(args.pfam_file)
