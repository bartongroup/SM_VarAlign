import argparse
from Bio import AlignIO, SeqIO
import logging
import sys

sys.path.insert(0, '/homes/smacgowan/lib/python/TBB_ProteoFAV')
import proteofav


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create variant counts column wise for an alignment.')
    parser.add_argument('alignment', type=str, help='Path to the alignment.')
    parser.add_argument('--format', type=str, default='stockholm',
                        help='Alignment format.')
    parser.add_argument('--filter', type=str,
                        help='An inclusive sequence ID filter to process only a subset of sequences.')
    args = parser.parse_args()
