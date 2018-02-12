#!/usr/bin/env python
import argparse

from varalign import pfam

if __name__ == '__main__':
    # CLI
    parser = argparse.ArgumentParser(description='Index Pfam alignments.')
    parser.add_argument('pfam_file', type=str, help='Path to the Pfam file.')
    args = parser.parse_args()

    # Index Pfam file
    pfam.index_pfam(args.pfam_file)
