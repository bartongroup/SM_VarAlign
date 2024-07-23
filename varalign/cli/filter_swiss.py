import argparse
from Bio import AlignIO
import os
import varalign
from varalign.core import pfam

#!/usr/bin/env python


sp_whitelist = os.path.join(os.path.dirname(varalign.__file__), 'data', 'swissprot_uniprot-all.tab.gz')


def main():
    # CLI
    parser = argparse.ArgumentParser(description='Filter swissprot.')
    parser.add_argument('alignment', type=str, help='Alignment file.')
    parser.add_argument('--whitelist', type=str, help='Path to sequence whitelist (e.g. SwissProt', default=sp_whitelist)
    args = parser.parse_args()

    alignment = AlignIO.read(args.alignment, format='stockholm')
    new_alignment = pfam.filter_non_swissprot(alignment, swissprot_id_file=args.whitelist)
    path, filename = os.path.split(args.alignment)
    output = os.path.join(path, 'swissprot_'+filename)
    AlignIO.write(new_alignment, output, 'stockholm')


if __name__ == '__main__':
    main()
