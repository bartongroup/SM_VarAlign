import argparse
from Bio import AlignIO, SeqIO
import logging

log = logging.getLogger(__name__)


def main(args):
    """

    :param args:
    :return:
    """
    swissprot = SeqIO.index(args.local_uniprot, 'swiss')
    pfam = AlignIO.parse(args.local_pfam, 'stockholm')

    for a in pfam:
        pass


if __name__ == '__main__':
    # CLI Arguments
    parser = argparse.ArgumentParser(description='Create variant counts column wise for an alignment.')
    parser.add_argument('--local_uniprot', type=str, help='Path to uniprot_sprot.dat')
    parser.add_argument('--local_pfam', type=str, help='Path to PFAM-A.seed')
    args = parser.parse_args()

    main(args)