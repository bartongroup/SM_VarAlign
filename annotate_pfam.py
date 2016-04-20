import annotate_alignment
import argparse
from Bio import AlignIO, SeqIO  # Needs PR #768 #769 patched
import logging

log = logging.getLogger(__name__)


def main(local_uniprot, local_pfam, seq_id_filter, use_local_alignment, downloads):
    """

    :param args:
    :return:
    """
    swissprot = SeqIO.index(local_uniprot, 'swiss')
    pfam = AlignIO.parse(local_pfam, 'stockholm')

    for alignment in pfam:
        alignment_name = alignment.annotations['GF']['AC'][0]
        if 'PF0000' in alignment_name:
            annotate_alignment.main(alignment, alignment_name, seq_id_filter, use_local_alignment, swissprot,
                                    False, downloads)


if __name__ == '__main__':
    # CLI Arguments
    parser = argparse.ArgumentParser(description='Create variant counts column wise for an alignment.')
    parser.add_argument('--local_uniprot_file', type=str, help='Path to uniprot_sprot.dat')
    parser.add_argument('--local_pfam', type=str, help='Path to PFAM-A.seed')
    parser.add_argument('--seq_id_filter', type=str,
                        help='An inclusive filter to process only a subset of sequences.')
    parser.add_argument('--use_local_alignment', action='store_true',
                        help='Align sub-sequences to UniProt rather than enforcing exact match.')
    parser.add_argument('--downloads', type=str, default='.VarAlign',
                        help='A directory to store downloaded files.')
    args = parser.parse_args()

    local_uniprot = args.local_uniprot_file
    local_pfam = args.local_pfam
    seq_id_filter = args.seq_id_filter
    use_local_alignment = args.use_local_alignment
    downloads = args.downloads

    main(local_uniprot, local_pfam, seq_id_filter, use_local_alignment, downloads)