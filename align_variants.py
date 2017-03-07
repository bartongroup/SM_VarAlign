import argparse
import alignments
from Bio import AlignIO, SeqIO
import ensembl
import gnomad
import itertools
import logging
from operator import itemgetter
import pandas as pd
import sys


log = logging.getLogger(__name__)
log.setLevel('INFO')


def _fetch_variants_for_uniprot(uniprot, canonical=True, consequences=('missense_variant', 'synonymous_variant')):
    """
    Retrieve variants from Gnomad given a UniProt ID.

    :param sequence:
    :return:
    """
    # Map UniProt to human genome
    ensembl_gene = ensembl.get_xrefs(uniprot, features='gene')
    assert(len(ensembl_gene) == 1)  ## Check for unique mapping
    ensembl_range = ensembl.get_genomic_range(ensembl_gene[0])
    log.info('Mapped {} to {} on chr: {}, {}-{}'.format(uniprot, ensembl_gene, *ensembl_range))
    # Lookup variants and parse VEP annotation
    variants = [x for x in gnomad.gnomad.fetch(*ensembl_range)]
    veps = [gnomad.get_vep_raw(x) for x in variants]
    variant_indices = [i for i, x in enumerate(veps) for _ in x]  # index each vep entry to the variant record
    vep_table = pd.DataFrame(list(itertools.chain.from_iterable(veps)),
                             columns=gnomad.CSQ_Format, index=variant_indices)
    log.info('Found {} variants.'.format(len(variants)))
    # filter variants on VEP
    query = 'SWISSPROT == @uniprot & Consequence in @consequences'
    if canonical:
        query += ' & CANONICAL == "YES"'
    vep_table.query(query, inplace=True)
    variants = itemgetter(*vep_table.index)(variants)
    log.info('Returning {} variants after filtering.'.format(len(variants)))
    return variants, vep_table



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create variant counts column wise for an alignment.')
    parser.add_argument('alignment', type=str, help='Path to the alignment.')
    parser.add_argument('--format', type=str, default='stockholm',
                        help='Alignment format.')
    parser.add_argument('--filter', type=str,
                        help='An inclusive sequence ID filter to process only a subset of sequences.')
    args = parser.parse_args()
