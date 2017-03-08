import argparse
import logging
from operator import itemgetter

import ensembl
import gnomad
from gnomad import tabulate_variant_effects

log = logging.getLogger(__name__)
log.setLevel('INFO')


def _fetch_variants_for_uniprot(uniprot, canonical=True, consequences=('missense_variant', 'synonymous_variant'),
                                only_snvs=True):
    """
    Retrieve variants from Gnomad given a UniProt ID.

    :param sequence:
    :return:
    """
    # Map UniProt to human genome
    ensembl_genes = ensembl.get_xrefs(uniprot, features='gene')
    assert(len(ensembl_genes) >= 1)  # Check for at least one mapping
    ensembl_ranges = [ensembl.get_genomic_range(x) for x in ensembl_genes]
    for i in range(len(ensembl_genes)):
        log.info('Mapped {} to {} on chr: {}, {}-{}'.format(uniprot, ensembl_genes[i], *ensembl_ranges[i]))
    # Identify and remove non-standard sequence regions
    non_standard_ranges = [i for i, x in enumerate(ensembl_ranges) if x[0] not in ensembl.standard_regions]
    if len(non_standard_ranges) > 0:
        non_standard_ranges.sort(reverse=True)
        [ensembl_genes.pop(i) for i in non_standard_ranges]
        [ensembl_ranges.pop(i) for i in non_standard_ranges]
        log.info('Removed %s non-standard sequence regions.', len(non_standard_ranges))
    if len(ensembl_ranges) == 0:
        raise ValueError('Could not map {} to the genome.'.format(uniprot))

    # Lookup variants
    log.info('Retrieving variants...')
    lookup_ranges = ensembl.merge_ranges(ensembl_ranges, min_gap=1000)
    variants = [x for _range in lookup_ranges for x in gnomad.gnomad.fetch(*_range)]  # TODO: Add progress bar?
    log.info('Found {} variants.'.format(len(variants)))

    # filter variants on VEP
    vep_table = tabulate_variant_effects(variants)
    query = 'SWISSPROT == @uniprot & Consequence in @consequences'
    # TODO: Refactor optional filter query to config?
    if canonical:
        query += ' & CANONICAL == "YES"'
    if only_snvs:
        query += ' & VARIANT_CLASS == "SNV"'
    vep_table.query(query, inplace=True)
    assert not any(vep_table['Allele'].reset_index().duplicated())  # Assume this filter gives one effect per variant allele
    variants = itemgetter(*vep_table.index)(variants)
    vep_table.reset_index(drop=True, inplace=True)  # Reset index to match filtered variant list
    assert len(variants) == len(vep_table)
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
