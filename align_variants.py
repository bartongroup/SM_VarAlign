import argparse
from config import defaults
import logging
from operator import itemgetter

import ensembl
import gnomad
from gnomad import tabulate_variant_effects

log = logging.getLogger(__name__)
log.setLevel('INFO')


def _build_vep_filter(canonical=eval(defaults.canonical), consequences=defaults.consequences,
                      additional=defaults.additional):
    """
    Build a VEP filter string suitable for DataFrame.query().

    :param canonical:
    :param consequences:
    :param additional:
    :return:
    """
    query = []
    if consequences != ['']:
        query.append('Consequence in {}'.format(consequences))
    if canonical:
        query.append('CANONICAL == "YES"')
    if additional != '':
        query.append(additional)
    return ' & '.join(query)


def _fetch_variants_for_uniprot(uniprot, filter_swissprot=True):
    """
    Retrieve variants from Gnomad given a UniProt ID.

    :param uniprot:
    :param filter_swissprot:
    :return:
    """
    # Map UniProt to genome
    lookup_ranges = map_uniprot_to_genome(uniprot)

    # Lookup variants
    log.info('Retrieving variants...')
    variants = [x for _range in lookup_ranges for x in gnomad.gnomad.fetch(*_range)]  # TODO: Add progress bar?
    log.info('Found {} variant sites.'.format(len(variants)))

    # filter variants on VEP
    vep_table = tabulate_variant_effects(variants)
    query = _build_vep_filter()
    if filter_swissprot:
        query += ' & SWISSPROT == "{}"'.format(uniprot)
    log.info('Keeping variants where: %s', query)
    vep_table.query(query, inplace=True)
    if vep_table.empty:
        raise ValueError('No variants pass filter.')
    # Assume this filter gives one effect per variant allele
    assert not any(vep_table['Allele'].reset_index().duplicated())
    variants = itemgetter(*vep_table.index)(variants)  # Filter variant record by VEP filter
    assert len(variants) == len(vep_table)
    vep_table.reset_index(drop=True, inplace=True)  # Reset index to match filtered variant list
    variants = [gnomad.split_variant(variant, allele)[0] for variant, allele in zip(variants, vep_table.Allele)]
    log.info('Returning {} variants after filtering.'.format(len(variants)))

    return variants, vep_table


def map_uniprot_to_genome(uniprot, collapse=True):
    """
    Map a UniProt entry to the genome.

    :param uniprot:
    :param collapse:
    :return:
    """
    # Map to Gene IDs
    ensembl_genes = ensembl.get_xrefs(uniprot, features='gene')  #TODO: This could be transcrpts or translations...
    assert (len(ensembl_genes) >= 1)  # Check for at least one mapping
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
    # Check for no mapping
    if len(ensembl_ranges) == 0:
        raise ValueError('Could not map {} to the genome.'.format(uniprot))
    # Collapse ranges if desired
    if collapse:
        ensembl_ranges = ensembl.merge_ranges(ensembl_ranges, min_gap=1000)
    return ensembl_ranges


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create variant counts column wise for an alignment.')
    parser.add_argument('alignment', type=str, help='Path to the alignment.')
    parser.add_argument('--format', type=str, default='stockholm',
                        help='Alignment format.')
    parser.add_argument('--filter', type=str,
                        help='An inclusive sequence ID filter to process only a subset of sequences.')
    args = parser.parse_args()
