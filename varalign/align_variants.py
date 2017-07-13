import argparse
import logging
from operator import itemgetter

import alignments
import ensembl
import gnomad
from config import defaults
from gnomad import tabulate_variant_effects
from numpy import vectorize
import pandas as pd
from tqdm import tqdm
from time import sleep

from Bio import AlignIO, SeqIO


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


def map_uniprot_to_genome(uniprot, species='homo_sapiens', collapse=True):
    """
    Map a UniProt entry to the genome.

    :param uniprot:
    :param species:
    :param collapse:
    :return:
    """
    # Map to Gene IDs
    ensembl_genes = ensembl.get_xrefs(uniprot, species=species, features='gene')  #TODO: This could be transcrpts or translations...
    if len(ensembl_genes) == 0:
        return None  # no mapping
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


def align_variants(alignment_info, species='HUMAN'):
    """

    :param alignment_info:
    :param species:
    :return:
    """
    # ----- Map sequences to genome -----
    data = alignment_info.loc[alignment_info['species'] == species, ['seq_id', 'uniprot_id']]
    # TODO: If get transcript ID can use to filter variant table
    genomic_ranges = [
        (row.uniprot_id, row.seq_id, map_uniprot_to_genome(row.uniprot_id, species=species))
        for row in tqdm(data.itertuples(), total=len(data))
    ]
    genomic_ranges = [x for x in genomic_ranges if x[2] is not None]  # Filter unmapped sequences
    log.info("Mapped {} sequences to genome.".format(len(genomic_ranges)))

    # Add ranges to alignment info
    genomic_mapping_table = pd.DataFrame(genomic_ranges, columns=['uniprot_id', 'seq_id', 'genomic_ranges'])
    alignment_info = alignment_info.merge(genomic_mapping_table, on=['uniprot_id', 'seq_id'], how='left')

    # ----- Fetch variants for the mapped genomic ranges -----
    sequence_variant_lists = []
    for uniprot_id, seq_id, ranges in tqdm(genomic_ranges):
        sequence_variant_lists.append((seq_id, [x for r in ranges for x in gnomad.gnomad.fetch(*r)]))

    # Parse Variants to table
    all_variants = [(variant, seq_id) for seq_id, l in sequence_variant_lists for variant in l]
    variants_table = gnomad.vcf_row_to_table(*zip(*all_variants))

    # ----- Add source UniProt identifiers to the table -----
    # Create UniProt ID series that shares an index with the variant table
    source_uniprot_ids = alignment_info.set_index('seq_id')['uniprot_id']
    source_uniprot_ids.name = ('External', 'SOURCE_ACCESSION')
    source_uniprot_ids.index.name = 'SOURCE_ID'

    # Add IDs to variant tables
    variants_table = variants_table.join(source_uniprot_ids)

    # return variants_table

    # ----- Filter variant table -----
    # See version 1 of notebook for other ideas (e.g. Protin_position in UniProt range...)
    # Reduce transcript duplication
    is_canonical = variants_table[('VEP', 'CANONICAL')] == 'YES'
    is_ccds = variants_table[('VEP', 'CCDS')] != ''

    # Only want those that can map to a residue
    is_protein_coding = variants_table[('VEP', 'BIOTYPE')] == 'protein_coding'
    at_protein_position = variants_table[('VEP', 'Protein_position')] != ''

    # Filter least useful effects
    is_not_modifier = variants_table[('VEP', 'IMPACT')] != 'MODIFIER'

    # Source protein filter
    swissprot_matches_source = (variants_table['External', 'SOURCE_ACCESSION'] == variants_table['VEP', 'SWISSPROT'])
    vcontains = vectorize(lambda x, y: x in y)
    trembl_matches_source = vcontains(variants_table[('External', 'SOURCE_ACCESSION')],
                                      variants_table[('VEP', 'TREMBL')])

    trembl_matches_source[:] = False  # OVERRIDE TREMBL TO KEEP ONLY SWISSPROT

    # Apply filter
    filtered_variants = variants_table[is_canonical & is_protein_coding & is_not_modifier & is_ccds &
                                       (swissprot_matches_source | trembl_matches_source) &
                                       at_protein_position]
    log.info('Redundant rows:\t{}'.format(sum(filtered_variants.reset_index('Feature').index.duplicated())))
    log.info('Total rows:\t{}'.format(len(filtered_variants)))

    # return filtered_variants

    # ----- Map variants to columns -----
    mapping_table = pd.DataFrame(alignment_info['mapping'].tolist(),
                                 index=[alignment_info.index, alignment_info['seq_id']])  # From top of notebook
    mapping_table.reset_index(inplace=True)
    mapping_table = pd.melt(mapping_table, id_vars=['level_0', 'seq_id'])
    mapping_table.dropna(subset=['value'], inplace=True)

    # Used to include alignment row in index, see **aligned_variants_dev1**
    indexed_map_table = pd.DataFrame(mapping_table['value'].tolist(),
                                     columns=['Column', 'Protein_position'],
                                     index=[mapping_table['seq_id']]).reset_index()
    indexed_map_table = indexed_map_table.set_index(['seq_id', 'Protein_position']).sort_index()
    indexed_map_table.index.rename(['SOURCE_ID', 'Protein_position'], inplace=True)
    indexed_map_table.columns = pd.MultiIndex.from_tuples([('Alignment', x) for x in indexed_map_table.columns],
                                                          names=['Type', 'Field'])
    # Coerce Protein_position
    filtered_variants.loc[:, ('VEP', 'Protein_position')] = pd.to_numeric(
        filtered_variants[('VEP', 'Protein_position')],
        errors='coerce')
    # Set index
    filtered_variants.reset_index(['SITE', 'ALLELE_NUM', 'Feature'], inplace=True)
    filtered_variants.set_index(('VEP', 'Protein_position'), append=True, inplace=True)
    filtered_variants.index.set_names(['SOURCE_ID', 'Protein_position'], inplace=True)
    filtered_variants.sort_index(inplace=True)

    # Merge
    alignment_variant_table = indexed_map_table.join(filtered_variants)  # Drops variants that map outside alignment
    alignment_variant_table.sort_index(inplace=True)
    alignment_variant_table.head()

    return alignment_variant_table


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create variant counts column wise for an alignment.')
    parser.add_argument('alignment', type=str, help='Path to the alignment.')
    parser.add_argument('--format', type=str, default='stockholm',
                        help='Alignment format.')
    parser.add_argument('--filter', type=str,
                        help='An inclusive sequence ID filter to process only a subset of sequences.')
    args = parser.parse_args()
