import argparse
import logging
from operator import itemgetter

import pandas as pd
from Bio import AlignIO
from numpy import vectorize
from tqdm import tqdm

import alignments
import ensembl
import gnomad
from config import defaults
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


def _map_uniprot_to_genome(uniprot, species='homo_sapiens', collapse=True):
    """
    Map a UniProt entry to the genome.

    :param uniprot:
    :param species:
    :param collapse:
    :return:
    """
    # Map to Gene IDs
    ensembl_genes = ensembl.get_xrefs(uniprot, species=species,
                                      features='gene')  # TODO: This could be transcrpts or translations...
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


def _default_variant_filter(variants_table):
    """
    Apply standard filters to a alignment derived variant table.

    :param variants_table:
    :return:
    """
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
    return filtered_variants


def _mapping_table(alignment_info):
    """
    Construct a alignment column to sequence residue mapping table.

    :param alignment_info:
    :return:
    """
    mapping_table = pd.DataFrame(alignment_info['mapping'].tolist(),
                                 index=[alignment_info.index, alignment_info['seq_id']])  # From top of notebook
    mapping_table.reset_index(inplace=True)
    mapping_table = pd.melt(mapping_table, id_vars=['level_0', 'seq_id'])
    mapping_table.dropna(subset=['value'], inplace=True)
    indexed_map_table = pd.DataFrame(mapping_table['value'].tolist(),
                                     columns=['Column', 'Protein_position'],
                                     index=[mapping_table['seq_id']]).reset_index()
    indexed_map_table = indexed_map_table.set_index(['seq_id', 'Protein_position']).sort_index()
    indexed_map_table.index.rename(['SOURCE_ID', 'Protein_position'], inplace=True)
    indexed_map_table.columns = pd.MultiIndex.from_tuples([('Alignment', x) for x in indexed_map_table.columns],
                                                          names=['Type', 'Field'])
    return indexed_map_table


def align_variants(alignment, species='HUMAN'):
    """

    :param alignment_info:
    :param species:
    :return:
    """
    # ----- Parse alignment info -----
    alignment_info = alignments.alignment_info_table(alignment, species)

    # ----- Map sequences to genome -----
    # TODO: If get transcript ID can use to filter variant table
    genomic_ranges = [
        (row.seq_id, _map_uniprot_to_genome(row.uniprot_id, species=species))
        for row in tqdm(alignment_info.itertuples(), total=len(alignment_info))
    ]
    log.info("Mapped {} sequences to genome.".format(len(genomic_ranges)))

    # Add ranges to alignment info
    genomic_mapping_table = pd.DataFrame(genomic_ranges, columns=['seq_id', 'genomic_ranges'])
    alignment_info = alignment_info.merge(genomic_mapping_table, on=['seq_id'], how='left')

    # ----- Fetch variants for the mapped genomic ranges -----
    sequence_variant_lists = [(row.seq_id, (x for r in row.genomic_ranges for x in gnomad.gnomad.fetch(*r)))
                              for row in alignment_info.dropna(subset=['genomic_ranges']).itertuples()]
    all_variants = ((variant, seq_id) for seq_id, range_reader in tqdm(sequence_variant_lists)
                    for variant in range_reader)
    variants_table = gnomad.vcf_row_to_table(*zip(*all_variants))

    # ----- Add source UniProt identifiers to the table -----
    # Create UniProt ID series that shares an index with the variant table
    source_uniprot_ids = alignment_info.set_index('seq_id')['uniprot_id']
    source_uniprot_ids.name = ('External', 'SOURCE_ACCESSION')
    source_uniprot_ids.index.name = 'SOURCE_ID'
    # Add IDs to variant tables
    variants_table = variants_table.join(source_uniprot_ids)

    # ----- Filter variant table -----
    filtered_variants = _default_variant_filter(variants_table)
    log.info('Redundant rows:\t{}'.format(sum(filtered_variants.reset_index('Feature').index.duplicated())))
    log.info('Total rows:\t{}'.format(len(filtered_variants)))

    # ----- Map variants to columns -----
    # Generate alignment column / sequence residue mapping table
    indexed_map_table = _mapping_table(alignment_info)
    # Coerce Protein_position to correct type
    filtered_variants.loc[:, ('VEP', 'Protein_position')] = pd.to_numeric(
        filtered_variants.loc[:, ('VEP', 'Protein_position')],
        errors='coerce')
    # Set index for merge
    filtered_variants.reset_index(['SITE', 'ALLELE_NUM', 'Feature'], inplace=True)
    filtered_variants.set_index(('VEP', 'Protein_position'), append=True, inplace=True)
    filtered_variants.index.set_names(['SOURCE_ID', 'Protein_position'], inplace=True)
    filtered_variants.sort_index(inplace=True)
    # Merge to map
    alignment_variant_table = indexed_map_table.join(filtered_variants)  # Drops variants that map outside alignment
    alignment_variant_table.sort_index(inplace=True)

    return alignment_info, alignment_variant_table


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create variant counts column wise for an alignment.')
    parser.add_argument('alignment', type=str, help='Path to the alignment.')
    parser.add_argument('--format', type=str, default='stockholm',
                        help='Alignment format.')
    parser.add_argument('--filter', type=str,
                        help='An inclusive sequence ID filter to process only a subset of sequences.')
    args = parser.parse_args()
