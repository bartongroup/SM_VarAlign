import itertools

from config import defaults
from copy import deepcopy
import logging
import pandas as pd
import vcf  ## Requires pysam functionality
from vcf.utils import trim_common_suffix


log = logging.getLogger(__name__)

# Load the VCF and log some key metadata
log.info('Loading %s...', defaults.gnomad)
gnomad = vcf.Reader(filename=defaults.gnomad)
log.info('File format: {}'.format(gnomad.metadata.get('fileformat')))
log.info('Reference: {}'.format(gnomad.metadata.get('reference')))

# Parse VEP annotation format
CSQ_Format = gnomad.infos['CSQ'].desc.split(' Format: ')[1].split('|')
log.info('CSQ Format: {}'.format(CSQ_Format))

# Parse INFO header
info_header = pd.DataFrame([v for k, v in gnomad.infos.items()], dtype='str')
info_flag_num = 0
info_value_num = 1
info_allele_num = -1
info_flag_fields = info_header.query('num == @info_flag_num').id.tolist()
info_value_fields = info_header.query('num == @info_value_num').id.tolist()
info_allele_fields = info_header.query('num == @info_allele_num').id.tolist()
log.info('INFO flags: {}'.format(str(info_flag_fields)))
log.info('INFO value: {}'.format(str(info_value_fields)))
log.info('INFO per allele: {}'.format(str(info_allele_fields)))

# Annotations that need special handling during variant allele expansion
standard_num_values = [info_flag_num, info_value_num, info_allele_num]
special_handling = {'INFO': list(info_header.query('num not in @standard_num_values').id)}
msg = 'The following INFO fields cannot be resolved to a variant and will be stripped when melting variants: {}'
log.info(msg.format(str(special_handling['INFO'])))

# Parse FORMAT header
format_header = pd.DataFrame([v for k, v in gnomad.formats.items()], dtype='str')
log.info('FORMAT fields: {}'.format(str(list(format_header.id))))

# Parse FILTER header
filter_header = pd.DataFrame([v for k, v in gnomad.filters.items()], dtype='str')
log.info('FILTER flags: {}'.format(str(list(filter_header.id))))


def get_vep_annotation(record, fields=CSQ_Format):
    """
    Retrieve and format the VEP annotation from a VCF record.
    """
    parsed = pd.DataFrame([x.split('|') for x in record.INFO['CSQ']],
                          columns=CSQ_Format)
    if len(fields) < len(CSQ_Format):
        parsed = parsed[fields]

    return parsed


def get_vep_raw(record):
    """
    Retrive VEP annotation from VCF record with minimal parsing.
    """
    return [x.split('|') for x in record.INFO['CSQ']]


def tabulate_variant_effects(variants):
    """
    Efficiently parse the VEP CSQ field from a list of VCF records into a table.

    :param variants: List of VCF records.
    :return: DataFrame of VEP annotations indexed to variant list.
    """
    veps = [get_vep_raw(x) for x in variants]
    variant_indices = [i for i, x in enumerate(veps) for _ in x]  # index each vep entry to the variant record
    vep_table = pd.DataFrame(list(itertools.chain.from_iterable(veps)),
                             columns=CSQ_Format, index=variant_indices)
    return vep_table


def tabulate_variant_info(variants, split=False):
    """
    Efficiently parse INFO from a list of VCF records into a table.

    :param variants: List of VCF records.
    :param split: Whether to return seperate tables for site and allele annotations.
    :return: DataFrame of INFO annotations indexed to variant list.
    """
    info_table = pd.DataFrame([x.INFO for x in variants])
    if split:
        return _split_info_table(info_table)

    else:
        return info_table


def _split_info_table(info_table):
    """
    Split an INFO table into separate tables for site and allele level annotations.
    :param info_table: 
    :return: 
    """
    available = info_table.columns
    allele_fields = [x for x in info_allele_fields if x in available]
    site_fields = [x for x in info_value_fields + info_flag_fields if x in available]
    return info_table[site_fields], info_table[allele_fields]


def split_variant(variant, alleles=[], exclude=special_handling['INFO'], value_only=False):
    """
    Split a multiallelic variant.

    :param variant:
    :param alleles: Alleles to be extracted, either by ALT allele or ALT allele index (default: all alleles).
    :param exclude: INFO fields to be dropped from the new variant records.
    :return:
    """
    def _extract_allele(variant, allele_index, exclude):
        new_variant = deepcopy(variant)
        # ALT field
        new_variant.ALT = [new_variant.ALT[allele_index], ]
        # INFO field
        info = new_variant.INFO
        [info.pop(x, 0) for x in exclude]  # Exclude ignored info fields
        if value_only:
            new_variant.INFO = {k: v[allele_index] if isinstance(v, list) else v for k, v in info.items()}
        else:
            new_variant.INFO = {k: [v[allele_index], ] if isinstance(v, list) else v for k, v in info.items()}
        return new_variant

    # Work out alleles to process, either all or selection by ALT allele base or index
    if alleles == ():
        alleles = range(len(variant.ALT))
    else:
        if isinstance(alleles, int):
            pass
        if isinstance(alleles[0], str):  # Will work with single string or list
            try:
                alleles = [variant.ALT.index(x) for x in alleles]
            except ValueError:
                # Happens when VEP allele doesn't match variant
                # Could occur because a multiallelic site contains an indel, MNP, and/or a SNP
                trimmed_alts = [trim_common_suffix(str(x), variant.REF)[0] for x in variant.ALT]
                alleles = [trimmed_alts.index(x) for x in alleles]

    # Extract desired alleles
    allelic_records = []
    for i in alleles:
        new_variant = _extract_allele(variant, i, exclude)
        allelic_records.append(new_variant)

    return allelic_records


def vcf_row_to_table(variants):
    """
    Parse PyVCF Variant to a DataFrame.
    
    
    :param variants: 
    :return: 
    """

    print 'Processing {} variants...'.format(len(variants))

    # Build record table
    print 'Constructing site record table...'
    records = []
    for site, variant in enumerate(variants):
        for allele, alt in enumerate(variant.ALT):
            records.append([site, allele, variant.ID, variant.CHROM, variant.POS, variant.REF, alt, variant.QUAL,
                            variant.FILTER])
    row_record = pd.DataFrame(records, columns=['SITE', 'ALLELE_NUM', 'ID', 'CHROM', 'POS', 'REF', 'ALT', 'QUAL',
                                                'FILTER'])
    row_record.set_index(['SITE', 'ALLELE_NUM'], inplace=True)

    # VEP table
    print 'Tabulating VEP records...'
    vep_table = tabulate_variant_effects(variants)
    vep_table['ALLELE_NUM'] = vep_table['ALLELE_NUM'].astype(int)
    vep_table['ALLELE_NUM'] = vep_table['ALLELE_NUM'] - 1
    vep_table.index.name = 'SITE'
    vep_table.set_index('ALLELE_NUM', inplace=True, append=True)  # NB. Exclude 'Feature' to join on index later

    # INFO tables
    print 'Processing INFO fields...'
    site_info, allele_info = tabulate_variant_info(variants, split=True)
    site_info.index.name = 'SITE'
    # Split allele INFO fields
    tables = []
    for site_num, variant in enumerate(variants):
        info_dict = variant.INFO
        allele_info_records = {k: info_dict[k] for k in info_allele_fields if k in info_dict}
        fields, both_allele_values = zip(*allele_info_records.items())
        for allele_num, single_allele_values in enumerate(zip(*both_allele_values)):
            allele_entry = dict(zip(fields, single_allele_values))
            allele_entry['ALLELE_NUM'] = allele_num
            allele_entry['SITE'] = site_num
            tables.append(allele_entry)
    split_allele_info = pd.DataFrame(tables).set_index(['SITE', 'ALLELE_NUM'])

    # Create MultiIndex for sub-table columns
    row_record.columns = pd.MultiIndex.from_tuples([('Row', x) for x in row_record.columns],
                                                   names=['Type', 'Field'])
    vep_table.columns = pd.MultiIndex.from_tuples([('VEP', x) for x in vep_table.columns],
                                                  names=['Type', 'Field'])
    site_info.columns = pd.MultiIndex.from_tuples([('Site_INFO', x) for x in site_info.columns],
                                                  names=['Type', 'Field'])
    split_allele_info.columns = pd.MultiIndex.from_tuples([('Allele_INFO', x) for x in split_allele_info.columns],
                                                          names=['Type', 'Field'])

    # Merge all the sub-tables
    print 'Joining sub-tables...'
    # 1. Merged at site level
    merged_variant_table = row_record.join(site_info)
    # 2. Merged at allele level
    merged_variant_table = merged_variant_table.join(split_allele_info)
    merged_variant_table = merged_variant_table.join(vep_table)
    # 3. Add Feature to index
    print 'Formatting result...'
    merged_variant_table.set_index(('VEP', 'Feature'), append=True, inplace=True)
    merged_variant_table.index.set_names(['SITE', 'ALLELE_NUM', 'Feature'], inplace=True)
    # merged_variant_table = merged_variant_table.reorder_levels(['SITE', 'ALLELE_NUM', 'Feature'])
    # merged_variant_table.sort_index(inplace=True)

    return merged_variant_table
