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


def tabulate_variant_info(variants):
    """
    Efficiently parse INFO from a list of VCF records into a table.

    :param variants: List of VCF records.
    :return: DataFrame of INFO annotations indexed to variant list.
    """
    return pd.DataFrame([x.INFO for x in variants])


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
