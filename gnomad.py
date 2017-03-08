import itertools

from config import defaults
from copy import deepcopy
import logging
import pandas as pd
import vcf  ## Requires pysam functionality


log = logging.getLogger(__name__)


log.info('Loading %s...', defaults.gnomad)
gnomad = vcf.Reader(filename=defaults.gnomad)

# Parse VEP annotation format
CSQ_Format = gnomad.infos['CSQ'].desc.split(' Format: ')[1].split('|')


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


def split_variant(variant, alleles=[], exclude=('CSQ',)):
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
        new_variant.INFO = {k: [v[allele_index], ] if isinstance(v, list) else v for k, v in info.items()}
        return new_variant

    # Work out alleles to process, either all or selection by ALT allele base or index
    if alleles == ():
        alleles = range(len(variant.ALT))
    else:
        if isinstance(alleles, int):
            pass
        if isinstance(alleles[0], str):  # Will work with single string or list
            alleles = [variant.ALT.index(x) for x in alleles]

    # Extract desired alleles
    allelic_records = []
    for i in alleles:
        new_variant = _extract_allele(variant, i, exclude)
        allelic_records.append(new_variant)

    return allelic_records
