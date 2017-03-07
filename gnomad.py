import itertools

from config import defaults
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