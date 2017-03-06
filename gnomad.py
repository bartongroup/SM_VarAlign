from config import defaults
import pandas as pd
import vcf


gnomad = vcf.Reader(filename=defaults.gnomad)

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
