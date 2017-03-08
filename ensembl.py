from config import defaults
import requests
import sys


default_server = defaults.api_ensembl

standard_regions = ('1', '2', '3', '4', '5', '6', '7', '8', '9', '10',
                    '11', '12', '13', '14', '15', '16', '17', '18', '19', '20',
                    '21', '22', 'X', 'Y')


def get_xrefs(query_id, species='homo_sapiens', features=('gene', 'transcript', 'translation'), server=default_server):
    """
    Lookup EnsEMBL xrefs for an external ID and get feature IDs.
    """
    endpoint = "/xrefs/symbol"
    ext = '/'.join([endpoint, species, query_id]) + "?"
    r = requests.get(server+ext, headers={"Content-Type": "application/json"})

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    return [x['id'] for x in r.json() if x['type'] in features]


def get_genomic_range(query_id, server=default_server):
    """
    Get the genomic range for an EnsEMBL gene or transcript.
    """
    endpoint = '/lookup/id'
    ext = '/'.join([endpoint, query_id]) + "?"
    r = requests.get(server+ext, headers={"Content-Type": "application/json"})

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    decoded = r.json()

    return str(decoded['seq_region_name']), decoded['start'], decoded['end']
