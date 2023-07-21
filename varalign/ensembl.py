import os
import sys
import time
import pandas as pd
import logging
import requests
import requests_cache

from varalign.config import defaults

default_server = defaults.api_ensembl

# Globals for rate-limiting
reqs_per_sec = 15
req_count = 0
last_req = 0

log = logging.getLogger(__name__)

standard_regions = ('1', '2', '3', '4', '5', '6', '7', '8', '9', '10',
                    '11', '12', '13', '14', '15', '16', '17', '18', '19', '20',
                    '21', '22', 'X', 'Y')


def ratelimit():
    """
    Check if we need to rate limit ourselves.

    :return:
    """
    global req_count
    global reqs_per_sec
    global last_req
    if req_count >= reqs_per_sec:
        delta = time.time() - last_req
        if delta < 1:
            time.sleep(1 - delta)
        last_req = time.time()
        req_count = 0


def update_ratelimit(response):
    """
    Update parameters used for ratelimiting after a request.

    :param response:
    :return:
    """
    # No need to increment if response taken from cache.
    if not response.from_cache:
        global req_count
        req_count += 1


def get_xrefs(query_id, species='homo_sapiens', features=('gene', 'transcript', 'translation'), server=default_server):
    """
    Lookup EnsEMBL xrefs for an external ID and get feature IDs.
    """
    ratelimit()

    endpoint = "/xrefs/symbol"
    ext = '/'.join([endpoint, species, query_id]) + "?"

    with requests_cache.CachedSession(os.path.join('.varalign', 'ensembl_cache')) as s:
        r = s.get(server+ext, headers={"Content-Type": "application/json"})

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    update_ratelimit(r)

    return [x['id'] for x in r.json() if x['type'] in features]

### ADDED BY JSU

def get_canonical_transcript(stable_id):
    """
    returns the canonical transctipt ID for a stable ENSEMBL gene id.
    """
    feats = pd.read_json("https://grch37.rest.ensembl.org/overlap/id/{}?content-type=application/json;feature=gene".format(stable_id), convert_axes = False)
    gene = feats.query('id == @stable_id').copy()
    canonical_id = gene.canonical_transcript.tolist()[0].split(".")[0] # keep id, not version
    prot_feats = pd.read_json("https://grch37.rest.ensembl.org/overlap/id/{}?content-type=application/json;feature=cds".format(stable_id), convert_axes = False) # new
    cannon_prot = prot_feats.query('Parent == @canonical_id').drop_duplicates("protein_id") # new
    if len(cannon_prot) > 1:
        log.warning("There are {} protein ids for {}".format(str(len(cannon_prot)), stable_id))
    prot_canon_id = cannon_prot.protein_id.tolist()[0] # new
    return prot_canon_id

def get_genomic_range_JSU(canonical_transcript, region, server=default_server):
    """
    Get the genomic range for an EnsEMBL gene or transcript.
    """
    #ratelimit()
    start, end = region
    endpoint = "/map/translation" # "/map/cds" is wrong
    ext = '/'.join([endpoint, canonical_transcript, "{}-{}".format(str(start), str(end))]) + "?"

    #print(server+ext)
    with requests_cache.CachedSession(os.path.join('.varalign', 'ensembl_cache')) as s:
        r = s.get(server+ext, headers={"Content-Type": "application/json"})

    if not r.ok:
        log.error("Could not get genomic ranges for CDS {} region {}-{}".format(canonical_transcript, str(region[0]), str(region[1])))
        return None, 0, 0
        #r.raise_for_status()
        #sys.exit()
        #return tuple()

    decoded = r.json()["mappings"]
    decoded_filt = []
    for el in decoded:
        if el["strand"] == 0:
            log.error("Strand is 0 for {}".format(el))
        else:
            decoded_filt.append(el)
    if len(decoded_filt) == 0:
        return ('0', 0, 0)
    #decoded = [el for el in decoded if el["strand"] != 0 else log.error("")]
    if decoded_filt[0]["strand"] == 1:
        prot_start = decoded_filt[0]["start"]
        prot_end = decoded_filt[-1]["end"]
    elif decoded_filt[0]["strand"] == -1:
        prot_start = decoded_filt[-1]["start"]
        prot_end = decoded_filt[0]["end"]
    if prot_start > prot_end:
        log.error("start {} >= end {} for {}".format(str(prot_start), str(prot_end), canonical_transcript))
    #update_ratelimit(r)
    # add condition to check all seq_region_names are the same
    #return str(decoded['seq_region_name']), decoded['start'], decoded['end']
    return str(decoded_filt[0]['seq_region_name']), prot_start, prot_end

### ADDED BY JSU



def get_genomic_range(query_id, server=default_server):
    """
    Get the genomic range for an EnsEMBL gene or transcript.
    """
    ratelimit()

    endpoint = '/lookup/id'
    ext = '/'.join([endpoint, query_id]) + "?"

    with requests_cache.CachedSession(os.path.join('.varalign', 'ensembl_cache')) as s:
        r = s.get(server+ext, headers={"Content-Type": "application/json"})

    if not r.ok:
        r.raise_for_status()
        sys.exit()

    decoded = r.json()

    update_ratelimit(r)

    return str(decoded['seq_region_name']), decoded['start'], decoded['end']


def merge_ranges(ranges, min_gap=150):
    """
    Merge a set of genomic ranges into non-overlapping sets.

    :param ranges:
    :param min_gap: Minimum gap to allow between consecutive ranges.
    :return:
    """
    ranges = ranges[:]
    ranges.sort(key=lambda a: a[2])  # Sort by end
    ranges.sort(key=lambda a: a[1])  # Sort by start
    ranges.sort(key=lambda a: a[0])  # Sort by region

    new_ranges = [list(ranges.pop(0))]
    for region, start, end in ranges:
        if region == new_ranges[-1][0]:
            if start <= new_ranges[-1][2] + min_gap:
                new_ranges[-1][2] = end  # Merge range
            else:
                new_ranges.append([region, start, end])  # New range on same region
        else:
            new_ranges.append([region, start, end])  # New range on different region

    return [tuple(x) for x in new_ranges]
