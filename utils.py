"""
Utility functions.
"""
import re
import urllib2
import requests
from retry import retry
import logging

log = logging.getLogger(__name__)

# http://stackoverflow.com/questions/9446387/how-to-retry-urllib2-request-when-fails
@retry(urllib2.URLError, tries=4, delay=3, backoff=2)
def urlopen_with_retry(url):
    return urllib2.urlopen(url)  #TODO: 404 should be handled differently


def query_uniprot(search_terms, first=False):
    """
    Query the UniProt API for proteins that have particualar characteristics.

    :param search_terms: A tuple of UniProt Query search terms.
    E.g. ('keyword:Disease', 'reviewed:yes', 'organism:human', 'database:(type:pdb)')
    :return: A list of UniProt IDs
    """
    url = 'http://www.uniprot.org/uniprot'
    params = {'query': ' AND '.join(search_terms),
              'format': 'tab', 'columns': 'id', 'sort':'score'}
    r = requests.get(url, params=params)
    uniprots = r.content.split('\n')[1:] # Drop column header
    if '' in uniprots:
        uniprots.remove('')
    if len(uniprots) == 0:
        return None
    log.info('Retreived {} UniProt IDs matching query.'.format(len(uniprots)))
    if first:
        return uniprots[0]
    else:
        return uniprots


def worse_than(SO_term):
    """
    Identify what variant effects are worse than the specified term.

    :param SO_term: Sequence Ontology term
    :return: list of SO terms 'worse than' and including the query
    """
    # http://www.ensembl.org/info/genome/variation/predicted_data.html#consequences
    ranked_terms = ('transcript_ablation', 'splice_acceptor_variant', 'splice_donor_variant',
                    'stop_gained', 'frameshift_variant', 'stop_lost', 'start_lost', 'transcript_amplification',
                    'inframe_insertion', 'inframe_deletion', 'missense_variant', 'protein_altering_variant',
                    'splice_region_variant', 'incomplete_terminal_codon_variant', 'stop_retained_variant',
                    'synonymous_variant')
    return ranked_terms[:ranked_terms.index(SO_term) + 1]


def parse_seq_name(seq_name):
    """
    Extract identifier portion of alignment sequence name.

    Sequence name strings often contain metadata (e.g., residue ranges P12345/1-21). This function returns the ID
    portion.

    :param seq_name: Alignment sequence identifier.
    :return:
    """
    return re.search('\w*', seq_name).group().strip()