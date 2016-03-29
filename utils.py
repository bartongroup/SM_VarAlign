import urllib2
import requests
from retry import retry
import logging

log = logging.getLogger(__name__)

# http://stackoverflow.com/questions/9446387/how-to-retry-urllib2-request-when-fails
@retry(urllib2.URLError, tries=4, delay=3, backoff=2)
def urlopen_with_retry(url):
    return urllib2.urlopen(url)

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