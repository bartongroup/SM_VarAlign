from Bio import SeqIO
from config import defaults
import pandas as pd
import io
import logging
import os
import requests
from utils import urlopen_with_retry

log = logging.getLogger(__name__)

def _has_version(uniprot_id):
    """
    Tell is a UniProt ID has a version or not.

    :param uniprot_id:
    :return:
    """
    return '.' in uniprot_id


def _get_version(uniprot_id):
    """
    Get version number from versioned UniProt ID.

    :param uniprot_id:
    :return:
    """
    return uniprot_id.split('.')[1]


def _strip_version(uniprot_id):
    """
    Strips the version number from a UniProt ID is present

    :param uniprot_id:
    :return:
    """
    return uniprot_id.split('.')[0]


def _get_entry_history(uniprot_id):
    """

    Example URL:
    http://www.uniprot.org/uniprot/A0A0B4J1V8.tab?version = *

    :param uniprot_id:
    :return:
    """
    url = defaults.api_uniprot + _strip_version(uniprot_id) + '.tab?version=*'
    s = requests.get(url).content
    return pd.read_table(io.StringIO(s.decode('utf-8')))


def _get_latest_revision_number(uniprot_id):
    """
    Get the number of the latest sequence revision for a UniProt entry.

    :param uniprot_id:
    :return:
    """
    history = _get_entry_history(_strip_version(uniprot_id))
    return history['Sequence version'].max()


def get_uniprot_fasta(uniprot_id):
    """
    Robust method to connect to a specific, possibly versioned, UniProt entry with caching.

    :param uniprot_id:
    :return:
    """
    local_path = os.path.join(defaults.uniprot_cache, 'sequences', uniprot_id + '.fasta')
    if os.path.isfile(local_path):
        # Reload from cache
        log.info('Loading {}...'.format(local_path))
        handle = open(local_path, 'r')
        return SeqIO.read(handle, "fasta")
    else:
        # Format URL for sequence version retrieval
        if _has_version(uniprot_id):
            sequence_version = _get_version(uniprot_id)
            uniprot_id = _strip_version(uniprot_id)
            # Match required sequence version to latest entry version that is equivalent
            history = _get_entry_history(uniprot_id)
            matches = history['Sequence version'] == int(sequence_version)
            matched_entry_version = history[matches]['Entry version'].max()
            url = defaults.api_uniprot + uniprot_id + '.fasta' + '?version=' + str(matched_entry_version)
        else:
            url = defaults.api_uniprot + uniprot_id + '.fasta'
        handle = urlopen_with_retry(url)
        sequence = SeqIO.read(handle, "fasta")
        # Write to cache
        if not os.path.exists(os.path.dirname(local_path)):
            os.makedirs(os.path.dirname(local_path))
        SeqIO.write(sequence, local_path, "fasta")
        return sequence
