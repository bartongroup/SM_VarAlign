from Bio import SeqIO
from config import defaults
import pandas as pd
import io
import os
import requests
from utils import urlopen_with_retry


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
        handle = open(local_path, 'r')
        return SeqIO.read(handle, "fasta")
    else:
        # Unrevised sequences don't allow version number
        if _has_version(uniprot_id):
            version = _get_version(uniprot_id)
            if version == '1':
                current_version = _get_latest_revision_number(uniprot_id)
                if current_version == 1:
                    uniprot_id = _strip_version(uniprot_id)
            else:
                uniprot_id = _strip_version(uniprot_id)
                history = _get_entry_history(uniprot_id)
                matches = history['Sequence version'] == int(version)
                most_recent = history[matches]['Entry version'].max()
                uniprot_id = _strip_version(uniprot_id) + '.' + str(most_recent)

        url = defaults.api_uniprot + uniprot_id + '.fasta'
        handle = urlopen_with_retry(url)
        return SeqIO.read(handle, "fasta")
