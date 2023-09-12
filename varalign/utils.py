"""
Utility functions.
"""
import copy
import logging
import os
import re

import requests
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from varalign.six.moves import urllib

from varalign.retry import retry

ALIGNMENT_CHARS = set('ACDEFGHIKLMNPQRSTVWY-')

log = logging.getLogger(__name__)
log.setLevel('INFO')


# http://stackoverflow.com/questions/9446387/how-to-retry-urllib2-request-when-fails
@retry(urllib.error.URLError, tries=4, delay=3, backoff=2)
def urlopen_with_retry(url):
    return urllib.request.urlopen(url)  #TODO: 404 should be handled differently


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


def filter_alignment(alignment, seq_id_filter):
    """

    :param alignment:
    :param seq_id_filter:
    :return:
    """
    passing_seqs = []
    for seq in alignment:
        if seq_id_filter is not None and seq_id_filter not in seq.id:
            log.info('Filtering sequence {}.'.format(seq.id))
        else:
            passing_seqs.append(seq)
    filtered_alignment = MultipleSeqAlignment(passing_seqs)
    filtered_alignment.annotations = alignment.annotations
    return filtered_alignment


def sanitise_alignment(aln):
    """

    :param alignment:
    :return:
    """
    alignment_copy = copy.deepcopy(aln)
    modified = {}
    for seqrec in alignment_copy:
        # Letter annotations have to be put aside before seq is mutated
        annots = seqrec.letter_annotations
        seqrec.letter_annotations = {}

        # Mark problem residues for the log
        modified['lowercase'] = [str(i) for i, c in enumerate(str(seqrec.seq)) if c.islower()]
        modified['X'] = [str(i) for i, c in enumerate(str(seqrec.seq)) if c == 'X']
        modified['.'] = [str(i) for i, c in enumerate(str(seqrec.seq)) if c == '.']
        modified['Z'] = [str(i) for i, c in enumerate(str(seqrec.seq)) if c == 'Z']
        modified['B'] = [str(i) for i, c in enumerate(str(seqrec.seq)) if c == 'B']

        # Sanitise seq string
        new_seq_str = str(seqrec.seq).upper()
        new_seq_str = new_seq_str.replace('X', 'G')  # Any AA
        new_seq_str = new_seq_str.replace('Z', 'E')  # Glutamine or Glutamic acid
        new_seq_str = new_seq_str.replace('B', 'D')  # Aspartic acid or Asparagine
        new_seq_str = new_seq_str.replace('.', '-')

        # Check if there's anything left weird
        unk_chars = set(new_seq_str).difference(ALIGNMENT_CHARS)
        if unk_chars:
            log.warning('Unrecognised characters ({}) remain in {}.'.format(''.join(unk_chars),
                                                                            seqrec.id))

        # Mutate seq
        seqrec.seq = Seq(new_seq_str)

        # Put back annotations
        seqrec.letter_annotations = annots

    # Log modified columns
    if modified['lowercase']:
        log.info('Fixed columns with lowercase letters: {}'.format(','.join(set(modified['lowercase']))))
    if modified['X']:
        log.info('Replaced X with G in columns: {}'.format(','.join(set(modified['X']))))
    if modified['.']:
        log.info('Replaced . with - in columns: {}'.format(','.join(set(modified['.']))))
    if modified['Z']:
        log.info('Replaced Z with E in columns: {}'.format(','.join(set(modified['Z']))))
    if modified['B']:
        log.info('Replaced B with D in columns: {}'.format(','.join(set(modified['B']))))

    return alignment_copy


def is_missense_variant(variants):
    """
    Identify missense variants.

    This is not automatically as trivial as `table.type == 'missense_variant'` because multiallelic variant are all
    given the same 'type' if bundled in a single record.

    :param variants: Variant table
    :return: Boolean mask
    """
    mask = (variants['type'] == 'missense_variant') & \
           (variants['from_aa'] != variants['to_aa_expanded'])
    return mask


def is_from_to_variant(native, mutant, variants):
    """
    Identify variants that are specific mutations.

    :param native: Wild-type residue
    :param mutant: Mutant resdiue
    :param variants: Variant table
    :return: Boolean mask
    """
    mask = (variants['from_aa'] == native) & (variants['to_aa_expanded'] == mutant)
    return mask


def is_worse_than_type(type, variants):
    """
    Filter variants annotated with SO term worse than specified.

    :param type: SO term
    :param variants: Variant table
    :return: Boolean mask
    """
    mask = variants.type.apply(lambda x: x in worse_than(type))
    return mask


def is_common_variant(variants, maf):
    """
    Classify variants according to population frequency.

    :param variants: Variant table
    :param maf: Minimum frequency to call common (float or None, such that non-singleton is common)
    :return: Boolean mask
    """
    if not maf:
        mask = variants['minor_allele_frequency'].notnull()
    elif isinstance(maf, float):
        mask = variants['minor_allele_frequency'] >= maf

    return mask


def is_non_synonomous(variants):
    """
    Identify non-synonymous variants.

    :param variants: Variant table
    :return: Boolean mask
    """
    mask = variants['from_aa'] != variants['to_aa_expanded']
    return mask


def make_dir_if_needed(path):
    try:
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise