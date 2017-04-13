import argparse
import logging
import os
from itertools import tee, izip
from subprocess import call

import pandas as pd
from Bio import AlignIO  # Needs PR #768 #769 patched
from Bio.Alphabet import IUPAC

import pfam
from utils import filter_alignment, sanitise_alignment

log = logging.getLogger(__name__)
logging.captureWarnings(True)
logging.basicConfig(level='INFO', format='%(asctime)s - %(levelname)s - %(message)s ')

def apply_column_mask(alignment, column_mask):
    """
    Apply a column mask to an family.

    :param alignment: Alignment to be masked.
    :param column_mask: Boolean mask.
    :type alignment: Bio.Align.MultipleSeqAlignment
    :type column_mask: list
    :return: Masked family.
    :rtype: Bio.Align.MultipleSeqAlignment
    """
    def pairwise(iterable):
        "s -> (s0,s1), (s1,s2), (s2, s3), ..."
        a, b = tee(iterable)
        next(b, None)
        return izip(a, b)

    # Identify boundaries between empty / non-empty and format into ranges
    is_mask_boundary = []
    for i, j in zip(column_mask[:-1], column_mask[1:]):
        is_mask_boundary.append(i != j)
    boundary_indexes = [None]  # list[None:i]
    for n, k in enumerate(is_mask_boundary):
        if k:
            boundary_indexes.append(n + 1)  # +1 for slicing
    boundary_indexes.append(None)  # list[i:None]
    ranges = pairwise(boundary_indexes)

    # Now use those block lengths to build new family
    # 1) Initialise new family
    start, stop = ranges.next()
    keep_range = not column_mask[0]
    if not keep_range:
        start, stop = ranges.next()
        keep_range = not keep_range
    masked_alignment = alignment[:, start:stop]
    keep_range = not keep_range

    # 2) Proceed through remaining ranges
    for start, stop in ranges:
        if keep_range:
            masked_alignment = masked_alignment + alignment[:, start:stop]
        keep_range = not keep_range

    return masked_alignment


if __name__ == '__main__':
    # CLI Arguments
    parser = argparse.ArgumentParser(description='Create variant counts column wise for an family.')
    parser.add_argument('--local_pfam', type=str, help='Path to PFAM-A.seed')
    parser.add_argument('--seq_id_filter', type=str,
                        help='An inclusive filter to process only a subset of sequences.')
    parser.add_argument('--families', type=str, help='Path to list of families.')
    parser.add_argument('--jabaws_jar', type=str, default='/Applications/min-jabaws-client-2.1.0.jar',
                        help='Path to JABAWS client.')
    parser.add_argument('--jabaws_server', type=str, default='http://www.compbio.dundee.ac.uk/jabaws',
                        help='URL of JABAWS server.')
    parser.add_argument('--aacons_preset', type=str, default='Quick conservation',
                        help='AACons preset.')
    args = parser.parse_args()

    # Arguments and parameters
    local_pfam = args.local_pfam
    seq_id_filter = args.seq_id_filter
    family_id_list = args.families

    # Identify array number and lookup Pfam list for corresponding AC
    job_number = int(os.getenv('SGE_TASK_ID', 1))
    log.info('Array job number: {}'.format(job_number))
    with open(family_id_list, 'r') as ids:
        for i in xrange(job_number):
            desired_family = ids.readline().strip()
    log.info('Searching Pfam for family: {}'.format(desired_family))

    # Find family with Pfam index
    index_path = local_pfam + '.idx'
    if not os.path.isfile(index_path):
        log.error('Pfam index file not found')
        raise SystemExit

    start = pfam.lookup_index(index_path, desired_family)
    if start is None:
        log.error('Pfam {} could not be found in the index'.format(desired_family))
        raise SystemExit
    family = pfam.read_family(local_pfam, start)
    family_name = family.annotations['GF']['AC'][0]  # This requires biopython patches #768 #769
    log.debug('Read family {}...'.format(family_name))


    # Store full family name
    log.info('Found family: {}'.format(family_name))
    family_name += '_' + os.path.basename(local_pfam)

    # Filter unwanted sequences
    if seq_id_filter is not None:
        log.info('Filtering family for sequences without "{}"'.format(seq_id_filter))
        family = filter_alignment(family, seq_id_filter)
        if len(family) == 0:
            log.warning('No sequences passed filter. Exiting.')
            raise SystemExit

    # TODO: Alternatively could sanitise the sequences completely
    # Some sanitisation is unambiguous
    log.info('Sanitising columns...')
    family = sanitise_alignment(family)

    # Remove empty columns and those with unknown characters
    log.info('Removing empty and unrecognised columns...')
    allowed_chars = IUPAC.IUPACProtein.letters + '-'
    is_empty_column = []
    contains_unk_chars = []
    for column in range(family.get_alignment_length()):
        is_empty_column.append(all([x == '-' for x in family[:, column]]))
        contains_unk_chars.append(any([x not in allowed_chars for x in family[:, column]]))
    log.info('Removed empty columns: {}'.format(','.join([str(i + 1) for i, x in enumerate(is_empty_column) if x])))
    log.info('Removed malformed columns: {}'.format(','.join([str(i + 1) for i, x in enumerate(contains_unk_chars) if x])))
    column_mask = [x | y for x, y in zip(is_empty_column, contains_unk_chars)]
    orig_col_nums = [ind + 1 for ind, x in enumerate(column_mask) if not x]
    masked_alignment = apply_column_mask(family, is_empty_column)

    # Write
    masked_fasta = family_name + '_filtered_no_empty.fa'
    log.info('Writing masked alignment to {}'.format(masked_fasta))
    with open(masked_fasta, 'wb') as output:
        AlignIO.write(masked_alignment, output, 'fasta')

    # Execute JABAWS
    log.info('Launching JABAWS...')
    jabaws_jar = args.jabaws_jar
    jabaws_server = args.jabaws_server
    AACons_presets = args.aacons_preset
    fasta_file = masked_fasta
    scores_out = family_name + '_jabaws_scores.txt'
    jabaws_cmd = ['java', '-jar', jabaws_jar,
                  '-h={}'.format(jabaws_server),
                  '-s=AAConWS',
                  '-i={}'.format(fasta_file),
                  '-r={}'.format(AACons_presets),
                  '-o={}'.format(scores_out)]
    #call(jabaws_cmd, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
    call(jabaws_cmd)

    # Read, parse and save results
    log.info('Formatting JABAWS results...')
    cons_scores = pd.read_table(scores_out, sep=' ', comment='>', index_col=0, header=None)
    cons_scores = cons_scores.transpose()
    cons_scores['orig_col_nums'] = orig_col_nums
    cons_scores_file = scores_out.replace('txt', 'csv')
    log.info('Writing formatted AACons results to {}'.format(cons_scores_file))
    cons_scores.to_csv(cons_scores_file)
    log.info('DONE.')