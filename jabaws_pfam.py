import argparse
from Bio import AlignIO, SeqIO  # Needs PR #768 #769 patched
import logging
from config import defaults
from utils import filter_alignment
from itertools import groupby, tee, izip
import os
from subprocess import call
import pandas as pd

log = logging.getLogger(__name__)
logging.captureWarnings(True)
logging.basicConfig(level=9,
                    format='%(asctime)s - %(levelname)s - %(message)s ')

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

    # Parse Pfam family
    log.info('Parsing PFAM: {}'.format(local_pfam))
    pfam = AlignIO.parse(local_pfam, 'stockholm')

    # Identify array number and lookup Pfam list for corresponding AC
    job_number = int(os.environ.get('SGE_TASK_ID'))
    log.info('Array job number: {}'.format(job_number))
    with open(family_id_list, 'r') as ids:
        for i in xrange(job_number):
            desired_family = ids.readline().strip()
    log.info('Searching Pfam for family: {}'.format(desired_family))

    # Iterate through Pfam until found desired family, and exit if out-of-range
    log.info('Iterating through Pfam generator...')
    for family in pfam:
        try:
            family = pfam.next()
            family_name = family.annotations['GF']['AC'][0]  # This requires biopython patches #768 #769
            log.debug('Read family {}...'.format(family_name))
            if family_name.startswith(desired_family):
                break
        except StopIteration:
            log.error('PFAM index out-of-range: could not find {} in {}.'.format(desired_family, local_pfam))
            raise SystemExit

    # Store full family name
    log.info('Found family: {}'.format(family_name))
    family_name += '_' + os.path.basename(local_pfam)

    # Filter unwanted sequences
    log.info('Filtering family for sequences without "{}"'.format(seq_id_filter))
    if seq_id_filter:
        family = filter_alignment(family, seq_id_filter)
        if len(family) == 0:
            log.warning('No sequences passed filter. Exiting.')
            raise SystemExit

    # Remove empty columns from family and write
    # Identify empty columns
    log.info('Removing empty columns...')
    is_empty_column = []
    for column in range(family.get_alignment_length()):
        is_empty_column.append(all([x == '-' for x in family[:, column]]))
    orig_col_nums = [ind + 1 for ind, x in enumerate(is_empty_column) if not x]
    masked_alignment = apply_column_mask(family, is_empty_column)
    masked_fasta = family_name + '_filtered_rmc.fa'
    log.info('Writing masked alignment to {}'.format(masked_fasta))
    with open(masked_fasta, 'wb') as output:
        AlignIO.write(masked_alignment, output, 'fasta')

    # Execute JABAWS
    log.info('Launching JABAWS...')
    jabaws_jar = args.jabaws_jar
    jabaws_server = args.jabaws_server
    AACons_presets = args.aacons_preset
    fasta_file = masked_fasta
    scores_out = 'alignment_name_jabaws_scores.txt'
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