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

def remove_empty_columns(aln):
    """
    Remove empty columns from an alignment.

    :param aln: Alignment.
    :return: Alignment with empty columns removed and RLE for blocks. (aln, empty_run_lengths)
    """
    def reverse_degapping(empty_run_lengths):
        """

        :param empty_run_lengths:
        :param degapped_column_properties:
        :return:
        """
        i = 0
        rebuilt = []
        for length, insert_gap in empty_run_lengths:
            run = 0
            while run < length:
                if insert_gap:
                    i += 1
                else:
                    rebuilt.append(i)
                    i += 1
                run += 1
        return rebuilt


    def pairwise(iterable):
        "s -> (s0,s1), (s1,s2), (s2, s3), ..."
        a, b = tee(iterable)
        next(b, None)
        return izip(a, b)

    # Identify empty columns
    is_empty_column = []
    for column in range(aln.get_alignment_length()):
        is_empty_column.append(all([x == '-' for x in aln[:, column]]))

    # Provide original column numbers
    empty_run_lengths = [(len(list(g)), k)
                         for k, g in groupby(is_empty_column)]
    orig_col_nums = reverse_degapping(empty_run_lengths)

    # Identify boundaries between empty / non-empty
    is_boundary = []
    for i, j in zip(is_empty_column[:-1], is_empty_column[1:]):
        is_boundary.append(i != j)

    # These indexes will be compiled into sequential, pairwise tuples for slicing
    boundary_indexes = [None]  # list[None:i]
    for n, k in enumerate(is_boundary):
        if k:
            boundary_indexes.append(n + 1) # +1 for slicing
    boundary_indexes.append(None) # list[i:None]

    # Boundaies as range blocks
    ranges = pairwise(boundary_indexes)

    # Initialise new alignment
    keep_range = not is_empty_column[0]
    start, stop = ranges.next()
    if keep_range:
        degapped = aln[:, start:stop]
    else:
        start, stop = ranges.next()
        keep_range = not keep_range
        degapped = aln[:, start:stop]
    keep_range = not keep_range

    # Proceed through ranges
    for start, stop in ranges:
        if keep_range:
            degapped = degapped + aln[:, start:stop]
        keep_range = not keep_range

    return degapped, orig_col_nums


if __name__ == '__main__':
    # CLI Arguments
    parser = argparse.ArgumentParser(description='Create variant counts column wise for an alignment.')
    parser.add_argument('--local_pfam', type=str, help='Path to PFAM-A.seed')
    parser.add_argument('--seq_id_filter', type=str,
                        help='An inclusive filter to process only a subset of sequences.')
    parser.add_argument('--families', type=str, help='Path to list of families.')
    parser.add_argument('--downloads', type=str, default=defaults.db_root,
                        help='A directory to store downloaded files.')
    args = parser.parse_args()

    # Arguments and parameters
    local_pfam = args.local_pfam
    seq_id_filter = args.seq_id_filter
    families_list = args.families
    downloads = args.downloads

    # Parse Pfam alignment
    log.info('Parsing PFAM: {}'.format(local_pfam))
    pfam = AlignIO.parse(local_pfam, 'stockholm')

    # Identify array number and corresponding Pfam
    job_number = int(os.environ.get('SGE_TASK_ID'))
    log.info('Array job number: {}'.format(job_number))
    with open(families_list, 'r') as family_acs:
        for i in xrange(job_number):
            desired_pfam_ac = family_acs.readline().strip()
    log.info('Searching Pfam for family: {}'.format(desired_pfam_ac))

    # Iterate through Pfam until found desired family, and exit if out-of-range
    log.info('Iterating through Pfam generator...')
    for alignment in pfam:
        try:
            alignment = pfam.next()
            alignment_name = alignment.annotations['GF']['AC'][0]
            log.debug('Read family {}...'.format(alignment_name))
            if alignment_name.startswith(desired_pfam_ac):
                break
        except StopIteration:
            log.error('PFAM index out-of-range: could not find {} in {}.'.format(desired_pfam_ac, local_pfam))
            raise SystemExit

    # Store full alignment name
    log.debug('Reading PFAM alignment ID.')
    alignment_name = alignment.annotations['GF']['AC'][0]  # This requires biopython patches #768 #769
    alignment_name += '_' + os.path.basename(local_pfam)
    log.info('Found alignment: {}'.format(alignment_name))

    # Filter unwanted sequences
    log.info('Filtering alignment for sequences without "{}"'.format(seq_id_filter))
    if seq_id_filter:
        alignment = filter_alignment(alignment, seq_id_filter)
        if len(alignment) == 0:
            log.warning('No sequences passed filter. Exiting.')
            raise SystemExit

    # Remove empty columns from alignment and write
    processed_alignment, orig_col_nums = remove_empty_columns(alignment)
    degapped_fasta = alignment_name + '_filtered_rmc.fa'
    with open(degapped_fasta, 'wb') as file_:
        AlignIO.write(processed_alignment, file_, 'fasta')

    # Run analysis: 1) Write alignment, 2) Execute JABAWS
    log.info('Processing alignment...')
    # Construct JABAWS command and execute
    jabaws_jar = '/Applications/min-jabaws-client-2.1.0.jar'
    jabaws_server = 'http://www.compbio.dundee.ac.uk/jabaws'
    AACons_presets = 'Quick conservation'
    fasta_file = degapped_fasta
    scores_out = 'alignment_name_jabaws_scores.txt'
    jabaws_cmd = ['java', '-jar', jabaws_jar,
                  '-h={}'.format(jabaws_server),
                  '-s=AAConWS',
                  '-i={}'.format(fasta_file),
                  '-r={}'.format(AACons_presets),
                  '-o={}'.format(scores_out)]
    #call(jabaws_cmd, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
    call(jabaws_cmd)

    # Write results: either read from stdout or set JABAWS argument output file
    log.info('Writing results...')
    # Read, parse and save results
    cons_scores = pd.read_table(scores_out,
                                sep=' ', comment='>', index_col=0, header=None)
    cons_scores = cons_scores.transpose()
    cons_scores['orig_col_nums'] = orig_col_nums
    cons_scores.to_csv(scores_out.replace('txt', 'csv'))
    log.info('DONE.')