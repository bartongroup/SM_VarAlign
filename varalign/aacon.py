import argparse
from Bio import AlignIO  # Needs PR #768 #769 patched
from Bio.Alphabet import IUPAC
from jabaws import apply_column_mask
import logging
import pandas as pd
import pfam
from utils import filter_alignment, sanitise_alignment
import subprocess


log = logging.getLogger(__name__)
log.setLevel('INFO')


def _reformat_alignment_for_aacon(aln):
    """
    Reformat an alignment so that it is compatible with AACon.

    :param aln: Alignment
    :return:
    """
    aacon_alignment = sanitise_alignment(aln)
    log.info('Removing empty and unrecognised columns...')

    # Identify empty columns and those with unknown characters
    allowed_chars = IUPAC.IUPACProtein.letters + '-'
    is_empty_column = []
    contains_unk_chars = []
    for column in range(aacon_alignment.get_alignment_length()):
        is_empty_column.append(all([x == '-' for x in aacon_alignment[:, column]]))
        contains_unk_chars.append(any([x not in allowed_chars for x in aacon_alignment[:, column]]))
    log.info('Empty columns: {}'.format(','.join([str(i + 1) for i, x in enumerate(is_empty_column) if x])))
    log.info('Malformed columns: {}'.format(','.join([str(i + 1) for i, x in enumerate(contains_unk_chars) if x])))

    # Remove those columns
    column_mask = [x | y for x, y in zip(is_empty_column, contains_unk_chars)]
    orig_col_nums = [ind + 1 for ind, x in enumerate(column_mask) if not x]
    aacon_alignment = apply_column_mask(aacon_alignment, is_empty_column)
    log.info('Empty and malformed columns removed.')

    return aacon_alignment, orig_col_nums


def _run_aacon(aln, aacon_jar_path='/homes/smacgowan/bin/compbio-conservation-1.1.jar'):
    """
    Run standalone AACon on an alignment.

    :param aln: Alignment.
    :param aacon_jar_path: Path to AACon jar.
    :return:
    """
    aacon_input = 'aacon_input.fa'
    aacon_output = 'aacon_scores.out'

    # Write alignment to disk as fasta
    log.info('Writing alignment to {}...'.format(aacon_input))
    with open(aacon_input, 'wb') as output:
        AlignIO.write(aln, output, 'fasta')

    # Run AACon on alignment
    log.info('Launcing AACon...')
    with open('aacon.log', 'wb') as aacon_log:
        subprocess.call(["java", "-jar", aacon_jar_path,
                         "-i={}".format(aacon_input),
                         "-o={}".format(aacon_output)],
                        stdout=aacon_log, stderr=aacon_log)


if __name__ == '__main__':
    # CLI
    parser = argparse.ArgumentParser(description='Run AACon a Pfam alignment.')
    parser.add_argument('alignment', type=str, help='Path to the alignment.')
    args = parser.parse_args()

    alignment = AlignIO.read(args.alignment, format='stockholm')

    # Format for AACon and run
    aacon_alignment, orig_col_nums = _reformat_alignment_for_aacon(alignment)
    _run_aacon(aacon_alignment)

    # Index to orginal column numbers and reformat
    log.info('Formatting JABAWS results...')
    cons_scores_table = pd.read_table('aacon_scores.out', sep=' ', comment='>', index_col=0, header=None)
    cons_scores_table = cons_scores_table.transpose().dropna()  # Needed to add `dropna` as extra row...
    cons_scores_table['orig_col_nums'] = orig_col_nums

    # Save result
    cons_scores_file = 'aacon_scores.out'.replace('out', 'csv')
    log.info('Writing formatted AACons results to {}'.format(cons_scores_file))
    cons_scores_table.to_csv(cons_scores_file)
