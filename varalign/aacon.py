import argparse
import logging
import os
import subprocess

import pandas as pd
from Bio import AlignIO  # Needs PR #768 #769 patched
from Bio.Alphabet import IUPAC

import varalign
from varalign.jabaws import apply_column_mask
from varalign.utils import sanitise_alignment, make_dir_if_needed

log = logging.getLogger(__name__)
log.setLevel('INFO')

aacon_path = os.path.join(os.path.dirname(varalign.__file__), 'lib', 'aacon', 'compbio-conservation-1.1.jar')


def _reformat_alignment_for_aacon(aln):
    """
    Reformat an alignment so that it is compatible with AACon.

    :param aln: Alignment
    :return:
    """
    # Perform general sanitisation
    aacon_alignment = sanitise_alignment(aln)

    # Identify empty columns and those with unknown characters
    allowed_chars = IUPAC.IUPACProtein.letters + '-'
    is_empty_column = []
    contains_unk_chars = []
    for column in range(aacon_alignment.get_alignment_length()):
        is_empty_column.append(all([x == '-' for x in aacon_alignment[:, column]]))
        contains_unk_chars.append(any([x not in allowed_chars for x in aacon_alignment[:, column]]))

    # Remove those columns
    column_mask = [x | y for x, y in zip(is_empty_column, contains_unk_chars)]
    orig_col_nums = [ind + 1 for ind, x in enumerate(column_mask) if not x]
    aacon_alignment = apply_column_mask(aacon_alignment, is_empty_column)

    # Log modifications
    if any(is_empty_column):
        log.info('Removed empty columns: {}'.format(','.join([str(i + 1) for i, x in enumerate(is_empty_column) if x])))
    if any(contains_unk_chars):
        log.info('Removed malformed columns: {}'.format(','.join([str(i + 1) for i, x in enumerate(contains_unk_chars) if x])))

    return aacon_alignment, orig_col_nums


def _run_aacon(aln, column_index=None, aacon_jar_path=aacon_path):
    """
    Run standalone AACon on an alignment.

    :param aln: Alignment.
    :param column_index: Additional column numbering to add to results table.
    :param aacon_jar_path: Path to AACon jar.
    :return:
    """
    tmp_dir = os.path.join('.varalign', 'aacon')
    make_dir_if_needed(tmp_dir)

    aacon_input = os.path.join(tmp_dir, 'aacon_input.fa')
    aacon_output = os.path.join(tmp_dir, 'aacon_scores.out')

    # Write alignment to disk as fasta
    with open(aacon_input, 'w') as output:
        AlignIO.write(aln, output, 'fasta')
    log.info('AACon formatted alignment saved to {}'.format(aacon_input))

    # Run AACon on alignment
    log.info('Launcing AACon...')
    with open(os.path.join(tmp_dir, 'aacon.log'), 'wb') as aacon_log:
        subprocess.call(["java", "-jar", aacon_jar_path,
                         "-i={}".format(aacon_input),
                         "-o={}".format(aacon_output)],
                        stdout=aacon_log, stderr=aacon_log)

    # Read results and format
    aacon_table = pd.read_table(aacon_output, sep=' ', comment='>', index_col=0, header=None)
    aacon_table = aacon_table.transpose().dropna()  # Needed to add `dropna` as extra row...
    aacon_table.columns = aacon_table.columns.str.replace('#', '')
    aacon_table.columns = aacon_table.columns.str.lower()

    if column_index:
        aacon_table['column_index'] = column_index
        aacon_table.set_index('column_index', inplace=True)

    return aacon_table


if __name__ == '__main__':
    # CLI
    parser = argparse.ArgumentParser(description='Run AACon a Pfam alignment.')
    parser.add_argument('alignment', type=str, help='Path to the alignment.')
    args = parser.parse_args()

    alignment = AlignIO.read(args.alignment, format='stockholm')

    # Format for AACon and run
    aacon_alignment, orig_col_nums = _reformat_alignment_for_aacon(alignment)
    alignment_conservation = _run_aacon(aacon_alignment, orig_col_nums)

    # Save result
    cons_scores_file = 'aacon_scores.csv'
    alignment_conservation.to_csv(cons_scores_file)
    log.info('Formatted AACons results saved to {}'.format(cons_scores_file))

