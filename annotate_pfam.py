import annotate_alignment
import argparse
from Bio import AlignIO, SeqIO  # Needs PR #768 #769 patched
import logging
from config import defaults
from utils import filter_alignment
import os
import pandas as pd
import pfam_index


log = logging.getLogger(__name__)


if __name__ == '__main__':
    # CLI Arguments
    parser = argparse.ArgumentParser(description='Create variant counts column wise for an alignment.')
    parser.add_argument('--local_uniprot_file', type=str, help='Path to uniprot_sprot.dat')
    parser.add_argument('--local_pfam', type=str, help='Path to PFAM-A.seed')
    parser.add_argument('--seq_id_filter', type=str,
                        help='An inclusive filter to process only a subset of sequences.')
    parser.add_argument('--downloads', type=str, default=defaults.db_root,
                        help='A directory to store downloaded files.')
    parser.add_argument('--families', type=str, help='Path to list of families.')
    args = parser.parse_args()

    # Arguments and parameters
    local_uniprot_path = args.local_uniprot_file
    local_pfam = args.local_pfam
    seq_id_filter = args.seq_id_filter
    use_local_alignment = False
    downloads = args.downloads
    family_id_list = args.families

    # Initialise local UniProt if provided
    local_uniprot_file = args.local_uniprot_file
    if local_uniprot_file is not None:
        if args.format != 'stockholm':
            log.error('Can only use local UniProt with Stockholm alignments that have AC annotations.')
            raise TypeError
        log.info('Indexing local UniProt database: {}'.format(local_uniprot_path))
        local_uniprot_index = SeqIO.index(local_uniprot_path, 'swiss')
    else:
        local_uniprot_index = None

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

    start = pfam_index.lookup_index(index_path, desired_family)
    if start is None:
        log.error('Pfam {} could not be found in the index'.format(desired_family))
        raise SystemExit
    alignment = pfam_index.read_family(local_pfam, start)
    alignment_name = alignment.annotations['GF']['AC'][0]  # This requires biopython patches #768 #769
    log.debug('Read family {}...'.format(alignment_name))

    # Store full family name
    log.info('Found family: {}'.format(alignment_name))
    alignment_name += '_' + os.path.basename(local_pfam)

    # Filter unwanted sequences
    log.info('Filtering alignment for sequences without "{}"'.format(seq_id_filter))
    if seq_id_filter:
        alignment = filter_alignment(alignment, seq_id_filter)
        if len(alignment) == 0:
            log.warning('No sequences passed filter. Exiting.')
            raise SystemExit

    # Run analysis
    log.info('Processing alignment...')
    merged_table, fisher_results, rvis_scores = annotate_alignment.main(alignment, alignment_name, use_local_alignment,
                                                                        local_uniprot_index, downloads)

    # Write results
    log.info('Writing results...')
    merged_table.to_csv(alignment_name + '_merged_table.csv')
    scores = pd.DataFrame(fisher_results)
    scores.columns = ['or', 'p']
    scores['rvis'] = pd.Series(rvis_scores)
    scores.to_csv(alignment_name + '_scores.csv')
    log.info('DONE.')
