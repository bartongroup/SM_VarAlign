import annotate_alignment
import argparse
from Bio import AlignIO, SeqIO  # Needs PR #768 #769 patched
import logging
from config import defaults
from utils import filter_alignment
import os
import pandas as pd


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
    args = parser.parse_args()

    # Arguments and parameters
    local_uniprot_path = args.local_uniprot_file
    local_pfam = args.local_pfam
    seq_id_filter = args.seq_id_filter
    use_local_alignment = False
    downloads = args.downloads

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

    # Parse Pfam alignment
    log.info('Parsing PFAM: {}'.format(local_pfam))
    pfam = AlignIO.parse(local_pfam, 'stockholm')

    # Identify array number
    job_number = int(os.environ.get('SGE_TASK_ID'))
    log.info('Array job number: {}'.format(job_number))

    # Iterate correct number of times through Pfam, and exit if out-of-range
    log.info('Iterating through PFAM generator.')
    for i in xrange(job_number):
        try:
            alignment = pfam.next()
        except StopIteration:
            log.error('PFAM index out-of-range: less than {} alignments in PFAM file.'.format(job_number))
            raise SystemExit

    log.debug('Reading PFAM alignment ID.')
    alignment_name = alignment.annotations['GF']['AC'][0]  # This requires biopython patches #768 #769
    alignment_name += '_' + local_pfam
    log.info('Found alignment: {}'.format(alignment_name))

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