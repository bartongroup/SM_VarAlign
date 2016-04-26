import argparse
from Bio import SeqIO

def batch_iterator(iterator, batch_size) :
    """Returns lists of length batch_size.
 
    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.
 
    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True #Make sure we loop once
    while entry :
        batch = []
        while len(batch) < batch_size :
            try :
                entry = iterator.next()
            except StopIteration :
                entry = None
            if entry is None :
                #End of file
                break
            batch.append(entry)
        if batch :
            yield batch


if __name__ == '__main__':
	# CLI
	parser = argparse.ArgumentParser(description='Script to split a fasta file.')
	parser.add_argument('path_to_fasta', type=str, help='Path to fasta file.')
	args = parser.parse_args()

	# Write subsets of the FASTA file
	record_iter = SeqIO.parse(open(args.path_to_fasta),"fasta")
	for i, batch in enumerate(batch_iterator(record_iter, 1000)) :
	    filename = "group_%i.fasta" % (i+1)
	    handle = open(filename, "w")
	    count = SeqIO.write(batch, handle, "fasta")
	    handle.close()
	    print "Wrote %i records to %s" % (count, filename)