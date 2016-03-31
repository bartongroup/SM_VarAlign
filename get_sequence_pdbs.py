"""
This module will take an alignment as input and provide PDB files stripped down to the regions covered by the sequence.
"""

import argparse
from Bio import AlignIO, SeqIO, pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.pairwise2 import format_alignment
from Bio.Alphabet import generic_protein
from Bio.PDB import *
from Bio.PDB.Dice import ChainSelector, extract
import code
import sys
sys.path.extend(['/Users/smacgowan/PycharmProjects/ProteoFAV'])
from proteofav.structures import sifts_best
from annotate_alignment import parse_seq_name
from utils import query_uniprot
import os
import re
import logging

log = logging.getLogger(__name__)


# Subclass chain selector so that it selects only C-alphas
class alpha_select(ChainSelector):
    def accept_atom(self, atom):
        name = atom.get_id()
        if name == 'CA':
            return 1
        else:
            return 0


def main(args):
    downloads = '.SequencePdbs'
    output_dir = 'sequence_pdbs'

    alignment = AlignIO.read(args.fasta_file, "fasta")

    output_pdb_files = []
    for seq in alignment:
        # Get PDBs through SIFTS
        seq_name = parse_seq_name(seq.id)
        if args.seq_name_ids == 'UniProt':
            pass
        elif args.seq_name_ids == 'UniProt_gene':
            query = re.search('[^_\W]+', seq_name).group().strip()
            seq_name = query_uniprot(('gene:' + query, 'reviewed:yes', 'organism:human'), first=True)
        elif args.seq_name_ids == 'mnemonic':
            query = seq_name
            seq_name = query_uniprot(('mnemonic:' + query, 'reviewed:yes', 'organism:human'), first=True)
        else:
            log.error('Sequence names must correspond to either UniProt IDs or gene names')
            raise NotImplemented
        sifts_pdb = sifts_best(seq_name)
        if sifts_pdb is None:
            log.info('No SIFTs mapping for {}.'.format(seq_name))
            continue

        # Pick out first X-ray structure
        # TODO: better structure selection
        pdb_id = None
        for prot, record in sifts_pdb.iteritems():
            for r_dict in record:
                if r_dict.get('experimental_method') == 'X-ray diffraction':
                    pdb_id = r_dict.get('pdb_id')
                    break
        if pdb_id is None:
            log.info('No X-ray structure for {}.'.format(seq_name))
            continue

        # Download and parse that PDB
        pdbl = PDBList()
        pdb_dir = os.path.join(downloads, 'PDBs')
        pdbl.retrieve_pdb_file(pdb_id, pdir=pdb_dir)
        parser = PDBParser()
        pdb_file_path = os.path.join(pdb_dir, 'pdb' + pdb_id + '.ent')
        structure = parser.get_structure(pdb_id, pdb_file_path)

        # Now need to identify PDB residue numbers to keep, can do this either through Biopython PPBuilder or
        # ProteoFAV structure tables

        # First build PDB chains
        ppb = PPBuilder()
        model = structure[0]
        pdb_chain_sequences = []
        for chain in model.get_list():
            chain_id = chain.get_id()
            chain_pp = Seq('', generic_protein)
            for pp in ppb.build_peptides(chain):
                chain_pp = chain_pp + pp.get_sequence()  # Concatenate all PPs from a single chain
            pdb_chain_sequences.append((structure.get_id(), chain_id, chain_pp, ppb.build_peptides(chain)))

        # Now compare to sequence and identify PDB selection
        seq_str = str(seq.seq).replace('-', '').upper()
        seq_matched = False
        for pdb, c, pdb_seq, pp_list in pdb_chain_sequences:
            # TODO: Pick best rather than first match?
            if not seq_matched and str(seq_str) in str(pdb_seq):
                # Find start
                length = len(str(seq_str))
                seq_start = str(pdb_seq).find(str(seq_str))
                seq_end = seq_start + length
                # Identify the pp that has the first residue and get the resseq
                pp_coverage = 0
                for pp in pp_list:
                    pp_coverage += len(pp)
                    if seq_start < pp_coverage:
                        start = pp[seq_start + pp_coverage - len(pp)].get_id()[1]
                        break
                # Identify the pp with the last residue
                pp_coverage = 0
                for pp in pp_list:
                    pp_coverage += len(pp)
                    if seq_end < pp_coverage:
                        end = pp[seq_end + pp_coverage - len(pp)].get_id()[1]
                        break

                select_residues = (pdb, c, start, end)
                seq_matched = True
        if not seq_matched:
            log.info('Could not match {} to {}.'.format(seq_name, pdb_id))
            continue

        pdb_id = select_residues[0]
        #if phd_resnums[1] > 0:  # Having some problems just now...
        #    continue
        chain_id = select_residues[1]  # Assuming this is the chain we got the sequence for earlier
        start = select_residues[2]
        end = select_residues[3]
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        filename = os.path.join(output_dir, seq_name + '_' + structure.get_id() + '.pdb')
        #extract(structure, chain_id, start, end, filename)
        sel = alpha_select(chain_id, start, end)
        io = PDBIO()
        io.set_structure(structure)
        io.save(filename, sel)

        output_pdb_files.append(filename)
        log.info('Successfully matched {} to {}.'.format(seq_name, pdb_id))

    # Lastly create MAMMOTH input file
    with open('MAMMOTH_input.pdb', 'w') as mammoth_file:
        mammoth_header = ('<pre>\n'
                          'REMARK MAMMOTH-mult format rules:\n'
                          'REMARK   (1) read only "CA" entries\n'
                          'REMARK   (2) duplicate residues are ignored (same residue number\n'
                          'REMARK       consecutevely)\n'
                          'REMARK   (3) each structure must be separated by "TER"\n'
                          'REMARK   (4) 20 residues or more.\n'
                          'REMARK')
        mammoth_file.write(mammoth_header)
    with open('MAMMOTH_input.pdb', 'a') as mammoth_file:
        for pdb_file_path in output_pdb_files:
            with open(pdb_file_path, 'r') as pdb_file:
                for line in pdb_file.readlines():
                    if 'END' not in line:
                        mammoth_file.write(line)
    with open('MAMMOTH_input.pdb', 'a') as mammoth_file:
        mammoth_file.write('</pre>\n')



if __name__ == '__main__':
    # CLI Arguments
    parser = argparse.ArgumentParser(description='Fetch and subset PDB files that correspond to sequences in an alignment.')
    parser.add_argument('fasta_file', type=str, help='Path to fasta file containing MSA.')
    parser.add_argument('--seq_name_ids', type=str, default='mnemonic',
                        help='Definitions the alignment sequences names come from ["UniProt", "mnemonic", "UniProt_gene"]')
    parser.add_argument('--interpreter', action='store_true',
                        help='Drop into interactive python session once analysis is complete.')
    args = parser.parse_args()

    main(args)

    if args.interpreter:
        code.interact(local=dict(globals(), **locals()))