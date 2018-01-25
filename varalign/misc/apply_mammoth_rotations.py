import argparse
from Bio.PDB import *
import numpy as np

if __name__ == '__main__':
    # CLI
    parser = argparse.ArgumentParser(description='Apply a mammoth rotation to a PDB file.')
    parser.add_argument('mammoth', type=str, help='Path to mammoth rotation file.')
    parser.add_argument('number', type=int, help='The rotation to use from the mammoth rotation file.')
    parser.add_argument('pdb', type=str, help='Path to pdb file')
    args = parser.parse_args()

    # read mammoth rotation and translation
    name_strings = []
    transform_strings = []
    with open(args.mammoth, 'r') as mammoth_file:
        lines = mammoth_file.readlines()
        for line in lines:
            if len(line) == 126:
                transform_strings.append(line.strip())
            else:
                name_strings.append(line.strip())

    # parse desired transformation
    # NB: first line is header
    rot = transform_strings[args.number].strip().split()[1:10]
    cn = transform_strings[args.number].strip().split()[10:13] # centre of mass
    tr = transform_strings[args.number].strip().split()[13:]

    # Convert to arrays
    rot = np.array([rot[:3], rot[3:6], rot[6:]], 'f')
    rot = np.transpose(rot)
    cn = np.array(cn, 'f')
    tr = np.array(tr, 'f')
    identity = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], 'f')

    print('Applying:\n{}\n{}'.format(name_strings[args.number - 1],
                                     transform_strings[args.number]))
    print('Rotation:\n{}'.format(rot))
    print('Centre:\n{}'.format(cn))
    print('Translation:\n{}'.format(tr))

    # Apply to all atoms in PDB file
    pdb_parser = PDBParser()
    structure = pdb_parser.get_structure(args.pdb.replace('.pdb', ''), args.pdb)
    for model in structure.get_iterator():
        for chain in model.get_iterator():
            for residue in chain.get_iterator():
                for atom in residue.get_iterator():
                    atom.transform(identity, -1 * cn)
                    atom.transform(rot, tr)

    # Write the new file
    io = PDBIO()
    io.set_structure(structure)
    io.save(args.pdb.replace('.pdb', '_rotated.pdb'))
