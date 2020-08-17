#!/usr/bin/env python3
from sys import argv
from Bio.PDB import PDBList

for name in argv[1:]:
    PDBList().retrieve_pdb_file(name, file_format='pdb')
