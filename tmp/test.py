# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP

pdb_path = "/home/sdv/m2bi/lxenard/Documents/PGP/CompBio/data/2n90.pdb"

# Ouverture et parsig du fichier PDB
p = PDBParser()
structure = p.get_structure("2N90", pdb_path)
model = structure[0]
dssp = DSSP(model, pdb_path)

# Liste des propriétés accessibles par un tuple (chaine, n° résidu)
print(dssp[('A', 1)])
# ASA
print(dssp[('A', 1)][3])

print(dssp[('A', 1)])
print(dssp[('B', 101)])
