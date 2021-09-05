# -*- coding: utf-8 -*-
"""
Created on Sun Sep  5 16:56:22 2021

@author: Laura Xénard
"""

from pathlib import Path
import time

from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP

import protein as ptn


if __name__ == '__main__':

    start = time.time()

    # TODO: Argument à sortir par argparse
    pdb_path = "G:/RAID/Fac/M2_BI/PGP/CompBio/data/2n90.pdb"
    model = 0
    chain = 'A'

    # Ouverture et parsing du fichier PDB
    p = PDBParser()
    ptn_id = Path(pdb_path).stem
    structure = p.get_structure(ptn_id, pdb_path)
    dssp = DSSP(structure[model], pdb_path)



    # Création des résidus
    residues = []
    for i_res, res in enumerate(structure[model][chain]):

        pt = ptn.Point(*res['CA'].coord)
        tmp = ptn.Residue(res.resname, pt, dssp[(chain, i_res+1)][3])
        residues.append(tmp)





# =============================================================================
#     try:
#         hydrophobic = residues[0].is_hydrophobic()
#     except ValueError:
#         print(f"Can't determine hydrophobicity of {res}: unknown amino acid.")
#     print(hydrophobic)
# =============================================================================






    end = time.time() - start
    #print('\nDONE in {:.0f} min {:.2f} s.'.format(end // 60, end % 60))