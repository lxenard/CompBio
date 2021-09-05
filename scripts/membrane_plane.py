# -*- coding: utf-8 -*-
"""
Created on Sun Sep  5 16:56:22 2021

@author: Laura Xénard
"""

import argparse
from pathlib import Path
import time

from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP

import protein as ptn


def barycenter(residues_list):
    """
    Compute the barycenter of a list of residues.

    Parameters
    ----------
    residues_list : list(ptn.Residues)
        The residues from which to find the barycenter.

    Returns
    -------
    x_bary, y_bary, z_bary : (float, float, float)
        Respectively x, y and z coordinates of the residues barycenter.

    """
    x_bary = sum([r.coord.x for r in residues_list]) / len(residues_list)
    y_bary = sum([r.coord.y for r in residues_list]) / len(residues_list)
    z_bary = sum([r.coord.z for r in residues_list]) / len(residues_list)
    return x_bary, y_bary, z_bary


if __name__ == '__main__':

    start = time.time()

    # TODO: Argument à sortir par argparse
    pdb_path = 'G:/RAID/Fac/M2_BI/PGP/CompBio/data/2n90.pdb'
    model = 0
    chain = 'A'
    IS_EXPOSED_THRESHOLD = 0.3
    DEBUG = False
    N_DIRECTIONS = 8 # nb de points pour échantillonner la demi-sphère

    # Opening and parsing of the PDB file.
    p = PDBParser()
    ptn_id = Path(pdb_path).stem
    structure = p.get_structure(ptn_id, pdb_path)
    dssp = DSSP(structure[model], pdb_path)

    # Selection of the exposed residues.
    # For better results, the burrowed residues are not taken into account
    # during the membrane detection.
    # TODO: à vérifier en comparant les résultats avec le traitement de tous
    # les résidus
    residues = []  # Store the exposed residues.
    for i_res, res in enumerate(structure[model][chain]):
        # For simplification, the position of a residue is defined as the
        # position of its Cα.
        pt = ptn.Point(*res['CA'].coord)
        asa = dssp[(chain, i_res+1)][3]  # Accessible surface area.
        tmp = ptn.Residue(res.id[1], res.resname, pt, asa)
        if tmp.is_exposed(IS_EXPOSED_THRESHOLD):
            residues.append(tmp)

    # Place the center of the coordinate system at the residues barycenter.
    bary = ptn.Point(*barycenter(residues))
    if DEBUG:
        print(f"Barycenter: {bary}")
    for res in residues:
        res.coord -= bary

    # Generate N equidistributed points on the surface of a sphere













# =============================================================================
#     try:
#         hydrophobic = residues[0].is_hydrophobic()
#     except ValueError:
#         print(f"Can't determine hydrophobicity of {res}: unknown amino acid.")
#     print(hydrophobic)
# =============================================================================






    end = time.time() - start
    print('\nDONE in {:.0f} min {:.2f} s.'.format(end // 60, end % 60))