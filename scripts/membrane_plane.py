# -*- coding: utf-8 -*-
"""
Created on Sun Sep  5 16:56:22 2021

@author: Laura Xénard
"""


from pathlib import Path
import time

from Bio.PDB import PDBParser, Dice
from Bio.PDB.Atom import Atom
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.Residue import Residue
import numpy as np

import protein as ptn
import settings as st

import warnings


if __name__ == '__main__':

    start = time.time()

    st.init()  # Set the global parameters.

    # Opening and parsing of the PDB file.
    p = PDBParser()
    ptn_id = Path(st.PDB).stem
    structure = p.get_structure(ptn_id, st.PDB)
    prot = ptn.Protein(structure, st.MODEL, st.CHAIN, st.FIRST_RESIDUE,
                       st.LAST_RESIDUE)

    # For better results, the burrowed residues are not taken into account
    # during the membrane detection, only the exposed ones.
    # TODO: à vérifier en comparant les résultats avec le traitement de tous
    # les résidus
    # Place the center of the coordinate system at the residues barycenter.
    bary = ptn.Point.barycenter(prot.residues_exposed)
    if st.DEBUG:
        print(f"Barycenter: {bary}")
    prot.move(-bary)

    centered_slices = []
    # obtenir les plans orthogonaux centrés sur l'origine
    for v in prot.vectors:
        s = ptn.Slice(prot, 0, v, st.SCORE_METHOD)
        centered_slices.append(s)

    # Translation of the centered slices along their normal vector to get
    # all possible slices.
    shift = 1
    slices = centered_slices.copy()
    for sli in centered_slices:

        # Plans translated toward the end of the normal vector ('up').
        center = shift
        new_sli_up = ptn.Slice(prot, center, sli.normal, st.SCORE_METHOD)
        while len(new_sli_up.residues) != 0:
            slices.append(new_sli_up)
            center += shift
            new_sli_up = ptn.Slice(prot, center, sli.normal, st.SCORE_METHOD)

        # Plans translated toward the start of the normal vector ('down').
        center = -shift
        new_sli_down = ptn.Slice(prot, center, sli.normal, st.SCORE_METHOD)
        while len(new_sli_down.residues) != 0:
            slices.append(new_sli_down)
            center -= shift
            new_sli_down = ptn.Slice(prot, center, sli.normal, st.SCORE_METHOD)

    slices.sort(reverse=True)
    best = next(((i, sli) for i, sli in enumerate(slices) if sli.score < 1), None)
    assert best is not None, "All slices have a score of 1."
    print(best)
    best_index, best_sli = best
    for res in best_sli.residues:
        print(res.num, res.aa)

    # Thickening the best found slice to maximise the score.
    base_score = best_sli.score
    increment = 1

    # Toward the end of the normal vector ('up').
    new_scores_up = [base_score]
    best_sli.thicken(increment, normal_direction=True)
    cpt = 0
    while best_sli.score != 0:
        new_scores_up.append(best_sli.score)
        best_sli.thicken(increment, normal_direction=True)
        cpt += 1
        if cpt >= 5 and new_scores_up[-5:].count(new_scores_up[-1]) == 5:
            break

    # Resetting the slice
    best_sli.thicken(-increment * cpt, normal_direction=False)
    print(best_sli)

    # Toward the start of the normal vector ('down').
    new_scores_down = [base_score]
    best_sli.thicken(increment, normal_direction=False)
    cpt = 0
    while best_sli.score != 0:
        new_scores_down.append(best_sli.score)
        best_sli.thicken(increment, normal_direction=False)
        cpt += 1
        if cpt >= 5 and new_scores_down[-5:].count(new_scores_down[-1]) == 5:
            break

    # Resetting the slice
    best_sli.thicken(-increment * cpt, normal_direction=True)
    print(best_sli)

    # Adding 2 dummy residues to represent the membrane delimiting planes.
    last_id = prot.structure[prot.model][prot.chain][prot.res_ids_pdb[-1]].id
    mem1_id = (last_id[0], last_id[1]+1, last_id[2])
    mem2_id = (last_id[0], last_id[1]+2, last_id[2])
    mem1 = Residue(mem1_id, 'MEM', '')
    mem2 = Residue(mem2_id, 'MEM', '')
    prot.structure[prot.model][prot.chain].add(mem1)
    prot.structure[prot.model][prot.chain].add(mem2)

    # Creating grids to place dummy atoms in order to represent the
    # membrane 2 delimiting planes.
    resolution = 1
    shift = best_sli.center
    thickness = best_sli.thickness
    a = best_sli.normal.end.x
    b = best_sli.normal.end.y
    c = best_sli.normal.end.z
    # Finding bounding coordinates of the protein.
    x_min, x_max, y_min, y_max, z_min, z_max = prot.find_bounding_coord()
    # Building the grids.
    xx, yy = np.mgrid[x_min:x_max:resolution, y_min:y_max:resolution]
    zz1 = (-a*xx - b*yy - thickness[0] + shift) / c
    zz2 = (-a*xx - b*yy + thickness[1] + shift) / c
    # Translation of the membrane planes to account for centering the 3D
    # space onto the protein barycenter.
    xx = xx + bary.x
    yy = yy + bary.y
    zz1 = zz1 + bary.z
    zz2 = zz2 + bary.z

    # Adding the dummy atoms to the 'membrane' residues.
    cpt = 1
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', PDBConstructionWarning)
        for x, y, z in zip(xx.flatten(), yy.flatten(), zz1.flatten()):
            new = Atom(f'D{cpt}', np.array([x, y, z]), 0, 1, 32,
                       f' D{cpt} ', cpt)
            mem1.add(new)
            cpt += 1
    cpt = 1
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', PDBConstructionWarning)
        for x, y, z in zip(xx.flatten(), yy.flatten(), zz2.flatten()):
            new = Atom(f'D{cpt}', np.array([x, y, z]), 0, 1, 32,
                       f' D{cpt} ', cpt)
            mem2.add(new)
            cpt += 1

    # Saving a PDB file that include the membrane position
    filename = 'G:/RAID/Fac/M2_BI/PGP/CompBio/tmp/2n90_extract.pdb'
    io = PDBIO()
    # It's not possible to extract something from another model than 0.
    # TODO: derivate Bio.PDB.Dice.ChainSelector to implement model
    # choice support.
    Dice.extract(prot.structure, prot.chain, prot.res_ids_pdb[0],
                 prot.res_ids_pdb[-1]+2, filename)

    end = time.time() - start
    print('\nDONE in {:.0f} min {:.2f} s.'.format(end // 60, end % 60))
