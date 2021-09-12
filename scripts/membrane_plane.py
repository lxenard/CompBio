# -*- coding: utf-8 -*-
"""
Created on Sun Sep  5 16:56:22 2021

@author: Laura Xénard
"""


from pathlib import Path
import time

from Bio.PDB import PDBParser

import protein as ptn
import settings as st



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
    best = next(((i, sli) for i, sli in enumerate(slices)
                 if sli.score < 1), None)
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

    # Adding atoms representing the membrane position to the protein
    # and saving the updated protein as a PDB file.
    prot.add_membrane(best_sli)
    filename = (f"{ptn_id}_membrane_m{prot.model}_c{prot.chain}_"
                f"fr{prot.res_ids_pdb[0]}_lr{prot.res_ids_pdb[-1]}.pdb")
    filepath = Path(st.PDB).parent.joinpath(filename)
    prot.save_pdb(str(filepath))

    end = time.time() - start
    print('\nDONE in {:.0f} min {:.2f} s.'.format(end // 60, end % 60))
