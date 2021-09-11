# -*- coding: utf-8 -*-
"""
Created on Sun Sep  5 16:56:22 2021

@author: Laura Xénard
"""

import itertools
from pathlib import Path
import time

from Bio.PDB import PDBParser
from mayavi import mlab
import numpy as np

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
    prot.move(bary)

    mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(600, 600))
    mlab.clf()

    mlab.points3d(0, 0, 0, scale_factor=0.8, color=(1, 0, 0))
# =============================================================================
#     vectors = []
#     for point in sphere.surf_pts:
#         #mlab.points3d(point.x, point.y, point.z, scale_factor=0.05)
#         vectors.append(ptn.Vector(point))
# =============================================================================
    # mlab.show()

    to_draw = 1

# =============================================================================
#     mlab.plot3d(prot.vectors[to_draw].get_xx(), prot.vectors[to_draw].get_yy(),
#                 prot.vectors[to_draw].get_zz(), color=(0, 1, 0),
#                 tube_radius=None)
# =============================================================================

    centered_slices = []
    # obtenir les plans orthogonaux centrés sur l'origine
    for v in prot.vectors:
        s = ptn.Slice(prot, 0, v, st.SCORE_METHOD)
        centered_slices.append(s)

    # Drawing all the exposed residues of the protein.
    for res in prot.residues_exposed:
        mlab.points3d(res.coord.x, res.coord.y, res.coord.z,
                      scale_factor=1, color=(0.5, 0, 0.5))

# =============================================================================
#     a = centered_slices[to_draw].normal.end.x
#     b = centered_slices[to_draw].normal.end.y
#     c = centered_slices[to_draw].normal.end.z
#     x, y = np.mgrid[-20:20:1000j, -20:20:1000j]
#     z = (-a*x - b*y + -7) / c
#     zz = (-a*x - b*y + 7) / c
#     mlab.surf(x, y, z)
#     mlab.surf(x, y, zz)
# =============================================================================

# =============================================================================
#     for res in centered_slices[to_draw].residues:
#         mlab.points3d(res.coord.x, res.coord.y, res.coord.z,
#                       scale_factor=1, color=(0, 1, 0))
# =============================================================================

    #mlab.show()

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

# =============================================================================
#     for sli in slices[:100]:
#         print(sli.score)
#
#     exit
# =============================================================================

    to_draw = best_index
    shift = slices[to_draw].center
    thickness = slices[to_draw].thickness
    a = slices[to_draw].normal.end.x
    b = slices[to_draw].normal.end.y
    c = slices[to_draw].normal.end.z
    x, y = np.mgrid[-20:20:1000j, -20:20:1000j]
    z = (-a*x - b*y - thickness[0] + shift) / c
    zz = (-a*x - b*y + thickness[1] + shift) / c
    mlab.surf(x, y, z, color=(0, 0, 0.8))
    mlab.surf(x, y, zz, color=(0, 0, 0.8))

    for res in slices[to_draw].residues:
        mlab.points3d(res.coord.x, res.coord.y, res.coord.z,
                      scale_factor=1, color=(0, 1, 1))

    mlab.plot3d(slices[to_draw].normal.get_xx(), slices[to_draw].normal.get_yy(),
                slices[to_draw].normal.get_zz(), color=(0, 1, 1),
                tube_radius=None)

    mlab.show()

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


    # TODO: pour le renvoi des résultats, ne pas oublier de translater
    # la membrane puisque le centre du repère a été déplacé sur le barycentre

    end = time.time() - start
    print('\nDONE in {:.0f} min {:.2f} s.'.format(end // 60, end % 60))
