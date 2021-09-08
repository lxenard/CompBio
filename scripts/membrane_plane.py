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
from mayavi import mlab
import numpy as np

import settings as st
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
    return ptn.Point(x_bary, y_bary, z_bary)


if __name__ == '__main__':

    start = time.time()

    st.init()  # Set the global parameters.

    # Opening and parsing of the PDB file.
    p = PDBParser()
    ptn_id = Path(st.PDB).stem
    structure = p.get_structure(ptn_id, st.PDB)
    prot = ptn.Protein(structure, st.MODEL, st.CHAIN)

    # For better results, the burrowed residues are not taken into account
    # during the membrane detection, only the exposed ones.
    # TODO: à vérifier en comparant les résultats avec le traitement de tous
    # les résidus
    # Place the center of the coordinate system at the residues barycenter.
    bary = barycenter(prot.residues_exposed)
    if st.DEBUG:
        print(f"Barycenter: {bary}")
    prot.move(bary)

    # Sample the space in roughly N_DIRECTIONS vectors all passing by the
    # center of the coordinate system.
    sphere = ptn.Sphere()
    sphere.sample_surface(st.N_DIRECTIONS*2)

    mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(600, 600))
    mlab.clf()

    mlab.points3d(0, 0, 0, scale_factor=0.1, color=(1, 0, 0))
    vectors = []
    for point in sphere.surf_pts:
        mlab.points3d(point.x, point.y, point.z, scale_factor=0.05)
        vectors.append(ptn.Vector(point))
    # mlab.show()

    to_draw = 8

    mlab.plot3d(vectors[to_draw].get_xx(), vectors[to_draw].get_yy(), vectors[to_draw].get_zz(),
                color=(0, 1, 0), tube_radius=None)

    slices = []
    # obtenir les plans orthogonaux
    for v in vectors:
        # TODO: a la création d'une slice les résidues de la prot sont
        # recalculées sans prendre en compte le shift
        # => revoir les classes slice / prot
        s = ptn.Slice(0, v, structure, st.MODEL, st.CHAIN)
        slices.append(s)

    for res in prot.residues_exposed:
        mlab.points3d(res.coord.x, res.coord.y, res.coord.z,
                      scale_factor=1, color=(0.5, 0, 0.5))

    a = slices[to_draw].normal.end.x
    b = slices[to_draw].normal.end.y
    c = slices[to_draw].normal.end.z
    x, y = np.mgrid[-20:20:1000j, -20:20:1000j]
    z = (-a*x - b*y + -7) / c
    zz = (-a*x - b*y + 7) / c
    mlab.surf(x, y, z)
    mlab.surf(x, y, zz)

    for s, sli in enumerate(slices):
        n = sli.find_residues()
        print(n)
# =============================================================================
#         try:
#             # traiter le cas où il n'y a aucun résidus dans la tranche
#             sli.compute_score()
#         except ValueError:
#             print("Method must be 'ASA' or 'simple'")
#         print(f"Score slice {s} : {sli.score}")
# =============================================================================

    for res in slices[to_draw].residues:
        mlab.points3d(res.coord.x, res.coord.y, res.coord.z,
                      scale_factor=1, color=(0, 1, 0))

    mlab.show()

# =============================================================================
#     for s, sli in enumerate(slices):
#         sli.find_residues(residues)
#         sli.compute_score()
#         print(f"Score slice {s} : {sli.score}")
# =============================================================================

    print(max(slices))

    # obtenir les plans translatés
# =============================================================================
#     slices[0].find_residues(residues)
#     for res in slices[0].residues:
#         mlab.points3d(res.coord.x, res.coord.y, res.coord.z,
#                       scale_factor=1.1, color=(0, 1, 0))
#
#     mlab.show()
#     slices[0].compute_score()
#     print(f"Score slice 0 : {slices[0].score}")
# =============================================================================

    # TODO: pour le renvoi des résultats, ne pas oublier de translater
    # la membrane puisuqe le centre du repère a été déplacé sur le barycentre

    end = time.time() - start
    print('\nDONE in {:.0f} min {:.2f} s.'.format(end // 60, end % 60))