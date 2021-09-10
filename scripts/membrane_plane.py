# -*- coding: utf-8 -*-
"""
Created on Sun Sep  5 16:56:22 2021

@author: Laura Xénard
"""


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
    prot = ptn.Protein(structure, st.MODEL, st.CHAIN)

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

    to_draw = 10

    mlab.plot3d(prot.vectors[to_draw].get_xx(), prot.vectors[to_draw].get_yy(),
                prot.vectors[to_draw].get_zz(), color=(0, 1, 0),
                tube_radius=None)

    slices = []
    # obtenir les plans orthogonaux centrés sur l'origine
    for v in prot.vectors:
        s = ptn.Slice(prot, 0, v, 'ASA')
        slices.append(s)

    # Drawing all the exposed residues of the protein.
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

    for res in slices[to_draw].residues:
        mlab.points3d(res.coord.x, res.coord.y, res.coord.z,
                      scale_factor=1, color=(0, 1, 0))

    mlab.show()

    print(max(slices))

    # obtenir les plans translatés

    # TODO: pour le renvoi des résultats, ne pas oublier de translater
    # la membrane puisuqe le centre du repère a été déplacé sur le barycentre

    end = time.time() - start
    print('\nDONE in {:.0f} min {:.2f} s.'.format(end // 60, end % 60))