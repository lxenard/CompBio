# -*- coding: utf-8 -*-
"""
Created on Sun Sep  5 23:00:28 2021

@author: haiba

Nécessite python 3.7
"""

import math

import numpy as np
from mayavi import mlab



nb = 100
radius = 1

points_x = []
points_y = []
points_z = []
a = 4.0 * math.pi * (radius**2.0 / nb)
d = math.sqrt(a)
m_theta = int(round(math.pi / d))
d_theta = math.pi / m_theta
d_phi = a / d_theta

for m in range(0, m_theta):
    theta = math.pi * (m + 0.5) / m_theta
    m_phi = int(round(2.0 * math.pi * math.sin(theta) / d_phi))
    for n in range(0, m_phi):
        phi = 2.0 * math.pi * n / m_phi
        x = radius * math.sin(theta) * math.cos(phi)
        y = radius * math.sin(theta) * math.sin(phi)
        z = radius * math.cos(theta)
# =============================================================================
#         points_x.append(x)
#         points_y.append(y)
#         points_z.append(z)
# =============================================================================

        if z >= 0:
            points_x.append(x)
            points_y.append(y)
            points_z.append(z)


points = np.zeros((len(points_x), 3))
print(points.shape[0])
for i in range(points.shape[0]):
    points[i, :] = [points_x[i], points_y[i], points_z[i]]
# =============================================================================
#     points[i, 0] = points_x[i]
#     points[i, 1] = points_y[i]
#     points[i, 2] = points_z[i]
# =============================================================================

xx, yy, zz = np.hsplit(points, 3)


# Create a sphere
pi = np.pi
cos = np.cos
sin = np.sin
phi, theta = np.mgrid[0:pi:101j, 0:2 * pi:101j]

x = radius*sin(phi)*cos(theta)
y = radius*sin(phi)*sin(theta)
z = radius*cos(phi)

mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(400, 300))
mlab.clf()




mlab.mesh(x , y , z, color=(0.0,0.5,0.5))
mlab.points3d(xx, yy, zz, scale_factor=0.05)


mlab.show()