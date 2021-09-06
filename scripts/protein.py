# -*- coding: utf-8 -*-
"""
Created on Sun Sep  5 16:55:30 2021

@author: Laura Xénard

https://www.cmu.edu/biolphys/deserno/pdf/sphere_equi.pdf
"""

import math

class Point:

    def __init__(self, x=0, y=0, z=0):
        self.x = x
        self.y = y
        self.z = z

    def __repr__(self):
        return f"[{self.x}, {self.y}, {self.z}]"

    def __add__(self, p):
        x = self.x + p.x
        y = self.y + p.y
        z = self.z + p.z
        return Point(x, y, z)

    def __sub__(self, p):
        x = self.x - p.x
        y = self.y - p.y
        z = self.z - p.z
        return Point(x, y, z)


class Vector:

    start = Point(0, 0, 0)

    def __init__(self, end=Point(0, 0, 0)):
        self.end = end

    def __repr__(self):
        return f"(start: {self.start}, end: {self.end})"

    def get_xx(self):
        return [self.start.x, self.end.x]

    def get_yy(self):
        return [self.start.y, self.end.y]

    def get_zz(self):
        return [self.start.z, self.end.z]

class Sphere:

    origin = Point(0, 0, 0)

    def __init__(self, radius=1):
        self.radius = radius
        self.surf_pts = []

    def sample_surface(self, nb):
        """
        Generate around 'nb' equidistributed Points on the surface of the
        demi z-positive Sphere.
        This code implements the Deserno algorithm:
        https://www.cmu.edu/biolphys/deserno/pdf/sphere_equi.pdf
        and is based on dinob0t's code:
        https://gist.github.com/dinob0t/9597525

        Parameters
        ----------
        nb : int
            Number of Points to be generated.

        Raises
        ------
        ValueError
            Raised when the number of Points to generate is inferior or equal
            to 0.

        Returns
        -------
        int
            The number of Points that have been generated. It may not be
            exactly equal to 'nb' but should be close to it.

        """
        if nb <= 0:
            raise ValueError

        a = 4.0 * math.pi * (self.radius**2.0 / nb)
        d = math.sqrt(a)
        m_theta = int(round(math.pi / d))
        d_theta = math.pi / m_theta
        d_phi = a / d_theta

        for m in range(0, m_theta):
            theta = math.pi * (m + 0.5) / m_theta
            m_phi = int(round(2.0 * math.pi * math.sin(theta) / d_phi))
            for n in range(0, m_phi):
                phi = 2.0 * math.pi * n / m_phi
                x = self.radius * math.sin(theta) * math.cos(phi)
                y = self.radius * math.sin(theta) * math.sin(phi)
                z = self.radius * math.cos(theta)
                if z>= 0:
                    self.surf_pts.append(Point(x, y, z))
        return len(self.surf_pts)


class Residue:

    def __init__(self, num=0, aa='', p=Point(), asa=0):
        self.num = num
        self.aa = aa
        self.coord = p
        self.asa = asa

    def __repr__(self):
        return f"({self.num}, {self.aa}, coord: {self.coord}, asa: {self.asa})"

    def is_hydrophobic(self):
        """
        Determine if the residue is hydrophobic or not.

        Raises
        ------
        ValueError
            Raised when the amino acid of the residue is not valid.

        Returns
        -------
        bool
            True if the residue is hydrophobic, False otherwise.

        """
        hydrophobic = ('PHE', 'GLY', 'ILE', 'LEU', 'MET', 'VAL', 'TRP', 'TYR')
        hydrophilic = ('ALA', 'CYS', 'ASP', 'GLU', 'HIS', 'LYS', 'ASN', 'PRO',
                       'GLN', 'ARG', 'SER', 'THR')
        if self.aa not in hydrophobic and self.aa not in hydrophilic:
            raise ValueError

        if self.aa in hydrophobic:
            return True
        else:
            return False

    def is_exposed(self, threshold=0.3):
        """
        Determine if the residue is exposed to solvent or membrane, or if
        it is burrowed in the protein.

        Parameters
        ----------
        threshold : float, optional
            Threshold defining is a residue is exposed or burrowed. The
            default is 0.3.

        Returns
        -------
        bool
            True if the residue is exposed to solvent or membrane, False if
            it is burrowed in the protein.

        """
        if self.asa >= threshold:
            return True
        else:
            return False


class Slice:

    def __init__(self, center, normal):
        self.center = center
        self.normal = normal
        self.thickness = 14
        self.score = 0
        self.residues = []

    def find_residues(self, residues):

        for res in residues:

            d = self.thickness / 2

            # Normal vector
            a = self.normal.end.x
            b = self.normal.end.y
            c = self.normal.end.z

            # Plane vector
            x = res.coord.x
            y = res.coord.y
            z = res.coord.z



            if a*x + b*y + c*z <= d and a*x + b*y + c*z >= -d:
                self.residues.append(res)
                print('True')
            else:
                print('False')

# =============================================================================
#         # Equation du plan supérieur :
#         a*x + b*y + z*c + d = 0
#         # Equation du plan inférieur :
#         a*x + b*y + z*c - d= 0
#
#         z = (-a*x - b*y + d) / c
# =============================================================================


    def compute_score(self):
        pass

    def thicken(self, increment=1, direction='top'):
        """


        Parameters
        ----------
        increment : float, optional
            How much to thicken the Slice. The default is 1.
        direction : {'top', 'bottom'}, optional
            In which direction to thicken the Slice. The default is 'top'.

        Returns
        -------
        None.

        """
        pass