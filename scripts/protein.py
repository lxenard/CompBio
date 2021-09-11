# -*- coding: utf-8 -*-
"""
Created on Sun Sep  5 16:55:30 2021

@author: Laura Xénard

https://www.cmu.edu/biolphys/deserno/pdf/sphere_equi.pdf
"""

import math

from Bio.PDB.DSSP import DSSP

import settings as st


class Point:
    """
    Represent a point in 3D space.

    Parameters
    ----------
    x : float, optional
        The x coordinate. The default is 0.
    y : float, optional
        The y coordinate. The default is 0.
    z : float, optional
        The z coordinate. The default is 0.

    Attributes
    ----------
    x : float
        The x coordinate.
    y : float
        The y coordinate.
    z : float
        The z coordinate.

    """
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

    @classmethod
    def barycenter(cls, residues_list):
        """
        Compute the barycenter of a list of residues.

        Parameters
        ----------
        residues_list : list(ptn.Residues)
            The residues from which to find the barycenter.

        Returns
        -------
        Point
            The barycenter of the Residues.

        """
        x_bary = sum([r.coord.x for r in residues_list]) / len(residues_list)
        y_bary = sum([r.coord.y for r in residues_list]) / len(residues_list)
        z_bary = sum([r.coord.z for r in residues_list]) / len(residues_list)
        return Point(x_bary, y_bary, z_bary)


class Vector:
    """
    Represent a vector in 3D space.

    Parameters
    ----------
    end : Point
        The end point of the vector.

    Attributes
    ----------
    start : Point
        The start point of the vector.
    end : Point
        The end point of the vector.

    """
    start = Point(0, 0, 0)

    def __init__(self, end):
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
    """
    Represent a sphere in 3D space.

    Parameters
    ----------
    radius : int, optional
            Radius of the sphere. The default is 1.

    Attributes
    ----------
    origin : Point
        The center of the sphere.
    radius : int, optional
            Radius of the sphere. The default is 1.
    surf_pts : list(Point)
            Points on the surface of the Sphere.

    """
    origin = Point(0, 0, 0)

    def __init__(self, radius=1):
        self.radius = radius
        self.surf_pts = []

    def sample_surface(self, nb):
        """
        Generate equidistributed Points on the surface of the demi Sphere.

        The Points are generated on the z-positive part of the Sphere.

        Parameters
        ----------
        nb : int
            The ideal number of Points to be generated. In some cases, less
            Points might be generated.

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

        Notes
        -----
        This code implements the Deserno algorithm:
        https://www.cmu.edu/biolphys/deserno/pdf/sphere_equi.pdf
        and is based on dinob0t's code:
        https://gist.github.com/dinob0t/9597525

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
                if z >= 0:
                    self.surf_pts.append(Point(x, y, z))
        return len(self.surf_pts)


class Residue:
    """
    TODO

    Parameters
    ----------

    Attributes
    ----------

    """
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
        Determine if the residue is exposed to solvent or burrowed.

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


class Protein():
    """
    TODO

    Parameters
    ----------

    Attributes
    ----------

    """
    def __init__(self, structure, model=0, chain='A', first_residue=None):
        self.structure = structure
        self.model = model
        self.chain = chain

        ids = [d.get_id()[1] for d in structure[0][chain]]
        if first_residue is None:
            self.res_ids_pdb = ids
        else:
            try:
                first_index = ids.index(first_residue)
                self.res_ids_pdb = ids[first_index:]
            except (ValueError, IndexError):
                print(f"WARNING: residue {first_residue} does not exist in "
                      f"model {self.model} chain {self.chain}. "
                      "Starting from the first existing residue instead"
                      f" which is {ids[0]}.")
                self.res_ids_pdb = ids

        self.vectors = []
        self.sample_space()

        self.residues_burrowed = []
        self.residues_exposed = []
        self.find_exposed_residues()

    def sample_space(self):
        """
        Sample the space in several vectors.

        All the vectors pass by the center of the coordinate system.

        Returns
        -------
        None.

        """
        sphere = Sphere()
        sphere.sample_surface(st.N_DIRECTIONS*2)
        for point in sphere.surf_pts:
            self.vectors.append(Vector(point))

    def find_exposed_residues(self):
        """
        Find the residues exposed to solvent or membrane.

        Returns
        -------
        None.

        """
        dssp = DSSP(self.structure[self.model], st.PDB)
        #i_res = self.res_ids_pdb[0]
        ids = self.res_ids_pdb
        for i_res, res in zip(ids, self.structure[self.model][self.chain]):
            # For simplification, the position of a residue is defined as the
            # position of its Cα.
            print(i_res, res)
            pt = Point(*res['CA'].coord)
            asa = dssp[(self.chain, i_res)][3]  # Accessible surface area.
            tmp = Residue(res.id[1], res.resname, pt, asa)
            if tmp.is_exposed(st.IS_EXPOSED_THRESHOLD):
                self.residues_exposed.append(tmp)
            else:
                # Also save burrowed residues in case they are needed later.
                self.residues_burrowed.append(tmp)

    def move(self, shift):
        """
        TODO

        Parameters
        ----------
        shift : Point
            DESCRIPTION.

        Returns
        -------
        None.

        """
        for res in self.residues_exposed:
            res.coord -= shift
        for res in self.residues_burrowed:
            res.coord -= shift


class Slice():
    """
    TODO

    Parameters
    ----------
    protein : TYPE
        DESCRIPTION.
    center : TYPE
        DESCRIPTION.
    normal : TYPE
        DESCRIPTION.
    method : TYPE, optional
        DESCRIPTION. The default is 'ASA'.

    Attributes
    ----------

    """
    def __init__(self, protein, center, normal, method='ASA'):
        self.protein = protein
        self.center = center
        self.normal = normal
        self.score_method = method
        self.thickness = [7, 7]
        self.residues = []
        self.score = 0
        n_res = self.find_residues()
        # If there's no residues in the slice, no need to update the score.
        if n_res != 0:
            try:
                self.compute_score()
            except ValueError:
                print("Method must be 'ASA' or 'simple'")

    def __repr__(self):
        thickness = sum(self.thickness)
        nb_residues = len(self.residues)
        return (f"(center: {self.center}, normal: {self.normal}, "
                f"thickness: {thickness}, nb_residues: {nb_residues}, "
                f"score: {self.score})")

    def __lt__(self, other):
        if self.score < other.score:
            return True
        else:
            return False

    def __gt__(self, other):
        if self.score > other.score:
            return True
        else:
            return False

    def __le__(self, other):
        if self.score <= other.score:
            return True
        else:
            return False

    def __ge__(self, other):
        if self.score >= other.score:
            return True
        else:
            return False

    def find_residues(self):
        """
        Find the exposed Residues that are inside the Slice.

        Returns
        -------
        int
            Number of exposed Residues inside the Slice.

        """
        for res in self.protein.residues_exposed:

            # Normal vector.
            a = self.normal.end.x
            b = self.normal.end.y
            c = self.normal.end.z

            # Plane vector.
            x = res.coord.x
            y = res.coord.y
            z = res.coord.z

            # Position of the planes along the normal vector.
            d1 = self.center - self.thickness[0]
            d2 = self.center + self.thickness[1]

            # The residues between the 2 planes are inside the Slice.
            if a*x + b*y + c*z >= d1 and a*x + b*y + c*z <= d2:
                self.residues.append(res)

        return len(self.residues)

    def compute_score(self, method='ASA'):
        """
        Compute the membrane score of the Slice.

        The higher the score the more probable it is that the Slice
        represents the position of a membrane.

        Parameters
        ----------
        method : {'ASA', 'simple'}, optional
            The method used to compute the Slice score. 'ASA' computes
            the ratio accessible surface area (ASA) of hydrophobic
            residues to ASA of all residues. 'simple' computes the ratio
            number of hydrophobic residues to total number of residues.
            The default is 'ASA'.

        Raises
        ------
        ValueError
            Raised when the method to use is neither 'ASA' nor 'simple'.

        Returns
        -------
        None.

        """
        if method != 'ASA' and method != 'simple':
            raise ValueError

        if method == 'ASA':
            hydrophobic_asa = 0
            total_asa_slice = 0
            total_asa_prot = 0

            for res in self.protein.residues_exposed:
                total_asa_prot += res.asa

            for res in self.residues:
                try:
                    if res.is_hydrophobic():
                        hydrophobic_asa += res.asa
                    total_asa_slice += res.asa
                except ValueError:
                    print(f"Can't determine hydrophobicity of {res}: "
                          f"unknown amino acid.")
            self.score = hydrophobic_asa / total_asa_slice

        elif method == 'simple':
            cpt = 0
            for res in self.residues:
                try:
                    if res.is_hydrophobic():
                        cpt += 1
                except ValueError:
                    print(f"Can't determine hydrophobicity of {res}: "
                          f"unknown amino acid.")
            self.score = cpt / len(self.residues)

    def thicken(self, increment=1, normal_direction=True):
        """
        Thicken the Slice along the direction of the normal vector.

        Parameters
        ----------
        increment : float, optional
            How much to thicken the Slice. The default is 1.
        normal_direction : bool, optional
            In which direction to thicken the Slice. The default is True.
            True corresponds to the normal vector direction, False to the
            opposite direction.

        Returns
        -------
        None.

        """
        if normal_direction:
            self.thickness[1] += increment
        else:
            self.thickness[0] += increment
        # TODO: gérer la maj automatique des résidus compris dans la membrane
        # et le calcul du score.
