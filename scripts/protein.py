# -*- coding: utf-8 -*-
"""
Created on Sun Sep  5 16:55:30 2021

@author: Laura Xénard
"""

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


class Residue:

    def __init__(self, num=0, aa='', p=Point(), asa=0):
        self.num=num
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

    def __init__(self):

        pass

    def find_residues(self):
        pass

    def compute_score(self):
        pass

    def thicken(self):
        pass