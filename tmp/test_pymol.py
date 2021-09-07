#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  7 15:14:42 2021

@author: lxenard
"""

import pymol
import __main__
__main__.pymol_argv = ['pymol', '-qc']  # Quiet and no GUI


pymol.finish_launching()

spath = '/home/sdv/m2bi/lxenard/Documents/PGP/CompBio/data/2n90.pdb'
sname = '2N90'

pymol.cmd.load(spath, sname)
pymol.cmd.disable("all")
pymol.cmd.enable(sname)
pymol.cmd.png("/home/sdv/m2bi/lxenard/Documents/PGP/my_image.png")

pymol.cmd.quit()
