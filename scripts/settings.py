# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 19:02:23 2021

@author: Laura Xénard
"""

import argparse

def init():

    parser = argparse.ArgumentParser()
    parser.add_argument("pdb", help="path of the pdb file to process")
    parser.add_argument("-m", "--model", type=int, default=0,
                        help="model of the protein to process")
    parser.add_argument("-c", "--chain", default='A',
                        help="chain of the protein to process")
    parser.add_argument("-t", "--threshold", type=float, default=0.3,
                        help=("ASA threshold from which residues are considered "
                              "as exposed to solvent or membrane"))
    parser.add_argument("-n", "--ndir", type=int, default=20,
                        help=("number of directions for sampling the space when "
                              "looking for membrane plane"))
    parser.add_argument("-d", "--debug", action="store_true",
                        help="activate debug mode")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="increase output verbosity")
    args = parser.parse_args()


    global PDB, MODEL, CHAIN, IS_EXPOSED_THRESHOLD, N_DIRECTIONS
    global DEBUG, VERBOSE
    PDB = args.pdb
    MODEL = args.model
    CHAIN = args.chain
    IS_EXPOSED_THRESHOLD = args.threshold
    N_DIRECTIONS = args.ndir
    DEBUG = args.debug
    VERBOSE = args.verbose


    print(args.verbose)

# =============================================================================
#     pdb_path = 'G:/RAID/Fac/M2_BI/PGP/CompBio/data/2n90.pdb'
#     #pdb_path = '/home/sdv/m2bi/lxenard/Documents/PGP/CompBio/data/2n90.pdb'
#     model = 0
#     chain = 'A'
#     IS_EXPOSED_THRESHOLD = 0.3
#     DEBUG = False
#     N_DIRECTIONS = 20  # nb de points pour échantillonner la demi-sphère
#
#     global myList
#     myList = []
#     mylist.append(arg)
# =============================================================================