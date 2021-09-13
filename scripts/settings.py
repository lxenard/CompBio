# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 19:02:23 2021

@author: Laura Xénard
"""

import argparse


parse = False


def init():

    if parse:
        parser = argparse.ArgumentParser()
        parser.add_argument("pdb", help="path of the pdb file to process")
        parser.add_argument("-m", "--model", type=int, default=0,
                            help="model of the protein to process")
        parser.add_argument("-c", "--chain", default='A',
                            help="chain of the protein to process")
        parser.add_argument("-f", "--first_residue", default=None,
                            help="id of the first residue to consider")
        parser.add_argument("-l", "--last_residue", default=None,
                            help="id of the last residue to consider")
        parser.add_argument("-t", "--threshold", type=float, default=0.3,
                            help=("ASA threshold from which residues are considered "
                                  "as exposed to solvent or membrane"))
        parser.add_argument("-n", "--ndir", type=int, default=100,
                            help=("number of directions for sampling the space when "
                                  "looking for membrane plane"))
        parser.add_argument("-s", "--score_method", default='ASA',
                            help="method to use to find the membrane position")
        parser.add_argument("-d", "--debug", action="store_true",
                            help="activate debug mode")
        parser.add_argument("-v", "--verbose", action="store_true",
                            help="increase output verbosity")
        args = parser.parse_args()

        global PDB, MODEL, CHAIN, FIRST_RESIDUE, LAST_RESIDUE
        global IS_EXPOSED_THRESHOLD, N_DIRECTIONS, SCORE_METHOD
        global DEBUG, VERBOSE
        PDB = args.pdb
        MODEL = args.model
        CHAIN = args.chain
        FIRST_RESIDUE = args.first_residue
        LAST_RESIDUE = args.last_residue
        IS_EXPOSED_THRESHOLD = args.threshold
        N_DIRECTIONS = args.ndir
        DEBUG = args.debug
        VERBOSE = args.verbose

        # print(args.verbose)
    else:
        PDB = 'G:/RAID/Fac/M2_BI/PGP/CompBio/data/2n90.pdb'
        #PDB = 'G:/RAID/Fac/M2_BI/PGP/CompBio/data/2lop.pdb'
        #PDB = 'G:/RAID/Fac/M2_BI/PGP/CompBio/data/6g79.pdb'
        #PDB = '/home/sdv/m2bi/lxenard/Documents/PGP/CompBio/data/2n90.pdb'
        MODEL = 0
        CHAIN = 'A'
        FIRST_RESIDUE = None
        LAST_RESIDUE = None
        SCORE_METHOD = 'simple'
        IS_EXPOSED_THRESHOLD = 0.4
        N_DIRECTIONS = 20
        DEBUG = True
        VERBOSE = False
