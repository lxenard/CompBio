# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 19:02:23 2021
@author: Laura Xénard

This module parses the option when using the program membrane_plane.py.

"""


import argparse



def init():
    """
    Parse the arguments and initialise the parameters.

    Returns
    -------
    None.

    """

    parser = argparse.ArgumentParser()
    parser.add_argument("pdb", help="path of the pdb file to process")
    parser.add_argument("-m", "--model", type=int, default=0,
                        help="model of the protein to process")
    parser.add_argument("-c", "--chain", default='A',
                        help="chain of the protein to process")
    parser.add_argument("-f", "--first_residue", type=int, default=None,
                        help="id of the first residue to consider")
    parser.add_argument("-l", "--last_residue", type=int, default=None,
                        help="id of the last residue to consider")
    parser.add_argument("-t", "--threshold", type=float, default=0.4,
                        help=("ASA threshold from which residues are "
                              "considered as exposed to solvent or membrane"))
    parser.add_argument("-n", "--ndir", type=int, default=100,
                        help=("number of directions for sampling the space "
                              "when looking for membrane plane"))
    parser.add_argument("-s", "--score_method", default='simple',
                        help="method to use to find the membrane position")
    parser.add_argument("-i", "--thicken_increment", default=1,
                        help="how much to increment the membrane thickness "
                        "when looking for the maximal scoring thickness")
    parser.add_argument("-d", "--debug", action="store_true",
                        help="activate debug mode")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="increase output verbosity")
    args = parser.parse_args()

    # Data parameters.
    global PDB, MODEL, CHAIN, FIRST_RESIDUE, LAST_RESIDUE
    PDB = args.pdb
    MODEL = args.model
    CHAIN = args.chain
    FIRST_RESIDUE = args.first_residue
    LAST_RESIDUE = args.last_residue

    # Computing parameters.
    global IS_EXPOSED_THRESHOLD, N_DIRECTIONS, SCORE_METHOD, THICKEN_INCREMENT
    IS_EXPOSED_THRESHOLD = args.threshold
    N_DIRECTIONS = args.ndir
    SCORE_METHOD = args.score_method
    THICKEN_INCREMENT = args.thicken_increment

    global DEBUG, VERBOSE
    DEBUG = args.debug
    VERBOSE = args.verbose
    if DEBUG:
        VERBOSE = True
