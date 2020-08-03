#!/usr/bin/env python

from __future__ import print_function, division, absolute_import
import numpy as np
import h5py
import argparse
import os
import sys
### here we need to add the path to the personal myfunc directory (location of some read functions)
sys.path.insert(0,sys.argv[0].replace('/InputOutputManipulators/generate_kanamori_umatrix.py','/General'))
### in order to import the functions available in input.py and interatction.py and readwrite and custom errors
import readwrite as rw
import custom_errors as err

parser = argparse.ArgumentParser(description='This script generate an umatrix file for the Kanamori \
                                              interaction in w2dynamics file format.')
parser.add_argument('ufile', help='Name of the file containing the umatrix.')
parser.add_argument('U', help='U value.')
parser.add_argument('J', help='J value.')
parser.add_argument('V', help='V value.')
parser.add_argument('norb', help='Number of orbitals.')

args = parser.parse_args()
destfile = args.ufile
U = args.U
J = args.J
V = args.V
norb = int(args.norb)

u_matrix = np.zeros(shape=(norb,2,norb,2,norb,2,norb,2), dtype=float)

for b1 in range(0,norb):
    u_matrix[b1,(0,1),b1,(1,0),b1,(0,1),b1,(1,0)] = U
    for b3 in range(0,norb):
        if b3 != b1:
            u_matrix[b1,(0,0,1,1),b3,(0,1,0,1),b1,(0,0,1,1),b3,(0,1,0,1)] = V
            u_matrix[b1,(0,1),b3,(0,1),b3,(0,1),b1,(0,1)] = J

for b1 in range(0,norb):
    for b2 in range(norb):
        if b1 != b2:
            u_matrix[b1,(0,1),b1,(1,0),b2,(0,1),b2,(1,0)] = J # pair hopping
            u_matrix[b1,(0,1),b2,(1,0),b2,(0,1),b1,(1,0)] = J # spin flip

rw.write_umatrix(filename=destfile, umatrix=u_matrix, spin=True)
