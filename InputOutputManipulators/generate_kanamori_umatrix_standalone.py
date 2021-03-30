#!/usr/bin/env python

from __future__ import print_function, division, absolute_import
import numpy as np
import h5py
import argparse
import sys
import re

class ShapeError(Exception):
    """Exeption raised fo error in the shaping.
    """

def write_umatrix(filename, umatrix, spin):
    '''Funtion to write the umatrix either with spin dependency or without'''
    try:
        if type(filename) is not str:
            raise TypeError("1. argument (filenname) must be of type str.")
        if type(umatrix) is not np.ndarray:
            raise TypeError("2. argument (umatrix) must be of type numpy.ndarray.")
        if type(spin) is not bool:
            raise TypeError("3. argument (spin) must be of type bool.")
        if len(umatrix.shape)!=4 and len(umatrix.shape)!=8:
            raise ShapeError("Shape of umatrix is invalid")
        if len(umatrix.shape)==4 and spin:
            raise ShapeError("Spin independent umatrix detected although spin"\
                                 " was set to True.")
        if len(umatrix.shape)==8 and not spin:
            raise ShapeError("Spin dependent umatrix detected although spin"\
                                 " was set to False.")
        if spin and (umatrix.shape[1]!=2 or umatrix.shape[3]!=2 or umatrix.shape[5]!=2\
                     or umatrix.shape[7]!=2):
            raise ShapeError("Unexpected shape of spin dependent umatrix. Expected"\
                             " shape: (Nbands,2,Nbands,2,Nbands,2,Nbands,2)")
    except TypeError:
        raise
        sys.exit()
    except ShapeError:
        raise
        sys.exit()

    f = open(filename, 'w')
    f.write("# non zero elements of interaction matrix U_ijkl\n")
    f.write("%i BANDS\n"%umatrix.shape[0])
    spindict={0:'u',1:'d'}
    if not spin:
        for index, value in np.ndenumerate(umatrix):
            if value != 0.:
                f.write("%i %i %i %i  %016.16f\n"%(index[0]+1,index[1]+1,index[2]+1,\
                                                   index[3]+1, value))
    else:
        for index, value in np.ndenumerate(umatrix):
            if value != 0:
                f.write("%s %s %s %s  %016.16f\n"%(str(index[0]+1)+spindict[index[1]],\
                                                   str(index[2]+1)+spindict[index[3]],\
                                                   str(index[4]+1)+spindict[index[5]],\
                                                   str(index[6]+1)+spindict[index[7]],\
                                                   value))

parser = argparse.ArgumentParser(description='This script generate an umatrix file for the Kanamori \
                                              interaction in w2dynamics file format.')
parser.add_argument('ufile', help='Name of the output file containing the umatrix'\
                                  ' without fileending.')
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

#force fileending to be .dat
fname = re.sub(".out","",re.sub(".txt","",re.sub(".dat","",destfile)))
fname = fname + ".dat"
write_umatrix(filename=fname, umatrix=u_matrix, spin=True)
