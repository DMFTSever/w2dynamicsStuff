#!/bin/python

from __future__ import print_function, division, absolute_import
import numpy as np
import h5py
import argparse
import os
import sys
### here we need to add the path to the personal myfunc directory (location of some read functions)
sys.path.insert(0,sys.argv[0].replace('/InputOutputManipulators/getOccDiagonals.py','/General'))
import hdf5_stuff as h5s

parser = argparse.ArgumentParser(description='This script print the diagonal occupancies, \
                                              and the differences for spin up and down per band.')
parser.add_argument('source', help='This is the source file. It will not be changed.')
parser.add_argument('-fg', '--fromGtau', help='Occs will be taken from Gtau(0), where 0 is the\
                                               smallest matsubara frequency.', \
                    action='store_true', default=False)
parser.add_argument('-i', '--iteration', type=int, default=0, help='Enforce to use iteration i. Overrides -ul')
parser.add_argument('-ul', '--useLast', help='Enforce to use last iteration found\
                                                        not second to last (default).', \
                    action='store_true', default=False)

args = parser.parse_args()
sourcefile = args.source
fromGtau = args.fromGtau
iteration = args.iteration
useLast = args.useLast

if sourcefile == 'latest':
    try:
        sourcefile = h5s.get_latest_hdf5()
    except:
        raise RuntimeError("Could not identify latest hdf5 file")
else:
    pass

with h5py.File(sourcefile,'r') as source:
        if iteration != 0:
            if '/dmft-{:03}'.format(iteration) in source:
                iteration = "{:03d}".format(iteration)
            else:
                raise RuntimeError("Iteration {} not in hdf5 file!".format(iteration))
        else:
            if useLast:
                iteration = "{:03d}".format(h5s.get_latest_dmft_iter(source, safe=False))
            else:
                iteration = "{:03d}".format(h5s.get_latest_dmft_iter(source, safe=True))
        print('Using dmft iteration: {}'.format(iteration))
        nat = source['/.config'].attrs.get('general.nat')
        nds = np.zeros((nat,),dtype=int)
        for a in range(0,nat):
            nds[a] = source['/.config'].attrs.get('atoms.'+str(a+1)+'.nd')
        print("\nDiagonal entries of rho1 with differences:\n")
        for a in range(0,nat):
            print("atom {}:\n".format(a+1))
            try:
                if fromGtau:
                    gtau = source['/dmft-{}/ineq-{:03}/gtau-full/value'.format(iteration,a+1)][:]
                    rho = 1-gtau[:,:,:,:,0]
                else:
                    rho = source['/dmft-{}/ineq-{:03}/rho1/value'.format(iteration,a+1)][:]
            except:
                raise RuntimeError("Could not extract density from file {} at iteration {}".format(sourcefile, iteration))
            for b in range(0,nds[a]):
                    print('  band {:2} up    : {:.5}'.format(b+1, rho[b,0,b,0].real))
                    print('  band {:2} down  : {:.5}'.format(b+1, rho[b,1,b,1].real))
                    print('  Del (up-down) : {:.5}\n'.format(rho[b,0,b,0].real-rho[b,1,b,1].real))
            print("  Mean Del      : {}\n".format(1/nds[a]*(rho[:,0,:,0].trace().real-rho[:,1,:,1].trace().real)))
            

