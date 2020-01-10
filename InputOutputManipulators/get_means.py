#!/bin/usr/env python3

import sys
import os
import argparse
import numpy as np
import h5py as hdf5

parser = argparse.ArgumentParser(description='This script copies all input and averages over the last\
                                              n iterations of the source file and places the it into \
                                              the destination file. If the source file contains not  \
                                              enough valid iterations you will receive an Error.')
parser.add_argument('source', help='This is the source file. It will not be changed.', type=str)
parser.add_argument('dest', help='This is the destination file which will be generated and contain \
                                  the input and averaged objects from the source file.', type=str)
parser.add_argument('iters', help='Iterations to be averaged on. <upper> <lower>', nargs=2, type=int)

args = parser.parse_args()
sourcefile = args.source
destfile = args.dest
limits = args.iters

if sourcefile == destfile:
    raise RuntimeError('Source and destination equal. Please specify another destination.')
if limits[1]<limits[0]:
    raise RuntimeError('Upper limit has to be large than lower limit!')
elif limits[0]<0 or limits[1]<0:
    raise RuntimeError('Limits have to be positive numbers!')

try:
    with hdf5.File(sourcefile,'r') as oldfile:
        with hdf5.File(destfile,'w') as newfile:
            
            maxiter = int(oldfile['/.config'].attrs.get('general.dmftsteps'))
            if limits[1]>maxiter:
                raise RuntimeError(('Upper limit higher then highest iteration in sourcefile. Specify '
                                    'an upper limit that is lower than the highest iteration in the '
                                    'sourcefiel'))
            for i in range(limits[0],limits[1]+1):
                if '/dmft-{:03d}'.format(i) in oldfile:
                    pass
                else:
                    raise RuntimeError(('Did not find iteration {} in source file. Specify valid '
                                        'boundaries!').format(i))
            nineq = 0
            for key in oldfile['/dmft-001/'].keys():
                if key.split('-')[0] == 'ineq':
                    nineq += 1

            for atr in oldfile.attrs.keys():
                newfile.attrs.create(atr,oldfile.attrs.get(atr))
        
            #we also take the start in finish group from the oldfile so hgrep accepts it as input
            for key in oldfile.keys():
                if 'dmft' in key: #or 'finish' in key or 'start' in key:
                    pass
                else:
                    try:
                        oldfile[key].copy(source=oldfile[key], dest=newfile)
                    except:
                        raise RuntimeError('Could not copy group or dataset: '+str(key))

            newfile['/.config'].attrs['z.mean.sourcefile'] = sourcefile
            newfile['/.config'].attrs['z.mean.iterations'] = '{}-{}'.format(limits[0],limits[1])

            newfile.create_group('/dmft-001')
            for atom in range(1,nineq+1):
                newfile.create_group('/dmft-001/ineq-{:03d}'.format(atom))
                for quantity in ['giw-full', 'siw-full', 'gtau-full', 'g0iw-full', 'fiw-full', 'ftau-full']:
                    newfile.create_group('/dmft-001/ineq-{:03d}/{}'.format(atom, quantity))
                    obj = np.array([oldfile['/dmft-{:03d}/ineq-{:03d}/'.format(i,atom)+quantity+'/value'][:] for i in range(limits[0],limits[1]+1)])
                    mean = np.mean(obj, axis=0)
                    err = np.std(obj, axis=0)
                    newfile.create_dataset('/dmft-001/ineq-{:03d}/{}/value'.format(atom, quantity), data=mean)
                    newfile.create_dataset('/dmft-001/ineq-{:03d}/{}/error'.format(atom, quantity), data=err)
except:
    os.remove(destfile)
    raise
    

