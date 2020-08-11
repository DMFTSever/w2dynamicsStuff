from __future__ import print_function, division, absolute_import
import numpy as np
import h5py
import os

def print_hdf5_tree(item):
    if isinstance(f[item], h5py.Group):
        print(item + " -> Group")
    elif isinstance(f[item], h5py.Dataset):
        print(item + " -> Dataset")
    else:
        pass

def get_latest_hdf5(directory=os.getcwd()):
    return max([os.path.join(directory, fn) for fn in os.listdir(directory) if h5py.is_hdf5(os.path.join(directory,fn))],\
               key=os.path.getctime)

def get_latest_dmft_iter(source, safe=True):
    maxiter = int(source['/.config'].attrs.get('general.dmftsteps'))
    if '/finish' in source:
        return maxiter
    else:
        for i in range(maxiter,0,-1):
            # return i-1 to ensure that dmft iteration has completted!
            if i==1 and safe:
                raise RuntimeError('Found no valid safe dmft iteration. Aborting...')
            elif '/dmft-%03i'%i in source:
                if safe:
                    return i-1
                else:
                    return i
            elif i==1:
                raise RuntimeError('Found no valid dmft iteration. Aborting...')
            else:
                pass
