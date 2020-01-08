import numpy as np
import h5py

def print_hdf5_tree(item):
    if isinstance(f[item], h5py.Group):
        print item + " -> Group"
    elif isinstance(f[item], h5py.Dataset):
        print item + " -> Dataset"
    else:
        pass
