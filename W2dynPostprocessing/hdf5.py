import h5py as hdf5
import numpy as np

class W2dynfile():

    def __init__(self,filename):
        self.filename = filename
        self.file = hdf5.File(filename,'r')

    def walk_file(self,iterations='all'):
        def to_int(string):
            try:
                iteration = int(string.strip("dmft-"))
            except:
                return -1
            return iteration

        if iterations != 'all' and type(iterations)!= list:
            raise RuntimeError("iterations of to be of type list!")
        maxiter = int(self.file['/.config'].attrs.get('general.dmftsteps'))
        if '/finish' in self.file:
            pass
        else:
            for i in range(maxiter,0,-1):
                print('/dmft-%3i'%i)
                if i==1:
                    raise RuntimeError('Found no valid dmft iteration to copy. Aborting...')
                elif '/dmft-%03i'%i in self.file:
                    maxiter = str('%03i'%(i-1))
                    print('Last dmft iteration in file: ', i )
                    break
                else:
                    pass
        for key_level1 in self.file.keys():
            if iterations!='all' and 'dmft' in key_level1 \
                                 and to_int(key_level1) not in iterations:
                pass
            else:
                print("{:20}".format(key_level1), self.file[key_level1])
                try:
                    for key_level2 in self.file[key_level1].keys():
                        print("  {:18}".format(key_level2), \
                              self.file[key_level1][key_level2])
                        try:
                            for key_level3 in self.file[key_level1][key_level2].keys():
                                print("    {:16}".format(key_level3),\
                                      self.file[key_level1][key_level2][key_level3])
                                try:
                                    for key_level4 in self.file[key_level1][key_level2]\
                                                               [key_level3].keys():
                                        print("      {:14}".format(key_level4),\
                                              self.file[key_level1][key_level2]\
                                                       [key_level3][key_level4])
                                except:
                                    pass
                        except:
                            pass
                except:
                    pass
                print("")
        for atr in self.file.attrs.keys():
            print(self.file.attrs.get(atr))
        for x in iterations:
            if x!='all' and int(x)>maxiter:
                print("ITERATION {} LARGER THAN MAXIMAL ITERATION IN FILE".format(x))
