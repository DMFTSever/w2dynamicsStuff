import h5py as hdf5
import numpy as np

class MissmatchError(Exception):

    def __init__(self,message,filename):
        self.message  = message
        self.filename = filename

    def __str__(self):
        string = "Concerning file: {}\n".format(self.filename)
        string += self.message
        return string

class W2dynfile():

    def __init__(self,filename):
        self.filename = filename
        self.file     = hdf5.File(filename,'r')
        self.maxiter  = int(self.file['/.config'].attrs.get('general.dmftsteps'))
        if '/finish' in self.file:
            self.lastiter = self.maxiter
        else:
            for i in range(self.maxiter,0,-1):
                if i==1:
                    print("Found no valid dmft iteration in file: {}"\
                          .format(self.filename))
                    print("Here is what is in the file:")
                    self.walk_file()
                    raise RuntimeError('Aborting...')
                elif '/dmft-%03i'%i in self.file:
                    self.lastiter = i-1 #str('%03i'%(i-1))
                    print('Last valid dmft iteration in file {}: {}'\
                          .format(self.filename,i-1))
                    break
                else:
                    pass
        self.natoms = int(self.file['/.config'].attrs.get('general.nat'))
        self.nineq  = 0
        for element in self.file['/dmft-001'].keys():
            if "ineq" in element:
                self.nineq += 1
        self.equiv = self.file['/.config'].attrs.get('general.equiv')
        if self.equiv == 'none' and self.natoms != self.nineq:
            raise MissmatchError("Unknown quivalence of atoms! Aborting...",self.filename)
        elif self.equiv == 'none':
            self.equiv = np.arrange(1,self.natoms+1)
        else:
            self.equiv = self.equiv.astype(int)+1
        self.niw = self.file['/.config'].attrs.get('qmc.niw')
        self.ntau  = self.file['/.config'].attrs.get('qmc.ntau')

    def walk_file(self,iterations='all'):
        def to_int(string):
            try:
                iteration = int(string.strip("dmft-"))
            except:
                return -1
            return iteration

        if iterations != 'all' and type(iterations)!= list:
            raise MissmatchError("iterations need of to be of type list!",self.filename)

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
            print(atr, self.file.attrs.get(atr))
        if iterations != 'all':
            for x in iterations:
                if int(x)> self.lastiter:
                    print("ITERATION {} LARGER THAN MAXIMAL ITERATION IN FILE".format(x))

    def get_tauax(self,ntau='all'):
        if ntau == 'all':
            ntau = self.ntau
        elif ntau > self.ntau or type(ntau) != int:
            raise MissmatchError("ntau hast to be of type integer < self.ntau",\
                                 self.filename)
        return self.file['/.axes/taubin'][0:ntau]

    def get_iwax(self,niw='all'):
        if niw == 'all':
            niw = self.niw
        elif niw > self.niw or type(niw) != int:
            raise MissmatchError("niw hast to be of type integer < self.niw",\
                                 self.filename)
        return self.file['/.axes/iw'][self.niw-niw:self.niw+niw]

    def get2dim_iw(self,quant,atom,iterations='all',niw='all'):
        if iterations == 'all':
            iterations = np.arange(1,self.lastiter+1)
        elif type(iterations) != list:
            raise MissmatchError("iterations have to be of type list",self.filename)
        if niw == 'all':
            niw = self.niw
        elif niw > self.niw or type(niw) != int:
            raise MissmatchError("niw hast to be of type integer < self.niw",self.filename)
        atstr = str('%03i'%self.equiv[atom-1])
        innershape = self.file['/dmft-001/ineq-{}/{}/value'.format(atstr,quant)].shape[:-1]
        val = np.zeros((len(iterations),)+innershape+(2*niw,),dtype=complex)
        for i in iterations:
            val[i-1,:,:,:]=self.file['/dmft-{:03d}/ineq-{}/{}/value'.format(i,atstr,quant)\
                           ][:,:,self.niw-niw:self.niw+niw]
        return val

    def get2dim_tau(self,quant,atom,iterations='all',ntau='all'):
        if iterations == 'all':
            iterations = np.arange(1,self.lastiter+1)
        elif type(iterations) != list:
            raise MissmatchError("iterations have to be of type list",self.filename)
        if ntau == 'all':
            ntau = self.ntau
        elif ntau > self.ntau or type(ntau) != int:
            raise MissmatchError("ntau hast to be of type integer < self.ntau",\
                                 self.filename)
        atstr = str('%03i'%self.equiv[atom-1])
        innershape = self.file['/dmft-001/ineq-{}/{}/value'.format(atstr,quant)].shape[:-1]
        val = np.zeros((len(iterations),)+innershape+(2*ntau,),dtype=complex)
        for i in iterations:
            val[i-1,:,:,:]=self.file['/dmft-{:03d}/ineq-{}/{}/value'.format(i,atstr,quant)\
                           ][:,:,0:ntau]
        return val
