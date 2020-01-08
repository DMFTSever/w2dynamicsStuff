##################################################################################
# This script rotates an object into or from the basis in which Hloc is diagonal #
##################################################################################

import numpy as np
#import matplotlib.pyplot as plt
import sys
import os
import argparse
import h5py as hdf5
sys.path.insert(0,sys.argv[0].replace('/InputOutputManipulators/rotate_objects.py','/General'))
import readwrite as rw 
import custom_errors as err 

class shift(argparse.Action):
    def __call__(self, parser, namespace, values, option_sting=None):
        newvalues = []
        for value in values:
            newvalues.append(value-1)
        setattr(namespace, self.dest, newvalues)

parser = argparse.ArgumentParser(description="This script rotates objects depending on tau either in the basis where\
                                              Hloc is diagonal or from the basis where Hloc is diagonal into the LS\
                                              Basis. On default it rotates gtau. ")
parser.add_argument('hkfile', help='PATH to the k dependent Hamiltonian. Needs to be spin dependend', type=str)
parser.add_argument('Natoms', help='Number of atoms', type=int)
parser.add_argument('Nbands', help='Number of bands per atom', type=int)
parser.add_argument('hdf5file', help='PATH to the HDF5 file with the w2dynamics results', type=str)
parser.add_argument('--toDIAG', default=False, action='store_true', help="Rotate from the LS to the diagonal basis. On\
                                                                          default it the scripts rotates from the \
                                                                          diagonal to the LS basis.")
parser.add_argument('-A', '--Atoms', default=[0], type=int, dest='Atoms', nargs='*', action=shift,\
                    help="Objects of this atoms will be rotated.")
parser.add_argument('--notgtau', default=True, action='store_false', help="Specify in order not to rotated gtau. In \
                                                                           default the script rotates gtau.") 
parser.add_argument('--ftau', default=False, action='store_true', help="Specify in order to rotated ftau.")
parser.add_argument('--giw', default=False, action='store_true', help="Specify in order to rotated giw.")
parser.add_argument('--siw', default=False, action='store_true', help="Specify in order to rotated siw.")
parser.add_argument('-s', '--spin', default=False, action='store_true', help="Specify that Hamiltonian has spin \
                                                           dependency.")
parser.add_argument('-i', '--iteration', default=None, type=int, help="Specify which iteration of hdf5 file to use.")
args = parser.parse_args()
hkfile = args.hkfile
Natoms = args.Natoms
Nbands = args.Nbands
hdf5file = args.hdf5file
atoms = args.Atoms
toDIAG = args.toDIAG
rotate_gtau = args.notgtau
rotate_ftau = args.ftau
rotate_giw = args.giw
rotate_siw = args.siw
iteration = args.iteration
Nspins = int(args.spin) + 1

for atom in atoms:
    if atom < 0 or atom > Natoms-1:
        raise err.InputError("Atom list out of range!")
        sys.exit()

if rotate_gtau == False:
    rot_objs = []
else:
    rot_objs = ["gtau"]
if rotate_ftau == True:
    rot_objs.append("ftau")
if rotate_giw == True:
    rot_objs.append("giw")
if rotate_siw == True:
    rot_objs.append("siw")

if Nspins == 2:
    print "Loading spin dependent hamiltonian on hk and kpoints on kpoints."
    hk, kpoints = rw.read_hk_wannier(hkfile, spin=True)
else:
    print "Loading spin independent hamiltonian on hk and kpoints on kpoints."
    hk, kpoints = rw.read_hk_wannier(hkfile, spin=False)
Nk=kpoints.shape[0]
print "hk.shape", hk.shape
print "kpoints.shape", kpoints.shape


print "Building hkmean i.e. averaging over all k points"
hkmean = 1./Nk * np.sum(hk, axis=0)
print "hkmean.shape", hkmean.shape
hkmean = hkmean.reshape(Nspins*Nbands*Natoms,Nspins*Nbands*Natoms)
print "hkmean.shape", hkmean.shape
print "Building hkmean_pa (hkmean per atom)"
hkmean_pa = np.zeros((Natoms,Nspins*Nbands,Nspins*Nbands),dtype=complex)
for i in range(0,Natoms):
    hkmean_pa[i] = hkmean[i*Nspins*Nbands:(i+1)*Nspins*Nbands,i*Nspins*Nbands:(i+1)*Nspins*Nbands]

print "hkmean_pa.shape", hkmean_pa.shape
    
if Nspins == 2:
    for atom in atoms:
        print "hkmean_pa[" + str(atom) + "]"
        for i in range(0,Nspins*Nbands):
            print '%+05.5f%+05.5fi' %(np.real(hkmean_pa[atom,i,0]), np.imag(hkmean_pa[atom,i,0])), \
                  '%+05.5f%+05.5fi' %(np.real(hkmean_pa[atom,i,1]), np.imag(hkmean_pa[atom,i,1])), \
                  '%+05.5f%+05.5fi' %(np.real(hkmean_pa[atom,i,2]), np.imag(hkmean_pa[atom,i,2])), \
                  '%+05.5f%+05.5fi' %(np.real(hkmean_pa[atom,i,3]), np.imag(hkmean_pa[atom,i,3])), \
                  '%+05.5f%+05.5fi' %(np.real(hkmean_pa[atom,i,4]), np.imag(hkmean_pa[atom,i,4])), \
                  '%+05.5f%+05.5fi' %(np.real(hkmean_pa[atom,i,5]), np.imag(hkmean_pa[atom,i,5])) 
else:
    for atom in atoms:
        print "hkmean_pa[" + str(atom) + "]"
        for i in range(0,Nspins*Nbands):
            print '%+05.5f%+05.5fi' %(np.real(hkmean_pa[atom,i,0]), np.imag(hkmean_pa[atom,i,0])), \
                  '%+05.5f%+05.5fi' %(np.real(hkmean_pa[atom,i,1]), np.imag(hkmean_pa[atom,i,1])), \
                  '%+05.5f%+05.5fi' %(np.real(hkmean_pa[atom,i,2]), np.imag(hkmean_pa[atom,i,2])), \

print "Solving the Eigenvalue problem for hkmean_pa"
EigVal = np.zeros((Natoms,Nspins*Nbands),dtype=complex)
EigVec = np.zeros_like(hkmean_pa)
for i in range(0,Natoms):
	EigVal[i], EigVec[i] = np.linalg.eigh(hkmean_pa[i])
	print "Eigenvalues atom ", i+1
	print EigVal[i]

#test for non spin rotation
if Nspins == 1:
    print "Blowing up Transformation matrix to two spins"
    temp = EigVec.copy()
    EigVec = np.zeros((Natoms,Nbands,2,Nbands,2), dtype=complex)
    EigVec[:,:,0,:,0] = temp[:]
    EigVec[:,:,1,:,1] = temp[:]
    print "EigVec.shape ", EigVec.shape
    EigVec = EigVec.reshape(Natoms,2*Nbands,2*Nbands)
    print "EigVec.shape ", EigVec.shape

print "Loading hdf5 file: " + hdf5file
f=hdf5.File(hdf5file,"r")
if iteration == None:
    print 'Using last iteration'
    iteration = "%03i" %int(f["/.config"].attrs.get("general.dmftsteps")) #using the last Iteration
else:
    print 'Using %03i iteration' %iteration
    iteration = "%03i" %iteration

#Checking for directory to put results
if toDIAG:
    basis = 'hloc'
else:
    basis = 'LS'
basedir = "./" + basis + "_basis"
if not os.path.exists(basedir):
    os.makedirs(basedir)

for obj in rot_objs:
    for atom in atoms:
        atomstr = "%03i" %(atom+1)
        print ""
        print "Loading object " + obj + " of atom " + atomstr
        values = f["dmft-" + iteration + "/ineq-" + atomstr + "/" + obj + "-full/value"][:]
        values = np.array(values,dtype=complex)
        if obj in ['gtau']:
            axis = f[".axes/taubin"][:]
            ax = "tau"
        elif obj in ['ftau']:
            axis = f[".axes/tauf"][:]
            ax = "tau"
        else:
            axis = f[".axes/iw"][:]
            ax = "iw"
        n = values.shape[-1]
        print values.shape, n
        print axis.shape
        print "Building basis transformation matrices"
        if toDIAG:
            Trafo = np.array(EigVec[atom])
            InvTrafo = np.conjugate(np.transpose(Trafo))
        else:
            InvTrafo = np.array(EigVec[atom])
            Trafo = np.conjugate(np.transpose(InvTrafo))
        print "Checking if transformation is unitary."
        try:
            if not np.allclose(np.dot(InvTrafo,Trafo),np.eye(2*Nbands)): raise RuntimeError("Full transformation \
               is not unitary! You may try to change np.linalg.eigh to np.linalg.eig at diagonalisation procedure")
        except RuntimeError:
            raise
            sys.exit()
        else:
            print "Transformation all clear."
            
        print "Rotating"
        values = values.reshape(2*Nbands,2*Nbands,n)
        print obj + ".shape", values.shape
        for i in range(0,n):
            values[:,:,i] = np.dot(InvTrafo,np.dot(values[:,:,i],Trafo))
        values = values.reshape(Nbands,2,Nbands,2,n)
        print obj + ".shape", values.shape

        print "Saving " + obj
        targetdir = basedir + "/atom-" + atomstr
        if not os.path.exists(targetdir):
            os.makedirs(targetdir)
        for b1 in range(0,Nbands):
            for s1 in range (0,2):
                for b2 in range(0,Nbands):
                    for s2 in range (0,2):
                        filename= targetdir + "/"  + obj + "_" + basis + "_" + str(b1+1) + str(s1+1) + str(b2+1) \
                                  + str(s2+1) +".dat"
                        targetf = open(filename, 'w')
                        targetf.write('#%14s' % ax)
                        targetf.write('%15s' % ("Re("+obj+")"))
                        targetf.write('%15s' % ("Im("+obj+")"))
                        targetf.write('\n')
                        data=np.column_stack((axis,np.real(values[b1,s1,b2,s2,:]),np.imag(values[b1,s1,b2,s2,:])))
                        np.savetxt(targetf,data,fmt='%15.10f')
                        targetf.close()

f.close()        

