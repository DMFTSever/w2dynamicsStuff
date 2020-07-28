###########################################################################################################################################################
#This script Rotates Hk+SOC into the eigenbasis of Hk+SOC. It also prints out some checks if the rotation is valid. Make sure that thoose checks are ok!!!!
############################################################################################################################################################

#imports
from __future__ import print_function, division, absolute_import
import h5py as hdf5
import numpy as np
import sys
import argparse
### here we need to add the path to the personal myfunc directory (location of some read functions)
sys.path.insert(0,sys.argv[0].replace('/InputOutputManipulators/rotate_U_ijkl_to_Hloc_basis.py','/General'))
### in order to import the functions available in input.py and interatction.py and readwrite and custom errors
import readwrite as rw
import custom_errors as err

parser = argparse.ArgumentParser(description="This script transfomrs U_ijkl of one atom to the local eigenbasis of the provided hamiltonian\
					      .The output ufile is spin dependent by default. A pure band transformation can be activated\
					      with the --bandsonly option and does only work with non spin dependent hamiltonians.")
parser.add_argument('hkfile', help='PATH to the k dependent Hamiltonian', type=str)
parser.add_argument('Natoms', help='Number of atoms', type=int)
parser.add_argument('Nbands', help='Number of bands per atom', type=int)
parser.add_argument('ufile', help='PATH to the U_ijkl file', type=str)
parser.add_argument('atom', help='The transformation matrix of this atom will be used to transform U_ijkl', type=int)
parser.add_argument('postfix', help='Postfix for the filename of the transformation matrix', type=str)
parser.add_argument('-s', '--spin', action='store_true', default=False,\
	            help="If this option is set script expects a Hamiltonian with spin dependency, spin indices varying slowest \
			 (Wannier90 Convention).", dest='spin')
parser.add_argument('-bo','--bandsonly', action='store_true', default=False,\
	            help="If this option is set script expects a ufile without spin dependency. It will then only transform the bands.\
			 Only works if option -s is NOT set!.", dest='bandsonly')
parser.add_argument('--corrector', action='store_true', default=False,\
	            help="If this option is set the umatrix will be corrected before it is transformed. Entries corresponding to the same band and\
			 spin in the first 4 or last 4 indices are set to zero. Has no influence if -bo (--bandsonly) is set.", dest='corrector')
args = parser.parse_args()
filename = args.hkfile
ufilename = args.ufile
spin = args.spin
Nspins = int(spin)+1
Natoms = args.Natoms
Nbands = args.Nbands
atom = args.atom
bandsonly = args.bandsonly #This bool is used throughout the script describing wether U_ijkl has spins or not!
			   #Usage may be counterintuitive. Name is choosen such that users only use it if they really want to!
Ntspins = 2-int(bandsonly)
corrector = args.corrector
postfix = args.postfix

try:
	if atom > Natoms: raise RuntimeError("Number supplied for atom must be smaller than number supplied for Natoms")
	if spin and bandsonly: raise RuntimeError("--bandsonly (-bo) only works without --spin (-s)!")
except RuntimeError:
	raise
	sys.exit()

### load hamiltonian
print("Loading hamiltonian on hk and kpoints on kpoints from file:"+filename)
hk, kpoints = rw.read_hk_wannier(filename,spin=spin)

print("hk.shape", hk.shape)
print("kpoints.shape", kpoints.shape)

### get number of k-points
Nk=kpoints.shape[0]
print("Number of k points: ", Nk)

try:
	if hk.shape[1] < Natoms*Nbands: raise err.InputError("Input hkfile does not have spin entries! Please use -s False.")
	if hk.shape[1] > Natoms*Nbands: raise err.InputError("Input hkfile does have spin entries! Please use -s True.")
except err.InputError:
	raise
	sys.exit()

#load umatrix
#checking if ufile is spin dependent
print("Loading Umatrix from file:"+ufilename)
ufile=open(ufilename,'r')
rw.readaline(ufile)
line = rw.readaline(ufile)

try:
	if (('u' in line) or ('d' in line)) and bandsonly: 
		raise err.InputError("Detected spin dependency in input ufile.\
				      Cannot perform bandsonly transformation!") 
	elif ('u' in line) or ('d' in line): 
		uspin = True
	else:
		uspin = False
except err.InputError:
	raise
	sys.exit()

ufile.close()

umatrix = rw.read_umatrix(ufilename,spin=uspin)
print("umatrix.shape", umatrix.shape)

#blow up umatrix if it is not already spin dependent
if (not uspin) and (not bandsonly):
	print("Blowing up U matrix")
	newumatrix = np.zeros((umatrix.shape[0],) + (2,) + (umatrix.shape[1],) + (2,) + (umatrix.shape[2],) + (2,) + (umatrix.shape[3],) + (2,),dtype=complex)
	newumatrix[:,0,:,0,:,0,:,0] = umatrix[:,:,:,:]
	newumatrix[:,1,:,0,:,1,:,0] = umatrix[:,:,:,:]
	newumatrix[:,0,:,1,:,0,:,1] = umatrix[:,:,:,:]
	newumatrix[:,1,:,1,:,1,:,1] = umatrix[:,:,:,:]
	umatrix = newumatrix
	if corrector:
		for index, x in np.ndenumerate(umatrix):
			if (index[0:2] == index[2:4]) or (index[4:6] == index[6:8]): 
				umatrix[index] = 0.
	print("umatrix.shape", umatrix.shape)

#Building hkmean
print("Building hkmean i.e. averaging over all k points")
hkmean = 1./Nk * np.sum(hk, axis=0)
print("hkmean.shape", hkmean.shape)

#reshape hkmean
hkmean = hkmean.reshape(Nspins*Nbands*Natoms,Nspins*Nbands*Natoms)
print("hkmean.shape", hkmean.shape)

#Creating array with mean hamiltonians per atom
print("Building hkmean_pa (hkmean per atom)")
hkmean_pa = np.zeros((Natoms,Nspins*Nbands,Nspins*Nbands),dtype=complex)
for i in range(0,Natoms):
	hkmean_pa[i] = hkmean[i*Nspins*Nbands:(i+1)*Nspins*Nbands,i*Nspins*Nbands:(i+1)*Nspins*Nbands]

print("hkmean_pa.shape", hkmean_pa.shape)

#Solving the Eigenvalue problem
if np.allclose(hkmean_pa[atom-1]-np.diag(np.diagonal(hkmean_pa[atom-1])), np.zeros_like(hkmean_pa[atom-1])):
    print("Atom "+str(atom) + " already diagonal. No diagonalisation necessary for this atom")
    EigVal = np.diagonal(hkmean_pa[atom-1])
    EigVec = np.identity(hkmean_pa[atom-1].shape[0], dtype=complex)
else:
    print("Solving the Eigenvalue problem of hkmean_pa for atom %i"%atom)
    EigVal, EigVec = np.linalg.eigh(hkmean_pa[atom-1])
    print("Eigvalues of atom " + str(atom), EigVal)

#Creating the Basis Transformation
print("Building basis transformation matrix of atom %i"%atom)
Trafo = EigVec
if (not spin) and (not bandsonly):
	print("Expanding transformation to both spin channels.")
	temptrafo = np.zeros((Nbands,2,Nbands,2),dtype=complex)
	temptrafo[:,0,:,0] = Trafo
	temptrafo[:,1,:,1] = Trafo
	Trafo = temptrafo.reshape(2*Nbands,2*Nbands)
InvTrafo = Trafo.conjugate().transpose()
print("Trafo.shape", Trafo.shape)

print("Checking if transformation is unitary.")
try:
	if not np.allclose(np.dot(InvTrafo,Trafo),np.eye(Ntspins*Nbands)): 
		raise RuntimeError("Transformation is not unitary! You may try to change np.linalg.eigh\
				    to np.linalg.eig at diagonalisation procedure")
except RuntimeError:
	raise
	sys.exit()
else:
	print("Transformation all clear.")

print("Writting Trafo (Q): Hdiag = Q^{-1}HQ")
np.savetxt("UTrafo_at"+str(atom)+"_"+postfix+".dat",Trafo)

#transforming umatrix
print("Transforming umatrix with transformation matrix of atom %i"%atom)
umatrix = umatrix.reshape(Ntspins*Nbands,Ntspins*Nbands,Ntspins*Nbands,Ntspins*Nbands)
umatrix = np.einsum('ai,bj,ijkl,kc,ld',InvTrafo,InvTrafo,umatrix,Trafo,Trafo,dtype=complex)
print("umatrix.shape", umatrix.shape)

if not bandsonly:
	umatrix = umatrix.reshape(Nbands,2,Nbands,2,Nbands,2,Nbands,2)

#Writting output
newfilename1 = ufilename[:-4] + "_rot_re.dat"
newfilename2 = ufilename[:-4] + "_rot_im.dat"
print("Writting Re(U_ijkl) to file " + newfilename1)
rw.write_umatrix(newfilename1, umatrix.real, spin=(not bandsonly))
print("Writting Im(U_ijkl) to file " + newfilename2)
rw.write_umatrix(newfilename2, umatrix.imag, spin=(not bandsonly))
