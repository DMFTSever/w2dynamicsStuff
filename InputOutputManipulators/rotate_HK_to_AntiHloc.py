###########################################################################################################################################################
#This script Rotates Hk+SOC into the eigenbasis of Hk+SOC. It also prints out some checks if the rotation is valid. Make sure that thoose checks are ok!!!!
############################################################################################################################################################

#imports
import h5py as hdf5
import numpy as np
import sys
import argparse
### here we need to add the path to the personal myfunc directory (location of some read functions)
auxdir="/dss/dsshome1/0D/di76rir/lib/MyPythonScripts/W2dynamics/General"
sys.path.insert(0,auxdir)
### in order to import the functions available in input.py and interatction.py and readwrite and custom errors
import readwrite as rw
import custom_errors as err

parser = argparse.ArgumentParser(description="This script transforms a k dependent Hamiltonian to its local eigenbasis")
parser.add_argument('hkfile', help='PATH to the k dependent Hamiltonian', type=str)
parser.add_argument('Natoms', help='Number of atoms', type=int)
parser.add_argument('Nbands', help='Number of bands per atom', type=int)
parser.add_argument('-s', '--spin', action='store_true', default=False,\
	            help="If this option is set script expects a Hamiltonian with spin dependency, spin indices varying slowest \
			 (Wannier90 Convention).", dest='spin')
args = parser.parse_args()
filename = args.hkfile
spin = args.spin
Nspins = int(spin)+1
Natoms = args.Natoms
Nbands = args.Nbands

### load hamiltonian
print "Loading hamiltonian on hk and kpoints on kpoints."
hkfile = file(filename)
hk, kpoints = rw.read_hk_wannier(filename, spin=spin)

print "hk.shape", hk.shape
print "kpoints.shape", kpoints.shape

### get number of k-points
Nk=kpoints.shape[0]
print "Number of k points: ", Nk

try:
	if hk.shape[1] < Natoms*Nbands: raise err.InputError("Input hkfile does not have spin entries! Please don't use -s.")
	if hk.shape[1] > Natoms*Nbands: raise err.InputError("Input hkfile does have spin entries! Please use -s.")
except err.InputError:
	raise
	sys.exit()

#Building hkmean
print "Building hkmean i.e. averaging over all k points"
hkmean = 1./Nk * np.sum(hk, axis=0)
print "hkmean.shape", hkmean.shape

#reshape hkmean
hkmean = hkmean.reshape(Nspins*Nbands*Natoms,Nspins*Nbands*Natoms)
print "hkmean.shape", hkmean.shape

#Creating array with mean hamiltonians per atom
print "Building hkmean_pa (hkmean per atom)"
hkmean_pa = np.zeros((Natoms,Nspins*Nbands,Nspins*Nbands),dtype=complex)
for i in range(0,Natoms):
	hkmean_pa[i] = hkmean[i*Nspins*Nbands:(i+1)*Nspins*Nbands,i*Nspins*Nbands:(i+1)*Nspins*Nbands]

print "hkmean_pa.shape", hkmean_pa.shape

for a in range(0,Natoms):
	print "hkmean_pa[" + str(a) + "]"
	if spin:
		for i in range(0,Nspins*Nbands):
			print '%+05.5f%+05.5fi' %(np.real(hkmean_pa[a,i,0]), np.imag(hkmean_pa[a,i,0])), '%+05.5f%+05.5fi' %(np.real(hkmean_pa[a,i,1]), np.imag(hkmean_pa[a,i,1])), '%+05.5f%+05.5fi' %(np.real(hkmean_pa[a,i,2]), np.imag(hkmean_pa[a,i,2])), '%+05.5f%+05.5fi' %(np.real(hkmean_pa[a,i,3]), np.imag(hkmean_pa[a,i,3])), '%+05.5f%+05.5fi' %(np.real(hkmean_pa[a,i,4]), np.imag(hkmean_pa[a,i,4])), '%+05.5f%+05.5fi' %(np.real(hkmean_pa[a,i,5]), np.imag(hkmean_pa[a,i,5])) 
	else:
		for i in range(0,Nspins*Nbands):
			print '%+05.5f%+05.5fi' %(np.real(hkmean_pa[a,i,0]), np.imag(hkmean_pa[a,i,0])), '%+05.5f%+05.5fi' %(np.real(hkmean_pa[a,i,1]), np.imag(hkmean_pa[a,i,1])), '%+05.5f%+05.5fi' %(np.real(hkmean_pa[a,i,2]), np.imag(hkmean_pa[a,i,2]))

#Solving the Eigenvalue problem
print "Solving the Eigenvalue problem for hkmean_pa"
EigVal = np.zeros((Natoms,Nspins*Nbands),dtype=complex)
EigVec = np.zeros_like(hkmean_pa)
for i in range(0,Natoms):
	EigVal[i], EigVec[i] = np.linalg.eigh(hkmean_pa[i])
	print "Eigvalues of atom " + str(i), EigVal[i]

print "Checking if transformations of each atom are unitary."
eye = np.eye(Nspins*Nbands)
for i in range(0,len(EigVal)):
	print "Checking transformation of atom %i" %(i+1)
	try:
		if not np.allclose(np.dot(EigVec[i].conj().T,EigVec[i]),eye): raise RuntimeError("Transformation of atom %i is not unitary! You may try to change np.linalg.eigh to np.linalg.eig at diagonalisation procedure"%(i+1))
	except RuntimeError:
		raise
		sys.exit()
	else:
		print "Transformation of atom %i all clear."%(i+1)

#Creating the Basis Transformation
print "Building full basis transformation matrix"
Trafo = np.zeros_like(hkmean)
for i in range(0,Natoms):
	Trafo[i*Nspins*Nbands:(i+1)*Nspins*Nbands,i*Nspins*Nbands:(i+1)*Nspins*Nbands] = EigVec[i] 
InvTrafo = Trafo.conj().T

print "Checking if full transformation is unitary."
try:
	if not np.allclose(np.dot(InvTrafo,Trafo),np.eye(Nspins*Nbands*Natoms)): raise RuntimeError("Full transformation is not unitary! You may try to change np.linalg.eigh to np.linalg.eig at diagonalisation procedure")
except RuntimeError:
	raise
	sys.exit()
else:
	print "Full transformation all clear."

#Rotating Hamiltonian
print "Transforming the Hamiltonian"
hk = hk.reshape(Nk,Nspins*Nbands*Natoms,Nspins*Nbands*Natoms)
print "hk.shape", hk.shape
for i in range(0,Nk):
	hk[i] = np.dot(Trafo,np.dot(hk[i],InvTrafo))

#Checking if rotated hk_mean is diagonal with entries equal to the calculated eigenvalues
#print "Checking if transformed hk_mean per atom is diagonal with entries equal to the calculated eigenvalues"
#oldhkmean = hkmean
#hkmean = 1./Nk * np.sum(hk, axis=0)
#for i in range(0,Natoms):
#	try:
#		temphkmean = hkmean[i*Nspins*Nbands:(i+1)*Nspins*Nbands,i*Nspins*Nbands:(i+1)*Nspins*Nbands] 
#		if not np.allclose(temphkmean,np.diag(EigVal[i])): raise RuntimeError("Transformed hamiltonian is locally not diagonal for atom %i. Numeric errors to big!"%(i+1))
#	except RuntimeError:
#		raise
#		sys.exit()
#	else:
#		print "Transformed local hamiltonian of atom %i all clear."%(i+1)

hk = hk.reshape(Nk,Nbands*Natoms,Nspins,Nbands*Natoms,Nspins)
print "hk.shape", hk.shape

#Writting output
newfilename = filename[:-4]+"_rot.dat"
print "Writting Hamiltonian to file "+newfilename
rw.write_hk_wannier(newfilename,hk,kpoints)
