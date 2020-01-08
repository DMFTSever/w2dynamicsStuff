# imports
import h5py as hdf5
import numpy as np
import sys
### here we need to add the path to the w2dynamics main directory
auxdir="/dss/dsshome1/0D/di76rir/lib/MyPythonScripts/W2dynamics/General"
sys.path.insert(0,auxdir)
### in order to import the functions available in input.py
import readwrite as rw
import argparse

parser = argparse.ArgumentParser(description="This script adds SOC to an existing t2g Hamiltonian. This will quadruple\
					      its size.")

parser.add_argument('hkfile', help='PATH to an existing file containing the hamiltonian', type=str)
parser.add_argument('Natoms', help='Number of atoms of the hamiltonian / brillouin zone', type=int)
parser.add_argument('epsilon', help='Epsilon Parameter of SOC', type=float)

args = parser.parse_args()
hkfile = args.hkfile
epsilon = args.epsilon
Na = args.Natoms

#setting the epsilon parameter for the SOC
print "epsilonSOC = ", epsilon

### load hamiltonian
print "Loading hamiltonian from file: ", hkfile
hk, kpoints = rw.read_hk_wannier(hkfile,spin=False)

print "hk.shape", hk.shape
print "kpoints.shape", kpoints.shape

### get number of k-points
Nk=kpoints.shape[0]
print "Nk", Nk

### number of d-orbitals
Nd=hk.shape[1]
print "Nd", Nd

#Setting atoms and orbitals per atom
Ndpa = Nd/Na       #Number of d-orbitals per atom
Ns = 2             #Number of spins
print "Na", Na 
print "Ndpa", Ndpa 
print "Ns", Ns

#Blowing the Hamiltonian to up for spin inclusion
hk = hk.transpose(0,1,3,2,4)
print "hk.shape", hk.shape
print "Blowing up Hamiltonian"
hk = np.append(np.append(hk,np.zeros_like(hk),axis=4),np.append(hk,np.zeros_like(hk),axis=4)[:,:,:,:,::-1],axis=3)
print "hk.shape", hk.shape

#build hSOC from pauli matricies
print "Building hSOC"
hSOC = np.zeros_like(hk[1])
sigma1 = np.array([[0,1],[1,0]])
sigma2 = np.array([[0,-1j],[1j,0]])
sigma3 = np.array([[1,0],[0,-1]])

#order xz, xz, yz
for i in range(0,Na):
		hSOC[0+3*i,1+3*i] = -1j*sigma1
		hSOC[0+3*i,2+3*i] = 1j*sigma2
		hSOC[1+3*i,0+3*i] = 1j*sigma1
		hSOC[1+3*i,2+3*i] = -1j*sigma3
		hSOC[2+3*i,0+3*i] = -1j*sigma2
		hSOC[2+3*i,1+3*i] = 1j*sigma3
hSOC = epsilon * hSOC
print "hSOC.shape", hSOC.shape

#build hkSOC
print "Building hkSOC"
hkSOC = np.zeros_like(hk)
for i in range(0,hkSOC.shape[0]):
	hkSOC[i] = hSOC
print "hkSOC.shape", hkSOC.shape

#Adding together hk and hkSOC
print "Adding hkSOC to hk"
hk += hkSOC

#Building hkmean as check
hkmean = 1./Nk*np.sum(hk, axis=0)
hkmean = hkmean.transpose(0,2,1,3)
hkmean = hkmean.reshape(Ns*Nd,Ns*Nd)

#writting hkmean of the first atom
print "hkmean including SOC for the first atom" 
for i in range(0,6):
	print '%+3.3f%+3.3fi' % (np.real(hkmean[i,0]), np.imag(hkmean[i,0])), '%+3.3f%+3.3fi' % (np.real(hkmean[i,1]), np.imag(hkmean[i,1])), '%+3.3f%+3.3fi' % (np.real(hkmean[i,2]), np.imag(hkmean[i,2])), '%+3.3f%+3.3fi' % (np.real(hkmean[i,3]), np.imag(hkmean[i,3])), '%+3.3f%+3.3fi' % (np.real(hkmean[i,4]), np.imag(hkmean[i,4])), '%+3.3f%+3.3fi' % (np.real(hkmean[i,5]), np.imag(hkmean[i,5]))

#Shaping it such that it can be written in wanner90 format
print "Shaping Hamiltonian for passing to write function"
print "hk.shape", hk.shape
hk = hk.transpose(0,1,3,2,4)
print "hk.shape", hk.shape

#Writting output
hkfilenew = hkfile[:-4] + "_SOC.dat"
print "Writting Hamiltonian to file: " + hkfilenew
rw.write_hk_wannier(hkfilenew,hk,kpoints)

