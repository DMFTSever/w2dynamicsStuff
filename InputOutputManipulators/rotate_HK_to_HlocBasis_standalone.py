#####################################################################################
# This script Rotates Hk+SOC into the eigenbasis of Hk+SOC. It also prints out some #
# checks if the rotation is valid. Make sure that thoose checks are ok!!!!          #
#####################################################################################

from __future__ import print_function, division, absolute_import
import h5py as hdf5
import numpy as np
import os
import sys
import argparse
import warnings

########################
# functions and errors #
########################

class InputError(Exception):
    """Exeption raised fo error in the input.
    """

def warnwrapper(message):
    print("\n------- !!! WARNING !!! -------")
    warnings.warn(message)
    print("------- !!! ------- !!! -------\n")

def warnwrapper_soft(message):
    print("\n----------- warning -----------")
    print(message)
    print("-------------------------------\n")

def readaline(readfile, comments='#'):
    '''
       This is a function to reading the next line which is not a comment from
       a file ignoring everthing after the comment symbol
    '''
    while True:
        line = readfile.readline()
        if line[0] == comments:
            pass
        else:
            break
    commentpos = 0
    for letter in line:
        if letter == comments:
            return line[0:commentpos] + '\n'
        commentpos += 1
    return line

def read_hk_wannier(filename,spin=False):
    '''
       Function to read hk from a file with name 'filename'. spin dependency is False
       by default. If spin is set to True it assumes wannier convention (orbital index
       running fastest, spin index slowest). It returns the hk of form of an numpy array
       of form (number of kpoints, number of orbitals, number of spins, number of
       orbitals, number of spins), the kpoints in form of a numpy array of form (number
       of kpoints, x, y, z)\nreturn kpoints, hk
    '''
    f         = open(filename, 'r')
    firstline = readaline(f)
    nk        = int(firstline.split()[0]) #Number of k points of hk
    nbands    = int(firstline.split()[1]) #Number of orbitals (inkluding spin) of hk
    hk        = np.zeros((nk,nbands,nbands),dtype=complex)
    kpoints   = np.zeros((nk,3))

    for i in range(0,nk):
        splitline = readaline(f).split()
        kpoints[i,0] = float(splitline[0])
        kpoints[i,1] = float(splitline[1])
        kpoints[i,2] = float(splitline[2])
        for j in range(0,nbands):
            splitline = readaline(f).split()
            for k in range(0,nbands):
                hk[i,j,k] = complex(float(splitline[2*k]),float(splitline[2*k+1]))
    if spin:
        hk = hk.reshape(nk, 2, int(nbands/2), 2, int(nbands/2))
        hk = hk.transpose(0,2,1,4,3)
    else:
        hk = hk.reshape(nk, 1, nbands, 1, nbands)
        hk = hk.transpose(0,2,1,4,3)
    f.close()
    return hk, kpoints

def write_hk_wannier(filename,hk,kpoints):
    '''
       Function to write hk to a file named 'filename' in wannier format (spin indices
       running slower than orbital indices). Assumes hk to be of the form
       (nk,nbands,spins,nbands,spins) and kpoint(nk,3). nd is the number of orbitals
       (excluding spin), nk the number of kpoints.
    '''
    f      = open(filename, 'w')
    nk     = hk.shape[0]
    nbands = hk.shape[1]
    nspin  = hk.shape[2]
    hk     = hk.transpose(0,2,1,4,3)
    hk     = hk.reshape(nk,nspin*nbands,nspin*nbands)

    ### write first line: number of k-points, number wannier orbitals, number of bands
    f.write(str(nk)+" "+str(nspin*nbands)+" "+str(nspin*nbands)+"\n")

    for k in range(0,nk):
        f.write(str(kpoints[k,0])+"  "+str(kpoints[k,1])+"  "+str(kpoints[k,2])+" \n")
        for i in range(0,nspin*nbands):
            for j in range(0,nspin*nbands):
                f.write('%+010.8f' % np.real(hk[k,i,j]))
                f.write("  ")
                f.write('%+10.8f' % np.imag(hk[k,i,j]))
                f.write("  ")
            f.write("\n")

###########################
#          MAIN           #
###########################

parser = argparse.ArgumentParser(description="""
                                             This script transforms a k dependent
                                             Hamiltonian to its local eigenbasis
                                             """)
parser.add_argument('hkfile',   help='PATH to the k dependent Hamiltonian', type=str)
parser.add_argument('Natoms',   help='Number of atoms', type=int)
parser.add_argument('Ndbands',  help='Number of d-bands per atom', type=int)
parser.add_argument('--Npbands', help="""Number of p-bands in total. Assumed to be
                                         appended to last atom""",\
                    type=int, default=0)
parser.add_argument('-s', '--spin', action='store_true', default=False,\
                    help="""
                         If this option is set script expects a Hamiltonian with spin
                         dependency, spin indices varying slowest (Wannier90 Convention).
                         """,
                    dest='spin')
parser.add_argument('--all-atoms-equiv', action='store_true', default=False, dest='equiv',\
                    help="All atoms are assumed to be equivalent")

args     = parser.parse_args()
filename = args.hkfile
spin     = args.spin
Nspins   = int(spin)+1
Natoms   = args.Natoms
Ndbands  = args.Ndbands
Npbands  = args.Npbands
equiv    = args.equiv

if Npbands != 0:
    warnwrapper_soft("You have activated p-bands. All p-bands are assumed"\
                     " to be appeneded to the last atom!")

### load hamiltonian
print("Loading hamiltonian on hk and kpoints on kpoints.")
hkfile      = filename
hk, kpoints = read_hk_wannier(filename, spin=spin)
print("hk.shape", hk.shape)
print("kpoints.shape", kpoints.shape)

### get number of k-points
Nk=kpoints.shape[0]
print("Number of k points: ", Nk)

try:
    if hk.shape[1] != Natoms*Ndbands+Npbands:
        raise InputError("Hk dimension do not match with Natom, Ndbands, Npbands, Nspins!")
except InputError:
    raise
    sys.exit()

#Building hkmean
print("Building hkmean i.e. averaging over all k points")
hkmean = 1./Nk * np.sum(hk, axis=0)

#reshape hkmean
hkmean = hkmean.reshape(Nspins*(Ndbands*Natoms+Npbands),Nspins*(Ndbands*Natoms+Npbands))

#Creating array with mean hamiltonians per atom
print("Building hkmean_pa (hkmean per atom, d-bands only)")
hkmean_pa = np.zeros((Natoms,Nspins*Ndbands,Nspins*Ndbands),dtype=complex)
for i in range(0,Natoms):
    hkmean_pa[i] = hkmean[i*Nspins*Ndbands:(i+1)*Nspins*Ndbands,\
                          i*Nspins*Ndbands:(i+1)*Nspins*Ndbands]

for a in range(0,Natoms):
    print("hkmean_pa[" + str(a) + "]")
    if spin:
        for i in range(0,Nspins*Ndbands):
            print('%+05.5f%+05.5fi'%(np.real(hkmean_pa[a,i,0]),np.imag(hkmean_pa[a,i,0])),\
                  '%+05.5f%+05.5fi'%(np.real(hkmean_pa[a,i,1]),np.imag(hkmean_pa[a,i,1])),\
                  '%+05.5f%+05.5fi'%(np.real(hkmean_pa[a,i,2]),np.imag(hkmean_pa[a,i,2])),\
                  '%+05.5f%+05.5fi'%(np.real(hkmean_pa[a,i,3]),np.imag(hkmean_pa[a,i,3])),\
                  '%+05.5f%+05.5fi'%(np.real(hkmean_pa[a,i,4]),np.imag(hkmean_pa[a,i,4])),\
                  '%+05.5f%+05.5fi'%(np.real(hkmean_pa[a,i,5]),np.imag(hkmean_pa[a,i,5])))
    else:
        for i in range(0,Nspins*Ndbands):
            print('%+05.5f%+05.5fi'%(np.real(hkmean_pa[a,i,0]),np.imag(hkmean_pa[a,i,0])),\
                  '%+05.5f%+05.5fi'%(np.real(hkmean_pa[a,i,1]),np.imag(hkmean_pa[a,i,1])),\
                  '%+05.5f%+05.5fi'%(np.real(hkmean_pa[a,i,2]),np.imag(hkmean_pa[a,i,2])))

#Solving the Eigenvalue problem
print("Solving the Eigenvalue problem for hkmean_pa")
EigVal = np.zeros((Natoms,Nspins*Ndbands),dtype=complex)
EigVec = np.zeros_like(hkmean_pa)
if equiv:
    print("All atoms assumed ot be equivalent. Only diagonalizing Atom 1")
    if np.allclose(hkmean_pa[0]-np.diag(np.diagonal(hkmean_pa[0])),\
                   np.zeros_like(hkmean_pa[0])):
        print("Atom 1 already diagonal. No diagonalisation necessary")
        print("Since all atoms are assumed equivalent nothing to do. Aborting...")
        sys.exit()
    else:
        EigVal[0], EigVec[0] = np.linalg.eigh(hkmean_pa[0])
        print("Eigvalues of atom 1", EigVal[0])
    print("Copying trafo of atom 1 to all other atoms")
    for i in range(1,Natoms):
        EigVal[i], EigVec[i] = EigVal[0], EigVec[0]
else:
    for i in range(0,Natoms):
        if np.allclose(hkmean_pa[i]-np.diag(np.diagonal(hkmean_pa[i])),\
                       np.zeros_like(hkmean_pa[i])):
            print("Atom "+str(i+1) + " already diagonal. No diagonalisation necessary")
            EigVal[i] = np.diagonal(hkmean_pa[i])
            EigVec[i] = np.identity(hkmean_pa[i].shape[0])
        else:
            EigVal[i], EigVec[i] = np.linalg.eigh(hkmean_pa[i])
            print("Eigvalues of atom " + str(i), EigVal[i])

print("Checking if transformations of each atom are unitary.")
eye = np.eye(Nspins*Ndbands)
for i in range(0,len(EigVal)):
    print("Checking transformation of atom %i" %(i+1))
    try:
        if not np.allclose(np.dot(EigVec[i].conj().T,EigVec[i]),eye):
            raise RuntimeError("""Transformation of atom %i is not unitary! You may try to
                                  change np.linalg.eigh to np.linalg.eig at diagonalisation
                                  procedure"%(i+1)""")
    except RuntimeError:
        raise
        sys.exit()
    else:
        print("Transformation of atom %i all clear."%(i+1))

#Creating the Basis Transformation
print("Building full basis transformation matrix")
Trafo = np.zeros_like(hkmean)
for i in range(0,Natoms):
    Trafo[i*Nspins*Ndbands:(i+1)*Nspins*Ndbands,\
          i*Nspins*Ndbands:(i+1)*Nspins*Ndbands] = EigVec[i]
if Ndbands != 0:
    Trafo[Natoms*Nspins*Ndbands:,Natoms*Nspins*Ndbands:] = np.eye(Npbands)
InvTrafo = Trafo.conjugate().transpose()

print("Checking if full transformation is unitary.")
try:
    if not np.allclose(np.dot(InvTrafo,Trafo),np.eye(Nspins*Ndbands*Natoms)):
        raise RuntimeError("Full transformation is not unitary! You may try to change"\
                           " np.linalg.eigh to np.linalg.eig at diagonalisation procedure")
except RuntimeError:
    raise
    sys.exit()
else:
    print("Full transformation all clear.")

print("Writting Trafos (Q): Hdiag = Q^{-1}HQ")
for i in range(0,Natoms):
    np.savetxt("Trafo_at"+str(i+1)+"_toHlocBasis.dat",EigVec[i])
    np.save("Trafo_at"+str(i+1)+"_toHlocBasis",EigVec[i])
np.savetxt("Trafo_Full_toHlocBasis.dat",Trafo)
np.save("Trafo_Full_toHlocBasis",Trafo)

#Rotating Hamiltonian
print("Transforming the Hamiltonian")
hk = hk.reshape(Nk,Nspins*(Ndbands*Natoms+Npbands),Nspins*(Ndbands*Natoms+Npbands))
for i in range(0,Nk):
    hk[i] = np.dot(InvTrafo,np.dot(hk[i],Trafo))

#Checking if rotated hk_mean is diagonal with entries equal to the calculated eigenvalues
print("Checking if transformed hk_mean per atom is diagonal with entries equal"\
      " to the calculated eigenvalues")
oldhkmean = hkmean
hkmean = 1./Nk * np.sum(hk, axis=0)
warn=False
for i in range(0,Natoms):
    try:
        temphkmean = hkmean[i*Nspins*Ndbands:(i+1)*Nspins*Ndbands,\
                            i*Nspins*Ndbands:(i+1)*Nspins*Ndbands]
        if not np.allclose(temphkmean,np.diag(EigVal[i])):
            if equiv:
                warn       = True
                diag       = np.diag(temphkmean)
                offdiag    = temphkmean-np.diag(diag)
                maxdiag    = np.amax(abs(diag))
                diffdiag    = np.amax(diag) - np.amin(diag)
                maxoffdiag = np.amax(abs(offdiag))
                warnwrapper("Transformed hamiltonian is locally not diagonal for"\
                            " atom {}. Continuing since --all-atoms-equiv is"\
                            " set. Check results!\n"\
                            "Largest diagonal entry (absolute): {}\n"\
                            "Largest offdiagonal entry (absolute): {}\n"\
                            "Largest difference of diagonal entries: {}\n"\
                            .format(i+1, maxdiag, maxoffdiag, diffdiag))
            else:
                raise RuntimeError("Transformed hamiltonian is locally not diagonal for"\
                                   " atom %i. Numeric errors to big!"%(i+1))
    except RuntimeError:
        raise
        sys.exit()
    else:
        if warn:
            continue
        else:
            print("Transformed local hamiltonian of atom %i all clear."%(i+1))

hk = hk.reshape(Nk,Ndbands*Natoms+Npbands,Nspins,Ndbands*Natoms+Npbands,Nspins)

#Writting output
newfilename = os.path.splitext(filename)[0]+"_inHlocBasis.dat"
print("Writting Hamiltonian to file "+newfilename)
write_hk_wannier(newfilename,hk,kpoints)
