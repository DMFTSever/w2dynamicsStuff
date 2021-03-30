#####################################################################################
# This script Rotates Hk+SOC into the eigenbasis of Hk+SOC. It also prints out some #
# checks if the rotation is valid. Make sure that thoose checks are ok!!!!          #
#####################################################################################

#imports
from __future__ import print_function, division, absolute_import
import h5py as hdf5
import numpy as np
import sys
import argparse
import warnings
import os

########################
# functions and errors #
########################

class InputError(Exception):
    """Exeption raised fo error in the input.
    """

class ShapeError(Exception):
    """Exeption raised fo error in the shaping.
    """

def warnwrapper(message):
    print("\n------- !!! WARNING !!! -------")
    warnings.warn(message)
    print("------- !!! ------- !!! -------\n")

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

def read_umatrix(filename, spin):
    '''Function to read a umatrix from a .dat file. 2 possible Formats
       1.: 1 1 1 1 value (spin==False)
       2.: 1u 1d 1u 1d value (spin==True)
        '''
    try:
        if type(filename) is not str:
            raise TypeError("1. argument (filename) must be of type str")
        if type(spin) is not bool:
            raise TypeError("2. argument (spin) must be of type bool:")
    except TypeError:
        raise
        sys.exit()

    f=open(filename, 'r')
    firstline = readaline(f).split()
    try:
        Nbands = int(firstline[0])
        if firstline[1].lower() != "bands":
            raise InputError("Expecting first non comment line to be of the form: # BANDS")
    except InputError:
        raise
        sys.exit()

    Nspin = int(spin)+1
    if spin:
        umatrix=np.zeros((Nbands,2,Nbands,2,Nbands,2,Nbands,2))
    else:
        umatrix=np.zeros((Nbands,Nbands,Nbands,Nbands))

    spindict = {'u':0,'d':1}
    spininfile = False
    try:
        for line in f:
            splitline = line.split()
            newsplitline = []
            for element in splitline[0:-1]:
                newsplitline.append(list(element))
            splitline = [y for x in newsplitline for y in x] + [splitline[-1]]
            index = np.zeros((Nspin*4,), dtype=int)
            for i in range(0,(len(splitline)-1)):
                element = splitline[i]
                if element in spindict:
                    index[i] = int(spindict[element])
                    spininfile = True
                else:
                    index[i] = int(element)-1
            umatrix[tuple(index)] = float(splitline[-1])
    except IndexError:
        raise InputError("Specifiefied spin dependency does not"\
                         " match ufile spin dependency")
        sys.exit()

    if spininfile!=spin:
        raise InputError("Specifiefied spin dependency does not"\
                         " match ufile spin dependency")
        sys.exit()

    return umatrix

def write_umatrix(filename, umatrix, spin):
    '''Funtion to write the umatrix either with spin dependency or without'''
    try:
        if type(filename) is not str:
            raise TypeError("1. argument (filenname) must be of type str.")
        if type(umatrix) is not np.ndarray:
            raise TypeError("2. argument (umatrix) must be of type numpy.ndarray.")
        if type(spin) is not bool:
            raise TypeError("3. argument (spin) must be of type bool.")
        if len(umatrix.shape)!=4 and len(umatrix.shape)!=8:
            raise ShapeError("Shape of umatrix is invalid")
        if len(umatrix.shape)==4 and spin:
            raise ShapeError("Spin independent umatrix detected although spin"\
                                 " was set to True.")
        if len(umatrix.shape)==8 and not spin:
            raise ShapeError("Spin dependent umatrix detected although spin"\
                                 " was set to False.")
        if spin and (umatrix.shape[1]!=2 or umatrix.shape[3]!=2 or umatrix.shape[5]!=2\
                     or umatrix.shape[7]!=2):
            raise ShapeError("Unexpected shape of spin dependent umatrix. Expected"\
                             " shape: (Nbands,2,Nbands,2,Nbands,2,Nbands,2)")
    except TypeError:
        raise
        sys.exit()
    except ShapeError:
        raise
        sys.exit()

    f = open(filename, 'w')
    f.write("# non zero elements of interaction matrix U_ijkl\n")
    f.write("%i BANDS\n"%umatrix.shape[0])
    spindict={0:'u',1:'d'}
    if not spin:
        for index, value in np.ndenumerate(umatrix):
            if value != 0.:
                f.write("%i %i %i %i  %016.16f\n"%(index[0]+1,index[1]+1,index[2]+1,\
                                                   index[3]+1, value))
    else:
        for index, value in np.ndenumerate(umatrix):
            if value != 0:
                f.write("%s %s %s %s  %016.16f\n"%(str(index[0]+1)+spindict[index[1]],\
                                                   str(index[2]+1)+spindict[index[3]],\
                                                   str(index[4]+1)+spindict[index[5]],\
                                                   str(index[6]+1)+spindict[index[7]],\
                                                   value))

########################
#       MAIN           #
########################

parser = argparse.ArgumentParser(description="""
                                             This script transfomrs U_ijkl according to
                                             the provided transformation. The output ufile
                                             is spin dependent by default. A pure band
                                             transformation can be activated with the
                                             --bandsonly option and does only work with
                                             non spin dependent hamiltonians.
                                             CURENTLY NO SUPPORT FOR COMPLEX UMATRIX AS
                                             INPUT
                                             """)
parser.add_argument('ufile', help='PATH to the U_ijkl file', type=str)
parser.add_argument('trafofile', help='PATH to the trafo file. Can be a .npy or a text '\
                                      'file.', type=str)
parser.add_argument('postfix', type=str,
                    help='Postfix for the filename of the transformed U_ijkl')
parser.add_argument('-s', '--spin', action='store_true', default=False, dest='spin',
                    help="If this option is set script expects transformation with spin"\
                         " dependency.")
parser.add_argument('-bo','--bandsonly', action='store_true', default=False,
                    dest='bandsonly',
                    help="If this option is set script expects a ufile without spin "\
                         "dependency. It will then only transform the bands. Only works "\
                         "if option -s is NOT set!.")

args      = parser.parse_args()
ufilename = args.ufile
trafofile = args.trafofile
spin      = args.spin
Nspins    = int(spin)+1
bandsonly = args.bandsonly # This bool is used throughout the script describing wether
                           # U_ijkl has spins or not!
Ntspins   = 2-int(bandsonly)
postfix   = args.postfix

try:
    if spin and bandsonly:
        raise RuntimeError("--bandsonly (-bo) only works without --spin (-s)!")
except RuntimeError:
    raise
    sys.exit()

#load umatrix
print("Loading Umatrix from file:"+ufilename)

#checking if ufile is spin dependent
ufile=open(ufilename,'r')
Nbands = int(readaline(ufile).split()[0])
line = readaline(ufile)
try:
    if (('u' in line) or ('d' in line)) and bandsonly:
        raise InputError("Detected spin dependency in input ufile."\
                         " Cannot perform bandsonly transformation!")
    elif ('u' in line) or ('d' in line):
        uspin = True
    else:
        uspin = False
except InputError:
    raise
    sys.exit()
ufile.close()

umatrix = read_umatrix(ufilename,spin=uspin)
print("umatrix.shape", umatrix.shape)

#blow up umatrix if it is not already spin dependent
if (not uspin) and (not bandsonly):
    print("Expanding U matrix to both spin channels")
    newumatrix = np.zeros(  (umatrix.shape[0],) + (2,) + (umatrix.shape[1],) + (2,)\
                          + (umatrix.shape[2],) + (2,) + (umatrix.shape[3],) + (2,)\
                          ,dtype=complex)
    newumatrix[:,0,:,0,:,0,:,0] = umatrix[:,:,:,:]
    newumatrix[:,1,:,0,:,1,:,0] = umatrix[:,:,:,:]
    newumatrix[:,0,:,1,:,0,:,1] = umatrix[:,:,:,:]
    newumatrix[:,1,:,1,:,1,:,1] = umatrix[:,:,:,:]
    umatrix = newumatrix

# reading basis transformation
print("Reading basis trafo from {}".format(trafofile))
fileending = os.path.splitext(trafofile)[1]
if fileending == '.npy':
    Trafo = np.load(trafofile)
else:
    try:
        Trafo = np.loadtxt(trafofile)
    except:
        raise RuntimeError("Could not load file: {}".format(trafofile))

if not spin and Trafo.shape[0]==2*Nbands:
    raise InputError("Looks like the supplied transformation has a spin dependency."\
                     " Use -s option!")

print("Trafo.shape", Trafo.shape)
if (not spin) and (not bandsonly):
    print("Expanding transformation to both spin channels.")
    temptrafo = np.zeros((Nbands,2,Nbands,2),dtype=complex)
    temptrafo[:,0,:,0] = Trafo
    temptrafo[:,1,:,1] = Trafo
    Trafo = temptrafo.reshape(2*Nbands,2*Nbands)
InvTrafo = Trafo.conjugate().transpose()

#transforming umatrix
print("Transforming umatrix")
umatrix = umatrix.reshape(Ntspins*Nbands,Ntspins*Nbands,Ntspins*Nbands,Ntspins*Nbands)
umatrix = np.einsum('ai,bj,ijkl,kc,ld',InvTrafo,InvTrafo,umatrix,Trafo,Trafo,dtype=complex)

if not bandsonly:
    umatrix = umatrix.reshape(Nbands,2,Nbands,2,Nbands,2,Nbands,2)

#Writting output
newfilename = os.path.splitext(ufilename)[0] + "_{}".format(postfix)
if np.allclose(umatrix.imag, np.zeros_like(umatrix.imag)):
    name = newfilename + ".dat"
    print("Writting U_ijkl to file " + name)
    write_umatrix(name, umatrix.real, spin=(not bandsonly))

else:
    namereal = newfilename + "_re.dat"
    nameimag = newfilename + "_im.dat"
    print("Writting Re(U_ijkl) to file " + namereal)
    write_umatrix(namereal, umatrix.real, spin=(not bandsonly))
    print("Writting Im(U_ijkl) to file " + nameimag)
    write_umatrix(nameimag, umatrix.imag, spin=(not bandsonly))
