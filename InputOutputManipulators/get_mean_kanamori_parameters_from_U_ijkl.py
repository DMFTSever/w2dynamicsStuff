###################################################################################################
# This script extracts the mean kanamori parameter of spin or non spin dependent U_ijkl (umatrix) #
###################################################################################################

# imports
import h5py as hdf5
import numpy as np
import sys
import argparse
### here we need to add the path to the personal myfunc directory (location of some read functions)
sys.path.insert(0,sys.argv[0].replace('/InputOutputManipulators/get_mean_kanamori_parameters_from_Uijkl.py','/General'))
### in order to import the functions available in input.py and interatction.py and readwrite and custom errors
import readwrite as rw
import custom_errors as err

parser = argparse.ArgumentParser(description="This script extracts the mean kanamori Parameter U,J,V from a given U_ijkl file.\
					      The file can be spin or not spin dependent, however this needs to be specified by\
					      as a second argument")
parser.add_argument('ufile', help='PATH to the U_ijkl file', type=str)
parser.add_argument('--nospin', action='store_false', default=True, help='This specifies that the ufile is not spin dependent.', dest='spin')

args = parser.parse_args()
filename = args.ufile
spin = args.spin

#load umatrix
print "Loading Umatrix"
umatrix = rw.read_umatrix(filename,spin=spin)
umatrix = np.real(umatrix)
print "umatrix.shape", umatrix.shape

if (not spin):
	print "Temporary blowing up U matrix"
	newumatrix = np.zeros((umatrix.shape[0],) + (2,) + (umatrix.shape[1],) + (2,) + (umatrix.shape[2],) + (2,) + (umatrix.shape[3],) + (2,))
	newumatrix[:,0,:,0,:,0,:,0] = umatrix[:,:,:,:]
	newumatrix[:,1,:,0,:,1,:,0] = umatrix[:,:,:,:]
	newumatrix[:,0,:,1,:,0,:,1] = umatrix[:,:,:,:]
	newumatrix[:,1,:,1,:,1,:,1] = umatrix[:,:,:,:]
	umatrix = newumatrix
	print "umatrix.shape", umatrix.shape

bands = umatrix.shape[0]
spins = umatrix.shape[1] 

Udd = []
Jdd = []
Vdd = []
rest = []
U = np.zeros_like(umatrix)
J = np.zeros_like(umatrix)
V = np.zeros_like(umatrix)
spindict = {0:'u',1:'d'}

for spin1 in range(0,spins):
	for spin2 in range(0,spins):
		for spin3 in range(0,spins):
			for spin4 in range(0,spins):
				if spin1 == spin3 and spin2 == spin4:
					for band1 in range(0,bands): 
						for band2 in range(0,bands): 
							for band3 in range(0,bands): 
								for band4 in range(0,bands): 
									if band1 == band2 == band3 == band4 and spin1 != spin2:
										Udd.append(umatrix[band1,spin1,band2,spin2,band3,spin3,band4,spin4])
										U[band1,spin1,band2,spin2,band3,spin3,band4,spin4] = Udd[-1]
									elif band1 == band3 and band2 == band4 and band1 != band2:
										Vdd.append(umatrix[band1,spin1,band2,spin2,band3,spin3,band4,spin4])
										V[band1,spin1,band2,spin2,band3,spin3,band4,spin4] = Vdd[-1]
									elif band1 == band4 and band2 == band3 and band1 != band2 and spin1 == spin2: 
										Jdd.append(umatrix[band1,spin1,band2,spin2,band3,spin3,band4,spin4])
										J[band1,spin1,band2,spin2,band3,spin3,band4,spin4] = Jdd[-1]
									elif band1 == band2 and band3 == band4 and band1 != band3 and spin1 != spin2: 
										Jdd.append(umatrix[band1,spin1,band2,spin2,band3,spin3,band4,spin4])
										J[band1,spin1,band2,spin2,band3,spin3,band4,spin4] = Jdd[-1]
									elif band1 == band4 and band2 == band3 and band1 != band3 and spin1 != spin2: 
										Jdd.append(umatrix[band1,spin1,band2,spin2,band3,spin3,band4,spin4])
										J[band1,spin1,band2,spin2,band3,spin3,band4,spin4] = Jdd[-1]
									
										
print ""
print "Nonzero U entries:"
for index, value in np.ndenumerate(U):
	if value != 0:
		print str(index[0]+1)+spindict[index[1]], str(index[2]+1)+spindict[index[3]], str(index[4]+1)+spindict[index[5]], str(index[6]+1)+spindict[index[7]], value
print ""
print "Nonzero V entries:"
for index , value in np.ndenumerate(V):
	if value != 0:
		print str(index[0]+1)+spindict[index[1]], str(index[2]+1)+spindict[index[3]], str(index[4]+1)+spindict[index[5]], str(index[6]+1)+spindict[index[7]], value
print ""
print "Nonzero J entries:"
for index, value in np.ndenumerate(J):
	if value != 0:
		print str(index[0]+1)+spindict[index[1]], str(index[2]+1)+spindict[index[3]], str(index[4]+1)+spindict[index[5]], str(index[6]+1)+spindict[index[7]], value
print ""
print "-----RESULTS-----"
print ""
#print "Udd: ", Udd                                        
print "mean(Udd) = ", np.mean(Udd)
print "std(Udd) = ", np.std(Udd)
NUdd = len(np.nonzero(Udd)[0])
print "Number of entries contributing: " + str(NUdd)
print ""
#print "Vdd: ", Vdd
print "mean(Vdd) = ", np.mean(Vdd)
print "std(Vdd) = ", np.std(Vdd)
NVdd =  len(np.nonzero(Vdd)[0])
print "Number of entries contributing: " + str(NVdd)
print ""
#print "Jdd: ", Jdd
print "mean(Jdd) = ", np.mean(Jdd)
print "std(Jdd) = ", np.std(Jdd)
NJdd = len(np.nonzero(Jdd)[0])
print "Number of entries contributing: " + str(NJdd)
print ""
print "Total number of entries contributing: " + str(NUdd + NVdd + NJdd)
