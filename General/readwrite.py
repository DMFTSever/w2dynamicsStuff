#########################################################################################################################################################
#Script containing functions: 	1. readaline: reads next not comment line from file, neclecting everything after a comment letter. Appends \n at the end#
#                               2. checkline: checks whether a line contains a comment and returns everthing without the comment. Appends \n at the end #
#				3. read_hk_wannier: reads a Hamiltonian in wannier format								#
#				4. write_hk_wannier: writes a Hamiltonian to wannier format								#
#                               + all lot more....
#Severino Adler 24.09.2018                                                                                                                              #
#########################################################################################################################################################

import numpy as np
import sys
import custom_errors as err

def readaline(readfile, comments='#'):
    '''This is a function to reading the next line which is not a comment from a file ignoring everthing after the comment symbol'''
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

def checkline(line, comments='#'):
    '''This function checks whether there is a comment in line. It returns everthing up to the comment including a \n.'''
    commentpos=0
    for letter in line:
        if letter == comments:
            return line[0:commentpos] + '\n'
        commentpos += 1
    return line 

def read_hk_wannier(filename,spin=False):
	'''Function to read hk from a file with name 'filename'. spin dependency is False by default. If spin is set to True it assumes wannier convention (orbital index running fastest, spin index slowest). It returns the hk of form of an numpy array of form (number of kpoints, number of orbitals, number of spins, number of orbitals, number of spins), the kpoints in form of a numpy array of form (number of kpoints, x, y, z)\nreturn kpoints, hk'''
	f = open(filename, 'r')
	firstline = readaline(f)
	nk = int(firstline.split()[0]) #Number of k points of hk
	nd = int(firstline.split()[1]) #Number of orbitals (inkluding spin) of hk
	hk = np.zeros((nk,nd,nd),dtype=complex)
	kpoints = np.zeros((nk,3))
	for i in range(0,nk):
		splitline = readaline(f).split()
		kpoints[i,0] = float(splitline[0])
		kpoints[i,1] = float(splitline[1])
		kpoints[i,2] = float(splitline[2])
		for j in range(0,nd):
			splitline = readaline(f).split()
			for k in range(0,nd):
				hk[i,j,k] = complex(float(splitline[2*k]),float(splitline[2*k+1]))
	if spin:
		hk = hk.reshape(nk, 2, nd/2, 2, nd/2)
		hk = hk.transpose(0,2,1,4,3)
	else:
		hk = hk.reshape(nk, 1, nd, 1, nd)
		hk = hk.transpose(0,2,1,4,3)
	
	f.close()
	return hk, kpoints

def write_hk_wannier(filename,hk,kpoints):
	'''Function to write hk to a file named 'filename' in wannier format (spin indices running slower than orbital indices).\
       Assumes hk to be of the form (nk,nbands,spins,nbands,spins) and kpoint(nk,3). nd is the number of orbitals \
       (excluding spin), nk the number of kpoints.'''
	f = open(filename, 'w')
	nk = hk.shape[0]
	nbands = hk.shape[1]
	nspin = hk.shape[2]
	hk = hk.transpose(0,2,1,4,3)			
	hk = hk.reshape(nk,nspin*nbands,nspin*nbands)

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

def read_umatrix(filename, spin):
	'''Function to read a umatrix from a .dat file. 2 possible Formats
	   1.: 1 1 1 1 value (spin==False)
	   2.: 1u 1d 1u 1d value (spin==True)
        '''
	try:
		if type(filename) is not str: raise TypeError("1. argument (filename) must be of type str")
		if type(spin) is not bool: raise TypeError("2. argument (spin) must be of type bool:")
	except TypeError:
		raise
		sys.exit()

	f=open(filename, 'r')
	firstline = readaline(f).split()
	try:
		Nbands = int(firstline[0])
		if firstline[1].lower() != "bands": raise err.InputError("Expecting first non comment line to be of the form: # BANDS")
	except err.InputError:
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
		raise err.InputError("Specifiefied spin dependency does not match ufile spin dependency")
		sys.exit()	

	try:
		if spininfile!=spin: raise err.InputError("Specifiefied spin dependency does not match ufile spin dependency")
	except err.InputError:
		raise
		sys.exit()
		
	return umatrix	

def write_umatrix(filename, umatrix, spin):
	'''Funtion to write the umatrix either with spin dependency or without'''
	try:
		if type(filename) is not str: raise TypeError("1. argument (filenname) must be of type str.")
		if type(umatrix) is not np.ndarray: raise TypeError("2. argument (umatrix) must be of type numpy.ndarray.")
		if type(spin) is not bool: raise TypeError("3. argument (spin) must be of type bool.")
		if len(umatrix.shape)!=4 and len(umatrix.shape)!=8: raise err.ShapeError("Shape of umatrix is invalid")
		if len(umatrix.shape)==4 and spin: raise err.ShapeError("Spin independent umatrix detected although spin was set to True.")
		if len(umatrix.shape)==8 and not spin: raise err.ShapeError("Spin dependent umatrix detected although spin was set to False.")
		if spin and (umatrix.shape[1]!=2 or umatrix.shape[3]!=2 or umatrix.shape[5]!=2 or umatrix.shape[7]!=2):
			raise err.ShapeError("Unexpected shape of spin dependent umatrix. Expected shape: (Nbands,2,Nbands,2,Nbands,2,Nbands,2)")
	except TypeError:
		raise
		sys.exit()
	except err.ShapeError:
		raise
		sys.exit()
	
	f = open(filename, 'w')
	f.write("# non zero elements of interaction matrix U_ijkl\n")
	f.write("%i BANDS\n"%umatrix.shape[0])
	spindict={0:'u',1:'d'}
	if not spin:
		for index, value in np.ndenumerate(umatrix):
			if value != 0.:
				f.write("%i %i %i %i  %010.10f\n"%(index[0]+1,index[1]+1,index[2]+1,index[3]+1, value))
	else:
		for index, value in np.ndenumerate(umatrix):
			if value != 0:
				f.write("%s %s %s %s  %010.10f\n"%(str(index[0]+1)+spindict[index[1]],str(index[2]+1)+spindict[index[3]],\
							    str(index[4]+1)+spindict[index[5]],str(index[6]+1)+spindict[index[7]],value))
		

def write_function(filename, axis, axisname, function, functionname, header=""):
	"""Function to write a variable living on axis to a file. Including a header, which will be included with a preceding #"""
	try:
		if type(filename) is not str: raise TypeError("1. argument (filename) must be of type str")
		if type(axisname) is not str: raise TypeError("3. argument (axisname) must be of type str")
		if type(functionname) is not str: raise TypeError("5. argument (functionname) must be of type str")
		if type(header) is not str: raise TypeError("6. argument (header) must be of type str")
		if type(axis) is not np.ndarray: raise TypeError("2. argument (axis) must be of type numpy.array")
		if type(function) is not np.ndarray: raise TypeError("4. argument (function) must be of type numpy.array")
		f = open(str(filename),'w')
		if len(axis.shape) != 1 or len(function.shape) != 1 or axis.shape != function.shape: raise err.ShapeError("Shape missmatch of axis and/or function, both must have one dimension and same shape!")
	except err.ShapeError as error:
		raise
		sys.exit()
	except TypeError as error:
		raise
		sys.exit()
	except:
		raise
		print "Error: Cannot open file: " + str(filename)
		sys.exit()
	
	formats = [10,10]
	if len(str(axis[0])) > 10: formats[0] = len(str(axis[0]))	
	if len(str(max(np.real(function[0]),np.imag(function[0])))) > 10: formats[1] = len(str(max(np.real(function[0]),np.imag(function[0]))))	
	strf = ["%+"+str(formats[0]-1)+"s","%+"+str(formats[1])+"s"]
	ff = ["%+"+str(formats[0])+".8f","%+"+str(formats[1])+".8f"]

	if header != "":
		f.write("# " + header)	
		f.write("\n")
	f.write("#")
	f.write(strf[0] % axisname)
	f.write("  ")
	f.write(strf[1] %("Re("+functionname+")"))
	f.write("  ")
	f.write(strf[1] %("Im("+functionname+")"))
	f.write("\n")
	for i in range(0, axis.shape[0]):
		f.write(ff[0] % axis[i])
		f.write("  ")
		f.write(ff[1] % np.real(function[i]))
		f.write("  ")
		f.write(ff[1] % np.imag(function[i]))
		f.write("\n")
	f.close()
		 
def write_np_array(filename, arr, indexshift=0):
	try:
		if type(filename) is not str: raise TypeError("1. argument (filename) must of type str!")
		if type(arr) is not np.ndarray: raise TypeError("2. argument (array) must be of type numpy.ndarray!")
		if type(indexshift) is not str: raise TypeError("3. argument (indexshift) must be of type int.")
		if True in (np.array(arr.shape) > 100000): raise err.ShapeError("Shape of input array must NOT have entries > 100000!")
	except TypeError:
		raise
		sys.exit()

	writefile = open(filename, 'w')
	arrshape = arr.shape
	writefile.write(str(arrshape))
	writefile.write('\n')
	iterlist = []
	for element in arrshape:
		if len(iterlist) == 0:
			for i in range(0,element):
				iterlist.append((i,))
		else:	
			tempiterlist = iterlist
			iterlist = []
			for iterelement in tempiterlist:
				for i in range(0,element):	
					iterlist.append(iterelement + (i,))
	for indices in iterlist:
		for index in indices:
			writefile.write("%5i" %(index+indexshift))
			writefile.write(" ")
		writefile.write("  ")
		writefile.write(str(arr[indices]))
		writefile.write('\n')
	writefile.close()

	
if __name__ == "__main__":
	### here we need to add the path to the w2dynamics main directory
	auxdir="/home/lv70961/sever2/programs/w2dynamics/"
	sys.path.insert(0,auxdir)
	### in order to import the functions available in input.py and interatction.py and readwrite and custom errors
	import auxiliaries.input as io
	import dmft.interaction as ia
	
	hfile = file("Hk_LiOsO3_SOC.dat")
	hk, kpoints = io.read_hamiltonian(hfile, spin_orbit=True)
	print hk.shape
	hk1, kpoints1 = read_hk_wannier("Hk_LiOsO3_SOC.dat", spin=True)
	print hk1.shape
	print np.allclose(hk, hk1)
	print np.allclose(kpoints, kpoints1)

