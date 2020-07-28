from __future__ import print_function, division, absolute_import
import numpy as np
import matplotlib.pyplot as plt
import sys
import argparse
import glob
import os
import fnmatch
import h5py 
import matplotlib.cm as colormap
import matplotlib.colors as mcolors 
from collections import OrderedDict ### here we need to add the path to the personal myfunc directory (location of some read functions)
sys.path.insert(0,sys.argv[0].replace('/InputOutputManipulators/generate_maxent_plots.py','/General'))
### in order to import the functions available in input.py and interatction.py and readwrite and custom errors
import readwrite as rw 
import custom_errors as err 

def check_xrange(strrange):
	try:
		strrange = str(strrange)
		if strrange == ':':
			ranges = [0,-1]
		else:
			maxval = float(strrange)
			if maxval <= 0: raise err.InputError("xrange value must be > 0")
			ranges = [-maxval,maxval]
	except:
		raise
		sys.exit()
	else:
		return ranges	

def check_yrange(strrange):
	try:
		strrange = str(strrange)
		if not ':' in strrange: raise err.InputError("Type of yrange must be min:max with max > min or only :") 
		elif strrange == ':':
			ranges = [0,-1]
		else:
			ranges = strrange.split(':')
			ranges = [float(ranges[0]),float(ranges[1])]
			if ranges[0] >= ranges[1]: 
				raise err.InputError("Type of yrange must be min:max with max > min > 0 or only :")
	except:
		raise
		sys.exit()
	else:
		return ranges	

def check_hamiltonian(hamiltonian):
	possibles = ["Density", "Kanamori", "ReadUmatrix", "ReadNormalUmatrix"]
	hamiltonian = str(hamiltonian)
	if hamiltonian not in possibles: raise err.InputError("Found not allowed Hamiltonian")
	return hamiltonian

def check_channel(strspin):
	spindict = {'up':0, 'down':1}
	try:
		if strspin.lower() in spindict.keys(): return spindict[strspin.lower()]	
		else: raise err.InputError("Specified spin unknown. Must be either 'up' oder 'down'")
	except: 
		raise
		sys.exit()

parser = argparse.ArgumentParser(description="This script plots A(omega) over beta. It searches recusivly through all \
					      directories starting from startpath searching for a Maxent directory")
parser.add_argument('path', help='PATH from where to start recursivly searching for Maxent directories.', type=str)
parser.add_argument('basetitle', help='Basename for the title of the plots.', type=str)
parser.add_argument('--xrange', help="range of w values which will be plotted. Input style: max or :. max > 0.\
		    		    Actual range will be [-max:max]. : is interpreted as full range which is also default"\
		  , type=check_xrange, default=[0,-1])
parser.add_argument('--yrange', help="range of y-axis which will be plotted. Input style: min:max. max > min >= 0. Actuall range\
				       will be [min:max]", default=[0,-1], type=check_yrange)
parser.add_argument('-H','--Hamiltonians', help='Hamiltonian types to plot. Choose from list: Density, Kanamori, ReadUmatrix, \
						 ReadNormalUmatrix', type=check_hamiltonian, nargs='*')
parser.add_argument('-B','--Betas', help='Beta values to plot. If Lambdas and Betas are specified one of both must be only one \
                                          value.', type=int, nargs='*', default=None)
parser.add_argument('-L','--Lambdas', help='Lambda values to plot. If Lambdas are specified Betas have to be secified as well \
                                            and one of both must be only one value.', type=float, nargs='*', default=None)
parser.add_argument('-A', '--Atoms', help='Inequivalent atoms to plot. Each entry must be smaller then Nat of found hdf5 file.\
					   If the inequivalent atom does not exist in one or many hdf5 files it will be skipped\
					   for thoose directories. Default = 1', type=int, default=[1], nargs='*')
parser.add_argument('-S','--Spin', help="Use followed 'by' up or 'down' in order to plot only one spin channel (up,down). Per default\
					 both spin channels get ploted", type=check_channel, default = 2)
parser.add_argument('--ylabels', help='Label y-axis', action = 'store_true', default = False)
parser.add_argument('--svg', help='Save all files as svg. Default is pdf', action = 'store_true', default = False)
parser.add_argument('--maxentdir', help='Set directory name for maxent calculation. Default is Maxent', type=str, default='Maxent')
parser.add_argument('--xi', help='use xi instead of lambda', action = 'store_true', default = False)
args = parser.parse_args()

startpath = args.path
basetitle = args.basetitle
xranges = args.xrange
yranges = args.yrange
basetitle = args.basetitle
atoms = set(args.Atoms)
maxatom = max(atoms)
searchpatterns = []
spin = args.Spin
svg = args.svg
maxentdir = args.maxentdir
lambdas = args.Lambdas
ylabels = args.ylabels
xi = args.xi

plotbetas=True
plotlambdas=False

if lambdas != None  and args.Betas == None:
	raise err.InputError('Please specify --Betas as well')
	sys.exit()

if lambdas != None and args.Betas != None:
        if len(lambdas)>1 and len(args.Betas)>1:
		raise err.InputError('If both Betas and Lambdas are specified one of both must be only one value')
		sys.exit()
	elif len(args.Betas)==1:
		plotlambdas=True
		plotbetas=False

if args.Betas is not None:
	searchpatterns += [['general.beta',args.Betas]]
for atom in atoms:
	if args.Hamiltonians is not None:
		searchpatterns += [['atoms.' + str(atom) +'.hamiltonian',args.Hamiltonians]]
hamiltonians = args.Hamiltonians

#Walking through the directories to obtain the newest hdf5 file
temppaths = []
filedirpaths = []
for dirpath, dirnames, filenames in os.walk(startpath):
	for filename in fnmatch.filter(filenames, '*.hdf5'):
		temppaths.append(os.path.join(dirpath, filename))
	else:
		if len(temppaths) > 0:
			filedirpaths.append([max(temppaths, key=os.path.getctime), dirpath])
		temppaths = []

usepaths = []

if len(searchpatterns) != 0:
	for filepath, dirpath in filedirpaths:
		tempfile = h5py.File(filepath)
		usefile = True
		for element in searchpatterns:
			atr = tempfile["/.config"].attrs.get(element[0])
			if atr not in element[1]: usefile = False
		try:
			if tempfile["/.config"].attrs.get('general.nat') < maxatom: 
				raise err.InputError("Atom Index list out of range for " + filepath)
		except err.InputError:
			raise
			sys.exit()
		if usefile:
			usepaths.append([filepath, dirpath])
		tempfile.close()
else:
	usepaths = filedirpaths

if lambdas is not None:
	temppaths = usepaths[:]
	usepaths = []
	for path, dirpath in temppaths:
		lambdaa = 0
		parts=path.split('/')[-1].split('_')
		for element in parts:
			if 'lambda' in element:
				lambdaa = float(element[6:])
		if lambdaa in lambdas:
			usepaths.append([path,dirpath])

if len(usepaths) == 0: raise RuntimeError("No matching files found")

#extracing number of bands from first file and first atom
tempfile = h5py.File(usepaths[0][0])	
Nbands = int(tempfile["/.config"].attrs.get("atoms.1.nd"))

plotbool = False
tempusepaths = []
for i in range(0,len(usepaths)):
	if maxentdir in os.listdir(usepaths[i][1]):
		tempusepaths.append([usepaths[i][0], os.path.join(usepaths[i][1], maxentdir)])
		plotbool = True
usepaths = tempusepaths
if plotbool == False: raise RuntimeError("No Maxent directory found")
	
for pair in usepaths:
	print(pair)

#PLOTTING
plt.rcParams["text.usetex"] = True
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.sans-serif'] = 'Computer Modern Sans Serif'

if plotbetas==True:
	cms = [colormap.jet, colormap.gnuplot, colormap.gnuplot2, colormap.brg, colormap.rainbow, colormap.ocean, colormap.brg]
	if args.Betas==None:
		norm = mcolors.Normalize(1./52.,1./4.6) #hardcoded beta values between 5 and 50 here!
	else:
		temp = np.sort(args.Betas)
		norm = mcolors.Normalize(1./(temp[-1]+0.4),1./(temp[0]-0.4))
elif plotlambdas == True:
	cms = [colormap.jet, colormap.gnuplot, colormap.gnuplot2, colormap.brg, colormap.rainbow, colormap.ocean, colormap.brg]
	temp = np.sort(lambdas)
	norm = mcolors.Normalize(1./(temp[-1]+1),1./(temp[0]+1))
spindict = {0:'up',1:'down'}
hamiltoniandict = {'Density':'Density', 'Kanamori':'Kanamori', 'ReadUmatrix':'full Coulomb', 'ReadNormalUmatrix':'full Coulomb'}

if hamiltonians == None: hamiltonians = ["Density", "Kanamori", "ReadUmatrix", "ReadNormalUmatrix"]
for hamiltonian in hamiltonians:
	
	print("")
	print("Plotting for Hamiltonian: ", hamiltonian)
	#Generating figures

	if spin == 2 and ('12' in os.listdir(usepaths[0][1] + "/atom1/")) :
		#only diagonal entries
		figaw, mainaw = plt.subplots(Nbands+1,2,figsize=(15,20))
		figaw.delaxes(mainaw[Nbands,1])
		box = mainaw[Nbands,0].get_position()
                box.x0 = box.x0 + 0.21
                box.x1 = box.x1 + 0.21
                mainaw[Nbands,0].set_position(box)
		print(type(mainaw))
		figaw.suptitle(r"Diagonal entries of $A(\omega)$" +'\nfor ' + hamiltoniandict[hamiltonian] + ' Hamiltonian', fontsize=20)
		for i in range(0,Nbands+1):
			if xranges[1] != -1:
				mainaw[i,0].set_xlim([xranges[0],xranges[1]])
				mainaw[i,1].set_xlim([xranges[0],xranges[1]])
			if yranges[1] != -1:
				mainaw[i,0].set_ylim([yranges[0],yranges[1]])
				mainaw[i,1].set_ylim([yranges[0],yranges[1]])
			if i != Nbands:
				mainaw[i,0].set_title(r"$A(\omega)$ of orbital " + str(i+1) + " spin up")
				mainaw[i,1].set_title(r"$A(\omega)$ of orbital " + str(i+1) + " spin down")
				mainaw[i,0].set_xlabel(r"$\omega$")
				mainaw[i,1].set_xlabel(r"$\omega$")
				if ylabels==True:
					mainaw[i,0].set_ylabel(r"$A(\omega)$")
					mainaw[i,1].set_ylabel(r"$A(\omega)$")
			else:
				mainaw[i,0].set_title(r"Average $A(\omega)$ of all orbitals and spins")
				mainaw[i,0].set_xlabel(r"$\omega$")
				if ylabels==True:
					mainaw[i,0].set_ylabel(r"$A(\omega)$")

	else:
		if spin==2: spin = 0
		#only diagonal entries
		figaw, mainaw = plt.subplots(Nbands+1,1,figsize=(7.5,20))
		figaw.suptitle(r"Diagonal entries of $A(\omega)$" +'\nfor ' + hamiltoniandict[hamiltonian] + ' Hamiltonian', fontsize=20)
		for i in range(0,Nbands):
			if xranges[1] != -1:
				mainaw[i].set_xlim([xranges[0],xranges[1]])
			if yranges[1] != -1:
				mainaw[i].set_ylim([yranges[0],yranges[1]])
			mainaw[i].set_title(r"$A(\omega)$ of orbital " + str(i+1) + " spin up")
			mainaw[i].set_xlabel(r"$\omega$")
		if xranges[1] != -1:
			mainaw[Nbands].set_xlim([xranges[0],xranges[1]])
		if yranges[1] != -1:
			mainaw[Nbands].set_ylim([yranges[0],3*yranges[1]])
		mainaw[Nbands].set_title(r"Average $A(\omega)$ of all orbitals spin up")
		mainaw[Nbands].set_xlabel(r"$\omega$")

	print("Entering for loop to extraxt and plot data")
	foundfile = False
	for filepath, dirpath in usepaths:
		lambdaa = 0
		parts=filepath.split('/')[-1].split('_')
		for element in parts:
			if 'lambda' in element:
				lambdaa = float(element[6:])
		print('')
		print(dirpath)

		#extracting A(omega) of atoms to plot
		print("Extracting metadata of directory ", dirpath)
		f = h5py.File(filepath)
		hamilton = str(f["/.config"].attrs.get("atoms.1.hamiltonian")).capitalize()
		if hamilton.lower() == hamiltonian.lower():
			foundfile = True
			beta = int(f["/.config"].attrs.get("general.beta"))
			ineqatoms = [x for x in range(1,len(os.listdir(dirpath))+1)]
			plotatoms = [x for x in atoms if x in ineqatoms]	
			f = open(dirpath + "/atom1/11/maxent.in") #using first atom first band here
			for line in f:
				line = rw.checkline(line)
				split=line.split("=")
				if split[0]=='NFREQ':
					nfreq = int(split[1])
			f.close()
			awshape = (4,2,nfreq)
			aw = np.zeros((len(plotatoms),)+tuple(awshape))
			waxis = np.zeros((nfreq,))

			f = open(dirpath + "/atom1/11/spec_111.dat")
			i=0
			for line in f:
				line = rw.checkline(line)
				if line == '\n': continue
				waxis[i] = float(line.split()[0])
				i+=1
			f.close()
			
			for i in range(1,len(plotatoms)+1):
				for j in range(1,Nbands+1):
					if spin == 2 and ('12' in os.listdir(dirpath + "/atom1/")) :
						for s in range(1,3):	
							f = open(dirpath + "/atom" + str(i) + "/" + str(j) + str(s) + "/spec_" + str(i) + str(j) + str(s) + ".dat")
							k=0
							for line in f:
								line = rw.checkline(line)
								if line == '\n': continue
								aw[i-1,j-1,s-1,k] = line.split()[1]
								k+=1
							f.close()
								
					else:
						if spin == 2: spin = 0
						f = open(dirpath + "/atom" + str(i) + "/" + str(j) + str(spin+1) + "/spec_" + str(i) + str(j) + str(spin+1) + ".dat")
						k=0
						for line in f:
							line = rw.checkline(line)
							if line == '\n': continue
							aw[i-1,j-1,spin,k] = line.split()[1]
							k+=1
						f.close()
				if spin == 2 and ('12' in os.listdir(dirpath + "/atom1/")) :
					for s in range(1,3):	
						f = open(dirpath + "/atom" + str(i) + "/sum/spec.dat")
						k=0
						for line in f:
							line = rw.checkline(line)
							if line == '\n': continue
							aw[i-1,3,s-1,k] = line.split()[1]
							k+=1
						f.close()
							
				else:
					if spin == 2: spin = 0
					f = open(dirpath + "/atom" + str(i) + "/sum/spec.dat")
					k=0
					for line in f:
						line = rw.checkline(line)
						if line == '\n': continue
						aw[i-1,3,spin,k] = line.split()[1]
						k+=1
					f.close()
		
			#add to plots of figures
			print("Adding data to plots")
			if plotbetas== True:
				label= r"$\beta$ " + str(beta)
				normval = 1./beta
			elif plotlambdas == True:
				label = r" $\xi$ " + str(lambdaa*2)
				normval = 1./(lambdaa+1)
			if spin == 2:
				for i in range(0,len(plotatoms)):
					cm = cms[i%5]
					lw = 2
					for j in range(0,Nbands):
						mainaw[j,0].plot(waxis[:], np.real(aw[i,j,0,:]), label=label + " atom " + str(plotatoms[i]), lw=lw, c=cm(norm(normval)))
						mainaw[j,1].plot(waxis[:], np.real(aw[i,j,1,:]), label=label + " atom " + str(plotatoms[i]), lw=lw, c=cm(norm(normval)))
					mainaw[Nbands,0].plot(waxis[:], np.real(aw[i,Nbands,0,:])/(Nbands*2), label=label + " atom " + str(plotatoms[i]), lw=lw, c=cm(norm(normval)))
					
			else:
				for i in range(0,len(plotatoms)):
					cm = cms[i%5]
					lw = 2
					for j in range(0,Nbands):
						mainaw[j].plot(waxis[:], np.real(aw[i,j,spin,:]), label=label + " atom " + str(plotatoms[i]), lw=lw, c=cm(norm(normval)))
					mainaw[Nbands].plot(waxis[:], np.real(aw[i,Nbands,spin,:])/Nbands, label=label + " atom " + str(plotatoms[i]), lw=lw, c=cm(norm(normval)))
	if foundfile:
		print("Found at least one file for this Hamiltonian. Generating output files")

		if spin==2:
			for i in range(0,Nbands):
				mainaw[i,0].legend(bbox_to_anchor=(0,1), loc='upper left')
				mainaw[i,1].legend(bbox_to_anchor=(0,1), loc='upper left')
			mainaw[Nbands,0].legend(bbox_to_anchor=(0,1), loc='upper left')
		else:
			for i in range(0,Nbands+1):
				mainaw[i].legend(bbox_to_anchor=(1,1), loc='upper right')

		if svg:
			filenameaw = basetitle + "_" + hamiltonian + "_spec.svg"
		else:
			filenameaw = basetitle + "_" + hamiltonian + "_spec.pdf"
		figaw.savefig(filenameaw)
	else:
		print("Found no file for this Hamiltonian. No output generated.")
				
