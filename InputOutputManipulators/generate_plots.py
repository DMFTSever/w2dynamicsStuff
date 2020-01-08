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
auxdir="/dss/dsshome1/0D/di76rir/lib/MyPythonScripts/W2dynamics/General"
sys.path.insert(0,auxdir)
### in order to import the functions available in input.py and interatction.py and readwrite and custom errors
import readwrite as rw 
import custom_errors as err 

def check_xrange(strrange):
	try:
		strrange = str(strrange)
		if not ':' in strrange: raise err.InputError("Type of xrange must be min:max with max > min >= 0 or only :") 
		elif strrange == ':':
			ranges = [0,-1]
		else:
			ranges = strrange.split(':')
			ranges = [float(ranges[0]),float(ranges[1])]
			if ranges[0] >= ranges[1] or ranges[0] < 0: 
				raise err.InputError("Type of xrange must be min:max with max > min >= 0 or only :")
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
			ranges = [0,1]
		else:
			ranges = strrange.split(':')
			ranges = [-float(ranges[1]),-float(ranges[0])]
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

parser = argparse.ArgumentParser(description="This script plots Tr(G(iw)), the diagonal channels of G(iw), S(iw) and some \
					      occupancies over beta. It searches recusivly through all directories starting\
					      from startpath and uses the newest file per directory. With the --allfiles\
					      option it uses all files per directory")
parser.add_argument('path', help='PATH from where to start recursivly searching for hdf5 files.', type=str)
parser.add_argument('xrange', help="range of iw values which will be plotted. Input style: min:max. max > min >= 0",\
		    type=check_xrange)
parser.add_argument('basetitle', help='Basename for the title of the plots.', type=str)
parser.add_argument('--syrange', help="range of y-axis for self energy values which will be plotted. Input style: min:max. max > min. Actuall range\
				       will be -max:-min", default=[0,1], type=check_yrange)
parser.add_argument('--gyrange', help="range of y-axis green function values which will be plotted. Input style: min:max. max > min. Actuall range\
				       will be -max:-min", default=[0,1], type=check_yrange)
parser.add_argument('-H','--Hamiltonians', help='Hamiltonian types to plot. Choose from list: Density, Kanamori, ReadUmatrix, \
						 ReadNormalUmatrix', type=check_hamiltonian, nargs='*')
parser.add_argument('-B','--Betas', help='Beta values to plot. If Lambdas and Betas are specified one of both must be only one \
                                          value.', type=int, nargs='*', default=None)
parser.add_argument('-L','--Lambdas', help='Lambda values to plot. If Lambdas are specified Betas have to be secified as well \
                                            and one of both must be only one value.', type=float, nargs='*', default=None)
parser.add_argument('-A', '--Atoms', help='Inequivalent atoms to plot. Each entry must be smaller then Nat of given hdf5 file.\
					   If the inequivalent atom does not exist in one or many hdf5 files it will be skipped\
					   for thoose files. Default = 1', type=int, default=[1], nargs='*')
parser.add_argument('-S','--Spin', help="Use followed 'by' up or 'down' in order to plot only one spin channel (up,down). Per default\
					 both spin channels get ploted", type=check_channel, default = 2)
parser.add_argument('--allfiles', help='Plot all files. Not only newest per directory', action = 'store_true', default = False)
parser.add_argument('--ylabels', help='Label y-axis', action = 'store_true', default = False)
parser.add_argument('--svg', help='Save all files as svg. Default is pdf', action = 'store_true', default = False)
parser.add_argument('--xi', help='use xi instead of lambda', action = 'store_true', default = False)
parser.add_argument('--useiter', help='Use this iteration if hdf5 file does not contain last iteration', type=int, default = 0)
args = parser.parse_args()

startpath = args.path
allfiles =args.allfiles
basetitle = args.basetitle
xranges = args.xrange
syranges = args.syrange
gyranges = args.gyrange
basetitle = args.basetitle
atoms = set(args.Atoms)
maxatom = max(atoms)
searchpatterns = []
spin = args.Spin
svg = args.svg
lambdas = args.Lambdas
ylabels = args.ylabels
xi = args.xi
useiter = args.useiter

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
filepaths = []
for dirpath, dirnames, filenames in os.walk(startpath):
	for filename in fnmatch.filter(filenames, '*.hdf5'):
		temppaths.append(os.path.join(dirpath, filename))
	if allfiles:
		for i in range(0,len(temppaths)):
			filepaths.append(temppaths[i])
	else:
		if len(temppaths) > 0:
			filepaths.append(max(temppaths, key=os.path.getctime))
		temppaths = []

usepaths = []
if len(searchpatterns) != 0:
	for filepath in filepaths:
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
			usepaths.append(filepath)
		tempfile.close()
else:
	usepaths = filepaths

if lambdas is not None:
	temppaths = usepaths[:]
	usepaths = []
	for path in temppaths:
		lambdaa = 0
		parts=path.split('/')[-1].split('_')
		for element in parts:
			if 'lambda' in element:
				lambdaa = float(element[6:])
		if lambdaa in lambdas:
			usepaths.append(path)

if len(usepaths) == 0: raise RuntimeError("No matching files found")

#extracing number of bands from first file and first atom
tempfile = h5py.File(usepaths[0])	
Nbands = int(tempfile["/.config"].attrs.get("atoms.1.nd"))

#PLOTTING
plt.rcParams["text.usetex"] = True
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.sans-serif'] = 'Computer Modern Sans Serif'

if plotbetas==True:
	if args.Betas==None:
		norm = mcolors.Normalize(1./52.,1./4.6) #hardcoded beta values between 5 and 50 here!
	else:
		temp = np.sort(args.Betas)
		norm = mcolors.Normalize(1./(temp[-1]+0.4),1./(temp[0]-0.4))
elif plotlambdas == True:
	temp = np.sort(lambdas)
	norm = mcolors.Normalize(1./(temp[-1]+1),1./(temp[0]+1))
	
if plotbetas == True:
	cms = [colormap.jet, colormap.gnuplot, colormap.gnuplot2, colormap.brg, colormap.rainbow, colormap.ocean, colormap.brg]
elif plotlambdas == True:
	cms = [colormap.jet, colormap.gnuplot, colormap.gnuplot2, colormap.brg, colormap.rainbow, colormap.ocean, colormap.brg]
	#cms = [colormap.gist_ncar, colormap.nipy_spectral, colormap.gist_ncar, colormap.nipy_spectral, colormap.gist_ncar, colormap.nipy_spectral, colormap.brg]
spindict = {0:'up',1:'down'}
hamiltoniandict = {'Density':'Density', 'Kanamori':'Kanamori', 'ReadUmatrix':'full Coulomb', 'ReadNormalUmatrix':'full Coulomb'}

if hamiltonians == None: hamiltonians = ["Density", "Kanamori", "ReadUmatrix", "ReadNormalUmatrix"]
for hamiltonian in hamiltonians:
	
	print ""
	print "Plotting for Hamiltonian: ", hamiltonian
	#Generating figures
	figsign, mainsign = plt.subplots(1,1,figsize=(6.5,5))
	figsign.suptitle("Average sign of " + hamiltoniandict[hamiltonian] + " Hamiltonian", fontsize=15)
	mainsign.set_ylabel("Average sign")
	signlabel= r"$\beta$"
	if plotlambdas == True:
		if xi == True:
			signlabel = r"$\xi$"
		else:
			signlabel = r"$\lambda$"
	mainsign.set_xlabel(signlabel)
	if spin==2:
		#only diagonal entries
		figsiw, mainsiw = plt.subplots(Nbands,2,figsize=(15,15))
		figsiw.suptitle(r"Diagonal entries of Im($\Sigma(\omega_n)$)" +'\nfor ' + hamiltoniandict[hamiltonian] + ' Hamiltonian', fontsize=20)
		for i in range(0,Nbands):
			if xranges[1] != -1:
				mainsiw[i,0].set_xlim([xranges[0],xranges[1]])
				mainsiw[i,1].set_xlim([xranges[0],xranges[1]])
			if syranges[1] != 1:
				mainsiw[i,0].set_ylim([syranges[0],syranges[1]])
				mainsiw[i,1].set_ylim([syranges[0],syranges[1]])
			mainsiw[i,0].set_title(r"Im($\Sigma(\omega_n)$) of orbital " + str(i+1) + " spin up")
			mainsiw[i,1].set_title(r"Im($\Sigma(\omega_n)$) of orbital " + str(i+1) + " spin down")
			mainsiw[i,0].set_xlabel(r"$\omega_n$")
			mainsiw[i,1].set_xlabel(r"$\omega_n$")
			if ylabels:
				mainsiw[i,0].set_ylabel(r"Im($\Sigma(\omega_n)$)")
				mainsiw[i,1].set_ylabel(r"Im($\Sigma(\omega_n)$)")
		figgiw, maingiw = plt.subplots(Nbands,2,figsize=(15,15))
		figgiw.suptitle(r'Diagonal entries of Im(G($\omega_n$)) ' + '\nfor ' + hamiltoniandict[hamiltonian] + ' Hamiltonian', fontsize=20)
		for i in range(0,Nbands):
			if xranges[1] != -1:
				maingiw[i,0].set_xlim([xranges[0],xranges[1]])
				maingiw[i,1].set_xlim([xranges[0],xranges[1]])
			if gyranges[1] != 1:
				maingiw[i,0].set_ylim([gyranges[0],gyranges[1]])
				maingiw[i,1].set_ylim([gyranges[0],gyranges[1]])
			maingiw[i,0].set_title(r"Im($G(\omega_n)$) of orbital " + str(i+1) + " spin up")
			maingiw[i,1].set_title(r"Im($G(\omega_n)$) of orbital " + str(i+1) + " spin down")
			maingiw[i,0].set_xlabel(r"$\omega_n$")
			maingiw[i,1].set_xlabel(r"$\omega_n$")
			if ylabels:
				maingiw[i,0].set_ylabel(r"Im(G($\omega_n$))")
				maingiw[i,1].set_ylabel(r"Im(G($\omega_n$))")

		figtrgiw, maintrgiw = plt.subplots()
		figtrgiw.suptitle(r'Im(Tr(G($\omega_n$))) ' + '\nof ' + hamiltoniandict[hamiltonian] + ' Hamiltonian', fontsize=15)
		if xranges[1] != -1:
			maintrgiw.set_xlim([xranges[0],xranges[1]])
		if gyranges[1] != 1:
			maintrgiw.set_ylim([6*gyranges[0],gyranges[1]])
		maintrgiw.set_xlabel(r"$\omega_n$")
		if ylabels:
		        maintrgiw.set_ylabel(r"Im(Tr[G($\omega_n$)])")

		occlabel= r"$\beta$"
		if plotlambdas == True:
			if xi == True:
				occlabel = r"$\xi$"
			else:
				occlabel = r"$\lambda$"
		figoccn, mainocc = plt.subplots(4,2,figsize=(15,20))
		figoccn.suptitle('Occupation numbers \nof ' + hamiltoniandict[hamiltonian] + ' Hamiltonian', fontsize=20)
		mainocc[0,0].set_title(r"Average occupation per orbital for spin up")
		mainocc[0,1].set_title(r"Average occupation per orbital for spin down")
		mainocc[1,0].set_title(r"Average double occupation per orbital")
		mainocc[1,1].set_title(r"Relative average double occupation per orbital")
		mainocc[2,0].set_title(r"Average Hund's like occupation for spin up (i,$\uparrow$,j,$\uparrow$) $i \neq j$")
		mainocc[2,1].set_title(r"Average Hund's like occupation for spin down (i,$\downarrow$,j,$\downarrow$) $i \neq j$")
		mainocc[3,0].set_title(r"Relative average Hund's like occupation for spin up (i,$\uparrow$,j,$\uparrow$) $i \neq j$")
		mainocc[3,1].set_title(r"Relative average Hund's like occupation for spin down (i,$\downarrow$,j,$\downarrow$) $i \neq j$")
		mainocc[0,0].set_xlabel(occlabel)
		mainocc[0,1].set_xlabel(occlabel)
		mainocc[1,0].set_xlabel(occlabel)
		mainocc[1,1].set_xlabel(occlabel)
		mainocc[2,0].set_xlabel(occlabel)
		mainocc[2,1].set_xlabel(occlabel)
		mainocc[3,0].set_xlabel(occlabel)
		mainocc[3,1].set_xlabel(occlabel)

	else:
		#only diagonal entries
		figsiw, mainsiw = plt.subplots(Nbands,1,figsize=(7.5,15))
		figsiw.suptitle(r"Diagonal entries of Im($\Sigma(iw)$)" + '\nfor ' + hamiltoniandict[hamiltonian] + ' Hamiltonian', fontsize=20)
		for i in range(0,Nbands):
			if xranges[1] != -1:
				mainsiw[i].set_xlim([xranges[0],xranges[1]])
			if syranges[1] != 1:
				mainsiw[i].set_ylim([syranges[0],syranges[1]])
			mainsiw[i].set_title(r"Im($\Sigma(\omega_n)$) of orbital " + str(i+1) + " spin "+spindict[spin])
			mainsiw[i].set_xlabel(r"$\omega_n$")
			if ylabels:
				mainsiw[i].set_ylabel(r"Im($\Sigma(\omega_n)$)")
		figgiw, maingiw = plt.subplots(Nbands,1,figsize=(7.5,15))
		figgiw.suptitle(r'Diagonal entries of Im(G($\omega_n$)) ' + '\nfor ' + hamiltoniandict[hamiltonian] + ' Hamiltonian', fontsize=20)
		for i in range(0,Nbands):
			if xranges[1] != -1:
				maingiw[i].set_xlim([xranges[0],xranges[1]])
			if gyranges[1] != 1:
				maingiw[i].set_ylim([gyranges[0],gyranges[1]])
			maingiw[i].set_title(r"Im($G(\omega_n)$) of orbital " + str(i+1) + " spin " + spindict[spin])
			maingiw[i].set_xlabel(r"$\omega_n$")
			if ylabels:
				maingiw[i].set_ylabel(r"Im(G($\omega_n$))")

		figtrgiw, maintrgiw = plt.subplots()
		figtrgiw.suptitle(r'Im(Tr(G($\omega_n$))) ' + '\nof ' + hamiltoniandict[hamiltonian] + ' Hamiltonian', fontsize=15)
		if xranges[1] != -1:
			maintrgiw.set_xlim([xranges[0],xranges[1]])
		if gyranges[1] != 1:
			maintrgiw.set_ylim([4*gyranges[0],gyranges[1]])
		maintrgiw.set_xlabel(r"$\omega_n$")
		if ylabels:
		        maintrgiw.set_ylabel(r"Im(Tr[G($\omega_n$)])")

		occlabel= r"$\beta$"
		if plotlambdas == True:
			if xi == True:
				occlabel = r"$\xi$"
			else:
				occlabel = r"$\lambda$"
		figoccn, mainocc = plt.subplots(4,1,figsize=(7.5,20))
		figoccn.suptitle('Occupation numbers \nof ' + hamiltoniandict[hamiltonian] + ' Hamiltonian', fontsize=20)
		mainocc[0].set_title(r"Average occupation per orbital for spin " +spindict[spin])
		mainocc[1].set_title(r"Average double occupation per orbital")
		mainocc[2].set_title(r"Average Hund's like occupation for spin " + spindict[spin] + r" (i,$\uparrow$,j,$\uparrow$) $i \neq j$")
		mainocc[3].set_title(r"Relative average Hund's like occupation for spin " + spindict[spin] + r" (i,$\uparrow$,j,$\uparrow$) $i \neq j$")
		mainocc[0].set_xlabel(occlabel)
		mainocc[1].set_xlabel(occlabel)
		mainocc[2].set_xlabel(occlabel)
		mainocc[3].set_xlabel(occlabel)

	print "Entering for loop to extraxt and plot data"
	foundfile = False
	for filepath in usepaths:

		#extracting siw, giw, occ of atoms to plot
		print "Extracting metadata of file ", filepath
		lambdaa = 0
		parts=filepath.split('/')[-1].split('_')
		for element in parts:
			if 'lambda' in element:
				lambdaa = float(element[6:])
		f = h5py.File(filepath)
		hamilton = str(f["/.config"].attrs.get("atoms.1.hamiltonian")).capitalize()
		if hamilton.lower() == hamiltonian.lower():
			foundfile = True
			beta = int(f["/.config"].attrs.get("general.beta"))
			intiteration = int(f["/.config"].attrs.get("general.dmftsteps")) #using the last Iteration
			ineqatoms_str= [str(x) for x in fnmatch.filter(f["/dmft-001"].keys(),"ineq-*")] #using 1. Iteration
			ineqatoms = [int(x[-3:]) for x in ineqatoms_str]
			plotatoms = [x for x in atoms if x in ineqatoms]	
			for i in range(0,100):
				iteration = "%03i" %(intiteration-i)
				print "Trying to plot iteration Nr. %s" %iteration
				try:
					giwshape = f["/dmft-" + iteration + "/ineq-001/giw-full/value"][:].shape
					siwshape = f["/dmft-" + iteration + "/ineq-001/siw-full/value"][:].shape
					occshape = f["/dmft-" + iteration + "/ineq-001/occ/value"][:].shape
				except:
					if useiter != 0:
						print "Did not find iteration Nr. %s. Using iteration Nr. "%iteration + str(useiter)
						iteration = "%03i" %useiter
						giwshape = f["/dmft-" + iteration + "/ineq-001/giw-full/value"][:].shape
						siwshape = f["/dmft-" + iteration + "/ineq-001/siw-full/value"][:].shape
						occshape = f["/dmft-" + iteration + "/ineq-001/occ/value"][:].shape
						break
					continue
				break
			giw = np.zeros((len(plotatoms),)+tuple(giwshape), dtype=complex)
			trgiw = np.zeros((len(plotatoms),giwshape[-1]), dtype=complex)
			siw = np.zeros((len(plotatoms),)+tuple(siwshape), dtype=complex)
			try: occ
			except: at1fullocc = np.zeros((1,) + tuple(occshape))
			try: fullsigns
			except: fullsigns = np.zeros((1,len(plotatoms)), dtype=complex)
			occ = np.zeros((len(plotatoms),)+tuple(occshape))
			iwaxis = f["/.axes/iw"][:]
			lb = len(iwaxis)/2
			tempsigns = []
			for i in range(0,len(plotatoms)):
				atom = "%03i"%int(plotatoms[i])
				giw[i] = f["/dmft-" + iteration + "/ineq-" + atom + "/giw-full/value"][:]
				trgiw[i] = np.trace(giw[i].reshape(giwshape[0]*giwshape[1],giwshape[2]*giwshape[3],giwshape[4]),axis1=0,axis2=1)
				siw[i] = f["/dmft-" + iteration + "/ineq-" + atom + "/siw-full/value"][:]
				occ[i] = np.real(f["/dmft-" + iteration + "/ineq-" + atom + "/occ/value"][:])
				tempsigns.append(f["/dmft-" + iteration + "/ineq-" + atom + "/sign/value"].value)
			if not np.all(at1fullocc):
				at1fullocc[0] = occ[0]
				fullsigns[0] = tempsigns[:] 
				betalist = [beta]
				lambdalist = [lambdaa]
			else:
				temp = np.zeros((1,) + occshape)
				temp[0] = occ[0]
				at1fullocc = np.concatenate((at1fullocc, temp), axis=0)
				temp = np.zeros((1,len(plotatoms)), dtype=complex)
				temp[0] = tempsigns[:]
				fullsigns = np.concatenate((fullsigns, temp), axis=0)
				betalist.append(beta)
				lambdalist.append(lambdaa)

			#add to plots of figures
			print "Adding data to plots"

			if plotbetas== True:
				label= r"$\beta$ " + str(beta)
				normval = 1./beta
			elif plotlambdas == True:
				if xi == True:
					label = r" $\xi$ " + str(lambdaa*2)
				else:
					label = r" $\lambda$ " + str(lambdaa)
				normval = 1./(lambdaa+1)
			if spin == 2:
				for i in range(0,len(plotatoms)):
					cm = cms[i%5]
					lw = 3 - i*0.5
					mainsiw[0,0].plot(iwaxis[lb:], np.imag(siw[i,0,0,0,0,lb:]), label=label  + " atom " + str(plotatoms[i]), lw=lw, marker='o', c=cm(norm(normval)))
					mainsiw[0,1].plot(iwaxis[lb:], np.imag(siw[i,0,1,0,1,lb:]), label=label  + " atom " + str(plotatoms[i]), lw=lw, marker='o', c=cm(norm(normval)))
					mainsiw[1,0].plot(iwaxis[lb:], np.imag(siw[i,1,0,1,0,lb:]), label=label  + " atom " + str(plotatoms[i]), lw=lw, marker='o', c=cm(norm(normval)))
					mainsiw[1,1].plot(iwaxis[lb:], np.imag(siw[i,1,1,1,1,lb:]), label=label  + " atom " + str(plotatoms[i]), lw=lw, marker='o', c=cm(norm(normval)))
					mainsiw[2,0].plot(iwaxis[lb:], np.imag(siw[i,2,0,2,0,lb:]), label=label  + " atom " + str(plotatoms[i]), lw=lw, marker='o', c=cm(norm(normval)))
					mainsiw[2,1].plot(iwaxis[lb:], np.imag(siw[i,2,1,2,1,lb:]), label=label  + " atom " + str(plotatoms[i]), lw=lw, marker='o', c=cm(norm(normval)))
				for i in range(0,len(plotatoms)):
					cm = cms[i%5]
					lw = 3 - i*0.5
					maingiw[0,0].plot(iwaxis[lb:], np.imag(giw[i,0,0,0,0,lb:]), label=label  + " atom " + str(plotatoms[i]), lw=lw, marker='o', c=cm(norm(normval)))
					maingiw[0,1].plot(iwaxis[lb:], np.imag(giw[i,0,1,0,1,lb:]), label=label  + " atom " + str(plotatoms[i]), lw=lw, marker='o', c=cm(norm(normval)))
					maingiw[1,0].plot(iwaxis[lb:], np.imag(giw[i,1,0,1,0,lb:]), label=label  + " atom " + str(plotatoms[i]), lw=lw, marker='o', c=cm(norm(normval)))
					maingiw[1,1].plot(iwaxis[lb:], np.imag(giw[i,1,1,1,1,lb:]), label=label  + " atom " + str(plotatoms[i]), lw=lw, marker='o', c=cm(norm(normval)))
					maingiw[2,0].plot(iwaxis[lb:], np.imag(giw[i,2,0,2,0,lb:]), label=label  + " atom " + str(plotatoms[i]), lw=lw, marker='o', c=cm(norm(normval)))
					maingiw[2,1].plot(iwaxis[lb:], np.imag(giw[i,2,1,2,1,lb:]), label=label  + " atom " + str(plotatoms[i]), lw=lw, marker='o', c=cm(norm(normval)))
			else:
				for i in range(0,len(plotatoms)):
					cm = cms[i%5]
					lw = 3 - i*0.5
					mainsiw[0].plot(iwaxis[lb:], np.imag(siw[i,0,spin,0,spin,lb:]), label=label  + " atom " + str(plotatoms[i]), lw=lw, marker='o', c=cm(norm(normval)))
					mainsiw[1].plot(iwaxis[lb:], np.imag(siw[i,1,spin,1,spin,lb:]), label=label  + " atom " + str(plotatoms[i]), lw=lw, marker='o', c=cm(norm(normval)))
					mainsiw[2].plot(iwaxis[lb:], np.imag(siw[i,2,spin,2,spin,lb:]), label=label  + " atom " + str(plotatoms[i]), lw=lw, marker='o', c=cm(norm(normval)))
				for i in range(0,len(plotatoms)):
					cm = cms[i%5]
					lw = 3 - i*0.5
					maingiw[0].plot(iwaxis[lb:], np.imag(giw[i,0,spin,0,spin,lb:]), label=label  + " atom " + str(plotatoms[i]), lw=lw, marker='o', c=cm(norm(normval)))
					maingiw[1].plot(iwaxis[lb:], np.imag(giw[i,1,spin,1,spin,lb:]), label=label  + " atom " + str(plotatoms[i]), lw=lw, marker='o', c=cm(norm(normval)))
					maingiw[2].plot(iwaxis[lb:], np.imag(giw[i,2,spin,2,spin,lb:]), label=label  + " atom " + str(plotatoms[i]), lw=lw, marker='o', c=cm(norm(normval)))
			for i in range(0,len(plotatoms)):
				cm = cms[i%5]
				lw = 3 - i*0.5
				maintrgiw.plot(iwaxis[lb:], np.imag(trgiw[i,lb:]), label=label  + " atom " + str(plotatoms[i]), lw=lw, marker='o', c=cm(norm(normval)))
	if foundfile:
		print "Found at least one file for this Hamiltonian. Generating output files"

		norm2 = mcolors.Normalize(0,Nbands)
		if plotbetas == True:
			orderlist = np.argsort(betalist,axis=0)
			ordplotlist = np.sort(betalist, axis=0)
			ordat1fullocc = np.zeros_like(at1fullocc)
			ordfullsigns = np.zeros_like(fullsigns)
			for i in range(0,len(ordplotlist)):
				ordat1fullocc[i] = at1fullocc[orderlist[i]]
				ordfullsigns[i,:] = fullsigns[orderlist[i],:]
		elif plotlambdas == True:
			orderlist = np.argsort(lambdalist,axis=0)
			ordplotlist = np.sort(lambdalist, axis=0)
			if xi == True:
				ordplotlist = ordplotlist*2
			ordat1fullocc = np.zeros_like(at1fullocc)
			ordfullsigns = np.zeros_like(fullsigns)
			for i in range(0,len(ordplotlist)):
				ordat1fullocc[i] = at1fullocc[orderlist[i]]
				ordfullsigns[i,:] = fullsigns[orderlist[i],:]
		markers = ['s', 'v', '^','x','plus','point']
		
		for i in range(0,len(plotatoms)):
			cm = cms[6]
			mainsign.plot(ordplotlist[:], abs(np.real(ordfullsigns[:,i])), label="Atom " + str(plotatoms[i]), lw=3, marker=markers[i], markersize=8, c=cm(norm2(i)))
		if spin==2:
			for i in range(0,Nbands):
				cm = cms[6]
				lw = 3
				mainocc[0,0].plot(ordplotlist[:], ordat1fullocc[:,i,0,i,0], label="orbital " + str(i+1), lw=lw, marker=markers[i], markersize=8, c=cm(norm2(i)))
				mainocc[0,1].plot(ordplotlist[:], ordat1fullocc[:,i,1,i,1], label="orbital " + str(i+1), lw=lw, marker=markers[i], markersize=8, c=cm(norm2(i)))
				mainocc[1,0].plot(ordplotlist[:], ordat1fullocc[:,i,0,i,1], label="orbital " + str(i+1), lw=lw, marker=markers[i], markersize=8, c=cm(norm2(i)))
				mainocc[1,1].plot(ordplotlist[:], ordat1fullocc[:,i,0,i,1]/(ordat1fullocc[:,i,0,i,0]*ordat1fullocc[:,i,1,i,1]), label="orbital " + str(i+1), lw=lw, marker=markers[i], markersize=8, c=cm(norm2(i)))
				mainocc[2,0].plot(ordplotlist[:], ordat1fullocc[:,i,0,(i+1)%Nbands,0], label=str(i+1) + "," + str((i+1)%Nbands+1), lw=lw, marker=markers[i], markersize=8, c=cm(norm2(i)))
				mainocc[2,1].plot(ordplotlist[:], ordat1fullocc[:,i,1,(i+1)%Nbands,1], label=str(i+1) + "," + str((i+1)%Nbands+1), lw=lw, marker=markers[i], markersize=8, c=cm(norm2(i)))
				mainocc[3,0].plot(ordplotlist[:], ordat1fullocc[:,i,0,(i+1)%Nbands,0]/(ordat1fullocc[:,i,0,i,0]*ordat1fullocc[:,(i+1)%Nbands,0,(i+1)%Nbands,0]),\
						  label=str(i+1) + "," + str((i+1)%Nbands+1), lw=lw, marker=markers[i], markersize=8, c=cm(norm2(i)))
				mainocc[3,1].plot(ordplotlist[:], ordat1fullocc[:,i,1,(i+1)%Nbands,1]/(ordat1fullocc[:,i,1,i,1]*ordat1fullocc[:,(i+1)%Nbands,1,(i+1)%Nbands,1]),\
						  label=str(i+1) + "," + str((i+1)%Nbands+1), lw=lw, marker=markers[i], markersize=8, c=cm(norm2(i)))
		else:
			for i in range(0,Nbands):
				cm = cms[6]
				lw = 3
				mainocc[0].plot(ordplotlist[:], ordat1fullocc[:,i,0,i,0], label="orbital " + str(i+1), lw=lw, marker=markers[i], markersize=8, c=cm(norm2(i)))
				mainocc[1].plot(ordplotlist[:], ordat1fullocc[:,i,0,i,1], label="orbital " + str(i+1), lw=lw, marker=markers[i], markersize=8, c=cm(norm2(i)))
				mainocc[2].plot(ordplotlist[:], ordat1fullocc[:,i,0,(i+1)%Nbands,0], label=str(i+1) + "," + str((i+1)%Nbands+1), lw=lw, marker=markers[i], markersize=8, c=cm(norm2(i)))
				mainocc[3].plot(ordplotlist[:], ordat1fullocc[:,i,0,(i+1)%Nbands,0]/(ordat1fullocc[:,i,0,i,0]*ordat1fullocc[:,(i+1)%Nbands,0,(i+1)%Nbands,0]),\
						label=str(i+1) + "," + str((i+1)%Nbands+1), lw=lw, marker=markers[i], markersize=8, c=cm(norm2(i)))

		if spin==2:
			for i in range(0,Nbands):
				mainsiw[i,0].legend(bbox_to_anchor=(1,0), loc='lower right')
				mainsiw[i,1].legend(bbox_to_anchor=(1,0), loc='lower right')
			for i in range(0,Nbands):
				maingiw[i,0].legend(bbox_to_anchor=(1,0), loc='lower right')
				maingiw[i,1].legend(bbox_to_anchor=(1,0), loc='lower right')
		else:
			for i in range(0,Nbands):
				mainsiw[i].legend(bbox_to_anchor=(1,0), loc='lower right')
			for i in range(0,Nbands):
				maingiw[i].legend(bbox_to_anchor=(1,0), loc='lower right')

		maintrgiw.legend(bbox_to_anchor=(1,0), loc='lower right')
		mainsign.legend(bbox_to_anchor=(1,1), loc='upper right')

		spanoccplot = ordplotlist[-1]-ordplotlist[0]
		offset = spanoccplot*0.1
		mainsign.set_xlim([ordplotlist[0]-offset,ordplotlist[-1]+offset])
		if spin==2:
			mainocc[0,0].set_xlim([ordplotlist[0]-offset,ordplotlist[-1]+offset])
			mainocc[0,1].set_xlim([ordplotlist[0]-offset,ordplotlist[-1]+offset])
			mainocc[1,0].set_xlim([ordplotlist[0]-offset,ordplotlist[-1]+offset])
			mainocc[1,1].set_xlim([ordplotlist[0]-offset,ordplotlist[-1]+offset])
			mainocc[2,0].set_xlim([ordplotlist[0]-offset,ordplotlist[-1]+offset])
			mainocc[2,1].set_xlim([ordplotlist[0]-offset,ordplotlist[-1]+offset])
			mainocc[3,0].set_xlim([ordplotlist[0]-offset,ordplotlist[-1]+offset])
			mainocc[3,1].set_xlim([ordplotlist[0]-offset,ordplotlist[-1]+offset])
			mainocc[0,0].legend(bbox_to_anchor=(1,0.7), loc='upper right')
			mainocc[0,1].legend(bbox_to_anchor=(1,0.7), loc='upper right')
			mainocc[1,0].legend(bbox_to_anchor=(1,0.7), loc='upper right')
			mainocc[1,1].legend(bbox_to_anchor=(1,0.7), loc='upper right')
			mainocc[2,0].legend(bbox_to_anchor=(1,0.7), loc='upper right', title='orb.1,orb.2')
			mainocc[2,1].legend(bbox_to_anchor=(1,0.7), loc='upper right', title='orb.1,orb.2')
			mainocc[3,0].legend(bbox_to_anchor=(1,0.7), loc='upper right', title='orb.1,orb.2')
			mainocc[3,1].legend(bbox_to_anchor=(1,0.7), loc='upper right', title='orb.1,orb.2')
		else:
			mainocc[0].set_xlim([ordplotlist[0]-offset,ordplotlist[-1]+offset])
			mainocc[1].set_xlim([ordplotlist[0]-offset,ordplotlist[-1]+offset])
			mainocc[2].set_xlim([ordplotlist[0]-offset,ordplotlist[-1]+offset])
			mainocc[3].set_xlim([ordplotlist[0]-offset,ordplotlist[-1]+offset])
			mainocc[0].legend(bbox_to_anchor=(1,0.7), loc='upper right')
			mainocc[1].legend(bbox_to_anchor=(1,0.7), loc='upper right')
			mainocc[2].legend(bbox_to_anchor=(1,0.7), loc='upper right', title='orb.1,orb.2')
			mainocc[3].legend(bbox_to_anchor=(1,0.7), loc='upper right', title='orb.1,orb.2')

			
		if svg:
			filenamesiw = basetitle + "_" + hamiltonian + "_siw.svg"
			filenamegiw = basetitle + "_" + hamiltonian + "_giw.svg"
			filenametrgiw = basetitle + "_" + hamiltonian + "_trgiw.svg"
			filenameocc = basetitle + "_" + hamiltonian + "_occ.svg"
			filenamesign = basetitle + "_" + hamiltonian + "_sign.svg"
		else:
			filenamesiw = basetitle + "_" + hamiltonian + "_siw.pdf"
			filenamegiw = basetitle + "_" + hamiltonian + "_giw.pdf"
			filenametrgiw = basetitle + "_" + hamiltonian + "_trgiw.pdf"
			filenameocc = basetitle + "_" + hamiltonian + "_occ.pdf"
			filenamesign = basetitle + "_" + hamiltonian + "_sign.pdf"
		
		figsiw.savefig(filenamesiw)
		figgiw.savefig(filenamegiw)
		figtrgiw.savefig(filenametrgiw)
		figoccn.savefig(filenameocc)
		figsign.savefig(filenamesign)
	else:
		print "Found no file for this Hamiltonian. No output generated."
				
#add argparse to give startpath, optionals to specify beta and hamiltonians and one string to append to the title, specify xrange
#plot specific components of (siw, giw), with options to plot all,  Tr(giw), plot occ
