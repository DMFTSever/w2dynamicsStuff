import numpy as np
from numpy.linalg import inv
import scipy.optimize as opt
from scipy.optimize import curve_fit
import argparse

parser = argparse.ArgumentParser(description="This script builds the spectral function from the Hamiltonian and self\
					      energy")
parser.add_argument('mu', help="Chemical potential of calculation", type=float)
parser.add_argument('dc', help="Double Counting value", type=float)
parser.add_argument('smom11', help="Self Energy Moment 1 of band 1", type=float)
parser.add_argument('smom12', help="Self Energy Moment 1 of band 2", type=float)
parser.add_argument('smom13', help="Self Energy Moment 1 of band 3", type=float)
parser.add_argument('hkfile', help="hkfile", type=str)

args=parser.parse_args()

mu=args.mu
dc=args.dc
smom11=args.smom11
smom12=args.smom12
smom13=args.smom13
hkfile=args.hkfile

# *********************# 
# Reading from Files 1
nfreq=0
iwn = []
Imsiw111 = []
f = open('11/ImSiw.dat',"r")
for line in f:
  columns = line.split()
  iwn_c = float(columns[0])
  ISiw_c = float(columns[1])
  iwn.append(iwn_c)
  Imsiw111.append(ISiw_c)
  nfreq+=1
f.close  
###
nfreq=0
Resiw111 = []
iwn_kk = []
f = open('11/ReSiw.dat',"r")
for line in f:
  columns = line.split()
  iwn_c = float(columns[0])
  RSiw_c = float(columns[1])
  iwn_kk.append(iwn_c)
  Resiw111.append(RSiw_c)
  nfreq+=1
f.close  
###

# *********************# 
# Reading from Files 2
nfreq=0
iwn = []
Imsiw121 = []
f = open('21/ImSiw.dat',"r")
for line in f: 
  columns = line.split() 
  iwn_c = float(columns[0]) 
  ISiw_c = float(columns[1]) 
  iwn.append(iwn_c)
  Imsiw121.append(ISiw_c)
  nfreq+=1
f.close  
###
nfreq=0
Resiw121 = []
iwn_kk = []
f = open('21/ReSiw.dat',"r")
for line in f:
  columns = line.split()
  iwn_c = float(columns[0])
  RSiw_c = float(columns[1])
  iwn_kk.append(iwn_c)
  Resiw121.append(RSiw_c)
  nfreq+=1
f.close  
###

# *********************# 
# Reading from Files 3
nfreq=0
iwn = []
Imsiw131 = []
f = open('31/ImSiw.dat',"r")
for line in f:
  columns = line.split()
  iwn_c = float(columns[0])
  ISiw_c = float(columns[1])
  iwn.append(iwn_c)
  Imsiw131.append(ISiw_c)
  nfreq+=1
f.close  
###
nfreq=0
Resiw131 = []
iwn_kk = []
f = open('31/ReSiw.dat',"r")
for line in f:
  columns = line.split()
  iwn_c = float(columns[0])
  RSiw_c = float(columns[1])
  iwn_kk.append(iwn_c)
  Resiw131.append(RSiw_c)
  nfreq+=1
f.close  
###

#### *********************# 
#### Reading from Files 4
###nfreq=0
###iwn = []
###Imsiw141 = []
###f = open('ImSigma141.dat',"r")
###for line in f:
  ###columns = line.split()
  ###iwn_c = float(columns[0])
  ###ISiw_c = float(columns[1])
  ###iwn.append(iwn_c)
  ###Imsiw141.append(ISiw_c)
  ###nfreq+=1
###f.close  
######
###nfreq=0
###Resiw141 = []
###iwn_kk = []
###f = open('ReSigma141.dat',"r")
###for line in f:
  ###columns = line.split()
  ###iwn_c = float(columns[0])
  ###RSiw_c = float(columns[1])
  ###iwn_kk.append(iwn_c)
  ###Resiw141.append(RSiw_c)
  ###nfreq+=1
###f.close  
######

#### *********************# 
#### Reading from Files 5
###nfreq=0
###iwn = []
###Imsiw151 = []
###f = open('ImSigma151.dat',"r")
###for line in f:
  ###columns = line.split()
  ###iwn_c = float(columns[0])
  ###ISiw_c = float(columns[1])
  ###iwn.append(iwn_c)
  ###Imsiw151.append(ISiw_c)
  ###nfreq+=1
###f.close  
######
###nfreq=0
###Resiw151 = []
###iwn_kk = []
###f = open('ReSigma151.dat',"r")
###for line in f:
  ###columns = line.split()
  ###iwn_c = float(columns[0])
  ###RSiw_c = float(columns[1])
  ###iwn_kk.append(iwn_c)
  ###Resiw151.append(RSiw_c)
  ###nfreq+=1
###f.close  
######

# *********************# 
# Reading Hamiltonian
f = open(hkfile,'r')
line1 = f.readline()
kpoints = int(line1.split()[0])
nbands = int(line1.split()[1])
Hk_read = np.zeros((kpoints,nbands,nbands),dtype=complex)
for k in range(0,kpoints):
  line = f.readline()
  for b1 in range(0,nbands):
    line = f.readline().split()
    for b2 in range(0,nbands):
      Hk_read[k,b1,b2] = float(line[2*b2+0])+ 1j*float(line[2*b2+1])
f.close
#Hk=Hk_read[:,0:11,0:11]
Hk=Hk_read

# ****************************************************************************** #
# ****************************************************************************** #
# ****************************************************************************** #


Resiw111 = np.asarray(Resiw111)
Imsiw111 = np.asarray(Imsiw111)
Resiw121 = np.asarray(Resiw121)
Imsiw121 = np.asarray(Imsiw121)
Resiw131 = np.asarray(Resiw131)
Imsiw131 = np.asarray(Imsiw131)

# KramerKronig has fixed output => Interpolate to iwn grid.
Resiw_int111 = np.interp(iwn,iwn_kk,Resiw111)
Resiw_int121 = np.interp(iwn,iwn_kk,Resiw121)
Resiw_int131 = np.interp(iwn,iwn_kk,Resiw131)

iwn = np.asarray(iwn)
iwn_kk = np.asarray(iwn_kk)

eye = np.eye((Hk.shape[1]))
# Make complex Siw field
d0 = int(iwn.shape[0])
d1 = int(Hk.shape[1])
d2 = int(Hk.shape[2])
Siw = np.zeros((d0,d1,d2),dtype=complex)


Siw[:,0,0] = Resiw_int111 + 1j*Imsiw111 + smom11 + dc   # The last value is the ReShift.
Siw[:,1,1] = Resiw_int121 + 1j*Imsiw121 + smom12 + dc
Siw[:,2,2] = Resiw_int131 + 1j*Imsiw131 + smom13+ dc

Siw[:,3,3] = Resiw_int111 + 1j*Imsiw111 + smom11 + dc   # The last value is the ReShift.
Siw[:,4,4] = Resiw_int121 + 1j*Imsiw121 + smom12 + dc
Siw[:,5,5] = Resiw_int131 + 1j*Imsiw131 + smom13+ dc

Siw[:,6,6] = Resiw_int111 + 1j*Imsiw111 + smom11 + dc   # The last value is the ReShift.
Siw[:,7,7] = Resiw_int121 + 1j*Imsiw121 + smom12 + dc
Siw[:,8,8] = Resiw_int131 + 1j*Imsiw131 + smom13+ dc

Siw[:,9,9] = Resiw_int111 + 1j*Imsiw111 + smom11 + dc   # The last value is the ReShift.
Siw[:,10,10] = Resiw_int121 + 1j*Imsiw121 + smom12 + dc
Siw[:,11,11] = Resiw_int131 + 1j*Imsiw131 + smom13+ dc

Siw[:,12,12] = Resiw_int111 + 1j*Imsiw111 + smom11 + dc   # The last value is the ReShift.
Siw[:,13,13] = Resiw_int121 + 1j*Imsiw121 + smom12 + dc
Siw[:,14,14] = Resiw_int131 + 1j*Imsiw131 + smom13+ dc

Siw[:,15,15] = Resiw_int111 + 1j*Imsiw111 + smom11 + dc   # The last value is the ReShift.
Siw[:,16,16] = Resiw_int121 + 1j*Imsiw121 + smom12 + dc
Siw[:,17,17] = Resiw_int131 + 1j*Imsiw131 + smom13+ dc

#tiny = -1j*0.1
#Siw[:,5,5] = tiny
#Siw[:,6,6] = tiny
#Siw[:,7,7] = tiny
#Siw[:,8,8] = tiny
#Siw[:,9,9] = tiny
#Siw[:,10,10] = tiny
#Siw[:,16,16] = tiny
#Siw[:,17,17] = tiny
#Siw[:,18,18] = tiny
#Siw[:,19,19] = tiny
#Siw[:,20,20] = tiny
#Siw[:,21,21] = tiny

### Inversion performed for each k-point to obtain k-dependent spectral function Ak
##k = 0
##Ak = np.zeros((Hk.shape[0],iwn.shape[0],Hk.shape[1],Hk.shape[2]),dtype=complex)
##exit()
##for Hk_idx in Hk:
  ##Ziw = (iwn[:,None,None] + mu)*eye - Siw - Hk_idx
  ##Ak[k,:,:,:] = inv(Ziw)
  ##k+=1
### Obtain local spectral function A
##A = np.zeros((iwn.shape[0],Hk.shape[1],Hk.shape[2]),dtype=complex)
##A = sum(Ak)/Hk.shape[0]/np.pi

# Inversion performed for each k-point to obtain k-dependent spectral function Ak
k = 0
Ak = np.zeros((iwn.shape[0],Hk.shape[1],Hk.shape[2]),dtype=complex)
for Hk_idx in Hk:
  if np.mod(k,100)==0:
    print k
  Ziw = (iwn[:,None,None] + mu)*eye - Siw - Hk_idx
  Ak[:,:,:] = Ak[:,:,:] + inv(Ziw)
  k+=1
# Obtain local spectral function A
A = np.zeros((iwn.shape[0],Hk.shape[1],Hk.shape[2]),dtype=complex)
A = Ak/Hk.shape[0]/np.pi


A_verage = np.zeros((iwn.shape[0]),dtype=complex)


# Write into orbital files
bb = 0
for i1 in range(0,Hk.shape[1]):
  for i2 in range(0,Hk.shape[2]):
    fname = "A_"
    fname += str(i1)
    fname += "_"
    fname += str(i2)
    print 'Writing:',fname
    fname += ".dat"

    if i1==i2:
      f = open(fname,'w')
      for j in range(0,iwn.shape[0]):
        f.write( str(iwn[j]) + ' ' + str(-np.imag(A[j,i1,i2]))  + ' ' + str(np.real(A[j,i1,i2])) + '\n')  
      f.close

    if i1==i2:
      A_verage = A_verage + A[:,i1,i2]
      bb+=1

A_verage=A_verage/bb
f = open('A_verage.dat','w')
for j in range(0,iwn.shape[0]):
  f.write( str(iwn[j]) + ' ' + str(-np.imag(A_verage[j]))  + ' ' + str(np.real(A_verage[j])) + '\n')  
f.close

print 'Success'
