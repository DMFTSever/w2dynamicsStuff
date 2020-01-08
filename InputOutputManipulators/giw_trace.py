#################################################
# Tracing of the giw_local of a given hdf5 file #
#################################################

##########################################################################################
#usage: python giw_trace.py <hdf5 file to use> <number of atom for which the trace should be calculated> <iteration for which to perform trace, default: last>
##########################################################################################
import numpy as np
import h5py
import sys 
auxdir="/dss/dsshome1/0D/di76rir/lib/MyPythonScripts/W2dynamics/General"
sys.path.insert(0,auxdir)
import hdf5_stuff as h5s
import custom_errors as err
import readwrite as rw

try:
    if len(sys.argv) < 3: raise err.InputError("Error reading input.", "Not enough arguments found.")
    filename = sys.argv[1]
    at = int(sys.argv[2])
    if at > 999: raise ValueError
    at = "%03i"%at
    if filename[-5:] != ".hdf5": raise err.InputError("Error reading input.", "No hdf5 file found as second argument")
except err.InputError as error:
    print error.expression
    print error.message
    print "usage: python giw_trace.py <hdf5 file to use> <number of atom for which the trace should be calculated> <iteration for which to perform trace, default: last>"
    sys.exit()
except ValueError:
    print "Error assigning atom number to " + sys.argv[2]
    print "Must be integer smaller 999"
except:
    print "Some unkown error occourd during opening input file"
    sys.exit()

try:
    f = h5py.File(filename, 'r')
except:
    print "file " + filename + "  cannot be opened."
    sys.exit()

print "Extracing some config parameters"
Nat = f["/.config"].attrs.get('general.nat')

try: 
    if int(at) > Nat: raise err.InputError("Error getting number of d bands", "Atom to use not in hdf5 file")
    Ndbands =f["/.config"].attrs.get(("atoms."+ str(int(at)) +".nd"))
except err.InputError as error:
    print error.expression
    print error.message
    sys.exit()
except:
    print "Unkown error while extraction number of d bands."
    sys.exit()

try:
    if sys.argv[3] == "last":
        iteration = str(sys.argv[3])
    elif len(sys.argv[3]) == 1:
        iteration = "00" + sys.arv[3]
        int(iteration)
    elif len(sys.argv[3]) == 2:
        iteration = "0" + sys.argv[3]
        int(iteration)
    elif len(sys.argv[3]) == 3:
        iteration == sys.argv[3]
        int(iteration)
    else: 
        raise err.InputError("Error extracting giw of iteration" + sys.argv[3], "Not a valid iteration number")
except err.InputError as error:
    print error.expression
    print error.message
    sys.exit()
except ValueError:
    print "Error extractin giw of iteration" + iteration
    print "Not a valid iteration number. Needs to be integer."
    sys.exit()
except:
    print "No iteration given. Defaulting to last"
    iteration = "last"

print "Extracting giw-full of atom " + at + " of " + iteration + "iteration from " + filename
try:
    giw = f["/dmft-" + iteration + "/ineq-" + at + "/giw-full/value"][:,:,:,:,:]
except:
    print "Cannot extract giw-full of atom " + at + " of iteration " + iteration + " from " + filename
    sys.exit()
print "giw.shape ", giw.shape

print "Extracting iw axis from " + filename
iw = np.array(f["/.axes/iw"])
print "iw.shape", iw.shape
iwpoints = iw.shape[0]

print "Calculating Trace"
giw = giw.reshape(2*Ndbands,2*Ndbands,iwpoints)
print "giw.shape ", giw.shape
giw_trace = giw.trace()
print "giw_trace.shape ", giw_trace.shape

print "Writting giw_trace to file: " + "giw_atom_"+at+"_iteration_"+iteration+"_trace.dat"
rw.write_function("giw_atom_"+at+"_iteration_"+iteration+"_trace.dat", iw, "iw", giw_trace, "Tr(Giw)", "This file contains the trace of giw-full of "+filename)

f.close()
