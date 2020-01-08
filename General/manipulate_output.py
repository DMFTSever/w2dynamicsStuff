import numpy as np
import sys
import readwrite
import custom_errors as err		

def trace_giw(giw)
	''' Routine to calculate the Trace of a local Green's function of shape (bands,2,bands,2)'''

	try:
		if type(giw) is not np.ndarray: raise TypeError("Input argument needs to be of type numpy.ndarray")
		if len(giw.shape) != 5: raise err.ShapeError("Input argument needs to be of shape (bands,2,bands,2,iwpoints)")
		if giw.shape[1] != 2 or giw.shape[3] !- 2: raise err.ShapeError("Input argument needs to be of shape (bands,2,bands,2,iwpoints)")
		if giw.shape[0] != giw.shape[2]: raise err.ShapeError("Input argument needs to be of shape (bands,2,bands,2,iwpoints)")
	except TypeError:
		raise
		sys.exit()
	except err.ShapeError:
		raise
		sys.exit()

	Nbands = giw.shape[0]
	Niwpoints = giw.shape[4]
	giw = giw.reshape(2*Nbands,2*Nbands,Niwpoints)
	return giw.trace()
