import cython
cimport corca
import numpy as np
cimport numpy as np

ctypedef long long int64

def python_getorbits(np.ndarray[int, ndim=2, mode="c"] input not None, int nodenumber,int edgenumber):
	orbits = np.zeros((nodenumber*73), dtype=np.longlong)
	cdef np.ndarray[int64, ndim=1, mode = 'c'] np_buff = np.ascontiguousarray(orbits, dtype = np.longlong)
	corca.getorbits(<int*>input.data,nodenumber,edgenumber,<int64*>np_buff.data)
	orbits = orbits.reshape(nodenumber,73)
	return orbits

