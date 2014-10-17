from cpython cimport bool
import numpy as np
cimport numpy as np

cdef extern from "corca.h":
	ctypedef long long int64
	void getorbits(int *edgelist, int nodenumber, int edgenumber, int64 *my_orbit)
	void count5()
	void writeResults()



