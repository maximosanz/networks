import re
import numpy as np
import random as rd
import networkx as nx
import loadFunctions as loadf
import pickle
import sys


Nparam = 0

# Default values
Np = 500
percentile=10
Estart = 1.0
Eend = 0.1
restartF = ''
outputFilename = ''
transform = False

# Command line options

leave=False
if len(sys.argv)>1:
	for opt in sys.argv[1::2]:
		if not re.match('--',opt):
			leave =True
			break
		else:
			opt = opt[2:]
		val = sys.argv[sys.argv.index('--'+opt)+1]
		if opt == 'nparam':
			Nparam = int(val)
		elif opt == 'restart':
			restartF = val
		elif opt == 'npart':
			Np = int(val)
		elif opt == 'percentile':
			percentile = int(val)
		elif opt == 'Estart':
			Estart = float(val)
		elif opt == 'Eend':
			Eend = float(val)
		elif opt == 'out':
			outputFilename = val
		elif opt == 'transform':
			transform = True
else:
	leave=True
if not Nparam or not outputFilename:
	leave =True
if leave:
	print ('\nUsage:'
	'\n\nRequired arguments:'
	'\n\t--nparam Number of parameters'
	'\n\t--out Output file name'
	'\n\nOptional arguments (default values):'
	'\n\t--npart Number of particles (500)'
	'\n\t--percentile Percentile in the progression of auomated tolerance (10)'
	'\n\t--Estart Initial tolerance value (1.0)'
	'\n\t--Eend Final tolerance value (0.1)'
	'\n\t--restart Name of pickle file name for previous run (no restart)'
	'\n\t--transform Your generative model should produce graphs with node IDs from 0 to N (graph length). Run this option if that is not the case (slightly slower)'
	'\n\tNote: nparam, npart cannot change when restarting.'
	)
	sys.exit()

if restartF:
	results = pickle.load(open(restartF))
else:
    results = []

prior = np.empty((Nparam,2))

# Parameter ranges: all from 0 to 1

prior[:,0].fill(0) #lower bound
prior[:,1].fill(1) #upper bound

##########################
# Include here the graph you wish to compare to, the GDDs will be computed.
Datagraph = nx.read_edgelist('homo_domainppiedges.txt')
Data = loadf.get_GDDs(Datagraph,transform=True)

#######################

results = loadf.abcsmc_graphs(Data,Np,Estart,Eend,Nparam,prior,outputFilename=outputFilename,percentile=percentile,results=results,transform=transform)
