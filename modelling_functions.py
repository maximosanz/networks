import sys
import random as rd
import re
import networkx as nx
import os
import copy
import subprocess
import numpy as np
import orca
import datetime
import scipy
import scipy.spatial

# Returns the seedgraph and the compatibility information of each superfamily, which will be extended if more superfamilies come in.

def create_seed(G,loops=True,ringSize=5):
	supfamlist = []
	compatibilitylist = []
	if not G:
		seedgraph = nx.cycle_graph(ringSize)
	else:
		if not loops:
			for n in G.nodes():
				neigh = G.neighbors(n)
				if neigh == [n]:
					G.remove_node(n)
				elif n in neigh:
					G.remove_edge(n,n)
		initialnodenumber_i=len(G)

		seedgraph=nx.Graph()
		seedgraph.add_nodes_from(range(initialnodenumber_i))
	
		#adding 'superfamily' attribute to each of the nodes.
		for i in (range(initialnodenumber_i)):
			seedgraph.node[i]['family']=G.nodes()[i] 
		#adding edges to all the nodes
		for eachnode_i in seedgraph.nodes():
			nodefamily_s=seedgraph.node[eachnode_i]['family'] #determining superfamily of the node
			compatiblefamily_l=G.edge[nodefamily_s].keys()#determining what that superfamily is compatible with.
			for eachcompatiblefamily_s in compatiblefamily_l: # for each compatiblefamily of the node

				for targetnode_i,targetnodefamily_d in seedgraph.node.items(): #checking which nodes have the compatible family
					if (eachcompatiblefamily_s==targetnodefamily_d['family']):
						seedgraph.add_edge(eachnode_i,targetnode_i)
	i=0
	while i<len(seedgraph):
		supfamlist.append(set([i]))
		seedgraph.node[i]['family'] = i
		compatibilitylist.append(seedgraph.neighbors(i))
		i += 1
	seed = [seedgraph,supfamlist,compatibilitylist]
	return seed
	   
def export_edgelist(G,o,as_numbers=False,noLoops=False,indexBase=0,orca=False):
	if noLoops:
		for node in G.nodes():
			neig = G.neighbors(node)
			if G.degree(node) == 0 or neig==[node]:
				G.remove_node(node)
			elif node in neig:
				G.remove_edge(node,node)
	if orca: #First line is the size of graph and number of edges
		o.write(str(len(G))+' '+str(len(G.edges()))+'\n')
	nodelist = G.nodes()
	for edge in G.edges():
		n1 = edge[0]
		n2 = edge[1]
		if as_numbers:
			n1 = nodelist.index(n1)+indexBase
			n2 = nodelist.index(n2)+indexBase
		if not noLoops or n1 != n2:
			o.write(str(n1)+'\t'+str(n2)+'\n')
	return
	
def dd_model(parameters,seed,desiredfinalnodecount_i,dda=False,domain=False,loops=False,denovoDomain=False,nodemax=0,edgemax=0,orphanmax=0):
	if denovoDomain and not domain:
		print 'Error: Incompatible options. "denovoDomain" set to True but "domain" is False'
		return None
	if domain and not dda:
		print 'Error: Incompatible options. "domain" set to True but "dda" is False'
		return None

	# Determining the value of each parameter according to the type of model selected
	# Parameters should be provided in the following order: [pLoss,pGain,pNewSuperfamily,pParentChild]
	# Only required parameters for the selected type of model should be included. e.g. do not include pGain if dda is False
	npar = 1
	pLoss = parameters[0]
	if dda:
		pGain=parameters[1]
		npar += 1
		if denovoDomain:
			pNewSuperfamily=parameters[2]
			npar += 1
			if loops:
				pNewLoop=parameters[3]
				npar += 1
	if not loops:
		pParentChild=parameters[-1]
		npar += 1
	if len(parameters) != npar:
		print 'Error: Incorrect number of parameters provided'
		return None

	# Retrieving information from the seed
	seedgraph = seed[0]
	if domain:
		supfamlist = copy.deepcopy(seed[1])
		compatibilitylist = copy.deepcopy(seed[2])

	initialnodenumber_i=len(seedgraph)
	supfam_number = initialnodenumber_i
	
	graph1=seedgraph.copy()
	nodenumber_i=initialnodenumber_i - 1
	# x keeps track of how many nodes are orphan or contain only a self loop (not considered for the final count)
	x=graph1.degree().values().count(0)
	for node in graph1.nodes():
		if graph1.neighbors(node) == [node]:
			x += 1
	edgeno = len(graph1.edges())
	while nodenumber_i-x < desiredfinalnodecount_i:
		nodenumber_i+=1 

		#When nodemax, edgemax or orphanmax are provided, the function returns None if the limit is reached

		if nodemax and nodenumber_i > nodemax:
			return None
		if edgemax and edgeno > edgemax:
			return None
		if orphanmax and x > orphanmax:
			return None
		rndpick_i=rd.choice(graph1.nodes()) 
		attachtoOrphan = False
		attachtoSelforDad = False
		selfconnected = False

		if graph1.degree(rndpick_i) == 0 or graph1.neighbors(rndpick_i) == [rndpick_i]:
			orphanDad = True
		else:
			orphanDad = False
		graph1.add_node(nodenumber_i)

		# Modelling the de-novo generation of new superfamilies - Only if denovoDomain is True and its probability is satisfied
		#-------------------------------------------------------------------------------
		if denovoDomain and rd.random() < pNewSuperfamily:
			supfamlist.append(set([nodenumber_i]))
			nodefamily = graph1.node[rndpick_i]['family']
			graph1.node[nodenumber_i]['family']=supfam_number #A new superfamily appears
			compatibilitylist.append([nodefamily])
			compatibilitylist[nodefamily].append(supfam_number)
			graph1.add_edge(nodenumber_i,rndpick_i)
			edgeno += 1
			if loops and rd.random() < pNewLoop:
				graph1.add_edge(nodenumber_i,nodenumber_i)
				edgeno += 1
				compatibilitylist[supfam_number].append(supfam_number)
			supfam_number += 1
		else:

			# Duplication step 
			#-------------------------------------------------------------------------------
			if domain:
				nodefamily = graph1.node[rndpick_i]['family']
				supfamlist[nodefamily].update({nodenumber_i})
				graph1.node[nodenumber_i]['family']=nodefamily
			neig_l=graph1.neighbors(rndpick_i)
		
			for i_i in neig_l:
				graph1.add_edge(nodenumber_i,i_i)
				edgeno += 1
				if i_i == rndpick_i:
					graph1.add_edge(nodenumber_i,nodenumber_i)
					edgeno += 1
					selfconnected = True
					if not loops:
						print 'Error: Loops are not allowed. If you wish to incorporate loops set the variable "loops" to True'
						return None


			# Determining if should add interaction with parent
			#------------------------------------------------------------------------
			if not loops:
				a=rd.random() 
				if a<pParentChild:
					graph1.add_edge(nodenumber_i,rndpick_i)
					edgeno += 1
			elif selfconnected:
				# Edge between parent and child is there (as well as child self-loop).
				# Out of both self loops and the parent-child edge, two may disappear
			
				threedges = [(nodenumber_i,rndpick_i),(nodenumber_i,nodenumber_i),(rndpick_i,rndpick_i)]
				threedges.remove(rd.choice(threedges))
				for edge in threedges:
					b=rd.random()
					if b<pLoss:
						graph1.remove_edge(*edge)
						edgeno -= 1
				neig_l.remove(rndpick_i)
	
			# Divergence Step
			#--------------------------------------------------------------

			for i_i in neig_l: 
				b=rd.random()
				if b<pLoss:
					rndedge= rd.choice([(rndpick_i,i_i),(nodenumber_i,i_i)])
					graph1.remove_edge(*rndedge) 
					edgeno -= 1
		
			# Attachment Step
			#------------------------------------------------------------- This section doesn't exist in the DD model.

			if dda and rd.random()<pGain:
				attachfrom=rd.choice([nodenumber_i,rndpick_i])
				if not domain:
					pool = set(graph1.nodes())
				else:
					compatible_families = compatibilitylist[nodefamily]
					pool = set()
					for f in compatible_families:
						pool = pool.union(supfamlist[f])
				if not loops and attachfrom in pool:
					pool.remove(attachfrom)
				pool = pool.difference(set(graph1.neighbors(attachfrom)))
				pool = list(pool)
				if pool:
					attachto = rd.choice(pool)
					if attachto in [nodenumber_i,rndpick_i]:
						attachtoSelforDad = True
					if graph1.degree(attachto) == 0 or graph1.neighbors(attachto) == [attachto]:
						attachtoOrphan = True
					graph1.add_edge(attachfrom,attachto)
					edgeno += 1
		
		## Check which orphan nodes have been generated or lost
		if graph1.degree(nodenumber_i) == 0 or graph1.neighbors(nodenumber_i) == [nodenumber_i]:
			x += 1
		if (graph1.degree(rndpick_i) == 0 or graph1.neighbors(rndpick_i) == [rndpick_i]):
			if not orphanDad:
				x += 1
		elif orphanDad:
				x -= 1
		if attachtoOrphan and not attachtoSelforDad:
			x -= 1
	return graph1
	   


# Gets the graphlet degree distribution from an orca output file

def get_GDDs_fromfile(filename):
	m = np.loadtxt(filename,dtype=int)
	orb=0
	GDDs = []
	maxs = m.max(0)
	while orb<73:
		d = np.bincount(m[:,orb])[1:]
		ks = np.arange(1,len(d)+1,dtype=float)
		s = np.divide(d,ks)
		t = s/np.sum(s)
		GDDs.append(t)
		orb += 1
	return GDDs

# From a NetworkX graph

def get_GDDs(graph1,transform=False):
	if transform:
		newgraph = nx.Graph()
		i = 0
		transd = {}
		for n in graph1.nodes():
			transd[n] = i
			newgraph.add_node(i)
			i += 1
		for e in graph1.edges():
			newgraph.add_edge(transd[e[0]],transd[e[1]])
		graph1 = newgraph
	edges = []
	for e in graph1.edges():
		if e[0] != e[1]:
			edges.append(e)
	m = len(edges)
	edgelist = np.zeros((m,2),dtype='intc')
	i=0
	while i<m:
		e = edges[i]
		edgelist[i,:] = np.asarray(e,dtype='intc')
		i += 1
	n = graph1.number_of_nodes()
	simOrbits = orca.python_getorbits(edgelist,n,m)
	orb=0
	gdd1 = []
	maxs = simOrbits.max(0)
	while orb<73:
		d = np.bincount(simOrbits[:,orb])[1:]
		ks = np.arange(1,len(d)+1,dtype=float)
		s = np.divide(d,ks)
		t = s/np.sum(s)
		gdd1.append(t)
		orb += 1
	return gdd1

def distance(simGDD,dataGDD):
	i = 0
	distances = []
	while i<73:
		a = simGDD[i]
		b = dataGDD[i]
		if len(a) < len(b):
			c = b.copy()
			c[:len(a)] -= a
		else:
			c = a.copy()
			c[:len(b)] -= b
		c2 = c**2
		d = (np.sum(c2))**(0.5)
		d = d/(2**(0.5))
		distances.append(d)
		i += 1
	d = np.mean(distances)
	return d

# Produces an N-dimensional geometric random graph (default nd = 3 dimensions)

def geometric_model(n,d,nd=3):
    G = nx.Graph()
    i=0
    coords = np.empty((n,nd))
    while i<n:
        pos = np.random.uniform(size=nd)
        G.add_node(i)
        coords[i,:] = pos
        i += 1
    dmat = scipy.spatial.distance.pdist(coords,'euclidean')
    dmatred = scipy.spatial.distance.squareform(dmat)
    a = np.where((dmatred<d) & (dmatred != 0))
    n = len(a[0])
    i = 0
    x = []
    while i<n:
        if a[0][i]<a[1][i]:
            idx = (a[0][i],a[1][i])
            x.append(idx)
        i += 1
    G.add_edges_from(x)
    return G
