import pickle
import numpy as np
import random as rd
import networkx as nx
import orca

def perturb(p0,prior,pop):
    # uniform component wise perturbation kernel
    s2 = np.amax(pop,axis=0)
    s1 = np.amin(pop,axis=0)
    
    scale = (s2-s1)/2.0
    
    mini = p0-scale
    ind = np.where(mini<prior[:,0])
    mini[ind] = prior[ind,0]
    
    maxi = p0+scale
    ind = np.where(maxi>prior[:,1])
    maxi[ind] = prior[ind,1]
    
    p = np.random.uniform(size=len(p0),low=mini,high=maxi)
    
    return p
   
def computeWeights(pop,prior,i,weights):
    n = pop[i].shape[0]
    w = np.empty(n)
    if i==0:
		w.fill(1)
    else:
		k = 0
		while k < n:
			denom = 0.0
			
			ty = 1.0/(prior[:,1]-prior[:,0])
			ind = np.where((pop[i][k,:] < prior[:,0]) | (pop[i][k,:] > prior[:,1]))
			ty[ind] = 0
			numer = np.prod(ty)
			
			l = 0
			while l < n:
				tx = 1.0
				s2  = np.amax(pop[i-1],axis=0)
				s1 = np.amin(pop[i-1],axis=0)
				
				scale = (s2-s1)/2.0
				
				tx = np.prod(1.0/(2.0*scale))
				
				ind = np.where(((pop[i-1][l,:]-scale) > pop[i][k,:]) | ((pop[i-1][l,:]+scale) < pop[i][k,:]))
				
				if len(ind[0]) > 0:
					tx = 0.0
				
				denom = denom+tx*weights[i-1][l]
				
				l += 1
			
			w[k] = numer/denom
			
			k += 1
    w = w/np.sum(w)
    
    return w
    
    
    # The graph provided must have integers as node IDs (from 0 to n-1)
    # Otherwise set transform to True
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

def abcsmc_graphs(data,Np,Estart,Eend,Nparam,prior,outputFilename='abcsmcresults.pickle',percentile=25,results=[],transform=False):
    
    if not results:
		results = [0,0,0]
		pop = []
		weights = []
		distances = []
		i = 0
		epsilon = Estart
    else:
		pop = results[0]
		weights = results[1]
		distances = results[2]
		i = len(pop)
		epsilon = np.percentile(distances[i-1],percentile)
    while epsilon > Eend:
		if i > 0:
			epsilon = np.percentile(distances[i-1],percentile)
			if epsilon < Eend:
				epsilon = Eend
		Distance = np.empty(Np)
		print 'Population:',i+1,'Epsilon:',epsilon
		accepted = 0
		pop.append(np.empty((Np,Nparam)))
			
		while accepted < Np:
				
				# if pop == first sample from prior
				# else: sample from pop-1
				
			if i ==0:
				p = np.random.uniform(size=Nparam,low=prior[:,0],high=prior[:,1])
			else:
				ind = np.random.choice(np.arange(Np),1,replace=False,p=weights[i-1])
				p0 = pop[i-1][ind[0],:]
			
				#perturb p0
				p = perturb(p0,prior,pop[i-1])
				
			####################################
			## INCORPORATE GRAPH GENERATOR HERE, a function that takes the parameters, p as an argument and returns the graph, x
			## x must be a networkx Graph that will be compared to the data
			x = my_graph_generator(p) 
			
			# Following if statement in case model returns None for invalid graphs
			if x is not None:
				
				# x must have integers as node IDs (from 0 to len(x)-1)
				# Otherwise provide the argument transform=True to get_GDDs *will take longer
				sim = get_GDDs(x,transform=transform)
				d = distance(sim,data)
				
			else:
				d = epsilon	# Automatic rejection
			#####################################
			
			if d < epsilon:
				pop[i][accepted,:] = p
				Distance[accepted] = d
				accepted += 1
				print 'Hits:', accepted,p,d
		weights.append(computeWeights(pop,prior,i,weights))
		distances.append(Distance)
		results[0] = pop
		results[1] = weights
		results[2] = distances
		print 'Writing results... DO NOT STOP THE SCRIPT NOW!'
		o = open(outputFilename,'w')
		pickle.dump(results,o)	
		o.close()
		i += 1
	
    return results
    
    
    
