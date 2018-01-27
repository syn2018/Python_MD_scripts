"""
Networks is a module that implements a network of nodes and edges,
and various algorithms and measurements implemented on networks.
(Mathematically this is usually called a graph, but we are calling
it a network to avoid confusion with xy graphs, which we will try
to call "plots". Mathematicians use "network" for graphs
with weights on the edges.)
This module contains those parts of the answers to the Small World and
Percolation problems in Sethna's text that are 
independent of what the form of the network is.
"""

import NetGraphics

# -----------------------------------------------------------------------
# This file contains general-purpose class and function definitions for
# networks represented as undirected graphs.
# The initial UndirectedGraph class definition is used both by the 
# Small World and Percolation exercises.
# -----------------------------------------------------------------------

# ***** See also the two test routines following the UndirectedGraph ***** #
# ***** class definition: they should work with the bare-bones       ***** #
# ***** definition provided in the Hints file.                       ***** #

# -----------------------------------------------------------------------
# Defining undirected graph class: used by both Small World and Percolation
# exercises.
# -----------------------------------------------------------------------

class UndirectedGraph:
    """An UndirectedGraph g contains a dictionary (g.neighbor_dict) which
    maps node identifiers (keys) to lists of neighboring nodes (values).
    g.neighbor_dict[node] returns a list [node2, node3, node4] of neighbors.
    Node identifiers can be any non-mutable Python type (e.g., integers,
    tuples, strings, but not lists)."""

    def __init__(self):
        """UndirectedGraph() creates an empty graph g.
	g.neighbor_dict starts as an empty dictionary.  When nodes are
	added, the corresponding values need to be inserted into lists."""
        self.neighbor_dict = {}

    def HasNode(self, node):
        """Does not use the dict.keys() method, which would generate a 
	new list of all nodes each time this is called."""
        # return self.neighbor_dict.has_key(node)
        return node in self.neighbor_dict

    def AddNode(self, node):
	"""Uses HasNode(node) to determine if node has already been added."""
        if not self.HasNode(node):
            self.neighbor_dict[node] = []

    def AddEdge(self, node1, node2):
        """
	Add node1 and node2 to network first
	Adds new edge 
	(appends node2 to neighbor_dict[node1] and vice-versa, since it's
	an undirected graph)
	Do so only if old edge does not already exist 
	(node2 not in neighbor_dict[node1])
	"""
	self.AddNode(node1)
	self.AddNode(node2)
	if node2 not in self.neighbor_dict[node1]:
		self.neighbor_dict[node1].append(node2)
	if node1 not in self.neighbor_dict[node2]:
		self.neighbor_dict[node2].append(node1)
        # alternatively, and more obfuscatedly...
        # self.neighbor_dict.setdefault(node1, []).append(node2)
        # self.neighbor_dict.setdefault(node2, []).append(node1)
        #
        # and if one wanted to prohibit edges from nodes to themselves,
        # and prohibit duplicate edges, then begin with this:
        # if node1 == node2:
        #    return   # do not allow edges from node to itself
        # if node1 in self.neighbor_dict and node2 in self.GetNeighbors(node1):
        #    return   # do not duplicate edges

    def GetNodes(self):
        """g.GetNodes() returns all nodes (keys) in neighbor_dict"""
        return self.neighbor_dict.keys()

    def GetNeighbors(self, node):
        """g.GetNeighbors(node) returns a copy of the list of neighbors of
        the specified node.  A copy is returned (using the [:] operator) so
	that the user does not inadvertently change the neighbor list."""
        return self.neighbor_dict[node]

# Simple test routines

def testgraph():
    """testgraph() creates a simple 5-node UndirectedGraph, and then uses
    the imported NetGraphics module to layout and display the graph.
    The graph is returned from the function.
    
    Notice that we *can* add entries directly to the 'private' member
    variable neighbor_dict without using the AddNode() and/or
    AddEdge() methods, but the drawback of this is that it makes it
    difficult to change the internal implementation later should we
    decide to do so.  Unlike C++ and Java, Python has only a very
    weakly enforced notion of private members, allowing questionable
    programming practices as demonstrated here (if programmers are so
    inclined).  On the other hand, providing a SetNeighbors() method
    for the graph class would be perfectly reasonable, and very useful
    for certain problems."""
    g = UndirectedGraph()
    g.neighbor_dict[0] = [1,2]
    g.neighbor_dict[1] = [0,3]
    g.neighbor_dict[2] = [0,3]
    g.neighbor_dict[3] = [1,2,4]
    g.neighbor_dict[4] = [3]
    NetGraphics.DisplayCircleGraph(g)
    return g

def testplot():
    """testplot() illustrates how to use pylab to create x-y plots,
    including log-log plots.  Each plot must be killed for control to be
    returned to the interpreter."""
    import pylab
    pylab.plot([0,1,2,4,9,16,25])
    # same as pylab.plot([0,1,2,3,4,5], [0,1,2,4,9,16,25])
    pylab.show() # display plot
    # for log-log plot
    x = [0.001, 0.01, 0.1, 1.]
    y = [1.0e-06, 1.0e-04, 1.0e-02, 1.0]
    pylab.loglog(x, y)
    pylab.show()

# ***** After building and testing your network class, build some    ***** #
# ***** Small World or Percolation networks, making use of the       ***** #
# ***** corresponding hints file. Then return to analyze these       ***** #
# ***** networks using the general-purpose routines you write below. ***** #


# ***** Small World exercise routines start here                     ***** #

# -----------------------------------------------------------------------
# Routines for finding path lengths on graphs: used by the Small World
# Network exercise only
# -----------------------------------------------------------------------

def FindPathLengthsFromNode(graph, node):
    """Breadth--first search. See "Six degrees of separation" exercise
    from Sethna's book."""
    dist = 0
    distances = {node: dist}
    currentShell = [node]
    while len(currentShell) > 0:
        nextShell = []
        for k in currentShell:
	    for l in graph.GetNeighbors(k):
            	if l not in distances:
                	nextShell.append(l)
                	distances[l] = dist+1
        dist += 1
        currentShell = nextShell
    return distances

def FindAllPathLengths(graph):
    """
    FindAllPathLengths returns a dictionary, indexed by node pairs,
    storing the shortest path length between those nodes, e.g. for
    small-world networks
    """
    full_dict = {}
    for node in graph.GetNodes():
        distances = FindPathLengthsFromNode(graph, node)
        for node2 in distances.keys():
            full_dict[(min(node, node2), max(node, node2))] = distances[node2]
    return full_dict

def FindAveragePathLength(graph):
    """Averages path length over all pairs of nodes"""
    nodes = graph.GetNodes()
    count = 0
    total = 0.
    for node in nodes:
        distances = FindPathLengthsFromNode(graph, node)
        for node2, dist in distances.items():
            count += 1
            total += dist
    return total/count

# -----------------------------------------------------------------------
# Routine for calculating the clustering coefficient of a graph.
# This was a piece of the original Small World paper, but is not part of the 
# assigned exercise.
# -----------------------------------------------------------------------

def ComputeClusteringCoefficient(graph):
    """Computes clustering coefficient of graph"""
    num_bonds = 0
    num_total_bonds = 0
    C = 0.
    for node in graph.GetNodes():
        neighbors = graph.GetNeighbors(node)
        kv = len(neighbors)
        for index, nei1 in enumerate(neighbors):
            for nei2 in neighbors[index+1:]:
                if nei2 in graph.GetNeighbors(nei1):
                    num_bonds += 1
        Cv = float(num_bonds)/( (kv * (kv-1.) / 2.) )
        C += Cv
    C /= len(graph.GetNodes())
    return C

# -----------------------------------------------------------------------
# Routines for calculating "betweenness", which measures how many shortest
# paths between pairs of nodes on a graph pass through a given node or edge.
# Used in the Small Worlds exercise.
#
# References: (on the Web site)
# Mark E. J. Newman, "Scientific collaboration networks. ii. shortest paths,
# weighted networks, and criticality", Physical Review E 64: 016132, 2002.
# Michelle Girvan and Mark E. J. Newman, "Community structure in social
# and biological networks. Proceedings of the National Academy of Sciences
# 12, 7821-7826, 2002.
# -----------------------------------------------------------------------

def EdgeAndNodeBetweennessFromNode(graph, node):
    """
    Newman's edge and node betweenness inner loop
    Returns partial sum of edge, node betweenness
    """
    dist = 0
    distances = {node: dist}
    currentShell = [node]
    pred = {node:[]}	# predecessor list
    while len(currentShell) > 0:
        nextShell = []
        for k in currentShell:
	    k_pred_l = 0
	    for l in graph.GetNeighbors(k):
            	if l not in distances:
                    nextShell.append(l)
                    distances[l] = dist+1
		if distances[l] == dist+1:
#		     pred.setdefault(l,[]).append(k)
		    if l not in pred: pred[l]=[k]
		    else: pred[l].append(k)
        dist += 1
        currentShell = nextShell
    # Sort nodes by distances
    pairs = [(d,n) for (n,d) in distances.items()]
    pairs.sort()
    pairs.reverse()
    sortedNodes = [n for (d,n) in pairs]
    # Initialize b_k=1 for each of the nodes
    # Use sorted nodes in case disconnected
    nodeBt = dict([(n,1.) for n in sortedNodes])
    # Initialize edgeBt to nothing
    edgeBt = {}
    # Start at most distant node k, split b_k among predecessors and among bonds
    for k in sortedNodes:
        npred = len(pred[k])
	for p in pred[k]:
	   nodeBt[p] += nodeBt[k]/npred
	   edgeBt[(p,k)] = edgeBt.get((p,k),0.) + nodeBt[k]/npred
    return edgeBt, nodeBt

def EdgeAndNodeBetweenness(graph):
    """Returns Newman's edge, node betweenness"""
    edgeBt = {}
    nodeBt = {}
    for node in graph.GetNodes():
        localEdgeBt, localNodeBt = \
                 EdgeAndNodeBetweennessFromNode(graph, node)
        for node, bt in localNodeBt.items():
            nodeBt[node] = nodeBt.get(node, 0) + bt
        for edge, bt in localEdgeBt.items():
            edgeBt[edge] = edgeBt.get(edge, 0) + bt
    return edgeBt, nodeBt

# -----------------------------------------------------------------------
# Sample routines for reading in external files defining networks. Used
# primarily for the later portions of the Small World exercise.
# -----------------------------------------------------------------------

def ReadGraphFromEdgeFile(filename, conversion=None):
    """Reads file with (node1,node2) for each edge in graph"""
    g = UndirectedGraph()
    for line in file(filename):
        node1, node2 = line.split()
        if conversion is not None:
            node1, node2 = conversion(node1), conversion(node2)
        if node1 <= node2:
            g.AddEdge(node1, node2)
    return g

def ReadGraphFromNeighborFile(filename, conversion=None):
    """
    Reads file with [node1,node2,node3] for completely interconnected
    group of nodes (as in actors in a movie)
    Should be read in as a bipartite graph!
    """
    g = UndirectedGraph()
    for line in file(filename):
        sline = line.split()
        if conversion is not None:
            sline = map(conversion, sline)
        allPairs = [(actor1, actor2) for actor1 in sline for actor2 in sline \
                     if actor1 <= actor2]
        for pair in allPairs:
            g.AddEdge(*pair)
    return g


# ***** Percolation exercise routines start here                     ***** #

# -----------------------------------------------------------------------
# Routines for finding clusters in networks. Used in the Percolation exercise.
# -----------------------------------------------------------------------

def FindClusterFromNode(graph, node, visited=None):
    """Breadth--first search
    The dictionary "visited" should be initialized to False for
    all the nodes in the cluster you wish to find
    It's used in two different ways.
    (1) It's passed back to the
        calling program with all the nodes in the current cluster set to
        visited[nodeInCluster]=True, so that the calling program can skip
        nodes in this cluster in searching for other clusters.
    (2) It's used internally in this algorithm to keep track of the
        sites in the cluster that have already been found and incorporated
    See "Building a Percolation Network" in text for algorithm"""
    if visited is None:
        visited = dict([(n, False) for n in graph.GetNodes()])        
    cluster = [node]
    visited[node] = True
    currentShell = graph.GetNeighbors(node)
    while len(currentShell) > 0:
        nextShell = []
        for node in currentShell:
            if not visited[node]:
                nextShell.extend(graph.GetNeighbors(node))
                visited[node] = True
                cluster.append(node)
        currentShell = nextShell
    return cluster

def FindAllClusters(graph):
    """For example, find percolation clusters
    Set up the dictionary "visited" for FindClusterFromNode
    Set up an empty list "clusters"
    Iterate over the nodes;
        if it haven't been visited,
            find the cluster containing it
            append it to the cluster list
        return clusters
    Check your answer using
    NetGraphics.DrawSquareNetworkBonds(g, cl) and
    NetGraphics.DrawSquareNetworkSites(g, cl)
					            
    Optional: You may wish to sort your list of clusters according to their
    lengths, biggest to smallest
    For a list ell, the built-in method ell.sort() will sort the list
    from smallest to biggest;
    ell.sort(cmp) will sort the list according to the comparison function
    cmp(x, y) returns -1 if x < y, returns 0 if x==y, and returns 1 if x>y
    Define ReverseLengthCompare to compare two lists according to the
    unusual definition of inequality, l1<l2 if # len(l1) > len(l2)!
    """
    visited = dict([(n, False) for n in graph.GetNodes()])
    clusters = []
    for node in graph.GetNodes():
        if not visited[node]:
            cluster = FindClusterFromNode(graph, node, visited)
            clusters.append(cluster)
#   clusters.sort(lambda x,y: cmp(len(y), len(x))) # reverse sort of len()
    def ReverseLengthCompare(l1,l2): 
        #   if len(l2)>len(l1): return -1
        #   if len(l2)==len(l1): return 0
        #   if len(l2)<len(l1): return 1
    	return cmp(len(l2),len(l1))
    clusters.sort(ReverseLengthCompare) # reverse sort of len()
    return clusters

def GetSizeDistribution(clusters):
    """Given the clusters, makes up a dictionary giving the number
    of clusters of a given size.
    """
    sizes = [len(cl) for cl in clusters]
    max_size = max(sizes)
    hist = {}
    for cl in clusters: 
    	if len(cl) in hist:
	    hist[len(cl)] += 1
	else:
	    hist[len(cl)] = 1
    return hist
