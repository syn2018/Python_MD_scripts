"""SmallWorldNetwork contains the parts of a solution to the 
small-world network problem in Sethna's book that relate specifically
to the Watts-Strogatz small world networks. The more general graph
algorithms for path lengths and betweenness are in Networks.py."""

# ***** Start by reading the exercise "SixDegreesOfSeparation.pdf"   ***** #
# ***** from SmallWorld.html in                                      ***** #
# ***** www.physics.cornell.edu/sethna/StatMech/ComputerExercises/   ***** #

# ***** Then define the general-purpose UndirectedGraph class        ***** #
# ***** using NetworksHints.py (renamed Networks.py), or import your ***** #
# ***** answers previously written for Percolation.                  ***** #

# ***** Then return here to build some small-world networks          ***** #

import random, os
import numpy, pylab
import NetGraphics, MultiPlot

# Import your network definitions
import Networks 
reload(Networks)	# Helps with ipython %run command

# Small world and ring graphs

def MakeRingGraph(num_nodes, Z):
    """
    Makes a ring graph with Z neighboring edges per node.
    """
    pass

def AddRandomEdges(graph, num_edges_tried):
    """Attempts to add num_edges_tried random bonds to a graph. It may add 
    fewer, if the bonds already exist."""
    pass

def MakeSmallWorldNetwork(L, Z, p):
    """
    Makes a Watts and Strogatz small--world network of size L and Z 
    neighbors, except that there are p*Z*L/2 added shortcuts rather
    than rewirings."""
    pass

def SmallWorldSimple(L, Z, p, dotsize=4):
    """
    Generate and display small world network. Creates a graph g using
    MakeSmallWorldNetwork, and uses the NetGraphics command 
    DisplayCircleGraph, with only the mandatory argument g. Returns g.
    """
    pass

# ***** After creating, displaying, and debugging your small world   ***** #
# ***** graphs, go to Networks.py and develop and debug the routines ***** #
# ***** for finding path lengths in graphs. (They are put in         ***** #
# ***** Networks.py because they in principle could be used for any  ***** #
# ***** graph. Then return here to study the scaling properties of   ***** #
# ***** path lengths in small-world networks.                        ***** #

def MakePathLengthHistograms(L=100, Z=4, p=0.1):
    """
    Plots path length histograms for small world networks.
    Find list of all lengths
    Use pylab.hist(lengths, bins=range(max(lengths)), normed=True) """
    pass

def FindAverageAveragePathLength(L, Z, p, numTries):
    """Finds mean and standard deviation for path length between nodes,
    for a small world network of L nodes, Z bonds to neighbors, 
    p*Z*L/2 shortcuts, averaging over numTries samples"""
    pass

def GetPathLength_vs_p(L, Z, numTries, parray):
    """Calculates array of mean pathlengths and sigmas for small
    world networks; returns pathlengths and sigmas"""
    pass

def PlotPathLength_vs_p(L, Z, numTries=2,
                        parray=10.**numpy.arange(-3., 0.001, 0.25)):
    """Plots path length versus p"""
    pass

def PlotScaledPathLength_vs_pZL(LZarray, numtries=2, 
                                pZLarray=10.**numpy.arange(-1., 2.001, 0.25)):
    """
    PlotScaledPathLength_vs_pZL(((L1,Z1),(L2,Z2),...), numtries,
    				   [pZLa,pZLb,pZLc...])
    will plot the scaled path length for small world networks of size Li and
    neighbors Zi, at scaled rewiring probabilities pZLa, pZLb, ...
    Uses either MultiPlot.py to do the scaling, or rescales by hand, depending
    on the implementation chosen.
    To rescale, p is multiplied by Z*L and the mean path length ell is
    multiplied by 2*Z/L.
    """
    pass

# ***** Clustering coefficient was calculated in the original small  ***** #
# ***** world paper, but is not assigned (optional) in this exercise.***** #

def FindAverageClusteringCoefficient(L, Z, p, numTries):
    """Finds clustering coefficient for small world graph"""
    pass

def GetClustering_vs_p(L, Z, numTries, parray):
    pass

def PlotClustering_vs_p(L, Z, numTries,
                        parray=10.**numpy.arange(-3., 0.001, 0.1)):
    pass

# ***** Again, go to Networks.py to generate and debug your          ***** #
# ***** algorithms for measuring Betweenness. (The algorithms are    ***** #
# ***** described not in the exercise writeup, but in the original   ***** #
# ***** papers by Mark Newman and Michelle Girvan.                   ***** #
# ***** Use them on small-world networks here.   ***** #

def TestBetweennessSimple():
    """
    Makes a simple graph for which one can calculate the betweenness by 
    hand, to check your algorithm.
    """
    g = Networks.UndirectedGraph()
    g.AddEdge(0,1)
    g.AddEdge(0,2)
    g.AddEdge(1,3)
    g.AddEdge(2,3)
    g.AddEdge(3,4)
    edgeBt, nodeBt = Networks.EdgeAndNodeBetweenness(g)
    return g, edgeBt, nodeBt

def SmallWorldBetweenness(L, Z, p, dotscale=4, linescale=2, windowMargin=0.02):
    """
    Display small world network with edge and node betweenness,
    using NetGraphics routine DisplayCircleGraph, passing in arguments
    for edge-weights and node_weights. Passes through the arguments for 
    dotscale, linescale, and windowMargin, to fine-tune the graph
    """
    pass
