import random, numpy, numpy.random, os, time
import NetGraphics
import MultiPlot, IntegerHistogram
pylab = MultiPlot.pylab

# ***** Start by reading the exercise "PercolationComputation.pdf"   ***** #
# ***** (and later "PercolationScaling.pdf") from Percolation.html   ***** #
# ***** in www.physics.cornell.edu/sethna/StatMech/ComputerExercises/***** #

# ***** Then define the general-purpose UndirectedGraph class        ***** #
# ***** using NetworksHints.py (renamed Networks.py), or import your ***** #
# ***** answers previously written for SmallWorldNetworks.           ***** #

# ***** Then return here to build some percolation networks          ***** #

# Import your network definitions
import Networks
reload(Networks) # for ipython %run to reload properly

# -----------------------------------------------------------------------
# Bond percolation on a square lattice
# -----------------------------------------------------------------------

def MakeSquareBondPercolation(L, p):
    """Instantiate empty graph g = Networks.UndirectedGraph()
    Add nodes on square grid using g.AddNode
    Add horizontal and vertical bonds to neighbors with probability
    p using g.AddEdge
    (1) random.random() will generate random float 0 <= r < 1,
        so "random.random() < p" is true with probability p
    (2) (i+1)%L will give periodic boundary conditions
    return g
    Check your answer using
    NetGraphics.DrawSquareNetworkBonds(g)
    """
    pass

# ***** After creating, displaying, and debugging your bond          ***** #
# ***** percolation networks, write the routines for finding         ***** #
# ***** clusters in Networks.py (put there because they are general) ***** #
# ***** and use the NetGraphics routines  below to plot the clusters ***** #
# ***** in different colors to aid debugging.                        ***** #

def PlotBondPercolationBonds(L=10, p=0.5, seed=1):
    """
    Uses DrawSquareNetworkBonds in NetGraphics to graph the percolation
    network made by MakeSquareBondPercolation and the clusters returned
    by Networks.FindAllClusters. Best for small networks to debug.
    """
    pass

def PlotBondPercolation(L=10, p=0.5, seed=1):
    """
    Uses DrawSquareNetworkSites in NetGraphics to graph the percolation
    network made by MakeSquareBondPercolation and the clusters returned
    by Networks.FindAllClusters. Best for large networks to explore
    universality.
    """
    pass

# -----------------------------------------------------------------------
# Site percolation on a triangular lattice
# -----------------------------------------------------------------------

def MakeTriangularSitePercolation(L, p):
    """Triangular lattice is implemented by bonds to neighbors separated
    by [0,1], [1,0], [-1, 1] and their negatives, so we need an edge
    connecting [i,j] to [i,j+1], [i+1,j], [i-1,j+1], for each point
    in the lattice, modulo the lattice size L.
    Either
    (a) add one node at a time, and fill in all the neighbors, or
    (b) use numpy.random.random((L,L)) to
        fill a whole matrix at once to determine which sites are occupied;
      (i) add nodes and edges as appropriate, or
      (ii) dispense with the dictionary and write GetNodes() and
      GetNeighbors functions directly from the array
    Check your answer using
    NetGraphics.DrawTriangularNetworkSites(g, cl)
    (For small L, the graphics may cut off your graph early: use
    NetGraphics.DrawTriangularNetworkSites(g, cl, L) to fix this)
    """
    pass


def PlotSitePercolation(L=10, p=0.5, seed=1, scale=0):
    """
    Uses DrawTriangularNetworkSites to draw clusters.
    """
    pass

def PlotSitePercolationBiggest(L=10, p=0.5, seed=1, scale=0):
    random.seed(seed)
    """
    Uses DrawTriangularNetworkSites to draw only the largest cluster,
    by setting cl to the result of Networks.FindAllClusters (presuming
    it sorts by size) and by passing in [cl[0]] as the cluster list.
    """
    pass

# ***** More advanced exercise: knowledge of scaling and             ***** #
# ***** renormalization group useful 	                             ***** #
# ***** Read "PercolationScaling.pdf") from Percolation.html.        ***** #

def PlotLogLogSizeDistributionHist(L=100, p=0.5, log10binsize=0.25,
				   min=1.e-10, repeats=1):
    """Make bond percolation graph
    Find clusters
    Make list of sizes of all the clusters
    Make two dictionaries, S and D for the different size and probability
    curves. (The dictionary keys are used by MultiPlot as curve labels).
    Make two lists, S['bond'] = [0,1,2,...] and D['bond']=[D(0), D(1), ...],
    where D(S) is the number of clusters of size S
    pylab.plot(S['bond'],D['bond']), then pylab.show()
    should generate a plot of D(S)
    If you ensure that D(S) > 0 for all points (say, by adding min=1.e-10 to
    all the entries), you can do
        pylab.loglog(S['bond'],D['bond'])
    to do a log-log plot
    #
    Make two more lists, S['site'] and D['site'] for site percolation
    #
    You can now either change colors by hand
    ("bo" for blue circles with lines, "ro" for red, ...),
        pylab.plot(S['bond'],D['bond'], "bo-")
        pylab.plot(S['site'],D['site'], "ro-")
        pylab.show()
    or you can use our package
    MultiPlot.MultiPlot(S, D, xlabel='S', ylabel='D(S)')
    which should give legends and axis labels.
    Do help(MultiPlot.MultiPlot) for more details, for log-log plots, etc.
    #
    Add the theory curve using the cool property that one
    can take powers of numpy arrays:
    D['power'] = S['power']**(-tau)!
    #
    You'll find that the data curves become unreliable as soon as there
    are sizes with only a few clusters. We can cure this by binning
    several sizes in one bin. The bookkeeping for doing this is a bit
    messy, so we provide a package, IntegerHistogram.
    (1) Take your data points (sizes, all integers)
        make a long list of them called "sizes"
    (2) Set up the bins you like. I recommend
        bins =10.**arange(0.0,log10(L*L),0.25)
    (3) S['bond'], D['bond'] = IntegerHistogram.IntegerHistogram(sizes, bins)
        will give the centers and average counts per integer for each bin"""
    pass

def PlotLogLogBeta(L=10, pc=0.5, n=10, repeats=1):
    pass

def PlotPofp(Ls = [5, 10, 20], repeats=1):
    """Here, because you'll want to study small systems as well as large
    ones, you'll need to average over several realizations of the
    disorder.
    (1) Remember that the variance of the mean is 1/(repeats-1) times
    the variance of the sample, so the standard deviation
    is sqrt(<(P-Pbar)^2>/(repeats-1)) =
        sqrt((P2Sum - PSum**2/repeats)/(repeats*(repeats-1))
    (2) A good range of p might be 2.5/L on each side of pc?
    Generate the array of p values using numpy.arange.
    (3) To use MultiPlot.
     (a) Make dictionaries pValues[L, type], Pbar[L, type], and sigma[L, type]
     where type = 'site' or 'bond'
     (b) Plotting the (unrescaled) curves can be done with
     MultiPlot.MultiPlot(pValues, Pbar,
    			 xform='p->p', yform='P->P',
     			 yerrdata=sigma, yerrform = 'sigma->sigma',
     			 loc='upper left',
     			 showIt=True)
    (c) Doing your scaling collapses, you'll need to add
    keyNames=('L',"type") and a scalingParams dictionary to Multiplot's
    argument list."""
    pass




# -----------------------------------------------------------------------
# For percolation on a protein interaction network
# -----------------------------------------------------------------------

def RandomPrune(graph, delta):
    pass

def Percolate(graph, deltas):
    pass

def MultiPercolate(graph, deltas, numruns):
    pass

def TestMultiPercolate():
    pass
