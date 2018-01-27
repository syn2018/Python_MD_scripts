import numpy
import pylab		# Plots; also imports array functions cumsum, transpose

def RandomWalk(N=100, d=2):
    """
    Use numpy.cumsum and numpy.random.uniform to generate
    a 2D random walk of length N, each of which has a random DeltaX and
    DeltaY between -1/2 and 1/2.  You'll want to generate an array of 
    shape (N,d), using, for example, random.uniform(min, max, shape).
    """
    pass

def PlotRandomWalkXT(N=100):
    """
    Plot X(t) for one-dimensional random walk 
    """
    pass

def PlotRandomWalkXY(N=100):
    """
    Plot X, Y coordinates of random walk where 
        X = numpy.transpose(walk)[0]
        Y = numpy.transpose(walk)[1]
    To make the X and Y axes the same length, 
    use pylab.figure(figsize=(8,8)) before pylab.plot(X,Y) and
    pylab.axis('equal') afterward.
    """
    pass

def Endpoints(W=10000, N=10, d=2):
    """
    Returns a list of endpoints of W random walks of length N.
    (In one dimension, this should return an array of one-element arrays,
    to be consistent with higher dimensions.)
    One can generate the random walks and then peel off the final positions,
    or one can generate the steps and sum them directly, for example: 
        sum(numpy.random.uniform(-0.5,0.5,(10,100,2))
    """
    pass

def PlotEndpoints(W=10000, N=10, d=2):
    """
    Plot endpoints of random walks.
    Use numpy.transpose to pull out X, Y. 
    To plot black points not joined by lines use pylab.plot(X, Y, 'k.')
    Again, use pylab.figure(figsize=(8,8)) before and
    pylab.axis('equal') afterward.
    """
    pass

def HistogramRandomWalk(W=10000, N=10, d=1, bins=50):
    """
    Compares the histogram of random walks with the normal distribution
    predicted by the central limit theorem.
    #
    (1) Plots a histogram rho(x) of the probability that a random walk
    with N has endpoint X-coordinate at position x. 
    Uses pylab.hist(X, bins=bins, normed=1) to produce the histogram
    #
    (2) Calculates the RMS stepsize sigma for a random walk of length N
    (with each step uniform in [-1/2,1/2]
    Plots rho = (1/(sqrt(2 pi) sigma)) exp(-x**2/(2 sigma**2)) 
    for -3 sigma < x < 3 sigma on the same plot (i.e., before pylab.show).
    Hint: Create x using arange. Squaring, exponentials, and other operations
    can be performed on whole arrays, so typing in the formula for rho will
    work without looping over indices, except sqrt, pi, and exp need to be
    from the appropriate library (pylab, numpy, ...)
    """
    pass


if __name__=="__main__":
   demo()


