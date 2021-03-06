# 
# See the exercise "InvariantMeasure.pdf" from InvariantMeasure.html 
# in  http://www.physics.cornell.edu/sethna/StatMech/ComputerExercises/
#
"""Invariant Measure exercise"""

from IterateLogistic import *

def PlotInvariantDensityWithBoundaries(g, x0, num_boundaries, 
                                       args=(), xMax=0.5):
    """Plots the invariant density, together with the first num_boundaries
    iterates of xMax = 0.5 (which should coincide with folds, and hence
    cusps, in the invariant density). Plots the iterates f^[n](xMax) as red
    circles 'ro' at rho = n."""
    xAttractor = Iterate(g, x0, 1000, args)
    n, bins, patches = \
       pylab.hist(IterateList(g, xAttractor, 1000000, args), \
                    bins = 2000, normed = 1)
    boundaries = IterateList(g,xMax, num_boundaries, args)[1:]
    pylab.plot(boundaries, 1.+scipy.arange(len(boundaries)), 'ro-',
                 linewidth=1, antialiased=True)
    pylab.show()

# Plot bifurcation diagram; explain boundaries as images of x=1/2
# BifurcationDiagram(f, 0.1, 500, 128, scipy.arange(0.8, 1.0001, 0.001),
#                    showPlot=False)
# PlotBoundaries(f, 8, scipy.arange(0.8, 1.0001, 0.001)) 

def PlotBoundaries(g, nImages, etaArray, xMax=0.5):
    """
    For each parameter value eta in etaArray,
    iterate the point xMax nImages times, and plot the result
    (not including xMax) versus eta. We recommend using
       pylab.plot(etas, halfImages, 'ro')
    where the 'ro' will draw red circles.

    Usually xMax will be the peak in the function g (as hinted at by 
    its name).
    
    This can be used in conjunction with BifurcationDiagram to explain
    the boundary structure in the chaotic region. If you remove 
    pylab.show() from BifurcationDiagram, this plot will be 
    superimposed on the other.
    """
    for n in range(nImages):
        halfImages = [Iterate(g, xMax, n+1, (eta,)) for eta in etaArray]
        pylab.plot(etaArray, halfImages, 'r-', linewidth=1, antialiased=True)
    pylab.axis([etaArray[0],etaArray[-1],0,1])
    pylab.show()

def demo():
    """Demonstrates solution for exercise: example of usage"""
    print "Invariant Measure Demo"
    print "  Creating Invariant Density at mu=0.9"
    PlotInvariantDensityWithBoundaries(f,0.1,20,(0.9,))
    print "  Creating Bifurcation Diagram With Boundaries"
    BifurcationDiagram(f, 0.1, 500, 500, scipy.arange(0.8, 1.00001, 0.001),
                           showPlot=False)
    PlotBoundaries(f, 6, scipy.arange(0.8, 1.00001, 0.002))

if __name__=="__main__":
   demo()

