"""Basic functionality for iterated maps"""

import scipy
import pylab

# ***** Shared software, 
# ***** used by ChaosLyapunov, InvariantMeasure, FractalDimension, and 
# ***** PeriodDoubling exercises.

def f(x, mu):
    """
    Logistic map f(x) = 4 mu x (1-x), which folds the unit interval (0,1)
    into itself.

    We could treat f a function of one variable x, with a global 
    parameter mu but 
    (1) it would be bad programming practice [global variables are nasty]
    (2) other maps will have other parameter names
    (3) it will make it hard in the Period Doubling problem later to 
        vary mu inside routines.
    So, we recommend making f a function of both variables, and passing 
    the extra argument mu to the various functions that f uses explicitly.

    If you're feeling very fancy, instead of using a definition 
        def f(x, mu): ...
    you may want to try the lambda construction 
        f = lambda x, mu: ...
    It's not a common construction in physics, but it has applications
    in computer science (for example, in self-modifying programs). 

    """
    pass


def Iterate(g, x0, N, args=()):
    """
    Iterate the function g N times, starting at x0, with extra parameters
    passed in as a tuple args. Return g(g(...(g(x))...)). Used to find a 
    point on the attractor starting from some arbitrary point x0.

    Calling Iterate for the Feigenbaum map at mu=0.9 would look like
        Iterate(f, 0.1, 1000, (0.9,))
    Notice that, for a tuple with one element, you NEED the comma after
    the element (otherwise Python thinks the parentheses are for grouping 
    and not denoting a tuple). 

    We'll later be using Iterate to study the sine map 
        fsin(x,B) = B sin(pi x)
    so passing in the function and arguments will be necessary for 
    comparing the logistic map f to fsin.
    
    Inside Iterate you'll want to apply g(x0, *args), where the star in 
    front of args "unpacks" the tuple into its constituents and appends 
    them to the arguments of g.

    You may wish to assume that g depends only on one extra parameter,
    say eta: then you could just pass mu instead of (mu,) and not mess 
    with the arguments. But passing in a whole tuple of extra arguments 
    is a generally useful Pythonic trick.

    """
    pass


def IterateList(g, x0, N, args=()):
    """
    Iterate the function g N-1 times, starting at x0, with extra parameters
    passed in as a tuple args. Returns the entire list 
    (x, g(x), g(g(x)), ... g(g(...(g(x))...))). 

    Can be used to explore the dynamics starting from an arbitrary point 
    x0, or to explore the attractor starting from a point x0 on the 
    attractor (say, initialized using Iterate).

    For example, you can use Iterate to find a point xAttractor on the 
    attractor and IterateList to create a long series of points attractorXs
    (thousands, or even millions long, if you're in the chaotic region), 
    and then use
        pylab.hist(attractorXs, bins=500, normed=1)
        pylab.show()
    to see the density of points.
    """
    pass


def BifurcationDiagram(g, x0, nTransient, nCycle, etaArray, showPlot=True):
    """
    For each parameter value eta in etaArray,
    iterate g nTransient times to find a point on the attractor, and then
    make a list nCycle long to explore the attractor.

    To generate etaArray, it's convenient to use arange: for example,
    BifurcationDiagram(f, 0.1, 500, 128, scipy.arange(0.8, 1.00001, 0.001))
    Notice that you'll want the upper limit to be slightly larger than
    the biggest etaMax you want to plot (but less than etaMax + dEta)
    both to avoid rounding errors and because Python 'range' syntax ends
    with the point before the maximum.

    pylab.plot allows one to plot an entire array of abscissa-values versus an
    array of ordinate-values of the same shape. Our vertical axis (ordinate)
    is a list of arrays of attractor points of length nCycle (created by 
    IterateList after Iterating), one list per value of eta in etaArray.
    Our horizontal axis (abscissa) should thus be a list of arrays 
        [eta, eta, eta, ...] = [eta]*nCycle = eta*scipy.ones(nCycle) 
    of length nCycle.
    Use
        pylab.plot(etaMatrix, xMatrix, 'k,')
        pylab.show()
    to visualize the resulting bifurcation diagram, where 'k,' denotes 
    black points.
    """
    pass
