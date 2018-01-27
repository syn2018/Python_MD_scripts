# 
# See the exercise "CardiacDynamics.pdf" from CardiacDynamics.html 
# in  http://www.physics.cornell.edu/sethna/StatMech/ComputerExercises/
#
import scipy, scipy.optimize, scipy.integrate
import pylab

# ----------------------------------------------------------------------
# (1.1) Cardiac dynamics
# ----------------------------------------------------------------------

# ----------------------------------------------------------------------
# Part (a)
# ----------------------------------------------------------------------

def V_nullcline(v):
    """
    Returns the value of w(v) giving the V-nullcline of the Fitzhugh-Nagumo
    equations
    """
    return v-(1./3.)*v**3

def W_nullcline(v, gamma, beta):
    """
    Returns the value of w(v) giving the W-nullcline of the Fitzhugh-Nagumo
    equations
    """
    return (1./gamma)*(v + beta)

def V_nullclineArray(v0=-2.0, v1=2.0):
    """Generate a curve in the (v,w) plane representing the V-nullcline.

    Return two arrays:
    the set of v points for the curve, and the set of w points.

    The v points should lie between the specified arguments v0 and v1,
    set by default to generate the nullcline between -2 and 2."""
    v = scipy.arange(v0, v1, 0.01)
    return v, V_nullcline(v)

def W_nullclineArray(gamma, beta, v0=-1.5, v1=0.5):
    """Generate a curve in the (v,w) plane representing the W-nullcline.

    Return two arrays:
    the set of v points for the curve, and the set of w points.

    The v points should lie between the specified arguments v0 and v1,
    set by default to generate the nullcline between -1.5 and 0.5."""
    v = scipy.arange(v0, v1, 0.01)
    return v, W_nullcline(v, gamma, beta)

def PlotNullclines(gamma, beta):
    """Plot the V- and W-nullclines for specified gamma and beta.
    Optionally also find the fixed point (vstar, wstar) for the
    specified gamma and beta, and plot the (vstar, wstar) as a point
    along with the nullclines.

    Using pylab, multiple calls to pylab.plot() can be used to
    plot the various curves in one figure.  To plot only the
    nullclines, you may want to include a call to pylab.show() in
    this function, but you may want to remove it later when plotting
    nullclines along with trajectories.
    """
    V_null = V_nullclineArray()
    W_null = W_nullclineArray(gamma, beta)
    vstar, wstar = FindFixedPoint(gamma, beta)
    pylab.plot(V_null[0], V_null[1], 'b')
    pylab.plot(W_null[0], W_null[1], 'g')
    pylab.plot([vstar], [wstar], 'ro')
    #pylab.show()

def FindFixedPoint(gamma, beta):
    """
    Find and return the fixed point (vstar, wstar) for specified values of
    gamma and beta, by finding the intersection of the V- and W-nullclines.
    
    Use scipy.optimize.brentq to use Brent's method for finding the
    root of a specified function.  In this case the function is the
    difference between the two nullclines:

    def NullclineIntersection(v, gamma, beta):
        return <value of V nullcline>-<value of W nullcline>

    where both <value of V nullcline> and <value of W nullcline> are
    both represented as functions of v.

    Since both nullclines are functions of v, finding a root of the
    function NullclineIntersection() will yield the value of vstar.
    Finding wstar from vstar is trivial.

    vstar = scipy.optimize.brentq(NullClineIntersection, -2, 2.,
               args=(gamma,beta))
    will attempt to find a root to the specified function between the
    points v=-2 and v=2, using gamma and beta as additional arguments
    in the evaluation of the NullclineIntersection() function
    """ 
#   f = lambda v, gamma, beta: (v-(v**3)/3.)-((1./gamma)*(v+beta))
    def f(v, gamma, beta):
        return V_nullcline(v)-W_nullcline(v, gamma, beta)
    vstar = scipy.optimize.brentq(f, -2., 2., args=(gamma, beta))
    wstar = ((1./gamma)*(vstar+beta))
    return vstar, wstar

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

def dydt(y, t, eps, gamma, beta):
    """dydt(y, t, eps, gamma, beta) returns the instantaneous time
    derivative of the FitzHugh-Nagumo v and w variables, as a two-element
    list or scipy array, for use with the scipy odeint function"""
    v,w = y
    return [(1./eps) * (v - (1./3.)*v**3 - w), \
            eps*(v - gamma*w + beta)]
        
def RunFitzNag(y0, t=50., eps=0.2, gamma=0.8, beta=0.7):
    """RunFitzNag(y0, t, eps, gamma, beta) integrates the FitzHugh-Nagumo
    equations starting at state y0, up to time t, with parameters
    eps, gamma, and beta, using the scipy.integrate.odeint function.
    Return the trajectory provided by odeint."""
    tarray = scipy.arange(0.,t, float(t)/500.)
    ytraj = scipy.integrate.odeint(dydt, y0, tarray, args=(eps, gamma, beta))
    return ytraj

def Stimulate(delta, eps=0.2, gamma=0.8, beta=0.7):
    """Stimulate(delta, eps, gamma, beta) stimulates the FitzHugh-Nagumo
    model, starting at the resting state (vstar, wstar), by increasing
    vstar by delta.  The model should subsequently be integrated for
    sufficient time to allow the resting state to be reached.  Stimulate
    should used pylab.plot to generate a plot of the resulting trajectory
    (although not necessarily executing pylab.show() since it is
    useful to plot many stimulations on the same graph).
    You can use scipy.transpose or slicing to extract the list of v's and
    w's from the trajectory returned by RunFitzNag: e.g.
    vs = ytraj[:,0] or vs = scipy.transpose(ytraj)[0]"""
    vstar, wstar = FindFixedPoint(gamma, beta)
    y0 = [vstar+delta, wstar]
    ytraj = RunFitzNag(y0, 50., eps, gamma, beta)
    pylab.plot(ytraj[:,0], ytraj[:,1])
    
def MultiStimulate(eps=0.2, gamma=0.8, beta=0.7):
    """Stimulate along a line delta = [0,0.01,...,1.0], 
    starting at (vstar + delta, wstar) with parameters eps, gamma, beta:
    plot the trajectories along with the nullclines 
    """
    vstar, wstar = FindFixedPoint(gamma, beta)
    for delta in scipy.arange(0.,1.,0.01):
        y0 = [vstar+delta, wstar]
        ytraj = RunFitzNag(y0, 20., eps, gamma, beta)
        pylab.plot(ytraj[:,0], ytraj[:,1])
    V_null = V_nullclineArray()
    W_null = W_nullclineArray(gamma, beta)
    pylab.plot(V_null[0], V_null[1], 'r')
    pylab.plot(W_null[0], W_null[1], 'g')
    pylab.show()

def demo():
    """Demonstrates solution for exercise: example of usage"""
    print "Fitzhugh-Nagumo excitable medium demo"
    MultiStimulate()

if __name__=="__main__":
   demo()


