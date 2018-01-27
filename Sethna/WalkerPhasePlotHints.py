#
# See the exercise "Walker.pdf" from Walker.html
# in  http://www.physics.cornell.edu/sethna/StatMech/ComputerExercises/
#
import Walker
reload(W)	# for ipython to reload properly after changes in Walker
import scipy, pylab

def WalkerAttractor(gamma, 
		    theta0=0.235, thetaDot0=-0.23, phi0=0.475, phiDot0=-0.026,
		    numPoints=25, numTransient=25):
    """
    Returns the theta-values on the attractor for a walker run with 
    particular slope "gamma"
    Initial condition theta=0.235,thetaDot=-0.23,phi=0.475,phiDot=-0.026
    was empirically found allow the walker to avoid falling down for 
    gamma in an interesting range (from 0.014 to 0.019).
    (1) Create a walker with the given initial conditions
    (2) Set walker.gamma
    (3) Start an empty list of thetas
    (3) Run for numTransient steps to get onto the attractor 
    	(probably you'll want 50 for gamma > 0.018)
	Presume the time between steps is less than 20 
	(so tnew = w.Step(t, t+20.) will typically give a new step)
	If Step returns a step, remember to execute the heel strike
    (4) Run for numPoints further steps, appending theta after each heel strike
    Return the list of thetas
    """
    pass

def MakeBifurcationDiagram():
    """
    Collects points on the attractor in the range 0.014 < gamma < 0.019
    Plots them vs. theta
    #
    Probably the easiest thing to do is to get a big list of gammas and thetas
    (numPoints copies of the same gammas):
    gammaList = []
    thetaList = []
    for gamma in scipy.arange(...):
        thetas = WalkerAttractor(gamma)
	thetaList += thetas
	gammaList += (list of gammas of same length as thetas)
    pylab.plot(gammaList, thetaList, 'k.')            # plots points
    ...
    """
    pass


def demo():
    MakeBifurcationDiagram()

if __name__ == '__main__':
    demo()
    
