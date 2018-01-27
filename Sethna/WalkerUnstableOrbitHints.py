"""
When we walk, we are constantly making tiny corrections in order to not fall 
down. We're biologically engineered to have a nearly stable walk: one that
maintains itself with minimal corrections. This happens because we are
at an unstable fixed point. (Try balancing a meter stick vertically on your
palm: small corrections can keep it upright even though it's unstable.
Try "balancing" it at a 45 degree angle...)

Our walker model, even in regimes where it falls down by itself (no periodic
stable attractor), also has an unstable walking mode (unstable periodic orbit)
that is likely biologically significant.
"""
#
# See the exercise "Walker.pdf" from Walker.html
# in  http://www.physics.cornell.edu/sethna/StatMech/ComputerExercises/
#
import Walker
reload(W)	# for ipython to reload properly after changes in Walker
import scipy, scipy.optimize, pylab

def HeelMap(y, w):
    """
    We find this orbit by finding the initial condition (theta, thetaDot) after
    a heel strike that goes to itself after the next heel strike. This is
    done by defining HeelMap(y, w), which for a walker w and a state
    y[n] = (theta, thetaDot) immediately after a heelstrike, returns the 
    difference between y[n] and y[n+1]. It should print out an error message
    if the heel strike doesn't happen by t=20.
    """
    pass

def FindPeriodicOrbit(w, gamma, theta0, thetaDot0):
    """
    We use HeelMap along with a root finder (like scipy.optimize.fsolve)
    to find a nearby periodic orbit, starting at y=[theta0,thetaDot0].
    fsolve will look in the two dimensional space for a zero of HeelMap,
    which thus is a fixed point of the walker after a step, 
    which thus is a periodic orbit for the walker.
    fsolve doesn't care whether it's stable or not...
    The syntax for fsolve is
        scipy.optimize.fsolve(function, y0, args)
    where function(y, args) returns a vector to be zeroed of length len(y),
    y0 is an initial guess for the root, and 
    "args" is a tuple of potential additional arguments for function.
    (For HeelMap, "args" should be (w,). This one-element tuple has a comma
    after the Walker to indicate that it is a tuple and not a parenthetical
    expression.)
    """
    pass

def StanceAngleIntelligentWalker():
    """
    Stance angle theta vs. slope gamma for stable and unstable orbits
    of walker (duplicating figure 3 of Garcia, Chatterjee, Ruina, and Coleman, 
    "The Simplest Walking Model: Stability, Complexity, and Scaling".)
    We discovered through trial and error that the heelstrike condition
    for gamma=0.009, used as an initial condition, put us on one branch
    or the other, depending at whether we started at gamma=0.01 or gamma=0.02.
    We work outward slowly from that initial condition.
    #
    Set up lists of gammas starting from 0.01->0.046, 0.01->0, 0.02->0.046,
    and 0.02->0, with a step size somewhere near 0.003. (Bigger step sizes will
    cause the orbit to get lost.)
    For each list of gammas, 
        make a list of fixed--point thetaStar's
        start the walker in its default state (heelstrike condition for 0.009)
	work outward through the gammas, 
	    use FindPeriodicOrbit to find periodic (theta, thetaDot) for gamma
            set the walker to that pair
	    execute the heel strike
            append theta to thetaStar
        plot thetaStar vs. gammas
    """
    pass

if __name__=='__main__':
    demo()
    
