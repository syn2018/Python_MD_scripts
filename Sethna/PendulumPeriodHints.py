"""Find period of pendulum.
Start at rest: run forward till past zero point, then backtrack.
Multiply time by four to get period."""
#
# See the exercise "Walker.pdf" from Walker.html
# in  http://www.physics.cornell.edu/sethna/StatMech/ComputerExercises/
#

import scipy, scipy.integrate
import pylab

""" 
For convenience, set constants g = 9.8 m/s^2,
L = 1 m (the length of the bar)
"""

def dydt(y,t):
    """
    Define pendulum evolution law ydot = dydt(y,t) for odeint
    The input "y" is defined to be [theta, thetaDot]
    The input "t" is time, which is irrelevant for us 
    (some differential equations are time dependent, but not ours)
    return [thetaDot, thetaDotDot]
    """
    pass

def dzdtheta(z,theta):
    """
    Define pendulum evolution law dz/dtheta = dzdt / dthetadt
    The input z = [theta,thetaDot,t]
    dzdt = [thetaDot, thetaDotDot, 1]
    Return dzdt / thetaDot
    """
    pass

def pendulumPeriod(theta0):
    """
    Find first quarter of period by starting at rest with angle theta0,
    running forward until you pass zero (using dydt), then back up until
    you hit zero (using dzdt).
    (1) Pack initial conditions into vector y = [theta0, thetaDot0]
        Velocity thetaDot0 starts at zero
    (2) While theta > 0, integrate forward in small steps dt = 0.1
        using scipy.integrate.odeint.
        (a) odeint takes an array of time points to output, starting with
	    the original time. We just want time dt, so the array is 
	    times=scipy.array([t,t+dt]). (odeint does not check to convert
	    lists into arrays on input, so we need to do it explicitly.)
	(b) odeint returns a 'trajectory', a vector y(t) for each time 
	in the array. 
            y_traj = scipy.integrate.odeint(dydt, y, times)
	(c) The new y is thus y_traj[1] (or, when using intermediate
	time points, it is, y_traj[-1], the last element in the trajectory 
	array)
    (3) When theta first < 0, 
        (a) Pack current values into z0 = [theta, thetaDot, t]
	(b) Run odeint on dzdtheta with current z0 initial condition, for 
	an array thetas = scipy.array([theta, 0.]), store as z_traj
	(c) Pull out final t from z_traj[-1]
	(d) Return 4.0 times t for zero crossing
    """
    pass

"""
(1) Build an array of initial thetas, perhaps starting at pi/100, going to pi,
in steps of pi/100:  scipy.arange is useful here. (Zero won't work: why?)
You'll need pi = scipy.pi. 
(2) Find the array of periods for each initial theta. (You can use a loop
and append, or build an array using scipy.zeros and assign, or use the
clever Python construction [XXX I've forgotten the name]
    periods = [pendulumPeriod(theta0) for theta0 in theta0s]
(3) use pylab.plot to plot periods vs. theta0s. Don't forget pylab.show()
"""
