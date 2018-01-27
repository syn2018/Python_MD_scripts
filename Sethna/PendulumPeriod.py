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
g = 9.8
L = 1.0 # physical length of bar

def dydt(y,t):
    """
    Define pendulum evolution law ydot = dydt(y,t) for odeint
    The input "y" is defined to be [theta, thetaDot]
    The input "t" is time, which is irrelevant for us 
    (some differential equations are time dependent, but not ours)
    return [thetaDot, thetaDotDot]
    """
    theta,thetaDot = y
    thetaDotDot = -(g/L) * scipy.sin(theta)
    return [thetaDot,thetaDotDot]

def dzdtheta(z,theta):
    """
    Define pendulum evolution law dz/dtheta = dzdt / dthetadt
    The input z = [theta,thetaDot,t]
    dzdt = [thetaDot, thetaDotDot, 1]
    Return dzdt / thetaDot
    """
    theta,thetaDot,t = z
    thetaDotDot = -(g/L) * scipy.sin(theta)
    return [1.,thetaDotDot/thetaDot,1./thetaDot]

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
    # Pack evolving variables into vector of initial conditions
    # Start at rest. Run to where theta crosses zero: 1/4 of a period
    thetaDot0 = 0.0
    y = [theta0, thetaDot0]
    t = 0.0
    dt = 0.1 # Must be much smaller than period!
    while y[0] > 0:
        times = scipy.array([t, t+dt])
        y_traj = scipy.integrate.odeint(dydt, y, times)
	y = y_traj[-1]
	t += dt
    # Run backwards to zero crossing
    z0 = [y[0],y[1],t]
    thetas = scipy.array([y[0],0.])
    z_traj = scipy.integrate.odeint(dzdtheta, z0, thetas)
    theta, thetadot, t = z_traj[-1]
    # Period is four times time to first zero crossing
    period = 4.0*t
    return period

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

def demo():
    pi = scipy.pi
    theta0s = scipy.arange(pi/100.,pi,pi/100.)
    periods = [pendulumPeriod(theta0) for theta0 in theta0s]
    pylab.plot(theta0s, periods)
    pylab.show()

if __name__ == '__main__':
    demo()
