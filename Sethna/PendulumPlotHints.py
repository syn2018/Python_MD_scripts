# 
# See the exercise "Pendulum.pdf" from Pendulum.html 
# in  http://www.physics.cornell.edu/sethna/StatMech/ComputerExercises/
#
"""Hints for plotting theta(t) for the pendulum"""

# We introduce the pylab plotting package
# We illustrate it using the equation yDot = y solved using the Euler 
# method.

import scipy, scipy.integrate
import numpy
import pylab

# You'll want to include the physical properties and initial conditions
# for the pendulum here

# Successively smaller time steps

dts = [0.1,0.01,0.001]
tMax = 3.

for dt in dts:
    # You'll want to initialize theta and thetaDot here
    # and set up empty arrays for theta_of_t
    t = 0.0
    yOft = 1
    times = []
    ys = []
    eps = 1.0e-10		# To run up to and including tMax
    while t < tMax*(1.+eps): 	# add eps to avoid rounding errors
        times.append(t)
        ys.append(yOft)
        # Solve for theta(t+dt) in terms of theta(t) as before
	t += dt
        yOft += yOft * dt
    # Add curve to plot
    pylab.plot(times, ys)
pylab.show()

# Zoom in at t=1.0, where y(t) = e^t = e = 2.71828... 
# Notice how the accuracy improves as dt gets smaller, but only linearly in t

# Example for odeint:

def dydt(y,t):
    # Note: y is a vector of length one, but dydt is also a vector
    return y

y0 = [1.0]
times = numpy.arange(0.0,3.0,0.1)
trajectory = scipy.integrate.odeint(dydt, y0, times)
ys = trajectory[:,0]
pylab.plot(times, ys)
pylab.show()

