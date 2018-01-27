"""
StochasticCells runs a simple homodimerization reaction for a small total number
of molecules N, comparing the approximate continuum dynamics with the true
stochastic dynamics.
"""
#
# See the exercise "StochasticCells.pdf" from StochasticCells.html
# in  http://www.physics.cornell.edu/sethna/StatMech/ComputerExercises/
#
import scipy, scipy.integrate
RA = scipy.random
import pylab

"""Global values of rate constants (probably should be members of a class)"""
kb = 1.0
ku = 2.0

def dydt(y,t):
    """
    Gives the time evolution law dydt for y = [M,D], in the form needed by 
    odeint
    """
    pass

def PlotODE(N, tMax=1.0, dt = 0.01):
    """Plots the continuum time evolution of M(t), given N total 
    monomer molecules and no dimer molecules at time t=0. 
    Uses scipy.arange to produce an array of times; calls
    scipy.integrate.odeint(dydt, y0, times) to get yTrajectory
    M(t) is first column of yTrajectory = yTrajectory[:,0]
    uses pylab.plot to plot times versus M(t)
    """
    pass

def StochasticTrajectory(N, tMax=10.0):
    """
    Implements the Gillespie algorithm, as described in the text. If 
        t1, t2, t3, ..., tFinal 
    are the resulting reaction times and 
        M1, M2, M3, ..., MFinal 
    are the number of monomers just after each reaction, the routine returns 
    an array 
        times = [0.0, t1, t1, t2, t2, ..., tFinal, tMax]
    and an array 
        Ms = [N, N, M1, M1, M2, M2, ..., MFinal, MFinal]
    (suitable for plotting, since the concentration M(t) = M_n between t_n and
    t_{n+1}, and then abruptly jumps to M_{n+1}). This is easy to do:
    initialize them at t=0, append just before and just after each
    reaction, and add a point at t=tMax.

    You can generate tWait in two ways:
    	(1) import random, generate a random number r0 uniformly in [0,1),
	and let tWait = -log(1-r0)/gammaTot. (Check that this generates
	an exponential distribution! We take 1-r0 to avoid rare
	errors when r0 happens to be exactly zero.)
	(2) for convenience, define RA=scipy.random, and generate a random 
	number with an
	exponential distribution of mean 1/gammaTot. RA insists
	on returning an array of random numbers, so you need to tell it
	to give an array of length one and pull it out:
	    tWait = RA.exponential(1.0/gammaTot, 1)[0]
        Of course, you could also generate an exponential random number
        of width one and then divide by gammaTot:
	    tWait = RA.exponential(1.0, 1)[0]/gammaTot
    Generating r is easy: call RA.random() or random.random()
    and multiply by gammaTot.
    Notice that, since there are only two reactions, you can just check if
    r > bindingRate to see if you want to bind or unbind.
    """
    pass

def PlotODEVsStochastic(N, tMax=1.0, dt = 0.01):
    """Plots the continuum time evolution of M(t) versus
    the stochastic trajectory given by the Gillespie algorithm.
    Again, N total monomer molecules and no dimer molecules at time t=0. 
    """
    pass

def PlotODEVsAverageStochastic(N, nAverage, tMax=1.0, dt = 0.001):
    """Plots the continuum time evolution of M(t) versus
    the average of nAverage stochastic trajectories. 
    Computes the stochastic averages at the same 
        times = (dt, 2 dt, ...)
    that odeint uses (except for zero).
    The stochastic simulation returns values at irregular times: how
    can we evaluate M at regular intervals? We can find the stochastic
    interval 
      [stochasticTimes[indexStochastic], stochasticTimes[indexStochastic+1])
    containing time[index] by searching forward from the previous interval,
    something like this:
        ...
        indexStochastic = 0
	for index in range(len(times)):
	    while (indexStochastic < len(stochasticTimes)-1) \
	    	   & (stochasticTimes[indexStochastic+1] < times[index]):
	        indexStochastic+=1
	    (add stochastic M[indexStochastic] to total M[index]...)
    Or, we can use the "searchsorted" method defined for arrays:
        t = scipy.searchsorted(stochasticTimes, times)
        m = [stochasticMs[e-1] for e in t]
	(add m to total M)
    Check that your method is working by plotting the stochastic
    trajectory (as lines) and the interpolated trajectory (as points: 'ro')
    """
    pass


if __name__=="__main__":
   demo()
