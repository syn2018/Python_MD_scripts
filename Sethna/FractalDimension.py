# 
# See the exercise "FractalDimensions.pdf" from FractalDimensions.html 
# in  http://www.physics.cornell.edu/sethna/StatMech/ComputerExercises/
#
"""Fractal Dimension exercise"""

from IterateLogistic import *

def GetPn(mu, epsilonList, nSampleMax, nTransient=10000):
    """
    Generates probability arrays P_n[epsilon]. 
    Specifically,
     finds a point on the attractor by iterating nTransient times
     collects points on the attractor of size nSample
     for each epsilon in epsilonList, 
      generates bins of size epsilon for the range (0,1) of the function
          bins = scipy.arange(0.0,1.0+eps,eps)
      finds the number of points from the sample in each bin, using
      the histogram function
          numbers, bins = pylab.mlab.hist(sample, bins=bins)
      and computes the probability P_n[epsilon] of being in each bin.
    In the period doubling region the sample should of size 2^n so that 
    it covers the attractor evenly.
    """
    if mu < 0.85:
       nSample = min(nSampleMax, 4)
    else:
       nSample = nSampleMax
    xAttractor = Iterate(f, 0.1, nTransient, (mu,))
    sample = IterateList(f, xAttractor, nSample, (mu,))
    P_n = {}
    for eps in epsilonList:
        bins = scipy.arange(0.0,1.0+eps,eps)
        numbers, bins = pylab.mlab.hist(sample, bins=bins)
        P_n[eps] = numbers / float(nSample)	# Probability       
    return P_n


def DimensionEstimates(mu, epsilonList, nSampleMax):
    """
    Estimates the capacity dimension and the information dimension 
    for a sample of points on the line.
    The capacity dimension is defined as
       D_capacity = lim {eps->0} (- log(Nboxes) / log(eps))
    but converges faster as
       D_capacity = - (log(Nboxes[i+1])-log(Nboxes[i])) 
       			/ (log(eps[i+1])-log(eps[i]))
    where Nboxes is the number of intervals of size eps needed to 
    cover the space. The information dimension is defined as
       D_inf = lim {eps->0} sum(P_n log P_n) / log(eps)
    but converges faster as
       S0 = sum(P_n log P_n)
       D_inf = - (S0[i+1]-S0[i]) 
       			/ (log(eps[i+1])-log(eps[i]))
    where P_n is the fraction of the list 'sample' that is in bin n,
    and the bins are of size epsilon. You'll need to add a small 
    increment delta to P_n before taking the log: delta = 1.e-100 will 
    not change any of the non-zero elements, and P_n log (P_n + delta)
    will be zero if P_n is zero.

    Returns three lists, with epsilonBar (geometric mean of neighboring
    epsilonList values), and D_inf, and D_capacity values for each
    epsilonBar
    """
    D_inf = []
    D_capacity = []
    epsilonBar = []
    delta = 1.e-100	# Add to make log finite

    P_n = GetPn(mu, epsilonList, nSampleMax)

    Nboxes = [] 	# Number of non-zero P_n
    S0 = []		# Zero-dimensional entropy -sum(P_n log(P_n))

    for eps in epsilonList:
	Nboxes.append(sum( P_n[eps] > 0 ))
	S0.append(-sum(P_n[eps] * scipy.log(P_n[eps]+delta)))

    epsBar = []
    D_capacity = []
    D_inf = []
    for i in range(len(epsilonList)-1):
        epsi = epsilonList[i]
        epsiP1 = epsilonList[i+1]
        epsBar.append(scipy.sqrt(epsiP1*epsi))
        D_capacity_estimate = -(scipy.log(Nboxes[i+1])-scipy.log(Nboxes[i])) \
				/ (scipy.log(epsiP1)-scipy.log(epsi))
	D_capacity.append(D_capacity_estimate)
        D_inf_estimate = -(S0[i+1]-S0[i]) \
				/ (scipy.log(epsiP1)-scipy.log(epsi))
	D_inf.append(D_inf_estimate)

    return epsBar, D_capacity, D_inf

def PlotDimensionEstimates(mu, nSampleMax = 2**18, \
                           epsilonList=2.0**scipy.arange(-4,-16,-1)):
    """
    Plots capacity and information dimension estimates versus 
    epsilon, to allow one to visually extrapolate to zero. 
    Uses pylab.semilogx.

    We found using from 16 to one million bins useful. In the chaotic 
    region, the sample should be larger than the number of bins.

    Try mu = 0.9 in the chaotic region; compare with
    mu=0.9, nSampleMax = 10000 to see what the 'finite sample' effects of
    having fewer points than bins looks like.

    Try mu = 0.8 in the periodic limit cycle region. Are the dimensions
    of the attractor different?

    Try muInfinity = 0.892486418. Does the graph extrapolate to 
    near D_inf = 0.517098 and D_capacity = 0.538, as measured in the
    literature? 

    Try muInfinity +- 0.001, to see how the dimension of the attractor
    looks fractal on long length scales (big epsilon), but becomes 
    homogeneous on short length scales. This is the reverse of what 
    happens in percolation. It's because short length scales correspond 
    to long time scales...
    """
    epsilon, D_capacity, D_inf = \
        DimensionEstimates(mu, epsilonList, nSampleMax)
    pylab.plot(epsilon, D_capacity, 'bo', label=str("capacity %s"%mu))
    pylab.plot(epsilon, D_inf, 'ro', label=str("info %s"%mu))
    # Information dimension from theory at muInfinity
    capacity_theory = 0.*scipy.array(epsilon) + 0.538
    information_theory = 0.*scipy.array(epsilon) + 0.5170976
    pylab.plot(epsilon, capacity_theory, 'b-', label='capacity theory')
    pylab.plot(epsilon, information_theory, 'r-', label='information theory')
    pylab.semilogx()
    pylab.legend()
    pylab.show()

    return epsilon, D_capacity, D_inf


def yesno():
    response = raw_input('    Continue? (y/n) ')
    if len(response)==0:        # [CR] returns true
        return True
    elif response[0] == 'n' or response[0] == 'N':
        return False
    else:                       # Default
        return True
    
def demo():
    """Demonstrates solution for exercise: example of usage"""
    nSampleMax = 2**16
    epsilonList=2.0**scipy.arange(-4,-16,-1)
    muInfinity = 0.892486418
    print "Fractal Dimensions Demo"
    print "  Attractor becomes fractal at muInfinity = 0.892486..."
    BifurcationDiagram(f, 0.1, 500, 128, scipy.arange(0.89, 0.894, 0.00003))
    if not yesno(): return
    print "  Fractal Dimension Estimates at mu=0.9 (should be one)"
    PlotDimensionEstimates(0.9, nSampleMax = nSampleMax, \
                           epsilonList=epsilonList)
    if not yesno(): return
    print "  Fractal Dimension Estimates at mu=0.89 (should be zero)"
    PlotDimensionEstimates(0.89, nSampleMax = nSampleMax, \
                           epsilonList=epsilonList)
    if not yesno(): return
    print "  Fractal Dimension Estimates at mu=muInfinity (should be theory)"
    PlotDimensionEstimates(muInfinity, nSampleMax = nSampleMax, \
                           epsilonList=epsilonList)

if __name__=="__main__":
   demo()

