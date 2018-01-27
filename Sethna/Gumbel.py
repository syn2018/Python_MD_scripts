import pylab
import scipy
import RandomArray

def rhoN(H, mu, beta):
    """
    Returns Gumbel form for extremal value distribution;
    uses scipy.exp
    """
    return scipy.exp(-scipy.exp((mu-H)/beta)+(mu-H)/beta)/beta

def makePeaks(N, M):
    """
    Generates M arrays of N Gaussian random numbers (mean 0, std 1)
    using RandomArray.standard_normal((M,N))
    Returns vector of maxima of each of the M arrays 
      max(ra[j]) will give maximum in that row;
      [max(rv) for rv in ra] will give maxima for each row as vector
    """
    randomArrays = RandomArray.standard_normal((M,N));
    return [max(randomArray) for randomArray in randomArrays]

def plotPeakHistogram(N, M, showPlot=True):
    """
    Make peaks vector using makePeaks
    Calls pylab.hist on peaks 
      with normed=1 (normalize histogram to one), bins=30
    If showPlot, pylab.show() to display histogram
    """
    peaks = makePeaks(N, M)
    pylab.hist(peaks, normed=1, bins=30)
    if showPlot:
    	pylab.show()

def GumbelTheoryGaussian(H):
    """
    Uses B=0.5, delta=2., and Hstar=Hstar[1000]=3.09023 for Gaussian 
    Calculates mu and beta for Gumbel form
    Calls rhoN(H, mu, beta)
    """
    Hstar = 3.09023
    B = 1./2.
    delta = 2.
    mu = Hstar
    beta = 1./(B*delta*Hstar**(delta-1.))
    return rhoN(H, mu, beta)

def plotGumbelDistribution():
    """
    Generates range of Hs spanning histogram (scipy.arange(2.,6.,0.05)
    Generates rhos from Hs using GumbelTheoryGaussian
    Calls pylab.plot on Hs and rhos
      draw black line with 'k-'
      set linewidth=3
    Calls pylab.show()
    """
    Hs = scipy.arange(2., 6., 0.05)
    rhos = GumbelTheoryGaussian(Hs)
    pylab.plot(Hs, rhos, 'k-', linewidth=3)
    pylab.show()


def demo():
    """Demonstrates solution for exercise: example of usage"""
    print "Gumbel Demo"
    # Use smaller number for N, M when debugging!
    plotPeakHistogram(1000,10000,showPlot=False)
    plotGumbelDistribution()

if __name__=="__main__":
   demo()

