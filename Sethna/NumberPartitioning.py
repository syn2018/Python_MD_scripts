"""
A fair split? Number partitioning

Duty: To divide a list of n integers into two sublists with equal sums.
More vividly, divide n players with skills S[i] into two teams of equal total
skill.
"""

import scipy
import pylab
import RandomArray

# Four sample lists to test
S1 = scipy.array([10,13,23,6,20])
S2 = scipy.array([6,4,9,14,12,3,15,15])
S3 = scipy.array([93,58,141,209,179,48,225,228])
S4 = scipy.array([2474,1129,1388,3752,821,2082,201,739])


def AllTeams(n):
    """
    Returns a scipy array of all 2^n vectors of length n composed of +-1,
    giving the possible teams for n total players (so AllTeams(2) should
    return [[-1,-1],[-1,1],[1,-1],[1,1]]). 
    One method is to build an n x 2^n array, with ij entry j/(2**i) mod 2
    (x%2 is x mod 2), and then multiply the array by two and subtract one.
    (If you can do this entirely with scipy arrays, it'll be significantly
    faster: you'll need to use scipy.arange to 2**n, scipy.mod, and 
    scipy.transpose (C. R. Myers, private communication.) The speedup,
    though, won't be important for the sizes in this exercise.)
    """
#    return scipy.array([[(j/2**i)%2 for i in range(n)] 
#    				for j in range(2**n)])*2-1
    return scipy.transpose([scipy.mod(scipy.arange(2**n)/2**i, 2) \
                                for i in range(n)])*2-1


# Build allTeamsArray[n] = AllTeams(n) for n from zero to twelve, so we don't
# need to rebuild the team arrays for each new problem instance.
# (allTeamsArray[0] will look weird.)

allTeamsArray = [AllTeams(n) for n in range(13)]


def ExhaustivePartition(S):
    """
    Returns the minimum size team, given players S. This is nicely done by
    taking scipy.dot of allTeamsArray[len(S)] with S, and then using 
    abs and min.
    """
    return min(abs(scipy.dot(allTeamsArray[len(S)], S)))

def MakeRandomPartitionProblem(N, M):
    """
    Returns a random series of N integers in the range 1 < p < 2**M, guaranteed
    to sum to an even number. Use RandomArray.randint to generate a length N
    vector S of the appropriate range. While sum(S) mod 2 is not zero,
    re-generate S.
    """
    intSize = 2**M
    S = RandomArray.randint(1,intSize+1,N)
    while sum(S)%2 != 0:
	S = RandomArray.randint(1,intSize+1,N)
    return S

def pPerf(N, M, trials=100):
    """
    Returns the fraction of perfect partitions among "trials" runs of 
    N variables of M bits each.
    """
    nPerf = 0
    for trial in range(trials):
	best = ExhaustivePartition(MakeRandomPartitionProblem(N,M))
	if best==0:
	    nPerf+=1
    return float(nPerf)/float(trials)

def Plot_pPerf_versus_kappa(Ns=[3,5,7,9], trials=100):
    """
    Plots pPerf(N,M,trials) versus kappa=M/N for M in the range 1-2*N
    (For N in Ns, sets up empty lists of pPerfs and kappas; for M in
    the range appends float(M)/float(N) to kappas and pPerf(...) to pPerfs;
    then pylab.plot(kappas, pPerfs) and pylab.show().)
    """
    for N in Ns:
        pPerfs = []
        kappas = []
        for M in range(1,2*N):
	    kappas.append(float(M)/float(N))
	    pPerfs.append(pPerf(N,M,trials=trials))
        pylab.plot(kappas, pPerfs)
    pylab.show()

def Plot_pPerf_Scaling():
    """
    Plots the scaling function, for x from -6, 6, in steps of 0.01.
    (You'll need to use scipy.exp, scipy.sqrt, and scipy.pi. Use linewidth=3
    to distinguish the scaling curve.)
    """
    xs =  scipy.arange(-6,6,0.01)
    pPerfTheory = (1-scipy.exp(-scipy.sqrt(3/(2*scipy.pi))*2**(-xs)))
    pylab.plot(xs, pPerfTheory, linewidth=3)
    pylab.show()

def Plot_pPerf_versus_x(Ns=[3,5,7,9], trials=100):
    """
    Plots pPerf(N,M,trials) versus x for M in the range 1-2*N,
    and then calls Plot_pPerf_Scaling to generate the theory curve.
    """
    for N in Ns:
        pPerfs = []
        kappas = []
        for M in range(1,2*N):
	    kappas.append(float(M)/float(N))
	    pPerfs.append(pPerf(N,M,trials=trials))
        xs = N*(scipy.array(kappas)-1)+0.5*scipy.log2(N)
        pylab.plot(xs, pPerfs)
    Plot_pPerf_Scaling()
