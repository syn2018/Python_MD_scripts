import scipy
import RandomArray
import IsingStaticGraphics
import DynamicLattice

class IsingModel:
    """Ising model class"""

    def __init__(self, N=10, T=2./scipy.log(1.+scipy.sqrt(2.)), 
                       H=0., seed=(1,1)):
        """
	Call RandomArray.seed with arguments seed[0], seed[1];
	Set self.N to be N
	Build self.lattice to be random 50/50 0/1, using random_integers
	Call self.SetTemperatureField with T, H
	"""
	pass

    def SetTemperatureField(self, T, H):
        """
	Sets self.privateTemperature, self.privateField
	#
	#
	After debugging, for faster performance, set up algorithm-dependent
	data structures to avoid taking exponentials for every flipped spin
	#
	Heat bath algorithm:
	Sets up self.heatBathProbUp to be array of doubles (scipy.Float) 
	of length Z+1=5 (the number of different neighbors up (nUp) possible)
	Figure out spin energies eUp, eDown given nUp; figure out Boltzmann
	relative probabilities boltzUp, boltzDown; 
	If T!= 0, set heatBathProbUp(nUp) to boltzUp/(boltzUp+boltzDown) 
	otherwise set it to (0, 0.5, 1) if eUp is (positive, zero, negative).
	#
	Metropolis algorithm:
	Sets up self.MetropolisProbUp to be a 2x5 matrix of Floats
	(first index is current state of spin, second index is nUp)
	Iterate over nUp;
	  If T != 0
	   If eDown > eUp, set probability to flip from down to up to one,
	      up to down to 1-exp(-(eDown-eUp)/T)
	   If eDown <= eUp, set probability to flip from up to down to one,
	      down to up to 1-exp(-(eUp-eDown)/T)
	  Otherwise (T=0) set appropriately
	#
	Wolff algorithm
	set p
	"""
        pass

    def GetTemperature(self):
        """
	Returns self.privateTemperature
	"""
	pass

    def GetField(self):
        """
	Returns self.privateField
	"""
        pass

    def NeighborsUp(self, i, j):
        """ Sums self.lattice at four neighbor sites, modulo self.N """
        pass

    def SweepHeatBath(self, nTimes=1):
        """
	Slow variant (for debugging):
	For each time in range(ntimes):
	    For n in range(N):
	        Pick a random spin (i,j)
		Find NeighborsUp(i,j)
		Find the probability heatBathProbUp that the spin will be up
		Create a random number r in (0,1]
		if rand < heatBathProbUp, set spin lattice[i][j]=1
	#
	Fast variant:
        For each time in range(ntimes):
            Creates N random (i,j) pairs
	    Creates N random numbers in (0,1] 
		(note: use 1-RandomArray.random() to avoid zero)
	    if rand < heatBathProbUp for NeighborsUp(i,j) 
	    	set spin lattice[i][j]=1
	    else set it to zero
        """
        pass

    def SweepMetropolis(self, nTimes=1):
        """
        For each time in range(ntimes):
            Creates N random (i,j) pairs
	    Creates N random numbers in (0,1] 
	    if rand < MetropolisProbUp for current spin, NeighborsUp(i,j) 
		    set spin lattice[i][j]=1
	    else set it to zero
        """
        pass

    def WolffMoveRecursive(self):
        """
        Slow, recursive variant of Wolff move
        #
        Pick a random spin; remember its direction
        Flip it
        Call FlipNeighbors
        """
        pass

    def FlipNeighbors(self, i, j, oldSpin):
        """
        Used by WolffMoveRecursive
        #
        Initialize spinsFlipped to zero
        For m, n in neighbors of i, j:
           if lattice[m][n]==oldSpin and random()<p
       	      flip spin; add one to spinsFlipped
	      Call FlipNeighbors on (m,n); add to spinsFlipped
        return spinsFlipped
        """
        pass

    def WolffMove(self):
        """
        Faster, list-based Wolff move.
        #
        Pick a random spin; remember its direction as oldSpin
        Push it onto a list "toFlip" of spins to flip
        Set spinsFlipped = 0
        While there are spins left in toFlip
           Remove the first spin
           If it has not been flipped in between
              Flip it
              Add one to spinsFlipped
              For each of its neighbors
                  if the neighbor is in the oldSpin direction
    	          with probability p, put it on the stack
        Return spinsFlipped
        """
        pass

    def SweepWolff(self, nTimes=1, partialSweep=0):
        """
        Print error message if field is not zero
        do Wolff flips until number of spins flipped 
		(minus previous partialSweep)
        is greater or equal to N*N
        """
        pass

def runHeatBath(N=100, nSweepsPerShow=1, nShow=500):
    """
    Set up ising as an IsingModel
    Set up dl as an NxN DynamicLattice (dl = DynamicLattice.DynamicLattice(N,N))
    for t in range nShow, display lattice with dl.display(ising.Lattice)
    and then sweep with ising.SweepHeatBath(nSweepsPerShow)
    """
    pass

def runMetropolis(N=100, nSweepsPerShow=1, nShow=500):
    """
    Same as runHeatBath except use SweepMetropolis
    """
    pass

def runWolff(N=100, nSweeps=500):
    """
    Same as runHeatBath except initialize partialSweep=0,
    and use SweepWolff(partialSweep=partialSweep)
    """
    pass
