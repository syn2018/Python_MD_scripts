import scipy
import IsingStaticGraphics
import DynamicLattice

class IsingModel:
    """Ising model class"""

    def __init__(self, N=10, T=2./scipy.log(1.+scipy.sqrt(2.)), 
                       H=0., seed=1):
        """
	Call scipy.random.seed with argument seed;
	Set self.N to be N
	Build self.lattice to be random 50/50 0/1, using random_integers
	Call self.SetTemperatureField with T, H
	"""
        if seed==None:
	    scipy.random.seed()
	else:
            scipy.random.seed(seed)
        self.lattice = scipy.random.random_integers(0,1,(N,N))
        self.SetTemperatureField(T, H)
        self.N = N
#	# Checkerboard update
#	if N%2 != 0:
#	    print "N must be even for staggered method"
#        self.red = scipy.zeros((N,N))
#        self.black = scipy.zeros((N,N), float)
#	for i in range(N):
#            for j in range(N):
#	        if (i+j)%2 == 0.:
#		    self.red[i,j]=1.
#	else:
#	    self.black[i,j]=1.

    def SetTemperatureField(self, T, H):
        """
	Sets self.privateTemperature, self.privateField
	#
	#
	After debugging, for faster performance, set up algorithm-dependent
	data structures to avoid taking exponentials for every flipped spin
	#
	Heat bath algorithm:
	Sets up self.heatBathProbUp to be array of doubles  
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
        self.privateTemperature = T
        self.privateField = H
        J = 1.		# Convention: also float to avoid int problems
        # Set up heat bath algorithm lookup table
        self.heatBathProbUp = scipy.zeros(5, float)
        for nUp in range(0,5):	# Four neighbors on square lattice
            sumNbrs = 2*(nUp-2)	# Sum of spins of neighbors
            eUp = -J * sumNbrs - H
            eDown = J * sumNbrs + H
            if T != 0:
                boltzUp = scipy.exp(-eUp/T) 
                boltzDown = scipy.exp(-eDown/T) 
                self.heatBathProbUp[nUp] =  boltzUp/(boltzUp+boltzDown)
            else:
                if eUp>0:
                    self.heatBathProbUp[nUp]=0.
                elif eUp<0:
                    self.heatBathProbUp[nUp]=1.
                else:
                    self.heatBathProbUp[nUp]=0.5
        # Set up Metropolis algorithm lookup table
        self.MetropolisProbUp = scipy.zeros((2,5), float)
        for nUp in range(0,5):	# Four neighbors on square lattice
            sumNbrs = 2*(nUp-2)	# Sum of spins of neighbors
            eUp = -J * sumNbrs - H
            eDown = J * sumNbrs + H
	    if T != 0:
	        if eDown > eUp: # Down spin unstable
		    # If current spin is down, flip up
		    self.MetropolisProbUp[0,nUp]=1.
		    # If current spin is up, flip down with prob e^(-|dE|/T)
		    self.MetropolisProbUp[1,nUp]= 1.-scipy.exp(-(eDown-eUp)/T)
		else: # Up spin unstable
		    # If current spin is down, flip up with prob e^(-|dE|/T)
		    self.MetropolisProbUp[0,nUp]= scipy.exp(-(eUp-eDown)/T)
		    # If current spin is up, flip down
		    self.MetropolisProbUp[1,nUp]= 0.
	    else:
	        if eDown > eUp: # Down spin unstable
		    # If current spin is down, flip up
		    self.MetropolisProbUp[0,nUp]=1.
		    # If current spin is up, leave alone
		    self.MetropolisProbUp[1,nUp]= 0.
		elif eDown < eUp: # Up spin unstable
		    # If current spin is down, leave alone
		    self.MetropolisProbUp[0,nUp]= 0.
		    # If current spin is up, flip down
		    self.MetropolisProbUp[1,nUp]= 1.
        # Set up Wolff algorithm
	if T==0:
	    self.p = 0.
        else:
	    self.p = 1.0 - scipy.exp(-2.*J/T)

    def GetTemperature(self):
        """
	Returns self.privateTemperature
	"""
        return self.privateTemperature

    def GetField(self):
        """
	Returns self.privateField
	"""
        return self.privateField

    def NeighborsUp(self, i, j):
        """ Sums self.lattice at four neighbor sites, modulo self.N """
        ip1 = (i+1)%self.N
        im1 = (i-1)%self.N
        jp1 = (j+1)%self.N
        jm1 = (j-1)%self.N
        return (self.lattice[ip1][j] + self.lattice[im1][j]
                 + self.lattice[i][jp1] + self.lattice[i][jm1])

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
		(note: use 1-scipy.random.random() to avoid zero)
	    if rand < heatBathProbUp for NeighborsUp(i,j) 
	    	set spin lattice[i][j]=1
	    else set it to zero
        """
        for time in range(nTimes):
            iArr = scipy.random.randint(0, self.N, self.N*self.N)
            jArr = scipy.random.randint(0, self.N, self.N*self.N)
            randomArr = 1.-scipy.random.random(self.N*self.N)
            for i, j, rand in zip(iArr, jArr, randomArr):
                if rand < self.heatBathProbUp[self.NeighborsUp(i,j)]:
                    self.lattice[i][j] = 1
                else:
                    self.lattice[i][j] = 0

    def SweepMetropolis(self, nTimes=1):
        """
        For each time in range(ntimes):
            Creates N random (i,j) pairs
	    Creates N random numbers in (0,1] 
	    if rand < MetropolisProbUp for current spin, NeighborsUp(i,j) 
		    set spin lattice[i][j]=1
	    else set it to zero
        """
        for time in range(nTimes):
            iArr = scipy.random.randint(0, self.N, self.N*self.N)
            jArr = scipy.random.randint(0, self.N, self.N*self.N)
            randomArr = 1.-scipy.random.random(self.N*self.N)
            for i, j, rand in zip(iArr, jArr, randomArr):
                if rand < self.MetropolisProbUp[self.lattice[i][j],
	       					self.NeighborsUp(i,j)]:
                    self.lattice[i][j] = 1
                else:
                    self.lattice[i][j] = 0

    def WolffMoveRecursive(self):
        """
        Slow, recursive variant of Wolff move
        #
        Pick a random spin; remember its direction
        Flip it
        Call FlipNeighbors
        """
        i = scipy.random.randint(0,self.N)
        j = scipy.random.randint(0,self.N)
	oldSpin = self.lattice[i,j]
	self.lattice[i,j] = (self.lattice[i,j]+1)%2
	spinsFlipped = 1 + self.FlipNeighbors(i,j,oldSpin)
	return spinsFlipped

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
        spinsFlipped=0
        ip1 = (i+1)%self.N
        im1 = (i-1)%self.N
        jp1 = (j+1)%self.N
        jm1 = (j-1)%self.N
	neighbors = [(ip1,j),(im1,j),(i,jp1),(i,jm1)]
	for m,n in neighbors:
            if self.lattice[m][n] == oldSpin:
	        if scipy.random.random() < self.p:
	           self.lattice[m][n] = (self.lattice[m][n]+1)%2
	           spinsFlipped += 1+self.FlipNeighbors(m,n,oldSpin)
        return spinsFlipped

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
        i = scipy.random.randint(0,self.N)
        j = scipy.random.randint(0,self.N)
        oldSpin = self.lattice[i,j]
        toFlip = [(i,j)]
        spinsFlipped = 0
        while len(toFlip) > 0:
	    i, j = toFlip.pop(0)
	    # Check if flipped in between
	    if self.lattice[i,j] == oldSpin:
	        self.lattice[i,j] = (self.lattice[i,j]+1)%2
	        spinsFlipped += 1 
                ip1 = (i+1)%self.N
                im1 = (i-1)%self.N
                jp1 = (j+1)%self.N
                jm1 = (j-1)%self.N
	        neighbors = [(ip1,j),(im1,j),(i,jp1),(i,jm1)]
                for m, n in neighbors:
		    if self.lattice[m,n] == oldSpin:
		        if scipy.random.random() < self.p:
		            toFlip.append((m,n))
        return spinsFlipped

    def SweepWolff(self, nTimes=1, partialSweep=0):
        """
        Print error message if field is not zero
        do Wolff flips until number of spins flipped 
		(minus previous partialSweep)
        is greater or equal to N*N
        """
        if self.GetField()!=0.:
	    print "Field will be ignored by Wolff algorithm"
	partialSweep = 0
        for time in range(nTimes):
	    while partialSweep < self.N * self.N:
	        partialSweep += self.WolffMove()
            partialSweep = partialSweep - (self.N * self.N)

def runHeatBath(N=100, nSweepsPerShow=1, nShow=500):
    """
    Set up ising as an IsingModel
    Set up dl as an NxN DynamicLattice (dl = DynamicLattice.DynamicLattice(N,N))
    for t in range nShow, display lattice with dl.display(ising.Lattice)
    and then sweep with ising.SweepHeatBath(nSweepsPerShow)
    """
    ising = IsingModel(N, T=2./scipy.log(1.+scipy.sqrt(2.)) -0.1)
    dl = DynamicLattice.DynamicLattice((N,N))
    for t in range(nShow):
        dl.display(ising.lattice)
        ising.SweepHeatBath(nSweepsPerShow)

def runMetropolis(N=100, nSweepsPerShow=1, nShow=500):
    """
    Same as runHeatBath except use SweepMetropolis
    """
    ising = IsingModel(N, T=2./scipy.log(1.+scipy.sqrt(2.)) -0.1)
    dl = DynamicLattice.DynamicLattice((N,N))
    for t in range(nShow):
        dl.display(ising.lattice)
        ising.SweepMetropolis(nSweepsPerShow)

def runWolff(N=100, nSweeps=500):
    """
    Same as runHeatBath except initialize partialSweep=0,
    and use SweepWolff(partialSweep=partialSweep)
    """
    ising = IsingModel(N)
    dl = DynamicLattice.DynamicLattice((N,N))
    partialSweep=0
    for t in range(nSweeps):
        dl.display(ising.lattice)
        ising.SweepWolff(partialSweep=partialSweep)

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
    print "Ising Demo, T=Tc-0.1, N=100, nSteps=50"
    Tc = T=2./scipy.log(1.+scipy.sqrt(2.))-0.1
    N = 100
    nSweeps = 20
    ising = IsingModel(100, T=Tc)
    dl = DynamicLattice.DynamicLattice((N,N))
    print "Heat Bath Algorithm"
    for t in range(nSweeps):
        dl.display(ising.lattice)
        ising.SweepHeatBath()
    if not yesno(): return
    print "Metropolis Algorithm"
    ising = IsingModel(100, T=Tc)
    dl = DynamicLattice.DynamicLattice((N,N))
    for t in range(nSweeps):
        dl.display(ising.lattice)
        ising.SweepMetropolis()
    if not yesno(): return
    print "Wolff Algorithm"
    ising = IsingModel(100, T=Tc)
    dl = DynamicLattice.DynamicLattice((N,N))
    for t in range(nSweeps):
        dl.display(ising.lattice)
	partialSweep=0
        partialSweep = ising.SweepWolff(partialSweep=partialSweep)
       
if __name__=="__main__":
    demo()

###################################################3
# TEXTBOOK FIGURE UTILITIES
###################################################3

def CoarsenByNine(lattice):
    """
    Take square lattice of Ising spins size nxn = 3m x 3m
    Return coarse-grained matrix (majority rule) size mxm
    """
    n = len(lattice)
    if n%3 != 0:
        print "Size of matrix not divisible by 3"
    m = n/3
    coarse = scipy.zeros((m,m))
    for i in range(m):
        for j in range(m):
	    sum = 0
	    for p in range(3):
	        for q in range(3):
		    sum += lattice[3*i+p, 3*j+q] 
            if sum > 4:
	        coarse[i,j]=1
            else:
	        coarse[i,j]=0
    return coarse

def SaveCoarsenedLattices(pow=3, nSweeps=10, filenamePrefix="CoarsenedBy",
                          seed=(1,1)):
    n = 3**pow
    ising = IsingModel(n, seed=seed)
    for t in range(nSweeps):
        ising.SweepWolff()
    c = ising.lattice
    for p in range(pow):
        m = len(c)
        output = open(filenamePrefix+"%d"%p+".dat", "w")
	for i in range(m):
	    for j in range(m):
	        output.write("%d "%c[i,j])
            output.write("\n")
        output.write("\n")
        output.close()
	c = CoarsenByNine(c)
    
def SaveCoarsenedImages(pow=3, nSweeps=10, filenamePrefix="CoarsenedBy",
                          seed=(1,1)):
    n = 3**pow
    ising = IsingModel(n, seed=seed)
    for t in range(nSweeps):
        ising.SweepWolff()
    c = ising.lattice
    for p in range(pow):
        m = len(c)
	IsingStaticGraphics.DrawIsingLattice(c, 
			imfile=filenamePrefix+"%d"%p+".gif",imsize=729)
	c = CoarsenByNine(c)

# Abortive staggered method

    def SlowNeighborsUpMatrix(self):
        nUp = scipy.zeros((self.N,self.N))
	for i in range(self.N):
	    for j in range(self.N):
	        nUp[i][j] = self.NeighborsUp(i,j)
        return nUp

    def NeighborsUpMatrix(self):
        nUp = scipy.zeros((self.N,self.N))
	# Neighbor below
        nUp[0:-1] += self.lattice[1:]
	nUp[-1] += self.lattice[0]
	# Neighbor above
        nUp[1:] += self.lattice[0:-1]
	nUp[0] += self.lattice[-1]
	# Neighbor to left
        nUp[:,0:-1] += self.lattice[:,1:]
	nUp[:,-1] += self.lattice[:,0]
	# Neighbor to right
        nUp[:,1:] += self.lattice[:,0:-1]
	nUp[:,0] += self.lattice[:,-1]
	return nUp

    def SweepHeatBathStaggered(self, nTimes=1):
       for time in range(nTimes):
           randomArr = 1.-scipy.random.random((self.N,self.N))
	   # Checkerboard moves
	   # Flip even sites
	   nUp = self.NeighborsUpMatrix()
	   ## Want 
	   ## putmask(a, mask, v) results in a = v where mask is true
	   #scipy.putmask(self.lattice, red, XXX)
	   #nUp = self.NeighborsUpMatrix()
	   #putmask(self.lattice, red, XXX)
           #self.lattice = red * self.lattice 
	   #		+ black * numpy.less(randomArr, ...)

