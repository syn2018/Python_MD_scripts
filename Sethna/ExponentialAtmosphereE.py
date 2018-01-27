################################
# Exponential Atmosphere exercise
################################

################################
# Import Digital Material 
#
# ListOfAtoms
# Initializers
# Transformers
# Observers
# NeighborLocators
# BoundaryConditions
# Potentials
# Movers
#
################################

import DigitalMaterial as DM
reload(DM)

################################
# Dimension of space
################################

DM.dim = 2

################################
# Transfer other libraries from DM
################################

numpy = DM.numpy
RandomArray = DM.RandomArray
vi = DM.vi 			# Visual Python
pylab = DM.pylab

################################
# Set up pieces of simulation
################################


class LennardJonesInGravity (DM.MDSystem):
    """
    Appropriate potentials and stuff for Lennard Jones in graviational field
    """
    #
    def __init__(self, atoms=None, R=3, T=0.0, L=8., g=1.0):
	self.gravityPotential = DM.GravityPotential(g=g)
	self.LennardJonesPotential = DM.LennardJonesCutPotential()
	potential = \
		DM.CompositePotential([self.gravityPotential,
					self.LennardJonesPotential])
	boundaryConditions = DM.ReflectiveBoundaryConditions(L)
	neighborLocator = DM.SimpleNeighborLocator(
		self.LennardJonesPotential.cutoff, boundaryConditions)
	if atoms is None:
	    atoms = DM.TriangularSphericalClusterListOfAtoms(
			R=R, center=[L/2., L/2.], temperature=T,
			radius=self.LennardJonesPotential.latticeSpacing/2.0)
	else:
	    atoms = atoms
	self.displayObserver = \
	    DM.VisualDisplayAtomsObserver(atoms,L)
	self.energyObserver = DM.EnergyObserver(potential, 
			neighborLocator, boundaryConditions)
	observers = [self.displayObserver, self.energyObserver]
	mover = DM.RunVelocityVerlet;
	DM.MDSystem.__init__(self, L, atoms, observers, neighborLocator,
			     boundaryConditions, potential, mover)
    #
    def Equilibrate(self, nCoolSteps=10, temperature=0.3,
		    nStepsPerCool=50, timeStep=0.01):
	"""
	Equilibrates system, measuring and plotting temperatures 
	"""
	self.energyObserver.Reset()
	ListOfAtoms.Equilibrate(self, nCoolSteps, temperature, nStepsPerCool,
				timeStep)
	self.PlotTemperatures()
    #
    def PlotTemperatures(self):
	pylab.plot(2.0*numpy.array(self.energyObserver.KEs)/
			self.energyObserver.DOF())
	pylab.show()

def LiquidDistributions(atoms, g, T, L):
    pylab.figure(1)
    # Atoms repelled from floor by distance equal to radius
    pylab.hist([pos[1]-atoms.radius for pos in atoms.positions], 
			bins=10, normed=1.0)
    zs = numpy.arange(0,L-2.*atoms.radius, L/100.)	
    rho =(atoms.mass*g/T) * numpy.exp(-atoms.mass*g*(zs)/T)
    pylab.plot(zs, rho, "k-", linewidth=2)
    pylab.title('Height distribution')
    pylab.xlabel('Height')
    pylab.ylabel('Probability density rho(v)')
    pylab.figure(2)
    pylab.hist(atoms.velocities, bins=10, normed=1.0)
    sigmaV=numpy.sqrt(T/atoms.mass)
    vs = numpy.arange(-3.0*sigmaV, 3.0*sigmaV, sigmaV/50.)
    rho = (1.0/(numpy.sqrt(2.0*numpy.pi)*sigmaV)) \
		* numpy.exp(-0.5*atoms.mass*vs*vs/T)
    pylab.plot(vs, rho, "k-", linewidth=2)
    pylab.title('Velocity distribution')
    pylab.xlabel('Velocity v')
    pylab.ylabel('Probability density rho(v)')
    pylab.show()


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
    print "Exponential Atmosphere Demo: Liquid State"
    print "  Lennard Jones under gravity"
    sys = LennardJonesInGravity()
    sys.Run(nSteps=2000)
    if not yesno(): return
    print "  Position and velocity distributions"
    LiquidDistributions(sys.atoms, g=1.0, T=.75, L=8.)

if __name__=="__main__":
    demo()

