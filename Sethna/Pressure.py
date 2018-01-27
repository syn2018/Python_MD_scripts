################################
# Pressure exercise
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

DM.dim = 3

################################
# Transfer other libraries from DM
################################

numpy = DM.numpy
vi = DM.vi 			# Visual Python
pylab = DM.pylab

################################
# Set up pieces of simulation
################################

# XXX Make box of size L+2*radius, so V = LxLxL

class PressureMeasurement (DM.MDSystem):
    """
    Appropriate potentials and stuff for pressure measurements from wall
    collisions
    """
    #
    def __init__(self, atoms=None, nAtoms=500, T=10.0, L=20.):
	potential = \
		DM.Potential()
	boundaryConditions = DM.ReflectiveBoundaryConditions(L,
							impulseRecording=True)
	neighborLocator = DM.NoNeighborLocator(potential.cutoff,
							boundaryConditions)
	if atoms is None:
    	    atoms = DM.RandomListOfAtoms(L, nAtoms=nAtoms, temperature=T)
	self.displayObserver = DM.VisualDisplayAtomsObserver(atoms,L)
	observers = [self.displayObserver]
	mover = DM.RunVelocityVerlet;
	DM.MDSystem.__init__(self, L, atoms, observers, neighborLocator,
			     boundaryConditions, potential, mover)
    #
    def PlotImpulses(self):
    	"""
	Plots histograms of impulses Delta p on the wall due to specular
	collisions of gas molecules.
	"""
	impulses = []
	for d in range(DM.dim):
	    negImpulses = -numpy.array(
	    			self.boundaryConditions.impulses[d][0])
	    impulses.extend(negImpulses)
	    impulses.extend(self.boundaryConditions.impulses[d][1])
	pylab.hist(impulses, bins=50)
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
    print "Pressure Measurement Demo"
    sys = PressureMeasurement()
    sys.Run()
    if not yesno(): return
    print "  Impulse Histogram"
    sys.PlotImpulses()
    

if __name__=="__main__":
    demo()

