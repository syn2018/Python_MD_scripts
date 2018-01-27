# 
# See the exercise "CardiacDynamics.pdf" from CardiacDynamics.html 
# in  http://www.physics.cornell.edu/sethna/StatMech/ComputerExercises/
#
import scipy, scipy.optimize
from DynamicLattice import *

DISPLAYINTERVAL=10

def del2_5(A, dx):
    """del2_5(A, dx) returns a finite-difference approximation of the
    laplacian of the array A, with lattice spacing dx, using the five-point
    stencil:
        0   1   0
        1  -4   1
        0   1   0
    The no-flow boundary conditions mean that the sites on the boundary
    (the "ghost cells") are not really part of the simulation: they exist
    only to set the normal derivatives to zero. Our second derivatives
    will be zero along these ghost-cells, for convenience.

    Compute the laplacian by shifting the array A appropriately (up, down,
    left, and right) to implement the stencil.  The resulting array returned
    should represent the laplacian of the array A on all of the interior
    points of the array A[1:-1, 1:-1]: your function should set del2 to 
    zeros with the shape of A, and then
       del2[1:-1,1:-1] = (A[1:-1,2:] + ...)/(dx*dx)
    """
    pass

def del2_9(A, dx):
    """del2_9(A, dx) returns a finite-difference approximation of the
    laplacian of the array A, with lattice spacing dx, using the nine-point
    stencil:
        1/6   2/3   1/6
        2/3 -10/3   2/3
        1/6   2/3   1/6
    Compute the laplacian by shifting the array A appropriately (up, down,
    left, and right) to implement the stencil.  The resulting array returned
    should represent the laplacian of the array A on all of the interior
    points of the array A[1:-1, 1:-1]
    """
    pass

del2 = del2_9

def FindFixedPoint(gamma, beta):
    """Shared with previous code FitzNag"""
    f = lambda v, gamma, beta: (v-(v**3)/3.)-((1./gamma)*(v+beta))
    vstar = scipy.optimize.brentq(f, -2., 2., args=(gamma, beta))
    wstar = ((1./gamma)*(vstar+beta))
    return vstar, wstar


def CopyGhost(A):
    """
    Given a 2D array A, create a zero-derivative boundary condition
    at all four edges of the array. This is done by changing the values
    along the boundary rows and columns to be equal to the 
    second-from-the-boundary rows and columns. For example, to ensure
    the last column of A (in python, A[:,-1]) has a zero derivative, you
    copy the second-to-the-last column: 
        A[:,-1] = A[:,-2]
    """
    pass

class FitzNag2D:
    """FitzNag2D is a class to support the simulation of the FitzHugh-Nagumo
    equations on a 2D square grid."""

    def __init__(self, N=100, dx=1.0,
                 eps=0.2, gamma=0.8, beta=0.7):
        """Initialize a FitzNag2D2D class, consisting of NxN points on a
        square grid with lattice spacing dx, and model parameters eps,
        gamma and beta; the FitzHugh-Nagumo equations are:

        dv_dt = laplacian(v) + (1./eps) * (v - (1./3.)*v**3 - w)
        dw_dt = eps*(v - gamma*w + beta)

        The __init__ method will need to store the values of the
        specified parameters, and create the fields v and w (as scipy
        arrays of the appropriate size).  v and w should be initialized
        uniformly to their values at the fixed point (resting state),
        vstar and wstar, except for the initialization of a pulse of
        height XX within a square of size 10 at the center of the
        simulation grid.

        The __init__ method should also set a member variable, e.g.,
        self.pulseheight, that will indicate the height of the pulse
        generated by interactive mouse events.

        Animated displays will be taken care of by the DynamicLattice
        class which has been imported.  The FitzNag2D class should
        initialize an instance of DynamicLattice within the __init__
        method:

        self.DL = DynamicLattice((N, N), zmin=-2.0, zmax=2.0)

        This particular instance will animate dynamics on an NxN grid,
        creating grayscale images of a specified field based on the field
        value in the interval (zmin=-2.0, zmax=2.0).
        """
        pass

    def rhs(self):
        """self.rhs() sets the instantaneous value of the right-hand-side of
        the FitzHugh-Nagumo equations, by creating two scipy arrays
        self.dv_dt and self.dw_dt, which describe the time evolution of the
        v and w fields, respectively. self.rhs() will be called by the 
	step() method as part of the Euler time-stepping scheme."""
        pass

    def step(self, dt):
        """self.step(dt) increments the fields v and w by dt.
        
	The first and last rows and columns of our arrays are "ghost cells",
	which are added to the boundaries in order to make it convenient
	to enforce the boundary conditions. Our no-flow boundary condition
	wants the normal derivative (derivative perpendicular to the 
	boundary) to equal zero. 

	step(dt) first uses CopyGhost to copy the boundary (ghost) cells 
	to ensure the zero--derivative (no-flow) boundary condition.
	It then implements an Euler step, by calling self.rhs()
	to set the current values of the fields dv_dt and dw_dt, and
	then incrementing v by dt*dv_dt and similarly for w. You
	may wish to make dv_dt and dw_dt member variables of the 
	FitzNag2D class.
	"""
        pass

    def run(self, T, dt):
        """self.run(T, dt) integrates the FitzHugh-Nagumo equations for
        a time interval T, by taking steps of size dt.   self.run() should
        also increment an overall time counter by dt.

        Every DISPLAYINTERVAL steps (e.g., 10), interaction with the
        self.DL object should be undertaken, to update the display of
        the v field, and to determine if a pulse box has been selected
        with the mouse.  The display can be updated with the
        self.DL.display() method, which takes an array to be displayed
        as an argument (e.g., self.v).  An additional nicety is to set
        the title of the self.DL window with the current simulation
        time, using the self.DL.setTitle() method.

        The self.DL object supports the query IsBoxSelected(), to ascertain
        whether a region (box) of the grid has been selected with the mouse.
        If a box has been selected, the self.DL.GetMouseBox() can be
        called to unpack a set of grid points (x0, y0, x1, y1) delineating
        the selected box.  The v field within the selected box should be
        set to the pre-selected pulseheight.
        """
        pass
