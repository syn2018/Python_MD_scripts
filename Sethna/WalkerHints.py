"""Implements "The Simplest Walking Model" by Garcia, Chatterjee, 
Ruina, and Coleman."""
#
# See the exercise "Walker.pdf" from Walker.html
# in  http://www.physics.cornell.edu/sethna/StatMech/ComputerExercises/
#
from scipy.integrate import odeint
from scipy import *

class Walker:
    def __init__(self, L = 1.0,
                 theta = 0.2, thetaDot = -0.2, 
                 phi = 0.4001, phiDot = 0.0,
                 gamma = 0.009,
                 stanceFootPos = [0,0,0]):
        """Sets internal versions of these variables
        self.blah = blah
        or self.blah1, self.blah2, ... = blah1, blah2, ...
	"""
        pass

    def GetStateVector(self):
        """Returns [theta, thetaDot, phi, phiDot]"""
        pass

    def SetStateVector(self,stateVector):
        """Sets [theta, thetaDot, phi, phiDot]"""
        pass

    def GetBodyPos(self):
        """Calculate body position from stance foot position and theta
        3D vector: z-component zero"""
        pass

    def GetStanceFootPos(self):
        """Returns stance foot position"""
        pass

    def SetStanceFootPos(self,stanceFootPos):
        """Sets stance foot position"""
        pass

    def GetSwingFootPos(self):
        """Calculates swing foot position from body position, phi, and theta
        3D vector: z-component zero"""
        pass

    def CollisionCondition(self):
        """Returns condition c for heel strike"""
        pass

    def RoundCollisionConditionPositive(self):
        """Adjust phi by -2.0*c after heel strike if c<0
        Check to see that new collision condition > 0
        Prints error if c<-1.0e-10 (c should be near zero)"""
        pass

    def dydt(self, y,t):
        """Defines evolution law yDot = dydt(y,t) for odeint 
        dydt (and dzdc below) do not change the current state of the
        Walker or refer to it (except for the parameter gamma), since they
        are passed from odeint a state vector at which a derivative is to
        be evaluated."""
        pass

    def dzdc(self, z, c):
        """z = (theta, thetadot, phi, phidot, t)
        c = phi-2 theta = 0 at collision
        Define evolution law dz/dc = dzdt / dc/dt
        dzdc (and dydt above) do not change the current state of the
        Walker or refer to it (except for the parameter gamma), since they
        are passed from odeint a state vector at which a derivative is to
        be evaluated."""
        pass

    def ExecuteHeelstrike(self):
        """ExecuteHeelstrike switches roles of swing and stance legs, and
	changes momenta of legs based on impulse from heel, according
	to equation (4) in Garcia et al. above.
	ExecuteHeelstrike assumes the current state of the Walker has
        been set, and updates that state. 
	After ExecuteHeelstrike is run, phi is reset by 
	RoundCollisionConditionPositive to be big enough so that 
	CollisionCondition is just barely greater than zero."""
        pass

    def Walk(self, t_initial, t_final, dt=0.1):
        """Walk takes whatever number of Steps is needed to get from
        t_initial to t_final, stepping by dt. If Step returns a time,
	it executes a heelstrike and uses that time; otherwise it assumes
	Step ran to t_final"""
        pass

    def Step(self, t_initial, t_final, dt=0.1):
        """Step walks until the next heelstrike or t_final, whichever comes
        first. Step returns None if it reaches t_final, otherwise it 
	returns the heelstrike time.
        Step does not execute the heelstrike, so that must be called
        separately if a heelstrike condition has been reached.
	Step runs forward in time in steps of dt until either t_final or
	it passes a heelstrike (condition condition c changes from negative
	to positive). When it passes a heelstrike, it changes variables 
	from time t to c and integrates backwards to c=0. It sets the state 
	of the walker before returning."""
        pass
