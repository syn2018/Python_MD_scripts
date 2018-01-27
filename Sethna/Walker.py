#
# See the exercise "Walker.pdf" from Walker.html
# in  http://www.physics.cornell.edu/sethna/StatMech/ComputerExercises/
#
"""Implements "The Simplest Walking Model" by Garcia, Chatterjee, 
Ruina, and Coleman."""
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
	self.L, self.theta, self.thetaDot, self.phi, self.phiDot, self.gamma \
		= L, theta, thetaDot, phi, phiDot, gamma
	self.stanceFootPos = array(stanceFootPos)

    def GetStateVector(self):
        """Returns [theta, thetaDot, phi, phiDot]"""
        return array([self.theta, self.thetaDot, self.phi, self.phiDot])

    def SetStateVector(self,stateVector):
        """Sets [theta, thetaDot, phi, phiDot]"""
        self.theta, self.thetaDot, self.phi, self.phiDot = stateVector

    def GetBodyPos(self):
        """Calculate body position from stance foot position and theta
        3D vector: z-component zero"""
    	return self.stanceFootPos \
	    + [-self.L*sin(self.theta), self.L*cos(self.theta), 0]

    def GetStanceFootPos(self):
        """Returns stance foot position"""
    	return self.stanceFootPos

    def SetStanceFootPos(self,stanceFootPos):
        """Sets stance foot position"""
        self.stanceFootPos = stanceFootPos

    def GetSwingFootPos(self):
        """Calculates swing foot position from body position, phi, and theta
        3D vector: z-component zero"""
        return self.GetBodyPos() \
	    + [-self.L*sin(self.phi-self.theta), 
	       -self.L*cos(self.phi-self.theta), 0]

    def CollisionCondition(self):
        """Returns condition c for heel strike"""
	return self.phi - 2.*self.theta

    def RoundCollisionConditionPositive(self):
        """Adjust phi by -2.0*c after heel strike if c<0
        Check to see that new collision condition > 0
        Prints error if c<-1.0e-10 (c should be near zero)"""
	c = self.CollisionCondition()
	if c<0:
	    if c < -1.0e-10:
	        print "Unexpectedly large negative collision condition = ", c, \
			" after heelstrike"
	    self.phi = self.phi - 2.0*c
	    if self.CollisionCondition() < 0:
	        print "Problems with collision condition"

    def dydt(self, y,t):
        """Defines evolution law yDot = dydt(y,t) for odeint 
        dydt (and dzdc below) do not change the current state of the
        Walker or refer to it (except for the parameter gamma), since they
        are passed from odeint a state vector at which a derivative is to
        be evaluated."""
        theta,thetaDot,phi,phiDot = y
        thetaDotdot = sin(theta-self.gamma)
        phiDotdot = thetaDotdot + \
                    (thetaDot**2)*sin(phi)-cos(theta-self.gamma)*sin(phi)
        return [thetaDot,thetaDotdot,phiDot,phiDotdot]

    def dzdc(self, z, c):
        """z = (theta, thetadot, phi, phidot, t)
        c = phi-2 theta = 0 at collision
        Define evolution law dz/dc = dzdt / dc/dt
        dzdc (and dydt above) do not change the current state of the
        Walker or refer to it (except for the parameter gamma), since they
        are passed from odeint a state vector at which a derivative is to
        be evaluated."""
        theta,thetaDot,phi,phiDot,t = z
        y = array([theta, thetaDot, phi, phiDot])
        thetaDot, thetaDotdot, phiDot, phiDotdot = self.dydt(y, t)
        cDot = phiDot - 2.*thetaDot
        return [thetaDot/cDot,thetaDotdot/cDot,phiDot/cDot,phiDotdot/cDot,
                1./cDot]

    def ExecuteHeelstrike(self):
        """ExecuteHeelstrike switches roles of swing and stance legs, and
	changes momenta of legs based on impulse from heel, according
	to equation (4) in Garcia et al. above.
	ExecuteHeelstrike assumes the current state of the Walker has
        been set, and updates that state. 
	After ExecuteHeelstrike is run, phi is reset by 
	RoundCollisionConditionPositive to be big enough so that 
	CollisionCondition is just barely greater than zero."""
        # Heelstrike exchanges swing and stance legs
        self.SetStanceFootPos(self.GetSwingFootPos())
	# Copy pre-heelstrike variables
        theta, thetaDot, phi, phiDot = \
               self.theta, self.thetaDot, self.phi, self.phiDot
        # Heelstrike changes momenta, exchanges legs
        self.theta = -theta # Stuck leg has angle theta-phi = -theta
        self.thetaDot = cos(2.*theta) * thetaDot # From eq. 4
        self.phi = -2.*theta # Exchange of swinging leg
        self.phiDot = cos(2.*theta) * (1.-cos(2.*theta)) * thetaDot # eq. 4
	self.RoundCollisionConditionPositive()
	# Alternative method: use matrix multiplication
        #strikeMatrix = array([[-1., 0., 0., 0.],
        #                     [0, cos(2.*z[0]), 0., 0.],
        #                     [-2., 0., 0., 0.],
        #                     [0., cos(2.*z[0])*(1.-cos(2.*z[0])), 0., 0.]])
        #ynew = dot(strikeMatrix, z[0:4])
        #w.SetStateVector(ynew)

    def Walk(self, t_initial, t_final, dt=0.1):
        """Walk takes whatever number of Steps is needed to get from
        t_initial to t_final, stepping by dt. If Step returns a time,
	it executes a heelstrike and uses that time; otherwise it assumes
	Step ran to t_final"""
        t_cur = t_initial
        while t_cur < t_final:
            t_step = self.Step(t_cur, t_final, dt)
            if t_step is None:
                t_cur = t_final   # non heelstrike return from Step
            else:
                t_cur = t_step
                self.ExecuteHeelstrike()

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
        old_CollisionCondition = self.CollisionCondition()
	t = t_initial
        while t < t_final:
            tf = min(t+dt, t_final)
            self.SetStateVector(odeint(self.dydt, self.GetStateVector(),
                                       array([t_initial,tf]))[-1])
            if old_CollisionCondition < 0 and self.CollisionCondition() >= 0:
                # Run backward to heelstrike
                y = self.GetStateVector()
                z0 = [y[0],y[1],y[2],y[3],t]
                z=odeint(self.dzdc,z0,array([self.CollisionCondition(),0.]))[-1]
                self.SetStateVector([z[0],z[1],z[2],z[3]]) # before heelstrike
                t = z[4]
                return t
            old_CollisionCondition = self.CollisionCondition()
	    t = tf
        return None

def test():
    w = Walker()
    dt = 0.1
    t_initial = 0.
    t_final = 20.
    w.Walk(t_initial, t_final, dt)
    return w

