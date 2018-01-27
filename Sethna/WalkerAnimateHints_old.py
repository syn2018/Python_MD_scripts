"""Animates Walker"""
#
# See the exercise "Walker.pdf" from Walker.html
# in  http://www.physics.cornell.edu/sethna/StatMech/ComputerExercises/
#
import visual as V
import Walker as W
reload(W)	# for ipython to reload properly after changes in Walker
import scipy

V.scene.title = 'Walker'
V.scene.height = V.scene.width = 600
# Follow walker as it moves forward
V.scene.autocenter = 1
# Scale set automatically by system, but don't keep rescaling as walker moves
V.scene.autoscale = 1
# Slow down walker to reasonable speed
framesPerSecond = 120

# Walker Display

class WalkerDisplay:
    """Defines VPython cylinders for stanceLeg, swingLeg."""
    def __init__(self, w):
        """Store walker being displayed.
	Set up stanceLeg and swingLeg to be cylinders.
	stanceLeg extends from stance foot position to body position.
	swingLeg extends from body position to swing foot position."""
    	self.w = w # the walker being displayed
	d = 0.06   # thickness of legs
        pass

class CenteredWalkerDisplay:
    """Defines display centered on midpoint between stance and swing feet."""
    def __init__(self, w):
        """Store walker being displayed.
	Set up stanceLeg and swingLeg to be cylinders.
	stanceLeg extends from stance foot position to body position.
	swingLeg extends from body position to swing foot position."""
        pass

    def update(self):
        """Updates pos and axis for stanceLeg and swingLeg."""
        pass


def SillyWalk(w, dt=0.02):
    """Integrates walker differential equations, but ignores heel strikes.
    (1) Set scene.autocenter=0 to avoid following walker center of mass.
    (2) Pick times for some reasonable range (say to 100), 
        find trajectory from odeint.
    (3) Run through trajectory, setting state of system and 
        then updating display. Run rate(framesPerSecond) inside loop
	to slow things down.
    (4) Set scene.autocenter back to 1."""
    pass


def SimpleWalk(w, dwc=None, gamma=0.009, t_initial=0., t_final=20.2, dt=0.02):
    """Set gamma for walker w. Walk forward in steps of dt, 
    and update display."""
    pass


def MultiWalkFromFile(w, filename="AttractorPointsPruned.dat"):
    """Reads gammas, initial conditions from file to start on attractor"""
    pass

