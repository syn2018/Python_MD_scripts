#
# See the exercise "Repressilator.pdf" from Repressilator.html
# in  http://www.physics.cornell.edu/sethna/StatMech/ComputerExercises/
#
import scipy, scipy.integrate
import random

class Chemical:

    """Chemical is a class describing chemicals and/or chemical complexes.
    Chemical.amount: the amount of the Chemical in the cell.
    Chemical.reactions: a list of reactions whose rates are changed if
    Chemical.amount is altered."""

    def __init__(self, amount):
        # Chemical.amount represents the amount of this Chemical in the cell
        self.amount = amount
	# Chemical.reactions is a list of reactions whose rates are changed if
	# Chemical.amount is altered.
	self.reactions = []	

class DegradationReaction:

    """DegradationReaction describes the removal of one molecule of
    the specified substrate, with specified rate_constant.
    Overall rate for this reaction to occur is
    substrate.amount * rate_constant."""

    def __init__(self, substrate, rate_constant):
        self.stoichiometry = {substrate: -1}
        self.substrate = substrate
        self.rate_constant = rate_constant
	substrate.reactions.append(self)

    def GetRate(self):
        return self.substrate.amount * self.rate_constant
    
class CatalyzedSynthesisReaction:

    """CatalyzedSynthesisReaction describes the synthesis of product in
    the presence of a catalyst:  catalyst -> catalyst + product,
    with specified rate_constant.
    Overall rate for this reaction to occur is
    catalyst.amount * rate_constant."""

    def __init__(self, catalyst, product, rate_constant):
        self.stoichiometry = {product: 1}
        self.catalyst = catalyst
        self.rate_constant = rate_constant
	product.reactions.append(self)

    def GetRate(self):
        return self.catalyst.amount * self.rate_constant
    
class HeterodimerBindingReaction:

    """HeterodimerBindingReaction describes the binding of two distinct
    types of chemicals, A and B, to form a product dimer: A + B -> dimer,
    with specified rate_constant.
    Overall rate for this reaction to occur is
    A.amount * B.amount * rate_constant."""

    def __init__(self, A, B, dimer, rate_constant):
        self.stoichiometry = {A: -1, B: -1, dimer: 1}
        self.A = A
        self.B = B
        self.rate_constant = rate_constant
	A.reactions.append(self) 
	B.reactions.append(self)
        dimer.reactions.append(self) 

    def GetRate(self):
        return self.A.amount * self.B.amount * self.rate_constant

class HeterodimerUnbindingReaction:

    """HeterodimerBindingReaction describes the unbinding of a
    heterodimer into two distinct types of chemicals, A and B:
    dimer -> A + B, with specified rate_constant.
    Overall rate for this reaction to occur is
    dimer.amount * rate_constant."""

    def __init__(self, dimer, A, B, rate_constant):
        self.stoichiometry = {A: 1, B: 1, dimer: -1}
        self.dimer = dimer
        self.rate_constant = rate_constant
	A.reactions.append(self) 
	B.reactions.append(self)
	dimer.reactions.append(self) 

    def GetRate(self):
        return self.dimer.amount * self.rate_constant

# ---------------------------------------------------------------------

class Repressilator:

    """
    Repressilator is a base class for the three-gene repressilator
    system introduced by Elowitz and Leibler, ref.  All reactions
    in the stochastic version of the Repressilator are implemented, with
    rates as given in the paper.

    The set of Chemical states (all 3 mRNAs, all 3 proteins, and all 9
    promoter states) have a defined ordering, as stored in the
    dictionary chemIndex, which maps Chemicals to integer indices in
    the ordering.  This ordering is necessary to communicate with utilities
    (such as scipy.integrate.odeint) that require chemical amounts to
    be stored as a 1D array.  The GetStateVector method returns the
    Chemical amounts as a scipy array in the appropriate order, and
    the SetStateVector method sets the Chemical amounts based on a
    supplied array.
    """

    def __init__(self,
                 mRNA_degradation_rate=scipy.log(2.)/120.,
                 protein_degradation_rate=scipy.log(2.)/600.,
                 translation_rate=0.167,
                 unocc_transcription_rate=0.5,
                 occ_transcription_rate=5.0e-04,
                 P_binding=1.0,
                 P1_unbinding=224.0,
                 P2_unbinding=9.0):
        self.t = 0.
        self.chemIndex = {}
        self.reactions = []
        self.mRNA_degradation_rate = mRNA_degradation_rate
        self.protein_degradation_rate = protein_degradation_rate
        self.translation_rate = translation_rate
        self.unocc_transcription_rate = unocc_transcription_rate
        self.occ_transcription_rate = occ_transcription_rate
        self.P_binding = P_binding
        self.P1_unbinding = P1_unbinding
        self.P2_unbinding = P2_unbinding
        self.lacI = Chemical(30.0)
        self.tetR = Chemical(20.0)
        self.cI = Chemical(10.0)
        self.LacI = Chemical(300.0)
        self.TetR = Chemical(200.0)
        self.CI = Chemical(100.0)
        self.PlacI = Chemical(1.0)
        self.PtetR = Chemical(1.0)
        self.PcI = Chemical(1.0)
        self.PlacI_CI = Chemical(0.0)
        self.PtetR_LacI = Chemical(0.0)
        self.PcI_TetR = Chemical(0.0)
        self.PlacI_CI_CI = Chemical(0.0)
        self.PtetR_LacI_LacI = Chemical(0.0)
        self.PcI_TetR_TetR = Chemical(0.0)
        self.mRNAs = [self.lacI, self.tetR, self.cI]
        self.proteins = [self.LacI, self.TetR, self.CI]
        self.P0 = [self.PlacI, self.PtetR, self.PcI]
        self.P1 = [self.PlacI_CI, self.PtetR_LacI, self.PcI_TetR]
        self.P2 = [self.PlacI_CI_CI, self.PtetR_LacI_LacI, self.PcI_TetR_TetR]
        for chem in self.mRNAs:
            self.AddChemical(chem)
        for chem in self.proteins:
            self.AddChemical(chem)
        for chem in self.P0:
            self.AddChemical(chem)
        for chem in self.P1:
            self.AddChemical(chem)
        for chem in self.P2:
            self.AddChemical(chem)
        for i in range(3):
            self.AddReaction(\
                DegradationReaction(self.mRNAs[i], mRNA_degradation_rate))
            self.AddReaction(\
                DegradationReaction(self.proteins[i],protein_degradation_rate))
            self.AddReaction(\
                CatalyzedSynthesisReaction(self.mRNAs[i], self.proteins[i], \
                                           translation_rate))
            self.AddReaction(\
                CatalyzedSynthesisReaction(self.P0[i],
                                           self.mRNAs[i], 
                                           unocc_transcription_rate))
            self.AddReaction(\
                CatalyzedSynthesisReaction(self.P1[i],
                                           self.mRNAs[i], 
                                           occ_transcription_rate))
            self.AddReaction(\
                CatalyzedSynthesisReaction(self.P2[i],
                                           self.mRNAs[i], 
                                           occ_transcription_rate))
            j = (i+1)%3
            self.AddReaction(\
                            HeterodimerBindingReaction(self.P0[i],
                                           self.proteins[j],
                                           self.P1[i],
                                           P_binding))
            self.AddReaction(\
                HeterodimerBindingReaction(self.P1[i],
                                           self.proteins[j],
                                           self.P2[i],
                                           P_binding))
            self.AddReaction(\
                HeterodimerUnbindingReaction(self.P1[i],
                                             self.P0[i],
                                             self.proteins[j],
                                             P1_unbinding))
            self.AddReaction(\
                HeterodimerUnbindingReaction(self.P2[i],
                                             self.P1[i],
                                             self.proteins[j],
                                             P2_unbinding))

        self.rates = scipy.zeros(len(self.reactions), float)
	for rIndex, r in enumerate(self.reactions):
	    self.rates[rIndex] = r.GetRate()

    def AddChemical(self,chemical):
	self.chemIndex[chemical] = len(self.chemIndex)

    def GetChemicalIndex(self, chemical):
        return self.chemIndex[chemical]

    def AddReaction(self,reaction):
        self.reactions.append(reaction)

    def GetStateVector(self):
        c = scipy.zeros(len(self.chemIndex), float)        
        for chem, index in self.chemIndex.items():
            c[index] = chem.amount
        return c

    def SetFromStateVector(self, c):
        for chem, index in self.chemIndex.items():
            chem.amount = c[index]


class StochasticRepressilator (Repressilator):

    """
    StochasticRepressilator is a stochastic implementation of
    the three-gene repressilator system, with time evolution implemented
    using Gillespie's Direct Method, as described in detail in Step.
    """

    def ComputeReactionRates(self):
        """ComputeReactionRates computes the current rate for every
        reaction defined in the network, and stores the rates in self.rates."""
        for index, r in enumerate(self.reactions):
            self.rates[index] = r.GetRate()

    def Step(self, dtmax):
        """Step(self, dtmax) implements Gillespie's Direct Simulation Method,
        executing at most one reaction and returning the time increment
        required for that reaction to take place.  If no reaction is executed,
        the specified maximal time increment dtmax is returned.

        (1) all reaction rates are computed
        (2) a total rate for all reactions is found
        (3) a random time is selected, to be drawn from an exponential
            distribution with mean value given by the inverse of the total
            rate, e.g., ran_dtime = -scipy.log(1.-random.random())/total_rate
        (4) if the random time is greater than the time interval under
            consideration (dtmax), then no reaction is executed and dtmax
            is returned
        (5) otherwise, a reaction is chosen at random with relative
            probabilities given by the relative reaction rates;
            this is done by 
            (5a) uniformly drawing a random rate from the interval from
                 [0., total rate)
            (5b) identifying which reaction rate interval corresponds to
                 the randomly drawn rate, e.g., 

                 |<-------------------total rate---------------------->|
                 |<----r0----->|<-r1->|<--r2-->|<-----r3----->|<--r4-->|
                 |                                 X                   |

                 Randomly drawn rate X lands in interval r3
        (6) the chosen reaction is executed
        (7) the time at which the reaction is executed is returned
        """
        self.ComputeReactionRates()
        total_rate = sum(self.rates)
        # get exponentially distributed time
        ran_time = -scipy.log(1.-random.random())/total_rate
        if ran_time > dtmax:
            return dtmax
        # get uniformly drawn rate in interval defined by total_rate
        ran_rate = total_rate*random.random()
        # find interval corresponding to random rate
        reac_index = len(self.rates) - sum(scipy.cumsum(self.rates) > ran_rate)
        reaction = self.reactions[reac_index]
        # execute specified reaction
        for chem, dchem in reaction.stoichiometry.items():
            chem.amount += dchem
        # return time at which reaction takes place
        return ran_time

    def Run(self, T, delta_t=0.0):
        """Run(self, T, delta_t) runs the StochasticRepressilator for
        a specified time interval T, returning the trajectory at
        specified time intervals delta_t (or, if delta_t == 0.0, after
        every reaction)
        """
        tfinal = self.t + T
        if delta_t == 0.:
            ts = [self.t]
            trajectory = [self.GetStateVector()]
            while self.t < tfinal:
                dt = self.Step(tfinal-self.t)
                self.t += dt
                ts.append(self.t)
                trajectory.append(self.GetStateVector())
            return scipy.array(ts), scipy.array(trajectory)
        else:
            eps=1.0e-06
            ts = scipy.arange(0., T+eps, delta_t)
            trajectory = scipy.zeros((len(ts), len(self.chemIndex)),
                                     float)
            trajectory[0] = self.GetStateVector()
            tindex = 0
            while self.t < tfinal:
                dt = self.Step(ts[tindex+1]-self.t)
                self.t += dt
                if self.t >= ts[tindex+1]:
                    tindex += 1
                    for chem, cindex in self.chemIndex.items():
                        trajectory[tindex][cindex] = chem.amount
            return ts, trajectory
    

class DeterministicRepressilator (Repressilator):

    """
    DeterministicRepressilator is a deterministic implementation of
    the three-gene repressilator system, with time evolution implemented
    by summing up all reactions as appropriate to form a differential
    equation describing the time rate of change of all chemical
    constituents in the model.
    """

    def dcdt(self, c, t):
        """dcdt(self, c, t) returns the instantaneous time rate of
        change of the DeterministicRepressilator system, given chemical
        concentration state vector c and current time t, for use in
        integration by scipy.integrate.odeint.

        dcdt loops through all reactions defined in the Repressilator
        system, computes the rates of those reactions, and increments
        those elements in a dc_dt array that are affected by the
        reaction under consideration.

        the fully assembled dc_dt array is returned by this method.
        """
        self.SetFromStateVector(c)
	dc_dt = scipy.zeros(len(self.chemIndex), float)
	for index, r in enumerate(self.reactions):
	    rate = r.GetRate()
	    for chem, dchem in r.stoichiometry.items():
	        dc_dt[self.chemIndex[chem]] += (dchem * rate)
        return dc_dt
        
    def Run(self, tmax, dt):
        """Run(self, tmax, dt) integrates the DeterministicRepressilator
        for a time tmax, returning the trajectory at time steps as
        specified by dt, by calling scipy.integrate.odeint with the
        self.dcdt method describing the time derivative of the system

        Run should return the time array on which the trajectory is computed,
        along with the trajectory corresponding to those time points.
        """

        eps = 1.0e-06
        ts = scipy.arange(0., tmax+eps, dt)
        c = self.GetStateVector()
        traj = scipy.integrate.odeint(self.dcdt, c, ts)
        self.SetFromStateVector(traj[-1])
        return ts, traj

def RunStochasticRepressilator(T=100., dt=1., plots=False, plotPromoter=False):
    """RunStochasticRepressilator(tmax, dt, plots=False, plotPromoter=False)
    creates and runs a StochasticRepressilator for the specified time
    interval T, returning the trajectory in time increments dt,
    optionally using pylab to make plots of mRNA, 
    protein, and promoter amounts along the trajectory.
    """
    sr = StochasticRepressilator()
    sts, straj = sr.Run(T, dt)
    curvetypes = ['r-', 'g-', 'b-']
    if plots:
	import pylab
        pylab.figure(1)
        for i in range(3):
            pylab.plot(sts, straj[:,sr.chemIndex[sr.mRNAs[i]]], curvetypes[i])
        pylab.figure(2)
        for i in range(3):
            pylab.plot(sts, straj[:,sr.chemIndex[sr.proteins[i]]], curvetypes[i])
        if plotPromoter:
            pylab.figure(3)
            for i in range(3):
                promoter_state = (0.99+0.01*i)*\
                                 (straj[:,sr.chemIndex[sr.P1[i]]]\
                                  +2.*straj[:,sr.chemIndex[sr.P2[i]]])
                pylab.plot(sts, promoter_state, curvetypes[i])
        pylab.show()
    return sr, sts, straj

def RunDeterministicRepressilator(T=100., dt=1., plots=False,
                                  plotPromoter=False):
    """RunDeterministicRepressilator(tmax, dt, plots=False, plotPromoter=False)
    creates and runs a DeterministicRepressilator for the specified time
    interval T, returning the trajectory in time increments dt,
    optionally using pylab to make plots of mRNA, protein,
    and promoter amounts along the trajectory.
    """
    dr = DeterministicRepressilator()
    dts, dtraj = dr.Run(T, dt)
    curvetypes = ['r-', 'g-', 'b-']
    if plots:
        import pylab
        pylab.figure(1)
        for i in range(3):
            pylab.plot(dts, dtraj[:,dr.chemIndex[dr.mRNAs[i]]], curvetypes[i])
        pylab.figure(2)
        for i in range(3):
            pylab.plot(dts, dtraj[:,dr.chemIndex[dr.proteins[i]]], curvetypes[i])
        if plotPromoter:
            pylab.figure(3)
            for i in range(3):
                promoter_state = (0.99+0.01*i)*\
                                 (dtraj[:,dr.chemIndex[dr.P1[i]]]\
                                  +2.*dtraj[:,dr.chemIndex[dr.P2[i]]])
                pylab.plot(dts, promoter_state, curvetypes[i])
        pylab.show()
    return dr, dts, dtraj

def CompareDetSto(tmax=100., dt=1., RNAFactor=1., TelegraphFactor=1.,
                  plots=True, plotPromoter=False):
    """
    CompareDetSto compares the Deterministic and Stochastic variants of
    the Repressilator by running each, and plotting the trajectories
    alongside each other.  It is useful to return both of the Repressilator
    instances as well as both of the trajectories created.

    CompareDetSto(tmax, dt) runs each Repressilator for tmax seconds,
    returning the trajectory for each in increments of dt.

    The optional arguments RNAFactor and TelegraphFactor are used to explore
    the effects of shot noise and telegraph noise on the dynamics of
    the StochasticRepressilator; in particular,

    - initial RNA concentrations and the transcription rate are multiplied by
    RNAFactor, and the translation rate is divided by RNAFactor, thereby
    increasing RNA concentrations without otherwise affecting the continuum
    equations
    - promoter binding and unbinding rates are multiplied by TelegraphFactor
    
    The argument plots is a boolean variable affecting whether or not
    the RNA and protein concentrations for each Repressilator are plotted.

    The argument plotPromoter is a boolean variable affecting whether or not
    the various promoter states are plotted.  Because the promoter states
    switch discretely (in the StochasticRepressilator) between 0, 1, and 2,
    it is useful to separate them from each other so what the state of
    each can be ascertained.  For a trajectory straj, such a separation
    might look like:

    if plotPromoter:
       pylab.figure(3)
       pylab.plot(sts, 0.99*(straj[:,9]+2.*straj[:,12]), 'r-')
       pylab.plot(sts, 1.0*(straj[:,10]+2.*straj[:,13]), 'b-')
       pylab.plot(sts, 1.01*(straj[:,11]+2.*straj[:,14]), 'g-')
    """

    r0 = Repressilator()
    translation_rate = r0.translation_rate / RNAFactor
    occ_transcription_rate = r0.occ_transcription_rate * RNAFactor
    unocc_transcription_rate = r0.unocc_transcription_rate * RNAFactor
    P_binding = r0.P_binding * TelegraphFactor
    P1_unbinding = r0.P1_unbinding * TelegraphFactor
    P2_unbinding = r0.P2_unbinding * TelegraphFactor
    sr = StochasticRepressilator(translation_rate=translation_rate,
         occ_transcription_rate = occ_transcription_rate,
         unocc_transcription_rate = unocc_transcription_rate,
         P_binding = P_binding,
         P1_unbinding = P1_unbinding,
         P2_unbinding = P2_unbinding)
    for rna in sr.mRNAs:
        rna.amount *= RNAFactor
    sts, straj = sr.Run(tmax, dt)
    dr = DeterministicRepressilator(translation_rate=translation_rate,
         occ_transcription_rate = occ_transcription_rate,
         unocc_transcription_rate = unocc_transcription_rate,
         P_binding = P_binding,
         P1_unbinding = P1_unbinding,
         P2_unbinding = P2_unbinding)
    for rna in dr.mRNAs:
        rna.amount *= RNAFactor
    dts, dtraj = dr.Run(tmax, dt)
    if plots:
        import pylab
        pylab.figure(1)
        pylab.plot(sts, (1./RNAFactor)*straj[:,0], 'r-')
        pylab.plot(sts, (1./RNAFactor)*straj[:,1], 'b-')
        pylab.plot(sts, (1./RNAFactor)*straj[:,2], 'g-')
        pylab.plot(dts, (1./RNAFactor)*dtraj[:,0], 'r-')
        pylab.plot(dts, (1./RNAFactor)*dtraj[:,1], 'b-')
        pylab.plot(dts, (1./RNAFactor)*dtraj[:,2], 'g-')
	pylab.figure(2)
        pylab.plot(sts, straj[:,3], 'r-')
        pylab.plot(sts, straj[:,4], 'b-')
        pylab.plot(sts, straj[:,5], 'g-')
        pylab.plot(dts, dtraj[:,3], 'r-')
        pylab.plot(dts, dtraj[:,4], 'b-')
        pylab.plot(dts, dtraj[:,5], 'g-')
        if plotPromoter:
            pylab.figure(3)
            pylab.plot(sts, 0.99*(straj[:,9]+2.*straj[:,12]), 'r-')
            pylab.plot(sts, 1.0*(straj[:,10]+2.*straj[:,13]), 'b-')
            pylab.plot(sts, 1.01*(straj[:,11]+2.*straj[:,14]), 'g-')
        pylab.show()
    return sr, dr, straj, dtraj


# ---------------------------------------------------------------------

def Animate(trajectory):
    """Animate(trajectory) will animate a Repressilator trajectory using
    objects from VPython/visual to represent protein, mRNA and promoter
    states.  A trajectory should be a scipy array consisting of a number
    of timesteps, where each time step is a 15-element state vector for
    the Repressilator, in the defined ordering.

    NOTE: VPython/visual uses threads to separate control of the
    graphics from the Python interpreter.  The VPython thread seems
    to conflict with threads initiated by pylab, and can lead to
    crashing of the Python interpreter.  Therefore it is generally not a
    good idea to do both VPython animations and pylab plots within
    the same session.

    Since trajectories can take a significant amount of time to generate,
    it is sometimes useful to save a generated trajectory to a file
    (especially given the propensity of thread-induced crashes as described
    above).  Trajectories are simply arrays, so they can be written to
    and read from disk using the ArrayIO module in the Scientific.IO package:

    from Scientific.IO import ArrayIO
    # write trajectory traj to file traj.dat
    ArrayIO.writeArray(traj, 'traj.dat')    
    # read file traj.dat and store trajectory as newtraj
    newtraj = ArrayIO.readArray('traj.dat')
    """

    import RepressilatorAnimation
    RepressilatorAnimation.ClearScene()
    animator = RepressilatorAnimation.RepressilatorAnimator()
    animator.AnimateTrajectory(trajectory)

def yesno():
    response = raw_input('    Continue? (y/n) ')
    if len(response)==0: 	# [CR] returns true
	return True
    elif response[0] == 'n' or response[0] == 'N':
        return False
    else:			# Default
    	return True

def demo():
    """Demonstrates solution for exercise: example of usage"""
    print "Repressilator Demo"
    dr, dts, dtraj = RunDeterministicRepressilator(T=10000., dt=10.0)
    print "Deterministic Repressilator"
    Animate(dtraj)
    if not yesno(): return
#    print "  Compare Deterministic to Stochastic"
#    sr, dr, straj, dtraj = CompareDetSto()
    sr, sts, straj = RunStochasticRepressilator(T=200.)
    Animate(straj)

if __name__=="__main__":
   demo()


