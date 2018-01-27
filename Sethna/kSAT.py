import numpy, numpy.random, pylab, random
# Graphics for animating solution process
import DynamicNetwork
reload(DynamicNetwork)

# Use to control whether you animate your solution
# Set to False if you want to run quickly
GRAPHICS = True

# We first generate some "real-world" logical satisfiability problems
# that emerge from graph coloring. This is mostly to give a feeling for
# what SAT problems are, and to give a couple of examples where we
# know the answer, to be used to debug our algorithms.

def MakeSATFrom3Colorability(colorGraphDict):
    """Transform a 3-colorability problem into a SAT exercise.
    Assume that the graph to be colored is a dictionary, with
    keys labeling the nodes of the graph, and lists of neighbors
    as edges for each node. We assume an undirected graph, so if 
    A has B in its list of neighbors, so B has A in its list.
    Assume colorGraphDict gives the neighbors of each edge:
	    ['A'] = ['B', 'C', 'D']
    Given an N-node graph, make up a 3N-variable SAT problem
    where node A has three logical variables A_red, A_green, and A_blue.
    # 
    Let's do this systematically. 
    (1) Define three dictionaries red, green, and blue, by iterating
	over the nodes ("for node in colorGraphDict"); 
	for the nth node [counting from zero] red=3*n+1, green=3*n+2, 
	and blue=3*n+3. (We'll use the negatives of these variables to denote
	the logical statement that the node is not colored that color:
	not-red = -(3*n+1), etc.)
	Each clause is a list of three integers (minus for "not"), so
	[6, -7, 3] says "node 2 is blue or node 3 is not red or node 1 is blue".
    (2) Then generate the clauses, appending to an empty list
	For each node, 
	 Append a clause saying that the node is colored at least one color
	 Append three clauses saying that the node is not colored two colors
	 For each edge (say between A and B), if B>A (to avoid double counting)
	   append three clauses saying that 
	   A and B are not colored the same color
    (3) Return the clauses
    """
    red = {}
    green = {}
    blue = {}
    nodeNumber = {}
    variables = []
    #
    n=0
    for node in colorGraphDict:
       red[node] = 3*n+1
       green[node] = 3*n+2
       blue[node] = 3*n+3
       n+=1
    # 
    clauses = []
    for node in colorGraphDict:
	# At least one color
	clauses.append([red[node], green[node], blue[node]]) 	
	clauses.append([-red[node], -green[node]])  # Not both red and green
	clauses.append([-red[node], -blue[node]])   # Not both red and blue
	clauses.append([-green[node], -blue[node]]) # Not both green and blue
	for node2 in colorGraphDict[node]:
	    if node2 > node:
		clauses.append([-red[node], -red[node2]])     # Not both red 
		clauses.append([-green[node], -green[node2]]) # Not both green 
		clauses.append([-blue[node], -blue[node2]])   # Not both blue 
    #
    return clauses
    

def MakeMapsFromProblemSet():
    """Generate, by hand, the two maps in the Graph Coloring figure
    in the Satisfactory Map Coloring exercise, as dictionaries:
    leftGraph={}
    leftGraph['A']=['B','C','D'], 
    ...
    return leftGraph, rightGraph
    """
    leftGraph={}
    leftGraph['A']=['B','C','D']
    leftGraph['B']=['A','C']
    leftGraph['C']=['A','B','D'] 
    leftGraph['D']=['A','C']
#
    rightGraph={}
    rightGraph['A']=['B','C','D'] 
    rightGraph['B']=['A','C','D'] 
    rightGraph['C']=['A','B','D'] 
    rightGraph['D']=['A','B','C'] 
#
    return leftGraph, rightGraph

def MakeRandom_kSATClauses_0(nVars, nClauses, k=3):
    """
    MakeRandom_kSATClauses(nVars, nClauses, k) generates a random k-SAT
    instance with nvariables and nClauses.
    We recommend using numpy.random.randint, making an array of shape 
    (nClauses,k) of integers in the range [1,nVars], and multiplying it
    by an array filled with +-1 (found by making an array of 0's and 1's,
    doubling it and subtracting one).
    Our clauses are lists and not arrays (because lists can have variable
    numbers of entries, mixing 2SAT and 3SAT, for example), so convert
    the result to a list using clauses.tolist().
    """
    vars = numpy.random.randint(1,nVars+1,(nClauses,k))
    signs = numpy.random.randint(0,2,(nClauses,k))*2-1
    return (signs * vars).tolist()

# We also want to generate random kSAT problems

def MakeRandom_kSATClauses(nvars, nclauses, k=3):
    """make_random_SAT_clauses(nvars, nclauses, k) generates a random k-SAT
    instance with nvariables and nclauses; each clause is guaranteed not
    to contain a given variable (or its negation) more than once"""
    signs = 2*numpy.random.randint(0,2, (nclauses,k))-1
    v = []
    for nc in range(nclauses):
        has_repeats = True
        while has_repeats:
            vc = numpy.random.randint(1, nvars+1, (k,))
            s = set(vc)
            if len(s) == k:
                has_repeats = False
        v.append(vc)
    if len(v) > 0:
        v = numpy.array(v)
    else:
        v = numpy.zeros((0,k), 'd')
    clauses = signs * numpy.array(v)
    return clauses.tolist()




# -------------------------------------------------------------------

class DPSATSolverBase:
	
    """
    SATSolverBase implements a Davis-Putnam algorithm for finding solutions
    to logical satisfiability problems stated in conjuctive normal form. 
    #
    This is a "base class", or virtual class: it does not implement
    one key method (TrySet), but awaits two different implementations
    (SlowDPSATSolver and FasterDPSATSolver) in classes derived from 
    SATSolverBase.
    """

    def __init__(self, clauseList=None):
	"""
	Constructs a SATSolverBase given a set of input clauses.
	It needs to create and store the following data members:
	#
	self.clauseList: the list of input clauses, each one consisting of
	a list of variables and negations of variables (stored as a list
	of signed integers)
	#
	self.vars: a list of the variables. This internal list will be 
	used also to determine the order in which the variables are 
	tried (given tentative values). Shuffling or reordering this 
	list will be an important diagnostic and optimization tool.
	Until we start the fancy methods, be sure to sort this list
	(self.vars.sort()) after it is assembled from the 
	clauses, so that your algorithm will set variables in the same
	order as the graphics program lays them out.
	#
	self.varStatus: a dictionary mapping variables to boolean
	values True or False, or None if variable is (so far) unset
	#
	self.nodeClauses: a dictionary mapping a variable
	to the indices of the clauses including that variable and maps 
	its negation to those for clauses including the negation of that
	variable. Specifically, if clauseList[i]=[n,-m,p] then 
	nodeClauses[n]=[..., i, ...] contains i, as should nodeClauses[-m]
	and nodeClauses[p].
	#	
	self.clauseStatus: a dictionary identifying
	whether the clause is not determined (None), or if it is True
	give the current state of the variables. (If a clause
	is False, for this algorithm, we immediately backtrack.)
	#
	if GRAPHICS is set to True, create 
	self.dynnet: an instance of a DynamicNetwork graphics window 
	for animating solver progress 
	(self.dynnet = DynamicNetwork.DynamicNetwork())
	"""
	self.clauseList = clauseList
	self.vars = []
	self.varStatus = {}
	self.nodeClauses = {}
	self.clauseStatus = {}
	if GRAPHICS: self.dynnet = DynamicNetwork.DynamicNetwork()
	#
	for clauseIndex, clause in enumerate(clauseList):
	    for literal in clause:
		var = abs(literal)
		if var not in self.vars:
		    self.vars.append(var)
		    self.varStatus[var]=None
		    self.nodeClauses[var] = []
		    self.nodeClauses[-var] = []
		self.nodeClauses[literal].append(clauseIndex)
	    self.clauseStatus[clauseIndex]=None
	#self.vars.sort()
	numVarOccur = numpy.zeros(max(self.vars)+1)
	for c in self.clauseList:
	    for lit in c:
	       numVarOccur[abs(lit)]+=1
        self.vars.sort(cmp = lambda x, y: cmp(numVarOccur[y],numVarOccur[x]))

    def SelectNextVariableToTry(self):
        """Gives next variable (index entry into list vars) for DP algorithm
	to try setting. Initially we explore the simple algorithm of
	just setting them in order: later we'll select them in MOMS order"""
        #return self.SelectNextVariableSimple()
	return self.SelectNextVariableMOMS()
    
    def SelectNextVariableSimple(self):
        """
	Search through self.vars for the first unset variable
	(the index for which varStatus[vars[varIndex]] is None). If all are
	set (we're done!) return None; otherwise return varIndex.
	"""
	varIndex = 0
	while self.varStatus[self.vars[varIndex]] is not None:
	    varIndex+=1
	    if varIndex == len(self.vars):
	        return None
	return self.vars[varIndex]

    def SelectNextVariableMOMS(self):
        """
	Select one of the MOMS literals (Most Occurences in unresolved
	clauses of Minimial Size) at random. Also, if there are no
	unresolved clauses, return None to signal that a SAT solution
	has been found.
	"Size" refers to the number of unset variables in the clause.
	Usually the minimal size will be two. (`Length one' literals will
	be reduced by FasterDPSATSolver.TrySet).
	(1) Assemble a list MSVars of the unset variables in unresolved clauses
	    of minimal size:
	    Set k=2; 
	    while no unresolved clauses of length k exist:
	      Loop over clauses,
	        if the clause is not True (resolved)
		    build a list of unset variables in the clause
		    if that list has k variables, 
		       note that length k clauses exist
		       extend the list of MSVars
	       If no unresolved clauses exist, return None
	       Increment k
	(2) Find the number of times each variable is contained in MSVars,
	and the maximum number
	(3) Collect the list of MOMS variables represented most often
	(4) Use random.choice to return random choice among MOMS variables
	"""
        MSVars = []
	unresolvedClausesExist = False
	length_k_ClausesExist = False
	k=2
	while not length_k_ClausesExist:
            for clauseIndex, clause in enumerate(self.clauseList):
	        if self.clauseStatus[clauseIndex]!=True:
	            unresolvedClausesExist=True
		    unresolvedVars = [abs(lit) 
		    			for lit in self.clauseList[clauseIndex]
					    if self.varStatus[abs(lit)]==None]
	            if len(unresolvedVars)==k:
		       length_k_ClausesExist = True
		       MSVars.extend(unresolvedVars)
	    if not unresolvedClausesExist:
	        # No unresolved clauses! Done!
	        return None
	    k+=1
	# Find number of times each variable is contained in MSVars    
	variableCounts = numpy.zeros(max(self.vars)+1)
	for var in MSVars:
	   variableCounts[var] += 1 
        mostOccurences = max(variableCounts)
	MOMSVariables = []
	for var in range(max(self.vars)+1):
	    if variableCounts[var] == mostOccurences:
	        MOMSVariables.append(var)
        return random.choice(MOMSVariables)

    def SelectNextVariableMOMS2(self):
        """
	Select one of the MOMS literals (Most Occurences in clusters of
	Minimial Size) at random.
	(1) Find the size of the smallest currently unsatisfied clauses 
	(usually 2 for 3SAT)
	(2) Find the number of times each unset literal is contained in those
	clauses
	(3) Collect the list of literals represented most often
	(4) Use random.choice to return random choice among most used literals
	"""

	unresolved = [(clauseIndex, 
			[self.varStatus[abs(lit)] for lit in clause])
			  for clauseIndex, clause in enumerate(self.clauseList)
			    if self.clauseStatus[clauseIndex]!=True]
	if len(unresolved)==0:
	    # No unresolved clauses! Done!
	    return None
	numNones = [lis.count(None) for i, lis in unresolved]
	minNones = min(numNones)
	MSVars = []
	for clauseIndex, vals in unresolved:
	    if vals.count(None)==minNones:
	        vars = [abs(lit) for lit in self.clauseList[clauseIndex]
			    if self.varStatus[abs(lit)]==None]
		MSVars.extend(vars)
	
	#clauseLengthBound = 100000000000	# Assume no very long clauses 
	#clauseLengthMinimum = clauseLengthBound
	#MSLiterals = []
	#for clauseIndex, clause in enumerate(self.clauseList):
	#    if self.clauseStatus[clauseIndex]!=True: # Ignore satisfied clauses
	#        clauseLength = 0
	#	unsetLiterals = []
	#	for literal in clause:
	#	    if self.varStatus[abs(literal)] == None: 
	#	        clauseLength+=1
	#	        unsetLiterals.append(literal)
	#	if clauseLength <= 1:
	#	    print "clauseLength too small in SelectNextVariableMOMS"
	#	elif clauseLength == clauseLengthMinimum:
	#	    MSLiterals.extend(unsetLiterals)
	#	elif clauseLength < clauseLengthMinimum:
	#	    clauseLengthMinimum=clauseLength
	#	    MSLiterals = unsetLiterals
	#if clauseLengthMinimum == clauseLengthBound:
	#    # No unresolved clauses! Done!
	#    return None
	## Now find maximum, and collect all variables with that maximum
	
	variableCounts = numpy.zeros(max(self.vars)+1)
	for literal in MSLiterals:
	   variableCounts[abs(literal)] += 1 
        mostOccurences = max(variableCounts)
	MOMSVariables = []
	for var in range(max(self.vars)+1):
	    if variableCounts[var] == mostOccurences:
	        MOMSVariables.append(var)
        return random.choice(MOMSVariables)

    def Solve(self):
	"""Solve() returns True if the clauses are satisfiable (SAT)
	and False if UNSAT, along with the number of tries required
	to find a solution. Solve() also leaves "self" in a solved
	state if it found one.
	#
	Solve works by 
	(0) initializing numTries = 0
	(1) selecting a literal litNext to try setting next. If all
	variables are set, declare victory: return True and numTries(=0)
	(2) loop over literal in [litNext, -litNext]
	      calling TrySet on that literal, incrementing numTries
	      If TrySet returns True, Solve recursively calls itself,
		and supplements numTries.
	    	  If this recursive call returns true, Solve declares 
		  victory and returns True and numTries
	      Otherwise, unwind: set the status of the variables in varsSet
	      and clauses in clausesSet to None
	(3) return False, numTries
	Also, if GRAPHICS==True, call Display() after each call to 
	TrySet or Solve.
	"""
	numTries = 0
	var = self.SelectNextVariableToTry()
	if var==None:
	    return True, numTries
	#
	for literal in (var, -var):
	    trySuccess, varsSet, clausesSet = self.TrySet(literal)
	    if trySuccess:
		solveSuccess, incNumTries = self.Solve()
		if GRAPHICS: self.Display()
		numTries += incNumTries
		if solveSuccess: 
		    # Success! Declare victory and leave.
	    	    if GRAPHICS: self.Display()
		    return True, numTries
	    # Unwind
	    numTries += 1
	    for vars in varsSet:
		self.varStatus[vars]=None
	    for clauseIndex in clausesSet:
		self.clauseStatus[clauseIndex]=None
	    if GRAPHICS: self.Display()
	# 
	return False, numTries 

    def TrySet(self, literal):
	"""Tries to set one literal to True, with various ramifications
	setting other variables depending on implementation. Returns 
	status (True or False), a list varsSet of variables that have been
	set, and a list clausesSetTrue of clauses that have been newly set
	to True.
	#
	TrySet is undefined in the base class: SlowDPSATSolver and 
	FasterDPSATSolver will implement it as subclasses of SATSolverBase.
	"""
	pass
	
    def Display(self):
	"""Display() displays the current status of the DPSATSolver by
	calling self.dynnet.displayFromLists with the following data:
	#
	nodelist: the list of variables
	edgelist: the list of edges between variables which occur in
		  clauses that have not been satisfied
	nodecolors: a dictionary mapping (positive) variable ids to colors
	clear = True: otherwise it won't redraw when edges are removed
	"""
	self.dynnet.displayFromLists(self.vars, self.GetUnresolvedEdges(),
	    self.GetNodeColors(), clear=True)

    def GetUnresolvedEdges(self):
	"""GetUnresolvedEdges() (used in Display) returns a list of all 
	pairs of variables that are unset, and are involved in a clause
	which is not yet resolved; 
	start with an empty list of unresolved edges,
	for each clauseIndex, clause in enumerate(clauseList)
	   if the clause status is not set
	      for each literal in the clause
		 if the variable is unset
		    for each second literal in the clause
		       if the second variable is > first variable and is unset
			   append to unresolved"""
	unresolved = []
	for clauseIndex, clause in enumerate(self.clauseList):
	   if self.clauseStatus[clauseIndex] is None:
	      for lit1 in clause:
		    v1 = abs(lit1)
		    if self.varStatus[v1] is None:
		       for lit2 in clause:
			  v2 = abs(lit2)
			  if v1>v2 and self.varStatus[v2] is None:
			     unresolved.append((v1,v2))
	return unresolved

    def GetNodeColors(self):
	"""GetNodeColors() returns a dictionary of colors for Display:
	green (0,255,0) if variable is True, red (255,0,0)
	if False, and white (255,255,255) if unset"""
	nodeColors = {}
	for var in self.vars:
	    if self.varStatus[var]==True: nodeColors[var]=(0,255,0)
	    if self.varStatus[var]==False: nodeColors[var]=(255,0,0)
	    if self.varStatus[var]==None: nodeColors[var]=(255,255,255)
	return nodeColors

    """You'll want to implement and test SlowDPSATSolver and Display
    on your colorability problems and some kSAT examples at this point, and
    then PlotSATFrac. Probably then you'll be motivated to improve on your
    algorithm using FasterDPSATSolver. At that point, you'll want to 
    start implementing these further class member functions."""

    def ShuffleVariableOrder(self):
	"""ShuffleVariableOrder() resets and then
	randomly shuffles the var_order list,
	using the random.shuffle method: used for diagnosis on 
	FasterDPSATSolver"""
	self.Reset()
	random.shuffle(self.vars)

    def Reset(self):
	"""Reset() sets each variable value to None and each element of
	clause_true to None. Used especially when repeated runs are needed,
	as in using ShuffleVariableOrder."""
	for var in self.vars:
	    self.varStatus[var]=None
	for clauseIndex in range(len(self.clauseList)):
	    self.clauseStatus[clauseIndex]=None

    def SetVariableOrderFromFrequency(self):
	"""SetVariableOrderFromFrequency() sorts the vars list
	by the frequency with which each variable (or its negation)
	occurs in the clauses, in order of decreasing frequency"""
	pass

    def Verify(self):
	"""
	Verify() returns True if the current settings of the variables
	satisfy all the SAT clauses; otherwise returns False. Used especially
	in debugging.
	"""
	for clause in self.clauseList:
	    SAT = False
	    for literal in clause:
		if literal > 0 and self.varStatus[literal]==True: SAT = True
		if literal < 0 and self.varStatus[-literal]==False: SAT = True
	    if SAT==False: 
		return False
	return True
    

class SlowDPSATSolver (DPSATSolverBase):

    def TrySet(self, literal):
	"""Try setting one literal true.
	(1) Sets self.varStatus[var] to True or False, 
	    where var = abs(literal) and the sign of the literal is 
	    positive for True and negative for False. 
	    Set varsSet=[var] and initialize clausesSet to [].
	(2) Then tests each clause involving -literal: 
	      if it is not already set True, 
		if all of its literals have variables that are set, 
		  (so the clause must be False) 
		  return False, varsSet=[var], and clausesSet = []
	(3) Then, for each clause involving the literal, if the clause
	    is unset, add the clause index to clausesSet and set the
	    clause status to True.
	    Return True, varsSet, and clausesSet"""
	#
	var = abs(literal)
	if literal>0:
	    self.varStatus[var] = True
	else:
	    self.varStatus[var] = False
	varsSet = [var]
	clausesSet = []
	# Negative clauses
	for clauseIndex in self.nodeClauses[-literal]:
	    if not self.clauseStatus[clauseIndex]: # Don't worry if already true
		varStatus = [self.varStatus[abs(lit)] 
	    			for lit in self.clauseList[clauseIndex]]
		if None not in varStatus:  # All variables set `wrong'
		    return False, varsSet, clausesSet
	# Positive clauses
	for clauseIndex in self.nodeClauses[literal]:
	    if self.clauseStatus[clauseIndex] is None:
		clausesSet.append(clauseIndex)
		self.clauseStatus[clauseIndex] = True
	return True, varsSet, clausesSet


class FasterDPSATSolver (DPSATSolverBase):

    def TrySet(self, literal):
	"""Improves on SlowDPSATSolver, by looking for a `length one'
	clause where all but one literal is set wrong, and using
	TrySet recursively to set that last literal true.
	"""
	# XXX Code duplication: should refactor?
	var = abs(literal)
	if literal>0:
	    self.varStatus[var] = True
	else:
	    self.varStatus[var] = False
	varsSet = [var]
	clausesSet = []
	# Negative clauses
	for clauseIndex in self.nodeClauses[-literal]:
	    if not self.clauseStatus[clauseIndex]: # Don't worry if already true
		varStatus = [self.varStatus[abs(lit)] 
	    			for lit in self.clauseList[clauseIndex]]
		if None not in varStatus:  # All variables set `wrong'
		    return False, varsSet, clausesSet
	# Positive clauses
	for clauseIndex in self.nodeClauses[literal]:
	    if self.clauseStatus[clauseIndex] is None:
		clausesSet.append(clauseIndex)
		self.clauseStatus[clauseIndex] = True
	# XXX Code duplication ends here
	# Test for 'length one' clauses
	for clauseIndex, clause in enumerate(self.clauseList):
	    if not self.clauseStatus[clauseIndex]: # Don't worry if already true
		unsetLiterals = [lit for lit in clause 
					if self.varStatus[abs(lit)]==None]
		if len(unsetLiterals)==1:
		   status, tryVarsSet, tryClausesSet = \
		   	 self.TrySet(unsetLiterals[0])
		   varsSet.extend(tryVarsSet)
		   clausesSet.extend(tryClausesSet)
		   if status==False:
		      return False, varsSet, clausesSet
		   else:	# Recursive: must be finished
		      return True, varsSet, clausesSet
	return True, varsSet, clausesSet

# -------------------------------------------------------------------
# Various SAT analyses and tests
# -------------------------------------------------------------------

def SATFrac(nVars, nClauses, nRuns, k=3, Solver=FasterDPSATSolver):
    """SATFrac(nVars, nClauses, nRuns, k=3, Solver=FasterDPSATSolver)
    calculates the fraction satisfiable and the times taken for 
    nRuns of a kSAT problem with nVars and nClauses, using the solver algorithm.
    Returns mean and standard deviation for SATFrac and time, and the list
    of times.
    #
    It sets GRAPHICS to False (after declaring it global).
    It sets up variables 'SATTot' to store the total numer of satisfiable 
    problems, 'times' a list to store the number of tries taken by the
    algorithm, 'timeTot' for the total time, and 'time2Tot' for the
    total squared time. [Remember that the variance 
    var=(t2Tot-tTot^2/N)/(N-1) and thus the standard deviation of 
    the mean sig=sqrt(var/N)].
    It then loops over nRuns
	makes clauses
	instantiates solver
	finds SAT, time for solution
	(option, verify solution: useful for debugging new solvers)
	increments SATTot, timeTot, time2Tot, and appends to times 
    finds SATFrac and timeAv
    if nRuns==1 returns them
    finds variance, sigma for the mean for SAT, time 
      (note SAT2Tot = SATTot: why?)
    returns SATFrac, SATFracSigmaMean, timeAv, timeSigmaMean
    """
    global GRAPHICS
    GRAPHICS=False
    SATTot = 0.
    timeTot = 0.
    time2Tot = 0.
    times = []
    #
    for run in range(nRuns):
	clauses = MakeRandom_kSATClauses(nVars, nClauses, k)
	sol = Solver(clauses)
	SAT, time = sol.Solve()
	if SAT != sol.Verify():
	   print "Error in SATFrac: solution doesn't check"
	if SAT: SATTot += 1
	timeTot += time
	time2Tot += time*time
	times.append(time)
    SATFrac = SATTot / nRuns
    timeAv = timeTot / nRuns
    if nRuns==1:	# No error bars
       return SATFrac, timeAv
    timeVar = (time2Tot - timeTot*timeTot/nRuns)/(nRuns-1)
    timeSigmaMean = numpy.sqrt(timeVar/nRuns)
    # SAT2Tot = SATTot
    SATFracVar = (SATFrac - SATFrac*SATFrac/nRuns)/(nRuns-1)
    SATFracSigmaMean = numpy.sqrt(SATFracVar/nRuns)
    return SATFrac, SATFracSigmaMean, timeAv, timeSigmaMean, times


def PlotSATFrac(nVars, nRuns, k=3, 
		ClausesOverVars=numpy.arange(0.5, 8.0, 0.5),
		Solver=FasterDPSATSolver):
    """
    PlotSATFrac(nVars, nRuns, k=3, ClausesOverVars=arange(0.5, 8.0, 0.5),
    Solver=FasterDPSATSolver)
    collects data from SATFrac for a variety of nClauses, and 
    makes three plots versus the ratio clauses/variables: one with the mean
    and standard deviation of SATFrac, one of mean and sigma for times,
    and one of the distribution of times found at that ratio. 
    Use pylab.errorbar to draw the first two plots, and assemble two long
    lists of ratios and times for the last plot. Use pylab.figure(n) to 
    switch to plot #n.
    """
    SAT = []
    sig = []
    tSAT = []
    tsig = []
    ratiosTogether = []
    timesTogether = []
    for ratio in ClausesOverVars:
	nClauses = int(ratio*nVars)
	SATFracValue, SATSig, time, timeSig, times = \
		SATFrac(nVars, nClauses, nRuns, k, Solver=Solver)
	SAT.append(SATFracValue)
	sig.append(SATSig)
	tSAT.append(time)
	tsig.append(timeSig)
	ratiosTogether.extend([ratio]*nRuns)
	timesTogether.extend(times)
    pylab.figure(1)
    pylab.errorbar(ClausesOverVars, SAT, yerr=sig)
    pylab.figure(2)
    pylab.errorbar(ClausesOverVars, tSAT, yerr=tsig)
    pylab.figure(3)
    pylab.plot(ratiosTogether, timesTogether, 'r.')
    pylab.show()

def TimeTails(SATInstance, nRuns=10, showPlot=True):
    """TimeTails(SATInstance, nRuns=10, showPlot=True)
    computes and returns the times for nRuns runs of a SAT instance 
    shuffling the order of the variable list vars (and hence the order
    in which they are attempted by the DP algorithm) and resetting the 
    instance in between. GRAPHICS is turned off.
    Optionally, a pylab.hist plot of the histogram can be generated and shown.
    """
    global GRAPHICS
    GRAPHICS=False
    times = []
    sat_old = None
    for run in range(nRuns):
        SATInstance.Reset()
        sat, time = SATInstance.Solve()
        if sat != sat_old and sat_old is not None:
            print "Error: inconsistent answers"
        sat_old = sat
	times.append(time)
	SATInstance.ShuffleVariableOrder()
    if showPlot:
        pylab.hist(times)
    return times

def NTimeTails(nInst, nRuns, nVars=10, nClauses=40, k=3, plot=True,
	       Solver=FasterDPSATSolver):
    """NTimeTails(nInst, nRuns, nVars, nClauses, plot=True,
    Solver=FasterDPSATSolver) computes the runtimes for nRuns runs each
    of nInst SAT instances (k-SAT with nVars variables and nClauses clauses),
    returning a list of the runtime (number of flips) of each run. Optionally,
    a pylab plot of the histogram can be generated and shown.
    """
    pass

def LogLogHist(seq, nBins=10):
    """LogLogHist(seq, nBins=10) makes a log-log histogram (with nBins bins)
    of a sequence of numbers contained in seq.  Bin widths are sized
    exponentially (i.e., uniformly in logs) and bin counts are normalized
    by bin width."""
    pass



def TimeSlowFast(nInst, nRuns, nVars, nClauses):
    """TimeSlowFast(nInst, nRuns, nVars, nClauses) compares the runtimes
    of the Slow and Faster DPSATSolvers.
    """
    pass

def TestSAT(nVars, nClauses, k, with_graphics = True,
	    Solver=FasterDPSATSolver):
    """
    TestSAT(nVars, nClauses, k, with_graphics=True,
    Solver=FasterDPSATSolver) generates a single instance of k-SAT
    with nVars variables and nClauses clauses, and solves it using the
    specified solver class.  If with_graphics is True, then a graphics
    window showing progress of the algorithm is shown.  The status of
    the solution (i.e., True=SAT or False=UNSAT), along with the number
    of variable flips, is returned.
    """
    global GRAPHICS
    GRAPHICS = with_graphics
    clauses = MakeRandom_kSATClauses(nVars, nClauses, k)
    sol = Solver(clauses)
    SAT, time = sol.Solve()
    if SAT != sol.Verify():
        print "Error in SATFrac: solution doesn't check"
    else:
        if SAT:
            print "There is a solution."
        else:
            print "There is not a solution."

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
    print "3SAT solver"
    while 1:
        print "  Test kSAT solver: k=3, N=30, M=120"
        TestSAT(30, 120, 3, with_graphics=True)
        if not yesno(): return

if __name__=="__main__":
   demo()
    
