import os

# No demos yet for NumberPartitioning.py, Pendulum.py, Primes.py, 
# RandText.py, SandPCPI.py, Walker.py

demos = ["ChaosLyapunov.py",
	 "Equilibration.py",
	 "ExponentialAtmosphereA.py",
	 "ExponentialAtmosphereE.py",
	 "FitzNag2D.py",
	 "FitzNag.py",
	 "FractalDimension.py",
	 "Gumbel.py",
	 "InvariantMeasure.py",
	 "Ising.py",
	 "IterateLogistic.py",
	 "kSAT.py",
	 "PairDistributionMD.py",
	 "PendulumPeriod.py",
	 "Percolation.py",
	 "PerfumeWalk.py",
	 "PeriodDoubling.py",
	 "Pressure.py",
	 "Primes.py",
	 "RandomMatrixTheory.py",
	 "RandomWalk.py",
	 "Repressilator.py",
	 "Repressilator_v0.py",
	 "SandP.py",
	 "SmallWorldNetworks.py",
	 "StochasticCells.py",
	 "WalkerAnimate.py",
	 "WalkerAnimate_v0.py",
	 "WalkerPhasePlot.py",
	 "WalkerUnstableOrbit.py"
	 ]


for demo in demos:
    response = raw_input(
              "Next demo: %s; run (type y), skip (n), or end (e)?"%demo)
    if len(response)==0:	# [CR] continues too
        os.system("python %s"%demo)
    elif (response[0]=="e" or response[0]=="E"):
            break
    elif (response[0]!="s" and response[0]!="S" 
    	   and response[0]!="n" and response[0]!="N"): # Default is continue
        os.system("python %s"%demo)
