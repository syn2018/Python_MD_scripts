#
# See the exercise "SandP.pdf" from SandP.html
# in  http://www.physics.cornell.edu/sethna/StatMech/ComputerExercises/
#
import pylab
import scipy
import RandomArray

CPITime = []
CPI = []
for line in file("CPIMonthly.dat"):
    t, cpi = map(float, line.split())   
    CPITime.append(t)
    CPI.append(cpi)

SandPTime = []
SandP = []
for line in file("SandPXY.dat"):
    t, sandp = map(float, line.split())   
    SandPTime.append(t)
    SandP.append(sandp)

SandPTime.reverse()
SandP.reverse()

SandPConstantDollar = SandP[:]
SandPTime = scipy.array(SandPTime)
SandPConstantDollar = scipy.array(SandPConstantDollar)

index = 0;

for SandPIndex, t in enumerate(SandPTime):
    while(t > CPITime[index+1]):
       index+=1
    tFraction = (t-CPITime[index]) / (CPITime[index+1]-CPITime[index])
    cpiInterpolated = CPI[index] + tFraction * (CPI[index+1]-CPI[index])
    SandPConstantDollar[SandPIndex] /= (cpiInterpolated/100.)

def RemoveMeanGrowth(traj, time):
    eta = scipy.log(traj[-1]/traj[0])/(time[-1]-time[0])
    trajRenormalized = traj / (traj[0]*scipy.exp(eta * (time-time[0])))
    return eta, trajRenormalized

def SandPConstantDollarOutput():
    outputSandP = file("SandPConstantDollar.dat", "w")
    dayLength = SandPTime[1]-SandPTime[0]
    for t, data in zip(SandPTime, SandPConstantDollar):
        days = int(scipy.round((t-SandPTime[0])/dayLength))
        outputSandP.write("%s %s\n"%(days, data))
    outputSandP.close()
    

# eta is 6.424% for SandP until November 2004

def PlotFit(traj=SandPConstantDollar, times=SandPTime, nRandomTrajs = 0,
            a=0.04, output=False):

    plotEm = []
    eta, trajFlattened = RemoveMeanGrowth(traj, times)
    plotEm.append(times)
    plotEm.append(trajFlattened)

    randomTimes = scipy.arange(times[0], times[-1]+1.e-10, (times[-1]-times[0])/(len(times)-1))
    randomTrajs = []
    for i in range(nRandomTrajs):
        randomTrajs.append(
	    scipy.exp(scipy.cumsum(a * (RandomArray.random(len(times))-0.5))))
    randomTrajsFlattened = randomTrajs[:]
    for rF in randomTrajsFlattened:
      eta, rF = RemoveMeanGrowth(rF, randomTimes)
      plotEm.append(randomTimes)
      plotEm.append(rF)
    
    pylab.plot(*plotEm)
    pylab.show()

    if output:
        outputSandP = file("SandPFlattened.dat", "w")
	for t, data in zip(times, trajFlattened):
            outputSandP.write("%s %s\n"%(t, data))
        outputSandP.close()
        outputRandom = file("OneDRandomFlattened.dat", "w")
	for rF in randomTrajs:
            eta, rF = RemoveMeanGrowth(rF, randomTimes)
	    for t, data in zip(randomTimes, rF):
                outputRandom.write("%s %s\n"%(t, data))
            outputRandom.write("\n")
        outputRandom.close()

