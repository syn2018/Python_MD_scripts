#
# See the exercise "SandP.pdf" from SandP.html
# in  http://www.physics.cornell.edu/sethna/StatMech/ComputerExercises/
#
# Import packages
# Old or new version of plotter?
import pylab
import scipy

# Read in the file of Standard and Poor's average versus time
t = []
SP = []
for line in file("SandPConstantDollars.dat"):
    day, sandp = map(float, line.split())   
    t.append(day)
    SP.append(sandp)

# Convert from "list" to "array" form (arrays can be multiplied and added, etc.)

t = scipy.array(t)
SP = scipy.array(SP)

# ***** Plot SP versus t
# ***** pylab.plot(t, SP)
# ***** pylab.show()
# ***** Note 9/11/01 is day 6903: is it the cause for the post-2000 drop?"
# ***** Zoom in and see: did the terrorist attach on the World Trade Center
# ***** trigger the stock market downturn, or was it just a small extra dip in
# ***** an overall pattern?

def P(lag):
   """
   Function which returns a list of percentage changes after a number
   "lag" of trading days. 
   P(1) gives the daily percentage changes, P(5) the weekly changes, etc.
   #
   Arrays add and subtract and multiply by scalars just like vectors
   SP[m,n] is the part of the array starting at m and ending at n
   Arrays and lists in python start at zero and end at len(list)-1.
   #
   So, the easy way to compute the vector of fractional changes after 
   a time "lag" is 
	ratios = SP[lag:N]/SP[0:N-lag]
   where N is len(SP)
   and the list of percentage changes is 
	P = 100.*(ratios-1.)
   You'll need also to
	return P
   """
   pass

def PlotPHistogram(lag):
   """
   Plotting a histogram using pylab is easy: make a list or array
   of data, and call
	(n,bins,patches) = pylab.hist(data,bins=100,normed=True)
   for 100 equal-sized bins, with the values rescaled so the area is one.
   You may need to do
	pylab.show()
   to get the graph to come up.
   #
   Usage:
       PlotPHistogram(1)  	# Day
       PlotPHistogram(5)  	# Week
       PlotPHistogram(252)	# Year
   """
   pass

def PlotLogPHistogram(lag):
   """
   The histogram function can be run without makeing a plot.
	(n,bins) = pylab.mlab.hist(P(lag), bins=100, normed=True)
   It returns the normalized counts in "n", the start of each bin
   in "bins", and some rectangle objects in "patches". To get the
   bin centers, you can add (bins[1]-bins[0])/2
   to each of the bins (adding a scalar to an array works element by
   element). Also you can take the log of an array as a whole, so you
   want to use pylab.plot to plot scipy.log(n+0.0001) versus binCenters. 
   (We need to add a small constant to avoid log(0) for empty bins.)
   #
   Usage:
	PlotLogPHistogram(5)	# Week
   """
   pass

# PlotLogPHistogram(5)

def V(lag):
   """
   Volatility: standard deviation of percentage change after time lag
   scipy.sum(A) adds the numbers in an array A
   A**2 is a new array whose entries are the entries of A squared
   len(A) gives the number of entries
   Find the mean percentage change after lag
   Find the variance var = (p - mean)**2/(len(p)-1) 
   [Note: N-1 in denom from stats]
   Return the volatility = scipy.sqrt(var)
   #
   Usage:
	lags = range(100)
	volatilities = [V[lag] for lag in lags]
	pylab.plot(...) 
   """
   pass

