import numpy

def IntegerHistogram(integers, bins):
    #
    # Histogram of numbers of integers in ranges bins[x] <= i < bins[x+1]
    # Returns two arrays:
    # centers of integer range ic
    # and normalized probability distribution D[ic]
    # 
    # Find out how many entries are in each bin 
    # (from scipy.stats.histogram2)
    n = numpy.searchsorted(numpy.sort(integers),bins)
    n = numpy.concatenate([n, [len(integers)]])
    Number = n[1:] - n[:-1]
    NumberOverflow = Number[-1]
    Number = Number[:-1]

    # Then, find out min and max integer in each range
    #
    # Inefficient way if range is enormous
    # intPossible = numpy.arange(bins[len(bins)-1])
    # n = numpy.searchsorted(numpy.sort(intPossible),bins)
    # n = numpy.concatenate([n, [len(intPossible)]])
    # iNumber = n[1:] - n[:-1]

    # Efficient method: gives array one shorter (missing zero at end)
    #iNumber = numpy.zeros(len(bins)-1, numpy.Float)
    #centers = numpy.zeros(len(bins)-1, numpy.Float)
    iNumber = numpy.zeros(len(bins)-1, float)
    centers = numpy.zeros(len(bins)-1, float)
    for i in range(len(bins)-1):
        if (numpy.floor(bins[i+1])==bins[i+1]):
	    iNumber[i] = numpy.floor(bins[i+1])-numpy.ceil(bins[i])
	else:
	    iNumber[i] = numpy.floor(bins[i+1])-numpy.ceil(bins[i])+1
	centers[i] = numpy.ceil(bins[i]) + 0.5*(iNumber[i]-1)

    # Eliminate 0/0 if bins have no integers
    # Is there a fancy method for this?
    outCenter = []
    outD = []
    for i in range(len(centers)):
        if iNumber[i] != 0:
	    outCenter.append(centers[i])
	    outD.append(Number[i]/iNumber[i])
	else:
	    assert(outD == [], 
	    	"integer_histogram has data in bin where no ints are!")

    return numpy.array(outCenter), numpy.array(outD)
