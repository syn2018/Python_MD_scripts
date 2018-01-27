"""
ScalingCollapse is a set of routines to implement rescalings of multiple
data sets according to a wide variety of functional forms.
"""

import numpy, pylab

def ScaleIt(xdata, ydata, xform="x -> x", yform="y -> y", 
            yerrdata={}, yerrform="yerr -> yerr", 
            keyNames=(), scalingParams={}, verbose=False):
    """
    xdata and ydata are two dictionaries with matched pairs of (x,y) arrays 
    [optionally labeled by some keys given (in order!) by keyNames, which 
    can be used in the rescaling functions and will be used in the legend].
    [Optional error bars for ydata can be given by yerrdata].
    xform and yform [and yerrform] are transformations or rescalings of the
    data, of the form 
      xform = "xDataVariable -> new ordinate", 
      yform = "yDataVariable -> new abscissa", 
      yerrform = "yErrorDataVariable -> new rescaled error".
    Any expression that can be evaluated using numpy, in terms of 
    the variables xDataVariable, yDataVariable and various other parameters 
    (name, value for which are items in scalingParams).
    Examples: 
      ScaleIt(xdata,ydata) will return them unchanged as x_sc, y_sc 
      ScaleIt(xdata,ydata,yerrdata=yerrdata) will return all three
      unchanged as x_sc, y_sc, yerr_sc
      ScaleIt(p_data,pathlength_data,yerrdata=pathlength_err_data,
              xform='p->p*Z*L', yform='ell->(2*ell*Z)/L',
              yerrform='sig->(2*sig*Z)/L', keyNames=('L','Z'),
      should scale x=p_data[L,Z] by Z*L, and y=path_length_data[L,Z] and the
      error bars by 2 Z/L.
      ScaleIt(S,D,yerrdata=dD, keyNames->('R'),
              xform='S->S**sigma * (R-Rc)', 
              yform='D->S**tau * D', yerrform='sig->S**tau * sig',
              scalingParams={'Rc':2.16,'tau':2.02})
      should return the scaling form S^tau D(S^sigma (R-Rc)) for cluster
      or avalanche sizes, for a series of curves D[R], S[R] giving the
      size distribution D(S,R).

    """
    x_sc = {}
    y_sc = {}
    yerr_sc = {}
    # incorporate numpy operations into the scalingParams dict so that
    #    numpy operations can appear in scaling forms
    scalingParams.update(numpy.__dict__)
    # xvar and yvar represent left-hand sides of transformations
    xvar = xform.split('->')[0].strip()
    yvar = yform.split('->')[0].strip()
    # xfunc and yfunc represent right-hand sides of transformations
    xfunc = xform.split('->')[1].strip()
    yfunc = yform.split('->')[1].strip()
    if yerrform:
        yerrvar = yerrform.split('->')[0].strip()
        yerrfunc = yerrform.split('->')[1].strip()
    for cP, x in xdata.items():
        #convert data to array (in case not already so)
        x = numpy.array(x)
        try:
            y = numpy.array(ydata[cP])
        except:
            raise KeyError, "parameters for xdata and ydata do not match"
        if len(x) != len(y):
            raise KeyError, \
	       "Length of x and y arrays for key %s do not agree"%str(cP)
        if keyNames != ():
            curveParams = dict(zip(keyNames, cP))
        else:
            curveParams = {}
        curveParams[xvar] = x
        curveParams[yvar] = y
        x_sc[cP] = eval(xfunc, scalingParams, curveParams)
        y_sc[cP] = eval(yfunc, scalingParams, curveParams)
        if cP in yerrdata:
            curveParams[yerrvar] = numpy.array(yerrdata[cP])
            yerr_sc[cP] = eval(yerrfunc, scalingParams, curveParams)
        if verbose:
	    if cP != ():
                for name, val in zip(keyNames, cP):
                    print '%s = %s' % (name, val)
            print '%s = %s' % (xvar, x)
            print '%s_sc = %s' % (xvar, x_sc[cP])
            print '%s = %s' % (yvar, y)
            print '%s_sc = %s' % (yvar, y_sc[cP])
            if yerrdata:
                print '%s = %s' % (yerrvar, yerr)
                print '%s_sc = %s' % (yerrvar, yerr_sc[cP])
    return x_sc, y_sc, yerr_sc

        
def ColorWheel():
    """Cycles through a selection of colors, symbols, and line styles
    for matlibplot.matlab.plot."""
    colors = ['b', 'g', 'r', 'c', 'm', 'k']
    symbols = ['o', 's', '^', 'v', '<', ">", 'x', 'D', 'h', 'p']
    lines = ['-', '--', '-.', ':']
    while 1:
        for l in lines:
           for s in symbols:
                for c in colors:
                   yield c + s + l

def MultiPlot(xdata, ydata, xlabel=None, ylabel=None,
              xform="x -> x", yform="y -> y", 
              yerrdata={}, yerrform="yerr -> yerr", 
              keyNames=(), scalingParams={},
              fmt={}, verbose=False,
              loc='upper right', showIt=True, log=None):
    """
    MultiPlot will plot, and optionally rescale,
    a series of x,y curves, with legends and symbols and different
    colors and line styles. 
    xdata and ydata are two dictionaries with matched pairs of (x,y) arrays,
    keyed on identifiers for the curves.
    Examples:
       xdata['curve1'] = [x0, x1, ...]
       ydata['curve1'] = [y0, y1, ...]
       ...
       and optionally,
       yerrdata['curve1'] = [yerr0, yerr1, ...]

       will generate the legend 'curve1' for the data described by
       (xdata['curve1'], ydata['curve1'])

       xdata[10, 4] = [x0, x1, ...]
       ydata[10, 4] = [y0, y1, ...]
       ...
       keyNames = ('L', 'Z')

       will generate the legend 'L=10, Z=4' for the data described by
       (xdata[10, 4], ydata[10, 4])

       xdata[2.16, 100] = [s0, s1, ...]
       ydata[2.16, 100] = [d0, d1, ...]
       ...
       keyNames = ('R', 'L')
       xform = 'S -> (S**sigma)*(R-Rc)'
       yform = 'D -> (S**tau)*D'
       scalingParameters = {'sigma': 0.3, 'tau': 1.5, 'Rc': 2.16}

       will scale the xdata S by the transformation in xform and the
       ydata D by the transformation in yform, using the value of the
       scaling parameters sigma, tau, and Rc as specified in the
       scalingParameters dictionary

    See ScaleIt for usage of most of the variables.
    See help(pylab.legend) for possible legend locations loc.
    Argument log='x', 'y', or 'xy' will generate log scales on the
    corresponding axes.
    """
    x_sc, y_sc, yerr_sc = ScaleIt(xdata, ydata, xform, yform, 
                                  yerrdata, yerrform,
                                  keyNames, scalingParams, verbose)
    xfunc = xform.split('->')[1].strip()
    yfunc = yform.split('->')[1].strip()
    # Could pull this out, to allow multiple curve sets
    colorWheel = ColorWheel() 
    lines = []
    labels = []
    for cP, x_scaled in x_sc.items():
        if cP in fmt:
            format = fmt[cP]
        else:
            format = colorWheel.next()
        label=""
        for i, v in enumerate(cP):
            if keyNames != ():
                label += keyNames[i] + "=" + str(v) + " "
            else:
                label += str(v) + " "
        if cP in yerrdata:
            line, errorbars = pylab.errorbar(x_scaled, y_sc[cP], 
                                               yerr=yerr_sc[cP], 
                                               fmt=format)
        else: 
            line = pylab.plot(x_scaled, y_sc[cP], format)
        lines.append(line)
        labels.append(label)
    labels_and_lines = zip(labels, lines)
    labels_and_lines.sort()
    labels = [ll[0] for ll in labels_and_lines]
    lines =  [ll[1] for ll in labels_and_lines]
    pylab.legend(lines, labels, loc)
    if xlabel:
        pylab.xlabel(xlabel)
    else:
        pylab.xlabel(xfunc)
    if ylabel:
        pylab.ylabel(ylabel)
    else:
        pylab.ylabel(yfunc)
    # Wait on pylab.show() by default so external program can set options
    if log == 'x':
        pylab.semilogx()
    elif log == 'y':
        pylab.semilogy()
    elif log == 'xy':
        pylab.loglog()
    if showIt: pylab.show()
    return x_sc, y_sc, yerr_sc
        
def testScalingCollapse():
    """Simple test of MultiPlot and ScaleIt"""
    p = {}
    ell = {}
    sig = {}
    p[(100, 2)] = numpy.array([4.,8.,12.])
    ell[(100, 2)] = numpy.array([10.,20.,30.])
    sig[(100, 2)] = numpy.array([0.5,0.5,0.5])
    p[(200, 4)] = numpy.array([1.,2.,3.])
    ell[(200, 4)] = numpy.array([11.,18.,33.])
    sig[(200, 4)] = numpy.array([1.0,2.0,3.0])
    pylab.figure(1)
    MultiPlot(p, ell,yerrdata=sig)
    pylab.figure(2)
    sig = {}
    sig[(100, 2)] = numpy.array([0.5,0.5,0.5])
    sig[(200, 4)] = numpy.array([1.0,2.0,3.0])
    p_sc, ell_sc, ellerr_sc = MultiPlot(p, ell, yerrdata=sig, 
            xform='p->p*Z*L', yform='ell->(2*ell*Z)/L', 
	    yerrform = 'sig->(2*sig*Z)/L', keyNames=('L', 'Z'), 
	    showIt=False, loc=3)
    pylab.semilogx()
    pylab.show()
    return p, ell


