# -*- coding: utf-8 -*-

"""
Small collection of robust statistical estimators based on functions from
Henry Freudenriech (Hughes STX) statistics library (called ROBLIB) that have
been incorporated into the AstroIDL User's Library.  Function included are:
  * biweightMean - biweighted mean estimator
  * mean - robust estimator of the mean of a data set
  * std - robust estimator of the standard deviation of a data set
  * checkfit - return the standard deviation and biweights for a fit in order 
    to determine its quality
  * linefit - outlier resistant fit of a line to data
  * polyfit - outlier resistant fit of a polynomial to data
  * reject_outliers - reject outlier value from input data (jmiguel@iaa.es)
  * r_division - scalar robust division
  * r_divisionN - Array robust division

For the fitting routines, the coefficients are returned in the same order as
numpy.polyfit, i.e., with the coefficient of the highest power listed first.

For additional information about the original IDL routines, see:
  http://idlastro.gsfc.nasa.gov/contents.html#C17
"""

from __future__ import division

import math
import numpy

__version__ = '0.4'
__revision__ = '$Rev$'
__all__ = ['biweightMean', 'mean', 'std', 'checkfit', 'linefit', 'polyfit', 
           '__version__', '__revision__', '__all__']

__iterMax = 25
__delta = 5.0e-7
__epsilon = 1.0e-20

def reject_outliers(data, m=2):
    """
    Reject outliers values from input data array
    """
    return data[abs(data - numpy.mean(data)) < m * numpy.std(data)]

def r_division(num, den):
    """
    Scalar robust division to avoid zero-division indefinition
    """
    if math.fabs(den)<__epsilon:
        return num
    else:
        return num/den

def r_divisionN(num, den):
    """
    Array robust division to avoid zero-division indefinition
    """
    den = numpy.where(numpy.fabs(den)<__epsilon, 1.0, den)
    return num/den
    
    
def biweightMean(inputData):
    """
    Calculate the mean of a data set using bisquare weighting.  
    
    Based on the biweight_mean routine from the AstroIDL User's 
    Library.
    """
    
    y = inputData.ravel()
    if type(y).__name__ == "MaskedArray":
        y = y.compressed()
    
    n = len(y)
    closeEnough = 0.03*numpy.sqrt(0.5/(n-1))
    
    diff = 1.0e30
    nIter = 0
    
    y0 = numpy.median(y)
    deviation = y - y0
    sigma = std(deviation)
    
    if sigma < __epsilon:
        diff = 0
    while diff > closeEnough:
        nIter = nIter + 1
        if nIter > __iterMax:
            break
        uu = ((y-y0)/(6.0*sigma))**2.0
        uu = numpy.where(uu > 1.0, 1.0, uu)
        weights = (1.0-uu)**2.0
        weights /= weights.sum()
        y0 = (weights*y).sum()
        deviation = y - y0
        prevSigma = sigma
        sigma = std(deviation, Zero=True)
        if sigma > __epsilon:
            diff = numpy.abs(prevSigma - sigma) / prevSigma
        else:
            diff = 0.0
            
    return y0


def r_nanmean(inputData, axis=None):
    """
    Compute the arithmetic mean of a data set ignoring the NaNs values
    and checking what numpy and scipy versions are installed.
    
    Note: Because numpy.nanmean() was implemented in v1.8.0,
    and scipy.nanmean() is deprecated in scipy 0.15.0 in favour of 
    numpy.nanmean(), we check the version and then decide the function
    to use.

    Note2: PAPI@panic22.caha.es has numpy 1.7.x and scipy 0.12.x
 
 
    Parameters:
    -----------
    inputDat : array_like
        Input array or object that can be converted to an array.
    axis : int, optional
        Axis along which the medians are computed. The default (axis=None) is 
        to compute the median along a flattened version of the array.
        axis: axis=None is to compute the median along a flattened version of 
        the array.
        
        
    """ 

    try:
        from numpy import nanmean as nanmean
        # fastest option
        return nanmean(inputData, axis=axis)
    except ImportError:
        try:
            from scipy.stats import nanmean as nanmean
            return nanmean(inputData, axis=axis)    
        except ImportError:
            # no built-it function, then compute mean masking NaNs
            # slower option: numpy.ma.masked_invalid(inputData).median()
            return numpy.mean(inputData[~numpy.isnan(inputData)], axis=axis)
          
def r_nanmedian(inputData, axis=None):
    """
    Compute the arithmetic median of a data set ignoring the NaNs values
    and checking what numpy and scipy versions are installed.
    
    Note: Because numpy.nanmedian() was implemented in v1.9.0,
    and scipy.nanmedian() is deprecated in scipy 0.15.0 in favour of 
    numpy.nanmedian(), we check the version and then decide the function
    to use.

    Note2: PAPI@panic22.caha.es has numpy 1.7.x and scipy 0.12.x
 
    Parameters:
    -----------
    inputDat : array_like
        Input array or object that can be converted to an array.
    axis : int, optional
        Axis along which the medians are computed. The default (axis=None) is 
        to compute the median along a flattened version of the array.
        axis: axis=None is to compute the median along a flattened version of 
        the array.
    
 
    """     
    
    try:
        from numpy import nanmedian as nanmedian
        # fastest option
        return nanmedian(inputData, axis=axis)
    except ImportError:
        try:
            from scipy.stats import nanmedian as nanmedian
            return nanmedian(inputData, axis=axis)    
        except ImportError:
            # no built-it function, then compute median masking NaNs
            # slower option: numpy.ma.masked_invalid(inputData).median()
            return numpy.median(inputData[~numpy.isnan(inputData)], axis=axis)


def mean(inputData, Cut=3.0):
    """
    Robust estimator of the mean of a data set.  Based on the 
    resistant_mean function from the AstroIDL User's Library.

    .. seealso::
        :func:`lsl.misc.mathutil.robustmean`
        
    2016-Feb-4: modified for use of r_nanmean
    
    """

    data = inputData.ravel()
    if type(data).__name__ == "MaskedArray":
        data = data.compressed()

    data0 = r_nanmedian(data)
    
    maxAbsDev = r_nanmedian(numpy.abs(data-data0)) / 0.6745
    if maxAbsDev < __epsilon:
        maxAbsDev = r_nanmean(numpy.abs(data-data0)) / 0.8000

    cutOff = Cut*maxAbsDev
    good = numpy.where( numpy.abs(data-data0) <= cutOff )
    good = good[0]
    dataMean = r_nanmean(data[good])
    dataSigma = math.sqrt( ((data[good] - dataMean)**2.0).sum() / len(good) )

    if Cut > 1.0:
        sigmaCut = Cut
    else:
        sigmaCut = 1.0
    if sigmaCut <= 4.5:
        dataSigma = dataSigma / (-0.15405 + 0.90723*sigmaCut - 0.23584*sigmaCut**2.0 + 0.020142*sigmaCut**3.0)

    cutOff = Cut*dataSigma
    good = numpy.where(  numpy.abs(data-data0) <= cutOff )
    good = good[0]
    dataMean = r_nanmean(data[good])
    if len(good) > 3:
        dataSigma = math.sqrt( ((data[good] - dataMean)**2.0).sum() / len(good) )

    if Cut > 1.0:
        sigmaCut = Cut
    else:
        sigmaCut = 1.0
    if sigmaCut <= 4.5:
        dataSigma = dataSigma / (-0.15405 + 0.90723*sigmaCut - 0.23584*sigmaCut**2.0 + 0.020142*sigmaCut**3.0)

    dataSigma = dataSigma / math.sqrt(len(good)-1)

    return dataMean


def std(inputData, Zero=False):
    """
    Robust estimator of the standard deviation of a data set.  
    
    Based on the robust_sigma function from the AstroIDL User's Library.
    
    2016-Feb-4: modified for use of r_nanmean
    
    """

    data = inputData.ravel()
    if type(data).__name__ == "MaskedArray":
        data = data.compressed()

    if Zero:
        data0 = 0.0
    else:
        data0 = r_nanmedian(data)
    maxAbsDev = r_nanmedian(numpy.abs(data-data0)) / 0.6745
    if maxAbsDev < __epsilon:
        maxAbsDev = r_nanmean(numpy.abs(data-data0)) / 0.8000
    if maxAbsDev < __epsilon:
        sigma = 0.0
        return sigma

    u = (data - data0) / 6.0 / maxAbsDev
    u2 = u**2.0
    good = numpy.where( u2 <= 1.0 )
    good = good[0]
    if len(good) < 3:
        print "WARNING:  Distribution is too strange to compute standard deviation"
        sigma = -1.0
        return sigma

    numerator = ((data[good]-data0)**2.0 * (1.0-u2[good])**2.0).sum()
    nElements = (data.ravel()).shape[0]
    denominator = ((1.0-u2[good])*(1.0-5.0*u2[good])).sum()
    sigma = nElements*numerator / (denominator*(denominator-1.0))
    if sigma > 0:
        sigma = math.sqrt(sigma)
    else:
        sigma = 0.0

    return sigma


def checkfit(inputData, inputFit, epsilon, delta, BisquareLimit=6.0):
    """
    Determine the quality of a fit and biweights.  Returns a tuple
    with elements:
      0. Robust standard deviation analog
      1. Fractional median absolute deviation of the residuals
      2. Number of input points given non-zero weight in the calculation
      3. Bisquare weights of the input points
      4. Residual values scaled by sigma
    
    This function is based on the rob_checkfit routine from the AstroIDL 
    User's Library.
    """
    
    data = inputData.ravel()
    fit = inputFit.ravel()
    if type(data).__name__ == "MaskedArray":
        data = data.compressed()
    if type(fit).__name__ == "MaskedArray":
        fit = fit.compressed()

    deviation = data - fit
    sigma = std(deviation, Zero=True)
    if sigma < epsilon:
        return (sigma, 0.0, 0, 0.0, 0.0)
    
    toUse = (numpy.where( numpy.abs(fit) > epsilon ))[0]
    if len(toUse) > 3:
        fracDev = 0.0
    else:
        fracDev = numpy.median(numpy.abs(deviation[toUse]/fit[toUse]))
    if fracDev < delta:
        return (sigma, fracDev, 0, 0.0, 0.0)
        
    biweights = numpy.abs(deviation)/(BisquareLimit*sigma)
    toUse = (numpy.where(biweights > 1))[0]
    if len(toUse) > 0:
        biweights[toUse] = 1.0
    nGood = len(data) - len(toUse)
    
    scaledResids = (1.0 - biweights**2.0)
    scaledResids = scaledResids / scaledResids.sum()
    
    return (sigma, fracDev, nGood, biweights, scaledResids)


def linefit(inputX, inputY, iterMax=25, Bisector=False, BisquareLimit=6.0, CloseFactor=0.03):
    """
    Outlier resistance two-variable linear regression function.
    
    Based on the robust_linefit routine in the AstroIDL User's Library.
    """
    
    xIn = inputX.ravel()
    yIn = inputY.ravel()
    if type(yIn).__name__ == "MaskedArray":
        xIn = xIn.compress(numpy.logical_not(yIn.mask))
        yIn = yIn.compressed()
    n = len(xIn)
    
    x0 = xIn.sum() / n
    y0 = yIn.sum() / n
    x = xIn - x0
    y = yIn - y0
    
    cc = numpy.zeros(2)
    ss = numpy.zeros(2)
    sigma = 0.0
    yFit = yIn
    badFit = 0
    nGood = n
    
    lsq = 0.0
    yp = y
    if n > 5:
        s = numpy.argsort(x)
        u = x[s]
        v = y[s]
        nHalf = n/2 -1
        x1 = numpy.median(u[0:nHalf])
        x2 = numpy.median(u[nHalf:])
        y1 = numpy.median(v[0:nHalf])
        y2 = numpy.median(v[nHalf:])
        if numpy.abs(x2-x1) < __epsilon:
            x1 = u[0]
            x2 = u[-1]
            y1 = v[0]
            y2 = v[-1]
        cc[1] = (y2-y1)/(x2-x1)
        cc[0] = y1 - cc[1]*x1
        yFit = cc[0] + cc[1]*x
        sigma, fracDev, nGood, biweights, scaledResids = checkfit(yp, yFit, __epsilon, __delta)
        if nGood < 2:
            lsq = 1.0
        
    if lsq == 1 or n < 6:
        sx = x.sum()
        sy = y.sum()
        sxy = (x*y).sum()
        sxx = (x*x).sum()
        d = sxx - sx*sx
        if numpy.abs(d) < __epsilon:
            return (0.0, 0.0)
        ySlope = (sxy - sx*sy) / d
        yYInt = (sxx*sy - sx*sxy) / d
        
        if Bisector:
            syy = (y*y).sum()
            d = syy - sy*sy
            if numpy.abs(d) < __epsilon:
                return (0.0, 0.0)
            tSlope = (sxy - sy*sx) / d
            tYInt = (syy*sx - sy*sxy) / d
            if numpy.abs(tSlope) < __epsilon:
                return (0.0, 0.0)
            xSlope = 1.0/tSlope
            xYInt = -tYInt / tSlope
            if ySlope > xSlope:
                a1 = yYInt
                b1 = ySlope
                r1 = numpy.sqrt(1.0+ySlope**2.0)
                a2 = xYInt
                b2 = xSlope
                r2 = numpy.sqrt(1.0+xSlope**2.0)
            else:
                a2 = yYInt
                b2 = ySlope
                r2 = numpy.sqrt(1.0+ySlope**2.0)
                a1 = xYInt
                b1 = xSlope
                r1 = numpy.sqrt(1.0+xSlope**2.0)
            yInt = (r1*a2 + r2*a1) / (r1 + r2)
            slope = (r1*b2 + r2*b1) / (r1 + r2)
            r = numpy.sqrt(1.0+slope**2.0)
            if yInt > 0:
                r = -r
            u1 = slope / r
            u2 = -1.0/r
            u3 = yInt / r
            yp = u1*x + u2*y + u3
            yFit = y*0.0
            ss = yp
        else:
            slope = ySlope
            yInt = yYInt
            yFit = yInt + slope*x
        cc[0] = yInt
        cc[1] = slope
        sigma, fracDev, nGood, biweights, scaledResids = checkfit(yp, yFit, __epsilon, __delta)
        
    if nGood < 2:
        cc[0] = cc[0] + y0 - cc[1]*x0
        return cc[::-1]
        
    sigma1 = (100.0*sigma)
    closeEnough = CloseFactor * numpy.sqrt(0.5/(n-1))
    if closeEnough < __delta:
        closeEnough = __delta
    diff = 1.0e20
    nIter = 0
    while diff > closeEnough:
        nIter = nIter + 1
        if nIter > iterMax:
            break
        sigma2 = sigma1
        sigma1 = sigma
        sx = (biweights*x).sum()
        sy = (biweights*y).sum()
        sxy = (biweights*x*y).sum()
        sxx = (biweights*x*x).sum()
        d = sxx - sx*sx
        if numpy.abs(d) < __epsilon:
            return (0.0, 0.0)
        ySlope = (sxy - sx*sy) / d
        yYInt = (sxx*sy - sx*sxy) / d
        slope = ySlope
        yInt = yYInt
        
        if Bisector:
            syy = (biweights*y*y).sum()
            d = syy - sy*sy
            if numpy.abs(d) < __epsilon:
                return (0.0, 0.0)
            tSlope = (sxy - sy*sx) / d
            tYInt = (syy*sx - sy*sxy) / d
            if numpy.abs(tSlope) < __epsilon:
                return (0.0, 0.0)
            xSlope = 1.0/tSlope
            xYInt = -tYInt / tSlope
            if ySlope > xSlope:
                a1 = yYInt
                b1 = ySlope
                r1 = numpy.sqrt(1.0+ySlope**2.0)
                a2 = xYInt
                b2 = xSlope
                r2 = numpy.sqrt(1.0+xSlope**2.0)
            else:
                a2 = yYInt
                b2 = ySlope
                r2 = numpy.sqrt(1.0+ySlope**2.0)
                a1 = xYInt
                b1 = xSlope
                r1 = numpy.sqrt(1.0+xSlope**2.0)
            yInt = (r1*a2 + r2*a1) / (r1 + r2)
            slope = (r1*b2 + r2*b1) / (r1 + r2)
            r = numpy.sqrt(1.0+slope**2.0)
            if yInt > 0:
                r = -r
            u1 = slope / r
            u2 = -1.0/r
            u3 = yInt / r
            yp = u1*x + u2*y + u3
            yFit = y*0.0
            ss = yp
        else:
            yFit = yInt + slope*x
        cc[0] = yInt
        cc[1] = slope
        sigma, fracDev, nGood, biweights, scaledResids = checkfit(yp, yFit, __epsilon, __delta)
        
        if nGood < 2:
            badFit = 1
            break
        diff1 = numpy.abs(sigma1 - sigma)/sigma
        diff2 = numpy.abs(sigma2 - sigma)/sigma
        if diff1 < diff2:
            diff = diff1
        else:
            diff = diff2
                
    cc[0] = cc[0] + y0 - cc[1]*x0
    return cc[::-1]


def polyfit(inputX, inputY, order, iterMax=25):
    """
    Outlier resistance two-variable polynomial function fitter.
    
    Based on the robust_poly_fit routine in the AstroIDL User's 
    Library.
    
    Unlike robust_poly_fit, two different polynomial fitters are used
    because numpy.polyfit does not support non-uniform weighting of the
    data.  For the weighted fitting, the SciPy Orthogonal Distance
    Regression module (scipy.odr) is used.
    """
    
    from scipy import odr
    
    def polyFunc(B, x, order=order):
        out = x*0.0
        for i in range(order+1):
            out = out + B[i]*x**i
    
    model = odr.Model(polyFunc)
    
    x = inputX.ravel()
    y = inputY.ravel()
    if type(y).__name__ == "MaskedArray":
        x = x.compress(numpy.logical_not(y.mask))
        y = y.compressed()
    n = len(x)
    
    x0 = x.sum() / n
    y0 = y.sum() / n
    u = x
    v = y
        
    nSeg = order + 2
    if (nSeg/2)*2 == nSeg:
        nSeg = nSeg + 1
    minPts = nSeg*3
    if n < 1000:
        lsqFit = 1
        cc = numpy.polyfit(u, v, order)
        yFit = numpy.polyval(cc, u)
    else:
        lsqfit = 0
        q = numpy.argsort(u)
        u = u[q]
        v = v[q]
        nPerSeg = numpy.zeros(nSeg) + n/nSeg
        nLeft = n - nPerSeg[0]*nSeg
        nPerSeg[nSeg/2] = nPerSeg[nSeg/2] + nLeft
        r = numpy.zeros(nSeg)
        s = numpy.zeros(nSeg)
        r[0] = numpy.median(u[0:nPerSeg[0]])
        s[0] = numpy.median(v[0:nPerSeg[0]])
        i2 = nPerSeg[0]-1
        for i in range(1,nSeg):
            i1 = i2
            i2 = i1 + nPerSeg[i]
            r[i] = numpy.median(u[i1:i2])
            s[i] = numpy.median(v[i1:i2])
        cc = numpy.polyfit(r, s, order)
        yFit = numpy.polyval(cc, u)
        
    sigma, fracDev, nGood, biweights, scaledResids = checkfit(v, yFit, __epsilon, __delta)
    if nGood == 0:
        return cc
    if nGood < minPts:
        if lsqFit == 0:
            cc = numpy.polyfit(u, v, order)
            yFit = numpy.polyval(cc, u)
            sigma, fracDev, nGood, biweights, scaledResids = checkfit(v, yFit, __epsilon, __delta)
            if nGood == 0:
                return __processPoly(x0, y0, order, cc)
            nGood = n - nGood
        if nGood < minPts:
            return 0
        
    closeEnough = 0.03*numpy.sqrt(0.5/(n-1))
    if closeEnough < __delta:
        closeEnough = __delta
    diff = 1.0e10
    sigma1 = 100.0*sigma
    nIter = 0
    while diff > closeEnough:
        nIter = nIter + 1
        if nIter > iterMax:
            break
        sigma2 = sigma1
        sigma1 = sigma
        g = (numpy.where(biweights > 0))[0]
        if len(g) < len(biweights):
            u = u[g]
            v = v[g]
            biweights = biweights[g]
        data = odr.RealData(u, v, sy=1.0/biweights)
        fit = odr.ODR(data, model, beta0=cc[::-1])
        out = fit.run()
        cc = out.beta[::-1]
        yFit = numpy.polyval(cc, u)
        sigma, fracDev, nGood, biweights, scaledResids = checkfit(v, yFit, __epsilon, __delta)
        if nGood < minPts:
            return cc
        diff1 = numpy.abs(sigma1 - sigma)/sigma
        diff2 = numpy.abs(sigma2 - sigma)/sigma
        if diff1 < diff2:
            diff = diff1
        else:
            diff = diff2
    return cc


##
## Other robust stats estimators
##

def robust_mean(x):
    y = x.flatten()
    n = len(y)
    y.sort()
    ind_qt1 = round((n+1)/4.)
    ind_qt3 = round((n+1)*3/4.)
    IQR = y[ind_qt3]- y[ind_qt1]
    lowFense = y[ind_qt1] - 1.5*IQR
    highFense = y[ind_qt3] + 1.5*IQR
    ok = (y>lowFense)*(y<highFense)
    yy=y[ok]
    return yy.mean(dtype='double')

#-------Robust Standard Deviation---

def robust_std(x):
    y = x.flatten()
    n = len(y)
    y.sort()
    ind_qt1 = round((n+1)/4.)
    ind_qt3 = round((n+1)*3/4.)
    IQR = y[ind_qt3]- y[ind_qt1]
    lowFense = y[ind_qt1] - 1.5*IQR
    highFense = y[ind_qt3] + 1.5*IQR
    ok = (y>lowFense)*(y<highFense)
    yy=y[ok]
    return yy.std(dtype='double')

#-------Robust variance---

def robust_var(x):
    y = x.flatten()
    n = len(y)
    y.sort()
    ind_qt1 = round((n+1)/4.)
    ind_qt3 = round((n+1)*3/4.)
    IQR = y[ind_qt3]- y[ind_qt1]
    lowFense = y[ind_qt1] - 1.5*IQR
    highFense = y[ind_qt3] + 1.5*IQR
    ok = (y>lowFense)*(y<highFense)
    yy=y[ok]
    return yy.var(dtype='double')

#----line fit------------------------
"""
y = a+bx.
input:  y, x
return: a,b, sd_a, sd_b, R^2
This is used where there is no measurement errors on x, y
http://mathworld.wolfram.com/LeastSquaresFitting.html
http://en.wikipedia.org/wiki/R-squared
"""

def linefit(x,y):
    """
    a,b,SEa,SEb,chi2 = linefit(x)
    """
    n=len(x)
    SSxx = n*np.var(x,dtype='double')
    SSyy = n*np.var(y,dtype='double')
    ybar = np.mean(y,dtype='double')
    xbar = np.mean(x,dtype='double')
    SSxy = np.sum(x*y,dtype='double') - n*xbar*ybar
    b = SSxy/SSxx
    a = ybar - b*xbar
    s = np.sqrt((SSyy-SSxy**2/SSxx)/(n-2))
    SEa = s*np.sqrt(1/n+xbar**2/SSxx)
    SEb = s/np.sqrt(SSxx)
    SStot = np.sum((y-ybar)**2,dtype='double')
    f = a + b*x
    SSerr = np.sum((y-f)**2,dtype='double')
    Rsq = 1-SSerr/SStot
    return a,b,SEa,SEb,Rsq


#---weighted LSQ fit to straightline based on NR---
""" 
y = a +b x, 
input:  x, y, yerr
return: a, b, sd_a, sd_b, chi2
See NR chapter 15
"""

def linfit(xx,yy,yyerr):
    if len(yy[np.abs(yy)>0]) < 10:
        b=0
        a=0
        SEa=0
        SEb=0
        chi2=999
        return a,b,SEa,SEb,chi2
    else:
        ok=yyerr > 0
        x = xx[ok]
        y = yy[ok]
        yerr = yyerr[ok]
        n=len(x)
        S=np.sum(1/yerr**2,dtype='double')
        Sx = np.sum(x/yerr**2,dtype='double')
        Sy = np.sum(y/yerr**2,dtype='double')
        Sxx = np.sum(x**2/yerr**2,dtype='double')
        Sxy = np.sum((x*y)/(yerr**2),dtype='double')
        delta = S*Sxx - Sx**2
        a = (Sxx*Sy-Sx*Sxy)/delta
        b = (S*Sxy-Sx*Sy)/delta
        SEa = np.sqrt(Sxx/delta)
        SEb = np.sqrt(S/delta)
        chi2 = np.sum((y-(a + b*x))**2/yerr**2,dtype='double')/(n-2.)
        return a,b,SEa,SEb,chi2

def MAD(a, c=0.6745, axis=None):
    """
    Median Absolute Deviation along given axis of an array:

    median(abs(a - median(a))) / c

    c = 0.6745 is the constant to convert from MAD to std; it is used by
    default

    """
    import numpy.ma as ma
    
    a = ma.masked_where(a!=a, a)
    if a.ndim == 1:
        d = ma.median(a)
        m = ma.median(ma.fabs(a - d) / c)
    else:
        d = ma.median(a, axis=axis)
        # I don't want the array to change so I have to copy it?
        if axis > 0:
            aswp = ma.swapaxes(a,0,axis)
        else:
            aswp = a
        m = ma.median(ma.fabs(aswp - d) / c, axis=0)

    return m
