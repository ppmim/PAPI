from __future__ import division

import numpy as np
from numpy import isnan
import warnings

from scipy.special import gammainc
from scipy.optimize import leastsq
import matplotlib.pyplot as pl



def findchisq(func, x, y, ysig, a):
    """ Calculate chi2 for function and data values

    func: model function, func(x, *a) models y data
    x,y: data points
    ysig: 1 sigma errors on y
    a: model parameters
    """
    resid = (y - func(x, *a)) / ysig
    return np.dot(resid, resid)


def polyfitr(x, y, order=2, clip=6, xlim=None, ylim=None,
             mask=None, debug=False):
    """ Fit a polynomial to data, rejecting outliers.

    Fits a polynomial f(x) to data, x,y.  Finds standard deviation of
    y - f(x) and removes points that differ from f(x) by more than
    clip*stddev, then refits.  This repeats until no points are
    removed.

    Inputs
    ------
    x,y:
        Data points to be fitted.  They must have the same length.
    order: int (2)
        Order of polynomial to be fitted.
    clip: float (6)
        After each iteration data further than this many standard
        deviations away from the fit will be discarded.
    xlim: tuple of maximum and minimum x values, optional
        Data outside these x limits will not be used in the fit.
    ylim: tuple of maximum and minimum y values, optional
        As for xlim, but for y data.
    mask: sequence of pairs, optional
        A list of minimum and maximum x values (e.g. [(3, 4), (8, 9)])
        giving regions to be excluded from the fit.
    debug: boolean, default False
        If True, plots the fit at each iteration in matplotlib.

    Returns
    -------
    coeff, x, y:
        x, y are the data points contributing to the final fit. coeff
        gives the coefficients of the final polynomial fit (use
        np.polyval(coeff,x)).

    Examples
    --------
    >>> x = np.linspace(0,4)
    >>> np.random.seed(13)
    >>> y = x**2 + np.random.randn(50)
    >>> coeff, x1, y1 = polyfitr(x, y)
    >>> np.allclose(coeff, [1.05228393, -0.31855442, 0.4957111])
    True
    >>> coeff, x1, y1 = polyfitr(x, y, order=1, xlim=(0.5,3.5), ylim=(1,10))
    >>> np.allclose(coeff, [3.23959627, -1.81635911])
    True
    >>> coeff, x1, y1 = polyfitr(x, y, mask=[(1, 2), (3, 3.5)])
    >>> np.allclose(coeff, [1.08044631, -0.37032771, 0.42847982])
    True
    """

    good = ~np.isnan(x) & ~np.isnan(y)
    x = np.asanyarray(x[good])
    y = np.asanyarray(y[good])
    isort = x.argsort()
    x, y = x[isort], y[isort]

    keep = np.ones(len(x), bool)
    if xlim is not None:
        keep &= (xlim[0] < x) & (x < xlim[1])
    if ylim is not None:
        keep &= (ylim[0] < y) & (y < ylim[1])
    if mask is not None:
        badpts = np.zeros(len(x), bool)
        for x0,x1 in mask:
            badpts |=  (x0 < x) & (x < x1)
        keep &= ~badpts

    x,y = x[keep], y[keep]
    if debug:
        fig = pl.figure()
        ax = fig.add_subplot(111)
        ax.plot(x,y,'.')
        ax.set_autoscale_on(0)
        pl.show()

    coeff = np.polyfit(x, y, order)
    if debug:
        pts, = ax.plot(x, y, '.')
        poly, = ax.plot(x, np.polyval(coeff, x), lw=2)
        pl.show()
        raw_input('Enter to continue')
    norm = np.abs(y - np.polyval(coeff, x))
    stdev = np.std(norm)
    condition =  norm < clip * stdev
    y = y[condition]
    x = x[condition]
    while norm.max() > clip * stdev:
        if len(y) < order + 1:
            raise Exception('Too few points left to fit!')
        coeff = np.polyfit(x, y, order)
        if debug:
            pts.set_data(x, y)
            poly.set_data(x, np.polyval(coeff, x))
            pl.show()
            raw_input('Enter to continue')
        norm = np.abs(y - np.polyval(coeff, x))
        stdev = norm.std()
        condition =  norm < clip * stdev
        y = y[condition]
        x = x[condition]

    return coeff,x,y

def wleastsq(x, y, ysig=None):
    """ Calculate the line of best fit with weights.

    Input
    -----
      x : sequence of floats
          input x data
      y : sequence of floats, len(x)
          input y data, f(x).
      ysig : sequence of floats, len(x), optional
         If the y none sigma errors are given, points are weighted
         by their inverse variance.

    Returns
    -------
      (b, a),(b_sigma, a_sigma)
        The fitted parameters and their one sigma errors.  The fitted
        line has equation y = a + b*x.
    """
    if ysig == None:
        ysig = np.ones(len(x))
    yvar = ysig * ysig   # variance

    s = np.sum(1. / yvar)
    sx = np.sum(x / yvar)
    sxx = np.sum(x*x / yvar)
    sy = np.sum(y / yvar)
    sxy = np.sum(x*y / yvar)

    # See NR section 15.2 for a derivation of the below solutions for
    # the best fit values of a and b.
    # 
    # y = a + b*x 

    temp = s*sxx - sx*sx
    a = (sxx*sy - sx*sxy) / temp
    b = (s*sxy - sx*sy) / temp
    sig_a = np.sqrt(sxx / temp)
    sig_b = np.sqrt(s / temp)

    return (b,a),(sig_b,sig_a)


def delta_chi2(ptarg, ndeg):
    """ Find delta chi2 corresponding to the given probability for n
    degrees of freedom (number of data points - number of free
    parameters).

    ptarg is the target probability (ptarg=0.9 correponds to 90%).

    ptarg is the probability that a model fitted to data with the
    given number of degrees of freedom will have a chisq val
    delta_chi2 larger than the minimum chisq value for the
    best-fitting parameters.
    """
    assert 0 <= ptarg <= 1
    d0, d1 = 0, 1e5
    a = 0.5 * ndeg
    while True:
        d = 0.5*(d0 + d1)
        p = gammainc(a, 0.5*d)
        #print ptarg, p, d, d0, d1
        if p > ptarg:
            d1 = d
        else:
            d0 = d
        if abs(p - ptarg) < 1e-4:
            break

    return d


