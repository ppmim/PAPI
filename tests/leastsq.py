'''
@file leastsq.py
@author Douglas Applegate
@date 9/7/07

@brief Functions to fit parameters of functions
'''

LEASTSQ_CVS_ID = "$Id: leastsq.py,v 1.3 2009-08-10 18:38:00 dapple Exp $"

##################################################

import numpy
from scipy import optimize
from scipy import linalg

import sys




################################################

def _prepData(points, vals, errs):
    if not isinstance(points,numpy.ndarray):
        points = numpy.array(points)

    if not isinstance(vals,numpy.ndarray):
        vals = numpy.array(vals)
    
    if errs is not None and not isinstance(errs,numpy.ndarray):
        errs = numpy.array(errs)
        
    return points, vals, errs


################################################

def leastsq(func, guess, points, vals, errs=None, fullOutput=False):
    ''' Performs non-linear least squares on a function & datasate

        @param func(x,params) where x is (possibly) a numpy array
                and params is a list. Assumes func can handle numpy
                vectors. (Use numpy.vectorize if in doubt)
        @param guess list of initial param values.
        @param points[i][j] coordinates of data points, where
                points[i] are individual datapoints and 
                points[i][j] are components of datapoint position
        @param vals value of data and each position in points
        @param errs error in vals at each position in points
        @param fullOutput selects how much info to return

        @returns params, {residual, covar}, isConverged 
                 where:
                  params is a list of best fit parameters
                  chisq is chisq or equivalent Student-t distribution
                  covar is the covariance matrix
                  isConvered is a bool describing convergence
    '''

    points, vals, errs = _prepData(points, vals, errs)

    def residuals(params,x,y,errs):
        predicted = func(x,params)
        if errs is None:
            return y - predicted
        else:
            return numpy.divide(y - predicted, errs)

    params, covar, info, mesg, ier = optimize.leastsq(residuals, guess,
                                             args = (points,vals,errs),
                                             full_output = True)

    if not hasattr(params,'__getitem__'):
        params = [params]
    
    isConverged = (ier == 1)


    if fullOutput:
        chisq = numpy.sum(residuals(params,points,vals,errs)**2)
        return params, chisq, covar, isConverged
    else:
        return params, isConverged

#########################################

def lin_leastsq(model, points, vals, errs = None, fullOutput = False,
                **keywords):
    ''' Performs linear least squares on a function & dataset
    
        @param model function to be fit to.
            Contains function basis(), which is an array of funcs
            that take points as arguments to form config matrix
            Or can be a list of said functions
        @param points[i][j] coordinates of data points, where
            points[i] are individual datapoints and 
            points[i][j] are components of datapoint position
        @param vals value of data and each position in points
        @param errs error in vals at each position in points
        @param fullOutput selects how much info to return
    
        @returns params, {residual, covar}, isConverged 
            where:
                params is a list of best fit parameters
                chisq is chisq or equivalent Student-t distribution
                covar is the covariance matrix
                isConvered is a bool describing convergence (always true)
    '''


    ###################

    def makeMatrixs(points, vals, errs):

        A = numpy.zeros((len(points),len(basis)))
        if errs is None:
            for i in xrange(len(points)):
                for j in xrange(len(basis)):
                    A[i,j] = basis[j](points[i])
            b = vals
        else:
            b = numpy.zeros((len(points),1))
            for i in xrange(len(points)):
                for j in xrange(len(basis)):
                    A[i,j] = (basis[j](points[i]))/errs[i]
                b[i] = vals[i]/errs[i]
        return A, b

    #####################

    def calcCovarMatrix(A, params):
        u,s,v = linalg.svd(A)

        covar = numpy.zeros((len(params),len(params)))
        for i in xrange(len(params)):
            for j in xrange(len(params)):
                for k in xrange(len(params)):
                    covar[i,j] += v[i,k]*v[j,k]/s[k]
                covar[j,i] = covar[i,j]
        return covar

    #####################

    def calcResiduals(params,x,y,errs):
        predicted = model(x,params)
        if errs is None:
            return y - predicted
        else:
            return numpy.divide(y - predicted, errs)

    #####################

    if hasattr(model, 'basis'):
        basis = model.basis()
    else:
        basis = model

    points, vals, errs = _prepData(points, vals, errs)


    A,b = makeMatrixs(points, vals, errs)

    (params, resids, rank, s) = linalg.lstsq(A, b)

    if fullOutput:
        
        covar = calcCovarMatrix(A, params)
        
        chisq = numpy.sum(calcResiduals(params,points,vals,errs)**2)

        return params, chisq, covar, True

    return params, True
        
########################################

def linear_leastsq(points, vals, errs = None, fullOutput = False,
                **keywords):
    ''' Performs linear least squares assuming a linear model
    
        @param points[i][j] coordinates of data points, where
            points[i] are individual datapoints and 
            points[i][j] are components of datapoint position
        @param vals value of data and each position in points
        @param errs error in vals at each position in points
        @param fullOutput selects how much info to return
    
        @returns params, {residual, covar}, isConverged 
            where:
                params is a list of best fit parameters
                chisq is chisq or equivalent Student-t distribution
                covar is the covariance matrix
                isConvered is a bool describing convergence (always true)
    '''

    ###################

    def makeMatrixs(points, vals, errs):

        npoints = len(points)
        if len(points.shape) == 1:
            nparams = 2
        else:
            nparams = points.shape[1] + 1

        A = numpy.ones( (npoints, nparams)  )
        if nparams == 2:
            for i in xrange(npoints):
                A[i,0] = points[i]
        else:
            A[:,:-1] = points

        b = vals
        if errs is not None:            
            for i in xrange(npoints):
                A[i,:] = A[i,:] / errs[i]
            b = b / errs
        return A, b

    #####################

    def calcCovarMatrix(A, params):

        nparams = len(params)

        u,s,v = linalg.svd(A)

        covar = numpy.zeros((nparams,nparams))
        for i in xrange(nparams):
            for j in xrange(nparams):
                for k in xrange(nparams):
                    covar[i,j] += v[i,k]*v[j,k]/s[k]
                covar[j,i] = covar[i,j]
        return covar

    #####################

    def calcResiduals(params,A,y):
        predicted = numpy.dot(A,params)
        return y - predicted

    #####################

    points, vals, errs = _prepData(points, vals, errs)

    A,b = makeMatrixs(points, vals, errs)

    (params, resids, rank, s) = linalg.lstsq(A, b)

    if fullOutput:
        
        covar = calcCovarMatrix(A, params)
        
        chisq = numpy.sum(calcResiduals(params,A,vals)**2)

        return params, chisq, covar, True

    return params, True
        

    
