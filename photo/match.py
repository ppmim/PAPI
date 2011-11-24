#!/usr/bin/python
# match.py
# Utilties for matching catalogs by RA & Dec or ID
# H. Ferguson 11/22/05

__version__ = 2.0 # Updated for Numpy 5/15/08 not tested yet
__author__ = "H. C. Ferguson"

from numpy import *
import angsep
import string
import sys
import copy

def matchsorted(ra,dec,ra1,dec1,tol):
    """ Find closest ra,dec within tol to a target in an ra-sorted list of ra,dec.
        Arguments:
          ra - Right Ascension decimal degrees (numpy sorted in ascending order)
          dec - Declination decimal degrees (numpy array)
          ra1 - RA to match (scalar, decimal degrees)
          ra1 - Dec to match (scalar, decimal degrees)
          tol - Matching tolerance in decimal degrees. 
        Returns:
          ibest - index of the best match within tol; -1 if no match within tol
          sep - separation (defaults to tol if no match within tol)
    """
    i1 = searchsorted(ra,ra1-tol)-1
    i2 = searchsorted(ra,ra1+tol)+1
    if i1 < 0:
        i1 = 0
    sep = angsep.angsep(ra[i1:i2],dec[i1:i2],ra1,dec1)
    # print "tolerance ",tol
    indices=argsort(sep)
    # print sep
    if sep[indices[0]] > tol:
        return -1, tol
    ibest = indices[0] + i1
    return ibest, sep[indices[0]]
    
def matchpos(ra1,dec1,ra2,dec2,tol):
    """ Match two sets of ra,dec within a tolerance.
        Longer catalog should go first
        Arguments:
          ra1 - Right Ascension decimal degrees (numpy array)
          dec1 - Declination decimal degrees (numpy array)
          ra2 - Right Ascension decimal degrees (numpy array)
          dec2 - Declination decimal degrees (numpy array)
          tol - Matching tolerance in decimal degrees. 
        Returns:
          ibest - indices of the best matches within tol; -1 if no match within tol
          sep - separations (defaults to tol if no match within tol)
    """
    indices = argsort(ra1)
    rasorted = ra1[indices]
    decsorted = dec1[indices]
    ibest = []
    sep = []
    for i in range(len(ra2)):
        j,s = matchsorted(rasorted,decsorted,ra2[i],dec2[i],tol)
        if j < 0:
            ibest += [j]
        else:
            ibest += [indices[j]]
        sep += [s]
    return array(ibest),array(sep)

def matchjoin(si1,si2,ibest,sep=[],dict1={},dict2={}):
    """ Keep only elements that match in both catalogs. 
        Arguments:
          si1 -- object with data arrays as attributes (e.g. si1.mag_auto, si1.id)
          si2 -- object with data arrays as attributes 
          ibest -- indices of si1 that match si2 (in order of si2)
        Keyword Arguments:
          sep -- array of separations, will be added as an attribute to second object
          dict1 -- dictionary of attributes for si1 (e.g. {'mag_auto':0,'id':1} )
          dict2 -- dictionary of attributes for si2 
        The objects si1 and si2 would normally be sextractor objects returned
        by sextutils.sextractor(). In that case, the column names are stored in
        si1._d and si2._d. If the objects are not sextractor objects, the user can
        provide separate dictionaries whose keys are the object attributes 
	that correspond to the numpy data arrays. 
       
        Returns:
          s1, s2 -- object arrays that include only the matching objects contained
              in both s1 and s2.
    """
    if len(dict1) == 0:
        dict1 = si1._d
    if len(dict2) == 0:
        dict2 = si2._d
    indices = compress(ibest > 0,ibest)
    # print indices
    s1 = copy.deepcopy(si1)
    s2 = copy.deepcopy(si2)
    for k in dict1.keys():
        if type(s1.__dict__[k]) == type([]):
             s1.__dict__[k] = char.array(s1.__dict__[k])
        s1.__dict__[k] = s1.__dict__[k][indices]
    flag = ibest > 0
    for k in dict2.keys():
        if type(s2.__dict__[k]) == type([]):
             s2.__dict__[k] = char.array(s2.__dict__[k])
        s2.__dict__[k] = compress(flag,s2.__dict__[k])
    if len(sep) > 0:
        s2.separation = compress(flag,sep)
    return s1,s2

def matchids(id1,id2):
    """ Match two sets of ids. 
        Returns: 
          ibest -- array of indices of i1 that match i2; -1 if no match
    """
    indices = argsort(id1)
    idsorted = id1[indices]
    ibest = []
    for i in range(len(id2)):
        j = matchidsorted(idsorted,id2[i])
        if j < 0:
            ibest += [j]
        else:
            ibest += [indices[j]]
    return array(ibest)

def matchidsorted(ids,targetid):
    """ Find id matches, return index in i1 that matches targetid; -1 if no match. """
    i1 = searchsorted(ids,targetid)
    if targetid == ids[i1]:
        ibest = i1
    else:
        ibest = -1 
    return ibest
