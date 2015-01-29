################################################################################
#
# PANICtool
#
# utils.py
#
# Last update 10/04/2008
#
################################################################################

"""
   Utils module containing some useful class and funtions used in PANICtool
"""

# Import necessary modules

# System modules
import os
import time


class clock:

    def __init__(self):
        self.start=0

    def tic(self):
        self.start=time.time()

    def tac(self):
        lap=time.time()-self.start
        return "Elapsed time(s): %f" %lap
    

def makemasterframe(list_or_array):
#from http://www.astro.ucla.edu/~ianc/python/_modules/ir.html#baseObject
    """
    If an array is passed, return it.  Otherwise, return a
    median stack of the input filename list.
    """
        
    if hasattr(list_or_array, 'shape') and len(list_or_array.shape)>1:
        masterframe = np.array(list_or_array, copy=False)
    else:
        masterframe = np.median(map(fits.getdata, list_or_array), axis=0)

    return masterframe
        