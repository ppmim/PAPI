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
import config
import time


class clock:

    def __init__(self):
        self.start=0

    def tic(self):
        self.start=time.time()

    def tac(self):
        lap=time.time()-self.start
        return "Elapsed time(s): %f" %lap
    
