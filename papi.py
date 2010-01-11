#!/usr/bin/env python

################################################################################
#
# papi (PAnic PIpeline)
#
# papi.py
#
# Last update 09/March/2009
#
################################################################################


    ####################################################################
    #                                                                  #
    #                   This is where it all starts.                   #
    #         The body of this program is at the end of this file      #
    #                                                                  #
    ####################################################################
    
from qt import *

import sys
import os
import os.path
import fnmatch
import time
import datahandler

#import tasks
# PANICtools threads
import reduce
import misc.fileUtils

# Interact with FITS files
import pyfits

# IRAF packages
from pyraf import iraf
from iraf import noao
from iraf import imred
from iraf import ccdred

# Math module for efficient array processing
import numpy

#Log
import misc.paLog
from misc.paLog import log    


    
    
################################################################################            
#  Testing
################################################################################
def usage():
    """Print help """
    print "Usage: %s filename" %(sys.argv[0])
    
if __name__ == "__main__": 
    log.debug("Starting the tests.....")
    
    if len(sys.argv)>1:
        try:
            print "HOLA!!"
        except:
            log.error("Error reading fits File")
    else:
        usage()