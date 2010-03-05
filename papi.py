#!/usr/bin/env python
__date__ = "$Date$"
__author__ = "$Author$"
__revision__ = "$Rev$"

################################################################################
#
# papi (PAnic PIpeline)
#
# papi.py
#
# Last update 04/March/2010
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
from optparse import OptionParser

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

#PAPI packages 
import datahandler
import reduce
import misc.fileUtils

class ReductionSet:
    def __init__(self, list_file, out_dir):
        """ Init function """
        
        self.list_file = list_file
        self.out_dir = out_dir
        
         
    def reduce(self):
        log.debug("Start data reduction")
        
        # Clean old files and copy new directories  ------------------------------
        #\rm *.fits
        $PAPI_HOME/cleanpapi
        $PAPI_HOME/copySource.py $1
        
        log.debub("End of data reduction")
        
    
################################################################################
# main
################################################################################
if __name__ == "__main__":
    log.debug( 'Start PAPI....')
    # Get and check command-line options
    
    
    usage = "usage: %prog [options] arg1 arg2 ..."
    parser = OptionParser(usage)
    
                  
    parser.add_option("-s", "--source",
                  action="store", dest="source_file_list",
                  help="Source file list of data frames. It can be a file or directory name.")
    
    parser.add_option("-o", "--output",
                  action="store", dest="output_filename", help="final reduced output image")
    
    parser.add_option("-t", "--type",
                  action="store", dest="type", default="quick", help="type of reduction (quick|science)")
                  
    parser.add_option("-m", "--obs_mode",
                  action="store", dest="obs_mode", default="dither", help="observing mode (dither|ext_dither)")
    
    parser.add_option("-d", "--outdir",
                  action="store_true", dest="out_dir", default=False,
                  help="output dir for intermidiate files")
    
    parser.add_option("-v", "--verbose",
                  action="store_true", dest="verbose", default=True,
                  help="verbose mode [default]")
    
    parser.add_option("-C", "--config_file",
                  action="store_true", dest="config_file", help="config file for the data reduction process")
                  
    parser.add_option("-S", "--show",
                  action="store_true", dest="show", default=False, help="show final reduced image")
                  
                  
    (options, args) = parser.parse_args()
    if not options.source_file_list or not options.output_filename or len(args)!=0: # args is the leftover positional arguments after all options have been processed
        parser.print_help()
        parser.error("incorrect number of arguments " )
    if options.verbose:
        print "reading %s ..." % options.source_file_list
    
    
    redSet = ReductionSet(options.source_file_list, options.output_filename)
    redSet.reduce()
    
    if options.show==True:
        redSet.show()    
