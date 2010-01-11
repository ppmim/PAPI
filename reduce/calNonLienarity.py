#!/usr/bin/env python
################################################################################
#
# PANICtool
#
# calNonLinearity.py
#
# Created    : 16/12/2009    jmiguel@iaa.es
# Last update: 
#
################################################################################
# Import necessary modules

import getopt
import sys
import os
import logging
import fileinput
import time
from optparse import OptionParser

import misc.fileUtils
import misc.utils as utils
import datahandler

# Interact with FITS files
import pyfits
import numpy

# Logging
from misc.paLog import log

class NonLinearityModel:
    """
    \brief Class used to compute the Non-linearity model for the detectors
    
    \par Class:
         NonLinearityModel   
    \par Purpose:
        Compute the non-linearity coefficients of the detector
    \par Description:
        Based on the algorithm described by T.H.Jarrett et all (1994) we compute the non-linearity coefficientes
        of the model
    \par Language:
        PyRaf
    \author
        JMIbannez, IAA-CSIC
    \todo 
           
    """
    def __init__(self, input_files, output_dir, output_filename="/tmp/NLC.fits"):
        """
        \brief Init the class
        
        \param input_data
            A list of sky files
        \param input_files list of FITS files as input for the computation
        \param output_filename File where coefficientes will be saved
        \retval coeff Array of values a0,a1,a2 where 
            a0 = bias level
            a1 = quantum efficiency or sensitivity
            a2 = the non-linearity or saturation
        """
        self.__input_files=input_files
        self.__output_file_dir=output_dir
        self.__output_filename=output_filename  # full filename (path+filename)
    
    def createModel(self):
          
        """
        \brief Create the Non-linearity correction model
        """   
        log.debug("Start createModel")
        start_time = time.time()
        t=utils.clock()
        t.tic()
        
        # Get the user-defined list of dark frames
        framelist = self.__input_files
        
        # STEP 0: Determine the number of frames 
        try:    
            nframes = len(framelist)
        except IndexError:
            log.error("No DARK frames defined")
            raise
        
        if not os.path.exists(os.path.dirname(self.__output_filename)):
            raise NameError, 'Wrong output path'
        if not self.__output_filename:
            log.error("Combined DARK frame not defined")
            raise "Wrong output filename"
    
        log.debug('Saved NLC Model to %s' , self.__output_filename)
        log.debug("createModel' finished %s", t.tac() )
        
        return True

################################################################################
# main
if __name__ == "__main__":
    
    
    usage = "usage: %prog [options]"
    parser = OptionParser(usage)
    
    parser.add_option("-s", "--source",
                  action="store", dest="source_file_list",
                  help="Source file list of FITS files")
    
    parser.add_option("-l", "--limit",
                  action="store", dest="satur_lim", default=40000, 
                  help="saturation limit")
    
    parser.add_option("-o", "--out_data",
                  action="store", dest="out_data", 
                  help="filename of out data file, contains (...)")
    
    parser.add_option("-c", "--coeff_file",
                  action="store", dest="out_coeff_file", 
                  help="filename of outputs coeffs, contains (a0, a1, a2)")
    
    parser.add_option("-t", "--ref_time",
                  action="store_true", dest="ref_time", default=False,
                  help="exptime used as reference")

    parser.add_option("-v", "--verbose",
                  action="store_true", dest="verbose", default=False,
                  help="verbose mode [default]")
    
    (options, args) = parser.parse_args()
    
    #Check required parameters
    if (not options.source_file_list or not options.satur_lim or not options.out_data or not options.out_coeff_file or not options.ref_time or len(args)!=0): # args is the leftover positional arguments after all options have been processed
        parser.print_help()
        parser.error("incorrect number of arguments " )
    if options.verbose:
        print "reading %s ..." % options.source_file_list
        
    filelist=[line.replace( "\n", "") for line in fileinput.input(options.source_file_list)]
    NLM = NonLinearityModel(filelist,"/tmp",options.out_coeff_file)
    NLM.createModel()
    
        