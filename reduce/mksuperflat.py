#!/usr/bin/env python
################################################################################
#
#
# PAPI (PANIC PIpeline)
#
# mksuperflat.py
#
# Compute a super sky flat using the dither frames (IRAF implementation)
#
# Created    : 13/03/2009    jmiguel@iaa.es -
# Last update: 15/04/2009    jmiguel@iaa.es - Created function and modified to accept command line arguments
#              03/03/2010    jmiguel@iaa.es - Big modification to convert to a class and make more checkings
#
# TODO
#  
################################################################################

################################################################################
# Import necessary modules

import getopt
import sys
import os
import logging
import fileinput
import time
from datetime import datetime
from optparse import OptionParser

import misc.fileUtils
import misc.utils as utils


# Interact with FITS files
import pyfits
import datahandler

# Import Pyro core
import Pyro.core
import Pyro.naming

# Logging
from misc.paLog import log

# Pyraf modules
from pyraf import iraf
from iraf import noao
from iraf import imred
from iraf import ccdred


class SuperSkyFlat:
    """
    \brief Class used to build a super sky Flat from a dither set of science frames containing objects 
    
    \par Class:
        SuperSkyFlat
    \par Purpose:
         Create a super flat field 
    \par Description:
            
    \par Language:
        PyRaf
    \param file_list
        A list FITS files or directory
    \param output_filename
        File where log will be written
    \retval 0
        If no error
    \author
        JMIbannez, IAA-CSIC
        
    """
    def __init__(self,  objects_file,  output_filename="/tmp/superFlat.fits",  bpm=None):
         
        self.objects_file = objects_file
        self.output_file_dir = os.path.dirname(output_filename)
        self.output_filename = output_filename  # full filename (path+filename)
        self.bpm = bpm
        
        self.m_MIN_N_GOOD=2
        self.m_min_flats=5
            
    def create(self):
      
        """
        \brief Create a log sheet text file from a set of FITS files
        """   
        log.debug("Start createSuperSkyFlat") 
        if os.path.exists(self.output_filename): os.remove(self.output_filename)
        
        # Check data integrity (all have the same properties)
        if os.path.exists( self.objects_file ):
            m_filelist=[line.replace( "\n", "") for line in fileinput.input(self.objects_file)]
        else:
            m_filelist=self.objects_file
            
        if not datahandler.checkDataProperties( m_filelist ):
            log.error("Data integrity ERROR, some files not having same properties")
            raise Exception("Found a data integrity error")  
        
        # Combine the images to find out the super Flat
        iraf.imcombine(input=("'"+"@"+self.objects_file+"'").replace('//','/'),
                    output=self.output_filename,
                    combine='median',
                    offset='none',
                    reject='sigclip',
                    lsigma=2.5,
                    hsigma=2.5,
                    scale='median',
                    zero='none'
                    #masktype='none'
                    #scale='exposure',
                    #expname='EXPTIME'
                    #ParList = _getparlistname ('flatcombine')
                )
        
        median = float(iraf.imstat (
                images=("'"+self.output_filename+"[100:900,100:900]'").replace('//','/'),
                fields='midpt',format='no',Stdout=1)[0])
                                                
        """iraf.imarith(operand1 = 'sflat.fits',
                        operand2 = median,
                        op = '/',
                        result ='sflatn.fits',
                        verbose = 'yes'
                        )                                              
        """
        flatframe = pyfits.open(self.output_filename,'update')
        #Add a new keyword-->DATAMODE
        flatframe[0].header.update('DATAMODE',median,'Data mode of the frame')
        flatframe.close(output_verify='ignore') # This ignore any FITS standar violation and allow write/update the FITS file
                
                           
                      
################################################################################
# main
if __name__ == "__main__":
    print 'Start SuperSkyFlat....'
    # Get and check command-line options
    
    
    usage = "usage: %prog [options] arg1 arg2"
    parser = OptionParser(usage)
    parser.add_option("-v", "--verbose",
                  action="store_true", dest="verbose", default=True,
                  help="verbose mode [default]")
    parser.add_option("-s", "--source",
                  action="store", dest="source_file_list",
                  help="Source file list of data frames. It has to be a fullpath file name")
    parser.add_option("-o", "--output",
                  action="store", dest="output_filename", help="output file to write SuperFlat")
    


    (options, args) = parser.parse_args()
    if not options.source_file_list or not options.output_filename or len(args)!=0: # args is the leftover positional arguments after all options have been processed
        parser.print_help()
        parser.error("incorrect number of arguments " )
    if options.verbose:
        print "reading %s ..." % options.source_file_list
    
    
    superflat = SuperSkyFlat(options.source_file_list, options.output_filename)
    superflat.create()
          
        
        