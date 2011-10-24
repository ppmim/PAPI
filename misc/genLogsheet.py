#!/usr/bin/env python
################################################################################
#
# PANICtool
#
# genLogsheet.py
#
# Created    : 09/12/2009    jmiguel@iaa.es
# Last update: 09/12/2009    
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

class LogSheet:
    """
    \brief Class used to build a log sheet from a set of FITS files 
    
    \par Class:
        LogSheet
    \par Purpose:
        Create a log sheet from a set of FITS files
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
    def __init__(self, file_list, output_filename="/tmp/logsheet.txt"):
        self.__file_list=file_list
        self.__output_filename=output_filename  # full filename (path+filename)
        
            
    def create(self):
      
        """
        @summary:  Create a log sheet text file from a set of FITS files
        """   
        log.debug("Start createLogSheet")
        
        
        # STEP 0:Get the user-defined list of frames
        if os.path.isdir(self.__file_list):
            # Input is a directory
            filelist = [os.path.join(self.__file_list, file) for file in os.listdir(self.__file_list)]
        else:    
            # Input is a file
            filelist = [line.replace( "\n", "") for line in fileinput.input(self.__file_list)]
            
        # STEP 1: Create the logsheet text file    
        logsheet = open(self.__output_filename,"w")
        logsheet.write("#-------------------------------------------------------------------------\n")
        logsheet.write("#LOG SHEET created on %s (sorted by DATE_OBS)\n" %(datetime.now()))
        logsheet.write("#-------------------------------------------------------------------------\n") 
        logsheet.write('#%4s  %-32s  %-12s    %-20s  %-12s  %-12s  %-10s  %-10s  %-10s  %-20s\n' % (" ID", "Filename", "Filter", "Type", "RA", "Dec", "TEXP", "NCOADDS","ITIME","DATE_OBS"))
        logsheet.write("#-------------------------------------------------------------------------------------------------------------------------------------------------------------------\n") 
        
        # STEP 2: Determine the number of darks frames to combine
        id=0
        sorted_list = {}
        for file in filelist:
            try:
                fitsf = datahandler.ClFits ( file )
                sorted_list[fitsf.getDateTimeObs()] = fitsf
            except:
                log.error( "Unexpected error reading FITS file %s" %file )
                
        keys = sorted_list.keys()
        keys.sort()
        for key in keys:
            f = sorted_list[key]
            logsheet.write('%4d  %-32s  %-12s  %-20s  %-12f  %-12f  %-10f  %-10d  %-10f  %-20s\n' % (id, os.path.basename(f.filename), f.getFilter(), f.getType(), f.ra, f.dec, f.expTime(), f.getNcoadds(), f.getItime(), f.getDateTimeObs()) )
            id = id+1    
                         
        logsheet.close()
        log.debug('Saved logsheet to %s' , self.__output_filename)
        
        return self.__output_filename
    
    def show (self):
        #TODO
        os.system("cat %s" %(self.__output_filename))
        pass
                
################################################################################
# main
################################################################################
if __name__ == "__main__":
    print 'Start Log Sheet generator....'
    # Get and check command-line options
    
    
    usage = "usage: %prog [options] arg1 arg2"
    parser = OptionParser(usage)
    parser.add_option("-v", "--verbose",
                  action="store_true", dest="verbose", default=True,
                  help="verbose mode [default]")
    parser.add_option("-s", "--source",
                  action="store", dest="source_file_list",
                  help="Source file list of data frames. It can be a file or directory name.")
    parser.add_option("-o", "--output file for logsheet",
                  action="store", dest="output_filename", help="write output logsheet to specified file")
    parser.add_option("-d", "--display",
                  action="store_true", dest="show", default=False,
                  help="show result on screen (stdout)")


    (options, args) = parser.parse_args()
    if not options.source_file_list or not options.output_filename or len(args)!=0: # args is the leftover positional arguments after all options have been processed
        parser.print_help()
        parser.error("incorrect number of arguments " )
    if options.verbose:
        print "reading %s ..." % options.source_file_list
    
    
    logsheet = LogSheet(options.source_file_list, options.output_filename)
    logsheet.create()
    if options.show==True:
        logsheet.show()
        
