#!/usr/bin/env python
################################################################################
#
# PANICtool
#
# applyDarkFlat.py
#
# Created    : 05/06/2009    jmiguel@iaa.es
# Last update: 05/06/2009    jmiguel@iaa.es
#              11/03/2010    jmiguel@iaa.es - Added out_dir for output files
#              18/03/2010    jmiguel@iaa.es - Modified to support only dark subtraction, only flatfielding or both 
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

import misc.fileUtils
import misc.utils as utils

# Pyraf modules
from pyraf import iraf
from iraf import noao
from iraf import imred
from iraf import ccdred

import numpy

# Interact with FITS files
import pyfits

# Import Pyro core
import Pyro.core
import Pyro.naming

# Logging
from misc.paLog import log

class ApplyDarkFlat:
    """
    \brief Class used to subtract a master dark and then divide by a master flat field
    
    \par Class:
        ApplyDarkFlat
    \par Purpose:
        Apply a master dark and master flat to a list of non-calibrated science files
    \par Description:
        For each file processed a new file it is generated with the same filename but with the suffix '_DF.fits'
    \par Language:
        PyRaf
    \param data
        A list of science files
    \param bpm
        Input bad pixel mask or NULL
    \param mdark
        Master dark to subtract
    \param mflat
        Master flat to divide by
    \retval file_list
        If no error, return the list of files generated as result of the current processing
    \author
        JMIbannez, IAA-CSIC
  
    """
    def __init__(self, sci_files, mdark=None, mflat=None, out_dir="/tmp/", bpm=None):
        self.__sci_files = sci_files
        self.__mdark = mdark
        self.__mflat = mflat
        self.__bpm = bpm
        self.__out_dir = out_dir
      
    
    def apply(self):
      
        """
        \brief Apply masters to sci file list
        """   
        log.debug("Start applyDarkFlat")
        
        start_time = time.time()
        t=utils.clock()
        t.tic()
    
        out_suffix=".fits"
        # Master DARK reading:
        if self.__mdark!=None:
            if not os.path.exists(self.__mdark):      # check whether input file exists
                log.error('File %s does not exist', self.__mdark)
                sys.exit(1)
            else:
                dark = pyfits.open(self.__mdark)
                dark_time=dark[0].header['EXPTIME']
                dark_data=dark[0].data
                out_suffix=out_suffix.replace(".fits","_D.fits")
        # No master dark privided, then no dark subtracted
        else:
            dark_data=0
            dark_time=None
            
        # Master FLAT reading
        if self.__mflat!=None:    
            if not os.path.exists(self.__mflat):      # check whether input file exists
                log.error('File %s does not exist', self.__mflat)
                sys.exit(2)
            else:
                flat = pyfits.open(self.__mflat)
                flat_time=flat[0].header['EXPTIME']
                flat_filter=flat[0].header['FILTER']
                flat_data=flat[0].data
                out_suffix=out_suffix.replace(".fits","_F.fits")
        else:
            flat_data = 1.0
            flat_time = None
            flat_filter = None 
               
               
        # Get the user-defined list of flat frames
        framelist = self.__sci_files
       
        
        # STEP 2: Check the  TYPE and FILTER of each science file
        # If any frame on list missmatch the FILTER, then the procedure will be aborted
        # EXPTIME do not need be the same, so EXPTIME scaling will be done
        n_removed=0
        result_file_list = [] # List of files generared as result of this procedure and that will be returned
            
        for iframe in framelist:
            if not os.path.exists(iframe):
                log.error("File '%s' does not exist", iframe)
                continue  
            f = pyfits.open(iframe)
            debug=False
            if (debug):
              print "Science frame %s EXPTIME= %f TYPE= %s FILTER= %s" %(iframe, f[0].header['EXPTIME'],f[0].header['OBJECT'], f[0].header['FILTER'])
            # Check FILTER
            if ( flat_filter!=None and f[0].header['FILTER']!=flat_filter ):
                log.error("Error: Task 'applyDarkFlat' found a frame with different FILTER. Skipping frame...")
                f.close()
                n_removed=n_removed+1
            else:
                # Delete old files
                (path,name)=os.path.split(iframe)
                newpathname=self.__out_dir+"/"+name.replace(".fits", out_suffix)
                misc.fileUtils.removefiles(newpathname)
                s_time=float(f[0].header['EXPTIME'])
                # STEP 2.1: Check EXPTIME and apply master DARK and master FLAT
                if dark_time!=None and s_time!=dark_time:
                    log.debug("Scaling master dark ...")
                    f[0].data = (f[0].data - dark_data*float(f[0].header['EXPTIME']/dark[0].header['EXPTIME']) )/flat_data
                    f[0].header.add_history('Dark subtracted (scaled) %s' %self.__mdark)
                    if flat_time!=None: f[0].header.add_history('Flat-Fielding with %s' %self.__mflat)
                else:
                    log.debug("Not master dark EXPTIME scaling required")
                    f[0].data = (f[0].data - dark_data)/ flat_data
                    print "Flat data=", flat_data
                    if dark_time!=None: f[0].header.add_history('Dark subtracted %s' %self.__mdark)
                    if flat_time!=None: f[0].header.add_history('Flat-Fielding with %s' %self.__mflat)        
                            
                # Write output to outframe (data object actually still points to input data)
                try:
                    f[0].scale('float32')
                    
                    f.writeto(newpathname, output_verify='ignore')
                except IOError:
                    raise ExError('Cannot write output to %s' % newpathname)
                     
                f.close()
                result_file_list.append(newpathname)            
                log.debug('Saved  new dark subtracted or/and flattened file  %s' ,  newpathname )
        
        log.debug(t.tac() )
        log.info("Successful end of applyDarkFlat !")
                
        return result_file_list
        
        
        
        
################################################################################
# Functions       
def usage ():
    print ''
    print 'NAME'
    print '       applyDarkFlat.py - Apply Dark, Flat or both\n'
    print 'SYNOPSIS'
    print '       applyDarkFlat.py [options]\n'
    print 'DESCRIPTION'
    print '       Subtract a master dark file or divide by a'
    print '       master flat field (or both) a list of science files. '
    print ' '
    print 'OPTIONS'
    print "-s / --source=      Source file list of data frames"
    print "-d / --dark=        Master dark frame to subtract (optional)"
    print "-f / --flat=        Master flat field to divide by (optional)"
    print "-v                  Verbose debugging output\n"
    print 'VERSION'
    print '       2009 June 05'
    print ''
    raise SystemExit

################################################################################
# main
if __name__ == "__main__":
    
    # Get and check command-line options
    args = sys.argv[1:]
    source_file_list =None
    dark_file =None
    flat_file =None
    out_dir="/tmp/"
    debug=False
    
    try:
        opts, args = getopt.getopt(args, "s:d:f:o:v", ['source=','dark=', 'flat=', 'out_dir='])
    except getopt.GetoptError:
        # print help information and exit:
        usage()
        sys.exit(1)

    
    for option, parameter in opts:
        if option in ("-s", "--source"):
            source_file_list = parameter
            print "Source file list =", source_file_list
            if not os.path.isfile(source_file_list):
                print 'Error, file list does not exists :', source_file_list
                sys.exit(1)
        if option in ("-d", "--dark"):
            dark_file = parameter
            print "Dark file =", dark_file
        if option in ("-f", "--flat"):
            flat_file = parameter
            print "Flat file =", flat_file
            
        if option in ("-o", "--out_dir"):
            out_dir = parameter
            print "Out dir =", out_dir
                
        if option in ("-v"):
            debug=True
            
    if  source_file_list==None or (flat_file==None and dark_file==None):
        usage()
        sys.exit(3)
    
    filelist=[line.replace( "\n", "") for line in fileinput.input(source_file_list)]
    #filelist=['/disk-a/caha/panic/DATA/ALHAMBRA_1/A0408060036.fits', '/disk-a/caha/panic/DATA/ALHAMBRA_1/A0408060037.fits']
    print "Files:",filelist
    res = ApplyDarkFlat(filelist, dark_file, flat_file, out_dir)
    res.apply()
    
