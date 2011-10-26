#!/usr/bin/env python
""" This module implement the class ApplyDarkFlat """

# Copyright (c) 2009 Jose M. Ibanez. All rights reserved.
# Institute of Astrophysics of Andalusia, IAA-CSIC
#
# This file is part of PAPI (PANIC Pipeline)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

################################################################################
#
# PANICtool
#
# applyDarkFlat.py
#
# Created    : 05/06/2009    jmiguel@iaa.es
# Last update: 05/06/2009    jmiguel@iaa.es
#              11/03/2010    jmiguel@iaa.es - Added out_dir for output files
#              18/03/2010    jmiguel@iaa.es - Modified to support only dark 
#                            subtraction, only flatfielding or both 
#              15/11/2010    jmiguel@iaa.es - Added normalization to flat-field
#              16/11/2010    jmiguel@iaa.es - Added support for MEF files
################################################################################

################################################################################
# Import necessary modules

import getopt
import sys
import os
import fileinput
import time

import misc.fileUtils
import misc.utils as utils
import datahandler

# Pyraf modules
#from pyraf import iraf
#from iraf import noao
#from iraf import imred
#from iraf import ccdred

import numpy as np

# Interact with FITS files
import pyfits


# Logging
from misc.paLog import log


class ExError (Exception): 
    """ Next class if for a general execution Exception """
    pass

class ApplyDarkFlat:
    """
    \brief Class used to subtract a master dark and then divide by a master 
            flat field
    
    \par Class:
        ApplyDarkFlat
    \par Purpose:
        Apply a master dark and master flat to a list of non-calibrated science files
    \par Description:
        For each file processed a new file it is generated with the same filename but 
        with the suffix '_DF.fits'
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
        If no error, return the list of files generated as result of the current 
        processing
    \author
        JMIbannez, IAA-CSIC
  
    """
    def __init__(self, sci_files, mdark = None, mflat = None, \
                 out_dir_ ="/tmp", bpm = None):
        self.__sci_files = sci_files  # list of files which apply dark and flat
        self.__mdark = mdark          # master dark to apply
        self.__mflat = mflat          # master flat to apply (dome, twlight) - not normalized !!
        self.__bpm = bpm              # not used at the moment
        self.__out_dir = out_dir_
        
    
    def apply(self):
      
        """
        @summary: Apply masters DARK and/or FLAT to science file list. 
        Both master DARK and FLAT are optional,i.e., each one can be applied 
        even the other is not present. 
        """   

        log.debug("Start applyDarkFlat")
        
        start_time = time.time()
        t = utils.clock()
        t.tic()
    
        # some variables 
        out_suffix = ".fits"
        dmef = False # flag to indicate if dark is a MEF file or not
        fmef = False # flag to indicate if flat is a MEF file or not
        n_ext = 1 # number of extension of the MEF file (1=simple FITS file)
        mode = 1 # mode of the flat frame (if mef, mode of chip 0)
        
        # #######################
        # Master DARK reading:
        if self.__mdark != None:
            if not os.path.exists(self.__mdark): # check whether input file exists
                log.error('File %s does not exist', self.__mdark)
                sys.exit(1)
            else:
                dark = pyfits.open(self.__mdark)
                cdark = datahandler.ClFits (self.__mdark)
                dark_time = cdark.expTime()
                #dark_data=dark[0].data
                if len(dark)>1:
                    dmef = True
                    n_ext = len(dark)-1
                    log.debug("Dark MEF file with %d extensions", n_ext)
                out_suffix = out_suffix.replace(".fits","_D.fits")
                log.debug("Master DARK to be subtracted : %s"%self.__mdark)
        # No master dark provided, then no dark to subtract
        else:
            log.debug("No master dark to be subtracted !")
            dark_data = 0
            dark_time = None
            
       # ######################
       # Master FLAT reading
        if self.__mflat != None:    
            if not os.path.exists(self.__mflat): # check whether input file exists
                log.error('File %s does not exist', self.__mflat)
                sys.exit(2)
            else:
                flat = pyfits.open(self.__mflat)
                cflat = datahandler.ClFits (self.__mflat)
                flat_time = cflat.expTime()
                flat_filter = cflat.getFilter()
                #flat_data=flat[0].data
                #MEF
                if len(flat)>1:
                    fmef = True
                    if len(dark) != len(flat):
                        raise Exception("Number of extensions does not match \
                        in Dark and Flat files!")
                    else: n_ext = len(flat)-1    
                    log.debug("Flat MEF file with %d extensions", n_ext)
                    
                # compute mode for n_ext normalization (in case of MEF, we normalize 
                # all extension wrt chip 0)
                naxis1 = cflat.getNaxis1()
                naxis2 = cflat.getNaxis2()
                if len(flat)>1: 
                    ext = 1
                else: 
                    ext = 0
                median = np.median(flat[ext].data[200:naxis1-200, 200:naxis2-200])
                mean = np.mean(flat[ext].data[200:naxis1-200, 200:naxis2-200])
                mode = 3*median-2*mean
                log.debug("MEDIAN= %f  MEAN=%f MODE(estimated)=%f ", \
                           median, mean, mode)
                log.debug("Flat-field will be normalized by MODE ( %f ) value", mode)
                out_suffix = out_suffix.replace(".fits","_F.fits")
                log.debug("Found master FLAT to divide by: %s"%self.__mflat)
        else:
            log.debug("No master flat to be divided by")
            flat_data = 1.0
            flat_time = None
            flat_filter = None
            mode = 1 
        
        # Get the user-defined list of flat frames
        framelist = self.__sci_files
        
        # STEP 2: Check the  TYPE and FILTER of each science file
        # If any frame on list missmatch the FILTER, then the procedure will be 
        # aborted
        # EXPTIME do not need be the same, so EXPTIME scaling will be done
        n_removed = 0
        
        # List of files generated as result of this procedure and that will be returned
        result_file_list = [] 
            
        for iframe in framelist:
            if not os.path.exists(iframe):
                log.error("File '%s' does not exist", iframe)
                continue  
            f = pyfits.open(iframe)
            cf = datahandler.ClFits (iframe)
            log.debug("Science frame %s EXPTIME= %f TYPE= %s FILTER= %s"\
                      %(iframe, cf.expTime(), cf.getType(), cf.getFilter()))
            
            # Check FILTER
            if ( flat_filter != None and cf.getFilter() != flat_filter):
                log.error("Error: Task 'applyDarkFlat' found a frame with \
                different FILTER. Skipping frame...")
                f.close()
                n_removed = n_removed+1
            else:
                # check Number of Extension 
                if (len(f)>1 and (len(f)-1) != n_ext) or len(f) != n_ext:
                    raise Exception("File %s does not match the number of \
                    extensions (%d)", iframe, n_ext)
                
                # Delete old files
                (path, name) = os.path.split(iframe)
                newpathname = (self.__out_dir + "/" + \
                             name.replace(".fits", out_suffix)).replace("//","/")
                misc.fileUtils.removefiles(newpathname)
                
                # Scale master DARK
                i_time = float(cf.expTime()) # all extension have the same TEXP
                if self.__mdark != None and dark_time!=None:
                    time_scale = float(i_time / dark_time)
                else: time_scale = 1
                for chip in range (0, n_ext):
                    if n_ext > 1: # it means, MEF
                        if self.__mdark != None: dark_data = dark[chip+1].data
                        else: dark_data = 0
                        if self.__mflat!=None: flat_data = flat[chip+1].data/mode # normalization wrt chip 0
                        else: flat_data = 1
                        sci_data = f[chip+1].data 
                    else:
                        if self.__mdark != None: dark_data = dark[0].data
                        else: dark_data = 0
                        if self.__mflat!=None: flat_data = flat[0].data/mode  # normalization
                        else: flat_time = 1     
                        sci_data = f[0].data 
                                                               
                    # STEP 2.1: Check EXPTIME and apply master DARK and master FLAT
                    if time_scale != 1.0: log.debug("Scaling master dark ...")
                    sci_data = (sci_data - dark_data*time_scale) / flat_data
                    #store back the new values
                    if n_ext > 1: # it means, MEF
                        f[chip+1].data = sci_data
                    else:
                        f[0].data = sci_data        
                if self.__mdark != None: 
                    f[0].header.add_history('Dark subtracted %s' %self.__mdark)
                if self.__mflat != None: 
                    f[0].header.add_history('Flat-Field with %s' %self.__mflat)        
                            
                # Write output to outframe (data object actually still points 
                # to input data)
                try:
                    f[0].scale('float32')
                    f.writeto(newpathname, output_verify='ignore')
                except IOError:
                    raise ExError('Cannot write output to %s' % newpathname)
                     
                f.close()
                result_file_list.append(newpathname)            
                log.debug('Saved new dark subtracted or/and flattened file  %s', \
                          newpathname )
        
        log.debug(t.tac() )
        log.info("Successful end of applyDarkFlat !")
                
        return result_file_list
        
        
        
        
################################################################################
# Functions       
def usage ():
    """ Usage function to be shown as help for the user """
    print ''
    print 'NAME'
    print '       applyDarkFlat.py - Apply Dark, Flat or both\n'
    print 'SYNOPSIS'
    print '       applyDarkFlat.py [options]\n'
    print 'DESCRIPTION'
    print '       Subtract a master dark file or divide by and a'
    print '       master flat field (or both) a list of science files.' 
    print '       Input Flat will be normalized before dividing by, so it must NOT be normalized'
    print ' '
    print 'OPTIONS'
    print "-s / --source=      Source file list of data frames"
    print "-d / --dark=        Master dark frame to subtract (optional)"
    print "-f / --flat=        Master (NOT normalized!) flat field to divide by (optional)"
    print "-v                  Verbose debugging output\n"
    print 'VERSION'
    print '       15 Nov 2010 '
    print ''
    raise SystemExit

################################################################################
# main
if __name__ == "__main__":
    
    # Get and check command-line options
    args = sys.argv[1:]
    source_file_list = None
    dark_file = None
    flat_file = None
    out_dir = "/tmp/"
    debug = False
    
    try:
        opts, args = getopt.getopt(args, "s:d:f:o:v", ['source=', 'dark=', \
                                                       'flat=', 'out_dir='])
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
            debug = True
            
    if source_file_list == None or (flat_file == None and dark_file == None):
        usage()
        sys.exit(3)
    
    filelist = [line.replace( "\n", "") for line in fileinput.input(source_file_list)]
    #filelist = ['/disk-a/caha/panic/DATA/ALHAMBRA_1/A0408060036.fits', 
    #'/disk-a/caha/panic/DATA/ALHAMBRA_1/A0408060037.fits']
    print "Files:", filelist
    res = ApplyDarkFlat(filelist, dark_file, flat_file, out_dir)
    res.apply()
    
