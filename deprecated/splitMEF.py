#!/usr/bin/env python
################################################################################
#
#
# PAPI (PANIC PIpeline)
#
# splitMEF.py
#
# Split a Multi Extension FITS file into N single FITS frames
#
# Created    : 10/09/2010    jmiguel@iaa.es -
# Last update: 
# TODO
#       - An alternative way is to use iraf.imcopy !
#       - Other alternative is iraf.mscsplit
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


# Logging
from misc.paLog import log



class SplitMEF:
    """
    \brief Class used to split a MEF into single FITS frames, coping all the 
    header information required
    
    \par Class:
        SpliMEF
    \par Purpose:
         Split MEF files into simple FITS files
    \par Description:
         All PANIC observational data will be recorded in the so-called 
         multi-extension FITS format. A MEF file is comprised of several 
         segments called Header/Data Units (HDUs). Every HDU consists of 
         an Header Unit (the well known FITS headers) in ASCII format followed 
         by an optional Data Unit. The first HDU is called the primary, and any 
         number of additional HDUs may follow. These additional HDUs are 
         referred to as FITS extensions.

         In the PANIC FITS, the primary HDU only contains ASCII header cards 
         describing the observation, but no data. The astronomical data arrays 
         are stored in additional image extensions. There are 4 image extensions 
         , 1 for each detector in the mosaic.
        
    \par Language:
        Python, PyFITS
    \param file_list
        A list FITS files or directory
    \param output_filename
        File where single FITS will be created
    \retval 0
        If no error
    \author
        JMIbannez, IAA-CSIC
        
    """
    def __init__(self,  input_files,  output_filename_suffix = ".%02d.fits",\
                  copy_keyword = ['DATE','OBJECT','DATE-OBS','RA','DEC',\
                                'EQUINOX','RADECSYS','UTC','LST','UT','ST',\
                                'AIRMASS','IMAGETYP','EXPTIME','TELESCOP',\
                                'INSTRUME','MJD-OBS','FILTER','FILTER1','FILTER2']):
         
        self.input_files = input_files
        self.out_filename_suffix = output_filename_suffix  # suffix 
        # List of FITS keywords to propagate from the primary HDU to the 
        # single FITS header
        self.copy_keyword=copy_keyword 
            
    def run(self):
      
        """
        \brief Run the splitting
        """   
        log.info("Starting SplitMEF")
         
        for file in self.input_files:        
            try:
                hdulist = pyfits.open(file)
            except IOError:
                print 'Error, can not open file %s' %(file)
                raise Exception("Error, can not open file %s"%file)
            
            try:
                if hdulist[0].header['EXTEND']!=True:
                    print 'Error, file %s is not a MEF file' %(file)
                    raise Exception("Error, file %s is not a MEF file"%(file))
            except KeyError:
                raise Exception("Error, file %s is not a MEF file"%(file))
            
            try:
                n_ext=hdulist[0].header['NEXTEND']
            except KeyError:
                print 'Warning, card NEXTEND not found. Counting number \
                of extensions...'
                n_ext=0
                while (1):
                    try:
                        if (hdulist[n_ext+1].header['XTENSION'] == 'IMAGE'): 
                            n_ext+= 1
                    except:
                        break
                print "->Found %d extensions" %(n_ext)
                         
            out_filenames=[]
            for i in range(1,n_ext+1):
                suffix=self.out_filename_suffix %i
                out_filenames.append(file.replace('.fits', suffix))
                out_hdulist = pyfits.HDUList([pyfits.PrimaryHDU (header = hdulist[i].header, data=hdulist[i].data)])
                out_hdulist.verify('silentfix')
                # now, copy extra keywords required
                for key in self.copy_keyword:
                    try:
                        value=hdulist[0].header.ascardlist()[key].value
                        comment=hdulist[0].header.ascardlist()[key].comment
                        out_hdulist[0].header.update(key,value,comment)
                    except KeyError:
                        print 'Warning, key %s cannot not be copied,\
                         is not in the header' %(key)
                # delete some keywords not required anymore
                del out_hdulist[0].header['EXTNAME']                
                out_hdulist.writeto(out_filenames[i-1], \
                                    output_verify = 'ignore', clobber = True)
                out_hdulist.close(output_verify = 'ignore')
                del out_hdulist
                print "File %s created " %(out_filenames[i-1])
            
        log.info("End of SplitMEF. %d files created", n_ext)
        return n_ext, out_filenames
            
                           
                      
################################################################################
# main
if __name__ == "__main__":
    
    # Get and check command-line options
    usage = "usage: %prog [options] arg1 arg2"
    parser = OptionParser(usage)
    
    parser.add_option("-f", "--file",
                  action="store", dest="file",
                  help="Input MEF file. It has to be a fullpath file name")
    
    parser.add_option("-l", "--input",
                  action = "store", dest="input_file_list",
                  help = "Source file list of data frames. It has\
                   to be a fullpath file name")
    
    parser.add_option("-s", "--suffix",
                  action = "store", dest = "out_suffix", \
                  help = "suffix to out files (default .%02d.fits)")
    

    (options, args) = parser.parse_args()
    
    if options.file:
        filelist = [options.file]
    elif options.input_file_list:
        filelist = [line.replace ("\n", "") for line in \
                    fileinput.input(options.input_file_list)]
    else:
        parser.print_help()
        parser.error("incorrect number of arguments " )
    if not options.out_suffix:
        options.out_suffix = ".%02d.fits"
        
    copy_keyword = ['DATE', 'OBJECT', 'DATE-OBS', 'RA', 'DEC', 'EQUINOX',\
                     'RADECSYS', 'UTC', 'LST', 'UT', 'ST', 'AIRMASS',\
                     'IMAGETYP', 'EXPTIME', 'TELESCOP', 'INSTRUME', 'MJD-OBS',\
                     'FILTER','FILTER2']    
    split = SplitMEF(filelist, options.out_suffix, copy_keyword)
    split.run()
        