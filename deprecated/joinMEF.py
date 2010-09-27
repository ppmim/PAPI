#!/usr/bin/env python
################################################################################
#
#
# PAPI (PANIC PIpeline)
#
# joinMEF.py
#
# Join a Multi Extension FITS file into 1 single FITS frame
#
# Created    : 13/09/2010    jmiguel@iaa.es -
# Last update: 
# TODO
#       -An alternative way is to use SWARP to create a single FITS frame ( swarp file_with_ext.fits )
#       However, this alternative takes longer.
#
#       - Other alternative is iraf.mscjoin
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
import numpy as np
import datahandler


# Logging
from misc.paLog import log



class JoinMEF:
    """
    \brief Class used to join a MEF into a single FITS frames, coping all the header information required
    
    \par Class:
        JoinMEF
    \par Purpose:
         Join MEF files into a single FITS frame
    \par Description:
         All PANIC observational data will be recorded in the so-called multi-extension FITS format. A MEF file is comprised of several segments called Header/Data Units (HDUs). Every HDU consists of an Header Unit (the well known FITS headers) in ASCII format followed by an optional Data Unit. The first HDU is called the primary, and any number of additional HDUs may follow. These additional HDUs are referred to as FITS extensions.

         In the PANIC FITS, the primary HDU only contains ASCII header cards describing the observation, but no data. The astronomical data arrays are stored in additional image extensions. There are 4 image extensions , 1 for each detector in the mosaic.
        
    \par Language:
        Python, PyFITS
    \param input_files
        A list FITS files or directory
    \param output_filename_suffix
        Suffix to add to the outfile
    \retval 0
        If no error
    \author
        JMIbannez, IAA-CSIC
        
    """
    def __init__(self,  input_files,  output_filename_suffix=".join.fits"):
         
        self.input_files = input_files
        self.out_filename_suffix = output_filename_suffix  # suffix 
            
    def run(self):
      
        """
        \brief Run the join of MEF 
        """   
        log.info("Starting JoinMEF")
         
        for file in self.input_files:        
            try:
                hdulist = pyfits.open(file)
            except IOError:
                print 'Error, can not open file %s' %(file)
                return 0
            
            try:
                if hdulist[0].header['EXTEND']!=True:
                    print 'Error, file %s is not a MEF file' %(file)
                    return 0
            except KeyError:
                print 'Error, file %s is not a MEF file' %(file)
                return 0
            
            try:
                next=hdulist[0].header['NEXTEND']
            except KeyError:
                print 'Warning, card NEXTEND not found. Counting number of extensions...'
                next=0
                while (1):
                    try:
                        if (hdulist[next+1].header['XTENSION'] == 'IMAGE'): next += 1
                    except:
                        break
            print "->Found %d extensions" %(next)
                         
            out_filename=file.replace(".fits",self.out_filename_suffix)
            width=naxis1=2048
            height=naxis2=2048
            temp12=np.zeros((height,width*2), dtype=np.float32)
            temp34=np.zeros((height,width*2), dtype=np.float32)
            for i in range(0, height):
                # Q1 i-row
                temp12[i,0 : width]=hdulist[1].data[i, 0 : width]
                # Q2 i-row
                temp12[i, width: 2*width]=hdulist[2].data[i, 0 : width]
                # Q3 i-row
                temp34[i, 0 : width]=hdulist[3].data[i, 0 : width]
                # Q4 i-row
                temp34[i, width : 2*width]=hdulist[4].data[i, 0 : width]

            joined_data = np.append(temp12, temp34).reshape(4096,4096)
            hdu = pyfits.HDUList([pyfits.PrimaryHDU(header=hdulist[0].header, data=joined_data)])
            #hdu.verify('silentfix')
            # now, copy extra keywords required
            try:
                hdu[0].header.update("BITPIX",-32)
                hdu[0].header.update("NAXIS1",4096)
                hdu[0].header.update("NAXIS2",4096)
            except KeyError:
                print 'Warning, some key cannot not be copied'
                            
            hdu.writeto(out_filename, output_verify='ignore', clobber=True)
            hdu.close(output_verify='ignore')
            del hdu
            print "File %s created " %(out_filename)
        
        log.info("End of JoinMEF. %d files created", next)
        return next 
            
                           
                      
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
                  action="store", dest="input_file_list",
                  help="Source file list of data frames. It has to be a fullpath file name")
    
    parser.add_option("-s", "--suffix",
                  action="store", dest="out_suffix", help="suffix to out files (default .%02d.fits)")
    

    (options, args) = parser.parse_args()
    
    if options.file:
        filelist=[options.file]
    elif options.input_file_list:
        filelist=[line.replace( "\n", "") for line in fileinput.input(options.input_file_list)]
        print filelist
    else:
        parser.print_help()
        parser.error("incorrect number of arguments " )
    if not options.out_suffix:
        options.out_suffix=".join.fits"
        
    join = JoinMEF(filelist, options.out_suffix)
    join.run()
          
        
        