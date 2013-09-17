#!/usr/bin/env python
################################################################################
#
#
# PAPI (PANIC PIpeline)
#
# modFITS.py
#
# Modify some keyword value
#
# Created    : 25/11/2010    jmiguel@iaa.es -
# Last update: 28/05/2013    jmiguel@iaa.es - 
#             Adapted to work with old version of PyFITS (i.e., arch.osn.iaa.es)
#              17/09/2013    jmiguel@iaa.es
#                            New name and now keywords can be added.
# TODO
#     - Give a map of keyword-value as input
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



# Interact with FITS files
import pyfits

def modFITS(files, keyword, value, ext=0):
    """
    Method used to modify a keyword in a list of FITS files.
    
    Parameters
    ----------
    files: sequence 
        A list FITS files
        
    keyword: str
        Header keyword to be modified. If not exists, it will be added. 
    
    value: str
        Value to be assigned to the keyword
         
    ext: int
        Optionally, extension number where header keyword will be looked for.
        By default, Primary HDU (ext=0) is used.

    Returns
    -------
        Number of files modified
        
    """
           
    print "Starting modFITS..."
    n = 0

    pyfitsversion = pyfits.__version__.split(".")
    pyfitsversion = [ a.split('-')[0] for a in pyfitsversion ]
    if map(int, pyfitsversion) < map(int, ['2','4']):
        old_version = True
        print "Warning, old version of PyFITS found. BITPIX could be modified !"
    else:
        old_version = False
        
    for file in files:        
        try:
            # To preserve image scale (BITPIX)--> do_not_scale_image_data 
            # (http://goo.gl/zYkc6)
            # Other option, is use pyfits.ImageHDU.scale_back
            if not old_version:
                hdulist = pyfits.open(file, "update", do_not_scale_image_data=True)
            else:
                hdulist = pyfits.open(file, "update")
                
        except IOError:
            print 'Error, can not open file %s' %(file)
            continue

        #Check if it is a MEF file 
        if ext>(len(hdulist)-1):
            print "[Error] Wrong Extension number for file: %s"%file
            continue
        
        try:
            if not keyword in hdulist[ext].header:
                print "Warning, keyword does not exists; it will be added."
                
            if old_version:
                # Store scaling information for later use, as this is 
                # deleted when the array is updated
                BZERO = hdulist[ext].header['BZERO']
                BSCALE = hdulist[ext].header['BSCALE']
                BITPIX = hdulist[ext].header['BITPIX']
                bitpix_designation = pyfits.ImageHDU.NumCode[BITPIX]

            hdulist[ext].header.update(keyword, value)
            if old_version:
                hdulist[ext].scale(bitpix_designation, 'old')
                hdulist[ext].header.update('BZERO', BZERO)
                hdulist[ext].header.update('BSCALE', BSCALE)
                hdulist[ext].header.update('BITPIX', BITPIX)

            hdulist.writeto(file, clobber=True, output_verify='ignore')
            n+=1
        except Exception,e:
            print "[Error] Cannot modify keyword %s: \n %s"%(keyword, str(e))
        hdulist.close()    
    
    print "End of modFITS. %d files modified"%n
    return n
    
                                  
################################################################################
# main
if __name__ == "__main__":
    
    # Get and check command-line options
    usage = "usage: %prog [options] arg1 arg2"
    parser = OptionParser(usage)
    
    parser.add_option("-f", "--fits",
                  action="store", dest="fits", type="str",
                  help="Input FITS file. It has to be a fullpath file name")
    
    parser.add_option("-l", "--input",
                  action="store", dest="input_file_list", type="str",
                  help="Source file list of data frames. It has to be a fullpath file name")
                  
    parser.add_option("-k", "--key_value",
                  action="store", dest="keyword", type="str", nargs=1,
                  help="Keyword space separated to be modified")
                  
    parser.add_option("-v", "--value", type="str",
                  action="store", dest="value",
                  help="Value to set to 'keyword'")
                                
    parser.add_option("-e", "--ext",
                  action="store", dest="extension_number", type="int", default=0,
                  help="Extension number in which to look for 'keyword' [0,N]")
    
    (options, args) = parser.parse_args()
    
    if options.fits:
        filelist = [options.fits]
    elif options.input_file_list:
        filelist = [line.replace( "\n", "") for line in fileinput.input(options.input_file_list)]
    else:
        parser.print_help()
        parser.error("incorrect number of arguments " )
    
    if not options.keyword or not options.value:
        parser.print_help()
        parser.error("incorrect number of arguments " )
        
    try:
        modFITS(filelist, options.keyword, options.value, options.extension_number)    
    except Exception, e:
        raise e
        
