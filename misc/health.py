#!/usr/bin/env python

# Copyright (c) 2012 IAA-CSIC  - All rights reserved. 
# Author: Jose M. Ibanez. 
# Instituto de Astrofisica de Andalucia, IAA-CSIC
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
# PAPI (PANIC PIpeline)
#
# health.py
#
# Tools used:
#
#    IRAF
#
# Created    : 29/10/2012    jmiguel@iaa.es -
# Last update: 
# TODO
#       
################################################################################

# Import necessary modules
from optparse import OptionParser
import sys
import os
import math

import atpy 
import numpy

import matplotlib

import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.cm as cm
import pylab
            

# Logging
from misc.paLog import log
import misc.utils
import datahandler


def run_health_check ( input_file, window='full-frame', out_filename="/tmp/hc_out.fits" ):
    """
    @summary: 
    Takes a input catalog (ascii file) listing all the files to be used in the
    check-health analysis and performs the computation required for it. 
    
    @note: It is based on JWF MIDAS routine. 
    
    @param cat1: Text file listing the files to use for checking 
    @param start: File number where start the packet 
    @param packet_size: Size of the packet size
    @param window: max. error for finding objects within (arcseconds)
    @param out_filename: filename where results will be saved
    @return: filename where results where saved
    """
    
    
    # Read the file list from input_file
    if os.path.exists(input_file)
        filelist = [line.replace( "\n", "") for line in fileinput.input(input_file)]
    else:
        raise Exception("Input file %s does not exist"%input_file)
    
    if window=='Q1':
        area = [10, 10, 2030, 2030]
    elif window=='Q2':
        area = [2100, 10, 4080, 2100]    
    elif window=='Q3':
        area = [2100, 2100, 4080, 4080]    
    elif window=='Q4':
        area = [10, 2100, 2100, 4080]
    elif window=='central':
        area = [512, 512, 3584, 3584]
    else: # or 'full-frame'
        area = [10,10, 4080, 4080] 
    
    return 0



################################################################################
# main
################################################################################
if __name__ == "__main__":

    log.debug( 'Health-Check routines for PANIC')
    
    # Get and check command-line options
        
    usage = "usage: %prog [options] arg1 arg2 ..."
    parser = OptionParser(usage)
    
    parser.add_option("-i", "--input_images",
                  action="store", dest="input_images", 
                  help="input image list (darks and flats")
    
    parser.add_option("-w", "--window",
                  action="store", dest="window", type=str,
                  help="Window to use for computations (default = %default)",
                  default='full-detector')
    
    parser.add_option("-s", "--start_packet",
                  action="store", dest="star_packet", type=int, default=1,
                  help="File where start the packet (default=%default)")
    
    parser.add_option("-p", "--packet-size",
                  action="store", dest="packet_size", type=int, default=10,
                  help="Number of files per packet (default=%default)")
                  
    parser.add_option("-o", "--output",
                  action="store", dest="output_filename", 
                  help="output plot filename (default = %default)",
                  default="health_check.pdf")
    
    parser.add_option("-v", "--verbose",
                  action="store_true", dest="verbose", default=True,
                  help="verbose mode [default=%default]")
    
                                
    (options, args) = parser.parse_args()
    
    
    if not options.input_images or len(args)!=0: 
    # args is the leftover positional arguments after all options have been processed
        parser.print_help()
        parser.error("Wrong number of arguments")

    if not os.path.exists(options.input_images):
        log.error ("Input image %s does not exist", options.input_images)
        sys.exit(0)
        
    try:
        run_health_check(options.input_image, options.start, options.packet_size, 
                     options.window, options.output_file)
    except Exception,e:
        log.info("Some error while running Health-Check routine: %s"%str(e))
        sys.exit(0)
        
