#! /usr/bin/env python
#encoding:UTF-8

# Copyright (c) 2010 Jose M. Ibanez All rights reserved.
# Institute of Astrophysics of Andalusia, IAA-CSIC
#
# This file is part of PAPI
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


import sys
import os
import time
import math
import shutil
import tempfile
import subprocess 
from optparse import OptionParser
import fileinput

# PAPI modules
import astromatic 
# Logging
from misc.paLog import log


class AstroWarp(object):
    """ Astrometric warping """

    def __init__(self, input_files, catalog='2MASS', coadded_file="/tmp/astrowarp.fits"):
        """ Instantiation method for AstroWarp class

        Keyword arguments:
        input_files - the set of overlapping reduced frames, i.e., dark, flat and sky subtracted 
        catalog     - the catalog to use for the astrometric calibration (by default, 2MASS)
        coadded_file - the output final coadded image
        
        """
        self.input_files = input_files
        self.catalog = catalog
        self.coadded_file = coadded_file
        
    def run(self):
        """ Start the computing of the coadded image, following the next steps:
      
            0. Initialize rought WCS header
            1. Call to SExtractor to create pixel object list (ldac catalog)
            2. Call to SCAMP to make the astrometric calibration of the overlapped set of frames (.head calibration)
            3. Call to SWARP to make the coadd, using the .head files generated by SCAMP and containing the distortion model parameters
            
        """
        
        # STEP 0: Run IRDR::initwcs to initialize rough WCS header, thus modify the file headers
        # initwcs also converts to J2000.0 EQUINOX
        initwcs_path=os.environ['PAPI_HOME']+'/irdr/bin/initwcs'
        #initwcs_path=os.environ['PAPI_HOME']+'/astrometry_scamp.pl'
        for file in self.input_files:
            args = [initwcs_path, "2mass","noregrid",file]
            args = [initwcs_path, file]
            print "ARGS=", args
            ret_code = subprocess.call(args)
            if ret_code!=0:
                raise RuntimeError("There was an error while running 'initwcs'")
               
        # STEP 1: Create SExtractor catalogs (.ldac)
        for file in self.input_files:
            sex = astromatic.SExtractor()
            #sex.config['CONFIG_FILE']="/disk-a/caha/panic/DEVELOP/PIPELINE/PANIC/trunk/config_files/sex.conf"
            sex.config['CHECKIMAGE_TYPE'] = "NONE"
            sex.config['CATALOG_TYPE'] = "FITS_LDAC"
            sex.config['CATALOG_NAME'] = file+".ldac"
            # Lauch SExtractor on a FITS file
            sex.run(file, updateconfig=True, clean=False)
                        
        # STEP 2: Make astrometric calibration (all overlapped-files together)
        scamp = astromatic.SCAMP()
        scamp.config['CONFIG_FILE']="/disk-a/caha/panic/DEVELOP/PIPELINE/PANIC/trunk/config_files/scamp.conf"
        scamp.config['ASTREF_CATALOG']="2MASS"
        cat_files = [f.replace( ".fits", ".fits.ldac") for f in self.input_files]
        #updateconfig=False means scamp will use the specified config file instead of the single config parameters
        scamp.run(cat_files, updateconfig=False, clean=False)
        
        # STEP 3: Make the coadding with SWARP, and using .head files created by SCAMP
        # It requires the files are overlapped, i.e., have an common sky-area
        swarp = astromatic.SWARP()
        swarp.config['CONFIG_FILE']="/disk-a/caha/panic/DEVELOP/PIPELINE/PANIC/trunk/config_files/swarp.conf"
        swarp.config['IMAGEOUT_NAME']=self.coadded_file
        swarp.config['HEADER_SUFFIX']='.fits.head'
        swarp.run(self.input_files, updateconfig=False, clean=False)
        
        
################################################################################
# main
################################################################################
if __name__ == "__main__":
    log.debug( 'Start AstroWarp....')
    # Get and check command-line options
    
    
    usage = "usage: %prog [options] arg1 arg2 ..."
    parser = OptionParser(usage)
    
                  
    parser.add_option("-s", "--source",
                  action="store", dest="source_file_list",
                  help="Source file list of data frames. It can be a file or directory name.")
    
    parser.add_option("-o", "--output",
                  action="store", dest="output_filename", help="final coadded output image")
    
    parser.add_option("-v", "--verbose",
                  action="store_true", dest="verbose", default=True,
                  help="verbose mode [default]")
    
                                
    (options, args) = parser.parse_args()
    
    if not options.source_file_list or not options.output_filename or len(args)!=0: # args is the leftover positional arguments after all options have been processed
        parser.print_help()
        parser.error("incorrect number of arguments " )
    if options.verbose:
        print "reading %s ..." % options.source_file_list
    
    filelist=[line.replace( "\n", "") for line in fileinput.input(options.source_file_list)]
    astrowarp = AstroWarp(filelist, catalog="2MASS", coadded_file=options.output_filename)
    
    try:
        astrowarp.run()
    except:
        log.error("Some error while running Astrowarp....")
        raise
        sys.exit()
    
    