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

def doAstrometry( input_image, output_image=None, catalog='2MASS'):
    """ Do the astrometric calibration to the input image
      
    The method does astrometry on the image entirely using Emmanuel Bertin's
    tools, namely SExtractor, SCAMP and SWarp. First, SExtractor is run on the
    image, saving the output catalog in the FITS_LDAC binary format that is then
    read by SCAMP, whose output is in turn a '.head' file. This file contains
    the updated astrometric information, which is finally merged with the 
    original image using the SWarp stacking tool. In case no output path is
    specified, the resulting, calibrated image is saved to a temporary
    directory, most likely /tmp.

    
    Note that the method does not receive any other parameters. How SExtractor,
    SCAMP and SWarp work is determined by their respective configuration files,
    located in the config file specified or with default values, and they are what 
    should be modified in case you need to adjust their operation so that they better suit
    your needs.

    Keyword arguments:
  
    input_image  - path to the FITS image whose astrometry is to be done.
    output_image - path to which the astrometically calibrated version of the
                  image. If no path is provided, it will be saved to a temporary
                  location.
    catalog      - external catalog to be used for the astrometric calibration. 
                   By default, 2MASS is used 

    """
    
    # Save the resulting image to a temporary file if no path was specified
    if output_image is None:
        output_fd, output_path = tempfile.mkstemp(suffix='.fits')
        os.close(output_fd)
    log.debug("*** Start Astrometric calibration ***")
    
    ## STEP 0: Run IRDR::initwcs to initialize rough WCS header, thus modify the file headers
    # initwcs also converts to J2000.0 EQUINOX
    log.debug("***Doing WCS-header initialization ...")
    initwcs_path=os.environ['PAPI_HOME']+'/irdr/bin/initwcs'
    args = [initwcs_path, input_image]
    print "ARGS=", args
    ret_code = subprocess.call(args)
    if ret_code!=0:
        raise RuntimeError("There was an error while running 'initwcs'")
            
    ## STEP 1: Create SExtractor catalog (.ldac)
    log.debug("*** Creating SExtractor catalog ....")
    sex = astromatic.SExtractor()
    #sex.config['CONFIG_FILE']="/disk-a/caha/panic/DEVELOP/PIPELINE/PANIC/trunk/config_files/sex.conf"
    #sex.ext_config['CHECKIMAGE_TYPE'] = "OBJECTS"
    sex.config['CATALOG_TYPE'] = "FITS_LDAC"
    sex.config['CATALOG_NAME'] = input_image + ".ldac"
    try:
        sex.run(input_image, updateconfig=True, clean=False)
    except: 
        raise
                    
    ## STEP 2: Make astrometric calibration 
    log.debug("Doing astrometric calibration....")
    scamp = astromatic.SCAMP()
    scamp.config['CONFIG_FILE']="/disk-a/caha/panic/DEVELOP/PIPELINE/PANIC/trunk/config_files/scamp.conf"
    scamp.ext_config['ASTREF_CATALOG']=catalog
    scamp.ext_config['SOLVE_PHOTOM']="N"
    scamp.ext_config['WRITE_XML']="N"
    cat_file = input_image.replace( ".fits", ".fits.ldac")
    #updateconfig=False means scamp will use the specified config file instead of the single config parameters
    #but, "ext_config" parameters will be used in any case
    try:
        scamp.run(cat_file, updateconfig=False, clean=False)
    except:
        raise
    
    ## STEP 3: Merge the astrometric parameters (.head keywords) with SWARP, and using .head files created by SCAMP
    log.debug("Merging astrometric calibration parameters and resampling ...")
    swarp = astromatic.SWARP()
    swarp.config['CONFIG_FILE']="/disk-a/caha/panic/DEVELOP/PIPELINE/PANIC/trunk/config_files/swarp.conf"
    swarp.ext_config['IMAGEOUT_NAME']=output_image
    try:
        swarp.run(input_image, updateconfig=False, clean=False)
    except:
        raise
    
    return output_image 

class AstroWarp(object):
    """ Astrometric warping """

    def __init__(self, input_files, catalog='2MASS', coadded_file="/tmp/astrowarp.fits"):
        """ Instantiation method for AstroWarp class

        Keyword arguments:
        input_files  - the set of overlapping reduced frames, i.e., dark subtracted, flatted and sky subtracted 
        catalog      - the catalog to use for the astrometric calibration (by default, 2MASS)
        coadded_file - the output final coadded image
        
        """
        self.input_files = input_files
        self.catalog = catalog
        self.coadded_file = coadded_file
        
    def run(self):
        """ Start the computing of the coadded image, following the next steps:
      
            0. Initialize rought WCS header
            1. Call up SExtractor to create pixel object list (ldac catalog)
            2. Call up SCAMP to make the astrometric calibration of the overlapped set of frames (.head calibration)
            3. Call up SWARP to make the coadd, using the .head files generated by SCAMP and containing the distortion model parameters
            4. Make the final astrometric calibration the the coadded frame (SExtractor+SCAMP)
            
        """
        
        log.debug("*** Start Astrowarp ***")
        ## STEP 0: Run IRDR::initwcs to initialize rough WCS header, thus modify the file headers
        # initwcs also converts to J2000.0 EQUINOX
        log.debug("***Doing WCS-header initialization ...")
        initwcs_path=os.environ['PAPI_HOME']+'/irdr/bin/initwcs'
        for file in self.input_files:
            args = [initwcs_path, file]
            print "ARGS=", args
            ret_code = subprocess.call(args)
            if ret_code!=0:
                raise RuntimeError("There was an error while running 'initwcs'")
               
        ## STEP 1: Create SExtractor catalogs (.ldac)
        log.debug("*** Creating SExtractor catalog ....")
        for file in self.input_files:
            sex = astromatic.SExtractor()
            #sex.config['CONFIG_FILE']="/disk-a/caha/panic/DEVELOP/PIPELINE/PANIC/trunk/config_files/sex.conf"
            #sex.ext_config['CHECKIMAGE_TYPE'] = "OBJECTS"
            sex.config['CATALOG_TYPE'] = "FITS_LDAC"
            sex.config['CATALOG_NAME'] = file + ".ldac"
            sex.run(file, updateconfig=True, clean=False)
                        
        ## STEP 2: Make the multi-astrometric calibration for each file (all overlapped-files together)
        log.debug("*** Doing multi-astrometric calibration ....")
        scamp = astromatic.SCAMP()
        scamp.config['CONFIG_FILE']="/disk-a/caha/panic/DEVELOP/PIPELINE/PANIC/trunk/config_files/scamp.conf"
        scamp.ext_config['ASTREF_CATALOG']=self.catalog
        scamp.ext_config['SOLVE_PHOTOM']="N"
        cat_files = [f.replace( ".fits", ".fits.ldac") for f in self.input_files]
        #updateconfig=False means scamp will use the specified config file instead of the single config parameters
        scamp.run(cat_files, updateconfig=False, clean=False)
        
        
        ## STEP 3: Make the coadding with SWARP, and using .head files created by SCAMP
        # It requires the files are overlapped, i.e., have an common sky-area
        log.debug("*** Coadding overlapped files....")
        swarp = astromatic.SWARP()
        swarp.config['CONFIG_FILE']="/disk-a/caha/panic/DEVELOP/PIPELINE/PANIC/trunk/config_files/swarp.conf"
        swarp.ext_config['COPY_KEYWORDS']='OBJECT,INSTRUME,TELESCOPE,IMAGETYP,FILTER,FILTER2,SCALE,MJD-OBS'
        swarp.ext_config['IMAGEOUT_NAME']=os.getcwd() + "/coadd_tmp.fits"
        swarp.run(self.input_files, updateconfig=False, clean=False)
        
        ## STEP 4: Make again the final astrometric calibration to the final coadd
        log.debug("*** Doing final astrometril calibration....")
        doAstrometry(os.getcwd() + "/coadd_tmp.fits", self.coadded_file, self.catalog)
        
        
        
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
    
    