#! /usr/bin/env python

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
# PAPI (PANIC PIpeline)
#
# makeobjmask.py
#
# Created    : 04/08/2008    jmiguel@iaa.es
# Last update: 04/08/2009    jmiguel@iaa.es
#              26/10/2010    jmiguel@iaa.es - Added new astromatic.sextractor
#                                             support
################################################################################
#
# Create object masks (SExtractor OBJECTS images) for a list of FITS images.
# Expects the command "sex" (SExtractor Version 2+) in path.  If weight maps
# exist they will be used (assume weight map filename given by replacing .fits
# with .weight.fits).
#

################################################################################
# Import necessary modules

import sys
import os
import glob
import fileinput
from optparse import OptionParser

import misc.utils as utils
from misc.paLog import log
import astromatic.sextractor


#-----------------------------------------------------------------------
def makeObjMask (inputfile, minarea=5,  threshold=2.0, outputfile="/tmp/out.txt"):
    """DESCRIPTION
                Create an object mask of the inputfile/s based on SExtractor
           
           INPUTS
                inputfile       it can be a regular expresion or a filename containing a filelist
                
           OPTIONAL INPUTS
                
                outputfile      filename of file list with the object file/s created by SExtractor
                
                p_min_area      SExtractor DETECT_MINAREA (minimun object area)
                           
                p_mask_thresh   SExtractor DETECT_THRESH
                
                
           OUTPUTS
                outputfile      Filepath containig the list of objects mask files created by SExtractor ending with '.objs' suffix
              
      """
         
    # Some pathname settings and check
    irdr_basedir=''
    try:
        irdr_basedir=os.environ['IRDR_BASEDIR']
    except KeyError:
        log.error("Please, setenv IRDR_BASEDIR")
        sys.exit(0)
        
    sex_config=irdr_basedir+"/src/config/default.sex"        
    if not os.path.exists(sex_config):      # check whether input file exists
        log.error( 'File %s does not exist', sex_config)
        sys.exit(1) 
           
    files = []       
    # Check if inputfile is a filename of a file list OR a regular expresion        
    if os.path.isfile(inputfile):
        files=[line.replace( "\n", "") for line in fileinput.input(inputfile)]
    else:
        files = glob.glob(inputfile)
        files.sort()
        
    f_out = open(outputfile,"w")
    
    for fn in files:
        if not os.path.exists(fn):      # check whether input file exists
            log.error( 'File %s does not exist', fn)
            raise Exception ("File %s does not exist",fn) 
     
        log.debug("*** Creating SExtractor object mask ....")
        sex = astromatic.SExtractor()
        #sex.config['CONFIG_FILE']=sex_config
        sex.config['CATALOG_TYPE'] = "FITS_LDAC"
        sex.config['CATALOG_NAME'] = fn + ".ldac"
        sex.config['DETECT_MINAREA'] = minarea
        sex.config['DETECT_THRESH'] = threshold
        sex.config['CHECKIMAGE_TYPE'] = "OBJECTS"
        sex.config['CHECKIMAGE_NAME'] = fn + ".objs"
        sex.config['SATUR_LEVEL'] = 300000
        if os.path.exists(fn.replace(".fits",".weight.fits")):
            sex.config['WEIGHT_TYPE']="MAP_WEIGHT"
            sex.config['WEIGHT_IMAGE']=fn.replace(".fits",".weight.fits")
            
        try:
            sex.run(fn, updateconfig=True, clean=False)
        except Exception,e: 
            log.debug("Some error while running SExtractor : %s", str(e))
            raise Exception("Some error while running SExtractor : \n %s",str(e))
        
        log.debug("Object mask file created for file : %s",fn)
        f_out.write(fn+".objs"+"\n")
            
    sex.clean(config=True, catalog=True, check=False)
    f_out.close()
    log.debug("Succesful ending of makeObjMask")
    
    
################################################################################
# main
if __name__ == "__main__":
    print 'Start makeObjMask....'
    
    usage = "usage: %prog [options] arg1 arg2 ..."
    parser = OptionParser(usage)
    
                  
    parser.add_option("-s", "--file",
                  action="store", dest="inputfile",
                  help="Source file list of data frames.")
    
    parser.add_option("-o", "--output",
                  action="store", dest="outputfile", help="output file with the list of catalogs")
    
    parser.add_option("-m", "--minarea",
                  action="store", dest="minarea", default=5, help="SExtractor DETECT_MINAREA (int)")
                  
    parser.add_option("-t", "--threshold",
                  action="store", dest="threshold", default=2, help="SExtractor DETECT_THRESH (int)")
                                   
    parser.add_option("-v", "--verbose",
                  action="store_true", dest="verbose", default=True,
                  help="verbose mode [default]")
    
                                
    (options, args) = parser.parse_args()
    
    if not options.inputfile or not options.outputfile or len(args)!=0: # args is the leftover positional arguments after all options have been processed
        parser.print_help()
        parser.error("Incorrect number of arguments " )
        

    makeObjMask( inputfile, minarea, threshold, outputfile)
    
