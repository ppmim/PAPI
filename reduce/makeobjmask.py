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
#              16/12/2010    jmiguel@iaa.es - Added single point object mask feature
#                                             (even for MEF files)
################################################################################
#
# Creates object masks (SExtractor OBJECTS images) for a list of FITS images.
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
import astropy.io.fits as fits

from misc.paLog import log
import astromatic.sextractor
import astromatic.ldac
import datahandler


#-----------------------------------------------------------------------
def makeObjMask (inputfile, minarea=5, maxarea=0,  threshold=2.0, 
                 saturlevel=300000, outputfile="/tmp/out.txt", 
                 single_point=False):
    """
    DESCRIPTION
      Create an object mask of the inputfile/s based on SExtractor
       
    INPUTS
      inputfile    it can be a regular expresion or a filename 
                          containing a filelist
            
    OPTIONAL INPUTS
            
      outputfile   filename of file list with the object file/s 
                   created by SExtractor
            
      minarea      SExtractor DETECT_MINAREA (min. # of pixels above threshold)
      
      maxarea      SExtractor DETECT_MAXAREA (max. # of pixels above threshold, 0=unlimited)
                       
      threshold    SExtractor DETECT_THRESH
            
      saturlevel   Pixel Saturation level 
            
      single_point  If true, means the image will be reduced to a 
                    single point object mask,i.e., a single pixel set to 1
                    for each detected object.
            
    OUTPUTS
      n             Number of object mask created        
    """
       
    log.debug("Start of [makeObjMask]")
    
    files = []       
    # Check if inputfile is FITS file 
    if datahandler.isaFITS(inputfile)==True:
        files=[inputfile]
    # or a text file having the list of files to be masked
    elif os.path.isfile(inputfile):
        files=[line.replace( "\n", "") for line in fileinput.input(inputfile)]
    # or must be a regular expresion
    else:
        files = glob.glob(inputfile)
        files.sort()
    
    f_out = open(outputfile,"w")
    
    sex = astromatic.SExtractor()
    n = 0
    for fn in files:
        if not os.path.exists(fn):      # check whether input file exists
            log.error( 'File %s does not exist', fn)
            raise Exception ("File %s does not exist"%fn) 
     
        log.debug("*** Creating SExtractor object mask for file %s....", fn)
        #sex.config['CONFIG_FILE']=sex_config
        sex.config['CATALOG_TYPE'] = "FITS_LDAC"
        sex.config['CATALOG_NAME'] = fn + ".ldac"
        sex.config['DETECT_MINAREA'] = minarea
        sex.config['DETECT_MAXAREA'] = maxarea
        sex.config['DETECT_THRESH'] = threshold
        sex.config['CHECKIMAGE_TYPE'] = "OBJECTS"
        sex.config['CHECKIMAGE_NAME'] = fn + ".objs"
        sex.config['SATUR_LEVEL'] = saturlevel
        if os.path.exists(fn.replace(".fits",".weight.fits")):
            sex.config['WEIGHT_TYPE']="MAP_WEIGHT"
            sex.config['WEIGHT_IMAGE']=fn.replace(".fits",".weight.fits")
        
        # Run SExtractor     
        try:
            sex.run(fn, updateconfig=True, clean=False)
        except Exception,e: 
            log.debug("Some error while running SExtractor : %s", str(e))
            raise Exception("Some error while running SExtractor : %s"%str(e))          
        
        # Check an output file was generated, otherwise an error happened !
        if not os.path.exists(fn+".objs"):
            log.error("Some error while running SExtractor, no object mask file found and expected %s"%(fn+".objs"))
            raise Exception("Some error while running SExtractor, no object mask file found")

        # Reduce the object mask to a single point mask, in which each object
        # is represented by a single, one-valued pixel, located at the
        # coordinates specified by its X_IMAGE and Y_IMAGE parameters in the
        # SExtractor catalog. The remaining pixels in the image
        # will be set to zero. Note that, therefore, the 'single-point' mask will
        # have as many non-zero pixels as objects are in the SExtractor catalog.
        # NOTE:This feature is used to compute the dither offsets, but not for 
        # while object masking in skysubtraction
        if single_point==True:
            log.debug("Single point mask reduction")
            # NOTE we update/overwrite the image and don't create a new one
            myfits = fits.open(fn+".objs", mode="update")
            if len(myfits)>1: # is a MEF file
                next = len(myfits)-1
            else: 
                next = 1
            for ext in range(next):
                if next==1: data = myfits[0].data
                else: data = myfits[ext+1].data
                data[:] = 0 # set to 0 all pixels
                x_size = len(data[0])
                y_size = len(data)
                #stars = read_stars(fn + ".ldac")
                try:
                    cat = astromatic.ldac.openObjectFile(fn+".ldac", 
                                                         table='LDAC_OBJECTS')
                    if len(cat)<=0:
                      log.warning("No object found in catalog %s"%(fn+".ldac")) 
                      continue
                    for star in cat:
                        if round(star['X_IMAGE'])<x_size and round(star['Y_IMAGE'])<y_size:
                            data[round(star['Y_IMAGE']),round(star['X_IMAGE'])] = 1 # Note: be careful with X,Y coordinates position
                except Exception,e:
                    print "Y_IMAGE=",star['Y_IMAGE']
                    print "X_IMAGE=",star['X_IMAGE']
                    myfits.close(output_verify='ignore')
                    raise Exception("Error while creating single point object mask :%s"%str(e))
                
            myfits.close(output_verify='ignore')
            log.debug("Object mask (single_point) file created for file : %s"%fn)
        else:
            log.debug("Object mask file created for file : %s",fn)    

        n+=1
        log.debug("Adding file: %s"%fn)
        f_out.write(fn+".objs"+"\n")
            
    sex.clean(config=False, catalog=True, check=False)

    f_out.close()
    log.debug("Successful ending of makeObjMask => %d object mask files created", n)
    return n
    
    
################################################################################
# main
if __name__ == "__main__":
    
    usage = "usage: %prog [options]"
    desc = """ Creates object masks (SExtractor OBJECTS images) for a list of FITS images.
Expects the command "sex" (SExtractor Version 2+) in path.  If weight maps
exist they will be used (assume weight map filename given by replacing .fits
with .weight.fits)."""

    parser = OptionParser(usage, description=desc)
    
                  
    parser.add_option("-s", "--file", type="str",
                  action="store", dest="inputfile",
                  help="It can be a source file listing data frames or a single FITS file to process.")
    
    parser.add_option("-o", "--output", type="str",
                  action="store", dest="outputfile", 
                  help="Output text file including the list of objects mask files created by SExtractor ending with '.objs' suffix")
    
    parser.add_option("-m", "--minarea", type="int", default=5,
                  action="store", dest="minarea", 
                  help="SExtractor DETECT_MINAREA (default=%default)")
                  
    parser.add_option("-t", "--threshold", type="float", default=2.0,
                  action="store", dest="threshold", 
                  help="SExtractor DETECT_THRESH (default=%default)")
    
    parser.add_option("-l", "--saturlevel", type="int", default=300000,
                  action="store", dest="saturlevel",
                  help="SExtractor SATUR_LEVEL (default=%default)")
    
    parser.add_option("-1", "--single_point", default=False,
                  action="store_true", dest="single_point",
                  help="Create a single point object mask (default=%default)")
    
    
    (options, args) = parser.parse_args()

    if len(sys.argv[1:])<1:
       parser.print_help()
       sys.exit(0)
    
    # args is the leftover positional arguments after all options have been 
    # processed
    if not options.inputfile or not options.outputfile or len(args)!=0: 
        parser.print_help()
        parser.error("Incorrect number of arguments " )
        
    try:
        makeObjMask( options.inputfile, options.minarea, options.threshold, 
                     options.saturlevel, options.outputfile, 
                     options.single_point)
    except Exception, e:
        log.error("Error while running 'makeobjectmask' %s"%str(e))
    else:
        log.debug("Well done!")
        