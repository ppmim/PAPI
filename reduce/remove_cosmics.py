#!/usr/bin/env python

# Copyright (c) 2011-2012 IAA-CSIC  - All rights reserved. 
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
#
# PAPI (PANIC PIpeline)
#
# remove_cosmics.py
#
# Created    : 09/03/2012    jmiguel@iaa.es -
#
# Last update: 
#
# TODO
#   - Include SATURATION_LEVEL and/or other option  
#   - Speed up ! (i.e.,number of iterations, .....)
# NOTE: cosmics module needs scipy module !!
################################################################################
"""
Next code is a wrapper to the cosmics python module of Malte Tewes:

cosmics.py is a small and simple python module to detect and clean cosmic ray 
hits on images (numpy arrays or FITS), using scipy, and based on Pieter van 
Dokkum's L.A.Cosmic algorithm.

L.A.Cosmic = Laplacian cosmic ray detection

U{http://www.astro.yale.edu/dokkum/lacosmic/}

(article : U{http://arxiv.org/abs/astro-ph/0108003})
"""

# Import necessary modules
from optparse import OptionParser
import sys

import astropy.io.fits as fits

# Logging
from misc.paLog import log
from misc.version import __version__

import thirdparty.cosmics.cosmics as cosmics


def remove_cr(in_image, out_image=None, overwrite=False, want_mask=False):
    """
    Remove cosmic rays in O2k or PANIC images
    
    Parameters
    ----------
    in_image : str
        Input filename to remove cosmic rays
    
    out_image : str
        Output filename of cosmics clean image
        
    overwrite: Boolean
        If true, the input file 'in_image' filename will be overwritten,
        otherwise, the 'out_image' filename will be used as output. 
    
    want_mask : Boolean
        If true, the mask with cosmics detected and removed is written into a
        FITS file.
        Otherwise, no mask file is created (default).
    
    Returns
    -------
    If all was successful, the name of the output file is returned
        
    """
    
    if overwrite:
        out_file = in_image
    else:   
        if not out_image:
            out_file = in_image.replace(".fits","_dcr.fits")
        else:
            out_file = out_image
            
    try:
        f_in = fits.open(in_image)
        if len(f_in)==1:
            data_in = f_in[0].data
        else:
            log.error("MEF files currently not supported !")
            raise Exception("MEF files currently not supported !")
            
    except Exception as e:
        log.error("Error opening FITS file : %s"%in_image)
        raise e
    

    try:
        # Read the FITS :
        array, header = cosmics.fromfits(in_image)
        
        # array is a 2D numpy array
        
        # Build the object :
        c = cosmics.cosmicsimage(array, gain=2.2, readnoise=10.0, sigclip = 5.0, 
                                 sigfrac = 0.3, objlim = 5.0)
        # There are other options, check the manual...
        
        # Run the full artillery :
        c.run(maxiter = 4)
        
        # Write the cleaned image into a new FITS file, conserving the original header :
        cosmics.tofits(out_file, c.cleanarray, header)
        
        fits.setval(out_file, keyword="PAPIVERS",
                        value= __version__, comment="PANIC Pipeline version")
        
        # If you want the mask, here it is :
        if want_mask:
            cosmics.tofits("mask.fits", c.mask, header)
            # (c.mask is a boolean numpy array, that gets converted here to an integer array)
    
    except Exception as e:
        log.error("Error removing cosmic rays in file : %s , Error %s:"%(in_image,str(e)))
        raise e
    
    return out_file
        
# main
if __name__ == "__main__":
    
    
    usage = "usage: %prog [options]"
    desc = "Remove the cosmic ray hits in the input image."
    parser = OptionParser(usage, description=desc)
    
    parser.add_option("-i", "--input_image",
                  action="store", dest="input_image", 
                  help="input image to remove cosmics")
                  
    parser.add_option("-o", "--output",
                  action="store", dest="output_image", 
                  help="output filename (default = %default)",
                  default="without_cosmics.fits")
    
    parser.add_option("-O", "--overwrite",
                  action="store_true", dest="overwrite", default=False,
                  help="overwrite the original image with the corrected one")

    parser.add_option("-m", "--mask",
                  action="store_true", dest="want_mask", default=False,
                  help="If true, the mask with cosmics detected and removed is written into a FITS file.")
                                
    (options, args) = parser.parse_args()
    
    if len(sys.argv[1:])<1:
       parser.print_help()
       sys.exit(0)
       
    if not options.input_image or len(args)!=0: 
        parser.print_help()
        parser.error("wrong number of arguments " )

    if not options.output_image:
        options.output_image = None

    try:    
        remove_cr(options.input_image, options.output_image, 
                         options.overwrite, options.want_mask)
    except Exception as e:
        log.error("Fail of remove_cr procedure: %s"%str(e))
    else:
        log.info("Well done!")
