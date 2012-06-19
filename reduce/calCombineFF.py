#! /usr/bin/env python
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
#
# PAPI (PANIC PIpeline)
#
# calCombineFF.py
#
# Created    : 29/03/2012    jmiguel@iaa.es -
# Last update: #! /usr/bin/env python

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
# TODO
################################################################################
"""
A trick to combine domeFF and skyFF :

Often for a run you have dome flats with an accumulated number of
electrons in the millions, but a poor match in illumination and color to the dark
sky. You also have a limited number of twilight flats or dark-sky images that can
be combined to make a dark-sky flat, but the total counts per pixel in either set of
flats is not very high. A fairly standard procedure is to 'median-smooth' dome
and twilight or dark-sky flat. A median smoothing replaces each pixel with the
median of the pixel values in a box of a given size on a side. The result is an image
that has been smoothed on the scale of the smoothing box size.
A procedure for taking advantage of the facts that the large-scale flat-field
variation of the dark-sky flat match that of the program frames and the dome flats
have very high S/N in each pixel goes as follows:
 
 (a) Median smooth the combined, dark-sky flat -this improves the S/N and
preserves the large-scale features of the flat.
 (b) Median smooth the combined dome flats using the same filter size as was
used for the dark-sky flat.
 (c) Divide the combined dome flat by it's median smoothed-version. The result is
a frame that is flat on large scales but contains all the high spatial frequency
flat-field information.
 (d) Now multiply the smoothed dark-sky frame and the result of the division in
the previous step. You now have a flat-field with the low spatial frequency
properties of the dark-sky flat combined with the high S/N, high spatial
frequency properties of the dome flat.

"""

# Import necessary modules
from optparse import OptionParser
import sys
import os

import pyfits
import numpy

import misc.fileUtils
import datahandler

# Logging
from misc.paLog import log

# IRAF
from pyraf import iraf
from iraf import noao
from iraf import mscred

def combineFF(domeFF, skyFF, combinedFF=None):
    """
    Combine a dome FF and a sky FF in order to have a FF with the low frequency
    proporties of the dark-sky flat combined with the high S/N, high spatial
    frequency properties of the dome FF.
    
    Basically:
    
     skyFF' = smooth(skyFF)
     domeFF' = smooth(domeFF)
     domeFF'' = domeFF / domeFF'
     combFF = skyFF' * modeFF''
     
    
    Parameters
    ----------
    domeFF : str
        Input filename of dome FF
    
    skyFF : str
        Input filename of the sky FF
        
    combinedFF: str
        Output filename of the combined FF generated
        
    Returns
    -------
    If all was successful, the name of the output file is returned
        
    """
    
    
    if not datahandler.isaFITS(domeFF) or not datahandler.isaFITS(skyFF):
        log.error("Some input FF is not a FITS file")
        raise Exception("Some input FF is not a FITS file")
        
    try:
        #smooth the domeFF
        log.debug("Doing Median smooth of domeFF ...")
        misc.fileUtils.removefiles(domeFF.replace(".fits", "_smooth.fits"))
        iraf.mscmedian(
                    input=domeFF,
                    output=domeFF.replace(".fits", "_smooth.fits"),
                    xwindow=20,
                    ywindow=20,
                    outtype="median")

        #Or using scipy ( a bit slower then iraf...)
        #from scipy import ndimage
        #filtered = ndimage.gaussian_filter(f[0].data, 20)
                    
        
        #smooth the skyFF
        log.debug("Doing Median smooth of skyFF ...")
        misc.fileUtils.removefiles(skyFF.replace(".fits", "_smooth.fits"))
        iraf.mscmedian(
                    input=skyFF,
                    output=skyFF.replace(".fits", "_smooth.fits"),
                    xwindow=20,
                    ywindow=20,
                    outtype="median"
                    )
                       
        # Divide domeFF by smoothed version
        misc.fileUtils.removefiles(domeFF.replace(".fits", "_div_smooth.fits"))
        iraf.imarith(operand1=domeFF,
                    operand2=domeFF.replace(".fits", "_smooth.fits"),
                    op='/',
                    pixtype='real',
                    result=domeFF.replace(".fits", "_div_smooth.fits"),
                    )

        # Combine skyFF with domeFF
        misc.fileUtils.removefiles(combinedFF)
        iraf.imarith(operand1=skyFF.replace(".fits", "_smooth.fits"),
                    operand2=domeFF.replace(".fits", "_div_smooth.fits"),
                    op='*',
                    pixtype='real',
                    result=combinedFF,
                    )
        
        # Remove all temporal
        misc.fileUtils.removefiles(domeFF.replace(".fits", "_smooth.fits"))
        misc.fileUtils.removefiles(skyFF.replace(".fits", "_smooth.fits"))
        misc.fileUtils.removefiles(domeFF.replace(".fits", "_div_smooth.fits"))
        
    except Exception,e:
        raise e
    
    log.info("Combined Flat-Field created : %s",combinedFF)
    
    return combinedFF
 
# main
if __name__ == "__main__":
    
    description = \
    "Module to combine a dome Flat-field and a sky Flat-field. "
    
    parser = OptionParser(description = description, 
                                   #formatter = wider_format, 
                                   usage = "usage: %prog [options] arg1 arg2 ...", 
                                   version = "%prog 1.0")
    
    parser.add_option("-d", "--domeFF",
                  action="store", dest="domeFF", 
                  help="input dome Flat-Field ")
    
    parser.add_option("-s", "--skyFF",
                  action="store", dest="skyFF", 
                  help="input sky Flat-Field ")
                  
    parser.add_option("-o", "--output",
                  action="store", dest="output_image", 
                  help="output filename of combined Flat-Field (default = %default)",
                  default="combinedFF.fits")
    
                                
    (options, args) = parser.parse_args()
    
    
    if not options.domeFF or not options.skyFF or len(args)!=0: 
    # args is the leftover positional arguments after all options have been processed
        parser.print_help()
        parser.error("wrong number of arguments " )
    if not options.output_image:
        options.output_image = None

    try:    
        combineFF(options.domeFF, options.skyFF, options.output_image)
    except Exception, e:
        log.error("Fail of combineFF procedure")
        raise e
    
    print "\nWell done !"
    