#! /usr/bin/env python
# Copyright (c) 2013 IAA-CSIC  - All rights reserved. 
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
# refpixcorr.py
#
# Created    : 21/Jan/2013    jmiguel@iaa.es -
# Last update: 
# TODO
################################################################################
"""
Reference pixels subtraction.

HAWAII-2RG detectors have an effective surface of 2040x2040 sensitive pixels. 
A 4-pixel wide border is used as reference to correct for relatively slow bias 
drifts, first seen during PANIC commissionning as column to column DC offsets 
on a scale of several tens of columns.

Here, it is assumed that each column is approximately constant and that any 
drift is seen as a level change between columns. In other words, we correct only
the presence of vertical bands. Only the rows of 4 reference pixels on top and 
at the bottom of images are used, those on the left and on the right are not used.

First, we correct the reference pixels themselves because a small DC level slope
is seen from line 0 to line 3 and from line 2047 to line 2044 (most probably due
to capacitive coupling effects). So, each of those 8 lines has its median 
subtracted from it.

Then, we consider the 4 reference pixels on top and at the bottom of each column. 
We consider they follow any DC change that affects the whole column. The idea 
is to subtract the DC level seen on those 8 reference pixels from the 2040 data 
pixels of the corresponding column.

We construct a 2048x8 matrix of reference pixels that we median along the 
vertical axis to obtain a 2048x1 reference line (that median of 8 reference 
pixels is to reduce noise in the measurement). Then, to further reduce noise 
and spurious effects, we smooth the 2048x1 line with a boxchar of 9. Then we 
subtract the corresponding scalar value for each data column (e.g. 
column 5 = column 5 minus 5th index of reference line).

Note: Taken from WIRCam Web http://bit.ly/YP5heF

"""

# Import necessary modules
from optparse import OptionParser
import sys

import pyfits
import numpy
import scipy.ndimage.filters

# Logging
from misc.paLog import log

def refPixelCorrection(in_image, out_image=None, overwrite=False):
    """
    Run reference pixel correction of PANIC images.
    
    Parameters
    ----------
    in_image : str
        Input filename to be corrected; it is assumed to be a 'GEIRS' format,
        i.e., 4kx4k in a single HDU.
    
    out_image : str
        Output filename of corrected image
        
    overwrite: Boolean
        If true, the input file 'in_image' filename will be overwritten,
        otherwise, the 'out_image' filename will be used as output. 
    
    Returns
    -------
    If all was successful, the name of the output file is returned
        
    """
    
    try:
        if (pyfits.getval(in_image, 'INSTRUME').lower()!='panic' or 
            pyfits.getval(in_image,'NAXIS1')!=4096 or 
            pyfits.getval(in_image, 'NAXIS2')!=4096):
            log.errro("Error, expected 4kx4k PANIC image.")
    except Exception, e:
        log.errro("Error, cannot read INSTRUME keyword.")
        raise e


    if overwrite:
        out_file = in_image
    else:   
        if not out_image:
            out_file = in_image.replace(".fits","_refcorr.fits")
        else:
            out_file = out_image
    
    # open the file
    try:
        f_in = pyfits.open(in_image)
        if len(f_in)==1:
            imraw = f_in[0].data
        else:
            log.errro("MEF files currently not supported !")
            raise Exception("MEF files currently not supported !")
    except Exception, e:
        log.error("Error opening FITS file : %s"%in_image)
        raise e
    
    
    #Set the lines to be used as reference pixels, those at the top and bottom of
    #the image.
    lineQ12 = [0,1,2,3,    # top Q1 and Q2
            2044,2045,2046,2047] # bottom Q1 and Q2
    lineQ34 = [2048,2049,2050,2051, # top Q3 y Q4
            4092,4093,4094,4095] # bottom Q3 y Q4
    
    
    nquads = 4 # number of quadrants, i.e.,detectors
    med = numpy.zeros([len(lineQ12), nquads], dtype=numpy.float)
    im = imraw
    
    #Try to correct the small slope seen in the reference pixels themselves, i.e.
    #from line 0 to 3 and from 2047 to 2044 by subtracting the median horizontal 
    #value from itself.
    reflines = numpy.zeros([len(lineQ12), 2048, nquads], dtype=numpy.float)

    #import pdb; pdb.set_trace()

    for i in range(len(lineQ12)):
        #Determine offset = median of each line
        #med[i] = numpy.median(im[line[i], :])
        med[i, 0] = numpy.median(im[lineQ12[i], 0:2048])
        med[i, 1] = numpy.median(im[lineQ12[i], 2048:4096])
        med[i, 2] = numpy.median(im[lineQ34[i], 0:2048])
        med[i, 3] = numpy.median(im[lineQ34[i], 2048:4096])
        print "MED_line Q1 -->",med[:,0]
        print "MED_line Q2 -->",med[:,1]
        print "MED_line Q3 -->",med[:,2]
        print "MED_line Q4 -->",med[:,3]
        #Apply offset so all columns have a median of zero
        reflines[i,:,0] = im[lineQ12[i], 0:2048] - med[i, 0]
        reflines[i,:,1] = im[lineQ12[i], 2048:4096] - med[i, 1]
        reflines[i,:,2] = im[lineQ34[i], 0:2048] - med[i, 2]
        reflines[i,:,3] = im[lineQ34[i], 2048:4096] - med[i, 3]


    #Now, for each column, use the median of the 8 reference pixels available.
    #Crunch 8 lines into 1. This is a vector of 2048 reference columns.
    averagecolumn = numpy.median(reflines, axis=0)
    #Smooth that column with box of 9. This is to prevent introducing noise but 
    #it also means that column to column variations on shorter timescales can not
    #be corrected for.
    box = 9
    #averagecolumn = smooth(averagecolumn, box)
    averagecolumn[:,0] = scipy.ndimage.filters.median_filter(averagecolumn[:,0], size=box)
    averagecolumn[:,1] = scipy.ndimage.filters.median_filter(averagecolumn[:,1], size=box)
    averagecolumn[:,2] = scipy.ndimage.filters.median_filter(averagecolumn[:,2], size=box)
    averagecolumn[:,3] = scipy.ndimage.filters.median_filter(averagecolumn[:,3], size=box)
    #Subtract the average, smoothed column from each data column

    #import pdb; pdb.set_trace()

    #Q1 and Q2
    for i in range(4, 2044):
        im[i,0:2048] = imraw[i,0:2048] - averagecolumn[:,0]
        im[i,2048:4096] = imraw[i,2048:4096] - averagecolumn[:,1]

    #Q3 and Q4
    for i in range(2052, 4096):
        im[i,0:2048] = imraw[i,0:2048] - averagecolumn[:,2]
        im[i,2048:4096] = imraw[i,2048:4096] - averagecolumn[:,3]
 
    ### write FITS ###
            
    hdu = pyfits.PrimaryHDU()
    hdu.scale('float32') # important to set first data type
    hdu.data = im
    hdulist = pyfits.HDUList([hdu])
    
    hdr0 = pyfits.getheader(in_image)
    hdr0.add_history('Reference pixel correction done.')
    hdu.header = hdr0
    
    try:
        hdulist.writeto(out_file, output_verify='ignore', clobber=overwrite)
        hdulist.close(output_verify='ignore')
    except Exception, e:
        log.error("Cannot write image %s"%out_file)
        raise e
  
    log.debug("End of refPixelCorrection")

    return out_file
    
# main
if __name__ == "__main__":
    
    
    usage = "usage: %prog [options]"
    desc = "Reference pixel correction of the input image."
    parser = OptionParser(usage, description=desc)
    
    parser.add_option("-i", "--input_image",
                  action="store", dest="input_image", 
                  help="input image to be corrected.")
                  
    parser.add_option("-o", "--output",
                  action="store", dest="output_image", 
                  help="output filename (default = %default)",
                  default="corrected.fits")
    
    parser.add_option("-O", "--overwrite",
                  action="store_true", dest="overwrite", default=False,
                  help="overwrite the original image with the corrected one")

    # TODO
    #parser.add_option("-S", "--check_stars",
    #              action="store_true", dest="check_stars", default=False,
    #              help="check if there are bright stars and take them into account for the cube median")  
                                
    (options, args) = parser.parse_args()
    
    if len(sys.argv[1:])<1:
       parser.print_help()
       sys.exit(0)

    if not options.input_image or len(args)!=0: 
    # args is the leftover positional arguments after all options have been processed
        parser.print_help()
        parser.error("wrong number of arguments " )
    
    if not options.output_image:
        options.output_image = None

    try:    
        refPixelCorrection(options.input_image, options.output_image, 
                         options.overwrite)
    except Exception, e:
        log.error("Fail of refPixelCorrection procedure: %s"%str(e))
    else:
        log.info("Well done!")
    