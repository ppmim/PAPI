#! /usr/bin/env python
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
# dxtalk.py
#
# Created    : 23/02/2012    jmiguel@iaa.es -
# Last update: 
# TODO
#   - object rejection in median cube computation
#   - smooth the cube ??    
################################################################################

# Import necessary modules
from optparse import OptionParser
import sys

import pyfits
import numpy

def de_crosstalk(in_image, out_image=None):
    """
    Remove cross-talk in O2k images
    """
 
    #log.debug("De-cross-talk procedure started !")
       
    if not out_image:
        out_file = in_image.replace(".fits","_dx.fits")
    else:
        out_file = out_image
        
    try:
        f_in = pyfits.open(in_image)
        if len(f_in)==1:
            data_in = f_in[0].data
        else:
            log.errro("MEF files currently not supported !")
    except Exception,e:
        log.error("Error openning FITS file : %s"%in_image)
        raise e
    
    #### Q1 ####
    n_stripes = 8
    width_st = 1024
    height_st = 128
    x_orig = 0
    y_orig = 0
    
    cube = numpy.zeros([n_stripes, height_st, width_st], dtype=numpy.float)
    data_out = numpy.zeros([n_stripes*height_st*2, width_st*2], dtype=numpy.float)
    
    for j in range (0,n_stripes):
        cube [j] = data_in[x_orig+j*height_st:x_orig+(j+1)*height_st, 
                           y_orig+0:y_orig+width_st]
        med_cube = numpy.median(cube,0)
        print "Shape1 = ",med_cube.shape
        
    for j in range(0,n_stripes):
        data_out[x_orig+j*height_st:x_orig+(j+1)*height_st, 
                 y_orig+0:y_orig+width_st] = cube[j]-med_cube
        
    #### Q3 ####
    x_orig = 1024
    y_orig = 1024
    width_st = 1024
    height_st = 128
    n_stripes = 8

    for j in range (0,n_stripes):
        cube [j] = data_in[x_orig+j*height_st:x_orig+(j+1)*height_st, 
                           y_orig+0:y_orig+width_st]
        med_cube = numpy.median(cube,0)
        
    for j in range(0,n_stripes):
        data_out[x_orig+j*height_st:x_orig+(j+1)*height_st, 
                 y_orig+0:y_orig+width_st] = cube[j]-med_cube
        
    
    #### Q2 ####
    n_stripes = 8
    width_st = 128
    height_st = 1024
    x_orig = 0
    y_orig = 1024
    
    cube = numpy.zeros([n_stripes, height_st, width_st], dtype=numpy.float)
    #cube.reshape((n_stripes, height_st, width_st))
    
    for j in range (0,n_stripes):
        cube [j] = data_in[x_orig:x_orig+height_st, 
                           y_orig+j*width_st:y_orig+(j+1)*width_st]
        med_cube = numpy.median(cube,0)
        
    for j in range(0,n_stripes):
        data_out[x_orig:x_orig+height_st, 
                 y_orig+j*width_st:y_orig+(j+1)*width_st] = cube[j]-med_cube

    #### Q4 ####
    n_stripes = 8
    width_st = 128
    height_st = 1024
    x_orig = 1024
    y_orig = 0
    
    for j in range (0,n_stripes):
        cube [j] = data_in[x_orig:x_orig+height_st, 
                           y_orig+j*width_st:y_orig+(j+1)*width_st]
        med_cube = numpy.median(cube,0)
        
    for j in range(0,n_stripes):
        data_out[x_orig:x_orig+height_st, 
                 y_orig+j*width_st:y_orig+(j+1)*width_st] = cube[j]-med_cube

    ### write FITS ###
        
    hdu = pyfits.PrimaryHDU()
    hdu.scale('float32') # important to set first data type
    hdu.data = data_out
    hdulist = pyfits.HDUList([hdu])
    
    hdr0 = pyfits.getheader(in_image)
    hdr0.update('HISTORY', 'De-crosstalk procedure executed ')
    
    try:
        hdulist.writeto(out_image, header=hdr0, clobber=clobber)
        hdulist.close(output_verify='ignore')
    except Exception,e:
        raise e
    
# main
if __name__ == "__main__":
    
    print "Probando DX-talk"
    
    usage = "usage: %prog [options] arg1 arg2 ..."
    parser = OptionParser(usage)
    
    parser.add_option("-i", "--input_image",
                  action="store", dest="input_image", 
                  help="input image to remove crosstalk ")
                  
    parser.add_option("-o", "--output",
                  action="store", dest="output_image", 
                  help="output filename (default = %default)",
                  default="dxtalk.fits")
                                
    (options, args) = parser.parse_args()
    
    
    if not options.input_image or len(args)!=0: 
    # args is the leftover positional arguments after all options have been processed
        parser.print_help()
        parser.error("wrong number of arguments " )
    if not options.output_image:
        options.output_image = None
        
    de_crosstalk (options.input_image, options.output_image)
    
    print "Se Acabo !"
    