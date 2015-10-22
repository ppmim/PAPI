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
#   - MEF support
#   - save memory using data_in.reshape for cube median computation
#   - In order not normalize wrt quadrant Q1, we should divide by the 
#     ratio median_Qn/median_Q1
################################################################################
"""
From "Characterization, Testing and Operation of Omega2000 Wide Field Infrared
Camera", Zoltan Kovacs et.al.

Although bright stars can saturate the detector, resetting of the full array
prevents this excess in the pixel values from causing any residual image 
effects in the following image of the dithering. Nevertheless, the satured
pixels generate a crosstalk between the data transfer lines of the different
channels of the quadrant in which they are situated. The data lines of the 
channels are organized in parallel and there might be an interference between 
the data lines transferring the high video signal and the neighbour ones. As a 
result of this crosstalk, a series of spots with the distances of 128 pixels 
from each other appeares in the whole quadrant, corresponding to each channel. 
The average values of the spots were lower than the background signal and their
difference was few percent, which is large enough to degrade the photometric
correctness at the place they are situated. These spots could not be measured
in the raw images but they were well discernible in the reduced frames (Fig. 9). 
This effect was a general feature of the operation of all the  HAWAII-2 detectors 
we tested and should be considered for the choice of pointing positions in any 
field of next observations.
"""

# Import necessary modules
from optparse import OptionParser
import sys

import astropy.io.fits as fits
import numpy

# Logging
from misc.paLog import log
from misc.version import __version__

def remove_crosstalk(in_image, out_image=None, overwrite=False):
    """
    Remove cross-talk in O2k or PANIC images.
    
    Parameters
    ----------
    in_image : str
        Input filename to be decrosstalk
    
    out_image : str
        Output filename of decrosstalked image
        
    overwrite: Boolean
        If true, the input file 'in_image' filename will be overwritten,
        otherwise, the 'out_image' filename will be used as output. 
    
    Returns
    -------
    If all was successful, the name of the output file is returned
        
    """
    
    try:
        if fits.getval(in_image, 'INSTRUME').lower()=='omega2000':
            return de_crosstalk_o2k(in_image, out_image, overwrite)
        elif fits.getval(in_image, 'INSTRUME').lower()=='panic':
            # we suppose an image of one single detector.
            return de_crosstalk_PANIC(in_image, out_image, overwrite)
        else:
            log.error("Instrument is not supported !")
            raise Exception("Instrument is not supported !")
    except Exception, e:
        raise e

    
 
def de_crosstalk_o2k(in_image, out_image=None, overwrite=False):
    """
    Remove cross-talk in O2k images (2kx2k).
    
    The image structure expected is as follow:
    
        I-----------------I
        I        I        I
        I   Q4   I   Q3   I
        I        I        I
        I-----------------I
        I        I        I
        I   Q1   I   Q2   I
        I        I        I
        I-----------------I
        
    where each quadrant (Qn) is 1kx1k and has 8 horizontal (Q1,Q3) or vertical
    (Q2,Q4) stripes of 128 pixels of length (width or heigh). 
    So, so quadrant pairs (Q1,Q3) and (Q2,Q4) are processed in the same way.
    """

    log.debug("Start remove_crosstalk (O2k)")

    if overwrite:
        out_file = in_image
    else:   
        if not out_image:
            out_file = in_image.replace(".fits","_dx.fits")
        else:
            out_file = out_image
            
    try:
        f_in = fits.open(in_image)
        if len(f_in)==1:
            data_in = f_in[0].data
        else:
            log.errro("MEF files currently not supported !")
            raise Exception("MEF files currently not supported !")
            
        if f_in[0].header['INSTRUME'].lower()!='omega2000':
            log.error("Only O2k instrument is supported !")
            raise Exception("Only O2k instrument is supported !")
    except Exception, e:
        log.error("Error opening FITS file : %s"%in_image)
        raise e
    
    background = numpy.median(data_in)
    #print "Image background estimation = ", background
    
    #### Q1 #### left-bottom, horizontal stripes 
    n_stripes = 8 # = no. channels
    width_st = 1024
    height_st = 128
    x_orig = 0
    y_orig = 0
    
    cube = numpy.zeros([n_stripes, height_st, width_st], dtype=numpy.float)
    data_out = numpy.zeros([n_stripes*height_st*2, width_st*2], dtype=numpy.float32)
    median = numpy.zeros([4], dtype=numpy.float32)
    
    for j in range (0, n_stripes):
        cube [j] = data_in[x_orig+j*height_st:x_orig+(j+1)*height_st, 
                           y_orig+0:y_orig+width_st]
    
    med_cube = numpy.median(cube, 0)
    median[0] = numpy.median(med_cube)
    #print "CUBE_MEDIAN[0] = ", median[0]
        
    for j in range(0, n_stripes):
        # subtract cube_median and add constant (skybkg) to preserve original count level
        data_out[x_orig+j*height_st:x_orig+(j+1)*height_st, 
                 y_orig+0:y_orig+width_st] = (cube[j] - med_cube) + background # median[0]
        
    #### Q3 #### right-top, horizontal stripes 
    x_orig = 1024
    y_orig = 1024
    width_st = 1024
    height_st = 128
    n_stripes = 8

    for j in range (0, n_stripes):
        cube [j] = data_in[x_orig+j*height_st:x_orig+(j+1)*height_st, 
                           y_orig+0:y_orig+width_st]

    med_cube = numpy.median(cube, 0)
    median[2] = numpy.median(med_cube)
    #print "CUBE_MEDIAN[2] = ", median[2]
        
    for j in range(0,n_stripes):
        # subtract cube_median and add constant (skybkg) to preserve original count level
        data_out[x_orig+j*height_st:x_orig+(j+1)*height_st, 
                 y_orig+0:y_orig+width_st] = (cube[j] - med_cube) + background #median[2]
        
    
    #### Q2 #### right-bottom, vertical stripes 
    n_stripes = 8
    width_st = 128
    height_st = 1024
    x_orig = 0
    y_orig = 1024
    
    cube = cube.reshape((n_stripes, height_st, width_st))
    
    for j in range (0,n_stripes):
        cube [j] = data_in[x_orig:x_orig+height_st, 
                           y_orig+j*width_st:y_orig+(j+1)*width_st]
    
    med_cube = numpy.median(cube, 0)
    median[1] = numpy.median(med_cube)
    #print "CUBE_MEDIAN[1] = ", median[1]
    
    for j in range(0, n_stripes):
        # subtract cube_median and add constant (skybkg) to preserve original count level
        data_out[x_orig:x_orig+height_st, 
                 y_orig+j*width_st:y_orig+(j+1)*width_st] = (cube[j] - med_cube) + background #median[1]

    #### Q4 #### left-top, vertical stripes  
    n_stripes = 8
    width_st = 128
    height_st = 1024
    x_orig = 1024
    y_orig = 0
    
    for j in range (0,n_stripes):
        cube [j] = data_in[x_orig:x_orig+height_st, 
                           y_orig+j*width_st:y_orig+(j+1)*width_st]

    med_cube = numpy.median(cube, 0)
    median[3] = numpy.median(med_cube)
    #print "CUBE_MEDIAN[3] = ", median[3]
    
    for j in range(0, n_stripes):
        # subtract cube_median and add constant (skybkg) to preserve original count level
        data_out[x_orig:x_orig+height_st, 
                 y_orig+j*width_st:y_orig+(j+1)*width_st] = (cube[j] - med_cube) + background #median[3]


    ##TODO: In order not normalize wrt quadrant Q1, we should divide by the 
    # ratio median_Qn/median_Q1
    """
    NOTA : Aunque mejora, NO termina de funcionar bien, se sigue notando una 
    leve diferencia en el nivel de fodo de los cuandrantes !!
    # Q1
    data_out[0:1024,0:1024] *=median[0]/median[0]
    # Q2
    data_out[0:1024,1024:2048] *=median[0]/median[1]
    # Q3
    data_out[1024:2048,1024:2048] *=median[0]/median[2]
    # Q4
    data_out[1024:2048,0:1024] *=median[0]/median[3]
    """  
    
    """
    Other test to do a zero offset normalization, instead of mult. scale normalization, much more 'dangerous'
    --> But this method does not work either....quadrant offset level is still on images !!!
    """
    """
    # Q1
    print "SCALE[0]=",(numpy.mean(median)-median[0])
    data_out[0:1024,0:1024] += (numpy.mean(median)-median[0])
    # Q2
    print "SCALE[1]=",(numpy.mean(median)-median[1])
    data_out[0:1024,1024:2048] += (numpy.mean(median)-median[1])
    # Q3
    print "SCALE[2]=",(numpy.mean(median)-median[2])
    data_out[1024:2048,1024:2048] += (numpy.mean(median)-median[2])
    # Q4
    print "SCALE[3]=",(numpy.mean(median)-median[3])
    data_out[1024:2048,0:1024] += (numpy.mean(median)-median[3])
    """
    
    #mode = 3*numpy.median(data_in[200:1800,200:1800])-2*numpy.mean(data_in[200:1800,200:1800])
    #data_out = data_out + mode
    #print "Global MODE = ",mode
    
    ### write FITS ###
    
            
    hdu = fits.PrimaryHDU()
    hdu.scale('float32') # important to set first data type
    hdu.data = data_out
    hdulist = fits.HDUList([hdu])
    
    hdr0 = fits.getheader(in_image)
    hdr0.add_history('De-crosstalk procedure executed ')
    hdr0.set('PAPIVERS', __version__, 'PANIC Pipeline version')
    hdu.header = hdr0
    
    try:
        hdulist.writeto(out_file, output_verify='ignore', clobber=overwrite)
        hdulist.close(output_verify='ignore')
    except Exception, e:
        raise e
  
    log.debug("End of remove_crosstalk (O2k)")

    return out_file

def de_crosstalk_PANIC_full_detector(in_image, out_image=None, overwrite=False):
    """
    ==> NOT USED !!
    
    Remove cross-talk in PANIC images (4kx4k).
    The frame structure expected is as follow:
    
        I-----------------I 
        I        I        I 
        I   Q4   I   Q3   I 
        I        I        I 
        I-----------------I 
        I        I        I 
        I   Q1   I   Q2   I 
        I        I        I 
        I-----------------I 
        
    where each quadrant (Qn) is 2kx2k and has 32 horizontal stripes of 64 pixels
    of height. So, all the quadrant are processed in the same way.  
    
    Parameters:
    ----------
    in_image: str
        Filename of the FITS image to remove the xtalk. It must be a single FITS,
        MEF are not supported yet.

    out_image: str
        Filaname of the output dxtalked image.

    overwrite: bool
        Flag to indicate if the out_image can be overwritten. 
    """
    
    log.debug("Start remove_crosstalk (PANIC)")

    if overwrite:
        out_file = in_image
    else:   
        if not out_image:
            out_file = in_image.replace(".fits", "_dx.fits")
        else:
            out_file = out_image
            
    try:
        f_in = fits.open(in_image)
        if len(f_in)==1:
            data_in = f_in[0].data
        else:
            log.errro("MEF files currently not supported !")
            raise Exception("MEF files currently not supported !")
            
        if f_in[0].header['INSTRUME'].lower()!='panic':
            log.error("Instrument %s is not supported !"%f_in[0].header['INSTRUME'])
            raise Exception("Instrument is not supported !")
    except Exception, e:
        log.error("Error openning FITS file : %s"%in_image)
        raise e
    
    background = numpy.median(data_in)
    #print "Image background estimation = ", background

    #### Q1 #### left-bottom, horizontal stripes 
    n_stripes = 32 # = no. channels
    width_st = 2048
    height_st = 64
    x_orig = 0
    y_orig = 0
    
    
    
    cube = numpy.zeros([n_stripes, height_st, width_st], dtype=numpy.float)
    data_out = numpy.zeros([n_stripes*height_st*2, width_st*2], dtype=numpy.float32)
    
    for j in range (0, n_stripes):
        cube [j] = data_in[x_orig+j*height_st:x_orig+(j+1)*height_st, 
                           y_orig+0:y_orig+width_st]
    
    med_cube = numpy.median(cube, 0)
    median = numpy.median(med_cube)
    #print "CUBE_MEDIAN = ", median
        
    for j in range(0, n_stripes):
        # subtract cube_median and add constant (skybkg) to preserve original count level
        data_out[x_orig+j*height_st:x_orig+(j+1)*height_st, 
                 y_orig+0:y_orig+width_st] = (cube[j]-med_cube) + background #median
        
    #### Q3 #### right-top, horizontal stripes 
    n_stripes = 32 # = no. channels
    width_st = 2048
    height_st = 64
    x_orig = 2048
    y_orig = 2048

    for j in range (0, n_stripes):
        cube [j] = data_in[x_orig+j*height_st:x_orig+(j+1)*height_st, 
                           y_orig+0:y_orig+width_st]

    med_cube = numpy.median(cube, 0)
    median = numpy.median(med_cube)
    #print "CUBE_MEDIAN = ", median
        
    for j in range(0, n_stripes):
        # subtract cube_median and add constant (skybkg) to preserve original count level
        data_out[x_orig+j*height_st:x_orig+(j+1)*height_st, 
                 y_orig+0:y_orig+width_st] = (cube[j]-med_cube) + background #median
        
    
    #### Q2 #### right-bottom, horizontal stripes 
    n_stripes = 32 # = no. channels
    width_st = 2048
    height_st = 64
    x_orig = 0
    y_orig = 2048
    
    cube = cube.reshape((n_stripes, height_st, width_st))
    
    for j in range (0, n_stripes):
        cube [j] = data_in[x_orig+j*height_st:x_orig+(j+1)*height_st, 
                           y_orig+0:y_orig+width_st]

    med_cube = numpy.median(cube, 0)
    median = numpy.median(med_cube)
    #print "CUBE_MEDIAN = ", median

        
    for j in range(0, n_stripes):
        # subtract cube_median and add constant (skybkg) to preserve original count level
        data_out[x_orig+j*height_st:x_orig+(j+1)*height_st, 
                 y_orig+0:y_orig+width_st] = (cube[j]-med_cube) + background #median

    #### Q4 #### left-top, horizontal stripes 
    n_stripes = 32 # = no. channels
    width_st = 2048
    height_st = 64
    x_orig = 2048
    y_orig = 0
    
    for j in range (0,n_stripes):
        cube [j] = data_in[x_orig+j*height_st:x_orig+(j+1)*height_st, 
                           y_orig+0:y_orig+width_st]

    med_cube = numpy.median(cube, 0)
    median = numpy.median(med_cube)
    #print "CUBE_MEDIAN = ", median

        
    for j in range(0, n_stripes):
        # subtract cube_median and add constant (skybkg) to preserve original count level
        data_out[x_orig+j*height_st:x_orig+(j+1)*height_st, 
                 y_orig+0:y_orig+width_st] = (cube[j]-med_cube) + background #median

    ### write FITS ###
    
            
    hdu = fits.PrimaryHDU()
    hdu.scale('float32') # important to set first data type
    hdu.data = data_out
    hdulist = fits.HDUList([hdu])
    
    hdr0 = fits.getheader(in_image)
    hdr0.add_history('De-crosstalk procedure executed ')
    hdr0.set('PAPIVERS', __version__, 'PANIC Pipeline version')
    hdu.header = hdr0
    
    try:
        hdulist.writeto(out_file, output_verify='ignore', clobber=overwrite)
        hdulist.close(output_verify='ignore')
    except Exception, e:
        raise e
      
    log.debug("End of remove_crosstalk (PANIC)")
    
    return out_file

def de_crosstalk_PANIC(in_image, out_image=None, overwrite=False):
    """
    Remove cross-talk in PANIC single detectors (2kx2k).
    The frame structure of a full frame of PANIC is as follow:
    
        I-----------------I 
        I        I        I 
        I   Q3   I   Q4   I 
        I        I        I 
        I-----------------I 
        I        I        I 
        I   Q1   I   Q2   I 
        I        I        I 
        I-----------------I 
        
    where each quadrant (Qn) is 2kx2k and has 32 horizontal stripes of 64 pixels
    of width. So, all the quadrant are processed in the same way.
    
    The procedure does not need to read the DET_ID keyword to know the detector or how
    to stripes are distributed along the detector. All detectors have the
    same orientation of the stipes. So, no matter which detector we are proceesing !
    
    
    DET_ID = SG1 -- Q1 - vertical stripes
    DET_ID = SG2 -- Q2 - vertical stripes
    DET_ID = SG3 -- Q3 - vertical stripes
    DET_ID = SG4 -- Q4 - vertical stripes
    
    """
    
    log.debug("Start remove_crosstalk (PANIC)")

    if overwrite:
        out_file = in_image
    else:   
        if not out_image:
            out_file = in_image.replace(".fits", "_dx.fits")
        else:
            out_file = out_image
            
    try:
        f_in = fits.open(in_image)
        if len(f_in)==1:
            data_in = f_in[0].data
        else:
            log.errro("MEF files currently not supported !")
            raise Exception("MEF files currently not supported !")
            
        if f_in[0].header['INSTRUME'].lower()!='panic':
            log.error("Instrument %s is not supported !"%f_in[0].header['INSTRUME'])
            raise Exception("Instrument is not supported !")
    except Exception, e:
        log.error("Error openning FITS file : %s"%in_image)
        raise e
    
    background = numpy.median(data_in)
    print "Image background estimation = ", background

    # All detectors have 32 vertical_stripes of 2048x64 (rows x columns) each one
    # NOTE: in python, x=rows and y=columns
    n_stripes = 32
    width_st = 64
    height_st = 2048
    x_orig = 0
    y_orig = 0
    
    cube = numpy.zeros([n_stripes, height_st, width_st], dtype=numpy.float32)
    data_out = numpy.zeros([height_st, n_stripes*width_st], dtype=numpy.float32)
    
    for j in range (0,n_stripes):
        cube [j] = data_in[x_orig:x_orig+height_st, 
                           y_orig+j*width_st:y_orig+(j+1)*width_st]

    med_cube = numpy.median(cube, 0)
    median = numpy.median(med_cube)
    
    for j in range(0, n_stripes):
        # subtract cube_median and add constant (skybkg) to preserve original count level
        data_out[x_orig:x_orig+height_st, 
                 y_orig+j*width_st:y_orig+(j+1)*width_st] = (cube[j] - med_cube) + background #median
        

    ### write FITS ###
    
    hdu = fits.PrimaryHDU()
    hdu.scale('float32') # important to set first data type
    hdu.data = data_out
    hdulist = fits.HDUList([hdu])
    
    hdr0 = fits.getheader(in_image)
    hdr0.add_history('De-crosstalk procedure executed ')
    hdr0.set('PAPIVERS', __version__, 'PANIC Pipeline version')
    hdu.header = hdr0
    
    try:
        hdulist.writeto(out_file, output_verify='ignore', clobber=overwrite)
        hdulist.close(output_verify='ignore')
    except Exception, e:
        raise e
      
    log.debug("End of remove_crosstalk (PANIC)")
    
    return out_file

# main
if __name__ == "__main__":
    
    
    usage = "usage: %prog [options]"
    desc = "Remove the cross-talk spots in the input image."
    parser = OptionParser(usage, description=desc)
    
    parser.add_option("-i", "--input_image",
                  action="store", dest="input_image", 
                  help="input image to remove crosstalk ")
                  
    parser.add_option("-o", "--output",
                  action="store", dest="output_image", 
                  help="output filename (default = %default)",
                  default="dxtalk.fits")
    
    parser.add_option("-O", "--overwrite",
                  action="store_true", dest="overwrite", default=False,
                  help="overwrite the original image with the corrected one")

    # TODO
    #parser.add_option("-S", "--check_stars",
    #              action="store_true", dest="check_stars", default=False,
    #              help="check if there are bright stars and take them into account for the cube median")  
                                
    (options, args) = parser.parse_args()
    
    if len(sys.argv[1:]) < 1:
       parser.print_help()
       sys.exit(0)
       
    if not options.input_image or len(args)!=0: 
    # args is the leftover positional arguments after all options have been processed
        parser.print_help()
        parser.error("wrong number of arguments " )
    if not options.output_image:
        options.output_image = None

    try:    
        remove_crosstalk(options.input_image, options.output_image, 
                         options.overwrite)
    except Exception, e:
        log.error("Fail of Crosstalk procedure: %s"%str(e))
    else:
        log.info("Well done!")
    
