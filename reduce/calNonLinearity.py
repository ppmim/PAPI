#!/usr/bin/env python

# Copyright (c) 2009-2012 IAA-CSIC  - All rights reserved. 
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
# nonLinearity.py
#    
# Created    : 16/12/2009    jmiguel@iaa.es -
# Last update: 15/02/2012    jmiguel@iaa.es
#
# TODO : computation of non-linearity based on a serie or darks ???
#       
################################################################################
# Import necessary modules
"""
Report on non-linearity measurements of the LUCI1 detector, 27.04.2011, by JDK, MPE.
------------------------------------------------------------------------------------

To obtain the raw data from the LBT archive, just search for FILENAME:
luci.20110407. This will result in a list of 228 fits files, the first
34 are dark frames, the others are flat fields with increasing
exposure times which can be used to measure the non-linearity.

The raw data was taken as follows, for the o2dcr readout mode: 2 flat
exposures with T_exp=t0, 2 flats with T_exp=2s, 2 flats with T_exp=t1,
2 flats with T_exp=2s, ..., T_exp=tn, T_exp=2s, T_exp=tn-1,
... T_exp=t0. A similar series is available for the mer readout
mode. Each exposure time is therefore sampled four times. All four
samples are used in my analysis. The mer mode indeed shows a small,
but significant, difference between the increasing and decreasing part
of the series (which is not taken into account in the analysis).

For each pixel, I have fitted a linear function to the counts as a
function of exposure time, for the first four exposure times only,
where the counts are below 20000, and the detections should still be
linear (but this can be debated). Then, I have fitted a third-order
polynomial to the data divided by the linear fit, as a function of
detected counts (in units of 10^4), up to deviations of 4%. The first
coefficient was fixed to 1, leaving only two free parameters. Attempts
to fit up to larger deviations resulted in erratic behaviour of the
fitted functions. These fits (and data) are shown for 36 pixels and
both modes in two pdf files available from the website specified below
(linearity_individual_mer.pdf and inearity_individual_o2dcr.pdf). Blue
data points are from the first half of the data set (i.e. increasing
exposure times), red data points are from the second half
(i.e. decreasing exposure times). Data points with a black circle are
not fitted. Pixels showing abnormal behaviour (e.g., pixel 400,400,
second plot in second row) are marked, and the fit parameters are set
to zero (resulting in a correction factor of one). As a byproduct,
therefore, a bad pixel map is created for each readout mode. Note that
pixel at the edges of the detector are also marked as bad, but the
cause is most probably absent illumination rather than a detector
problem.

The two free parameters are saved in fits files (one for each readout
mode): nonlin_coeffs_mer.fits and nonlin_coeffs_o2dcr.fits, which can
be used to correct LUCI data for non-linearity (up to 4%). In
addition, the bad pixel maps are called: badpixel_mer.fits and
badpixel_o2dcr.fits. These files are available from the MPE LUCIFER
website, at: www.mpe.mpg.de/ir/lucifer/links.php

A python script is also available from the website, which uses the
fits files to correct a LUCIFER fits file:
luci_correct_linearity.py. It is, however, also very simple to use in
iraf or IDL, as shown below.

In IRAF:
imarith myfile / 1e4 myfile4
imfunc  myfile4 square square
imarith nonlin_coeffs_o2dcr[*,*,1] * square fac2
imarith nonlin_coeffs_o2dcr[*,*,2] * myfile4 fac1
imarith fac1 + fac2 fac
imarith fac + 1 fac
imarith myfile / fac myfile_lincor
imdelete myfile4, fac, fac1, fac2

In IDL:
arr = readfits('myfile.fits') / 1e4
poly = readfits('nonlin_coeffs_o2dcr.fits')
fac = 1 + poly[*,*,1] * arr + poly[*,*,0] * arr^2
writefits,'myfile_lincor.fits',arr/fac*1e4

Note the factor 10^4 and the backward order of the polynomial coefficients.

"""

import getopt
import sys
import os
import logging
import fileinput
import time
from optparse import OptionParser

import misc.fileUtils
import misc.utils as utils
import datahandler

# Interact with FITS files
import pyfits
import numpy

# Logging
from misc.paLog import log

class NonLinearityModel(object):
    """
    Class used to compute the Non-linearity model for the detectors
    Based on the algorithm described by T.H.Jarrett et all (1994) we compute the non-linearity coefficientes
    of the model
    """
    def __init__(self, input_files=None, output_dir="/tmp", output_filename="/tmp/NLC.fits"):
        """
        Init the object
        
        Parameters
        ----------
        
        input_data: list
            A list of sky files
        input_files: list
            A list of FITS files as input for the computation
        output_filename: str
            File where coefficientes will be saved
        
        Returns
        -------
        
        coeff: array
            Array of values a0,a1,a2 where 
            a0 = bias level
            a1 = quantum efficiency or sensitivity
            a2 = the non-linearity or saturation
        """
        self.__input_files = input_files
        self.__output_file_dir = output_dir
        self.__output_filename = output_filename  # full filename (path+filename)
    
    def createModel(self):
          
        """
        Compute the Non-linearity correction model based on a serie of FITS files
        taken in order to get check the detector linearity.
         
        Parameters
        ----------
        
        Return
        ------
        
        
        """   
        log.debug("Start createModel")
        log.warn("Non-Linearity model computation is NOT yet implemented !!")
        
        return
    
    
        start_time = time.time()
        t=utils.clock()
        t.tic()
        
        # Get the user-defined list of dark frames
        framelist = self.__input_files
        
        # STEP 0: Determine the number of frames 
        try:    
            nframes = len(framelist)
        except IndexError:
            log.error("No DARK frames defined")
            raise
        
        if not os.path.exists(os.path.dirname(self.__output_filename)):
            raise NameError, 'Wrong output path'
        if not self.__output_filename:
            log.error("Combined DARK frame not defined")
            raise "Wrong output filename"
    
        log.debug('Saved NLC Model to %s' , self.__output_filename)
        log.debug("createModel' finished %s", t.tac() )
        
        return True
    
    def applyModel(self, source, model):
          
        """
        Do the Non-linearity correction using the supplied model
        
        Parameters
        ----------
        
        source : str
            FITS filename to be corrected
        
        model : str 
            FITS filename of the Non-Linearity model, ie., containing polynomial 
            coeffs for correction that must has been previously computed
        
        
        Returns
        -------
        
        
        TODO
        ----
        - adjust the read modes availables on PANIC (lir, mer, o2dcr, ...)
        - read the proper BadPixelMask ??
        
        """   
        log.debug("Start applyModel")
        
        if not os.path.exists(source):
            log.error("Cannot read source file %s"%source)
            raise Exception("Cannot read source file")
        
        if not os.path.exists(model):
            log.error("Cannot read non-linearity model file %s"%model)
            raise Exception("Cannot read non-linearity model file")
        
        # open the source file
        myfits_hdu = pyfits.open(source)
        myfits_hdr = myfits_hdu[0].header
        
        # store scaling information for later use, as this is deleted when the array is updated
        BZERO = myfits_hdr['BZERO']
        BSCALE = myfits_hdr['BSCALE']
        BITPIX = myfits_hdr['BITPIX']
        bitpix_designation = pyfits.ImageHDU.NumCode[BITPIX]
        # print BSCALE, BZERO, BITPIX, bitpix_designation
        # read data array
        myfits_data = myfits_hdu[0].data
        myfits_hdu.close()
        
        """
        # determine read mode (TO BE COMPLETED !)
        if myfits_hdr['READMODE'] == 'multiple.endpoints':
            rmode = 'mer'
        elif myfits_hdr['READMODE'] == 'line.interlaced.read':
            rmode = 'lir'
            rmode = 'o2dcr' #only for testing
        else:
            rmode = 'o2dcr' # to be confirmed !
        
        # read fits file containing polynomial for correction
        polyfile = model.replace(".fits", "_"+rmode+".fits")
        """
        polyfile = model
        
        print 'Reading '+polyfile
        # for some weird reason I don't manage to read the cube using the standard functions
        # but have to use the convenience function getdata
        poly = pyfits.getdata(polyfile)
        
        
        ### read bad pixel file (currently not used and commented out)
        # badfile = 'badpixel_'+rmode+'.fits'
        # print 'Reading '+badfile
        # bad =z pyfits.getdata(badfile)
        ###
        
        # correct fits file
        #Note the factor 10^4 and the backward order of the polynomial coefficients.
        # you may want to constrain this only to "positive" corrections
        fac = 1 / (1 + poly[1,:,:] * myfits_data/1e4 + poly[0,:,:] * (myfits_data/1e4)**2)
        
        # save new fits file
        mfnp = source.partition('.fits')
        # add _lincor before .fits extension, or at the end if no such extension present
        outfitsname = mfnp[0] + '_lincor' + mfnp[1] + mfnp[2]
        # keep same data type
        myfits_hdu[0].data = (myfits_data * fac).astype(myfits_data.dtype)
        if os.path.exists(outfitsname):
            os.unlink(outfitsname)
            print 'Overwriting '+outfitsname
        else:
            print 'Writing '+outfitsname
        # scale back data to original values
        myfits_hdu[0].scale(bitpix_designation,'old')
        myfits_hdu.writeto(outfitsname, output_verify='ignore')
        myfits_hdu.close()

        log.info("Linearity correction done >> %s"%outfitsname)
        
        return outfitsname
        
################################################################################
# main
if __name__ == "__main__":
    
    
    usage = "usage: %prog [options]"
    parser = OptionParser(usage)
    
    # Basic inputs
    parser.add_option("-s", "--source",
                  action="store", dest="source_file_list",
                  help="Source file list of FITS files")
    
    parser.add_option("-o", "--out_data", type="str",
                  action="store", dest="out_data", default="source_corr.fits",
                  help="filename of out data file (default=%default)")
    
    parser.add_option("-c", "--coeff_file", type="str",
                  action="store", dest="out_coeff_file", 
                  help="filename of outputs coeffs, contains (a0, a1, a2)")
    
    parser.add_option("-a", "--apply_model",
                  action="store_true", dest="apply_model", default=False,
                  help="Apply to input the given non-linearity model [default],\
                  otherwise, the model will be computed with input files if they \
                  fit to requirements for model computation")
    
    ## not sure if required
    """
    parser.add_option("-l", "--limit",
                  action="store", dest="satur_lim", default=40000, 
                  help="saturation limit")
    
    parser.add_option("-t", "--ref_time",
                  action="store_true", dest="ref_time", default=False,
                  help="exptime used as reference")
    
    """
    parser.add_option("-v", "--verbose",
                  action="store_true", dest="verbose", default=False,
                  help="verbose mode [default]")
    
    (options, args) = parser.parse_args()
    
    #Check required parameters
    if (not options.source_file_list or not options.out_data 
        or not options.out_coeff_file  or len(args)!=0): # args is the leftover positional arguments after all options have been processed
        parser.print_help()
        parser.error("incorrect number of arguments " )
    if options.verbose:
        print "reading %s ..." % options.source_file_list
    
    # Two purposes, apply the model to correct non-linearity     
    if options.apply_model:
        NLM = NonLinearityModel()
        NLM.applyModel(options.source_file_list, options.out_coeff_file)
    # or, compute the non-linearity model
    else:
        filelist = [line.replace( "\n", "") for line in fileinput.input(options.source_file_list)]
        NLM = NonLinearityModel(filelist, "/tmp", options.out_coeff_file)
        NLM.createModel()    
        
