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
# TODO : computation of non-linearity based of a set of dard subtracted flats,
#        so, dark subtrction must be done previously.
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
import misc.robust as robust

# Interact with FITS files
import pyfits
import numpy

# Logging
from misc.paLog import log

class NonLinearityModel(object):
    """
    Class used to compute the Non-linearity model for the detectors
    Based on the algorithm described by T.H.Jarrett et all (1994) we compute 
    the non-linearity coefficientes of the model.
    """
    def __init__(self, input_files=None, output_filename="/tmp/NLC.fits"):
        """
        Init the object.
        
        Parameters
        ----------
        
        input_data: list
            A list of sky files
        input_files: list
            A list of FITS files as input for the computation
        output_filename: str
            File where coefficientes will be saved
        
        
        coeff: array
            Array of values a0,a1,a2 where 
            a0 = bias level
            a1 = quantum efficiency or sensitivity
            a2 = the non-linearity or saturation
        """
        self.__input_files = input_files
        self.__output_filename = output_filename  # full filename (path+filename)
    
    def createModel(self):
          
        """
        NOT COMPLETED AT ALL !!!
        ************************

        Compute the Non-linearity correction model based on a serie of FITS files
        taken in order to get check the detector linearity.
        
        1) Determine the coefficients:

          1.1 Take a series of "flat" exposures with the same illumination but
              increasing exp.time
          1.2 Calculate a polynomial fit (deg=3) to these data pixel by pixel (dark subtracted!)
          1.3 Store the coefficients of these polynomials in frames, so you have
              a frame/plane for each coeff:

                plane_0 = coeff_0 (intercept, bias level)
                plane_1 = coeff_1 (quantum efficiency or sensitivity)
                plane_2 = coeff_2 (the non-linearity or saturation)
                plane_3 = coeff_3 (???)


        2) Apply the corrections:

          2.1 For each pixel in the frame to be corrected compute the polynomial, 
              add the difference between poly. and linear relation to measured flux.

        Parameters
        ----------
        
        Returns
        -------
        The filename of the file created having a cube N=4 planes, one for each
        coeff of the polynomial fitted. 

        
        """   
        
        log.debug("Start createModel")
        
        start_time = time.time()
        t = utils.clock()
        t.tic()
        
        # Get the user-defined list of flat frames
        framelist = self.__input_files
        
        # STEP 0: Determine the number of flat frames to combine
        try:    
            nframes = len(framelist)
        except IndexError:
            log.error("No FLAT frames defined")
            raise
        
        if not os.path.exists(os.path.dirname(self.__output_filename)):
            raise NameError, 'Wrong output path'
        if not self.__output_filename:
            log.error("Output FLAT frame not defined")
            raise Exception("Wrong output filename")
    
        # Change to the source directory
        base, infile = os.path.split(self.__output_filename) 
        
        flats = numpy.zeros(nframes, dtype=numpy.int)
         
        # STEP 1: Check TYPE(flat), READMODE and read the EXPTIME of each frame
        # print "FRAMELIST= %s" %framelist
        i = 0
        f_readmode = -1
        f_n_extensions = -1
        for iframe in framelist:
            fits = datahandler.ClFits(iframe)
            log.debug("Frame %s EXPTIME= %f TYPE= %s " %(iframe, fits.exptime, fits.type)) 
            # Check TYPE (flat)
            if not fits.isDomeFlat():
                log.warning("Warning: Task 'createDarkModel' found a non DomeFlat frame. Skipping %s", iframe)
                flats[i] = 0
            else:
                # Check READMODE
                if ( f_readmode!=-1 and (f_readmode!= fits.getReadMode() )):
                    log.error("Error: Task 'createNLModel' finished. Found a FLAT frame with different  READMODE")
                    flats[i] = 0  
                    #continue
                    raise Exception("Found a FLAT frame with different READMODE") 
                else: 
                    f_readmode = fits.getReadMode()
                    f_n_extensions = fits.getNExt()
                    log.debug("NEXT= %s"%(f_n_extensions))
                    flats[i] = 1
                
            i = i+1
        log.debug('All frames checked')   
        
        naxis1 = fits.naxis1
        naxis2 = fits.naxis2            
        nflats = (flats==1).sum()
        
        if nflats<10:
            log.error('Input dataset doesnt have enough frames. At least 10 Flat frames are needed')
            raise Exception("Flat sequence is too short. Al least 10 Flat frames are needed !")
        
        # Initialize some storage arrays
        times = numpy.zeros(nflats, dtype=numpy.float32)
        temp = numpy.zeros([nflats, f_n_extensions, naxis1, naxis2], dtype=numpy.float32)
        

        # loop the images
        counter = 0
        for i in range(0, nframes):
            file = pyfits.open(framelist[i])
            f = datahandler.ClFits ( framelist[i] )
            if flats[i]==1:
                for i_ext in range(0, f_n_extensions):
                    if f_n_extensions==1:
                        temp[counter, 0, :,:] = file[0].data
                    else:
                        log.debug("Found MEF file")
                        temp[counter, i_ext, :,:] = file[i_ext+1].data
                _mean = numpy.mean(temp[counter])
                _robust_mean = robust.mean(temp[counter].reshape(naxis1*naxis2*f_n_extensions))
                _median = numpy.median(temp[counter])
                _mode = 3*_median - 2*_mean
                log.info("Flat frame TEXP=%s , ITIME=%s ,MEAN_VALUE=%s , MEDIAN=%s ROBUST_MEAN=%s"%(f.expTime(), f.getItime(), _mean, _median, _robust_mean))
                times[counter] = float(f.expTime())
                counter = counter+1
                file.close()

        log.debug("Now fitting the linearity model...")
        # now collapse and fit the data
        # polyfit returns polynomial coefficients ordered from low to high.
        pol_degree = 3 # third-order polynomial, ie., 4 coeffs
        fit = numpy.polynomial.polynomial.polyfit(times, 
                            temp.reshape(len(times), naxis1*naxis2*f_n_extensions), deg=pol_degree)

        # Get the median value of the coeffs                
        coeff_0 = numpy.median(fit[0])
        coeff_1 = numpy.median(fit[1])
        coeff_2 = numpy.median(fit[2])
        coeff_3 = numpy.median(fit[3])

        log.info("Coeff_0 = %s"%coeff_0)
        log.info("Coeff_1 = %s"%coeff_1)
        log.info("Coeff_2 = %s"%coeff_2)
        log.info("Coeff_3 = %s"%coeff_3)
        
        misc.fileUtils.removefiles(self.__output_filename)               

        # Write result in a FITS
        hdulist = pyfits.HDUList()
        hdr0 = pyfits.getheader(framelist[numpy.where(flats==1)[0][0]])
        prihdu = pyfits.PrimaryHDU (data = None, header = None)
        try:
            prihdu.header.set('INSTRUME', hdr0['INSTRUME'])
            prihdu.header.set('TELESCOP', hdr0['TELESCOP'])
            prihdu.header.set('CAMERA', hdr0['CAMERA'])
            prihdu.header.set('MJD-OBS', hdr0['MJD-OBS'])
            prihdu.header.set('DATE-OBS', hdr0['DATE-OBS'])
            prihdu.header.set('DATE', hdr0['DATE'])
            prihdu.header.set('UT', hdr0['UT'])
            prihdu.header.set('LST', hdr0['LST'])
            prihdu.header.set('ORIGIN', hdr0['ORIGIN'])
            prihdu.header.set('OBSERVER', hdr0['OBSERVER'])
        except Exception,e:
            log.warning("%s"%str(e))

        prihdu.header.set('PAPITYPE','MASTER_LIN_MODEL')
        prihdu.header.add_history('Linearity model based on %s' % framelist)
        
        if f_n_extensions>1:
            prihdu.header.set('EXTEND', pyfits.TRUE, after = 'NAXIS')
            prihdu.header.set('NEXTEND', f_n_extensions)
            prihdu.header.set('FILENAME', self.__output_filename)
            hdulist.append(prihdu)
            for i_ext in range(0, f_n_extensions):
                hdu = pyfits.PrimaryHDU()
                hdu.scale('float32') # important to set first data type
                hdu.data = fit.reshape(pol_degree+1, f_n_extensions, naxis1, naxis2)[:, i_ext, :, :]
                hdulist.append(hdu)
                del hdu
        else:
            prihdu.scale('float32') # important to set first data type
            prihdu.data = fit.reshape(pol_degree+1, 1, naxis1, naxis2)[:, 0, :, :]
            hdulist.append(prihdu)
         
        
        # write FITS
        try:
            hdulist.writeto(self.__output_filename)
            hdulist.close(output_verify='ignore')
        except Exception,e:
            log.error("Error writing linearity model %s"%self.__output_filename)
            raise e
        
        log.debug('Saved Linearity Model to %s' , self.__output_filename)
        log.debug("createModel' finished %s", t.tac() )

        
        return self.__output_filename
    

    def applyModel(self, source, model, suffix='_LC', out_dir='/tmp'):
        """
        Do the Non-linearity correction using the supplied model. In principle,
        it should be applied to all raw images (darks, flats, science, ...).
        
        Parameters
        ----------
        source : str
            List of FITS file names to be corrected.
        
        model : str 
            FITS filename of the Non-Linearity model, ie., containing polynomial 
            coeffs (3rd order) for correction that must has been previously 
            computed. It must be a cube with 4 planes (a plane for each coeff), 
            and N extensions (one per detector). Planes definitions:

                plane_0 = coeff_0 (intercept, bias level)
                plane_1 = coeff_1 (quantum efficiency or sensitivity)
                plane_2 = coeff_2 (the non-linearity or saturation)
                plane_3 = coeff_3 (?)
        
        suffix: str
            Suffix to add to the input filename to generate the output filename.


        out_dir: str
            Directory where new corredted files will be saved

        Returns
        -------
        outfitsname: list
            The list of new corrected files created.
        
        TODO
        ----
        - adjust the read modes availables on PANIC (lir, mer, o2dcr, ...)
        - read the proper BadPixelMask ??
        - compute properly the correction to do !!! 
        - Proccess in parallel each extension.
        
        """   
        log.debug("Start applyModel")


        if len(source)<1:
            log.error("Found empty list of files")
            raise Exception("Found empty list of files")
        
        if not os.path.exists(model):
            log.error("Cannot read non-linearity model file %s"%model)
            raise Exception("Cannot read non-linearity model file")


        # Get number of planes, ie. number of coeffs of the polynomial
        # We suppose a 3rd degree (4 coeffs) polynomial fit of the linearity model 
        fits_model = pyfits.open(model)
        model_n_extensions = 1 if len(fits_model)==1 else len(fits_model)-1
        if model_n_extensions==1:
            data_model = fits_model[0].data
        else:
            data_model = fits_model[1].data

        #if data_model.shape[0]!=4:
        if data_model.shape[0]!=3:
            log.error("Linearity model does not match a 4-plane cube image")
            raise Exception("Linearity model does not match a 4-plane cube image")

        # loop the images
        new_filenames = []
        for i in range(0, len(source)):
            i_file = pyfits.open(source[i])
            f_n_extensions = 1 if len(i_file)==1 else len(i_file)-1
            
            log.debug("Raw Number of extensions = %s"%f_n_extensions)
            log.debug("Model Number of extensions = %s"%model_n_extensions)
            
            if f_n_extensions!=model_n_extensions:
                log.error("Model and Raw source do not match number of extensions")
                raise Exception("Model and Raw source do not match number of extensions")
            
            for i_ext in range(0, f_n_extensions):
                offset = 0 if f_n_extensions==1 else 1
                data_model = fits_model[i_ext+offset].data.astype(float)
                raw = i_file[i_ext+offset].data.astype(float)
                
                if not raw.shape==data_model[0].shape:
                    log.error("Shape/size of lin_model and source_data does not match")
                    raise Exception("Shape/size of lin_model and source_data does not match")
                
                log.info("Median value before correction = %s"%(numpy.median(raw)))
                corr = data_model[0] + data_model[1]*raw + data_model[2]*raw**2
                #corr = data_model[0] + data_model[1]*raw + data_model[2]*raw**2 + data_model[3]*raw**3
                # For simulation/test
                # Note: skyfilter could fail if corrected image has wrong values 
                # (e.g.: negative mean, very big/small float values)
                #corr = data_model[0] + 0.000004*raw + 0.0000001*raw**2 + 0.000000000001*raw**3
                diff = corr - raw
                log.info("Median value of difference = %s"%(numpy.median(diff)))
                i_file[i_ext+offset].data = corr
                log.info("Median value after correction = %s"%(numpy.median(i_file[i_ext+offset].data)))
            
            # Write output to outframe (data object actually still points 
            # to input data)
            # save new fits file
            mfnp = os.path.basename(source[i]).partition('.fits')
            # add suffix before .fits extension, or at the end if no such extension present
            outfitsname = out_dir + '/' + mfnp[0] + suffix + mfnp[1] + mfnp[2]
            outfitsname = os.path.normpath(outfitsname)

            # keep same data type
            if os.path.exists(outfitsname):
                os.unlink(outfitsname)
                print 'Overwriting ' + outfitsname
            else:
                print 'Writing ' + outfitsname
            
            try:
                i_file[0].scale('float32')
                i_file[0].header.add_history('Applied Non-linearity correction with %s'%model)
                i_file.writeto(outfitsname, output_verify='ignore')
            except IOError:
                raise ExError('Cannot write output to %s' % outfitsname)
            i_file.close()

            new_filenames.append(outfitsname)
           
        return new_filenames 

    def applyModel_LBT(self, source, model):
          
        """
        Do the Non-linearity correction using the supplied model
        
        Parameters
        ----------
        
        source : str
            FITS filename to be corrected
        
        model : str 
            FITS filename of the Non-Linearity model, ie., containing polynomial 
            coeffs for correction that must has been previously computed.
        
        
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
        
        # Correct fits file
        # Note the factor 10^4 and the backward order of the polynomial coefficients.
        # you may want to constrain this only to "positive" corrections
        fac = 1 / (1 + poly[1,:,:] * myfits_data/1e4 + poly[0,:,:] * (myfits_data/1e4)**2)
        
        # save new fits file
        mfnp = source.partition('.fits')
        # add suffix before .fits extension, or at the end if no such extension present
        outfitsname = mfnp[0] + suffix + mfnp[1] + mfnp[2]
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
    desc= """Compute the non-linearity of the detector from a set of darks and
flats frames.
Optionaly an already computed non-linearity model is applied.

NOTE: Not yet implemented !!
"""
    parser = OptionParser(usage, description=desc)
    
    # Basic inputs
    parser.add_option("-s", "--source",
                  action="store", dest="source_file_list",
                  help="Source file list of FITS files")
    
    parser.add_option("-o", "--out_data", type="str",
                  action="store", dest="out_data", default="source_corr.fits",
                  help="filename of out data file (default=%default)")
    
    parser.add_option("-c", "--coeff_file", type="str",
                  action="store", dest="out_coeff_file", 
                  help="FITS cube file of outputs coeffs, contains (a0, a1, a2, a3) planes")
    
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
    
    if len(sys.argv[1:])<1:
       parser.print_help()
       sys.exit(0)

    # Check required parameters
    if (not options.source_file_list or not options.out_data 
        or not options.out_coeff_file  or len(args)!=0): # args is the leftover positional arguments after all options have been processed
        parser.print_help()
        parser.error("incorrect number of arguments " )
    if options.verbose:
        print "reading %s ..." % options.source_file_list
    
    # Two purposes, apply the model to correct non-linearity     
    filelist = [line.replace( "\n", "") for line in fileinput.input(options.source_file_list)]
    if options.apply_model:
        NLM = NonLinearityModel(filelist)
        NLM.applyModel(filelist, options.out_coeff_file)
    # or, compute the non-linearity model
    else:
        NLM = NonLinearityModel(filelist, options.out_coeff_file)
        NLM.createModel()    
        
