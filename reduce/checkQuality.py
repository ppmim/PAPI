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
# PANICtool
#
# checkQuality.py
#
# Created    : 07/11/2008    jmiguel@iaa.es
# Last update: 23/06/2009    jmiguel@iaa.es
#              17/09/2010    jmiguel@iaa.es - Added support to MEF files
#              22/10/2010    jmiguel@iaa.es - Updated the way to find out MEF files
#              17/01/2013    jmiguel@iaa.es - Modified the call to Sextractor,
#                                             now using astromatic package.
################################################################################
# Import necessary modules

from optparse import OptionParser
import numpy
import numpy.ma as ma 
import os
import pyfits
import sys


import misc.utils as utils
import astromatic

from misc.paLog import log

class CheckQuality(object):
    """
    Class used to estimate the image quality values using SExtractor 
    
    Parameters
    ----------
    data: 
        A input image
    
    Returns
    -------
       If no error, a seeing estimation value
    """
    def __init__(self, input_file, isomin=10.0, ellipmax=0.3, edge_x=200, edge_y=200, 
                 pixsize=0.45, gain = 4.15, sat_level=1500000, write=False):
        
        self.input_file = input_file
        # Default parameters values
        self.isomin = float(isomin)
        self.ellipmax = float(ellipmax)
        self.edge_x = int(edge_x)
        self.edge_y = int(edge_y)
        self.pixsize = float(pixsize)
        self.gain = float(gain)
        self.satur_level = sat_level
        self.write = False
        self.verbose = False
        self.MIN_NUMBER_GOOD_STARS = 5
        
    
    def estimateFWHM(self):
        """ 
        A FWHM of the current image is estimated using the 'best' stars on it.
        Generating an ascii text catalog with Sextractor, we can read the FWHM 
        values and give an estimation of the FWHM computing the median of the 
        'best' values/stars than fulfill some requirements, ie., ellipticity, 
        snr, location, etc. 
        
        It is very important that sextractor config file has the SATUR_LEVEL 
        parameter with a suitable value. In other case, we won't get any value 
        for FWHM. 
        
        SNR estimation as FLUX_AUTO/FLUXERR_AUTO or FLUX_APER/FLUXERR_APER
        
        Returns
        -------
        A couple of values (efwhm, std):
        
        efwhm : float
            Estimated FWHM.
        std: float
            Standard deviation of the FWHM.
        """

        # SExtractor configuration
        
        catalog_file = "test.cat"
        try:
            sex_cnf = os.environ['PAPI_HOME'] + "/config_files/sextractor.sex"
        except Exception,e:
            log.error("Error, variable PAPI_HOME not defined.")
            raise e
        
        sex = astromatic.SExtractor()
        sex.config['CONFIG_FILE']= sex_cnf
        #sex.config['PARAMETERS_NAME'] = os.environ['PAPI_HOME'] + "/irdr/src/config/default.param"
        sex.ext_config['CATALOG_TYPE'] = "ASCII"
        sex.ext_config['CHECKIMAGE_TYPE'] = "NONE"
        sex.ext_config['PIXEL_SCALE'] = self.pixsize
        sex.ext_config['GAIN'] = self.gain
        sex.ext_config['SATUR_LEVEL'] = self.satur_level
        sex.ext_config['CATALOG_NAME'] = catalog_file
        #sex.ext_config['DETECT_THRESH'] = config_dict['astrometry']['mask_thresh']
        #sex.ext_config['DETECT_MINAREA'] = config_dict['astrometry']['mask_minarea']
        
        # SExtractor execution
        try:
            sex.run(self.input_file, updateconfig=False, clean=False)
        except Exception,e:
            log.error("Error running SExtractor: %s"%str(e))  
            raise e
        
        ## SExtractor Catalog columns required and expected (sextractor.param)
        #0  NUMBER
        #1  X_IMAGE
        #2  Y_IMAGE
        #3  ALPHA_J2000
        #4  DELTA_J2000
        #5  MAG_BEST
        #6  ISOAREA_IMAGE
        #7  ELLIPTICITY
        #8  FWHM_IMAGE
        #9  FLUX_RADIUS
        #10 FLUX_AUTO
        #11 FLUXERR_AUTO
        #12 FLAGS
        #13 FLUX_APER
        #14 FLUXERR_APER
        
        source_file = catalog_file

        try:
            if self.write: fits_file = pyfits.open(self.input_file, 'update')
            else: fits_file = pyfits.open(self.input_file, 'readonly')
        except Exception,e:
            log.error("Error while openning file %s",self.input_file)
            raise e
        
        try:
            if len(fits_file)>1: # is a MEF
                naxis1 = fits_file[1].header['NAXIS1']
                naxis2 = fits_file[1].header['NAXIS2']
            else:  # is a simple FITS      
                naxis1 = fits_file[0].header['NAXIS1']
                naxis2 = fits_file[0].header['NAXIS2']
        except KeyError,e:
            log.error("Error while reading FITS header NAXIS keywords :%s",str(e))
            raise Exception("Error while reading FITS header NAXIS keywords")
            
        
        # Now, read the SEx catalog
        #fwhm_world=[float(line.split()[7]) for line in fileinput.input(source_file)]
        #matrix=[line.split() for line in fileinput.input(source_file)]
        #b=numpy.array(matrix)
        
        a = numpy.loadtxt(source_file)
        
        good_stars = []
        # Select 'best' stars for the estimation
        std = numpy.std(a[:,8])
        print "Initial STD of FWHM=",std
        for i in range(0, a.shape[0]):
            x = a[i,1]
            y = a[i,2]
            isoarea = a[i,6]
            ellipticity = a[i,7]
            fwhm = a[i,8]
            flux = a[i,10]
            flux_err = a[i,11]
            flags = a[i,12]
            #fa=a[i,13]
            #fea=a[i,14]
            if flux_err!=0:
                snr = flux/flux_err
                #print "SNR=",snr
            else:
                continue
            if (x>self.edge_x and x<naxis1-self.edge_x and 
                y>self.edge_y and y<naxis2-self.edge_y and
                ellipticity<self.ellipmax and fwhm>0.1 and 
                fwhm<20 and flags==0 and   
                isoarea>float(self.isomin) and snr>20.0): 
                # and fwhm<5*std it does not work many times
                good_stars.append(a[i,:])
                #print "%s SNR_APER= %s " %(i, snr)
            else:
                """
                print "START #%s"%i
                print "  SNR=",snr
                print "  FWHM=",fwhm
                print "  AREA=",isoarea
                print "  ELLIP=", ellipticity
                """
        
        m_good_stars = numpy.array(good_stars)
        
        print "Found <%d> GOOD stars"%len(m_good_stars)
        
        if len(m_good_stars)>self.MIN_NUMBER_GOOD_STARS:
            std = numpy.std(m_good_stars[:,8])
            print "STD2 = ",std
            efwhm = numpy.median(m_good_stars[:,8])
            print "best-FWHM-median(pixels) = ", efwhm
            print "FLUX_RADIUS (as mentioned in Terapix T0004 explanatory table) =", numpy.median(m_good_stars[:,9])
            print "Masked-mean = ", ma.masked_outside(m_good_stars[:,8], 0.01, 3*std).mean()
            
            if self.write:
                fits_file[0].header.update('hierarch PAPI.SEEING', efwhm*self.pixsize)
                print "Fits keyword updated " 
            
        else:
            print "Not enough stars found !!"
            fits_file.close(output_verify='ignore')    
            return -1,-1
        
        #print "FWHM-median(pixels)= ", numpy.median(fwhm_world), numpy.amin(fwhm_world), numpy.amax(fwhm_world)
        #print "FWHM-mean(pixels)= ", numpy.mean(fwhm_world)
        
        fits_file.close(output_verify='ignore')
        
        # cleanup files
        os.unlink(catalog_file)
        
        return efwhm,std
      
    def estimateBackground(self, output_file):
        """ 
        Runs SExtractor to estimate the image background.
        
        Parameters
        ----------
        output_file: str
            Filename of the image background to be generated.
            
        Returns
        ------- 
        The image background file obtained by SEextractor,ie, output_file. 
        """
        
        # Sextractor config
        try:
            sex_cnf = os.environ['PAPI_HOME'] + "/config_files/sextractor.sex"
        except Exception,e:
            log.error("Error, variable PAPI_HOME not defined.")
            raise e
        sex = astromatic.SExtractor()
        #sex.config['CONFIG_FILE']= "/disk-a/caha/panic/DEVELOP/PIPELINE/PANIC/trunk/config_files/sex.conf"
        sex.config['CONFIG_FILE']= sex_cnf
        #sex.config['CATALOG_TYPE'] = "ASCII"
        sex.ext_config['CHECKIMAGE_TYPE'] = "BACKGROUND"
        sex.ext_config['PIXEL_SCALE'] = self.pixsize
        sex.ext_config['GAIN'] = self.gain
        sex.ext_config['SATUR_LEVEL'] = self.satur_level
        sex.ext_config['CHECKIMAGE_NAME'] = output_file
        
        # SExtractor execution
        try:
            sex.run(self.input_file, updateconfig=False, clean=False)
        except Exception,e:
            log.error("Error running SExtractor: %s"%str(e))  
            raise e
        else:
            return output_file     

            

################################################################################
# main
if __name__ == "__main__":
    
    # Get and check command-line options
        
    usage = "usage: %prog [options]"
    desc = """This module gives an estimation of the FWHM of the input image 
using best stars of its SExtractor catalog
"""
 
    parser = OptionParser(usage, description=desc)
    
    parser.add_option("-f", "--input_image",
                  action="store", dest="input_image", 
                  help="Input image to computer FWHM estimation.")
                  
    parser.add_option("-i", "--isoarea_min",
                  action="store", dest="isoarea_min",type=int,
                  help="Minimum value of ISOAREA (default = %default)",
                  default=10)
    
    parser.add_option("-S", "--snr",
                  action="store", dest="snr", type=int,
                  help="Min SNR of stars to use (default = %default)",
                  default=10)
    
    parser.add_option("-e", "--ellipmax",
                  action="store", dest="ellipmax", type=float, default=0.2,
                  help="Maximum SExtractor ELLIPTICITY (default = %default)")
                  
    parser.add_option("-x", "--edge_x",
                  action="store", dest="edge_x", type=int, 
                  help="Consider sources out of image borders on X axis (default = %default)",
                  default=200)
    
    parser.add_option("-y", "--edge_y",
                  action="store", dest="edge_y", type=int, 
                  help="Consider sources out of image borders on Y axis (default = %default)",
                  default=200)
    
    parser.add_option("-p", "--pixsize",
                  action="store", dest="pixsize", type=float, 
                  help="Pixel scale of the input image (default = %default)",
                  default=0.45)
    
    parser.add_option("-g", "--gain",
                  action="store", dest="gain", type=float, 
                  help="Detector gain (default = %default)",
                  default=4.15)
    
    parser.add_option("-l", "--satur_level",
                  action="store", dest="satur_level", type=float, 
                  help="Saturation level (default = %default)",
                  default=1500000)
    
    parser.add_option("-w", "--write",
                  action="store_true", dest="write", default=True,
                  help="Update header with PA_SEEING keyword [default=%default]")
                                
    (options, args) = parser.parse_args()
    
    if len(sys.argv[1:])<1:
       parser.print_help()
       sys.exit(0)
    
    if not options.input_image or len(args)!=0: 
    # args is the leftover positional arguments after all options have been processed
        parser.print_help()
        parser.error("wrong number of arguments " )
    

    if not os.path.exists(options.input_image):
        log.error ("Input image %s does not exist", options.input_image)
        sys.exit(0)
        
    try:
        cq = CheckQuality(options.input_image, options.isoarea_min, 
                          options.ellipmax, options.edge_x, options.edge_y, 
                          options.pixsize, options.gain, options.satur_level , 
                          options.write)
        cq.estimateFWHM()
    except Exception,e:
        log.error("There was some error: %s "%str(e))
        sys.exit(0)
        
    
