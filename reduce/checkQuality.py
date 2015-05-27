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
#              20/07/2014    jmiguel@iaa.es - 
################################################################################
# Import necessary modules

from optparse import OptionParser
import fileinput
import numpy
import numpy.ma as ma 
import os
import astropy.io.fits as fits
import sys

from pyraf import iraf
import re

import misc.utils as utils
import astromatic

from misc.paLog import log

class CheckQuality(object):
    """
    Class used to estimate the image quality values using SExtractor.
    Although MEF files are supported, the routine does not distinguish the
    values of each extesion/detector, what would be a good thing.

    
    Parameters
    ----------
    data: 
        A input image
    
    Returns
    -------
       If no error, a seeing estimation value
    """
    def __init__(self, input_file, isomin=32.0, ellipmax=0.3, edge_x=2, edge_y=2, 
                 pixsize=0.45, gain = 4.15, sat_level=1500000, write=False,
                 min_snr=5.0, window='all'):
        
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
        self.min_snr = min_snr
        self.MIN_NUMBER_GOOD_STARS = 0
        self.window = window


        if self.window == 'Q1':
            self.sex_input_file = input_file + '[%d]'%1 # ext1
        elif self.window == 'Q2':
            self.sex_input_file = input_file + '[%d]'%2 # ext2
        elif self.window == 'Q3':
            self.sex_input_file = input_file + '[%d]'%3 # ext3
        elif self.window == 'Q4':
            self.sex_input_file = input_file + '[%d]'%4 # ext4
        else:
            self.sex_input_file = input_file
             
        
    
    def estimateFWHM(self, psfmeasure=False):
        """ 
        A FWHM of the current image is estimated using the 'best' stars on it.
        Generating an ascii text catalog with Sextractor, we can read the FWHM 
        values and give an estimation of the FWHM computing the median of the 
        'best' values/stars than fulfill some requirements, ie., ellipticity, 
        snr, location, etc. 
        
        It is very important that sextractor config file has the SATUR_LEVEL 
        parameter with a suitable value. In other case, we won't get any value 
        for FWHM. 
        
        SNR estimation as FLUX_AUTO/FLUXERR_AUTO or FLUX_APER/FLUXERR_APER,
        that is, the signal divided by the noise.
        
        Returns
        -------
        The values (efwhm, std, x, y):
        
        efwhm : float
            Estimated FWHM (in pixels)
        std: float
            Standard deviation of the FWHM.
        x,y: coordinates of the last star found by sextractor. It is thought
        for imagenes (subwindows) with a single star on the field.

        """
        
        # Check whether detector selection can be done
        if self.window != 'all':
            f = fits.open(self.input_file)
            if len(f) != 5:
                raise Exception("Error, expected a MEF file with 4 extensions")

        
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
            sex.run(self.sex_input_file, updateconfig=False, clean=False)
        except Exception,e:
            log.error("Error running SExtractor: %s"%str(e))  
            raise e
        
        ## SExtractor Catalog columns required and expected (sextractor.param)
        # 0 NUMBER           # Running object number
        # --------- Position Parameters --------------
        # 1 X_IMAGE          # Object position along x
        # 2 Y_IMAGE          # Object position along y

        # 3 ALPHA_J2000      # Right ascension of barycenter (J2000) [deg]
        # 4 DELTA_J2000      # Declination of barycenter (J2000) [deg]


        # --------- Photometric Parameters -----------
        # 5 MAG_BEST         # Best of MAG_AUTO and MAG_ISOCOR
        # 6 ISOAREA_IMAGE    # Isophotal area above Analysis threshold [pixel**2]
        # 7 ELLIPTICITY      # 1 - B_IMAGE/A_IMAGE
        # 8 FWHM_IMAGE       # FWHM assuming a gaussian core
        # 9 FLUX_RADIUS      # Fraction-of-light radii 
        # 10 FLUX_APER        # Flux vector within fixed circular aperture(s)
        # 11 FLUXERR_APER     # RMS error vector for aperture flux(es) 

        # --------- Flags  -----------
        # 12 FLAGS            # Extraction flags

        # -- Test (not read by PAPI) --
        # 13 FLUX_AUTO        # Flux within a Kron-like elliptical aperture [counts] 
        # 14 FLUXERR_AUTO     # RMS error for AUTO flux [counts]
        # 15 XWIN_IMAGE       # Windowed position estimate along x [pix]
        # 16 YWIN_IMAGE       # Windowed position estimate along y [pix]
      
        source_file = catalog_file

        try:
            if self.write: fits_file = fits.open(self.input_file, 'update')
            else: fits_file = fits.open(self.input_file, 'readonly')
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
        
        a = numpy.loadtxt(source_file, ndmin=2)

        if len(a)==0:
            raise Exception("Empy catalog, No stars found.")
        

        good_stars = []
        # Select 'best' stars for the estimation
        std = numpy.std(a[:,8])
        print "Initial STD of FWHM=",std
        for i in range(0, a.shape[0]):
            x = a[i,1]
            y = a[i,2]
            xwin = a[i,15]
            ywin = a[i,16]
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
            if (x > self.edge_x and x < naxis1 - self.edge_x and 
                y > self.edge_y and y < naxis2 - self.edge_y and
                ellipticity < self.ellipmax and fwhm > 0.1 and 
                fwhm < 20 and flags <= 31 and   
                isoarea > float(self.isomin) and snr > self.min_snr): 
                # and fwhm<5*std it does not work many times
                good_stars.append(a[i,:])
                #print "%s SNR_APER= %s " %(i, snr)
                #print "ISO_AREA= %s  ISO_MIN=%s"%(isoarea,self.isomin)
            else:
                """
                print "STAR #%s"%i
                print "X=%s  Y=%s"%(x,y)
                print "  SNR=",snr
                print "  FWHM=",fwhm
                print "  AREA=",isoarea
                print "  ELLIP=", ellipticity
                """
                pass
        
        m_good_stars = numpy.array(good_stars)
        
        print "Found <%d> GOOD stars" % len(m_good_stars)
        
        if len(m_good_stars)>self.MIN_NUMBER_GOOD_STARS:
            std = numpy.std(m_good_stars[:,8])
            print "STD2 = ",std
            efwhm = numpy.median(m_good_stars[:,8])
            print "best-FWHM-median(pixels) = ", efwhm
            print "Mean-FWHM(px) = ", numpy.mean(m_good_stars[:,8])
            print "FLUX_RADIUS (as mentioned in Terapix T0004 explanatory table) =", numpy.median(m_good_stars[:,9])
            print "Masked-mean = ", ma.masked_outside(m_good_stars[:,8], 0.01, 3*std).mean()
            
            if self.write:
                fits_file[0].header.update('hierarch PAPI.SEEING', efwhm*self.pixsize)
                print "Fits keyword updated " 

            # 2nd Estimation Method (psfmeasure)
            if psfmeasure:
                coord_text_file = "/tmp/coord_file.txt"
                # Build the coord_file
                with open(coord_text_file, "w") as text_file:
                    for source in m_good_stars:
                        text_file.write("%s   %s\n"%(source[15], source[16]))
                try:        
                    pfwhm = self.getAverageFWHMfromPsfmeasure(self.input_file, coord_text_file)
                    log.debug("Average FWHM (psfmeasure-Moffat): %s"%pfwhm)
                except Exception,e:
                    log.error("Cannot run properly iraf.psfmeasure")
                    log.error("%s"%str(e))
        else:
            print "Not enough good stars found !!"
            fits_file.close(output_verify='ignore')    
            return -1,-1,-1,-1
        
        #print "FWHM-median(pixels)= ", numpy.median(fwhm_world), numpy.amin(fwhm_world), numpy.amax(fwhm_world)
        #print "FWHM-mean(pixels)= ", numpy.mean(fwhm_world)
        
        fits_file.close(output_verify='ignore')
        
        # cleanup files
        # os.unlink(catalog_file)
        
        return efwhm, std, xwin, ywin

    
    def getAverageFWHMfromPsfmeasure(self, image, coord_file):
        """
        Calculate the average Full Width Half Max for the objects in image
        at the coords specified in coord_file calling iraf.obsutil.psfmeasure.

        The coordinates in coord_file should be in the same world coordiantes
        as the WCS applied to the image.
        
        Exam. of coord_file:
        1024  1024
        1100  1200
        ...

        Note: iraf.obsutil.psfmeasure does not work with multiprocessing, as it 
        opens a window with the graphic representation of the measurements.
        So, in PAPI cannot be used with parallel reduction enabled.

        """

        iraf.noao(_doprint=0)
        iraf.obsutil(_doprint=0)
        iraf.unlearn("psfmeasure")
      
        psfmeasure = iraf.obsutil.psfmeasure
        # setup all paramaters
        psfmeasure.coords = "mark1"
        # 'world' gives problems with SIP headers (ie. Astrometry.net)
        # 'phycal' gives too many 'Warning: Invalid flux profile ' 
        psfmeasure.wcs = "logical" 
        psfmeasure.display = "no"
        psfmeasure.size = "GFWHM"  # Moffat profile
        psfmeasure.radius = 5
        psfmeasure.iterations = 2
        psfmeasure.imagecur = coord_file
        psfmeasure.graphcur = '/dev/null' #file that is empty by definition
        psfmeasure.logfile = "fwhm.log"
        res = psfmeasure(image, Stdout=1)[-1] # get last linet of output
        numMatch = re.compile(r'(\d+(\.\d+)?)')
        match = numMatch.search(res)

        return float (match.group(1))
    
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
using best stars of its SExtractor catalog. If input file is a MEF, the routine
gives a mean estimation of the FWHM of all extensions/detectors, not distinguishing
between them. However, the -W --window flag can be used to specify a certain
detector (Q1,Q2,Q3,Q4). 
"""
 
    parser = OptionParser(usage, description=desc)
    
    parser.add_option("-f", "--input_image",
                  action="store", dest="input_image", 
                  help="Input FITS image to run FWHM estimation.")

    parser.add_option("-F", "--input_file_list",
                  action="store", dest="input_file_list", 
                  help="Text file with all FITS images to be FWHM-estimated. "
                  "It will produce a Output.txt text file with the results.")

    parser.add_option("-i", "--isoarea_min",
                  action="store", dest="isoarea_min",type=int,
                  help="Minimum value of ISOAREA (default = %default)",
                  default=32)
    
    parser.add_option("-S", "--snr",
                  action="store", dest="snr", type=int,
                  help="Min SNR of stars to use (default = %default)",
                  default=5)
    
    parser.add_option("-e", "--ellipmax",
                  action="store", dest="ellipmax", type=float, default=0.3,
                  help="Maximum SExtractor ELLIPTICITY (default = %default)")
                  
    parser.add_option("-x", "--edge_x",
                  action="store", dest="edge_x", type=int, 
                  help="Consider sources out of image borders on X axis (default = %default)",
                  default=2)
    
    parser.add_option("-y", "--edge_y",
                  action="store", dest="edge_y", type=int, 
                  help="Consider sources out of image borders on Y axis (default = %default)",
                  default=2)
    
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
    
    parser.add_option('-W', '--window',
                      type='choice',
                      action='store',
                      dest='window',
                      choices=['Q1', 'Q2', 'Q3', 'Q4', 'all'],
                      default='all',
                      help="When input is a MEF, it means the "
                      "window/dectector/extension to process: "
                      "Q1, Q2, Q3, Q4, full [default: %default]")

    parser.add_option("-P", "--psfmeasure",
                  action="store_true", dest="psfmeasure", default=False,
                  help="Show iraf.obsutil.psfmeasure FWHM measurements [default=%default]")

    (options, args) = parser.parse_args()
    
    if len(sys.argv[1:])<1:
       parser.print_help()
       sys.exit(0)
    
    if (not options.input_image and not options.input_file_list) or len(args)!=0: 
    # args is the leftover positional arguments after all options have been processed
        parser.print_help()
        parser.error("Wrong number of arguments " )
    

    #if not os.path.exists(options.input_image) and not os.path.exists(options.input_file_list):
    #    log.error ("Input image %s does not exist", options.input_image)
    #    sys.exit(0)
    
    if options.input_file_list:
        if not os.path.exists(options.input_file_list):
            log.error ("Input image %s does not exist", options.input_file_list)
            sys.exit(0)

        filelist = [line.replace( "\n", "") 
                for line in fileinput.input(options.input_file_list)]
        
        text_file = open("Output.txt", "w")
        text_file.write("#  Filename \t FWHM \t STD \t X \t Y\n")
        text_file.write("#  X,Y =spatial coordinates of **last** star found. \n")

        for m_file in filelist:
            try:
                cq = CheckQuality(m_file, options.isoarea_min, 
                                  options.ellipmax, options.edge_x, options.edge_y, 
                                  options.pixsize, options.gain, options.satur_level , 
                                  options.write, options.snr, options.window)
                efwhm, std, x, y = cq.estimateFWHM()
                text_file.write("%s    %s    %s    %s    %s\n"%(m_file, efwhm, std, x, y))
            except Exception,e:
                log.error("There was some error with image %s : %s "%(m_file, str(e)))
                text_file.write("%s    %s    %s\n"%(m_file, -1, -1, -1, -1))
        
        text_file.close()
    
    elif options.input_image:
        if not os.path.exists(options.input_image):
            log.error ("Input image %s does not exist", options.input_image)
            sys.exit(0)

        try:
            cq = CheckQuality(options.input_image, options.isoarea_min, 
                                options.ellipmax, options.edge_x, options.edge_y, 
                                options.pixsize, options.gain, options.satur_level , 
                                options.write, options.snr, options.window)
            efwhm, std, k, k = cq.estimateFWHM(options.psfmeasure)
            log.info("FWHM= %s  STD= %s"%(efwhm, std))
        except Exception,e:
            log.error("There was some error with image %s : %s "%(options.input_image, str(e)))
    else:
        log.error("No input file given.")
        sys.exit(0)

