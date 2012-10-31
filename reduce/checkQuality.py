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
################################################################################
# Import necessary modules

import datetime
import getopt
import numpy
import numpy.ma as ma 
import os
import pyfits
import sys
import subprocess

import numpy
import fileinput


import misc.utils as utils

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
    def __init__(self, input_file, isomin=10.0, ellipmax=0.3, edge=200, 
                 pixsize=0.45, write=False, verbose=False, bpm=None):
        
        self.input_file = input_file
        self.bpm = bpm
        # Default parameters values
        self.isomin = float(isomin)
        self.ellipmax = float(ellipmax)
        self.edge = int(edge)
        self.pixsize = float(pixsize)
        self.write = False
        self.verbose = False
        self.MIN_NUMBER_GOOD_STARS = 5
        
    
    def estimateFWHM(self):
        """ 
         FWHM estimation
         --------------- 
         Generating a ascii text catalog with Sextractor, we can read the FWHM values 
         and give an estimation of the FWHM computing the median of all the values
        
         It is very important that sextractor config file has the SATUR_LEVEL parameter
         with a suitable value. In other case, we won't get any value for FWHM. 
        
         SNR estimation as FLUX_AUTO/FLUXERR_AUTO or FLUX_APER/FLUXERR_APER
       
        Returns:
	-------
        A couple of values (fwhm,std) 
    
        """
        
        # Sextractor config

        
        sex_exe = os.environ['TERAPIX']+"/sex"
        sex_cnf = os.environ['PAPI_HOME']+"/irdr/src/config/default.sex"
        
        ## Sextractor Catalog columns
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
        
        catalog_file = "test.cat"
        
        sex_cmd = sex_exe + " " + self.input_file + " -c " + sex_cnf + \
        " -PIXEL_SCALE 0.45 -GAIN 4.15 -SATUR_LEVEL 1500000 " + "-CHECKIMAGE_TYPE NONE" + \
        " -CATALOG_TYPE ASCII -CATALOG_NAME  " + catalog_file
        
        # SExtractor execution
        
        try:
            if utils.runCmd( sex_cmd )==0:
                raise Exception("Some error happended while running SExtractor")
        except Exception,e:
            log.error("Some error while runCmd")
            raise e
        
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
            
        
        #fwhm_world=[float(line.split()[7]) for line in fileinput.input(source_file)]
        #matrix=[line.split() for line in fileinput.input(source_file)]
        #b=numpy.array(matrix)
        a = numpy.loadtxt(source_file)
        
        good_stars = []
        # Select 'best' stars for the estimation
        std = numpy.std(a[:,8])
        print "STD=",std
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
            snr = flux/flux_err
            if x>self.edge and x<naxis1-self.edge and y>self.edge and y<naxis2-self.edge \
               and ellipticity<self.ellipmax and fwhm>0.1 and fwhm<20 and flags==0  \
               and isoarea>float(self.isomin) and snr>20.0: # and fwhm<5*std it does not work many times
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
            print "STD2=",std
            efwhm = numpy.median(m_good_stars[:,8])
            print "best-FWHM-median(pixels)", efwhm
            print "FLUX_RADIUS (as mentioned in Terapix T0004 explanatory table) =", numpy.median(m_good_stars[:,9])
            print "Masked-mean=", ma.masked_outside(m_good_stars[:,8], 0.01, 3*std).mean()
            
            if self.write:
                fits_file[0].header.update('hierarch PAPI.SEEING', efwhm*self.pixsize)
                print "Fits keyword updated " 
            
        else:
            print "Not enough stars found !!"
            fits_file.close(output_verify='ignore')    
            return -1,-1
        
        #print "FWHM-median(pixels)= ", numpy.median(fwhm_world), numpy.amin(fwhm_world), numpy.amax(fwhm_world)
        #print "FWHM-mean(pixels)= ", numpy.mean(fwhm_world)
        print "\n----------"
        
        fits_file.close(output_verify='ignore')
        
        # cleanup files
        os.unlink(catalog_file)
        
        return efwhm,std
      
    def estimateBackground(self, output_file):
        """ 
         Background estimation
         ---------------------- 
         Run SExtractor to esimte the image background 
         
        """
        
        # Sextractor config

        sex_exe = os.environ['TERAPIX']+"/sex "
        sex_cnf = os.environ['PAPI_HOME']+"/irdr/src/config/default.sex"
        background_image = output_file
        
        sex_cmd = sex_exe + " " + self.input_file + " -c " + sex_cnf + " -PIXEL_SCALE 0.45 -GAIN 4.15 -SATUR_LEVEL 1500000 " +\
        " -CHECKIMAGE_TYPE  BACKGROUND -CHECKIMAGE_NAME  " + background_image
        
        # SExtractor execution
        #os.chdir("/disk-a/caha/panic/DEVELOP/PIPELINE/PAPI/")
        
        if utils.runCmd( sex_cmd )==0:
            raise Exception("Some error happended while running SExtractor")
        else:
            return  background_image     

            
#-----------------------------------------------------------------------

def usage():
    print ''
    print 'NAME'
    print '       checkQuality.py - FWHM estimation\n'
    print 'SYNOPSIS'
    print '       checkQuality.py [options] -f file.fits\n'
    print 'DESCRIPTION'
    print '       Give an estimation of the FWHM of the input image using'
    print '       values computed with SExtractor'
    print ' '
    print 'OPTIONS'
    print '       -v : verbose debugging output\n'
    print '       -i --isomin 10      minimun SExtractor ISOAREA_IMAGE (float)'
    print '       -e --ellipmax 0.2   maximun SExtractor ELLIPTICITY (float)'
    print '       -d --edge 200       consider sources out of image borders (int)'
    print '       -p --pixsize 0.45   pixel size (float)'
    print '       -w --write          update header with PA_SEEING keyword (bool)'
    print 'VERSION'
    print '       23 June 2009'
    print ''
    raise SystemExit

#-----------------------------------------------------------------------




################################################################################
# main
if __name__ == "__main__":
    print 'Start CheckQuality....'
    
    # Read command line parameters
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'f:i:e:d:p:wv',['file=','isomin=','ellipmax=','edge=','pixsize=','write','verbose'])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    
    nargs = len(sys.argv[1:])
    nopts = len(opts)
      
    
    isomin = 10
    ellipmax = 0.3
    edge = 200
    pixsize =0.45
    write = False
    verbose = False
    inputfile =''
            
    for option, par in opts:
        if option in ('-v','--verbose'):      # verbose debugging output
            verbose = True
            print "Verbose true"
        if option in ("-f", "--file"):
            inputfile = par
            print "inputfile=", inputfile
        if option in ("-i", "--isomin"):
            isomin = par
            print "isomin=", isomin
        if option in ("-e", "--ellipmax"):
            ellipmax = par
            print "ellipmax=",ellipmax
        if option in ("-d", "--edge"):
            edge = par
            print "edge=", edge
        if option in ("-p", "--pixsize"):
            pixsize = par
            print "pixsize=",pixsize
        if option in ("-w", "--write"):
            write = True
            print "write=", write
                
    
    # Error checking:
    if not os.path.exists(inputfile):      # check whether input file exists
        print inputfile, 'Input file does not exist'
        usage()
    
    print '...reading', inputfile
    try:
        cq = CheckQuality(inputfile, isomin, ellipmax, edge, pixsize, write, verbose)
        cq.estimateFWHM()
    except:
        log.error("There was some error !!")
        raise
    
    print 'ending application....'

    
