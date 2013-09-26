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
# PANICtool
#
# calBPM.py
#
# Created    : 25/06/2009    jmiguel@iaa.es
# Last update: 25/06/2009    jmiguel@iaa.es
#              03/03/2010    jmiguel@iaa.es Added READMODE checking 
# 
# TODO:
#  NOT USED !!!! 
# - include master dark subtraction !!!
# - NO FUNCIONA BIEN, muy conservador ??! no da mascaras buenas !! 
#   da pocos pixeles malos ????
################################################################################


# Import necessary modules

import datetime
import os
import sys
import fileinput
from optparse import OptionParser


# Pyraf modules
import pyraf
from pyraf import iraf
from iraf import noao
from iraf import imred
from iraf import ccdred

import pyfits
import numpy

# Logging
from misc.paLog import log

import datahandler
import misc.fileUtils
import misc.utils as utils

class BadPixelMask(object):
    """
        Generate a bad pixel mask from a list of dark corrected dome flat images
        (extracted from VIRCAM pipeline, vircam_genbpm)
    
        A list of dark corrected dome flat images is given. A master flat 
        is created from all the input flats in the list. Each input flat is 
        then divided by the master. Bad pixels are marked on the new image as 
        those that are above or below the threshold (in sigma) in the new image. 
        Any pixel which has been marked as bad for more than a quarter of the 
        input images is defined as bad in the output mask. 
        
    """
    def __init__(self, input_file, outputfile, lthr=3.0, hthr=5.0, temp_dir="/tmp"):
        
        self.input_file = input_file # file with the list of files to read and process 
        # Default parameters values
        self.lthr = float(lthr)
        self.hthr = float(hthr)
        self.output = outputfile
        self.temp_dir = temp_dir
        
    
    def create(self):
        """ 
         Algorith to create the BPM
         -------------------------- 
         1. Combine all of the dome flats into a master
         2. Divide the resulting image by its median -->normalized MASTER_FLAT
         3. Create and zero the rejection mask
         4. Loop for all input images and divide each by the master flat
            
            4.1 Divide the resulting image by its median
            4.2 Get the standard deviation of the image
            4.3 Define the bad pixels
            
         5. Go through the rejection mask and if a pixel has been marked bad 
            more than a set number of times, then it is defined as bad
        """
        
        t = utils.clock()
        t.tic()
        
        
        # Read the file list
        filelist = [line.replace( "\n", "") for line in fileinput.input(self.input_file)]
        
        # Here we could check if each frame is a good dome flat !!!
        good_flats = []
        f_readmode = -1
        for iframe in filelist:
            fits=datahandler.ClFits(iframe)
            log.debug("Frame %s EXPTIME= %f TYPE= %s " %(iframe, fits.exptime, fits.type)) 
            # Check EXPTIME, TYPE (flat) and FILTER
            if not fits.isDomeFlat():
                log.warning("Warning: Task 'create BPM' found a non dome flat frame")
            else:
                # Check READMODE
                if ( f_readmode!=-1 and (f_readmode!= fits.getReadMode() )):
                    log.error("Error: Task 'calBPM' finished. Found a FLAT frame with different  READMODE")
                    raise Exception("Found a FLAT frame with different  READMODE") 
                else: 
                    f_readmode  =fits.getReadMode()
                    good_flats.append(iframe)
                
            
        if len(good_flats)<2:
            log.error('Not enough dome flats provided. At least 2 good flat frames are required')
            sys.exit(1)
        
        # Due a bug in PyRAF that does not allow a long list of files separated with commas as 'input' argument
        # we need to build again a text file with the good_files
        if len(good_flats)!=len(filelist):
            flats = self.temp_dir + "/flats.txt"
            ftemp = open(flats,"w")
            for flat in good_flats:
                ftemp.write(flat+"\n")
            ftemp.close() 
        else: flats = self.input_file
        
             
        # STEP 1: Make the combine of dome Flat frames
        # - Build the frame list for IRAF
        log.debug("Combining Flat frames...")
        flat_comb = self.temp_dir + '/flatcomb.fits'
        misc.fileUtils.removefiles(flat_comb)
        # Call IRAF task
        iraf.flatcombine(input='@'+flats.replace('//','/'), 
                        output=flat_comb, 
                        combine='median', 
                        ccdtype='none', 
                        process='no', 
                        reject='sigclip', 
                        subsets='no', 
                        scale='mode'
                        )
                        
        # STEP 2: Divide the resulting combined flat by its median (robust estimator)
        # Compute the mean of the image
        log.debug("Divide the resulting combined flat by its median...")
        median = float(iraf.imstat (
            images=flat_comb,
            fields='midpt',Stdout=1)[1])
        
        # Compute normalized flat
        iraf.imarith(operand1=flat_comb,
                    operand2=median,
                    op='/',
                    result=flat_comb,
                    )
                    
        f = pyfits.open(flat_comb)
        nx1 = f[0].header['NAXIS1']
        nx2 = f[0].header['NAXIS1']
        f.close()
                    
        # STEP 3: Create and zero the rejection mask
        bpm = numpy.zeros([nx1, nx2])
        
        # STEP 4: Loop for all input images and divide each by the master
        mflat = pyfits.open(flat_comb)
        tmpf = numpy.zeros([nx1, nx2])
        for flat in  good_flats:
            f_i = pyfits.open(flat)
            #ceros=(mflat[0].data==0).sum()
            #print "CEROS=", ceros
            mydata = numpy.where(mflat[0].data==0, 0.0001, mflat[0].data) # to avoid zero division error
            tmpf = f_i[0].data/mydata
            std = numpy.std(tmpf)
            tmpf.shape = nx1*nx2
            median = numpy.median(tmpf)
            mad = numpy.median(numpy.abs(tmpf - median))
            mad*=1.4826
            
            print ">>MEDIAN=",median
            print ">>STD=",std
            print ">>MAD=",mad
            
            #log.debug("Divide each flatted flat by its median")
            tmpf = tmpf/median
            
            low = 1.0 - self.lthr*mad/median
            high = 1.0 + self.hthr*mad/median
            
            print ">>LOW=", low
            print ">>HIGH=", high
            
            #STEP 4.3 Define the bad pixels
            tmpf.shape = nx1,nx2
            bpm[ (tmpf < low) | (tmpf > high)]+=1
            log.debug("BPM updated with current flat %s", flat)
            f_i.close()
            
        mflat.close()
        
        # STEP 5: Go through the rejection mask and if a pixel has been marked bad 
        # more than a set number of times (a quarter of number of images), 
        # then it is defined as bad.
        nbmax = numpy.maximum(2, len(good_flats)/4)
            
        bpm = numpy.where(bpm>nbmax,1,0) # bad pixel set to 1
        nbad = (bpm==1).sum()
        
        badfrac = float(nbad)/float(nx1*nx2)
        print ">>#Bad pixels", nbad
        print ">>BAD_FRAC", badfrac
        
        # Save the BPM
        misc.fileUtils.removefiles( self.output )               
        hdu = pyfits.PrimaryHDU()
        hdu.scale('int16') # importat to set first data type
        hdu.data = bpm     
        hdulist = pyfits.HDUList([hdu])
        hdu.header.update('PAPITYPE','MASTER_BPM','TYPE of PANIC Pipeline generated file')
        hdu.header.add_history('BPM created from %s' % good_flats)
        hdulist.writeto(self.output)
        hdulist.close(output_verify='ignore')
        
        # Remove temp files
        misc.fileUtils.removefiles(flat_comb)
        
        log.debug('Saved Bad Pixel Mask  to %s' , self.output)
        log.debug("createBPM' finished %s", t.tac() )
        

###############################################################################
usage = "usage: %prog [options] "
parser = OptionParser(usage)
 
               
parser.add_option("-s", "--source",
              action="store", dest="source_file_list",
              help="list of input (optionally dark corrected) dome flat images..")

parser.add_option("-o", "--output",
              action="store", dest="output_filename", 
              help="The output bad pixel mask.")

parser.add_option("-L", "--lthr",
              action="store", dest="lthr", type='float', default=3.0,
              help="The low rejection threshold in units of sigma [default 3]")

parser.add_option("-H", "--hthr",
              action="store", dest="hthr", type='float', default=5.0,
              help="The high rejection threshold in units of sigma [default 5]")

parser.add_option("-D", "--master_dark",
              action="store", dest="master_dark", type='str',
              help="[Optional] Master dark frame to subtract")    

parser.add_option("-S", "--show_stats",
              action="store_true", dest="show_stats", default=False,
              help="Show statistics [default False]")    

parser.add_option("-v", "--verbose",
              action="store_true", dest="verbose", default=True,
              help="verbose mode [default]")


################################################################################
# main
def main(arguments=None):
    
    if arguments is None:
        arguments = sys.argv[1:] # argv[0] is the script name
    (options, args) = parser.parse_args(args = arguments)

    if len(sys.argv[1:])<1:
       parser.print_help()
       return 2

    if len(args) !=0:
        parser.print_help()
        return 2 # used for command line syntax errors
    
    # Check mandatory arguments
    if not options.output_filename or not options.source_file_list:
        parser.print_help()
        parser.error("incorrect number of arguments " )
        return 2
        
    # Make sure we are not overwriting an existing file 
    if os.path.exists(options.output_filename):
        print "Error. The output file '%s' already exists."  % \
              (options.output_filename)
        return 1
    if options.master_dark:
        print "Sorry, dark subtraction not yet implemented."
        return 1
    
    try:
        bpm = BadPixelMask(options.source_file_list, options.output_filename, 
                       options.lthr, options.hthr)
        bpm.create()
    except Exception, e:
        log.error("Error creating BPM")
        return 0
        
###############################################################################
if __name__ == "__main__":
    print 'Start BadPixelMap....'
    sys.exit(main())
        
