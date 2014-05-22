#!/usr/bin/env python

# Copyright (c) 2011-2014 IAA-CSIC  - All rights reserved. 
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
# Created    : 22/05/2014    jmiguel@iaa.es
# Last update: 22/05/2014    jmiguel@iaa.es
# 
# TODO:
#  -- Check dark and flat integrity (type, expT, readMode, size, etc.)
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
from iraf import mscred

import pyfits
import numpy

# Logging
from misc.paLog import log

import datahandler
import misc.fileUtils
import misc.utils as utils
import misc.robust as robust


class BadPixelMask(object):
    """
    Generate a bad pixel mask from a set of darks with fixed exp. time, and
    and set of dome flat images.

    """
    

    def __init__(self, dark_list, flat_list, outputfile=None, lthr=4.0, hthr=4.0, 
                 temp_dir="/tmp"):
        
        """
        Init method

        Parameters:
        -----------
        dark_file: str
            File with list of darks files to read and process for HOT pixels.
        flat_file: str
            File with list of flats files to read and process for COLD pixels. 
        outputfile: str
            Name of the output BPM generated.


        """
        self.dark_list = dark_list 
        self.flat_list = flat_list

        # Default parameters values
        self.lthr = float(lthr)
        self.hthr = float(hthr)
        self.temp_dir = temp_dir
        self.raw_flag = raw_flag

        if outputfile==None:
            dt = datetime.datetime.now()
            self.output = self.temp_dir + 'BadPixMask'+dt.strftime("-%Y%m%d%H%M%S")
        else:
            self.output = outputfile
        
    
    def create(self):
        """ 
         Algorith to create the BPM
         --------------------------      
         1. Combine all of the darks into a master --> HOT pixels
            1.1 Define the threshold as :  75% of (mean-3*sigma)
            1.2 Mask (==1) pixels **above** the the threshold.
         2. Combine all of the dome flats into a master --> COLD pixels
            2.1 Define the threshold as :  50% of (mean)
            2.2 Mask (==1) pixels **below** the the threshold.
         3. Combine the HOT and COLD masks
        """
        
        t = utils.clock()
        t.tic()

        __epsilon = 1.0e-20
        

        if len(self.dark_list)<2:
            log.error('Not enough darks provided. At least 2 darks frames are required')
            raise Exception("Not enough darks provided. At least 2 darks frames are required")
        if len(self.dome_list)<2:
            log.error('Not enough dome flats provided. At least 2 flat frames are required')
            raise Exception("Not enough dome flats provided. At least 2 flat frames are required")
        

        # STEP 1: Make the combine of dark frames
        log.debug("Combining DARKS frames...")
        dark_comb = self.temp_dir + '/darkcomb.fits'
        misc.fileUtils.removefiles(dark_comb)
        # Call IRAF task (it works with MEF or simple images)
        iraf.mscred.darkcombine(input=("'"+"@"+self.dark_list+"'").replace('//','/'), 
                        output=dark_comb, 
                        combine='median', 
                        ccdtype='', 
                        process='no', 
                        reject='sigclip', 
                        subsets='no', 
                        scale='exposure'
                        )
        log.debug("Created combined dark %s"%dark_comb)

        # STEP 1.1: Define the threshold as :  75% of (mean-3*sigma)
        dark = pyfits.open(dark_comb)
        mean = robust.mean(dark[0].data)
        std = robust.std(dark[0].data)
        dark_threshold = (mean + 3*std)*0.75
        log.debug("Dark MEAN (robust) of detector = %s"%mean)
        log.debug("Dark STD (robust) of detector = %s"%std)
        log.debug("Dark Threshold = %s"%dark_threshold)

        # STEP 1.2: Mask (==1) pixels **above** the the threshold.
        bpm = numpy.zeros([dark[0].data.shape[0], dark[0].data.shape[1]], 
                    dtype=numpy.uint8)
        bpm[(dark > dark_threshold) | numpy.isnan(dark)]=1



        # STEP 2: Make the combine of dome Flat frames
        # - Build the frame list for IRAF
        log.debug("Combining Flat frames...")
        flat_comb = self.temp_dir + '/flatcomb.fits'
        misc.fileUtils.removefiles(flat_comb)
        # Call IRAF task (it works with MEF or simple images)
        iraf.mscred.flatcombine(input=("'"+"@"+self.flat_list+"'").replace('//','/'), 
                        output=flat_comb, 
                        combine='median', 
                        ccdtype='', 
                        process='no', 
                        reject='sigclip', 
                        subsets='no', 
                        scale='mode'
                        )
        log.debug("Created combined flat %s"%flat_comb)

        # STEP 2.1: Define the threshold as :  75% of (mean-3*sigma)
        flat = pyfits.open(flat_comb)
        mean = robust.mean(flat[0].data)
        std = robust.std(flat[0].data)
        flat_threshold = (mean + 3*std)*0.75
        log.debug("Flat MEAN (robust) of detector = %s"%mean)
        log.debug("Flat STD (robust) of detector = %s"%std)
        log.debug("Flat Threshold = %s"%flat_threshold)

        # STEP 2.2: Mask (==1) pixels **below** the the threshold.
        bpm[(flat < flat_threshold) | numpy.isnan(flat)]=1


        # Show stats
        nbad = (bpm==1).sum()
        badfrac = float(nbad/float(dark[0].data.size)
        log.info("# Bad pixels : %f"%(nbad))
        log.info("Fraction Bad pixel : %f"%(badfrac))

        # STEP 6: Save the BPM ---
        misc.fileUtils.removefiles(self.output)
        hdulist = pyfits.HDUList()     
        hdr0 = pyfits.getheader(flat_list[0])
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

        prihdu.header.set('PAPITYPE','MASTER_BPM','TYPE of PANIC Pipeline generated file')
        prihdu.header.add_history('BPM created from %s' % good_flats)

        if nExt>1:
            prihdu.header.set('EXTEND', pyfits.TRUE, after = 'NAXIS')
            prihdu.header.set('NEXTEND', nExt)
            prihdu.header.set('FILENAME', self.output)
            hdulist.append(prihdu)
            for i_ext in range(0, nExt):
                hdu = pyfits.PrimaryHDU()
                hdu.scale('int16') # important to set first data type
                hdu.data = bpm[i_ext]
                hdulist.append(hdu)
                del hdu
        else:
            prihdu.scale('int16') # important to set first data type
            prihdu.data = bpm[0]
            hdulist.append(prihdu)
         
        
        # write FITS
        try:
            hdulist.writeto(self.output)
            hdulist.close(output_verify='ignore')
        except Exception,e:
            log.error("Error writing linearity model %s"%self.output)
            raise e

        # Remove temp files
        #misc.fileUtils.removefiles(flat_comb)
        
        log.debug('Saved Bad Pixel Mask  to %s' , self.output)
        log.debug("createBPM' finished %s", t.tac() )
        
    

###############################################################################
usage = "usage: %prog [options] "
desc = """
Generate a bad pixel mask from a set of darks with fixed exp. time, and
and set of dome flat images.
"""

parser = OptionParser(usage, description=desc)
 
               
parser.add_option("-d", "--darks_source",
              action="store", dest="darks_file_list",
              help="list of input (optionally dark corrected) darks images..")
parser.add_option("-f", "--flats_source",
              action="store", dest="flats_file_list",
              help="list of input (optionally flats corrected) dome flat images..")

parser.add_option("-o", "--output",
              action="store", dest="output_filename", 
              help="The output bad pixel mask.")

parser.add_option("-L", "--lthr",
              action="store", dest="lthr", type='float', default=4.0,
              help="The low rejection threshold in units of sigma [default=%default]")

parser.add_option("-H", "--hthr",
              action="store", dest="hthr", type='float', default=4.0,
              help="The high rejection threshold in units of sigma [default=%default]")


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
    #if options.master_dark:
    #    print "Sorry, dark subtraction not yet implemented."
    #    return 1
    
    try:
        bpm = BadPixelMask(options.dark_list, options.flat_list,
                        options.output_filename, 
                        options.lthr, options.hthr)
        bpm.create()
    except Exception, e:
        log.error("Error running BPM: %s"%str(e))
        return 0
        
###############################################################################
if __name__ == "__main__":
    print 'Starting BadPixelMap....'
    sys.exit(main())
        
