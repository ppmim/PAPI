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

import astropy.io.fits as fits
import numpy

# Logging
from misc.paLog import log

import datahandler
import misc.fileUtils
import misc.utils as utils
import misc.robust as robust


class BadPixelMask(object):
    """
    Generates a Bad Pixel Mask (BPM) from a set of darks with fixed exp. time,
    and set of dome flat images. The output is a FITS file with the bad pixels 
    coded with 1 and good pixels as 0.

    - At least one set (>2) of darks or flats are required as input. 
    - Inputs files can be MEF files.

    """
    

    def __init__(self, dark_list, flat_list, outputfile=None, dthr=75.0, 
                fthr=15.0, temp_dir="/tmp"):
        
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
        self.dthr = float(dthr)
        self.fthr = float(fthr)
        self.temp_dir = temp_dir

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
            1.1 Define the threshold as :  75% of (median-3*sigma)
            1.2 Mask (==1) pixels **above** the threshold.
         2. Combine all of the dome flats into a master --> COLD pixels
            2.1 Define the threshold as :  15% of (median)
            2.2 Mask (==1) pixels **below** the threshold.
         3. Combine the HOT and COLD masks
        """
        
        t = utils.clock()
        t.tic()

        __epsilon = 1.0e-20
        

        if self.dark_list==None and self.flat_list==None:
            msg = "Neither Darks nor Flats images provided !"
            log.error(msg)
            raise Exception(msg)

        if self.dark_list!=None and len(self.dark_list)<3:
            log.error('Not enough darks provided. At least 3 darks frames are required')
            raise Exception("Not enough darks provided. At least 3 darks frames are required")
        
        if self.flat_list!=None and len(self.flat_list)<3:
            log.error('Not enough dome flats provided. At least 3 flat frames are required')
            raise Exception("Not enough dome flats provided. At least 3 flat frames are required")
        
        dark = None
        if self.dark_list!=None:    
            # STEP 1: Make the combine of dark frames
            log.debug("Combining DARKS frames...")
            dark_comb = self.temp_dir + '/darkcomb.fits'
            misc.fileUtils.removefiles(dark_comb)
            # Call IRAF task (it works with MEF or simple images)
            # With next combine, cosmic rays are rejected.
            iraf.mscred.darkcombine(input=("'"+"@"+self.dark_list+"'").replace('//','/'), 
                            output=dark_comb, 
                            combine='median', 
                            ccdtype='', 
                            process='no', 
                            reject='sigclip', 
                            scale='exposure'
                            )
            log.debug("Created combined Dark %s"%dark_comb)

            # STEP 1.1: Define the threshold as: 75% of (mean-3*sigma)
            dark = fits.open(dark_comb)
            nExt = 1 if len(dark)==1 else len(dark)-1
            if nExt==1: nx1,nx2 = dark[0].data.shape
            else: nx1,nx2 = dark[1].data.shape
            bpm = numpy.zeros([nExt, nx1, nx2], dtype=numpy.uint8)
            nbad_hot = numpy.zeros(nExt)

            # Loop over extensions
            for i_nExt in range(0, nExt):
                ext = i_nExt + int(nExt>1)
                if 'DET_ID' in dark[ext].header: 
                    det_id = dark[ext].header['DET_ID']
                else: 
                    det_id = i_nExt+1
                log.info("*** Detector %s"%det_id)
                median = numpy.mean(dark[ext].data)
                std = robust.std(dark[ext].data)
                dark_threshold = (median + self.dthr*std)
                #dark_threshold = self.dthr
                #dark_threshold = (mean + 3*std)*(self.dthr/100.0)
                log.debug("   Dark Median = %s"%median)
                log.debug("   Dark STD = %s"%std)
                log.debug("   Dark Threshold (median + D*sigma) = %s"%dark_threshold)
 
                # STEP 1.2: Mask (==1) pixels **above** the threshold.
                bpm[i_nExt, (dark[ext].data > dark_threshold) | numpy.isnan(dark[ext].data)] = 1
                nbad_hot[i_nExt] = (bpm[i_nExt]==1).sum()
                log.info("   # Hot-Bad pixels from Dark of detector %d : %d"
                    %(i_nExt+1, nbad_hot[i_nExt]))
        else:
            nbad_hot = 0
            log.info("# Hot-Bad pixels from Dark : No Darks provided !")


        # STEP 2: Make the combine of dome Flat frames
        # - Build the frame list for IRAF
        if self.flat_list!=None:    
            log.debug("Combining Flat frames...")
            flat_comb = self.temp_dir + '/flatcomb.fits'
            misc.fileUtils.removefiles(flat_comb)
            # Call IRAF task (it works with MEF or simple images)
            try:
                # With next combine, cosmic rays are rejected.
                iraf.mscred.flatcombine(input=("'"+"@"+self.flat_list+"'").replace('//','/'), 
                                output=flat_comb, 
                                combine='median', 
                                ccdtype='', 
                                process='no', 
                                reject='sigclip', 
                                scale='mode',
                                subsets='no'
                                )
                log.debug("Created combined Flat %s"%flat_comb)
            except Exception,e:
                raise e

            flat = fits.open(flat_comb)
            
            nExt = 1 if len(flat)==1 else len(flat)-1
            if nExt==1: nx1,nx2 = flat[0].data.shape
            else: nx1,nx2 = flat[1].data.shape
            nbad_cold = numpy.zeros(nExt)

            if self.dark_list==None:
                bpm = numpy.zeros([nExt, nx1, nx2], dtype=numpy.uint8)
                nbad_hot = numpy.zeros(nExt)

            # Loop over extensions
            for i_nExt in range(0, nExt):
                ext = i_nExt + int(nExt>1)
                if 'DET_ID' in flat[ext].header: 
                    det_id = flat[ext].header['DET_ID']
                else: 
                    det_id = i_nExt+1
                log.info("*** Detector %s"%det_id)
                ## Note: robust mean is very similar to median
                median = numpy.median(flat[ext].data)
                std = robust.std(flat[ext].data)
                # STEP 2.1: Define the thresholds as : +-%-ile of median, i.e., 
                # (tails of the distribution).
                low_flat_threshold = median*(self.fthr/100.0)  # COLD pixels
                # High_flat_threshold is not sure to have a good value...
                high_flat_threshold = median*(1+(100-self.fthr)/100.0) # HOT_Sat pixels
                log.info("    Flat Median = %s"%median)
                log.info("    Flat STD = %s"%std)
                log.info("    Flat Low  Threshold = %s"%low_flat_threshold)
                log.info("    Flat High Threshold = %s"%high_flat_threshold)


                # STEP 2.2: Mask (==1) pixels below (cold) and above (sat) the 
                # threshold. The resulted BPM is combined with bad pixel 
                # from dark frames.
                bpm[i_nExt, (flat[ext].data < low_flat_threshold) | 
                            (flat[ext].data > high_flat_threshold) |
                            numpy.isnan(flat[ext].data)] = 1

                nbad_cold[i_nExt] = (bpm[i_nExt]==1).sum() - nbad_hot[i_nExt]
                
                log.info("    # Cold-Bad pixels from Flats of detector %d: %d"
                        %(i_nExt+1, nbad_cold[i_nExt]))
        else:
            nbad_cold = 0
            log.info("# Cold-Bad pixels from Flats : No Flats provided !")

        # Show stats
        for i_nExt in range(0, nExt):
            nbad = (bpm[i_nExt]==1).sum()
            badfrac = float(nbad/float(bpm[i_nExt].size))
            ext = i_nExt + int(nExt>1)
            if dark!=None and 'DET_ID' in dark[ext].header: 
                det_id = dark[ext].header['DET_ID']
            else: 
                det_id = i_nExt+1
            log.info("# Bad pixels of detector %s: %d"%(det_id, nbad))
            log.info("Fraction Bad pixel of detector %s: %f"%(det_id, badfrac))

        # STEP 6: Save the BPM ---
        misc.fileUtils.removefiles(self.output)
        hdulist = fits.HDUList()     
        if self.flat_list!=None: hdr0 = flat[0].header
        else: hdr0 = dark[0].header

        prihdu = fits.PrimaryHDU (data = None, header = None)
        try:
            if 'INSTRUME' in hdr0: prihdu.header.set('INSTRUME', hdr0['INSTRUME'])
            if 'TELESCOP' in hdr0: prihdu.header.set('TELESCOP', hdr0['TELESCOP'])
            if 'CAMERA' in hdr0: prihdu.header.set('CAMERA', hdr0['CAMERA'])
            if 'MJD-OBS' in hdr0: prihdu.header.set('MJD-OBS', hdr0['MJD-OBS'])
            if 'DATE-OBS' in hdr0: prihdu.header.set('DATE-OBS', hdr0['DATE-OBS'])
            if 'DATE' in hdr0: prihdu.header.set('DATE', hdr0['DATE'])
            if 'UT' in hdr0: prihdu.header.set('UT', hdr0['UT'])
            if 'LST' in hdr0: prihdu.header.set('LST', hdr0['LST'])
            if 'ORIGIN' in hdr0: prihdu.header.set('ORIGIN', hdr0['ORIGIN'])
            if 'OBSERVER' in hdr0: prihdu.header.set('OBSERVER', hdr0['OBSERVER'])
        except Exception,e:
            log.warning("%s"%str(e))

        prihdu.header.set('PAPITYPE','MASTER_BPM','TYPE of PANIC Pipeline generated file')
        
        src_files = []
        if self.dark_list:
            src_files = [line.replace( "\n", "") for line in 
                                        fileinput.input(self.dark_list)]
        if self.flat_list:
            src_files+= [line.replace( "\n", "") for line in 
                                        fileinput.input(self.flat_list)]
        
        prihdu.header.add_history('BPM created from %s' % src_files)

        if nExt>1:
            prihdu.header.set('EXTEND', True, after = 'NAXIS')
            prihdu.header.set('NEXTEND', nExt)
            prihdu.header.set('FILENAME', self.output)
            hdulist.append(prihdu)
            for i_ext in range(0, nExt):
                hdu = fits.PrimaryHDU()
                hdu.scale('int16') # important to set first data type
                hdu.data = bpm[i_ext]
                ext = i_nExt + int(nExt>1)
                if dark!=None and 'DET_ID' in dark[ext].header:
                    hdu.header.set('DET_ID', dark[ext].header['DET_ID'])
                hdulist.append(hdu)
                del hdu
        else:
            prihdu.scale('int16') # important to set first data type
            prihdu.data = bpm
            hdulist.append(prihdu)
         
        
        # write bpm as FITS file (0==good pixels, 1==bad pixels)
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
Creates a bad pixel mask (BPM) from a set of darks with fixed exp. time, and
and set of dome flat images. The output is a FITS file with the bad pixels coded 
with 1 and good pixels as 0.
At least one set (>2) of darks or flats are required as input. Inputs files can
be MEF files or simple FITS, but prefered MEF files to distinguish the detectors.
"""

parser = OptionParser(usage, description=desc)
 
               
parser.add_option("-d", "--darks_source",
              action="store", dest="darks_file_list",
              help="list of input Darks images (at least 3)")

parser.add_option("-f", "--flats_source",
              action="store", dest="flats_file_list",
              help="list of input (optionally flats corrected) Dome Flat"
              " images (at least 3)")

parser.add_option("-o", "--output",
              action="store", dest="output_filename", type='str',   
              help="The output bad pixel mask (optional)")

parser.add_option("-D", "--dark_threshold",
              action="store", dest="dthr", type='float', default=20.0,
              help="The Dark sigma rejection threshold (above cut "
                "detector + 3*sigma counts) [default=%default]")

parser.add_option("-F", "--flat_threshold",
              action="store", dest="fthr", type='float', default=15.0,
              help="The Flat rejection threshold (below % cut of mean "
                "counts) [default=%default]")


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
    if (not options.darks_file_list and not options.flats_file_list):
        parser.print_help()
        parser.error("incorrect number of arguments " )
        return 2
        
    # Make sure we are not overwriting an existing file 
    if os.path.exists(options.output_filename):
        print "Error. The output file '%s' already exists."  % \
              (options.output_filename)
        return 1

    try:
        bpm = BadPixelMask(options.darks_file_list, options.flats_file_list,
                        options.output_filename, 
                        options.dthr, options.fthr)
        bpm.create()
    except Exception, e:
        log.error("Error running BPM: %s"%str(e))
        return 0
        
###############################################################################
if __name__ == "__main__":
    print 'Starting BadPixelMap....'
    sys.exit(main())
        
