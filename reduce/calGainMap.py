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
#
# PAPI (PANIC PIpeline)
#
# calGainMap.py
#
# Compute a super sky flat using the dither frames (IRAF implementation)
#
# Created    : 23/09/2010    jmiguel@iaa.es -
# Last update: 23/09/2009    jmiguel@iaa.es - 
#              21/10/2009    jmiguel@iaa.es - Add normalization wrt chip 1 and 
#                                             support for MEF
#              17/09/2013    jmiguel@iaa.es - prevent zero-division
# TODO:
#  - include bpm 
################################################################################

################################################################################
# Import necessary modules

import getopt
import sys
import os
import tempfile
import logging
from optparse import OptionParser

import misc.fileUtils
import misc.utils as utils


# Interact with FITS files
import astropy.io.fits as fits
import numpy as np
import datahandler
import reduce.calSuperFlat


# Logging
from misc.paLog import log
from misc.version import __version__

class SkyGainMap(object):
    """ Compute the gain map from a list of sky frames """
    def __init__(self, filelist,  output_filename="/tmp/gainmap.fits",  
                 bpm=None, temp_dir="/tmp/"):

        self.framelist = filelist
        self.output = output_filename
        self.bpm = bpm
        self.temp_dir = temp_dir
        
    def create(self):
        """ Creation of the Gain map"""
        
        #First, we create the Super Sky Flat
        try:
            output_fd, tmp_output_path = tempfile.mkstemp(suffix='.fits')
            os.close(output_fd)
            superflat = reduce.calSuperFlat.SuperSkyFlat(self.framelist, 
                                                         tmp_output_path, 
                                                         self.bpm, norm=False,
                                                         temp_dir=self.temp_dir)
            superflat.create()
        except Exception,e:
            log.error("Error while creating super sky flat: %s", str(e))
            raise e
        
        # Secondly, we create the proper gainmap
        try:
            g = GainMap(tmp_output_path, self.output)
            g.create()   
        except Exception,e:
            log.error("Error while creating gain map: %s", str(e))
            raise e
        
        os.remove(tmp_output_path)
        return self.output      

class TwlightGainMap(object):
    """ Compute the gain map from a list of twlight frames """
    
    def __init__(self, flats_filelist, master_dark, 
                 output_filename="/tmp/gainmap.fits",  
                 bpm=None, temp_dir="/tmp/"):

        self.framelist = flats_filelist
        self.master_dark = master_dark
        self.output = output_filename
        self.bpm = bpm
        self.temp_dir = temp_dir
        
    def create(self):
        """ Creation of the Gain map"""
        
        #First, we create the Twlight Sky Flat
        try:
            output_fd, tmp_output_path = tempfile.mkstemp(suffix='.fits')
            os.close(output_fd)
                
            twflat = reduce.calTwFlat.MasterTwilightFlat (self.framelist,
                                                            self.master_dark, 
                                                            tmp_output_path)
            twflat.createMaster()
        except Exception,e:
            log.error("Error while creating twlight flat: %s", str(e))
            raise e
        
        # Secondly, we create the gainmap
        try:
            g = GainMap(tmp_output_path, self.output)
            g.create()   
        except Exception,e:
            log.error("Error while creating gain map: %s", str(e))
            raise e
        
        # Clean-up
        os.remove(tmp_output_path)
        
        return self.output      
               

               
class DomeGainMap(object):
    """ Compute the gain map from a list of dome (lamp-on,lamp-off) frames """
    
    def __init_(self, filelist,  output_filename="/tmp/domeFlat.fits", bpm=None):
        """ """
        self.framelist = filelist
        self.output = output_filename
        self.bpm = bpm
        
    def create(self):
        """ Creation of the Gain map"""
        
        #First, we create the Dome Flat
        try:
            output_fd, tmp_output_path = tempfile.mkstemp(suffix='.fits')
            os.close(output_fd)
            domeflat = reduce.calDomeFlat.MasterDomeFlat(self.framelist, 
                                                         temp_dir="/tmp", 
                                                         output_filename=tmp_output_path)
            domeflat.create()
        except Exception,e:
            log.error("Error while creating master dome flat: %s", str(e))
            raise e
        
        # Secondly, we create the proper gainmap
        try:
            g = GainMap(tmp_output_path, self.output)
            g.create()   
        except Exception,e:
            log.error("Error while creating gain map: %s", str(e))
            raise e 
        os.remove(tmp_output_path)
        return self.output 
                 
class GainMap(object):
    """
    Description
    -----------
    Build a Gain Map from a Flat Field image (dome, twilight, science-sky)  
            
    Author
    ------
    JMIbanez, IAA-CSIC
        
    """
    def __init__(self,  flatfield,  output_filename="/tmp/gainmap.fits",  
                 bpm=None, do_normalization=True, mingain=0.5, maxgain=1.5, 
                 nxblock=16, nyblock=16, nsigma=5):
        """
        Initialization method
        
        Parameters
        ----------
        
        flatfield: str
            Master flat field (dome, sky) from which the gainmap will be created.
        
        output_filename: str
            Output filename of the gainmap to build
        
        bpm: str (optional) -- NOT USED YET
            Bad pixel mask to use optionally for the gainmap build 
        
        do_normalization: bool
            If true, a normalization will be done; and if the input master
            flat is a multi-chip frame, the gainmap will be normalized wrt chip1
            However, be aware that the master flat might already have the 
            normalization done.
        
        mingain : int (0.5)
            Minimal gain; pixels below this gain value are considered bad and 
            set to 0 
        
        maxgain: int (1.5)
            Maximal gain; pixels above this gain value are considered bad and 
            set to 0. 
        
        nxblock: int (16)
            X-size (pixels) of box used to compute local bkg (even)
        
        nyblock: int (16)
            Y-size (pixels) of box used to compute local bkg (even)
        
        nsigma: int (10)
            Number of (+|-) stddev from local bkg to be bad pixel (default=10)
            
        """
         
        self.flat = flatfield  # Flat-field image (normalized or not, because optionaly, normalization can be done here)
        self.output_file_dir = os.path.dirname(output_filename)
        self.output_filename = output_filename  # full filename (path+filename)
        self.bpm = bpm
        self.do_norm = do_normalization
        
        # Some default parameter values
        self.m_MINGAIN = mingain #pixels with sensitivity < MINGAIN are assumed bad 
        self.m_MAXGAIN = maxgain #pixels with sensitivity > MAXGAIN are assumed bad 
        self.m_NXBLOCK = nxblock #image size should be multiple of block size 
        self.m_NYBLOCK = nyblock
        self.m_NSIG = nsigma  #badpix if sensitivity > NSIG sigma from local bkg
        self.m_BPM = bpm   #external BadPixelMap to take into account   (TODO) 
                
               
    def create(self):
        
        """
        Given a NOT normalized flat field, compute the gain map taking 
        into account the input parameters and an optional Bad Pixel Map (bpm)
        
        Return
        ------
        If success, return the output filename image of the Gain Map generated, 
        where Bad Pixels = 0.0

        """
        
        log.debug("Start creating Gain Map for file: %s", self.flat) 
        if os.path.exists(self.output_filename): os.remove(self.output_filename)

        # Check if we have a MEF file
        f = datahandler.ClFits ( self.flat )
        isMEF = f.mef
        if (not isMEF): nExt = 1
        else: nExt = f.next
        
        naxis1 = f.naxis1
        naxis2 = f.naxis2
        offset1 = int(naxis1*0.1)
        offset2 = int(naxis2*0.1)
        nbad = 0

        if f.getInstrument()=='panic' and naxis1==4096 and naxis2==4096:
            is_a_panic_full_frame = True # i.e., a single extension (GEIRS) full frame
        else: is_a_panic_full_frame = False

        # Check if normalization is already done to FF or otherwise it must be
        # done here
         
        
        gain = np.zeros([nExt, naxis1, naxis2], dtype=np.float32)
        myflat = fits.open(self.flat)
        
        # Check if normalization is already done to FF or otherwise it must be
        # done here
        if isMEF: extN = 1
        else: extN = 0
        if np.median(myflat[extN].data)> 100:
            self.do_norm = True
            log.info("Normalization will be done !")
        else:
            self.do_norm = False
            log.info("**No** normalization will be done !")
            
        for chip in range(0, nExt):
            log.debug("Operating in CHIP %d", chip+1)
            if isMEF:
                #flatM=np.reshape(myflat[chip+1].data, naxis1*naxis2)
                flatM = myflat[chip+1].data                 
            else:
                flatM = myflat[0].data
            
            # ##############################################################################
            # Normalize the flat (if MEF, all extension is normalized wrt extension/chip 1)#
            # ##############################################################################
            if chip==0:
                if self.do_norm and not is_a_panic_full_frame:
                    median = np.median(flatM[offset1:naxis1-offset1, offset2:naxis2-offset2])
                    mean = np.mean(flatM[offset1:naxis1-offset1, offset2:naxis2-offset2])
                    mode = 3*median-2*mean
                    log.debug("MEDIAN= %f  MEAN=%f MODE(estimated)=%f ", median, mean, mode)
                    log.debug("Normalizing flat-field by MEDIAN ( %f ) value", median)
                
                elif self.do_norm and is_a_panic_full_frame:
                    median = np.median(flatM[offset1/2:naxis1/2-offset1/2, offset2/2:naxis2/2-offset2/2])
                    mean = np.mean(flatM[offset1/2:naxis1/2-offset1/2, offset2/2:naxis2/2-offset2/2])
                    mode = 3*median-2*mean
                    log.debug("MEDIAN= %f  MEAN=%f MODE(estimated)=%f ", median, mean, mode)
                    log.debug("Normalizing (PANIC full-frame) flat-field by MEDIAN ( %f ) value", median)

                else: median = 1.0 # normalization not required (but must be already done !!)
                
            # To avoid zero-division
            __epsilon = 1.0e-20
            if np.fabs(median) > __epsilon:
                flatM = flatM/median
            else:
                flatM = flatM
                
            # Check for bad pixel 
            gain[chip] = np.where( (flatM < self.m_MINGAIN) | (flatM > self.m_MAXGAIN), 0.0, flatM)
            m_bpm = np.where(gain[chip] == 0.0, 1, 0) # bad pixel set to 1
            nbad = (m_bpm == 1).sum()
            log.debug("Initial number of Bad Pixels : %d ", nbad)
                        
            # local dev map to find out pixel deviating > NSIGMA from local median
            dev = np.zeros((naxis1, naxis2), dtype=np.float32)
            buf = np.zeros((self.m_NXBLOCK, self.m_NYBLOCK), dtype=np.float32)
            
            # Foreach image block
            for i in range(0, naxis1, self.m_NYBLOCK):
                for j in range(0, naxis2, self.m_NXBLOCK):
                    box = gain[chip][i: i + self.m_NXBLOCK, j: j + self.m_NYBLOCK]
                    p = np.where(box > 0.0)
                    buf = box[p]
                    if len(buf) > 0: med = np.median(buf)
                    else: med = 0.0
                    dev[i:i + self.m_NXBLOCK, j: j + self.m_NYBLOCK] = np.where(box > 0, (box - med), 0)
                            
            """                
            # Foreach image block
            for i in range(0,naxis1, self.m_NYBLOCK):
                for j in range(0,naxis2, self.m_NXBLOCK ):
                    n=0
                    for k in range(i, i+self.m_NYBLOCK):                   # foreach pix in block
                        for l in range(j, j+self.m_NXBLOCK):
                            if gain[chip][k*naxis1+l]>0.0:                 # if good pixel  
                                buf[n] = gain[chip][k*naxis1+l]
                                n+=1
                    #block median
                    if n>0: med = np.median(buf)
                    else: med = 0.0
                                
                    for k in range(i, i+self.m_NYBLOCK):                   # foreach pix in block
                        for l in range(j, j+self.m_NXBLOCK):
                            if gain[chip][k*naxis1+l]>0.0:                 # if good pixel, subtract median  
                                dev[k*naxis1+l] = gain[chip][k*naxis1+l] - med
                            else:
                                dev[k*naxis1+l] = 0.0                       # already known badpix
            """                    
            med = np.median(dev)
            sig = np.median(np.abs(dev - med)) / 0.6745
            lo  = med - self.m_NSIG * sig
            hi  = med + self.m_NSIG * sig
                                
            #log.debug("MED=%f LO=%f HI=%f SIGMA=%f", med, lo, hi, sig)                    
                                
            # Find more badpix by local dev
            p = np.where( (dev<lo) | (dev>hi))
            gain[chip][p] = 0.0 # badpix
            log.debug("Final number of Bad Pixel = %d", (gain[chip]==0.0).sum())
            
                 
        # Now, write result in a (MEF/single)-FITS file             
        output = self.output_filename
        log.debug('Generating output file: %s', output)
        prihdr = myflat[0].header.copy()
        prihdr.set('PAPITYPE','MASTER_GAINMAP', 'TYPE of PANIC Pipeline generated file')
        prihdr.set('PAPIVERS', __version__, 'PANIC Pipeline version')
        prihdr.add_history('Gain map based on %s' % self.flat)
        if self.do_norm:
            prihdr.add_history('Normalization wrt chip 0 done.')
        
        fo = fits.HDUList()
        # Add primary header to output file...
        if isMEF:
            prihdu = fits.PrimaryHDU(None,prihdr)
            fo.append(prihdu)
            # Add each extension
            for chip in range(0, nExt):
                hdu = fits.ImageHDU(data=gain[chip], header=myflat[chip+1].header)
                hdu.scale('float32') # important to set first data type ??
                #hdu.header.update('EXTVER',1)
                fo.append(hdu)
                del hdu
        else: 
            prihdu = fits.PrimaryHDU(gain[0],prihdr)
            fo.append(prihdu)
        
        fo.writeto(output,output_verify='ignore')
        fo.close(output_verify='ignore')
        del fo
                
        return output
                                    
#############################################################################
# main
# ###########################################################################
if __name__ == "__main__":
    # Get and check command-line options
    
    usage = "usage: %prog [options]"
    desc = """Creates a master gain map from a given master flat field (dome, twilight or superflat)
optionally normalized. The bad pixels are set to 0"""
    
    parser = OptionParser(usage, description=desc)
    
    parser.add_option("-s", "--source", type="str",
                  action="store", dest="source_file",
                  help="Flat Field image optionally normalized. It has to be a fullpath file name (required)")
                  
    parser.add_option("-o", "--output", type="str",
                  action="store", dest="output_filename", 
                  help="output file to write the Gain Map")
    
    #TODO
    #parser.add_option("-b", "--bpm", type="str",
    #              action="store", dest="bpm", 
    #              help="Input Bad pixel map file (optional)")
   
    parser.add_option("-L", "--low", type="float", default=0.5,
                  action="store", dest="mingain", 
                  help="pixel below this gain value  are considered bad (default=0.5)")
    
    parser.add_option("-H", "--high", type="float", default=1.5,
                  action="store", dest="maxgain", 
                  help="pixel above this gain value are considered bad (default=1.5)")
                  
    parser.add_option("-x", "--nx", type="int", default=16,
                  action="store", dest="nxblock", 
                  help="X dimen. (pixels) to compute local bkg (even) (default=16)")
    
    parser.add_option("-y", "--ny", type="int", default=16,
                  action="store", dest="nyblock", 
                  help="Y dimen. (pixels) to compute local bkg (even) (default=16)")
                  
    parser.add_option("-n", "--nsigma", type="int", default=10,
                  action="store", dest="nsigma", 
                  help="number of (+|-)stddev from local bkg to be bad pixel (default=10)")
    
    parser.add_option("-N", "--normal",  default=True,
                  action="store_true", dest="normal", 
                  help="if true, the input flat-field will be normalized before build the gainmap (default=True)")                            
                  
    (options, args) = parser.parse_args()
    
    if len(sys.argv[1:]) < 1:
       parser.print_help()
       sys.exit(0)
        
    # args is the leftover positional arguments after all options have been processed 
    if not options.source_file or not options.output_filename or len(args)!=0: 
        parser.print_help()
        parser.error("incorrect number of arguments " )
   
    try:
        gainmap = GainMap(options.source_file, options.output_filename, bpm=None,
                      do_normalization=options.normal,
                      mingain=options.mingain, maxgain=options.maxgain, 
                      nxblock=options.nxblock, nyblock=options.nyblock, 
                      nsigma=options.nsigma )
                      
        gainmap.create()
    except Exception,e:
        log.error("Some kind of problem happened %s"%str(e))
          
        
        
