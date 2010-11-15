#!/usr/bin/env python
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
#              21/10/2019    jmiguel@iaa.es - Add normalization wrt chip 1 and support for MEF
# TODO
#  
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
import pyfits
import numpy as np
import datahandler
import reduce.calSuperFlat


# Logging
from misc.paLog import log

class SkyGainMap:
    """ Compute the gain map from a list of sky frames """
    def __init__(self, filelist,  output_filename="/tmp/superFlat.fits",  bpm=None):
        """ """
        self.framelist = filelist
        self.output = output_filename
        self.bpm = bpm
        
    def create(self):
        """ Creation of the Gain map"""
        
        #First, we create the Super Sky Flat
        try:
            output_fd, tmp_output_path = tempfile.mkstemp(suffix='.fits')
            os.close(output_fd)
            superflat = reduce.calSuperFlat.SuperSkyFlat(self.framelist, tmp_output_path, self.bpm, norm=False)
            superflat.create()
        except Exception,e:
            log.error("Error while creating super sky flat: %s", str(e))
            raise e
        
        # Secondly, we create the proper gainmap
        try:
            g=GainMap(tmp_output_path, self.output)
            g.create()   
        except Exception,e:
            log.error("Error while creating gain map: %s", str(e))
            raise e
        
        os.remove(tmp_output_path)
        return self.output      
               
class DomeGainMap:
    """ Compute the gain map from a list of dome (lamp-on,lamp-off) frames """
    
    def __init_(self, filelist,  output_filename="/tmp/domeFlat.fits",  bpm=None):
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
            domeflat = reduce.calDomeFlat.MasterDomeFlat(self.framelist, output_dir="/tmp", output_filename=tmp_output_path)
            domeflat.create()
        except Exception,e:
            log.error("Error while creating master dome flat: %s", str(e))
            raise e
        
        # Secondly, we create the proper gainmap
        try:
            g=GainMap(tmp_output_path, self.output)
            g.create()   
        except Exception,e:
            log.error("Error while creating gain map: %s", str(e))
            raise e 
        os.remove(tmp_output_path)
        return self.output 
                 
class GainMap:
    """
    \brief Class used to build a Gain Map from a Flat Field image (dome, twilight, science-sky)  
    
    \par Class:
        GainMap
    \par Purpose:
         Create a Gain Map from a master flat-field (dome, tw, sky)
    \par Description:
            
    \par Language:
        PyRaf
    \param file_list
        A list FITS files or directory
    \param output_filename
        File where log will be written
    \retval 0
        If no error
    \author
        JMIbanez, IAA-CSIC
        
    """
    def __init__(self,  flatfield,  output_filename="/tmp/gainmap.fits",  bpm=None):
         
        self.flat = flatfield  # Flat-field image (NOT normalized !!!!, because normalization is done here)
        self.output_file_dir = os.path.dirname(output_filename)
        self.output_filename = output_filename  # full filename (path+filename)
        self.bpm = bpm
        
        # Some default parameter values
        self.m_MINGAIN = 0.5   #pixels with sensitivity < MINGAIN are assumed bad 
        self.m_MAXGAIN = 1.5   #pixels with sensitivity > MAXGAIN are assumed bad 
        self.m_NXBLOCK = 16    #image size should be multiple of block size 
        self.m_NYBLOCK = 16
        self.m_NSIG    = 5.0   #badpix if sensitivity > NSIG sigma from local bkg
        self.m_BPM     = bpm   #external BadPixelMap to take into account   
                
                
    def create(self):
        # flat, nsig, nxblock, nyblock, mingain, maxgain, bpm=None ):
        
        """
        \brief Given a NOT normalized flat field, compute the gain map taking into account the input 
               parameters and an optional Bad Pixel Map (bpm)
        """
        
        log.debug("Start creating Gain Map") 
        if os.path.exists(self.output_filename): os.remove(self.output_filename)

        # Check if we have a MEF file
        f=datahandler.ClFits ( self.flat )
        isMEF=f.mef
        if (not isMEF): nExt=1
        else: nExt=f.next
        
        naxis1=f.naxis1
        naxis2=f.naxis2
        nbad=0
        
        gain=np.zeros([nExt, naxis1, naxis2], dtype=np.float32)
        myflat=pyfits.open(self.flat)
        for chip in range(0,nExt):
            log.debug("Operating in CHIP %d", chip+1)
            if isMEF:
                #flatM=np.reshape(myflat[chip+1].data, naxis1*naxis2)
                flatM=myflat[chip+1].data                            
            else:
                flatM=myflat[0].data
            
            # ##############################################################################
            # Normalize the flat (if MEF, all extension is normlized wrt extension/chip 1) #
            # ##############################################################################
            if chip==0:
                median=np.median(flatM[200:naxis1-200, 200:naxis2-200])
                mean=np.mean(flatM[200:naxis1-200, 200:naxis2-200])
                mode=3*median-2*mean
                log.debug("MEDIAN= %f  MEAN=%f MODE(estimated)=%f ", median, mean, mode)
                log.debug("Normalizing flat-field by MEDIAN ( %f ) value", median)
            flatM=flatM/median   
            
            # Check for bad pixel 
            gain[chip]=np.where( (flatM<self.m_MINGAIN) | (flatM>self.m_MAXGAIN), 0.0, flatM)
            m_bpm=np.where(gain[chip]==0.0, 1, 0) # bad pixel set to 1
            nbad=(m_bpm==1).sum()
            log.debug("Initial number of Bad Pixels : %d ", nbad)
                        
            # local dev map to find out pixel deviating > NSIGMA from local median
            dev = np.zeros((naxis1,naxis2), dtype=np.float32)
            buf = np.zeros((self.m_NXBLOCK,self.m_NYBLOCK), dtype=np.float32)
            
            # Foreach image block
            for i in range(0, naxis1, self.m_NYBLOCK):
                for j in range(0, naxis2, self.m_NXBLOCK):
                    box=gain[chip][i:i+self.m_NXBLOCK, j:j+self.m_NYBLOCK]
                    p=np.where(box>0.0)
                    buf=box[p]
                    if len(buf)>0: med = np.median(buf)
                    else: med = 0.0
                    dev[i:i+self.m_NXBLOCK, j:j+self.m_NYBLOCK] = np.where(box>0, (box - med), 0)
                            
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
            sig = np.median(np.abs(dev-med)) / 0.6745
            lo  = med - self.m_NSIG * sig
            hi  = med + self.m_NSIG * sig
                                
            #log.debug("MED=%f LO=%f HI=%f SIGMA=%f", med, lo, hi, sig)                    
                                
            # Find more badpix by local dev
            p=np.where( (dev<lo) | (dev>hi))
            gain[chip][p]=0.0
            log.debug("Final number of Bad Pixel = %d", (gain[chip]==0.0).sum())
            
                 
        # Now, write result in a (MEF/single)-FITS file             
        output=self.output_filename
        log.debug('Generating output file: %s', output)
        prihdr = myflat[0].header.copy()
        prihdr.update('PAPITYPE','MASTER_GAINMAP','TYPE of PANIC Pipeline generated file')
        prihdr.add_history('Gain map based on %s' % self.flat)
        fo = pyfits.HDUList()
        # Add primary header to output file...
        if isMEF: 
            prihdu=pyfits.PrimaryHDU(None,prihdr)
            fo.append(prihdu)
            # Add each extension
            for chip in range(0, nExt):
                hdu = pyfits.ImageHDU(data=gain[chip], header=myflat[chip+1].header)
                hdu.scale('float32') # importat to set first data type ??
                #hdu.header.update('EXTVER',1)
                fo.append(hdu)
                del hdu
        else: 
            prihdu=pyfits.PrimaryHDU(gain[0],prihdr)
            fo.append(prihdu)
        
        fo.writeto(output,output_verify='ignore')
        fo.close(output_verify='ignore')
        del fo        
                                    
################################################################################
# main
if __name__ == "__main__":
    # Get and check command-line options
    
    usage = "usage: %prog [options] arg1 arg2"
    parser = OptionParser(usage)
    
    parser.add_option("-v", "--verbose",
                  action="store_true", dest="verbose", default=True,
                  help="verbose mode [default]")
                  
    parser.add_option("-s", "--source",
                  action="store", dest="source_file",
                  help="Flat Field image NOT normalized. It has to be a fullpath file name (required)")
                  
    parser.add_option("-o", "--output",
                  action="store", dest="output_filename", help="output file to write the Gain Map")
    
    parser.add_option("-b", "--bpm",
                  action="store", dest="bpm", help="Bad pixel map (optional)")
   

    (options, args) = parser.parse_args()
    if not options.source_file or not options.output_filename or len(args)!=0: # args is the leftover positional arguments after all options have been processed
        parser.print_help()
        parser.error("incorrect number of arguments " )
    
    gainmap = GainMap(options.source_file, options.output_filename, options.bpm)
    gainmap.create()
          
        
        
