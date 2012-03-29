#!/usr/bin/env python
################################################################################
#
#
# PAPI (PANIC PIpeline)
#
# calSuperFlat.py
#
# Compute a super sky flat using the dither frames (IRAF implementation)
#
# Created    : 13/03/2009    jmiguel@iaa.es -
# Last update: 15/04/2009    jmiguel@iaa.es - Created function and modified to accept command line arguments
#              03/03/2010    jmiguel@iaa.es - Big modification to convert to a class and make more checkings
#              16/09/2010    jmiguel@iaa.es - Renamed to calSuperFlat and added support to MEF files
#              23/09/2010    jmiguel@iaa.es - Added (optional) gain map creation and/or normaliced flat field
#              21/10/2010    jmiguel@iaa.es - Removed gain map computation
#
# TODO
#    - dark subtraction ??? (SIMPLE doesn't it)
################################################################################

################################################################################
# Import necessary modules

import getopt
import sys
import os
import logging
import fileinput
import shutil
from optparse import OptionParser

import misc.fileUtils
import misc.utils as utils
import calGainMap 

# Interact with FITS files
import pyfits
import numpy
import datahandler


# Logging
from misc.paLog import log

# Pyraf modules
from pyraf import iraf
from iraf import noao
#from iraf import imred
#from iraf import ccdred
from iraf import mscred

class SuperSkyFlat(object):
    """
    \brief Class used to build a super sky Flat from a dither set of science 
           frames containing objects.
    
    \par Class:
        SuperSkyFlat
    \par Purpose:
         Create a super flat field 
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
        JMIbannez, IAA-CSIC
        
    """
    def __init__(self,  filelist,  output_filename="/tmp/superFlat.fits",  
                 bpm=None, norm=True, temp_dir="/tmp/", median_smooth=False):
        """
        Initialization method.
        
        Parameters:
        ----------
        
        filelist : list 
            It can be a file or a python-list having the list of files to use for the super-flat
            
        output_filename (optional): string
            Filename for the super flat created
            
        bpm (optional) : string 
            Bad pixel map to be used
            
        norm (optional): bool
            Flag to indicate if the super flat must be normalized
            
        median_smooth: bool
            If true, median smooth filter is applied to the combined Flat-Field

        """
                    
            
        if type(filelist)==type(list()): 
            self.filelist = filelist  # list of sources files to be used in sky-flat computation
        elif os.path.isfile(filelist):
            self.filelist= [line.replace( "\n", "") for line in fileinput.input(filelist)]
        else:
            raise Exception("Cannot read source files")
           
        self.temp_dir = os.path.dirname(output_filename)
        self.output_filename = output_filename  # full filename (path+filename)
        self.bpm = bpm
        self.norm = norm # if true, the flat field will be normalized
        self.__median_smooth = median_smooth
        
        # Some default parameter values
        self.m_MIN_N_GOOD = 2
        self.m_min_flats = 5
        self.m_MINGAIN = 0.5    #pixels with sensitivity < MINGAIN are assumed bad 
        self.m_MAXGAIN = 1.5    #pixels with sensitivity > MAXGAIN are assumed bad 
        self.m_NXBLOCK = 16     #image size should be multiple of block size 
        self.m_NYBLOCK = 16
        self.m_NSIG = 5.0       #badpix if sensitivity > NSIG sigma from local bkg

            
    def create(self):
      
        """
        \brief Create the super sky flat using sigma-clipping algorithm (and supporting MEF)
        """
        # del old files   
        log.debug("Start createSuperSkyFlat") 
        if os.path.exists(self.output_filename): os.remove(self.output_filename)
        
        # Check data integrity (all have the same properties)
        m_filelist=self.filelist
        #print "LISTA =", m_filelist
            
        if not datahandler.checkDataProperties( m_filelist ):
            log.error("Data integrity ERROR, some files not having same properties (FILTER, EXPTIME, NCOADDS or READMODE)")
            raise Exception("Found a data integrity error")
        
  
        tmp1 = (self.temp_dir + "/tmp_sf.fits").replace('//','/')
        misc.fileUtils.removefiles(tmp1)
        log.info("Combining images...(images are scaled to have the same median)")
        misc.utils.listToFile(m_filelist, self.temp_dir + "/files.txt") 
        # Combine the images to find out the super Flat using sigma-clip algorithm;
        # the input images are scaled to have the same median, the pixels containing 
        # objects are rejected by an algorithm based on the measured noise (sigclip),
        # and the flat-field is obtained by a median.
        iraf.mscred.combine(input=("'"+"@"+self.temp_dir+"/files.txt"+"'").replace('//','/'),
                    output=tmp1,
                    combine='median',
                    offset='none',
                    reject='sigclip',
                    lsigma=1.5,
                    hsigma=1.5,
                    scale='median',
                    zero='none'
                    #masktype='none'
                    #scale='exposure',
                    #expname='EXPTIME'
                    #ParList = _getparlistname ('flatcombine')
                )
        
        
        #Median smooth the superFlat
        ## Median smooth the master (normalized) flat
        if self.__median_smooth:
            log.debug("Doing Median smooth of FF ...")
            iraf.mscmedian(
                    input=tmp1,
                    output=tmp1.replace(".fits","_smooth.fits"),
                    xwindow=20,
                    ywindow=20,
                    outtype="median"
                    )
            shutil.move(tmp1.replace(".fits","_smooth.fits"), tmp1)
            
            #Or using scipy ( a bit slower then iraf...)
            #from scipy import ndimage
            #filtered = ndimage.gaussian_filter(f[0].data, 20)


        
        # (optional) Normalize the wrt chip 1
        if self.norm:
            f = pyfits.open(tmp1, ignore_missing_end=True)
            #MEF frame
            if len(f)>1:
                chip = 1 # normalize wrt to mode of chip 1
                naxis1 = f[0].header['NAXIS1']
                naxis2 = f[0].header['NAXIS2']
                offset1 = int(naxis1*0.1)
                offset2 = int(naxis2*0.1)
                mode = (3*numpy.median(f[chip].data[offset1:naxis1-offset1,
                                                    offset2:naxis2-offset2])-
                        2*numpy.mean(f[chip].data[offset1:naxis1-offset1, 
                                                  offset2:naxis2-offset2]))
                msg = "Normalization of MEF master flat frame wrt chip 1. (MODE=%d)"%mode
                
            # PANIC multi-chip full frame
            elif ('INSTRUME' in f[0].header and f[0].header['INSTRUME']=='panic'
                  and f[0].header['NAXIS1']==4096 and f[0].header['NAXIS2']==4096):
                # It supposed to have a full frame of PANIC in one single 
                # extension (GEIRS default)
                mode = (3*numpy.median(f[0].data[200:2048-200,200:2048-200])- 
                       2*numpy.mean(f[0].data[200:2048-200,200:2048-200]) )
                msg = "Normalization of (full) PANIC master flat frame wrt chip 1. (MODE=%d)"%mode
                
            else:
                naxis1 = f[0].header['NAXIS1']
                naxis2 = f[0].header['NAXIS2']
                offset1 = int(naxis1*0.1)
                offset2 = int(naxis2*0.1)
                mode = (3*numpy.median(f[0].data[offset1:naxis1-offset1,
                                                    offset2:naxis2-offset2])-
                        2*numpy.mean(f[0].data[offset1:naxis1-offset1, 
                                                  offset2:naxis2-offset2]))
                msg = "Normalization of master (O2k?) flat frame. (MODE=%d)"%mode 

            log.debug(msg)

            f.close()
            
            
            # Cleanup: Remove temporary files
            #misc.fileUtils.removefiles(self.output_filename)
            # Compute normalized flat
            iraf.mscred.mscarith(operand1=tmp1,
                    operand2=mode,
                    op='/',
                    pixtype='real',
                    result=self.output_filename.replace("//","/")
                    )
        else:
            os.rename(tmp1, self.output_filename) 
        
        
        # Update FITS header 
        f = pyfits.open(self.output_filename,'update')
        if self.norm: 
            f[0].header.add_history("[calSuperFlat] Normalized Super-Flat created from : %s"%str(m_filelist))
            f[0].header.add_history(msg)
        else:
            f[0].header.add_history("[calSuperFlat] Non-Normalized Super-Flat created from : %s"%str(m_filelist))
        f[0].header.update('PAPITYPE','MASTER_SKY_FLAT','TYPE of PANIC Pipeline generated file')
        
        #
        if 'PAT_NEXP' in f[0].header:
		          f[0].header.update('PAT_NEXP', 1, 'Number of Positions into the dither pattern')
	
        f[0].header.update('IMAGETYP','MASTER_SKY_FLAT','TYPE of PANIC Pipeline generated file')
        f.close(output_verify='ignore')
        log.debug("Image created : %s", self.output_filename)

        return self.output_filename
                                    
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
                  action="store", dest="source_file_list",
                  help="Source file list of data frames. It has to be a fullpath file name")
                  
    parser.add_option("-o", "--output",
                  action="store", dest="output_filename", help="output file to write SuperFlat")
    
    parser.add_option("-b", "--bpm",
                  action="store", dest="bpm", help="bad pixel map file (optional)", default=None)

    
    parser.add_option("-N", "--norm",
                  action="store_true", dest="norm", help="normalize output \
                  SuperFlat. If image is multi-chip, normalization wrt chip 1 \
                  is done. (optional)", default=False)
    
    parser.add_option("-m", "--median_smooth",
                  action="store_true", dest="median_smooth", default=False,
                  help="Median smooth the combined flat-field [default False]")    

    (options, args) = parser.parse_args()
    if not options.source_file_list or not options.output_filename or len(args)!=0: 
        # args is the leftover positional arguments after all options have been 
        # processed
        parser.print_help()
        parser.error("incorrect number of arguments " )
    if options.verbose:
        print "reading %s ..." % options.source_file_list
    
    filelist=[line.replace( "\n", "") for line in fileinput.input(options.source_file_list)]
    try:
	superflat = SuperSkyFlat(filelist, options.output_filename, 
                             options.bpm, options.norm, "/tmp/",
                             options.median_smooth)
    	superflat.create()
    except Exception,e:
	log.error("Error: %s"%str(e))
          
