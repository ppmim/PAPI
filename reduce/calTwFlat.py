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
# calTwilightFlat.py
#
# Created    : 19/05/2009    jmiguel@iaa.es
# Last update: 29/09/2009    jmiguel@iaa.es
#              14/12/2009    jmiguel@iaa.es  - Check NCOADDS; use ClFits class; 
#                                              Skip non TW flats and continue 
#                                              working with the good ones
#              03/03/2010    jmiguel@iaa.es  - Added READMODE checking 
#              20/09/2010    jmiguel@iaa.es  - Added support of MEF files
#              18/11/2010    jmiguel@iaa.es  - Added optional normalization by 
#                                              mode
#              22/02/2012    jmiguel@iaa.es  - Added Dark model use to build the 
#                                              required scaled dark 
#              27/03/2012    jmiguel@iaa.es  - Fixed bug wrt chip 1 normalization
#              03/12/2012    jmiguel@iaa.es  - Modified normalization by median instead of mode
#              05/08/2014    jmiguel@iaa.es  - Added support for DOME_FLAT series (not lamp on/off)
#
# TODO:
#   - take into account BPM !!!
#   - compute automatically the level of counts the twilight flats should 
#     have (lthr,hthr)
#   - median smooth
################################################################################

################################################################################
# Import necessary modules

import sys
import os
import logging
import fileinput
import time
import shutil
from optparse import OptionParser

import misc.fileUtils
import misc.utils as utils
from misc.version import __version__
import datahandler


# Pyraf modules
from pyraf import iraf
from iraf import noao
from iraf import mscred

import numpy

# Interact with FITS files
import astropy.io.fits as fits

# Logging
from misc.paLog import log

class ExError(Exception):
    pass

class MasterTwilightFlat (object):
    """
    Description:
    ------------
    Class used to build and manage a master calibration twilight flat.
    
    Twilight flats are quite good for low-spatial frequency QE variation 
    across the chip (large scale variation), but not for high-spatial 
    frequency (small scale  variations).

        
        1. Check the  TYPE(twilight) and FILTER of each Flat frame
        If any frame on list missmatch the FILTER, then the master 
        twflat will skip this frame and contiune with then next ones.
        EXPTIME do not need be the same, so EXPTIME scaling with 'mode' will be 
        done
           
           1.1: Check either over or under exposed frames
        
        2. We subtract a proper MASTER_DARK, it is required for TWILIGHT FLATS 
        because they might have diff EXPTIMEs
        
        3. Make the combine (with sigclip rejection) of dark subtracted Flat 
        frames scaling by 'mode'
        
        4. Normalize the tw-flat dividing by the mean value
        
    Author:
        JMIbannez, IAA-CSIC
  
    """
    def __init__(self, flat_files, master_dark_model,
                 output_filename="/tmp/mtwflat.fits", lthr=1000, hthr=100000,
                 bpm=None, normal=True, temp_dir="/tmp/", median_smooth=False):
        
        """
        Initialization method.
        
        Parameters
        ----------
        
        flat_files: list
        
        master_dark_model : string
            Master dark model to subtract (required) 
        output_filename: string
            Filename of the output master tw-flat to build
        lthr: int
            Low threshold to identify good twilight flats (default 100)
        hthr: int
            High threshold to identify good twilight flats (default 100000)
        bpm: string
            Bad Pixel mask to use (optional)
        normal: bool
            If true, the normalization will be done.
        temp_dir: string
            Temporal directory for temporal files needed.
        median_smooth: bool
            If true, median smooth filter is applied to the combined Flat-Field
            
        
        """
        
        self.__input_files = flat_files
        self.__master_dark = master_dark_model # if fact, it is a dark model
        self.__output_filename = output_filename  # full filename (path+filename)
        self.__bpm = bpm
        self.__normal = normal
        self.__temp_dir = temp_dir #temporal dir used for temporal/intermediate files
        self.__median_smooth = median_smooth
        
        self.m_MIN_N_GOOD = 3
        self.m_lthr = lthr
        self.m_hthr = hthr
        self.m_min_flats = 5
        
    
    def createMaster(self):
      
        """
        Create a master Tw FLAT from the flat file list
        """   
        log.debug("Start createMasterTwilightFlat")
        
        start_time = time.time()
        t=utils.clock()
        t.tic()
    
        # Cleanup old files
        
        misc.fileUtils.removefiles(self.__output_filename)
        
        # Get the user-defined list of flat frames
        framelist = self.__input_files
        
        # Check exist master DARK
        if not os.path.exists( self.__master_dark  ):
            log.error('Cannot find frame : "%s"' % self.__master_dark)
            raise Exception("Master Dark not found")
        
        # Determine the number of Flats frames to combine
        try:
            nframes = len(framelist[0])
        except IndExError:
            log.error("No FLAT frames defined")
            raise ExError('No FLAT frames defined')
        
        if nframes<self.m_min_flats:
            log.error("Not enough number (%s) of flat frames (>%s) to compute \
            master tw-flat",nframes, self.m_min_flats)
            raise Error("Flat sequence is too short, at least %s frames are required"%self.m_min_flats)
        
        if not os.path.exists(os.path.dirname(self.__output_filename)):
            log.error('Directory of combined FLAT frame does not exist')
            raise ExError('Directory of combined FLAT frame does not exist')
        if not self.__output_filename :
            log.error('Combined FLAT frame not defined')
            raise ExError('Combined FLAT frame not defined')
    
        
        # Change to the source directory
        base, infile   = os.path.split(self.__output_filename)
        iraf.chdir(base)
    
        
        # STEP 1: Check the  TYPE(twilight) and FILTER,READEMODE of each Flat frame
        # If any frame on list missmatch the FILTER, then the master twflat will be aborted
        # EXPTIME do not need be the same, so EXPTIME scaling will be done
        f_expt = -1
        f_type = ''
        f_filter = ''
        f_ncoadds = -1
        f_readmode = ''
        good_frames = []
        
        for iframe in framelist:
            f = datahandler.ClFits(iframe)
            log.debug("Checking data compatibility (filter, texp, type)")
            log.debug("Flat frame '%s' - EXPTIME= %f TYPE= %s FILTER= %s" %(iframe, f.expTime(),f.getType(), f.getFilter()))
            #Compute the mean count value in chip to find out good frames (enough check ??)
            mean = 0
            myfits = fits.open(iframe, ignore_missing_end=True)
            if f.mef==True:
                log.debug("Found a MEF file")
                #log.error("Sorry, MEF files are not supported yet !")
                #raise Exception('Sorry, MEF files are not supported yet !')
                try:
                    for i in range(1,f.next+1):
                        mean+=numpy.mean(myfits[i].data)
                    mean/=f.next
                    log.debug("MEAN value of MEF = %d", mean)
                except Exception,e:
                    raise e
            else:
                myfits = fits.open(iframe, ignore_missing_end=True)
                mean = numpy.mean(myfits[0].data)
                log.debug("MEAN value of MEF = %d", mean)
                
            myfits.close()            
            if ( f_expt!=-1 and (f.getFilter()!=f_filter or f.getType()!=f_type or f.getReadMode()!=f_readmode)) :
                log.error("Task 'createMasterTwFlat' Found a FLAT frame with different FILTER or TYPE.")
                raise Exception("Error, frame %s has different FILTER or TYPE" %(iframe))
                #continue
            else: 
                f_expt = f.expTime()
                f_filter = f.getFilter()
                f_readmode = f.getReadMode()
                if f.isTwFlat():
                    f_type = f.getType()
                else:
                    # Added (temporal) support for DOME_FLAT series (not lamp on/off)
                    log.warning("Frame %s does not look like a TwiLight Flat field" %(iframe))
                    f_type = f.getType(False)
                    #log.error("Error, frame %s does not look a TwiLight Flat field" %(iframe))
                    #raise Exception("Error, frame %s does not look a TwiLight Flat field" %(iframe))
            
            # STEP 1.1: Check either over or under exposed frames
            log.debug("Flat frame '%s' - FILTER=%s EXPTIME= %f TYPE= %s Mean= %f" %(iframe, f_filter, f_expt, f_type, mean))
            if mean>self.m_lthr and mean<self.m_hthr:
                good_frames.append(iframe)
            else:
                log.error("Frame %s skipped, either over or under exposed" %(iframe))
            
            
        if len(good_frames)>self.m_MIN_N_GOOD:
            log.info('Found %d flat frames with same filter [%s] and type:\n', len(good_frames), f_filter)
            for e in good_frames:
                log.info("--->%s",e)            
        else:
            log.error("Error, not enough good frames, exiting....")
            raise Exception("Error, not enough good flat frames")
                
        #Clobber existing output images
        iraf.clobber = 'yes'
        
        # STEP 2: We subtract a proper MASTER_DARK, it is required for TWILIGHT
        # FLATS because they might have diff EXPTIMEs
        # Prepare input list on IRAF string format
            
        log.debug("Start Dark subtraction. Master Dark -> %s"%self.__master_dark)   
        # Read Master Dark Model
        try:
            cdark = datahandler.ClFits( self.__master_dark )
            mdark = fits.open(self.__master_dark, ignore_missing_end=True)
            mdark_textp = cdark.expTime()
            # MASTER_DARK_MODEL is required !!!
            if not cdark.isMasterDarkModel() and (mdark_textp < f_expt-0.05 or mdark_textp>mdark_textp+0.05):
                log.error("File %s does not look a Master Dark Model"%self.__master_dark)
                raise Exception("Cannot find a scaled dark to apply")
        except Exception, e:
            mdark.close()
            raise e
        
        #Check MEF compatibility
        if not f.mef==cdark.mef:
            log.error("Type mismatch with MEF files")
            mdark.close()
            raise Exception("Type mismatch with MEF files")
        if f.mef:
            next = f.next # number of extension
        else:
            next = 0
                
        #t_dark = cdark.expTime()
        fileList = []
        for iframe in good_frames:
            # Remove old dark subtracted flat frames
            my_frame = self.__temp_dir + "/" + os.path.basename(iframe.replace(".fits","_D.fits"))
            misc.fileUtils.removefiles(my_frame)
            
            log.debug("Scaling master dark")
            # Build master dark with proper (scaled) EXPTIME and subtract ( I don't know how good is this method of scaling !!!)
            f = fits.open(iframe, ignore_missing_end=True)
            t_flat = datahandler.ClFits ( iframe ).expTime()
            #pr_mdark = (numpy.array(mdark[0].data, dtype=numpy.double)/float(mdark[0].header['EXPTIME']))*float(f[0].header['EXPTIME'])
            if next>0:
                for i in range(1,next+1):
                    if cdark.isMasterDarkModel():
                        log.info("Scaling MASTER_DARK_MODEL")
                        scaled_dark = mdark[i].data[1]*t_flat + mdark[i].data[0]
                    else:
                        log.info("Scaling MASTER_DARK")
                        scaled_dark = (mdark[i].data / mdark_textp ) * t_flat
                        
                    log.info("AVG(scaled_dark)=%s"%numpy.mean(scaled_dark))
                    f[i].data = f[i].data - scaled_dark
                    #f[i].data = f[i].data - mdark[i].data*float(t_flat/t_dark)
                    f[i].header.add_history('Dark subtracted %s (scaled)'
                                             %os.path.basename(self.__master_dark))
            else:
                scaled_dark = mdark[0].data[1]*t_flat + mdark[0].data[0]
                log.info("AVG(scaled_dark)=%s"%numpy.mean(scaled_dark))
                f[0].data = f[0].data - scaled_dark
                #f[0].data = f[0].data - mdark[0].data*float(t_flat/t_dark)
                f[0].header.add_history('Dark subtracted %s (scaled)' 
                                        %os.path.basename(self.__master_dark))    
                log.debug("Dark subtraction (scaled) done")
                
            #a=numpy.reshape(f[0].data, (2048*2048,))
            #print "MODE=", 3*numpy.median(a)-2*numpy.mean(a)
            #print "MEAN=" , numpy.mean(a)
            
            # Write output to outframe (data object actually still points to input data)
            try:
                f.writeto(my_frame, output_verify='fix')#'ignore')
                log.debug("Writing %s"%my_frame)
            except IOError:
                raise ExError('Cannot write output to %s' % my_frame)
                     
            f.close()
            fileList.append(my_frame)
                    
        # STEP 3: Make the combine of dark subtracted Flat frames scaling by 'mode'
        # - Build the frame list for IRAF
        log.debug("Combining dark subtracted Twilight flat frames...")
        comb_flat_frame = (self.__temp_dir + "/comb_tw_flats.fits").replace("//","/")
        misc.fileUtils.removefiles(comb_flat_frame)
        misc.utils.listToFile(fileList, self.__temp_dir + "/twflat_d.list")
        # - Call IRAF task
        # Combine the images to find out the Tw-Flat using sigma-clip algorithm;
        # the input images are scaled to have a common mode, the pixels containing 
        # objects are rejected by an algorithm based on the measured noise (sigclip),
        # and the flat-field is obtained by a median.
        iraf.mscred.flatcombine(input=("'"+"@"+self.__temp_dir+"/twflat_d.list"+"'").replace('//','/'),
                        output=comb_flat_frame,
                        combine='median',
                        ccdtype='',
                        process='no',
                        reject='sigclip',
                        subset='no',
                        scale='mode')
                        #verbose='yes'
                        #scale='exposure',
                        #expname='EXPTIME'
                        #ParList = _getparlistname ('flatcombine')
                        #)
        log.debug("Dark subtracted Twilight flat frames COMBINED")
        # Remove the dark subtracted frames
        for ftr in fileList: misc.fileUtils.removefiles(ftr)
        
        # STEP 3b (optional)
        #Median smooth the master flat
        if self.__median_smooth:
            log.debug("Doing Median smooth of FF ...")
            iraf.mscmedian(
                    input=comb_flat_frame,
                    output=comb_flat_frame.replace(".fits","_smooth.fits"),
                    xwindow=20,
                    ywindow=20,
                    outtype="median"
                    )
            shutil.move(comb_flat_frame.replace(".fits","_smooth.fits"), 
                        comb_flat_frame)
            
        #Or using scipy ( a bit slower then iraf...)
        #from scipy import ndimage
        #filtered = ndimage.gaussian_filter(f[0].data, 20)
        
        # STEP 4: Normalize the flat-field (if MEF, normalize wrt chip 1)
        # Compute the mean of the image
        if self.__normal:
            f = fits.open(comb_flat_frame, ignore_missing_end=True)
            if next>0:
                ##chip = 1 # normalize wrt to mode of chip 1
                # Because PANIC Chip1 is bad (>25% bad pixels), we use chip 2 for normalization
                chip = 2 # normalize wrt to mode of chip 1
                naxis1 = f[chip].header['NAXIS1']
                naxis2 = f[chip].header['NAXIS2']
                offset1 = int(naxis1*0.1)
                offset2 = int(naxis2*0.1)
                median = numpy.median(f[chip].data[offset1:naxis1-offset1,
                                                    offset2:naxis2-offset2])
                msg = "Normalization of MEF master flat frame wrt chip 1. (MEDIAN=%d)"%median
            elif ('INSTRUME' in f[0].header and f[0].header['INSTRUME'].lower()=='panic'
                  and f[0].header['NAXIS1']==4096 and f[0].header['NAXIS2']==4096):
                # It supposed to have a full frame of PANIC in one single 
                # extension (GEIRS default)
                median = numpy.median(f[0].data[2048-200:4096-200,200:2048-200])
                msg = "Normalization of (full) PANIC master flat frame wrt chip 1. (MEDIAN=%d)"%median
            else:
                # Not MEF, not PANIC full-frame, but could be a PANIC subwindow
                naxis1 = f[0].header['NAXIS1']
                naxis2 = f[0].header['NAXIS2']
                offset1 = int(naxis1*0.1)
                offset2 = int(naxis2*0.1)
                median = numpy.median(f[0].data[offset1:naxis1-offset1,
                                                    offset2:naxis2-offset2])
                msg = "Normalization of master (O2k?) flat frame. (MEDIAN=%d)"%median 
 

            f.close()
            log.debug(msg)
            
            # Cleanup: Remove temporary files
            misc.fileUtils.removefiles(self.__output_filename)
            # Compute normalized flat
            iraf.mscred.mscarith(operand1=comb_flat_frame,
                    operand2=median,
                    op='/',
                    pixtype='real',
                    result=self.__output_filename.replace("//","/"),
                    )
        else:
            shutil.move(comb_flat_frame, self.__output_filename)
        
        # Change back to the original working directory
        iraf.chdir()
        
        flatframe = fits.open(self.__output_filename,'update', 
                                ignore_missing_end=True)
        if self.__normal: 
            flatframe[0].header.add_history('Computed normalized master twilight flat')
            flatframe[0].header.add_history(msg)
        else: 
            flatframe[0].header.add_history('Computed master (not normalized) twilight flat')
        
        # Combined files is already added by IRAF:imcombine
        #flatframe[0].header.add_history('Twilight files: %s' %good_frames)
        
        # Add new keywords (PAPITYPE, PAPIVERS, IMAGETYP)
        flatframe[0].header.set('PAPITYPE',
                                   'MASTER_TW_FLAT',
                                   'TYPE of PANIC Pipeline generated file')
        flatframe[0].header.set('PAPIVERS',
                                   __version__,
                                   'PANIC Pipeline version')
        flatframe[0].header.set('IMAGETYP',
                                   'MASTER_TW_FLAT',
                                   'TYPE of PANIC Pipeline generated file') 
        if 'PAT_NEXP' in flatframe[0].header:
            flatframe[0].header.set('PAT_NEXP',
                                   1,
                                   'Number of positions into the current dither pattern')
        flatframe.close(output_verify='ignore') # This ignore any FITS standar violation and allow write/update the FITS file
        
        log.debug(t.tac())
        log.debug('Saved master TW_FLAT to %s', self.__output_filename )
    
        return self.__output_filename
        
        


def makemasterframe(list_or_array):
#from http://www.astro.ucla.edu/~ianc/python/_modules/ir.html#baseObject
    """
    If an array is passed, return it.  Otherwise, return a
    median stack of the input filename list.
    """
        
    if hasattr(list_or_array, 'shape') and len(list_or_array.shape)>1:
        masterframe = np.array(list_or_array, copy=False)
    else:
        masterframe = np.median(map(fits.getdata, list_or_array), axis=0)

    return masterframe
        
        

################################################################################
# main
if __name__ == "__main__":
    # Get and check command-line options
    usage = "usage: %prog [options]"  
    desc = """This module receives a series of Twilight Flats and
and a Master Dark Model and then creates a Master Twilight Flat-Field.
Note: Dome Flats series (not lamp ON/OFF) are also supported.
"""
    parser = OptionParser(usage, description = desc)
    
                  
    parser.add_option("-s", "--source",
                  action="store", dest="source_file_list", type='str',
                  help="Source file list of data frames."
                  " It can be a file or directory name.")
    
    parser.add_option("-d", "--master_dark_model", type='str',
                  action="store", dest="master_dark",
                  help="Master dark model to subtract each raw flat"
                  " (it will be scaled by TEXP)")
    
    parser.add_option("-o", "--output", type='str',
                  action="store", dest="output_filename", 
                  help="Final coadded output image")
    
    ## -optional
    
    parser.add_option("-b", "--master_bpm", type='str',
                  action="store", dest="master_bpm",
                  help="Bad pixel mask to be used (optional)", default=None)
    
    parser.add_option("-N", "--normalize",
                  action="store_true", dest="normalize", default=False,
                  help="Normalize master flat by median. If image is "
                  "multi-detector, then normalization wrt chip 1 is done)"
                  " [default=%default]")
    
    parser.add_option("-m", "--median_smooth",
                  action="store_true", dest="median_smooth", default=False,
                  help="Median smooth the combined flat-field [default=%default]")
    
    parser.add_option("-L", "--low", type='float', default=1000,
                  action="store", dest="minlevel", 
                  help="Flats with median level bellow are rejected "
                  "[default=%default].")
    
    parser.add_option("-H", "--high", type='float', default=100000,
                  action="store", dest="maxlevel", 
                  help="Flats with median level above are rejected "
                  "[default=%default].")
    
    (options, args) = parser.parse_args()
    
    
    if len(sys.argv[1:])<1:
       parser.print_help()
       sys.exit(0)
       
    if not options.source_file_list or not options.output_filename or not options.master_dark:
        parser.print_help()
        parser.error("incorrect number of arguments.")
    
    
    # Start proceduce
    
    filelist=[line.replace( "\n", "") 
              for line in fileinput.input(options.source_file_list )]
    try:
        mTwFlat = MasterTwilightFlat(filelist, options.master_dark,
                                     options.output_filename, options.minlevel,
                                     options.maxlevel,
                                     options.master_bpm,
                                     options.normalize,
                                     "/tmp",
                                     median_smooth=options.median_smooth)
        mTwFlat.createMaster()
    except Exception,e:
        log.error("Unexpected error: %s", str(e))
        sys.exit(1)
    
        
