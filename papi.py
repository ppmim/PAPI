#! /usr/bin/env python

# Copyright (c) 2009 Jose M. Ibanez. All rights reserved.
# Institute of Astrophysics of Andalusia, IAA-CSIC
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
# PAPI (PAnic PIpeline)
#
# papi.py
#
# Last update 04/March/2010
#
################################################################################


    ####################################################################
    #                                                                  #
    #                   This is where it all starts.                   #
    #         The body of this program is at the end of this file      #
    #                                                                  #
    ####################################################################

__date__ = "$Date: 2010-03-05 11:48:08 +0100 (Fri, 05 Mar 2010) $"
__author__ = "$Author: panic $"
__revision__ = "$Rev: 21 $"

    
#From system
import sys
import os
import os.path
from optparse import OptionParser
import fileinput
import glob
import shutil


# IRAF packages
from pyraf import iraf
from iraf import noao
from iraf import imred
from iraf import ccdred

# Math module for efficient array processing
import numpy

#Log
import misc.paLog
from misc.paLog import log    

#PAPI packages 
import datahandler
import reduce
import misc.fileUtils
import PAPI.linkSources as papi
import misc.utils
from reduce.makeobjmask import *
import reduce.imtrim
import reduce.astrowarp
import misc.mef 
import astromatic


class MEF_ReductionSet:
    """ 
    This class implement a reduction set of multi-extension-FITS and allow 
    all the basic data reduction operations over them. It also works with simple
    FITS files.
    
    TODO: parallel data reduction
    """
    
    def __init__(self, sci_list, out_dir=None, out_file="out.fits", dark=None, flat=None, bpm=None, file_n=0):
        """ Initialization """
        
        # Data set definition values
        self.science_list= sci_list  # the list containing the MEF science data filenames to reduce
        self.out_dir = out_dir       # directory where all output will be written
        self.out_file = out_file     # output filename for result frame (a fits file)
        self.master_dark = dark      # the master dark (filename) to use (input)
        self.master_flat = flat      # the master flat (filename) to use (input)
        self.bpm = bpm               # the master Bad Pixel Mask (filename) to use (input)
        self.file_n = file_n         # the file position of file to which run a specific funtion (nearSkySubtraction,....)
        
        # MEF-concerning 
        self.nExt = -1               # number of extensions (= number of ReductionSet's)
        self.l_red_set = []          # list of ReductionSet's. Filled in by the split() method
        
        # Checking output directory 
        if out_dir==None:
            try:
                os.mkdir(os.getcwd()+"/papi_out")
            except:
                log.error("Error creating output directory %s:", os.getcwd()+"/papi_out")
                raise
            self.out_dir=os.getcwd()+"/papi_out"
        elif not os.path.isdir(out_dir):
            try:
                os.mkdir(out_dir)
            except:
                log.error("Directory %s doesn't exist, error while creating.", out_dir)
                raise
                  
    def split(self):
        """ 
        Split the data (science & calibration) into N 'ReductionSet', where N is the number 
        of extension of the Multi-Extension FITS
        """
              
        #Source files 
        #load file list from file
        #sci_files=[line.replace( "\n", "") for line in fileinput.input(self.science_list)]
        sci_files = self.science_list
        
        sources=[]
        darks=[]
        flats=[]
        bpms=[]    
        # First, we need to check if we have MEF files
        if not datahandler.ClFits( sci_files[0] ).mef:
            self.nExt=1
            sources = sci_files
            darks.append(self.master_dark)
            flats.append(self.master_flat)
            bpms.append(self.bpm)
            rs = ReductionSet(sci_files, out_dir=self.out_dir, out_file=self.out_file, obs_mode="dither", \
                                dark=self.master_dark, flat=self.master_flat, bpm=self.bpm)          
            self.l_red_set.append( rs )
        else:            
            #Suppose we have MEF files ...
            kws=['DATE','OBJECT','DATE-OBS','RA','DEC','EQUINOX','RADECSYS','UTC','LST','UT','ST','AIRMASS','IMAGETYP','EXPTIME','TELESCOP','INSTRUME','MJD-OBS','FILTER2']    
            try:
                mef = misc.mef.MEF(sci_files)
                (self.nExt,new_sci_files)=mef.doSplit(".Q%02d.fits", out_dir=self.out_dir, copy_keyword=kws)
                n=0
                if self.master_dark!=None:
                    mef = misc.mef.MEF([self.master_dark])
                    (n,new_dark_files)=mef.doSplit(".Q%02d.fits", out_dir=self.out_dir, copy_keyword=kws)
                if self.master_flat!=None:
                    mef = misc.mef.MEF([self.master_flat])
                    (n,new_flat_files)=mef.doSplit(".Q%02d.fits", out_dir=self.out_dir, copy_keyword=kws)
                if self.bpm!=None:
                    mef = misc.mef.MEF([self.bpm])
                    (n,new_bpm_files)=mef.doSplit(".Q%02d.fits", out_dir=self.out_dir, copy_keyword=kws)
            except Exception,e:
                log.debug("Some error while splitting data set ...%s",str(e))
                raise
            # now, generate the new output filenames        
            for n in range(1,self.nExt+1):
                sources.append([self.out_dir+"/"+os.path.basename(file.replace(".fits",".Q%02d.fits"%n)) for file in sci_files])
                if self.master_dark: darks.append(self.out_dir+"/"+os.path.basename(self.master_dark.replace(".fits",".Q%02d.fits"%n)))
                else: darks.append(None)
                if self.master_flat: flats.append(self.out_dir+"/"+os.path.basename(self.master_flat.replace(".fits",".Q%02d.fits"%n)))
                else: flats.append(None)
                if self.bpm: bpms.append(self.out_dir+"/"+os.path.basename(self.bpm.replace(".fits",".Q%02d.fits"%n)))
                else: bpms.append(None)
                """
                for f in new_file_names:
                    #if re.search(".*(\.Q01)(.fits)$", f):
                        sources.append(f)
                """
                # Create the ReductionSet
                rs = ReductionSet(sources[n-1], out_dir=self.out_dir, out_file=self.out_dir+"/out_Q%02d.fits"%n, obs_mode="dither", \
                                  dark=darks[n-1], flat=flats[n-1], bpm=bpms[n-1], red_mode="single")
                self.l_red_set.append( rs )
                    
        return self.nExt,sources,darks,flats,bpms
    
    def subtractNearSky(self, fn=-1, out_filename="/tmp/skysub.fits"):
        """ run a near-sky subtraction to a give frame (fn) from the sci frame list"""
        
        ############
        self.split()
        ############
        
        outs=[]
        if fn==-1: m_fn = self.file_n
        else: m_fn = fn
        
        log.debug("# RedSet = %d",len(self.l_red_set))
        log.debug("# File pos = %d", m_fn)
        for rs in self.l_red_set:
            try:
                outs.append(rs.subtractNearSky(m_fn))
            except:
                log.error("Error while MEF near sky subtraction")
                raise
        
        #Package results from each extension into a MEF file (only if nExt>1)
        if len(self.l_red_set)>1:
            mef=misc.mef.MEF(outs)
            mef.createMEF(out_filename)
        else:
            shutil.move(outs[0], out_filename) 
        
        return out_filename
        
    def doReduction(self, red_mode="single"):
        """ Do the data reduction of all MEF's , in principle sequencially """
                                
        ############
        self.split()
        ############
    
        outs=[]
        for rs in self.l_red_set:
            try:
                # Here, we should programm the parallel data reduction 
                outs.append(rs.reduce(red_mode))
            except Exception, e:
                log.error("Error while MEF data reduction: %s",str(e))
                raise e
            
        #Package results from each extension into a MEF file (only if nExt>1)
        if len(self.l_red_set)>1:
            #log.debug("*** Creating MEF file with all outputs from reduction....***")
            #mef=misc.mef.MEF(outs)
            #mef.createMEF(self.out_file)
            # other option, do a SWARP to register the N-extension into one wide-single extension
            log.debug("*** Coadding overlapped files....")
            swarp = astromatic.SWARP()
            swarp.config['CONFIG_FILE']="/disk-a/caha/panic/DEVELOP/PIPELINE/PANIC/trunk/config_files/swarp.conf"
            swarp.ext_config['COPY_KEYWORDS']='OBJECT,INSTRUME,TELESCOPE,IMAGETYP,FILTER,FILTER2,SCALE,MJD-OBS'
            swarp.ext_config['IMAGEOUT_NAME']= self.out_file
            swarp.ext_config['WEIGHTOUT_NAME']=self.out_file.replace(".fits",".weight.fits")
            swarp.ext_config['WEIGHT_TYPE']='MAP_WEIGHT'
            swarp.ext_config['WEIGHT_SUFFIX']='.weight.fits'
            swarp.run(outs, updateconfig=False, clean=False)
        else:
            shutil.move(outs[0], self.out_file)
        
        log.info("*** Final reduced file %s created. Congratulations !! ***",self.out_file)
        
        return self.out_file
                   
    def doQuickReduction(self):
        """ Do the data quick reduction of all MEF's , in principle sequencially """
                                
        return self.doReduction(red_mode="single")
                
class ReductionSet:
    def __init__(self, sci_filelist, out_dir, out_file, obs_mode, dark=None, flat=None, bpm=None, red_mode="single"):
        """ Init function """
        
        # Input values
        self.sci_filelist = sci_filelist # list containing the science data filenames to reduce
        self.out_dir   = out_dir   # directory where all output will be written
        self.out_file  = out_file  # final reduced data file (out)
        self.obs_mode  = obs_mode  # observing mode (dither, dither_on_off, dither_off_on....)
        self.master_dark = dark    # master dark to use (input)
        self.master_flat = flat    # master flat to use (input)
        self.bpm = bpm             # master Bad Pixel Mask to use (input)
        self.red_mode = red_mode   # reduction mode (single=for QL, full=for science) 
        
        # Environment variables
        self.m_terapix_path = os.environ['TERAPIX']
        self.m_papi_path = os.environ['PAPI_HOME']
        
        
        self.m_LAST_FILES = []   # Contain the files as result of the last processing step (science processed frames)
        self.m_rawFiles = []     # Raw files (originals in the working directory)
        self.m_filter = ""       # Filter of the current data set (m_LAST_FILES)
        self.m_type = ""         # Type (dark, flat, object, ...) of the current data set; should be always object !
        self.m_expt = 0.0        # Exposition Time of the current data set files
        self.m_ncoadd = 0        # Number of coadds of the current data set files
        self.m_itime  = 0.0      # Integration Time of the currenct data set files
        self.m_n_files = ""      # Number of file in the data set

        self.MAX_MJD_DIFF = 6.95e-3  # Maximun seconds (600secs aprox) of temporal distant allowed between two consecutive frames 
        self.MIN_SKY_FRAMES = 5  # minimun number of sky frames required in the sliding window for the sky subtraction
        
          
    def checkData(self, chk_filter=True, chk_type=True, chk_expt=True, chk_itime=True, \
                  chk_ncoadd=True, chk_cont=True):
        """
        Return true is all files in file have the same filter and/or type and/or
        expt; false otherwise.
        Also check the temporal continuity (distant between two consecutive frames),
        if exceeded return False 
        
        \return True or False
        """
        
        f = datahandler.ClFits(self.m_LAST_FILES[0])
        
        filter_0 = f.getFilter()
        type_0 = f.getType()
        expt_0 = f.expTime()
        itime_0 = f.getItime()
        ncoadd_0 = f.getNcoadds()
        
        self.m_filter = filter_0
        self.m_type = type_0
        self.m_expt = expt_0
        self.m_itime = itime_0
        self.m_ncoadd = ncoadd_0
        
        mismatch_filter=False
        mismatch_type=False
        mismatch_expt=False
        mismatch_itime=False
        mismatch_ncoadd=False
        mismatch_cont=False
        
        prev_MJD=-1
        for file in self.m_LAST_FILES:
            fi=datahandler.ClFits( file )
            if chk_filter and not mismatch_filter: 
                if fi.getFilter() != filter_0:
                    log.debug("File %s does not match file FILTER", file)
                    mismatch_filter=True
            if chk_type and not mismatch_type: 
                if fi.getType() != type_0:
                    log.debug("File %s does not match file TYPE", file)
                    mismatch_type=True
            if chk_expt and not mismatch_expt: 
                if fi.expTime() != expt_0:
                    log.debug("File %s does not match file EXPTIME", file)
                    mismatch_expt=True
            if chk_itime and not mismatch_itime: 
                if fi.getItime() != itime_0:
                    log.debug("File %s does not match file ITIME", file)
                    mismatch_itime=True
            if chk_ncoadd and not mismatch_ncoadd: 
                if fi.getNcoadds() != ncoadd_0:
                    log.debug("File %s does not match file NCOADD", file)
                    mismatch_ncoadd=True
            if chk_cont and not mismatch_cont:
                if prev_MJD!=-1 and (fi.getMJD()-prev_MJD)>self.MAX_MJD_DIFF:
                    log.error("Maximmun time distant between two consecutives frames exceeded !!")
                    mismatch_cont=True
                else:
                    prev_MJD=fi.getMJD()
                    
        if mismatch_filter or mismatch_type or mismatch_expt or  mismatch_itime or mismatch_ncoadd or mismatch_cont:
            log.error("Data checking found a mismatch ....check your data files....")
            raise Exception("Error while checking data")
            #return False             
        else:    
            log.debug("All files match same file filter")
            return True
            
    def checkFilter(self):
        """
        Return true is all files in file have the same filter type, false otherwise
        
        \return True or False
        """
        f = datahandler.ClFits(self.m_LAST_FILES[0])
        filter_0 = f.getFilter()
        self.m_filter = filter_0
        for file in self.m_LAST_FILES:
            fi=datahandler.ClFits( file )
            if fi.getFilter() != filter_0:
                log.debug("File %s does not match file filter", file)
                return False
            
        log.debug("All files match same file filter")
        return True
        
        
    def checkType(self, type_to_check=None):
        """
        Return true is all files in file have the same type(science, dark, flat, ...), false otherwise
        
        \param type_to_check (\c string) type to check to
        \return True or False
        """
        
        if type_to_check==None:
            f = datahandler.ClFits(self.m_LAST_FILES[0])
            type_0 = f.getType()
        else:
            type_0 = type_to_check
             
        for file in self.m_LAST_FILES:
            fi = datahandler.ClFits(file)
            if fi.getType() != type_0:
                log.debug("File %s does not match file type %s", file, type_0)
                return False
            
        log.debug("All files match same file type")
        return True
                 
    def sortOutData(self, list=None):
        """
        Sort out input data files by MJD
        """
        
        dataset=[]
        if list==None:
            m_list=self.m_LAST_FILES
        else:
            m_list=list
            
        for file in m_list:
            fits=datahandler.ClFits(file) 
            dataset.append((file, fits.getMJD()))
            
        dataset=sorted(dataset, key=lambda data_file: data_file[1])          
        sorted_files=[]
        for tuple in dataset:
            sorted_files.append(tuple[0])
        
        return sorted_files                                
                     
    def getObsMode(self):
        """
        Return the type of dither sequece followed in the currect  'm_LAST_FILES' list. It could be:
            - dither (T-T-T-T-T- ....)
            - dither_on_off (T-S-T-S-T-S-T-S-T-....)
            - dither_off_on (S-T-S-T-S-T-S-T-S-....)
            - other  (non defined sequence,unknown)
            
            NOTE: To find out the dither/observing sequence then OBJECT keyword will be checked (see ClFits class)
        """
                   
        mode='other' # default
                   
        #Step 1: we suppose list file is sorted out  by MJD
        #Step 2: get the data type of the first file and check sequence starting from this file
        fits_0=datahandler.ClFits(self.m_LAST_FILES[0])
        fits_1=datahandler.ClFits(self.m_LAST_FILES[1])
        i=0
        if fits_0.isSky() and fits_1.isObject():
            mode='dither_off_on'
            # Then, we are going to suppose the sequence S-T-S-T-S- .... (dither_off_on)
            for file in self.m_LAST_FILES:
                if not i%2: #even
                    fits=datahandler.ClFits(file)
                    if not fits.isSky():
                        return 'other'
                elif i%2: #odd
                    fits=datahandler.ClFits(file)
                    if not fits.isObject():
                        return 'other'
                i=i+1         
        elif fits_0.isObject() and fits_1.isSky():
            # Then, we are going to suppose the sequence T-S-T-S-T- .... (dither_on_off)
            mode='dither_on_off'
            for file in self.m_LAST_FILES:
                if not i%2: #even
                    fits=datahandler.ClFits(file)
                    if not fits.isObject():
                        return 'other'
                elif i%2: #odd
                    fits=datahandler.ClFits(file)
                    if not fits.isSky():
                        return 'other'
                i=i+1                 
        elif fits_0.isObject() and fits_1.isObject():
            # Then, we are going to suppose the sequence T-T-T-T-T- .... (dither)
            mode='dither'
            for file in self.m_LAST_FILES:
                fits=datahandler.ClFits(file)
                if not fits.isObject():
                    return 'other'
        else:
            mode='other'
                        
        return mode
                        
    def getSkyFrames(self, list=None):
        """
        Given a list of files(data set), return the files identified as 'sky' frames in the m_LAST_FILES
        """
       
        sky_list=[]
        if list==None:
            m_list=self.m_LAST_FILES
        else:
            m_list=list
            
        for file in m_list:
            fits=datahandler.ClFits(file)
            if fits.isSky():
                sky_list.append(file)
                
        return sky_list
    
    def getDarkFrames(self, list=None):
        """
        Given a list of files (data set), return the files matching as dark frames
        sorted by TEXP.
        """
        
        
        if list==None:
            m_list=self.sci_filelist
        else: m_list=list
        
        match_list = []
        for file in m_list:
            fits = datahandler.ClFits(file)
            if fits.isDark():
                match_list.append((file, fits.expTime()))
        
        # Sort out frames by TEXP
        match_list=sorted(match_list, key=lambda data_file: data_file[1]) 
        
        return match_list
    
                         
    def skyFilter( self, list_file, gain_file, mask='nomask', obs_mode='dither' ):
        """
            For each input image, a sky frame is computed by combining a certain number of the closest images, 
            then this sky frame is subtracted to the image and the result is divided by the master flat; 
                         
            This function is a wrapper for skyfilter.c (IRDR)              
            
            INPUT
                list_file : a text file containing the suited structure in function of the observing_mode (the list shall be sorted by obs-date)
            
            OUTPUT
                The function generate a set of sky subtrated images (*.skysub.fits) and
                Return ONLY filtered images, when extended-source-mode sky frames are not included in the returned file list. 
                The out-file-list is previously ordered by obs-data.             
            
            VERSION
                1.0, 20090909 by jmiguel@iaa.es
        
            TODO: extended objects !!!!
        """               
        
        # Skyfilter parameters
        halfnsky=3
        destripe='none'
        out_files=[]
            
        if obs_mode=='dither':
            skyfilter_cmd=self.m_papi_path+'/irdr/bin/skyfilter '+ list_file + '  ' + gain_file +' '+ str(halfnsky)+' '+ mask + '  ' + destripe 
        elif obs_mode=='dither_on_off':
            skyfilter_cmd=self.m_papi_path+'/irdr/bin/skyfilteronoff '+ list_file + '  ' + gain_file +' '+ str(halfnsky)+' '+ mask + '  ' + destripe
        #elif obs_mode=='dither_on_off':
        #    skyfilter_cmd=self.m_papi_path+'/irdr/bin/skyfilteronoff '+ '/tmp/skylist_prueba.list' + '  ' + gain_file +' '+ str(halfnsky)+' '+ 'mask' + '  ' + destripe
        elif obs_mode=='dither_off_on':
            skyfilter_cmd=self.m_papi_path+'/irdr/bin/skyfilteroffon '+ list_file + '  ' + gain_file +' '+ str(halfnsky)+' '+ mask + '  ' + destripe
        elif obs_mode=='other':
            skyfilter_cmd=self.m_papi_path+'/irdr/bin/skyfilter_general '+ list_file + '  ' + gain_file +' '+ str(halfnsky)+' '+ mask + '  ' + destripe
        else:
            log.error("Observing mode not supported")
            raise
                  
        if misc.utils.runCmd( skyfilter_cmd )==1: # All was OK
            # Rename output sky-subtracted files
            for file in glob.glob(self.out_dir+'/*.fits.skysub'):
                shutil.move(file, file.replace('.fits.skysub', '.skysub.fits'))
                #out_files.append(file.replace('.fits.skysub', '.skysub.fits'))
            # Compose the output file list
            if obs_mode=='dither':
               out_files=[line.split()[0].replace('.fits', '.skysub.fits') for line in fileinput.input(list_file)]
            elif (obs_mode=='dither_on_off' or obs_mode=='dither_off_on' or obs_mode=='other'):
                out_files=glob.glob(self.out_dir+'/*.skysub.fits')
            # Sort-out data files by obs-data (actually not required when obs_mode='dither')
            out_files=self.sortOutData(out_files) 
            return out_files
        else:
            log.error("Some problem while running command %s", skyfilter_cmd) 
            return []
                                  
    def subtractNearSky (self, fn= 0):
        """
            Compute and subtract the nearest sky to the image in position 'fn' in the frame list (sci_filelist) 
                         
            This function make use of skyfilter_single.c (IRDR)              
            
            INPUT
                fn : file position in sci file list (the list is supposed to be sorted by obs-date)
            
            OUTPUT
                The function generate a sky subtrated image (*.skysub.fits)
            
            VERSION
                1.0, 20101011 by jmiguel@iaa.es
        
            TODO: extended objects !!!!
            
        """
                                                  
        near_list = self.sci_filelist
        if fn<0 or fn>len(near_list):
            log.error("Wrong frame number selected in near-sky subtraction")
            return None
        
        if len(near_list)<self.MIN_SKY_FRAMES:
            log.error("Wrong number of sky frames provided. Min number of sky frame is %d", self.MIN_SKY_FRAMES)
            return None
        
        # Create the temp list file of nearest (ar,dec,mjd) from current selected science file
        listfile=self.out_dir+"/nearfiles.list"
        utils.listToFile(near_list, listfile) 

        # Get the gain map
        gain=self.master_flat
        if not os.path.exists( gain ):
            raise Exception("Error, no gain map %s found"%gain)
            #TODO: try to compute GainMap !!!
                    
        hwidth=2
        cmd=self.m_papi_path+"/irdr/bin/skyfilter_single %s %s %d nomask none %d" %(listfile, gain, hwidth, fn)
        #Call external app skyfilter (papi)
        e=utils.runCmd( cmd )
        if e==1: # success
            return near_list[fn-1].replace(".fits", (".fits.skysub"))  
        else:
            log.error("Error while subtracting near sky")
            return None
                                                  
    def getPointingOffsets (self, images_in=None,  p_min_area=5, p_mask_thresh=2, p_offsets_file='/tmp/offsets.pap'):
        """DESCRIPTION
                Derive pointing offsets between each image using SExtractor OBJECTS (makeObjMask) and offsets (IRDR)
                
           INPUTS
                images_in        filename of list file (if = None, then use all .skysub.fits files in the out directory --> TB removed !) 
                 
                p_min_area       minimun area for a detected object (SExtractor DETECT_MINAREA)
                 
                p_mask_thresh    threshold for object detection (SExtractor DETECT_THRESH)
                      
           OUTPUTS
                offsets          two dimensional array with offsets
                
           NOTE
                Assumes that North is up and East is left
        """
            
        log.info("Starting getPointingOffsets....")
        offsets_mat=None
           
        # STEP 1: Create SExtractor OBJECTS images
        suffix='_'+self.m_filter+'.skysub.fits'
        output_list_file=self.out_dir+"/gpo_objs.pap"
        
        log.debug("Creaing OBJECTS images (SExtractor)....")

        if images_in==None:
            makeObjMask( self.out_dir+'*'+suffix , p_min_area, p_mask_thresh, output_list_file)
        elif os.path.isfile(images_in):
            makeObjMask( images_in , p_min_area, p_mask_thresh, output_list_file)
        else:
            log.error("Option not recognized !!!")    
                           
        # STEP 2: Compute dither offsets (in pixles) using cross-correlation technique ==> offsets
        #>mosaic objfiles.nip $off_err > offsets1.nip
        search_box=10 # half_width of search box in arcsec (default 10)
        offsets_cmd=self.m_papi_path+'/irdr/bin/offsets '+ output_list_file + '  ' + str(search_box) + ' >' + p_offsets_file
        if misc.utils.runCmd( offsets_cmd )==0:
            log.error ("Some error while computing dither offsets")
            return 0
        else:
            try:
                offsets_mat = numpy.loadtxt(p_offsets_file, usecols = (1,2,3)) # columns => (xoffset, yoffset, match fraction) in PIXELS  
            except IOError:
                log.debug("Any offsets read. Check images ....")   
        
        log.debug("END of getPointingOffsets")                        
        return offsets_mat
                            
    def coaddStackImages(self, input='/tmp/stack.pap', gain='/tmp/gain.fits', output='/tmp/coadd.fits', type_comb='average'):
        """
            Register the aligned-stack of dithered FITS images listed in 'input'
            file using offset specified inside, and calculate the mean plane.
            
            INPUTS:
                input   : file listing the file to coadd and the offsets between
                          each one
                          
                gain    : gain map file to use for the coaddition (it take into account the BPM)
                
                type_comb : type of combination to use (currently, only average available)
                
            OUTPUTS:
                output : coadded image (and the weight map .weight.fits)
                
            TODO: allow other kind of combination (median, ...see iraf.imcombine)
            
        """
                                                  
        log.info("Start coaddStackImages ...")                                          
        # STEP 1: Define parameters                                          
        input_file=input
        if input_file==None:
            log.error("Bad input file provided !")
            return
        
        if gain==None:
            gain_file=self.out_dir+"/gain_"+self.m_filter+".fits"
        else:
            gain_file=gain
        
        if output==None:
            output_file=self.out_dir+"/coadd_"+self.m_filter+".fits"     
        else:
            output_file=output
            
        weight_file=output_file.replace(".fits",".weight.fits")
        
        # STEP 2: Run the coadd                                           
        if type_comb=='average': # (use IRDR::dithercubemean)
            prog = self.m_papi_path+"/irdr/bin/dithercubemean "
            cmd  = prog + " " + input_file + " " + gain_file + " " + output_file + " " + weight_file 
            e=utils.runCmd( cmd )
            if e==0:
                log.debug("Some error while running command %s", cmd)
                return (None,None)
            else:
                log.debug("Succesful ending of coaddStackImages")
                return (output, weight_file)
        #elif type_comb=='median': # (use IRAF::imcombine)
                
        else: return (None,None)
    
                                              
    def createMasterObjMask( self, input_file, output_master_obj_mask ):
        """
        Create a master object mask from a input file using SExtractor and
        then dilate the object file by a certain factor to remove also the
        undetected object tails (dilate.c)
        """
                                                             
        log.info("Start createMasterObjMask....")
                                                             
        # STEP 1: create mask                                                     
        makeObjMask( input_file+"*", 5, 2.0)
        if os.path.exists(input_file+".objs"): 
            shutil.move(input_file+".objs", output_master_obj_mask)
            log.debug("New Object mask created : %s", output_master_obj_mask)
            
        # STEP 2: dilate mask (NOT DONE)
        """
        log.info("Dilating image ....(NOT DONE by the moment)")
        prog = self.m_papi_path+"/irdr/bin/dilate "
        scale = 0.5 #mult. scale factor to expand object regions; default is 0.5 (ie, make 50%% larger)
        cmd  = prog + " " + input_file + " " + str(scale)
        
        e=utils.runCmd( cmd )
        if e==0:
            log.debug("Some error while running command %s", cmd)
        else:
            log.debug("Succesful ending of createMasterObjMask")
        """        

        return output_master_obj_mask
                                        
    def makeAstrometry( self, input_file, catalog='2mass', re_grid=False):
        """
        Compute the astrometry of the given input_file
        (not used, deprecated)
        """
                                            
        if re_grid: regrid_str='regrid'
        else: regrid_str='noregrid'
                                            
        astrometry_cmd=self.m_papi_path+'/astrometry_scamp.pl '+ catalog + '  ' + regrid_str + ' ' + input_file
        if misc.utils.runCmd( astrometry_cmd )==0:
            log.error ("Some error while computing Astrometry")
            return 0
        else:
            return 1
                                            
        
    
    def cleanUpFiles(self):
        """Clean up files from the working directory, probably from the last execution"""
        """ TB reviewed """
        
        #misc.fileUtils.removefiles(self.out_dir+"/*.fits")
        misc.fileUtils.removefiles(self.out_dir+"/c_*", self.out_dir+"/dc_*",
                                   self.out_dir+"/*.nip", self.out_dir+"/*.pap" )
        misc.fileUtils.removefiles(self.out_dir+"/coadd*", self.out_dir+"/*.objs",
                                   self.out_dir+"/uparm*", self.out_dir+"/*.skysub*")
        misc.fileUtils.removefiles(self.out_dir+"/*.head", self.out_dir+"/*.list",
                                   self.out_dir+"/*.xml", self.out_dir+"/*.ldac",
                                   self.out_dir+"/*.png" )

    def prepareData(self):
        """Prepara input data files to be reduced, doing some FITS header modification and data """
        pass
    
    ############# Calibration Stuff ############################################
    def buildMasterDarks(self):
        """
        Look for dark files in the data set, group them by DIT,NCOADD and create
        the master darks, as many as found groups
        """
        
        log.debug("Creating Master darks...")
        # 1. Look for dark frames
        full_dark_list = self.getDarkFrames(self.sci_filelist)
        
        # 2. take the first group having the same TEXP
        last_texp = full_dark_list[0][1]
        outfile = self.out_dir+"/master_dark.fits" # class MasterDark will add as suffix (TEXP,NCOADD)
        group=[]
        l_mdarks=[]
        for tupla in full_dark_list:
            if tupla[1]==last_texp:
                group.append(tupla[0])
            else:
                try:
                    task=reduce.calDark.MasterDark (group, self.out_dir, outfile, texp_scale=False)
                    out=task.createMaster()
                    l_mdarks.append(out) # out must be equal to outfile
                except Exception,e:
                    log.error("Some error while creating master dark: %s",str(e))
                    raise e
                group=[]
                group.append(tupla[0])
                last_texp=tupla[1]
            if tupla==full_dark_list[-1]: # the last file in the list
                try:
                    task=reduce.calDark.MasterDark (group, self.out_dir, outfile, texp_scale=False)
                    out=task.createMaster()
                    l_mdarks.append(out) # out must be equal to outfile
                except Exception,e:
                    log.error("Some error while creating master dark: %s",str(e))
                    raise e
                
        return l_mdarks        
        
    def buildMasterDFlats(self):
        """
        Look for dark files in the data set, group them by DIT,NCOADD and create
        the master darks, as many as found groups
        """
        
    def buildMasterTwFlats(self):
        """
        Look for TwFlats files in the data set, group them by FILTER and create
        the master twflat, as many as found groups
        """
        
        log.debug("Creating Master TwFlat...")
        # 1. Look for twflat frames
        full_flat_list = self.getTwFlatFrames(self.sci_filelist)
        
        # 2. take the first group having the same TEXP
        last_texp = full_dark_list[0][1]
        outfile = self.out_dir+"/master_dark.fits" # class MasterDark will add as suffix (TEXP,NCOADD)
        group=[]
        l_mdarks=[]
        for tupla in full_dark_list:
            if tupla[1]==last_texp:
                group.append(tupla[0])
            else:
                try:
                    task=reduce.calDark.MasterDark (group, self.out_dir, outfile, texp_scale=False)
                    out=task.createMaster()
                    l_mdarks.append(out) # out must be equal to outfile
                except Exception,e:
                    log.error("Some error while creating master dark: %s",str(e))
                    raise e
                group=[]
                group.append(tupla[0])
                last_texp=tupla[1]
            if tupla==full_dark_list[-1]: # the last file in the list
                try:
                    task=reduce.calDark.MasterDark (group, self.out_dir, outfile, texp_scale=False)
                    out=task.createMaster()
                    l_mdarks.append(out) # out must be equal to outfile
                except Exception,e:
                    log.error("Some error while creating master dark: %s",str(e))
                    raise e
                
        return l_mdarks
        
    def reduce(self, red_mode):
        """ Main procedure for data reduction """
        
        log.info("###############################")
        log.info("#### Start data reduction #####")
        log.info("## MODE = %s  ##", self.red_mode)
        log.info("###############################")
        
        # set the reduction mode
        if red_mode!=None: self.red_mode=red_mode
        
        dark_flat = False
        
        # Clean old files 
        self.cleanUpFiles()
        
        # Change cwd to self.out_dir
        old_cwd=os.getcwd()
        os.chdir(self.out_dir)
        
        # Copy/link source files (file or directory) to reduce to the working directory
        if not os.path.dirname(self.sci_filelist[0])==self.out_dir:
            papi.linkSourceFiles(self.sci_filelist, self.out_dir)
            #files1=[line.replace( "\n", "") for line in fileinput.input(self.list_file)]
            self.m_LAST_FILES=[self.out_dir+"/"+os.path.basename(file_i) for file_i in self.sci_filelist]
        else:
            print "Input files already in output directory!"
            self.m_LAST_FILES=self.sci_filelist
            
        print "SOURCES=\n",self.m_LAST_FILES
        
        ######################################################
        # 0 - Some checks (filter, ....) 
        ######################################################
        log.info("**** Data ckecking ****")
        try:
            self.checkData(chk_filter=True, chk_type=False, chk_expt=True, chk_itime=True, chk_ncoadd=True)
        except:
            raise 
        ######################################################
        # 00 - Sort out data by MJD (self.m_LAST_FILES)
        ######################################################
        try:
            self.m_LAST_FILES=self.sortOutData()
        except:
            raise
        
        ######################################################
        # 000 - Find out dither mode
        ######################################################
        try:
            self.obs_mode=self.getObsMode()  # overwrite initial given observing mode
        except:
            raise
        
        log.info( "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" )
        log.info( "OBSERVING SEQUENCE DETECTED = %s", self.obs_mode)
        log.info( "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" )
        
        # keep a copy of the original file names
        self.m_rawFiles = self.m_LAST_FILES        
        
        ######################################
        # 1 - Apply dark, flat to ALL files 
        ######################################
        if self.master_dark!=None and self.master_flat!=None:
            log.info("**** Applying dark and Flat ****")
            dark_flat=True
            res = reduce.ApplyDarkFlat(self.m_LAST_FILES, self.master_dark, self.master_flat, self.out_dir)
            self.m_LAST_FILES = res.apply()
            
        ######################################    
        # 2 - Compute Super Sky Flat-Field 
        ######################################
        # - Find out what kind of observing mode we have (dither, ext_dither, ...)
        log.info('**** Computing Super-Sky Flat-Field ****')
        if self.obs_mode=="dither":
            log.debug("---> dither sequece <----")
            misc.utils.listToFile(self.m_LAST_FILES, self.out_dir+"/files.list")
            superflat = reduce.SuperSkyFlat(self.out_dir+"/files.list", self.out_dir+"/superFlat.fits", bpm=None, norm=False)
            superflat.create()
        elif self.obs_mode=="dither_on_off" or self.obs_mode=="dither_off_on" or self.obs_mode=="other":
            log.debug("----> EXTENDED SOURCE !!! <----")
            sky_list=self.getSkyFrames()
            misc.utils.listToFile(sky_list, self.out_dir+"/files.list")
            superflat = reduce.SuperSkyFlat(self.out_dir+"/files.list", self.out_dir+"/superFlat.fits", bpm=None, norm=False)
            superflat.create()                            
        else:
            log.error("Dither mode not supported")
            raise Exception("Error, dither mode not supported") 
              
              
        ######################################    
        # 3 - Compute Gain map and apply BPM
        ######################################
        log.info("**** Computing gain-map ****")
        gainfile = self.out_dir+'/gain_'+self.m_filter+'.fits'
        #nxblock=16
        #nyblock=16
        #The next values are to find out bad pixels 
        #nsig=5
        #mingain=0.7
        #maxgain=1.3
        g=reduce.calGainMap.GainMap(self.out_dir+"/superFlat.fits", gainfile)
        g.create() 
           
        ########################################
        # Add external Bad Pixel Map to gainmap
        ########################################     
        if self.bpm !=None:
            if not os.path.exists( self.bpm ):
                print('No external Bad Pixel Mask found. Cannot find file : "%s"' %self.bpm)
            else:
                iraf.imarith(operand1=gainfile,
                  operand2=self.bpm,
                  op='*',
                  result=gainfile.replace(".fits","_bpm.fits"),
                  verbose='yes'
                  )
                #replace old gain file
                os.rename(gainfile.replace(".fits","_bpm.fits"), gainfile)
                        
        #########################################
        # 4 - First Sky subtraction (IRDR) - sliding window technique
        #########################################
        log.info("**** 1st Sky subtraction (without object mask) ****")
        misc.utils.listToFile(self.m_LAST_FILES, self.out_dir+"/skylist1.list")
        # return only filtered images; in exteded-sources, sky frames  are not included 
        self.m_LAST_FILES=self.skyFilter(self.out_dir+"/skylist1.list", gainfile, 'nomask', self.obs_mode)      
                        
        #########################################
        # 5 - Quality assessment (FWHM, background, ellipticity, PSF quality)  
        #########################################
                            
        log.info("**** Quality Assessment: TODO ****")                   
                                                  
        #########################################
        # 6 - Compute dither offsets from the first sky subtracted/filtered images using cross-correlation
        #########################################
        misc.utils.listToFile(self.m_LAST_FILES, self.out_dir+"/files_skysub.list")
        offset_mat=self.getPointingOffsets (self.out_dir+"/files_skysub.list", 5, 3, self.out_dir+'/offsets1.pap')                
                        
        
        #########################################
        # 7 - First pass coaddition using offsets
        #########################################
        log.info("**** Initial coaddition of sky subtracted frames ****")
        fo=open(self.out_dir+'/offsets1.pap',"r")
        fs=open(self.out_dir+'/stack1.pap','w+')
        for line in fo:
          n_line = line.replace(".skysub.fits.objs", ".skysub.fits") 
          fs.write( n_line )
        fo.close()
        fs.close()    
        self.coaddStackImages(self.out_dir+'/stack1.pap', None, self.out_dir+'/coadd1.fits','average')
    
        ## END OF SINGLE REDUCTION  ##
        if self.obs_mode!='dither' or self.red_mode=="single":
            log.info("**** Doing Astrometric calibration of coadded result frame ****")
            reduce.astrowarp.doAstrometry(self.out_dir+'/coadd1.fits', self.out_file, "2MASS" ) 
            log.info("Generated output file ==>%s", self.out_file)
            log.info("#########################################")
            log.info("##### End of SINGLE data reduction ######")
            log.info("#########################################")
            return self.out_file 
        
        log.info("************************")
        log.info(" START SECOND PASS      ")
        log.info("************************")
        
        #########################################
        # 8 - Create master object mask
        #########################################
        log.info("**** Master object mask creation ****")
        obj_mask=self.createMasterObjMask(self.out_dir+'/coadd1.fits', self.out_dir+'/masterObjMask.fits')  
        
        ###########################################################
        # 9 - Second Sky subtraction (IRDR) using then OBJECT MASK
        ###########################################################
        log.info("**** Sky subtraction with 2nd object mask ****")
        #Compound masked sky list
        fs=open(self.out_dir+"/skylist2.pap","w+")
        i=0
        j=0
        for file in self.m_rawFiles:
            if dark_flat: 
                line = file.replace(".fits","_D_F.fits") + " " + obj_mask + " " + str(offset_mat[j][0]) + " " + str(offset_mat[j][1])
            else:
                line = file + " " + obj_mask + " " + str(offset_mat[j][0]) + " " + str(offset_mat[j][1])
            fs.write(line+"\n")
            if (self.obs_mode=='dither_on_off' or self.obs_mode=='dither_off_on') and i%2:
                j=j+1
            elif self.obs_mode=='dither':
                j=j+1
            i=i+1
                
        fs.close()
        self.m_LAST_FILES=self.skyFilter( self.out_dir+"/skylist2.pap", gainfile, 'mask', self.obs_mode)      
    
        #### EXIT ########
        #log.info("Sucessful end of Pipeline (I hope!)")
        #sys.exit()
        ##################
        
        #########################################
        # X1 - Compute field distortion (SCAMP internal stats)
        #########################################
        _astrowarp=False
        if _astrowarp:
            log.info("**** Astrometric calibration and stack of individual frames to field distortion correction ****")
            aw = reduce.astrowarp.AstroWarp(self.m_LAST_FILES, catalog="2MASS", coadded_file=options.output_filename)
            try:
                aw.run()
            except:
                log.error("Some error while running Astrowarp....")
                raise
        
        #### EXIT ########
        #log.info("Sucessful end of Pipeline (I hope!)")
        #sys.exit()
        ##################
            
        #########################################
        # X2 - Remove field distortion from individual images (SWARP)-regriding 
        #########################################
        
        #########################################
        # X3 - Coaddition of field distortion removed images (SWARP) 
        #########################################
    
        #########################################
        # X4 - Final Astrometric calibration (SCAMP) of the coadded image 
        #########################################
    
        #########################################################################################
        # 9 - Create second coadded image of the dithered stack using new sky subtracted frames (using the same offsets)
        #########################################################################################
        log.info("**** Coadding image free distorion frames ****")
        self.coaddStackImages(self.out_dir+'/stack1.pap', None, self.out_dir+'/coadd2.fits')
        reduce.imtrim.imgTrim(self.out_dir+'/coadd2.fits')
        
        #########################################
        # 10 - Make Astrometry
        #########################################
        log.info("**** Computing Astrometric calibration of coadded(2) result frame ****")
        reduce.astrowarp.doAstrometry(self.out_dir+'/coadd2.fits', self.out_file, "2MASS" )
        
        log.info("Generated output file ==>%s", self.out_file)
        
        os.chdir(old_cwd)
        
        log.info("##################################")
        log.info("##### End of data reduction ######")
        log.info("##################################")
        
        return self.out_file 
        
                            
################################################################################
# main
################################################################################
if __name__ == "__main__":
    log.debug( 'Start PAPI....')
    # Get and check command-line options
    
    
    usage = "usage: %prog [options] arg1 arg2 ..."
    parser = OptionParser(usage)
    
                  
    parser.add_option("-s", "--source",
                  action="store", dest="source_file_list",
                  help="Source file list of data frames. It can be a file or directory name.")
    
    parser.add_option("-o", "--output_file",
                  action="store", dest="output_filename", help="final reduced output image")
    
    parser.add_option("-d", "--outdir",
                  action="store", dest="out_dir", default="/tmp",
                  help="output dir for intermidiate files")
    
    parser.add_option("-t", "--type",
                  action="store", dest="type", default="quick", help="type of reduction (quick|science)")
                  
    parser.add_option("-m", "--obs_mode",
                  action="store", dest="obs_mode", default="dither", help="observing mode (dither|ext_dither)")
    
    
    parser.add_option("-D", "--dark",
                  action="store", dest="dark",
                  help="master dark to subtract")
    
    parser.add_option("-F", "--flat",
                  action="store", dest="flat",
                  help="master flat to divide by")

    parser.add_option("-b", "--bpm",
                  action="store", dest="bpm", help="bad pixel mask")

    
    parser.add_option("-C", "--config_file",
                  action="store_true", dest="config_file", help="config file for the data reduction process")
                  
    parser.add_option("-S", "--show",
                  action="store_true", dest="show", help="show final reduced image", default=False)
                  
    parser.add_option("-1", "--single",
                  action="store_true", dest="single", help="make a single reduction", default=False)              

    parser.add_option("-v", "--verbose",
                  action="store_true", dest="verbose", default=True,
                  help="verbose mode [default]")

                                
    (options, args) = parser.parse_args()
    
    if not options.source_file_list or not options.output_filename or len(args)!=0: # args is the leftover positional arguments after all options have been processed
        parser.print_help()
        parser.error("incorrect number of arguments " )
    if options.verbose:
        print "reading %s ..." % options.source_file_list
    
    
    # Create the MEF_RS (it works both simple FITS as MEF files)
    sci_files=[line.replace( "\n", "") for line in fileinput.input(options.source_file_list)]
    #mef_rs = MEF_ReductionSet( sci_files, options.out_dir, out_file=options.output_filename, \
    #                            dark=options.dark, flat=options.flat, bpm=options.bpm)
    #out=mef_rs.doQuickReduction()

    rs = ReductionSet( sci_files, options.out_dir, out_file=options.output_filename, obs_mode="dither", \
                                dark=options.dark, flat=options.flat, bpm=options.bpm, red_mode="single")
    out = rs.buildMasterDarks()
     
    if options.show:
        out.show()    
