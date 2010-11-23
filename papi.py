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
import tempfile
import dircache


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
import misc.utils
from reduce.makeobjmask import *
import reduce.imtrim
import reduce.astrowarp
import misc.mef 
import astromatic
import datahandler.dataset

class ReductionSet:
    def __init__(self, rs_filelist, out_dir, out_file, obs_mode, dark=None, flat=None, bpm=None, red_mode="single", classf=True):
        """ Init function """
        
        # Input values
        self.rs_filelist = rs_filelist # list containing the science data filenames to reduce
        self.out_dir   = out_dir   # directory where all output will be written
        self.out_file  = out_file  # final reduced data file (out)
        self.obs_mode  = obs_mode  # observing mode (dither, dither_on_off, dither_off_on....)
        self.master_dark = dark    # master dark to use (input)
        self.master_flat = flat    # master flat to use (input)
        self.bpm = bpm             # master Bad Pixel Mask to use (input)
        self.red_mode = red_mode   # reduction mode (single=for QL, full=for science) 
        self.classf = classf       # flag to decide if classification will be done (by OB_ID, OB_PAT, FILTER) or not (only classf by FILTER)
        self.gr_sci = []           # list of group of observing object/sequence (see getObjectSequences() ) 
        
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
        
        #DataBase (in memory)
        self.db=None
        
        print "CLASS=",classf
    
    def __initDB(self):
        """
        Initialize the data base, loading the full frame list
        
        NOTE: If we init the DB in the constructor __init__(), we will have problems
        if we call use the DB from any method of ReductionSet and from other thread.
        For more info, see http://stackoverflow.com/questions/393554/python-sqlite3-and-concurrency
        """
        #DataBase (in memory)
        self.db=datahandler.dataset.DataSet(source="file")
        self.db.createDB()
        self.db.load(self.rs_filelist)
        if self.master_dark!=None: self.db.insert(self.master_dark)
        if self.master_flat!=None: self.db.insert(self.master_flat)
        if self.bpm !=None: self.db.insert(self.bpm)
        
        self.db.ListDataSet()
        
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
    
    def split(self, frame_list):
        """ 
        Split the data in frame list (any kind, science or calibration) into N 'sub-list',
        where N is the number of extension of the Multi-Extension FITS.
        """
        
        log.debug("Starting split() method ....")
        
        new_frame_list=[] # a list of N list, where N=number of extension of the MEF 
        nExt=0
        if frame_list==None or len(frame_list)==0:
            return [],0
            
        # First, we need to check if we have MEF files
        if not datahandler.ClFits( frame_list[0] ).mef:
            nExt=1
            new_frame_list.append(frame_list)
        else:            
            #Suppose we have MEF files ...
            kws=['DATE','OBJECT','DATE-OBS','RA','DEC','EQUINOX','RADECSYS','UTC','LST','UT','ST','AIRMASS','IMAGETYP','EXPTIME','TELESCOP','INSTRUME','MJD-OBS','FILTER2']    
            try:
                mef = misc.mef.MEF(frame_list)
                (nExt, sp_frame_list)=mef.doSplit(".Q%02d.fits", out_dir=self.out_dir, copy_keyword=kws)
            except Exception,e:
                log.debug("Some error while splitting data set ...%s",str(e))
                raise
            # now, generate the new output filenames        
            for n in range(1,nExt+1):
                new_frame_list.append([self.out_dir+"/"+os.path.basename(file.replace(".fits",".Q%02d.fits"%n)) for file in frame_list])
                """
                for f in new_file_names:
                    #if re.search(".*(\.Q01)(.fits)$", f):
                        sources.append(f)
                """
        
        print "NEW_FRAME_LIST=",new_frame_list
        
        return new_frame_list, nExt
    
    def getObsMode(self):
        """
        Return the type of dither sequece followed in the currect  'm_LAST_FILES' list. It could be:
            - dither (T-T-T-T-T- ....)
            - dither_on_off (T-S-T-S-T-S-T-S-T-....)
            - dither_off_on (S-T-S-T-S-T-S-T-S-....)
            - other  (non defined sequence,unknown) (i.e,  T-S-T-T-S-T-T-S-T-....)
            
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

    def getCalibFor(self, sci_obj_list):
        """
        Given a list of frames belonging to a observing sequence for a given object (star, galaxi, whatever),
        return the most appropiate calibration files (master dark,flat,bpm) to
        reduce the sequence.
        
        Reduce 3-list of calibration files (dark, flat, bpm), even is each one has only 1 file 
        """
        obj_frame = datahandler.ClFits(sci_obj_list[0])
        # We take as sample, the first frame in the list, but all frames must
        # have the same features (expT,filter,ncoadda, readout-mode, ...)
        expTime = obj_frame.expTime()
        filter = obj_frame.getFilter()
        
        master_dark=self.db.GetFilesT('MASTER_DARK', expTime)
        master_flat=self.db.GetFilesT('MASTER_DOME_FLAT', -1, filter)
        if master_flat==[]:
            master_flat=self.db.GetFilesT('MASTER_TW_FLAT', -1, filter)
            
        master_bpm =self.db.GetFilesT('MASTER_BPM')
        
        log.debug("Found master dark %s", master_dark)
        log.debug("Found master flat %s", master_flat)
        log.debug("Found master bpm %s", master_bpm)
        
        return master_dark, master_flat, master_bpm
        
    def getDarkFrames(self):
        """
        Given a list of files (data set), return the files matching as dark frames
        sorted by TEXP.
        Return a list of tuples with (filename,tExp)
        """
        
        """
        if list==None:
            m_list=self.rs_filelist
        else: m_list=list
        
        match_list = []
        for file in m_list:
            fits = datahandler.ClFits(file)
            if fits.isDark():
                match_list.append((file, fits.expTime()))
        
        # Sort out frames by TEXP
        match_list=sorted(match_list, key=lambda data_file: data_file[1]) 
        
        return match_list
        """
        # query the DB for DARK frames 
        m_list = self.db.GetFilesT('DARK')
        
        match_list = []
        for file in m_list:
            fits = datahandler.ClFits(file)
            if fits.isDark():
                match_list.append((file, fits.expTime()))
        
        # Sort out frames by TEXP
        match_list=sorted(match_list, key=lambda data_file: data_file[1]) 
        
        return match_list # a list of tuples as (file,tExp)
    
    def getTwFlatFrames(self):
        """
        Query the DB (date set) to look for TwFlats sorted by FILTER.
        Return a list of tuples with (filename,filter)
        """

        # query the DB for TW_FLAT frames for any TEXP,FILTER 
        m_list = self.db.GetFilesT('SKY_FLAT')
        match_list = []
        for file in m_list:
            fits = datahandler.ClFits(file)
            if fits.isTwFlat():
                match_list.append((file, fits.getFilter()))
        
        
        # Sort out frames by FILTER
        match_list=sorted(match_list, key=lambda data_file: data_file[1]) 
        
        return match_list # a list of tuples as (file,filter)
    
    def getObjectSequences(self):
        """
        Query the DB (data set) to look for object/pointing sequences sorted by
        MJD; they can be grouped by:
        
            a) FILTER (it is used when no date set classification can be done)
            b) OB_ID,OB_PAT,FILTER.
            
        Return : a list of list of sequence files belonging to
        
        TBD: implement other alternative way of grouping : sort out by MJD and group by FILTER
        """
        
        if not self.classf:
            (filters, seq_list)=self.db.GetFilterFiles() # only SCIENCE FILES or SKY_FOR
            k=0
            for filter in filters:
                print "\nFILTER=%s"%filter
                print "-----------------------------------------------------------------------------------\n"
                for file in seq_list[k]:
                    print file
                k+=1
            log.debug("Found %d groups of SCI files: %s", len(filters), seq_list)
            
        else:
            (seq_par, seq_list)=self.db.GetSeqFiles(filter=None,type='SCIENCE')
            k=0
            for par in seq_par:
                print "\nSEQUENCE PARAMETERS - OB_ID=%s,  OB_PAT=%s, FILTER=%s"%(par[0],par[1],par[2])
                print "-----------------------------------------------------------------------------------\n"
                for file in seq_list[k]:
                    print file
                k+=1
        
            log.debug("Found %d groups of SCI files: %s", len(seq_par), seq_list)
        
        return seq_list
    
    def getDomeFlatFrames(self):
        """
        Query the DB (date set) to look for DomeFlats (on and off) sorted by FILTER.
        Return a list of tuples with (filename,filter)
        """

        # query the DB for DOME_FLAT frames for any TEXP,FILTER 
        m_list = self.db.GetFilesT('DOME_FLAT_LAMP_ON')
        m_list2 = self.db.GetFilesT('DOME_FLAT_LAMP_OFF')
        
        match_list = []
        for file in m_list+m_list2:
            fits = datahandler.ClFits(file)
            if fits.isDomeFlat():
                match_list.append((file, fits.getFilter()))
        
        
        # Sort out frames by FILTER
        match_list=sorted(match_list, key=lambda data_file: data_file[1]) 
        
        return match_list # a list of tuples as (file,filter)
    
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
                                  
    
    def subtractNearSky(self, near_list=None, file_pos=0, out_filename=None):
        """
            Compute and subtract the nearest sky to the image in position 'fn' in the given frame list (near_list) 
                         
            This function make use of skyfilter_single.c (IRDR)              
            
            INPUT
                file_pos : file position in sci file list (the list is supposed to be sorted by obs-date)
            
            OUTPUT
                The function generate a sky subtrated image (*.skysub.fits)
            
            VERSION
                1.1, 20101103 by jmiguel@iaa.es
        
            TODO: extended objects !!!!
            
        """
        log.debug("Start subtractNearSky")
        
        # default values
        if near_list==None: near_list = self.rs_filelist
        if out_filename==None: out_filename=self.out_file
        
        # Some previus checks
        if file_pos<0 or file_pos>len(near_list):
            log.error("Wrong frame number selected in near-sky subtraction")
            return None
        
        if len(near_list)<self.MIN_SKY_FRAMES:
            log.error("Wrong number of sky frames provided. Min number of sky frame is %d", self.MIN_SKY_FRAMES)
            return None
        
        # Get the gain map
        gain=self.master_flat
        if not os.path.exists( gain ):
            raise Exception("Error, no gain map file %s found"%gain)
            #TODO: try to compute GainMap using the given images !!!
        
        print "NEAR_LIST=",near_list    
        # 1. Split MEF file
        obj_ext, next = self.split(near_list) # it must return a list of list (one per each extension)
        out_ext = []
        
        # 2. Process each extension
        for n in range(next):
            log.debug("===> Processing extension %d", n+1)
            # Create the temp list file of nearest (ar,dec,mjd) from current selected science file
            listfile=self.out_dir+"/nearfiles.list"
            utils.listToFile(obj_ext[n], listfile)
            print "NEAR_FILES=", obj_ext[n]
            #Call external app skyfilter (papi)
            hwidth=2 ## TODO is it right ?????
            cmd=self.m_papi_path+"/irdr/bin/skyfilter_single %s %s %d nomask none %d" %(listfile, gain, hwidth, file_pos)
            e=utils.runCmd( cmd )
            if e==1: # success
                out_ext.append(obj_ext[n][file_pos-1].replace(".fits", (".fits.skysub")))  
            else:
                log.error("Some error while subtracting sky in extension #%d# ", n+1)
                raise Exception("Some error while subtracting sky in extension #%d# ", n+1)
        
        # 3. Package results back from each extension into a MEF file (only if nExt>1)
        if len(out_ext)>1:
            mef=misc.mef.MEF(out_ext)
            mef.createMEF(out_filename)
        elif len(out_ext)==1:
            shutil.move(out_ext[0], out_filename) 
        else:
            log.error("Some error while subtracting sky. Any output produced.")
            
        return out_filename
    
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
                            
    def coaddStackImages(self, input='/tmp/stack.pap', gain=None, output=None, type_comb='average'):
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
    
                                              
    def __createMasterObjMask( self, input_file, output_master_obj_mask ):
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
        
        misc.fileUtils.removefiles(self.out_dir+"/*.fits")
        misc.fileUtils.removefiles(self.out_dir+"/c_*", self.out_dir+"/dc_*",
                                   self.out_dir+"/*.nip", self.out_dir+"/*.pap" )
        misc.fileUtils.removefiles(self.out_dir+"/coadd*", self.out_dir+"/*.objs",
                                   self.out_dir+"/uparm*", self.out_dir+"/*.skysub*")
        misc.fileUtils.removefiles(self.out_dir+"/*.head", self.out_dir+"/*.list",
                                   self.out_dir+"/*.xml", self.out_dir+"/*.ldac",
                                   self.out_dir+"/*.png" )

    
    ############# Calibration Stuff ############################################
    def buildCalibrations(self):
        """
        Build the whole master calibration files from the currect calibrations files
        found in the data set (darks, flats)
        """
        
        log.debug("Start builing the whole calibration files ...")
        # If not initialized, Init DB
        if self.db==None: self.__initDB()
        master_files=[]
        try:
            master_files+=self.buildMasterDarks()
            master_files+=self.buildMasterDomeFlats()
            master_files+=self.buildMasterTwFlats()
            master_files+=self.buildMasterSuperFlats()
            master_files+=self.buildGainMaps(type="all")
        except Exception,e:
            log.error("Some error while builing master calibration files...: %s", str(e))
            raise e
        
        finally:
            log.debug("Calibration files created : %s", master_files)
            return master_files
        
    def buildGainMaps(self, type="all"):
        """
        Look for master flats (sky, twlight,dome or all) files in the data set, group them by
        FILTER and create the master gainmap, as many as found groups
        
        Return the list of gainmap created
        
        TODO: take into account the possibility to found several master flats and then combine them
        to build a gain map; at the moment only the first one found is used
        
        """
        
        log.debug("Building GainMap for %s Flats", type)
        l_gainmaps=[]
        # 1. Look for master flat frames
        full_flat_list=[]
        if type=="all":
            full_flat_list = self.db.GetFilesT(type="MASTER_SKY_FLAT", texp=-1, filter="ANY")
            full_flat_list+= self.db.GetFilesT(type="MASTER_TW_FLAT", texp=-1, filter="ANY")
            full_flat_list+= self.db.GetFilesT(type="MASTER_DOME_FLAT", texp=-1, filter="ANY")
        elif type=="sky":
            full_flat_list=self.db.GetFilesT(type="MASTER_SKY_FLAT", texp=-1, filter="ANY")
        elif type=="twlight":
            full_flat_list=self.db.GetFilesT(type="MASTER_TW_FLAT", texp=-1, filter="ANY")
        elif type=="dome":
            full_flat_list=self.db.GetFilesT(type="MASTER_DOME_FLAT", texp=-1, filter="ANY")
        else:
            log.error("Wrong type of master flat specified")
            raise Exception("Wrong type of master flat specified")
        
        if full_flat_list==None or len(full_flat_list)==0:
            log.warning("Any master flat field found")
            return []

        # 2. Group by filter
        sorted_list = []
        for file in full_flat_list:
            fits = datahandler.ClFits(file)
            sorted_list.append((file, fits.getFilter()))
        
        # Sort out frames by FILTER
        sorted_list=sorted(sorted_list, key=lambda data_file: data_file[1])
        
        # 3. take the first group having the same FILTER
        last_filter = sorted_list[0][1]
        print "SORTED LIST=", sorted_list
        group=[]
        build_master=False
        
        for tupla in sorted_list:
            if tupla[1]==last_filter:
                group.append(tupla[0])
                if tupla==sorted_list[-1]: # the last file in the list
                    build_master=True
            else: # we have found the end of the group 
                build_master=True
                
            if build_master:
                try:
                    #TODO: here we should check if we have more that one master flat, and if have, then combine them ... 
                    # generate a random filename for the master, to ensure we do not overwrite any file
                    output_fd, outfile = tempfile.mkstemp(suffix='.fits', dir=self.out_dir)
                    os.close(output_fd)
                    print "FLAT CANDIDATE=", group[0]
                    task=reduce.calGainMap.GainMap(group[0], outfile, bpm=None, do_normalization=True)
                    out=None
                    out=task.create()
                    l_gainmaps.append(out) # out must be equal to outfile
                except Exception,e:
                    log.error("Some error while creating gainmap: %s",str(e))
                    raise e
                
                # reset for new group
                group=[]
                group.append(tupla[0])
                last_filter=tupla[1]
                build_master=False
                
                if tupla==sorted_list[-1]: # the last file in the list
                    
        # insert products (gainmaps) into DB
        for f in l_gainmaps: self.db.insert(f)
        self.db.ListDataSet()  
        return l_gainmaps
    

    def buildMasterDarks(self):
        """
        Look for dark files in the data set, group them by DIT,NCOADD and create
        the master darks, as many as found groups
        """
        
        log.debug("Creating Master darks...")
        l_mdarks=[]
        # 1. Look for dark frames
        full_dark_list = self.getDarkFrames()
        if len(full_dark_list)<=0:
            log.warning("Any Dark frames found")
            return []
    
        # 2. take the first group having the same TEXP
        last_texp = full_dark_list[0][1]
        #outfile = self.out_dir+"/master_dark.fits" # class MasterDark will add as suffix (TEXP_NCOADDS)
        group=[]
        build_master=False
        for tupla in full_dark_list:
            if tupla[1]==last_texp:
                group.append(tupla[0])
                if tupla==full_dark_list[-1]: # the last file in the list
                    build_master=True
            else:
                build_master=True
            
            if build_master and len(group)>1:    
                try:
                    # generate a random filename for the master, to ensure we do not overwrite any file
                    output_fd, outfile = tempfile.mkstemp(suffix='.fits', dir=self.out_dir)
                    os.close(output_fd)
                    task=reduce.calDark.MasterDark (group, self.out_dir, outfile, texp_scale=False)
                    out=task.createMaster()
                    l_mdarks.append(out) # out must be equal to outfile
                except Exception,e:
                    log.error("Some error while creating master dark: %s",str(e))
                    raise e
                # reset for new group
                group=[]
                group.append(tupla[0])
                last_texp=tupla[1]
                build_master=False
                
        # insert products (master darks) into DB
        for f in l_mdarks: self.db.insert(f)
        self.db.ListDataSet()         
        return l_mdarks        
        
    def buildMasterDomeFlats(self):
        """
        Look for DOME FLATS files in the data set, group them by FILTER and create
        the master DomeFlat, as many as found groups (i.e., filters)
        """
        
        log.debug("Building Master DomeFlats...")
        l_mflats=[]
        # 1. Look for domeflat frames
        full_flat_list = self.getDomeFlatFrames()
        if len(full_flat_list)<=0:
            log.warning("Any DomeFlatField frames found")
            return []
            
        # 2. take the first group having the same FILTER
        last_filter = full_flat_list[0][1]
        group=[]
        build_master=False
        
        for tupla in full_flat_list:
            if tupla[1]==last_filter:
                group.append(tupla[0])
                if tupla==full_flat_list[-1]: # the last file in the list
                    build_master=True
            else: # we have found the end of the group 
                build_master=True
                
            if build_master and len(group)>1:
                try:
                    # generate a random filename for the master, to ensure we do not overwrite any file
                    output_fd, outfile = tempfile.mkstemp(suffix='.fits', dir=self.out_dir)
                    os.close(output_fd)
                    #outfile = self.out_dir+"/master_domeflat_%s.fits"%last_filter # added as suffix (FILTER)
                    task=reduce.calDomeFlat.MasterDomeFlat(group, os.path.dirname(outfile), outfile, None)
                    out=task.createMaster()
                    l_mflats.append(out) # out must be equal to outfile
                except Exception,e:
                    log.error("Some error while creating master DomeFlat: %s",str(e))
                    raise e
                
                # reset for new group
                group=[]
                group.append(tupla[0])
                last_filter=tupla[1]
                build_master=False            
                
        # insert products (master dome flats) into DB
        for f in l_mflats: self.db.insert(f)
        self.db.ListDataSet()  
        return l_mflats
    
    def buildMasterTwFlats(self):
        """
        Look for TwFlats files in the data set, group them by FILTER and create
        the master twflat, as many as found groups
        """
        
        log.debug("Building Master TwFlats...")
        l_mflats=[]
        # 1. Look for twflat frames
        full_flat_list = self.getTwFlatFrames()
        if len(full_flat_list)<=0:
            log.warning("Any TwFlatField frames found")
            return []
            
        # 2. take the first group having the same FILTER
        last_filter = full_flat_list[0][1]
        group=[]
        build_master=False
        
        for tupla in full_flat_list:
            if tupla[1]==last_filter:
                group.append(tupla[0])
                if tupla==full_flat_list[-1]: # the last file in the list
                    build_master=True
            else: # we have found the end of the group 
                build_master=True
                
            if build_master and len(group)>1:
                try:
                    master_dark=self.db.GetFilesT('MASTER_DARK') # could be > 1 master darks, then use the last(mjd sorted)
                    # if required, master_dark will be scaled in MasterTwilightFlat class
                    if len(master_dark)>0:
                        # generate a random filename for the master, to ensure we do not overwrite any file
                        output_fd, outfile = tempfile.mkstemp(suffix='.fits', dir=self.out_dir)
                        os.close(output_fd)
                        #outfile = self.out_dir+"/master_twflat_%s.fits"%last_filter # added as suffix (FILTER)
                        task=reduce.calTwFlat.MasterTwilightFlat(group, master_dark[-1], outfile, lthr=1000, hthr=100000, bpm=None)
                        out=task.createMaster()
                        l_mflats.append(out) # out must be equal to outfile
                    else:
                        # should we create master dark ??
                        log.error("MASTER_DARK not found")
                        raise Exception("MASTER_DARK not found")
                except Exception,e:
                    log.error("Some error while creating master TwFlat: %s",str(e))
                    raise e
                
                # reset for new group
                group=[]
                group.append(tupla[0])
                last_filter=tupla[1]
                build_master=False            
                
        # insert products (master twflats) into DB
        for f in l_mflats: self.db.insert(f)
        self.db.ListDataSet()  
        return l_mflats
    
    def buildMasterSuperFlats(self):
        """
        Look for science/sky files in the data set, group them by FILTER and
        temporal proximity and then create the master superFlats, as many as
        found groups.
        
        Return: the list of master superflat created
        
        TODO: need to be tested 
        """
        
        log.debug("Building Master SuperFlats...")
        
        # 1. Look for SCIENCE/SKY frames
        full_file_list = self.getObjectSequences()
        if len(full_file_list)<=0:
            log.warning("Any sequence science frames found")
            return []
            
        l_mflats=[]       
        for seq in full_file_list:
            if len(seq)>1:
                try:
                    # generate a random filename for the master super flat, to ensure we do not overwrite any file
                    output_fd, output_path = tempfile.mkstemp(suffix='.fits', dir=self.out_dir)
                    os.close(output_fd)
                    #outfile = self.out_dir+"/master_superflat_%s.fits"%last_filter # added as suffix (FILTER)
                    superflat = reduce.SuperSkyFlat(seq, output_path, bpm=None, norm=False)
                    out=superflat.create()
                    l_mflats.append(out)
                except Exception,e:
                    log.error("Some error while creating master SuperFlat: %s",str(e))
                    raise e
                
        # insert products (master SuperFlats) into DB
        for f in l_mflats: self.db.insert(f)
        self.db.ListDataSet()  
        return l_mflats # a list of master super flats created
    
    def reduceSet(self, red_mode="single"):
        """
        The main method for full set data reduction.
        
        Main steps:
        
         1. Build master calibration files
        
         2. Find out Objects/pointings/filter
        
         3. For each object/pointing/filter do :
         
            3.1 c=get_calibration_files (depend on filter/text)
            3.2 o=reduce_obj()
            3.3 exts=split(c,o)
            3.4 for each extension e do:
                out+=reduce_obj(o_e, c_e)
            3.5 create_joined_mef(out)
            3.6 Add objects to catalog (?)
            
         
        Return a list of N files produced as result of the data reduction of
        the N sequecend found.
        
        """
       
        log.info("Starting Reduction of data set ...")
        
        # Clean old files 
        self.cleanUpFiles()
        
        # Init DB
        if self.db==None: self.__initDB()
        
        if red_mode=="single": # we suppose it's Quick-Look mode, thus not use calibration files
            log.debug("Single/quick reduction mode, no calibration files will be built")
        else:
            log.debug("Building calibration for the whole files ...")
            
            if self.master_dark==None:
                log.debug("Building Master darks ...")
                try:
                    self.buildMasterDarks()
                except Exception,e:
                    log.error("Cannot build Master Darks !: %s",str(e))
                    
            if self.master_flat==None:
                log.debug("Building Master flats and gainmaps ...")
                try:     
                    self.buildMasterDomeFlats()
                except Exception,e:
                    log.error("Cannot build Master Dome Flats !: %s",str(e))
                
                try:    
                    self.buildMasterTwFlats()
                except Exception,e:
                    log.error("Cannot build Master Tw Flats !: %s",str(e))
                
                try:    
                    self.buildGainMaps(type="sky")
                except Exception,e:
                    log.error("Cannot build Sky Gain Maps !: %s",str(e))
                
                try:    
                    self.buildGainMaps(type="twlight")
                except Exception,e:
                    log.error("Cannot build TwLight Gain Maps !: %s",str(e))
                    
                try:    
                    self.buildGainMaps(type="dome")
                except Exception,e:
                    log.error("Cannot build  Gain Maps !: %s",str(e))    
        
            # TODO: and BPM ???
            
        return [] # PRUEBNA !!!!
        sequences=self.getObjectSequences()
        i=0
        out_ext=[] # it will store the partial extension reduced output filenames
        seq_result_outfile="" # will store the full-frame out filename of the current reduced sequence  
        seq_outfile_list=[] # will store the full-frame result for each sequence-data-reduction 
        # For each object/sequece, split and reduce de sequence
        for obj_seq in sequences:  #sequences is a list of list, and obj_seq is a list
            log.debug("===> Reduction of obj_seq %d",i)
            if len(obj_seq)<4:
                log.debug("Found a short Obs. object sequence. Only %d frames found. Required >4 frames",len(obj_seq))
                raise Exception("Found a short Obs. object sequence. Only %d frames found. Required >4 frames",len(obj_seq))
                #return []
            else:
                #avoid call getCalibFor() when red_mode="single"
                if red_mode=="single":
                    dark,flat,bpm=[],[],[]
                else:
                    dark, flat, bpm = self.getCalibFor(obj_seq)
                    # return 3 list of calibration frames (dark, flat, bpm), because there might be more than one master dark/flat/bpm
                obj_ext, next = self.split(obj_seq) # it must return a list of list (one per each extension)
                dark_ext, cext = self.split(dark)
                flat_ext, cext = self.split(flat)
                bpm_ext, cext = self.split(bpm)
                for n in range(next):
                    log.debug("===> Reducting extension group %d", n+1)
                    ## At the moment, we have the first calibration file for each extension; what rule could we follow ?
                    if dark_ext==[]: mdark=None
                    else: mdark=dark_ext[n][0]  # At the moment, we have the first calibration file for each extension
                    if flat_ext==[]: mflat=None
                    else: mflat=flat_ext[n][0]  # At the moment, we have the first calibration file for each extension
                    if bpm_ext==[]: mbpm=None
                    else: mbpm=bpm_ext[n][0]   # At the moment, we have the first calibration file for each extension
                    try:
                        out_ext.append(self.reduceObj(obj_ext[n], mdark, mflat, mbpm, red_mode, \
                                                      output_file=self.out_dir+"/out_Q%02d.fits"%(n+1)))
                    except Exception,e:
                        log.error("Some error while reduction of extension group %d of object sequence %d", n+1, i)
                        raise e
        
            # if all reduction were fine, now join/stich back the extensions in a widther frame
            seq_result_outfile=self.out_file.replace(".fits","_SEQ%02d.fits"%(i))
            if len(out_ext)>1:
                log.debug("*** Creating final output file joining/stiching single output frames....***")
                #mef=misc.mef.MEF(outs)
                #mef.createMEF(self.out_file)
                # other option, do a SWARP to register the N-extension into one wide-single extension
                log.debug("*** Coadding overlapped files....")
                swarp = astromatic.SWARP()
                swarp.config['CONFIG_FILE']="/disk-a/caha/panic/DEVELOP/PIPELINE/PANIC/trunk/config_files/swarp.conf"
                swarp.ext_config['COPY_KEYWORDS']='OBJECT,INSTRUME,TELESCOPE,IMAGETYP,FILTER,FILTER2,SCALE,MJD-OBS'
                swarp.ext_config['IMAGEOUT_NAME']= seq_result_outfile
                swarp.ext_config['WEIGHTOUT_NAME']=self.out_file.replace(".fits",".weight.fits")
                swarp.ext_config['WEIGHT_TYPE']='MAP_WEIGHT'
                swarp.ext_config['WEIGHT_SUFFIX']='.weight.fits'
                swarp.run(out_ext, updateconfig=False, clean=False)
            else:
                shutil.move(out_ext[0], seq_result_outfile)
        
            i+=1
            seq_outfile_list.append(seq_result_outfile)
            log.info("*** Obs. Sequence reduced. File %s created.  ***", seq_result_outfile)
                    
            
        #log.critical(" TBD : for each object sequence, we need a different final output file !!!!")
        log.critical("*** All Obs.Sequences (%d) has been reduced. Congratulations !!", i)
        
        for f in seq_outfile_list: self.db.insert(f)
        self.db.ListDataSet()
        
        return seq_outfile_list
        
    def reduceObj(self, obj_frames, master_dark, master_flat, master_bpm, red_mode, output_file):
        """ Given a set of object(sci) frames and (optional) master calibration files,
            run the data reduction of the observing object sequence, producing an reduced
            ouput frame if no error; otherwise return None or raise exception.
            
            NOTE: Currently this method only accepts single FITS files (not MEF), it means
                  the splitting must be done previusly to call this method.
        """
        
        log.info("##################################")
        log.info("#### Starting data reduction #####")
        log.info("#### MODE = %s  ##", self.red_mode)
        log.info("##################################")
        #print "OBJS =",obj_frames
        
        # set the reduction mode
        if red_mode!=None: self.red_mode=red_mode
        
        dark_flat = False
        
        # Clean old files 
        #self.cleanUpFiles()
        
        # Change cwd to self.out_dir
        old_cwd=os.getcwd()
        os.chdir(self.out_dir) 
        
        # Copy/link source files (file or directory) to reduce to the working directory
        if not os.path.dirname(obj_frames[0])==self.out_dir:
            misc.fileUtils.linkSourceFiles(obj_frames, self.out_dir)
            #files1=[line.replace( "\n", "") for line in fileinput.input(self.list_file)]
            self.m_LAST_FILES=[self.out_dir+"/"+os.path.basename(file_i) for file_i in obj_frames]
        else:
            print "Input files already in output directory!"
            self.m_LAST_FILES=obj_frames
            
        print "SOURCES=\n",self.m_LAST_FILES
        
        ######################################################
        # 0 - Some checks (filter, ....) 
        ######################################################
        # TODO : it could/should be done in reduceSet, to avoid the spliting ...??
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
        if master_dark!=None and master_flat!=None:
            log.info("**** Applying dark and Flat ****")
            dark_flat=True
            res = reduce.ApplyDarkFlat(self.m_LAST_FILES, master_dark, master_flat, self.out_dir)
            self.m_LAST_FILES = res.apply()
        
        ######################################    
        # 2 - Compute Super Sky Flat-Field 
        ######################################    
        if master_flat==None:
            # - Find out what kind of observing mode we have (dither, ext_dither, ...)
            log.info('**** Computing Super-Sky Flat-Field ****')
            master_flat=self.out_dir+"/superFlat.fits"
            if self.obs_mode=="dither":
                log.debug("---> dither sequece <----")
                misc.utils.listToFile(self.m_LAST_FILES, self.out_dir+"/files.list")
                superflat = reduce.SuperSkyFlat(self.out_dir+"/files.list", master_flat, bpm=None, norm=False)
                superflat.create()
            elif self.obs_mode=="dither_on_off" or self.obs_mode=="dither_off_on" or self.obs_mode=="other":
                log.debug("----> EXTENDED SOURCE !!! <----")
                sky_list=self.getSkyFrames()
                misc.utils.listToFile(sky_list, self.out_dir+"/files.list")
                superflat = reduce.SuperSkyFlat(self.out_dir+"/files.list", master_flat, bpm=None, norm=False)
                superflat.create()                            
            else:
                log.error("Dither mode not supported")
                raise Exception("Error, dither mode not supported") 
        else:
            log.info("Using the given (dome or twlight) master flat")
              
        ######################################    
        # 3 - Compute Gain map and apply BPM
        ######################################
        log.info("**** Computing gain-map ****")
        gainmap = self.out_dir+'/gain_'+self.m_filter+'.fits'
        #nxblock=16
        #nyblock=16
        #The next values are to find out bad pixels 
        #nsig=5
        #mingain=0.7
        #maxgain=1.3
        g=reduce.calGainMap.GainMap(master_flat, gainmap)
        g.create() 
           
        ########################################
        # Add external Bad Pixel Map to gainmap
        ########################################     
        if master_bpm !=None:
            if not os.path.exists( master_bpm ):
                print('No external Bad Pixel Mask found. Cannot find file : "%s"' %master_bpm)
            else:
                iraf.imarith(operand1=gainmap,
                  operand2=master_bpm,
                  op='*',
                  result=gainmap.replace(".fits","_bpm.fits"),
                  verbose='yes'
                  )
                #replace old gain file
                os.rename(gainmap.replace(".fits","_bpm.fits"), gainmap)
                        
        #########################################
        # 4 - First Sky subtraction (IRDR) - sliding window technique
        #########################################
        log.info("**** 1st Sky subtraction (without object mask) ****")
        misc.utils.listToFile(self.m_LAST_FILES, self.out_dir+"/skylist1.list")
        # return only filtered images; in exteded-sources, sky frames  are not included 
        self.m_LAST_FILES=self.skyFilter(self.out_dir+"/skylist1.list", gainmap, 'nomask', self.obs_mode)      
                        
        #########################################
        # 5 - Quality assessment (FWHM, background, ellipticity, PSF quality)  
        #########################################
                            
        log.info("**** Quality Assessment: TODO ****")                   

        ## -- una prueba con astrowarp ---
        """
        if self.obs_mode!='dither' or self.red_mode=="single":
            log.info("**** Doing Astrometric calibration and  coaddition result frame ****")
            #misc.utils.listToFile(self.m_LAST_FILES, self.out_dir+"/files_skysub.list")
            aw = reduce.astrowarp.AstroWarp(self.m_LAST_FILES, catalog="2MASS", coadded_file=output_file)
            try:
                aw.run()
            except Exception,e:
                log.error("Some error while running Astrowarp....")
                raise e
            log.info("Generated output file ==>%s", output_file)
            log.info("#########################################")
            log.info("##### End of SINGLE data reduction ######")
            log.info("#########################################")
            return output_file
        """
        #########################################
        # 6 - Compute dither offsets from the first sky subtracted/filtered images using cross-correlation
        #########################################
        misc.utils.listToFile(self.m_LAST_FILES, self.out_dir+"/files_skysub.list")
        offset_mat=self.getPointingOffsets (self.out_dir+"/files_skysub.list", 15, 5, self.out_dir+'/offsets1.pap')                
                        
        
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
        self.coaddStackImages(self.out_dir+'/stack1.pap', gainmap, self.out_dir+'/coadd1.fits','average')
    
        ## END OF SINGLE REDUCTION  ##
        if self.obs_mode!='dither' or self.red_mode=="single":
            log.info("**** Doing Astrometric calibration of coadded result frame ****")
            reduce.astrowarp.doAstrometry(self.out_dir+'/coadd1.fits', output_file, "2MASS" ) 
            log.info("Generated output file ==>%s", output_file)
            log.info("#########################################")
            log.info("##### End of SINGLE data reduction ######")
            log.info("#########################################")
            return output_file 
        
        log.info("************************")
        log.info(" START SECOND PASS      ")
        log.info("************************")
        
        #########################################
        # 8 - Create master object mask
        #########################################
        log.info("**** Master object mask creation ****")
        obj_mask=self.__createMasterObjMask(self.out_dir+'/coadd1.fits', self.out_dir+'/masterObjMask.fits')  
        
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
        self.m_LAST_FILES=self.skyFilter( self.out_dir+"/skylist2.pap", gainmap, 'mask', self.obs_mode)      
    
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
            except Exception,e:
                log.error("Some error while running Astrowarp....")
                raise e
        
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
        self.coaddStackImages(self.out_dir+'/stack1.pap', gainmap, self.out_dir+'/coadd2.fits')
        reduce.imtrim.imgTrim(self.out_dir+'/coadd2.fits')
        
        #########################################
        # 10 - Make Astrometry
        #########################################
        log.info("**** Computing Astrometric calibration of coadded(2) result frame ****")
        reduce.astrowarp.doAstrometry(self.out_dir+'/coadd2.fits', output_file, "2MASS" )
        
        log.info("Generated output file ==>%s", output_file)
        
        os.chdir(old_cwd)
        
        log.info("##################################")
        log.info("##### End of data reduction ######")
        log.info("##################################")
        
        return output_file 
        
    def test(self):
        log.debug("Una prueba")
        self.__initDB()
        self.db.GetOBFiles()
        (seq_par, seq_list)=self.db.GetSeqFiles(filter=None,type='SCIENCE')
        k=0
        for par in seq_par:
            print "\nSEQUENCE PARAMETERS - OB_ID=%s,  OB_PAT=%s, FILTER=%s"%(par[0],par[1],par[2])
            print "-----------------------------------------------------------------------------------\n"
            for file in seq_list[k]:
                print file
            k+=1
        
################################################################################
# main
################################################################################
if __name__ == "__main__":
    log.debug("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    log.debug(">Start PAPI....")
    log.debug("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    # Get and check command-line options
    
    
    usage = "usage: %prog [options] arg1 arg2 ..."
    parser = OptionParser(usage)
    
                  
    parser.add_option("-s", "--source",
                  action="store", dest="source",
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

    parser.add_option("-n", "--no_class",
                  action="store_true", dest="no_class", help="not try to do dataset classification", default=False)
    
    parser.add_option("-v", "--verbose",
                  action="store_true", dest="verbose", default=True,
                  help="verbose mode [default]")

                                
    (options, args) = parser.parse_args()
    
    if not options.source or not options.output_filename or len(args)!=0: # args is the leftover positional arguments after all options have been processed
        parser.print_help()
        parser.error("incorrect number of arguments " )
    if options.verbose:
        print "reading %s ..." % options.source
    
    
    sci_files=[]
    
    # Read file or list-directory 
    if os.path.isfile(options.source):
        sci_files=[line.replace("\n", "").replace('//','/') for line in fileinput.input(options.source)]
    elif os.path.isdir(options.source):
        for file in dircache.listdir(options.source):
            if file.endswith(".fits") or file.endswith(".fit"):
                sci_files.append((options.source+"/"+file).replace('//','/'))
                
    if options.single: red_m="single"
    else: red_m="full"
    
    # Create the RS (it works both simple FITS as MEF files)            
    rs = ReductionSet( sci_files, options.out_dir, out_file=options.output_filename, obs_mode="dither", \
                                dark=options.dark, flat=options.flat, bpm=options.bpm, red_mode=red_m, classf=not(options.no_class))
    
    """try:
        rs.test()
    except Exception,e:
        print "some error !!!"
        raise e
    
    #sys.exit(0)
    """
    
    try:
        rs.reduceSet(red_mode=red_m)
    except Exception,e:
        print "Cannot reduce the Data Set, check error log...."
        raise e
    else:
        print "Yuppiiii !!!"
        
    """
    try:
        out = rs.buildMasterDarks()
    except:
        print "Cannot build Master Darks !"
    
    try:
        out2 = rs.buildMasterTwFlats()
    except:
        print "Cannot build TwFlats !"
        
    try:
        out3 = rs.buildMasterDomeFlats()
    except:
        print "Cannot build Master Dome Flats"
    """
    
    if options.show:
        out.show()    
