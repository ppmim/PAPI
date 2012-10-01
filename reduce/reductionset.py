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
# PAPI (PAnic PIpeline)
#
# reductionset.py
#
# Last update 30/Nov/2010
#
################################################################################
    
#From system
#import sys
import os
import os.path
import fileinput
#import glob
import shutil
import tempfile
import dircache
#import math
#import pprocess # to be removed  !!
import multiprocessing
import itertools

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
import reduce.remove_cosmics
import reduce.astrowarp
import misc.mef 
import astromatic
from astromatic.swarp import *
import datahandler.dataset
import misc.collapse

#
# Next functions are needed to allow the use of  multiprocessing.Pool() with 
# class methods, that need to be picklable (at least 
# Solution obtained from :
# http://www.frozentux.net/2010/05/python-multiprocessing/  
# http://stackoverflow.com/questions/3288595/multiprocessing-using-pool-map-on-a-function-defined-in-a-class
# In addition, there is other version of the above functions:
# http://stackoverflow.com/questions/5429584/handling-the-classmethod-pickling-issue-with-copy-reg
# but it does not work (at least for me!) 
#

def _pickle_method(method):
    """
    Pickle methods properly, including class methods.
    """
    func_name = method.im_func.__name__
    obj = method.im_self
    cls = method.im_class
    return _unpickle_method, (func_name, obj, cls)

def _unpickle_method(func_name, obj, cls):
    """
    Unpickle methods properly, including class methods.
    """
    for cls in cls.mro():
        try:
            func = cls.__dict__[func_name]
        except KeyError:
            pass
        else:
            break
    return func.__get__(obj, cls)        

import copy_reg 
import types 

copy_reg.pickle(types.MethodType,  
    _pickle_method,  
    _unpickle_method)  



class ReductionSetException(Exception):
    pass


class ReductionSet(object):
    """
    This class implement a set of FITS files from an PANIC observation. There can
    be calibration or/and science files of any kind, althougt is expected they are
    compliant with a well-defined observation run with the PANIC OT.
    
    Then a Reduction Set (RS) can have:
        - dark sequences
        - flat sequences
        - focus sequences
        - science sequences
        - 
        
    """
     
    def __init__(self, rs_filelist, out_dir=None, out_file=None, obs_mode="dither", 
                 dark=None, flat=None, bpm=None, red_mode="quick", 
                 group_by="ot", check_data=True, config_dict=None, 
                 external_db_files=None, temp_dir = None,
                 *a, **k):
        """ 
        Init object 
        
        Parameters
        ----------
        rs_filelist : list
            list of files of the reduction set
        
        out_dir : str
            output dir where output files will be generated
        
        obs_mode : str
            'dither' - files belong to a dither pattern observation
            'other' - 
        
        dark : str
            if given, the filename of the master dark to be used for the reduction
        
        flat : str
            if given, the filename of the master flat to be used for the reduction
        
        external_db_files : str or list
            File list or Directory used as an external calibration database.
            Then, if during the reduction of a ReductionSet(RS) no calibration (dark, flat) 
            are found in the current RS, then PAPI will look for them into this directory.
            If the directory does not exists, of no calibration are found, then no calibrations
            will be used for the data reduction.
            Note that the calibrations into the current RS have always higher priority than
            the ones in the external calibration DB.   
        
        ....
        
        """
        
        super (ReductionSet, self).__init__ (*a, **k)
        
        # CONFIG dictionary
        if not config_dict:
            raise Exception("Config dictionary not provided ...")
        else:
            # dictionary with all sections (general, darks, dflats, 
            # twflats, skysub, fits, keywords, config_dicts) and their values
            self.config_dict = config_dict 

        # Input files
        if len(rs_filelist) <= 0:
            log.error("Empy file list, no files to reduce ...")
            raise ReductionSetException ("Empy file list, no files to reduce ...")
        else:
            for f in rs_filelist:
                if not os.path.exists(f) or not os.path.splitext(f)[1] == ".fits":
                    log.error("File %s does not exists or not has '.fits' extension", f)
                    raise ReductionSetException ("File %s does not exist or \
                    has not '.fits' extension"%f)
                    
        # Main directories
        self.rs_filelist = rs_filelist # list containing the science data filenames to reduce
        # Temporal directory for temporal files created during the data reduction process
        if temp_dir == None:
            if self.config_dict !=None:
                self.temp_dir = self.config_dict['general']['temp_dir']
            else:
                self.temp_dir = "/tmp/"
        else:
            self.temp_dir = temp_dir

        # Output directory for the results of the data reduction process            
        if out_dir == None:
            if self.config_dict !=None:
                self.out_dir = self.config_dict['general']['output_dir']
            else:
                self.out_dir = "/tmp/"
        else:
            self.out_dir = out_dir
            
        
        if not os.path.exists(self.out_dir):
            raise Exception("Output dir does not exist.")
        
        # Main output file resulted from the data reduction process
        if out_file==None: # if no filename, we choose a random filename
            output_fd, self.out_file = tempfile.mkstemp(suffix='.fits', prefix='red_set_', dir=self.out_dir)
            os.close(output_fd)
            os.unlink(self.out_file) # we only need the name
        else:    
            self.out_file = out_file  # final reduced data file (out)
            
        
        self.obs_mode = obs_mode  # observing mode (dither, dither_on_off, dither_off_on....)
        self.master_dark = dark # master dark to use (input)
        self.master_flat = flat # master flat to use (input)
        self.master_bpm = bpm # master Bad Pixel Mask to use (input)
        self.apply_dark_flat = self.config_dict['general']['apply_dark_flat'] # 0=no, 1=before, 2=after
        self.red_mode = red_mode # reduction mode (lemon=for LEMON pipeline, 
                                 # quick=for QL, science=for science) 
        log.debug("GROUP_BY = %s"%group_by)
        self.group_by = group_by.lower() # flag to decide how classification 
                                         # will be done, by OT (OB_ID, OB_PAT, 
                                         # FILTER and TEXP kws) or by 
                                         # FILTER (FILTER kw)
        self.check_data = check_data # flat to indicate if data checking need 
                                     # to be done (see checkData() method)
         
        
        if self.config_dict:
            self.m_terapix_path = self.config_dict['config_files']['terapix_bin']
            self.m_irdr_path = self.config_dict['config_files']['irdr_bin']
            # update config values from config_dict
            self.HWIDTH = self.config_dict['skysub']['hwidth'] # half width of sky filter window in frames
            self.MIN_SKY_FRAMES = self.config_dict['skysub']['min_frames']  # minimun number of sky frames required in the sliding window for the sky subtraction
            self.MAX_MJD_DIFF = self.config_dict['general']['max_mjd_diff'] # Maximun exposition time difference (seconds) between frames
            self.MIN_CORR_FRAC = self.config_dict['offsets']['min_corr_frac'] # Minimun overlap correlation fraction between offset translated images (from irdr::offset.c)
        else:
            print "Program should not enter here  !!!"
            # some "default" config values (see below how they are updated from the config_dict)
            # Environment variables
            """
            self.m_terapix_path = os.environ['TERAPIX']
            self.m_irdr_path = os.environ['IRDR_BIN']
            self.MAX_MJD_DIFF = (1/86400.0)*10*60 #6.95e-3  # Maximun seconds (10min=600secs aprox) of temporal distant allowed between two consecutive frames 
            self.HWIDTH = 2 #half width of sky filter window in frames
            self.MIN_SKY_FRAMES = 5  # minimun number of sky frames required in the sliding window for the sky subtraction
            self.MIN_CORR_FRAC = 0.1 # Minimun overlap correlation fraction between offset translated images (from irdr::offset.c)
            """
        
        # real "variables" holding current reduction status
        self.m_LAST_FILES = []   # Contain the files as result of the last 
                                 # processing step (science processed frames). 
                                 # Properly initialized in reduceSingleObj()
        self.m_rawFiles = []     # Raw files (originals in the working directory)
                                 # of the current sequence being reduced.
        self.m_filter = ""       # Filter of the current data set (m_LAST_FILES)
        self.m_type = ""         # Type (dark, flat, object, ...) of the current data set; should be always object !
        self.m_expt = 0.0        # Exposition Time of the current data set files
        self.m_ncoadd = 0        # Number of coadds of the current data set files
        self.m_itime  = 0.0      # Integration Time of the current data set files
        self.m_readmode = ""     # readout mode (must be the same for all data files)
        
        
        ##Local DataBase (in memory)
        self.db = None
        
        ##External DataBase (in memory) - optional
        # It is an optional external database provided when creating a ReductionSet 
        # instance and mainly can have calibration files required for the reduction 
        # of the current data set. It will be used during the on-line QL 
        # data reduction,e.g., tw_flats that require a master dark than can be 
        # in other RS. Mainly used in Quick-Look !!    
        self.ext_db = None
        
        # (optional) file list to build the external DB. We proceed this way, because
        # if we give a DB connection to the ReductionSet class instead of a list of files,
        # we can have problems because SQLite3 does not support access from multiple 
        # threads, and the RS.reduceSet() can be executed from other thread than 
        # it was created.
        # --> If can be a list of files or a directory name having the files
        # Further info: http://stackoverflow.com/questions/393554/python-sqlite3-and-concurrency  
        #               http://www.sqlite.org/cvstrac/wiki?p=MultiThreading
        self.ext_db_files = []
        if external_db_files == None:
            if self.config_dict !=None:
                cal_dir = self.config_dict['general']['ext_calibration_db']
                if os.path.isdir(cal_dir):
                    for file in dircache.listdir(cal_dir):
                        if file.endswith(".fits") or file.endswith(".fit"):
                            self.ext_db_files.append((cal_dir+"/"+file).replace('//','/'))
        else:
            self.ext_db_files = external_db_files    

        # Print config_dictonary values in log file for debugging
        log.info("[ReductionSet] CONFIGURATION VALUES OF RS ------------------")
        log.info(printDict(self.config_dict))
        log.info("------------------------------------------------------------")
  
    def __initDB(self):
        """
        @summary: Initialize the Data Bases (local and external), loading the full frame list.
        
        @attention: 
        If we initialize the DB in the constructor __init__(), we will have problems
        if we use the DB from any method of ReductionSet and then from other thread,
        i.e.,sqlite3 do not support access from multiple threads.
        For more info, see http://stackoverflow.com/questions/393554/python-sqlite3-and-concurrency
                           http://www.sqlite.org/cvstrac/wiki?p=MultiThreading
        """
        #Local DataBase (in memory)
        log.info("Initializing Local DataBase")
        try:
            self.db = datahandler.dataset.DataSet(self.rs_filelist)
            self.db.createDB()
            self.db.load()
        except Exception,e:
            log.error("Error while LOCAL data base initialization: \n %s"%str(e))
            raise Exception("Error while LOCAL data base initialization")
            
        
        self.m_LAST_FILES = self.rs_filelist # Later properly initialized in reduceSingleObj()
        
        if self.master_dark!=None: self.db.insert(self.master_dark)
        if self.master_flat!=None: self.db.insert(self.master_flat)
        if self.master_bpm !=None: self.db.insert(self.master_bpm)
        
        self.db.ListDataSet()
        
        #External DataBase (in memory)
        if len(self.ext_db_files)>0:
            log.info("Initializing External DataBase")
            try:
                self.ext_db = datahandler.dataset.DataSet(self.ext_db_files)
                self.ext_db.createDB()
                self.ext_db.load()
            except Exception,e:
                log.error("Error while EXTERNAL data base initialization: \n %s"%str(e))
                raise Exception("Error while EXTERNAL data base initialization")

            self.ext_db.ListDataSet()
        else:
            self.ext_db = None
                
        
    def checkData(self, chk_shape=True, chk_filter=True, chk_type=True, chk_expt=True,
                  chk_itime=True, chk_ncoadd=True, chk_cont=True,
                  chk_readmode=True, file_list=None):
        """
        Check if data properties match in the current file list.
        
        Parameters
        ---------- 
        
        check_shape: bool
            Flag to indicate if shape (naxis1, naxis2, naxis3) must be checked
            
        chk_filter: bool
            Flag to indicate if filter property must be checked
        
        chk_type: bool
            Flag to indicate if type (dark, flat, science, ...) must be checked
        
        chk_expt: bool
            Flag to indicate if type Exposition Time must be checked

        chk_itime: bool
            Flag to indicate if type Integration Time (ITIME) must be checked

        chk_ncoadd: bool
            Flag to indicate if type number of coadds (NCOADD) must be checked

        chk_cont: bool
            Flag to indicate if type temporal continuity (MJD) must be checked

        chk_readmode: bool
            Flag to indicate if type Readout mode  must be checked

        file_list : list
            A list having the FITS files to check
            
        Return: 
        -------
            - True whether all files in 'file_list' have the same properties 
            - False otherwise.
            
        Note:
        -----
        Also is checked the temporal continuity (distant between two 
        consecutive frames),if exceeded raise an exception 
        """
        
        if file_list and len(file_list)>0:
            files_to_check = file_list
        else:
            files_to_check = self.m_LAST_FILES
            
        f = datahandler.ClFits(files_to_check[0])
        
        filter_0 = f.getFilter()
        type_0 = f.getType()
        expt_0 = f.expTime()
        itime_0 = f.getItime()
        ncoadd_0 = f.getNcoadds()
        readmode_0 = f.getReadMode()
        shape_0 = f.shape
        
        self.m_filter = filter_0
        self.m_type = type_0
        self.m_expt = expt_0
        self.m_itime = itime_0
        self.m_ncoadd = ncoadd_0
        self.m_readmode = readmode_0
        
        
        mismatch_filter = False
        mismatch_type = False
        mismatch_expt = False
        mismatch_itime = False
        mismatch_ncoadd = False
        mismatch_cont = False
        mismatch_readmode = False
        mismatch_shape = False
        
        prev_MJD = -1
        for file in files_to_check:
            fi = datahandler.ClFits( file )
            if chk_shape and not mismatch_shape: 
                if fi.shape != shape_0:
                    log.debug("File %s does not match SHAPE", file)
                    mismatch_shape = True
                    break
            if chk_filter and not mismatch_filter: 
                if fi.getFilter() != filter_0:
                    log.debug("File %s does not match FILTER", file)
                    mismatch_filter = True
                    break
            if chk_type and not mismatch_type: 
                if fi.getType() != type_0:
                    log.debug("File %s does not match TYPE", file)
                    mismatch_type = True
                    break
            if chk_expt and not mismatch_expt: 
                if fi.expTime() != expt_0:
                #if prev_MJD!=-1 and ((fi.expTime()+self.MAX_MJD_DIFF)<expt_0 or
                #    (fi.expTime()-self.MAX_MJD_DIFF)>expt_0):   # more relaxed situation
                    log.debug("File %s does not match EXPTIME", file)
                    mismatch_expt = True
                    break
            if chk_itime and not mismatch_itime: 
                if fi.getItime() != itime_0:
                    log.debug("File %s does not match ITIME", file)
                    mismatch_itime = True
                    break
            if chk_ncoadd and not mismatch_ncoadd: 
                if fi.getNcoadds() != ncoadd_0:
                    log.debug("File %s does not match NCOADD", file)
                    mismatch_ncoadd = True
                    break
            if chk_readmode and not mismatch_readmode: 
                if fi.getReadMode() != readmode_0:
                    log.debug("File %s does not match READMODE", file)
                    mismatch_readmode = True
                    break
            if chk_cont and not mismatch_cont:
                if prev_MJD!=-1 and (fi.getMJD()-prev_MJD)>self.MAX_MJD_DIFF:
                    log.error("Maximmun time distant between two consecutives frames exceeded. File= %s",file)
                    mismatch_cont = True
                    break
                else:
                    prev_MJD=fi.getMJD()
                    
        if mismatch_shape or mismatch_filter or mismatch_type or mismatch_expt or  \
            mismatch_itime or mismatch_ncoadd or mismatch_cont:
            log.error("Data checking found a mismatch....check your data files....")
            #raise Exception("Error while checking data (filter, type, ExpT, Itime, NCOADDs, MJD)")
            return False             
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
        Ascending sort of input data files by MJD
        """
        
        dataset = []
        if list==None:
            m_list = self.m_LAST_FILES
        else:
            m_list = list
            
        for file in m_list:
            fits = datahandler.ClFits(file) 
            dataset.append((file, fits.getMJD()))
            
        dataset = sorted(dataset, key=lambda data_file: data_file[1])          
        sorted_files = []
        for tuple in dataset:
            sorted_files.append(tuple[0])
        
        return sorted_files

    def split(self, frame_list):
        """ 
        Split the data from the given frame list (any kind, science or 
        calibration) into N 'sub-list', where N is the number of extension of 
        the Multi-Extension FITS.
        
        Return a list with the split frames and the number of extensions.
        """
        
        log.debug("Starting split() method ....")
        
        new_frame_list=[] # a list of N list, where N=number of extension of the MEF 
        nExt=0
        if frame_list==None or len(frame_list)==0 or frame_list[0]==None:
            return [],0
            
        # First, we need to check if we have MEF files
        if not datahandler.ClFits( frame_list[0] ).isMEF() and \
            not datahandler.ClFits( frame_list[0] ).isFromGEIRS():
            nExt = 1
            new_frame_list.append(frame_list)
        else:            
            #Suppose we have MEF files ...
            if datahandler.ClFits(frame_list[0]).getInstrument()=="HAWKI":
                
                kws_to_cp = ['DATE','OBJECT','DATE-OBS','RA','DEC','EQUINOX','RADECSYS','UTC','LST',\
		  'UT','ST','AIRMASS','IMAGETYP','EXPTIME','TELESCOP','INSTRUME','MJD-OBS',\
		  'FILTER', 'FILTER1', 'FILTER2', "HIERARCH ESO TPL ID", "HIERARCH ESO TPL EXPNO", "HIERARCH ESO TPL NEXP"\
            ]   
            else: # PANIC
                kws_to_cp = ['DATE','OBJECT','DATE-OBS','RA','DEC','EQUINOX','RADECSYS','UTC','LST',\
          'UT','ST','AIRMASS','IMAGETYP','EXPTIME','TELESCOP','INSTRUME','MJD-OBS',\
          'FILTER', 'FILTER1', 'FILTER2', "OBS_TOOL", "OB_ID", "OB_PAT",\
          "PAT_EXPN", "PAT_NEXP", "END_SEQ"\
            ]
            
            if datahandler.ClFits(frame_list[0]).isFromGEIRS():
                try:
                    mef = misc.mef.MEF(frame_list)
                    (nExt, sp_frame_list) = mef.splitGEIRSToSimple(".Q%02d.fits", 
                                                                   out_dir=self.out_dir)
                except Exception,e:
                    log.debug("Some error while splitting GEIRS data set. %s",str(e))
                    raise
            else:
                log.debug("Splitting data files")
                try:
                    mef = misc.mef.MEF(frame_list)
                    (nExt, sp_frame_list) = mef.doSplit(".Q%02d.fits", 
                                                        out_dir=self.out_dir, 
                                                        copy_keyword=kws_to_cp)
                except Exception,e:
                    log.debug("Some error while splitting data set. %s",str(e))
                    raise
            
            # now, generate the new output filenames        
            for n in range(1,nExt+1):
                new_frame_list.append([self.out_dir+"/"+os.path.basename(file.replace(".fits",".Q%02d.fits"%n)) for file in frame_list])
                """
                for f in new_file_names:
                    #if re.search(".*(\.Q01)(.fits)$", f):
                        sources.append(f)
                """
        
        #print "NEW_FRAME_LIST=",new_frame_list
        
        return new_frame_list, nExt
    
    def getObsMode(self):
        """
        Return the type of dither sequece followed in the currect  'm_LAST_FILES' list. It could be:
            - dither (T-T-T-T-T- ....)
            - dither_on_off (T-S-T-S-T-S-T-S-T-....)
            - dither_off_on (S-T-S-T-S-T-S-T-S-....)
            - other  (non defined sequence,unknown) (e.g,  T-S-T-T-S-T-T-S-T-....)
            
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
            # check if all are objects ...
            # and if are, then, we are going to suppose the sequence T-T-T-T-T- .... (dither)
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
        @summary: Given a list of frames belonging to a observing sequence for 
        a given object (star, galaxy, whatever),return the most recently created 
        calibration files (master dark,flat,bpm) in order to reduce the sequence.
        The search of the calibration files is done, firstly in the local DB, but
        if no results, then in the external DB if it was provided.
          
        @return:  3 calibration files (dark, flat, bpm); If more than one master
        were found, the most recently created (according to MJD) is returned.
        If some master were not found, None is returned.
    
        """
        
        log.debug("Looking for calibration files into DB")
        
        master_dark = [] # we'll get a list of master dark candidates
        master_flat = [] # we'll get a list of master flat candidates
        master_bpm = [] # we'll get a list of master flat candidates
        
        obj_frame = datahandler.ClFits(sci_obj_list[0])
        # We take as sample, the first frame in the list, but all frames must
        # have the same features (expT,filter,ncoadd, readout-mode, ...)
        expTime = obj_frame.expTime()
        filter = obj_frame.getFilter()
        
        #DARK - Do NOT require equal EXPTIME Master Dark ???
        master_dark = self.db.GetFilesT('MASTER_DARK_MODEL', -1) 
        if len(master_dark)==0 and self.ext_db!=None:
            master_dark = self.ext_db.GetFilesT('MASTER_DARK_MODEL', -1) 
        #FLATS - Do NOT require equal EXPTIME, but FILTER
        master_flat = self.db.GetFilesT('MASTER_DOME_FLAT', -1, filter)
        if master_flat==[]:
            master_flat = self.db.GetFilesT('MASTER_TW_FLAT', -1, filter)
        if len(master_flat)==0 and self.ext_db!=None:
            master_flat = self.ext_db.GetFilesT('MASTER_DOME_FLAT', -1, filter)
            if len(master_flat)==0:
                master_flat=self.ext_db.GetFilesT('MASTER_TW_FLAT', -1, filter)

        #BPM                
        master_bpm = self.db.GetFilesT('MASTER_BPM')
        if len(master_bpm)==0 and self.ext_db!=None:
            master_bpm = self.ext_db.GetFilesT('MASTER_BPM')

        log.debug("Master Darks found %s", master_dark)
        log.debug("Master Flats found %s", master_flat)
        log.debug("Master BPMs  found %s", master_bpm)
        
        # Return the most recently created (according to MJD order)
        if len(master_dark)>0: r_dark = master_dark[-1]
        else: r_dark = None
        if len(master_flat)>0: r_flat = master_flat[-1]
        else: r_flat = None
        if len(master_bpm)>0: r_bpm = master_bpm[-1]
        else: r_bpm = None
        
        return r_dark, r_flat, r_bpm
        
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
        match_list = sorted(match_list, key=lambda data_file: data_file[1]) 
        
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
    
    def getSequences(self, show=True):
        """
        Look for sequences (calib,science) in the current data set following the
        grouping criteria given by self.group_by and orderted by MJD
        
        @param show: if True, print out in std output the found sequences
        
        @return: a list of lists of sequence files and their types (DARK, TW_FLAT, 
                DOME_FLAT, SCIENCE) ; see ClFits class for further details)
        
                In case of group_by, only 'SCIENCE' type is returned
        
        @attention: if group_by FILTER, only science sequences will be found
        
        """
        
        log.debug("[getSequences] Looking for Data Sequences into the DataSet")
        seqs = []
        seq_types = []
        
        if self.group_by=='ot':
            seqs, seq_types = self.getOTSequences(show)
        #elif self.group_by=='filter':
        #    seqs = self.getObjectSequences() # only look for science sequences !
        #    seq_types = ['SCIENCE']*len(seqs)
        #elif self.group_by=='none':
        elif self.group_by=='filter':
            if self.db==None: self.__initDB()
            seqs, seq_types = self.db.GetSequences(group_by='filter')
        else:
            log.error("[getSequences] Wrong data grouping criteria")
            raise Exception("[getSequences] Found a not valid data grouping criteria %s"%(self.group_by))

        # Print out the groups
        if show:
            k=0
            print "\n ========================================================="
            print " =========== GROUPED SEQUENCES (by %s) =============="%self.group_by
            print " ========================================================="
            for type in seq_types:
                print "\nSEQUENCE #[%d]  - TYPE= %s   FILTER= %s  TEXP= %f  #files = %d " \
                        %(k, type, self.db.GetFileInfo(seqs[k][0])[3], 
                          self.db.GetFileInfo(seqs[k][0])[4], 
                          len(seqs[k]))
                print "-------------------------------------------------------------------"
                for file in seqs[k]:
                    print file + " type= %s"%self.db.GetFileInfo(file)[2]
                k+=1
            log.debug("Found %d groups of files", len(seq_types))
        

        #To print the sequences into the log file
        debug = 1
        if debug:    
            k=0
            log.debug("=========================================================")
            log.debug("=========== GROUPED SEQUENCES (by %s) =============="%self.group_by)
            log.debug("=========================================================")
            for type in seq_types:
                log.debug("SEQUENCE #[%d]  - TYPE= %s   FILTER= %s  TEXP= %f  #files = %d " \
                        %(k, type, self.db.GetFileInfo(seqs[k][0])[3], 
                          self.db.GetFileInfo(seqs[k][0])[4], 
                          len(seqs[k])))
                log.debug("-------------------------------------------------------------------")
                for file in seqs[k]:
                    log.debug("%s type = %s" %(file, self.db.GetFileInfo(file)[2]))
                k+=1
            
        return seqs, seq_types
    
    
    def getOTSequences(self, show=True):
        """
        Look for sequences (calib, science) in the current data set and print the results
        The sequences must follow the PANIC Observing Tool schema, and basically are 
        detected by the pair {PAT_EXPN, PAT_NEXP} :
        
        OBS_TOOL, OB_ID, OB_PAT, PAT_EXPN, PAT_NEXP, END_SEQ
         
        
        The algorith followed is the next:
        
            1. Sort out the data files by MJD
            2. Do
                2.1 Look for the first file with OBS_TOOL=True and PAT_EXPN=1
                    2.1.1 Group all next files up to PAT_EXPN=PAT_NEXP and/or END_SEQ=True
            3. while (files in DataSet)
            4. Print files of each Group found 
            
        @param show: if True, print out the sequences found in the std output
             
        @return: two lists:
                - a list of lists of sequence files belonging to
                - a list with the types of the sequences (DARK, TW_FLAT, 
                DOME_FLAT, SCIENCE) ; see ClFits class for further details)
        
        @attention: this method is an (better) alertenative to getObjectSequences()        
        """
        
        log.debug("[getOTSequences] Looking for OT generated Data Sequences into the DataSet")
        
        seq_list = [] # list of list of sequence filenames
        seq_types = [] # list of types of sequence (dark, dflat, sflat, focus, science, ....)
        
        if self.db==None: self.__initDB()
        seq_list, seq_types = self.db.GetSequences(group_by='ot') # much more quick than read again the FITS files
        
        return seq_list, seq_types
        
    
    def getObjectSequences(self):
        """
        Query the DB (data set) to look for (science) object/pointing  sequences 
        sorted by MJD; they can be grouped by:
        
            a) FILTER, TEXP (it is used when no data set classification can be done)
            b) OB_ID, OB_PAT, FILTER, TEXP (not sure to be good way to look for ObjSequences ???? 
            
        @return: a list of list of sequence files belonging to
        
        @attention: other alternative way of grouping following OT schema is
        implemented in getOTSequences()
        
        @note: Only SCIENCE sequences are looked for !
        @todo: Look for CALIBRATION sequences as well !
        """
        
        log.debug("[getObjectSequences] Looking for GEIRS/MIDAS generated Data Sequences into the DataSet")
        
        # Init the DB    
        if self.db==None: self.__initDB()
        
        # group data file (only science) by Filter,TExp 
        seq_list = []
        seq_par = []
        
        if self.group_by=="filter":
            (seq_list, seq_par) = self.db.GetSequences(group_by='filter') # return a list of SCIENCE (or SKY) frames grouped by FILTER
                
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
        match_list = sorted(match_list, key=lambda data_file: data_file[1]) 
        
        return match_list # a list of tuples as (file,filter)
    
    def isaCalibSet(self):
        """Check all the current files in the ReductionSet list to find out if they are all calibatrion frames"""
        
        """
        filelist = self.db.GetFilesT(type="SCIENCE", texp=-1, filter="ANY")
        filelist+= self.db.GetFilesT(type="SKY", texp=-1, filter="ANY")
        filelist+= self.db.GetFilesT(type="SKY", texp=-1, filter="ANY")
        
        if len(filelist)>0:
            return False
        else:
            return True
        
        """
        for file in self.rs_filelist:
            fi = datahandler.ClFits(file)
            if fi.isScience():
                log.debug("File %s does not match as calibration file type %s", file)
                return False
            
        log.debug("All files seem to be calibration files")
        return True
        
        
    def skyFilter( self, list_file, gain_file, mask='nomask', 
                   obs_mode='dither', skymodel=None ):
        """
        @summary: For 'each' (really not each, depend on dither pattern, e.g., 
        extended sources) input image, a sky frame is computed by combining a 
        certain number of the closest images, then this sky frame is subtracted
        to the image and the result is divided by the master flat; 
                         
        This function is a wrapper for skyfilter.c (IRDR)              
        
        @param list_file: a text file containing the suited structure in 
        function of the observing_mode (the list shall be sorted by obs-date)

        @param gain_file: gain map file used as bad pixel mask
        
        @param mask: [mask|nomask] flag used to indicate if a object mask was 
        especified into the 'list_file'
        
        @param obs_mode: [dither|other] dither or other(e.g., nodding) pattern 
        to process
          
        @param skymodel: [median|min] sky model used for sky subtraction. Only 
        required if obs_mode=dither. (median=coarse fields, min=crowded fields)
           
        @return: 
        The function generate a set of sky subtrated images (*.skysub.fits) and
        Return ONLY filtered images; when extended-source-mode ,sky 
        frames are not included in the returned file list. 
        The out-file-list is previously ordered by obs-data.             
    
        @version: 1.0, 20090909 by jmiguel@iaa.es
    
        @todo: extended objects !!
        """               
        
        # Skyfilter parameters
        halfnsky = self.HWIDTH  # value from config file [skysub.hwidth]
        destripe = 'none'
        out_files = []
        
        # get the skymodel
        if skymodel==None:
            skymodel = self.config_dict['skysub']['skymodel']
            
        if obs_mode=='dither':
            # It comprises any dithering pattern for point-like objects
            skyfilter_cmd = self.m_irdr_path + '/skyfilter '+ list_file + \
            '  ' + gain_file + ' ' + str(halfnsky) + ' ' + mask + '  ' +  \
            destripe + '  ' + skymodel
        elif obs_mode=='other':
            # It comprises nodding pattern for extended objects 
            skyfilter_cmd = self.m_irdr_path+'/skyfilter_general '+ list_file \
            + '  ' + gain_file +' '+ str(halfnsky)+' '+ mask + '  ' + destripe
        elif obs_mode=='dither_on_off': 
            # actually not used
            skyfilter_cmd = self.m_irdr_path + '/skyfilteronoff ' + list_file + \
            '  ' + gain_file + ' ' + str(halfnsky) + ' ' + mask + '  ' + destripe 
        elif obs_mode=='dither_off_on':
            # actually not used 
            skyfilter_cmd = self.m_irdr_path + '/skyfilteroffon ' + list_file +\
            '  ' + gain_file + ' ' + str(halfnsky) + ' ' + mask + '  ' + destripe
        else:
            log.error("Observing mode not supported")
            raise
                  
        if misc.utils.runCmd( skyfilter_cmd )==1: # All was OK
            
            # Rename output sky-subtracted files
            """for file in glob.glob(self.out_dir+'/*.fits.skysub'):
                shutil.move(file, file.replace('.fits.skysub', '.skysub.fits'))
                out_files.append(file.replace('.fits.skysub', '.skysub.fits'))
            """
            # look for sky subtracted images created by irdr::skyfilter            
            files = [line.split(" ")[0].replace("\n","") for line in fileinput.input(list_file)] # it takes into account the two kind of possible inputs files to skyfilter
            for file in files:
                if os.path.exists(file+".skysub"): # it takes into acount dither_on_off and other extended obs. patterns
                    shutil.move(file.replace(".fits", ".fits.skysub"), file.replace(".fits", ".skysub.fits"))
                    out_files.append(file.replace(".fits", ".skysub.fits"))
            
            ## Compose the output file list
            #if obs_mode=='dither':
            #   out_files=[line.split()[0].replace('.fits', '.skysub.fits') for line in fileinput.input(list_file)]
            #elif (obs_mode=='dither_on_off' or obs_mode=='dither_off_on' or obs_mode=='other'):
            #    out_files=glob.glob(self.out_dir+'/*.skysub.fits')
            ##
            
            # Sort-out data files by obs-data (actually not required when obs_mode='dither')
            out_files = self.sortOutData(out_files) 
            
            return out_files
        else:
            log.error("Some problem while running command %s", skyfilter_cmd) 
            return []
                                  
    
    def subtractNearSky(self, near_list=None, file_pos=0, out_filename=None):
        """
            Compute and subtract the nearest sky to the image in position 'fn' 
            in the given frame list (near_list). 
                         
            This function make use of skyfilter_single.c (IRDR)              
            
            INPUT
                file_pos : file position in sci file list (the list is supposed 
                to be sorted by obs-date)
            
            OUTPUT
                The function generate a sky subtrated image (*.skysub.fits)
            
            VERSION
                1.1, 20101103 by jmiguel@iaa.es
                1.2, 20110297 added out_dir parameter to skyfilter_single
        
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
        
        if len(near_list)<(self.HWIDTH+1):
            log.error("Wrong number of sky frames provided. Min number of sky frame is %d", 
                      self.MIN_SKY_FRAMES)
            return None
        
        # 0.1 Get the gain map
        if not os.path.exists( self.master_flat ):
            #raise Exception("Error, gain map file <%s> not found"%gain)
            #TODO: --> DONE try to compute GainMap using the given images !!!
            log.debug("---> creating gain map <----")
            output_fd, l_gainMap = tempfile.mkstemp(suffix='.fits', dir=self.out_dir)
            os.close(output_fd)
            os.unlink(l_gainMap) # we only need the name
            output_fd, files_list= tempfile.mkstemp(suffix='.list', dir=self.out_dir)
            os.close(output_fd)
            os.unlink(files_list) # we only need the name
            try:
                misc.utils.listToFile(near_list, files_list)
                # Note: we must normalize wrt chip 1, so norm=True
                superflat = reduce.SuperSkyFlat(files_list, l_gainMap, bpm=None, 
                                                norm=True, 
                                                temp_dir=self.temp_dir)
                superflat.create()
            except Exception,e:
                log.error("Error while creating gain map : %s", str(e))
                raise
        else: l_gainMap = self.master_flat
        
        #0.2 Check if GainMap need to be split
        gain_ext, g_next = self.split([l_gainMap])
        if g_next==1: gain_ext*=4
        print "\nGAINEXT=",gain_ext
        
        # 1. Split MEF file (into the self.out_dir)
        obj_ext, next = self.split(near_list)
        # it must return a list of list (one per each extension) 
        out_ext = []
        
        # 2. Process each extension
        for n in range(next):
            log.debug("===> Processing extension %d", n+1)
            # Create the temp list file of nearest (ar,dec,mjd) from current 
            # selected science file
            listfile = self.out_dir+"/nearfiles.list"
            misc.utils.listToFile(obj_ext[n], listfile)
            print "NEAR_FILES=", obj_ext[n]
            print "GAIN_EXT_N=",gain_ext[n][0]
            #Call external app skyfilter (irdr)
            hwidth = self.HWIDTH
            
            cmd = self.m_irdr_path+"/skyfilter_single %s %s %d nomask none %d %s"\
                        %(listfile, gain_ext[n][0], hwidth, file_pos, self.out_dir)
            print "CMD=",cmd
            e = misc.utils.runCmd( cmd )
            if e==1: # success
                fname = self.out_dir+"/"+os.path.basename(obj_ext[n][file_pos-1].replace(".fits", (".fits.skysub")))
                out_ext.append(fname)  
            else:
                log.error("Some error while subtracting sky in extension #%d# ", n+1)
                raise Exception("Some error while subtracting sky in extension #%d# "%(n+1))
        
        # 3. Package results back from each extension into a MEF file (only if nExt>1)
        if len(out_ext)>1:
            mef = misc.mef.MEF(out_ext)
            mef.createMEF(out_filename)
        elif len(out_ext)==1:
            shutil.move(out_ext[0], out_filename) 
        else:
            log.error("Some error while subtracting sky. No output produced.")
            
        return out_filename
    
    def getPointingOffsets (self, images_in=None,  p_offsets_file='/tmp/offsets.pap'):
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
        #output_list_file=self.out_dir+"/gpo_objs.pap"
        output_fd, output_list_file = tempfile.mkstemp(suffix='.pap', dir=self.out_dir)
        os.close(output_fd)
        os.unlink(output_list_file) # we only need the name
        
        log.debug("Creating OBJECTS images (SExtractor)....")
        
        if self.config_dict:
            mask_minarea = self.config_dict['offsets']['mask_minarea']
            mask_thresh = self.config_dict['offsets']['mask_thresh']
            satur_level = self.config_dict['offsets']['satur_level']
            single_p= self.config_dict['offsets']['single_point']
        else:
            mask_minarea = 5
            mask_thresh = 1.5
            satur_level = 300000
            single_p = True
            
        if images_in==None: # then we use the images ending with suffing in the output directory
            makeObjMask( self.out_dir+'*'+suffix , mask_minarea, mask_thresh, satur_level, \
                        output_list_file, single_point=single_p)
        elif os.path.isfile(images_in): # we use the given list of images
            makeObjMask( images_in , mask_minarea, mask_thresh, satur_level, \
                        output_list_file, single_point=single_p)
        else:
            log.error("Option not recognized !!!")
            raise Exception("Wrong input frames given")
        
        # STEP 2: Compute dither offsets (in pixles) using cross-correlation technique ==> offsets
        #>mosaic objfiles.nip $off_err > offsets1.nip
        search_box=10 # half_width of search box in arcsec (default 10)
        offsets_cmd=self.m_irdr_path+'/offsets '+ output_list_file + '  ' + str(search_box) + ' >' + p_offsets_file
        if misc.utils.runCmd( offsets_cmd )==0:
            log.critical("Some error while computing dither offsets")
            raise Exception("Some error while computing dither offsets")
        else:
            try:
                offsets_mat = numpy.loadtxt(p_offsets_file, usecols = (1,2,3)) # columns => (xoffset, yoffset, match fraction) in PIXELS
                # check if correlation overlap fraction is good enought
                if (offsets_mat[:,2]<self.MIN_CORR_FRAC).sum()>1:
                    log.critical("Some error while computing dither offsets. Overlap correlation fraction is < %f",self.MIN_CORR_FRAC)
                    raise Exception("Wrong overlap correlation fraction for translation offsets")
                    
            except IOError:
                log.critical("No offsets read. There may be some problem while computing translation offsets ....")
                raise Exception("No offsets read")
        
        log.debug("END of getPointingOffsets")                        
        return offsets_mat
                            
    def coaddStackImages(self, input='/tmp/stack.pap', gain=None, 
                         output=None, type_comb='average'):
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
        input_file = input
        if input_file==None:
            log.error("Bad input file provided !")
            return
        
        if gain==None:
            gain_file = self.out_dir+"/gain_"+self.m_filter+".fits"
        else:
            gain_file = gain
        
        if output==None:
            output_file = self.out_dir+"/coadd_"+self.m_filter+".fits"     
        else:
            output_file = output
            
        weight_file = output_file.replace(".fits",".weight.fits")
        
        # STEP 2: Run the coadd                                           
        if type_comb=='average': # (use IRDR::dithercubemean)
            prog = self.m_irdr_path+"/dithercubemean "
            cmd  = prog + " " + input_file + " " + gain_file + " " + output_file + " " + weight_file 
            e = misc.utils.runCmd( cmd )
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
        if self.config_dict:
            mask_minarea = self.config_dict['skysub']['mask_minarea']
            mask_thresh = self.config_dict['skysub']['mask_thresh']
            satur_level = self.config_dict['skysub']['satur_level']
        else:
            print "Program should never enter here !!!"
            #mask_minarea = 5
            #mask_thresh = 1.5
            #satur_level = 300000
                               
        # BUG ! -> input_file+"*" as first parameter to makeObjMask ! (2011-09-23)                                                               
        makeObjMask( input_file, mask_minarea, mask_thresh, satur_level,
                    outputfile=self.out_dir+"/objmask_file.txt", single_point=False)
        if os.path.exists(input_file+".objs"): 
            shutil.move(input_file+".objs", output_master_obj_mask)
            log.debug("New Object mask created : %s", output_master_obj_mask)
            
        # STEP 2: dilate mask (NOT DONE)
        dilate = False                                                     
        if dilate:
            log.info("Dilating image ....")
            prog = self.m_irdr_path+"/dilate "
            scale = 1.5 #mult. scale factor to expand object regions; default is 0.5 (ie, make 50%% larger)
            cmd  = prog + " " + output_master_obj_mask + " " + str(scale)
            # dilate will overwrite the master object mask
            
            e = misc.utils.runCmd( cmd )
            if e==0:
                log.debug("Some error while running command %s", cmd)
            else:
                log.debug("Succesful ending of createMasterObjMask")
                

        return output_master_obj_mask
                                        
    
    def cleanUpFiles(self, list_dirs):
        """Clean up files from the working directory, probably from the last execution"""
        """ TB reviewed """
        
        for out_dir in list_dirs:
            misc.fileUtils.removefiles(out_dir+"/*.fits")
            misc.fileUtils.removefiles(out_dir+"/c_*", out_dir+"/dc_*",
                                       out_dir+"/*.nip", out_dir+"/*.pap" )
            misc.fileUtils.removefiles(out_dir+"/coadd*", out_dir+"/*.objs",
                                       out_dir+"/uparm*", out_dir+"/*.skysub*")
            misc.fileUtils.removefiles(out_dir+"/*.head", out_dir+"/*.list",
                                       out_dir+"/*.xml", out_dir+"/*.ldac",
                                       out_dir+"/*.png" )

    def purgeOutput(self):
        """
        Purge the output dir in order to remove all the intermediate files
        """
        
        log.info("Purging the output dir ...")
        
        out_dir = self.out_dir
                 
        misc.fileUtils.removefiles(out_dir+"/*.ldac",out_dir+"/py-sex*",
                                   out_dir+"/*.objs")
        misc.fileUtils.removefiles(out_dir+"/coadd1*", out_dir+"/*_D.fits",
                                       out_dir+"/*_F.fits", out_dir+"/*_D_F.fits" )
        misc.fileUtils.removefiles(out_dir+"/gain*.fits", out_dir+"/masterObjMask.fits",
                                       out_dir+"/*.pap", out_dir+"/*.list", out_dir+"/superFlat.fits")
        misc.fileUtils.removefiles(out_dir+"/*.head", out_dir+"/*.txt",
                                       out_dir+"/*.xml")#, out_dir+"/*.png")
        
    ############# Calibration Stuff ############################################
    def buildCalibrations(self):
        """
        Build the whole master calibration files from the currect calibrations files
        found in the data set (darks, flats)
        """
        
        log.debug("Start builing the whole calibration files ...")
        # If not initialized, Init DB
        if self.db==None: self.__initDB()
        master_files = []
        try:
            master_files += self.buildMasterDarks()
            master_files += self.buildMasterDomeFlats()
            master_files += self.buildMasterTwFlats()
            master_files += self.buildMasterSuperFlats()
            master_files += self.buildGainMaps(type="all")
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
        l_gainmaps = []
        # 1. Look for master flat frames
        full_flat_list = []
        if type=="all":
            full_flat_list = self.db.GetFilesT(type="MASTER_SKY_FLAT", texp=-1, filter="ANY")
            full_flat_list+= self.db.GetFilesT(type="MASTER_TW_FLAT", texp=-1, filter="ANY")
            full_flat_list+= self.db.GetFilesT(type="MASTER_DOME_FLAT", texp=-1, filter="ANY")
        elif type=="sky":
            full_flat_list = self.db.GetFilesT(type="MASTER_SKY_FLAT", texp=-1, filter="ANY")
        elif type=="twlight":
            full_flat_list = self.db.GetFilesT(type="MASTER_TW_FLAT", texp=-1, filter="ANY")
        elif type=="dome":
            full_flat_list = self.db.GetFilesT(type="MASTER_DOME_FLAT", texp=-1, filter="ANY")
        else:
            log.error("Wrong type of master flat specified")
            raise Exception("Wrong type of master flat specified")
        
        if full_flat_list==None or len(full_flat_list)==0:
            log.warning("No master flat field found")
            return []

        # 2. Group by filter
        sorted_list = []
        for file in full_flat_list:
            fits = datahandler.ClFits(file)
            sorted_list.append((file, fits.getFilter()))
        
        # Sort out frames by FILTER
        sorted_list = sorted(sorted_list, key=lambda data_file: data_file[1])
        
        # 3. take the first group having the same FILTER
        last_filter = sorted_list[0][1]
        group = []
        build_master = False
        k = 0
        while k<len(sorted_list):
            while k<len(sorted_list) and sorted_list[k][1]==last_filter:
                group.append(sorted_list[k][0])
                k+=1
            #create the new master
            try:
                #TODO: here we should check if we have more that one master flat, and if have, then combine them ... 
                # generate a random filename for the master, to ensure we do not overwrite any file
                output_fd, outfile = tempfile.mkstemp(suffix='.fits', dir=self.out_dir)
                os.close(output_fd)
                os.unlink(outfile) # we only need the name
                # get gainmap parameters
                if self.config_dict:
                    mingain = self.config_dict['gainmap']['mingain']
                    maxgain = self.config_dict['gainmap']['maxgain']
                    nxblock = self.config_dict['gainmap']['nxblock']
                    nyblock = self.config_dict['gainmap']['nyblock']
                    nsigma = self.config_dict['gainmap']['nsigma']
                else:
                    print "Program should never enter here !"
                    mingain = 0.5
                    maxgain = 1.5
                    nxblock = 16
                    nyblock = 16
                    nsigma = 5    
                
                task=reduce.calGainMap.GainMap(group[0], outfile, bpm=None, do_normalization=True,
                                               mingain=mingain, maxgain=maxgain, nxblock=nxblock,
                                               nyblock=nyblock, nsigma=nsigma)
                
                out=None
                out=task.create()
                l_gainmaps.append(out) # out must be equal to outfile
            except Exception,e:
                log.error("Some error while creating gainmap: %s",str(e))
                raise e
            if k<len(sorted_list):
                # reset the new group
                group=[]
                last_filter=sorted_list[k][1]

        # insert products (gainmaps) into DB
        for f in l_gainmaps: self.db.insert(f)
        self.db.ListDataSet()  
        return l_gainmaps
    

    def buildMasterDarks(self):
        """
        Look for dark files in the data set, group them by DIT,NCOADD and create
        the master darks, as many as found groups
        
        @return: A list of the master darks created
         
        @todo: how to process dark series with an increasing TExp ??
        """
        
        log.debug("Creating Master darks...")
        l_mdarks = []
        
        # 1. Look for dark frames
        full_dark_list = self.getDarkFrames() # a list of tuples as (file,tExp)
        if len(full_dark_list)<=0:
            log.warning("No Dark frames found")
            return []
    
        # 2. take the first group having the same TEXP
        last_texp = full_dark_list[0][1]
        #outfile = self.out_dir+"/master_dark.fits" # class MasterDark will add as suffix (TEXP_NCOADDS)
        group=[]
        k=0
        while k<len(full_dark_list):
            while k<len(full_dark_list) and full_dark_list[k][1]==last_texp:
                group.append(full_dark_list[k][0])
                k+=1
            #create the new master
            try:
                # Generate (and create the file) a random filename for the master, 
                # to ensure we do not overwrite any file
                output_fd, outfile = tempfile.mkstemp(suffix='.fits', dir=self.out_dir)
                os.close(output_fd)
                os.unlink(outfile) # we only need the name
                task = reduce.calDark.MasterDark (group, self.temp_dir, outfile, texp_scale=False)
                out = task.createMaster()
                l_mdarks.append(out) # out must be equal to outfile
            except Exception,e:
                log.error("Some error while creating master dark: %s",str(e))
                log.error("but, proceding with next dark group ...")
                #raise e
            if k<len(full_dark_list):
                # reset the new group
                group = []
                last_texp = full_dark_list[k][1]
                
        # insert products (master darks) into DB
        for f in l_mdarks: self.db.insert(f)
        self.db.ListDataSet()         
        return l_mdarks        
        
    def buildMasterDomeFlats(self):
        """
        Look for DOME FLATS files in the data set, group them by FILTER and create
        the master DomeFlat, as many as found groups,i.e., filters.
        
        @return: A list of the master Dome Flats created
        
        """
        
        log.debug("Building Master DomeFlats...")
        l_mflats=[]
        # 1. Look for domeflat frames
        full_flat_list = self.getDomeFlatFrames()
        if len(full_flat_list)<=0:
            log.warning("No DomeFlatField frames found")
            return []
            
        # 2. take the first group having the same FILTER
        last_filter = full_flat_list[0][1]
        group=[]
        k=0
        while k<len(full_flat_list):
            while k<len(full_flat_list) and full_flat_list[k][1]==last_filter:
                group.append(full_flat_list[k][0])
                k+=1
            #create the new master
            try:
                # generate a random filename for the master, to ensure we do not overwrite any file
                output_fd, outfile = tempfile.mkstemp(suffix='.fits', dir=self.out_dir)
                os.close(output_fd)
                os.unlink(outfile) # we only need the name
                #outfile = self.out_dir+"/master_domeflat_%s.fits"%last_filter # added as suffix (FILTER)
                task = reduce.calDomeFlat.MasterDomeFlat(group, self.temp_dir, outfile, None)
                out=task.createMaster()
                l_mflats.append(out) # out must be equal to outfile
            except Exception,e:
                log.error("Some error while creating master DomeFlat: %s",str(e))
                log.error("but, proceding with next flat group ...")
                #raise e
            if k<len(full_flat_list):
                # reset the new group
                group=[]
                last_filter=full_flat_list[k][1]
                
        # insert products (master dome flats) into DB
        for f in l_mflats: self.db.insert(f)
        self.db.ListDataSet()  
        return l_mflats
    
    def buildMasterTwFlats(self):
        """
        Look for TwFlats files in the data set, group them by FILTER and create
        the master twflat, as many as found groups.
        
        @return: A list of the master TwFlats created
        """
        
        log.debug("Building Master TwFlats...")
        l_mflats=[]
        # 1. Look for twflat frames
        full_flat_list = self.getTwFlatFrames()
        if len(full_flat_list)<=0:
            log.warning("No TwFlatField frames found")
            return []
            
        # 2. take the first group having the same FILTER
        last_filter = full_flat_list[0][1]
        group=[]
        k=0
        while k<len(full_flat_list):
            while k<len(full_flat_list) and full_flat_list[k][1]==last_filter:
                group.append(full_flat_list[k][0])
                k+=1
            try:
                master_dark = self.db.GetFilesT('MASTER_DARK_MODEL') # could be > 1 master darks, then use the last(mjd sorted)
                # if required, master_dark will be scaled in MasterTwilightFlat class
                if len(master_dark)>0:
                    # generate a random filename for the master, to ensure we do not overwrite any file
                    output_fd, outfile = tempfile.mkstemp(suffix='.fits', dir=self.out_dir)
                    os.close(output_fd)
                    os.unlink(outfile) # we only need the name
                    #outfile = self.out_dir+"/master_twflat_%s.fits"%last_filter # added as suffix (FILTER)
                    
                    task = reduce.calTwFlat.MasterTwilightFlat(group, 
                                                               master_dark[-1], 
                                                               outfile, lthr=1000, 
                                                               hthr=100000, 
                                                               bpm=None)
                    out = task.createMaster()
                    
                    l_mflats.append(out) # out must be equal to outfile
                else:
                    log.error("MASTER_DARK_MODEL not found. Cannot build master TwFlat")
                    raise Exception("MASTER_DARK_MODEL not found. Cannot build master TwFlat.")
            except Exception,e:
                log.error("Some error while creating master TwFlat: %s",str(e))
                log.error("but, proceding with next twflat group ...")
                raise Exception("Cannot build master TwFlat: %s"%(str(e)))
            
            if k<len(full_flat_list):
                # reset the new group
                group = []
                last_filter = full_flat_list[k][1]
                
        # insert products (master twflats) into DB
        for f in l_mflats: self.db.insert(f)
        self.db.ListDataSet()
          
        return l_mflats
    
    def buildMasterSuperFlats(self):
        """
        Look for science/sky files in the data set, group them by FILTER and
        temporal proximity and then create the master superFlats, as many as
        found groups.
        
        Return: A list of master SuperFlat created
        
        @todo: need to be tested 
        """
        
        log.debug("Building Master SuperFlats...")
        
        # 1. Look for SCIENCE/SKY frames
        full_file_list = self.getObjectSequences()
        if len(full_file_list)<=0:
            log.warning("No sequence science frames found")
            return []
            
        l_mflats=[]       
        for seq in full_file_list:
            if len(seq)>1:
                try:
                    # generate a random filename for the master super flat, to ensure we do not overwrite any file
                    output_fd, output_path = tempfile.mkstemp(suffix='.fits', dir=self.out_dir)
                    os.close(output_fd)
                    os.unlink(output_path) # we only need the name
                    #outfile = self.out_dir+"/master_superflat_%s.fits"%last_filter # added as suffix (FILTER)
                    superflat = reduce.SuperSkyFlat(seq, output_path, bpm=None, norm=False, temp_dir=self.temp_dir)
                    out=superflat.create()
                    l_mflats.append(out)
                except Exception,e:
                    log.error("Some error while creating master SuperFlat: %s",str(e))
                    log.error("but, proceding with next group ...")
                    #raise e
                
        # insert products (master SuperFlats) into DB
        for f in l_mflats: self.db.insert(f)
        self.db.ListDataSet()  
        return l_mflats # a list of master super flats created
    
    def reduceSet(self, red_mode=None, seqs_to_reduce=None):
        """
        @summary: This is the main method for full DataSet reduction supposed 
        it was obtained with the PANIC OT. 
        
        Main steps:
        
         1. Get all OT sequences
         
         2. For seq in Sequence
    
            ReduceSeq(seq)
            
         3. Insert the results into the local DB 
         
        @param red_mode: reduction mode (lemon, quick, science); default mode 
        is 'quick'.
        
        @param seqs_to_reduce: list of sequence number [0,N-1] to be reduced;
        default (None), all sequences found will be reduced.
           
        @return: the number of sequences successfully reduced.
        
        """
        
        log.debug("[reduceSet] Dataset reduction process...")

        reduced_sequences = 0
        files_created = []
        failed_sequences = 0
        
        # set the reduction mode
        if red_mode is not None:
            self.red_mode = red_mode
        #else, keep the red_mode from the object initialization
            
        # Look for the sequences     
        sequences, seq_types = self.getSequences()
        
        
        # Check which sequences are required to reduce (-S command line param)
        # If no sequence number was specified, all seqs will be re-ordered and
        # processed 
        if seqs_to_reduce==None:
            seqs_to_reduce = range(len(sequences))    
            # Re-order the sequences by type: DARK, DOME_FLAT, TW_FLAT, SCIENCE
            # This is required because some calibration sequence could be created 
            # after the science sequence, and might happen no other calibration is 
            # available to process the current sequence.
            sequences, seq_types = self.reorder_sequences( sequences, seq_types)
        
        if len(sequences)==0:
            raise Exception("No well-defined sequence to process was found")
        
        k = 0
        for seq,type in zip(sequences, seq_types):
            if k in seqs_to_reduce:
                log.debug("A Sequence is going to be reduced ... ")
                try:
                    files_created += self.reduceSeq(seq, type)
                    reduced_sequences+=1
                except Exception,e:
                    # If an error happen while proecessing a sequence, we 
                    # do NOT STOP, but continue with the next ones. 
                    # However, if there is only one sequence, raise the exception,
                    # what it is very useful for the QL
                    failed_sequences +=1
                    log.error("[reduceSet] Cannot reduce sequence : \n %s \n %s"%(str(seq),str(e)))
                    log.debug("[reduceSet] Procceding to next sequence...")
                    if len(sequences)==1:
                        raise e
            k = k + 1
    
        # print out the results
        #failed_sequences = len(seqs_to_reduce)-reduced_sequences
        log.debug("[reduceSet] All sequences processed.")
        log.debug("[reduceSet] Files generated # %d #: ***"%len(files_created))
        for r_file in files_created: log.debug("\t    - %s"%r_file)
        log.debug("\t    Sequences failed  # %d #: ***"%failed_sequences)

        # WARNING : Purging output !!
        if self.config_dict['general']['purge_output']:
            self.purgeOutput()
        
        # In order to have a complete track of the processing, copy log to the 
        # output directory.    
        shutil.copy(misc.paLog.file_hd.baseFilename, self.out_dir)
        

        return files_created
   
    def reorder_sequences (self, sequences, seq_types):
        """
        Re-order the sequences by type: DARK, DOME_FLAT, TW_FLAT, SCIENCE
        This is required because some calibration sequence could be created 
        after the science sequence, and might happen no other calibration is 
        available to process the current sequence.
        
        @param sequences: list of filename list for each sequence
        @param seq_types: list of types (DARK, DOME_FLAT, TW_FLAT, SCIENCE) of
        each sequence of the first parameter.
        
        @return: two lists:
                - a list of lists of sequence files belonging to
                - a list with the types of the sequences (DARK, TW_FLAT, 
                DOME_FLAT, SCIENCE) ; see ClFits class for further details)
        
        """

        log.debug("[reorder_sequences] start ...")
        
        req_types_order = ['DARK', 'DOME_FLAT', 'TW_FLAT', 'SKY_FLAT', 'SCIENCE']
        new_sequences = []
        new_seq_types = []
        
        for r_type in req_types_order:
            for i,type in enumerate(seq_types):
                log.debug("TYPE=%s  R_TYPE=%s"%(type,r_type))
                if type == r_type:
                    new_sequences.append(sequences[i])
                    new_seq_types.append(type)
                    log.info("[reorder_sequences]: Sequence is of type [%s]"%(r_type.lower()))
                #else:
                #    log.debug("[reorder_sequences]: Sequence is NOT of type [%s]"%(r_type.lower()))
            

        return new_sequences, new_seq_types
    
    def calc(self, args):
        """
        Method used only to use with Pool.map_asycn() function
        """
        return self.reduceSingleObj(*args)
      
    def reduceSeq(self, sequence, type):
        """
        Reduce/process a produced OT-sequence of files (calibration, science).
    
       
        @param sequence: list of files of the sequence to be reduced
        @param type: type of sequence (see ClFits.type) (currently not used !) 
        
        @return: list of filenames created by the reduction proccess

        NOTE: 
        Calibration files will not be splited for the building of the master 
        calibration file, but science files will be splitted. 
        In principle, not too much time could be saved if we split the calibration
        files during the building of master calibrations.
        However, the master calibration files will be splitted at the stage of
        the processing of scince files. 
        
        """
        
        
        log.debug("[reduceSeq] Starting ...")
        
        # start creating the pool of process
        n_cpus = self.config_dict['general']['ncpus']
        #n_cpus = multiprocessing.cpu_count()

        #Finally, it was not possible to create 'pool' as a global variable becasue
        #there are problems with sharing a global pool (or whatever variable) with
        #the Process created in the QL module for dispatching the on-line reduction.
        #So, what we do now is create a new pool everytime we need to do any iraf
        #depending task. In general, we create the pool for whatever task in 
        #this function, but in serial reduction of SCI seqs, the pool will not
        #be used.
        pool = multiprocessing.Pool(processes=n_cpus)
        
        files_created = []
        
        # Take the first file of the sequence in order to find out the type of 
        # the sequence
        fits = datahandler.ClFits(sequence[0])
        
        
        if fits.isDark():
            log.debug("[reduceSeq] A Dark sequence is going to be reduced: \n%s"%str(sequence))
            try:
                # Generate (and create the file) a random filename for the master, 
                # to ensure we do not overwrite any file
                output_fd, outfile = tempfile.mkstemp(suffix='.fits', prefix='mDark_', dir=self.out_dir)
                os.close(output_fd)
                os.unlink(outfile) # we only need the name
                
                #check and collapse if required (cube images)
                sequence = misc.collapse.collapse(sequence)
                
                # Check for EXPT in order to know how to create the master dark (dark model or fixed EXPT)     
                if (self.checkData(chk_shape=True, chk_filter=True, chk_type=True, chk_expt=True, 
                                   chk_itime=True, chk_ncoadd=True, chk_cont=True, 
                                   chk_readmode=True, file_list=sequence)==True):
                    # Orthodox master dark -- same EXPTIME & NCOADDS
                    log.debug("Found dark series with equal EXPTIME. Master dark is going to be created")
                    task = reduce.calDark.MasterDark (sequence, self.temp_dir, 
                                                      outfile, texp_scale=False,
                                                      bpm=None, normalize=False)
                    #out = task.createMaster()
                    
                    red_parameters = ()
                    result = pool.apply_async(task.createMaster, 
                                                 red_parameters)
                    result.wait()
                    out = result.get()
                    
                    log.critical("OUTPUT file generated %s"%out)
                    pool.close()
                    pool.join()

                    
                else:
                    log.info("Found a dark series of frames with different EXPTIME: Dark model will be created")
                    use_dark_model = True
                    if use_dark_model==True and self.red_mode !="quick":
                        # Build master dark model from a dark serie
                        task = reduce.calDarkModel.MasterDarkModel (sequence, 
                                                                    self.temp_dir, 
                                                                    outfile)
                        #out = task.createDarkModel()
                        
                        red_parameters = ()
                        result = pool.apply_async(task.createDarkModel, 
                                                 red_parameters)
                        result.wait()
                        out = result.get()
                    
                        pool.close()
                        pool.join()
                        
                    else:
                        log.warning("Found a dark serie with diff EXPTIME, but dark model processing not activated")
                        #raise Exception("Dark series are not processed in quick mode")
                        out = None
                 
                if out!=None: files_created.append(out) # out must be equal to outfile
            except Exception,e:
                log.error("[reduceSeq] Some error while creating master DARK: %s",str(e))
                raise e
        elif fits.isDomeFlat():
            log.debug("[reduceSeq] A DomeFlat sequence is going to be reduced: \n%s"%str(sequence))
            try:
                # generate a random filename for the master, to ensure we do not overwrite any file
                output_fd, outfile = tempfile.mkstemp(suffix='.fits', prefix='mDFlat_', dir=self.out_dir)
                os.close(output_fd)
                os.unlink(outfile) # we only need the name

                #check and collapse if required (cube images)
                sequence = misc.collapse.collapse(sequence)

                m_smooth = self.config_dict['dflats']['median_smooth']
                task = reduce.calDomeFlat.MasterDomeFlat(sequence, 
                                                         self.temp_dir, outfile, 
                                                         normal=True, # it is also done in calGainMap
                                                         median_smooth=m_smooth)
                #out = task.createMaster()
                
                red_parameters = ()
                result = pool.apply_async(task.createMaster, 
                                                 red_parameters)
                result.wait()
                out = result.get()
                    
                pool.close()
                pool.join()
                
                if out!=None: files_created.append(out) # out must be equal to outfile
            except Exception,e:
                log.error("[reduceSeq] Some error while creating master DomeFlat: %s",str(e))
                raise e
        elif fits.isTwFlat():
            log.debug("[reduceSeq] A TwFlat sequence is going to be reduced: \n%s"%str(sequence))
            try:
                #Look for the required MasterDark (any ExpTime);first in the Local DB (current RS), 
                #and if anyone found, then in the External DB
                #Local (ExpTime is not a constraint)
                master_dark = self.db.GetFilesT('MASTER_DARK_MODEL') # could there be > 1 master darks, then use the last(mjd sorted)
                
                #External (ExpTime is not a constraint)
                if len(master_dark)==0 and self.ext_db!=None:
                    master_dark = self.ext_db.GetFilesT('MASTER_DARK_MODEL') # could there be > 1 master darks, then use the last(mjd sorted)
                
                # if required, master_dark will be scaled in MasterTwilightFlat class
                if len(master_dark)>0:
                    log.debug("[reduceSeq] MASTER_DARK_MODEL found --> %s"%master_dark[-1])
                    # generate a random filename for the masterTw, to ensure we do not overwrite any file
                    output_fd, outfile = tempfile.mkstemp(suffix='.fits', prefix='mTwFlat_', dir=self.out_dir)
                    os.close(output_fd)
                    os.unlink(outfile) # we only need the name

                    #check and collapse if required (cube images)
                    sequence = misc.collapse.collapse(sequence)

                    m_smooth = self.config_dict['twflats']['median_smooth']
                    
                    task = reduce.calTwFlat.MasterTwilightFlat(sequence, 
                                                               master_dark[-1], 
                                                               outfile, 
                                                               lthr=1000, 
                                                               hthr=100000, 
                                                               bpm=None,
                                                               normal=True, # it is also done in calGainMap
                                                               temp_dir=self.temp_dir,
                                                               median_smooth=m_smooth)
                    #out = task.createMaster()
                    red_parameters = ()
                    result = pool.apply_async(task.createMaster, 
                                                 red_parameters)
                    result.wait()
                    out = result.get()
                    
                    pool.close()
                    pool.join()
                    
                    if out!=None: files_created.append(out) # out must be equal to outfile
                else:
                    # should we create master dark ??
                    log.error("[reduceSeq] MASTER_DARK_MODEL not found. Cannot build master TwFlat")
                    raise Exception("[reduceSeq] MASTER_DARK_MODEL not found")
            except Exception,e:
                log.error("[reduceSeq] Some error while creating master TwFlat: %s",str(e))
                raise e
        elif fits.isFocusSerie():
            log.debug("[reduceSeq] Focus Serie reduction is not yet implemented")
            # TODO
            pass
        elif fits.isScience():
            l_out_dir = ''
            results = None
            out_ext = []
            log.debug("[reduceSeq] Reduction of SCIENCE Sequence: \n%s"%str(sequence))
            if len(sequence) < self.config_dict['general']['min_frames']:
                log.info("[reduceSeq] Found a too SHORT Obs. object sequence.\n\
                 Only %d frames found. Required >%d frames"%(len(sequence),
                                                             self.config_dict['general']['min_frames']))
                raise Exception("Found a short Obs. object sequence. \n\
                Only %d frames found. Required >%d frames" %(len(sequence),
                                                            self.config_dict['general']['min_frames']))
            else:
                #check and collapse if required (cube images)
                sequence = misc.collapse.collapse(sequence)

                # Get calibration files
                dark, flat, bpm = None, None, None
                if self.red_mode == "quick":
                    if self.apply_dark_flat==1 or self.apply_dark_flat==2: 
                        dark, flat, bpm = self.getCalibFor(sequence)
                else:
                    dark, flat, bpm = self.getCalibFor(sequence)
                    # Return 3 filenames of master calibration frames (dark, flat, bpm), 


                # Check and split files if required
                obj_ext, next = self.split(sequence) # it must return a list of list (one per each extension)
                dark_ext, cext = self.split([dark])
                flat_ext, cext = self.split([flat])
                bpm_ext, cext = self.split([bpm])
                parallel = self.config_dict['general']['parallel']
                
                if parallel==True:
                    ######## Parallel #########
                    log.info("[reduceSeq] Entering PARALLEL data reduction ...")
                    try:
                        # Map the parallel process
                        #n_cpus = self.config_dict['general']['ncpus'] #anymore used, instead cpu_count()
                        # pool is defined at the end of the file, as a global variable

                        ##results = pprocess.Map(limit=n_cpus, reuse=1) 
                        # IF reuse=0, it block the application !! I don't know why ?? 
                        # though in pprocess examples it works! 
                        ##calc = results.manage(pprocess.MakeReusable(self.reduceSingleObj))

                        results = []
                          
                        for n in range(next):
                            ## only a test to reduce Q01
                            #log.critical("only a test to reduce Q01")
                            #if n!=0 and n!=1: continue
                            ## end-of-test
                            log.info("[reduceSeq] ===> (PARALLEL) Reducting extension %d", n+1)
                            ## At the moment, we have the first calibration file for each extension; what rule could we follow ?
                            if dark_ext==[]: mdark = None
                            else: mdark = dark_ext[n][0] # At the moment, we take the first calibration file found for each extension
                            if flat_ext==[]: mflat = None
                            else: mflat = flat_ext[n][0] # At the moment, we take the first calibration file found for each extension
                            if bpm_ext==[]: mbpm = None
                            else: mbpm = bpm_ext[n][0] # At the moment, we take the first calibration file found for each extension
                            
                            l_out_dir = self.out_dir + "/Q%02d" % (n+1)
                            if not os.path.isdir(l_out_dir):
                                try:
                                    os.mkdir(l_out_dir)
                                except OSError:
                                    log.error("[reduceSeq] Cannot create output directory %s",l_out_dir)
                            else: self.cleanUpFiles([l_out_dir])
                            
                            # async call to procedure
                            extension_outfilename = l_out_dir + "/" + os.path.basename(self.out_file.replace(".fits",".Q%02d.fits"% (n+1)))
                            ##calc( obj_ext[n], mdark, mflat, mbpm, self.red_mode, l_out_dir, extension_outfilename)
                            
                            red_parameters = (obj_ext[n], mdark, mflat, mbpm, 
                                              self.red_mode, l_out_dir, 
                                              extension_outfilename)
                            log.debug("Calling parameters ---> %s"%(str(red_parameters)))
                            
                            # Notice that the results will probably not come out 
                            # of the output queue in the same in the same order 
                            # as the corresponding tasks were put on the input 
                            # Pool.  If it is important to get the results back
                            # in the original order then consider using `Pool.map()` 
                            # or `Pool.imap()` (which will save on the amount of 
                            # code needed anyway).
                            #results += [pool.apply_async(self.reduceSingleObj, 
                            #                             red_parameters)]
                            pool = multiprocessing.Pool(2)
                            results += [pool.map_async(self.calc,
                                                      [red_parameters])]         
                        # Here is where we WAIT (BLOCKING) for the results 
                        # (result.get() is a blocking call).
                        # If the remote call raised an exception then that 
                        # exception will be reraised by get().

                        for result in results:
                            result.wait()
                            out_ext.append(result.get()[0]) # the 0 index is *ONLY* required if map_async is used !!!

                        ##for result in results:
                        ##    out_ext.append(result)
                            
                        #Prevents any more tasks from being submitted to the pool. 
                        #Once all the tasks have been completed the worker 
                        #processes will exit.
                        pool.close()
                        #Wait for the worker processes to exit. One must call 
                        #close() or terminate() before using join().
                        pool.join()
                        log.critical("[reduceSeq] DONE PARALLEL REDUCTION ")
                            
                    except Exception,e:
                        log.error("[reduceSeq] Error while parallel data reduction ! --> %s",str(e))
                        raise e
                    
                else:
                    ######## Serial #########
                    for n in range(next):
                        log.info("[reduceSeq] Entering SERIAL science data reduction ...")    
                        log.info("[reduceSeq] ===> (SERIAL) Reducting extension %d", n+1)
                        
                        ## At the moment, we have the first calibration file for each extension; what rule could we follow ?
                        if dark_ext==[]: mdark = None
                        else: mdark = dark_ext[n][0]  # At the moment, we have the first calibration file for each extension
                        
                        if flat_ext==[]: mflat = None
                        else: mflat = flat_ext[n][0]  # At the moment, we have the first calibration file for each extension
                        
                        if bpm_ext==[]: mbpm = None
                        else: mbpm = bpm_ext[n][0]    # At the moment, we have the first calibration file for each extension
                        
                        # Debug&Test
                        #if n==1: return None,None# only for a TEST !!!
                        #if n!=1: continue
                        # 
                        
                        try:
                            out_ext.append(self.reduceSingleObj(obj_ext[n],
                                                                mdark, mflat,
                                                                mbpm, self.red_mode,
                                                                out_dir=self.out_dir,
                                                                output_file = self.out_dir + \
                                                                "/out_Q%02d.fits"%(n+1)))
                        except Exception,e:
                            log.error("[reduceSeq] Error while serial data reduction of extension %d of object sequence", n+1)
                            raise e
        
            # If red_mode is 'lemon',then no warping of frames is required
            if self.red_mode=='lemon':
                files_created = out_ext
                log.info("*** Obs. Sequence LEMON-reduced. ***") 
                return files_created
            
            # if all reduction were fine, now join/stich back the extensions in a wider frame
            seq_result_outfile = self.out_file.replace(".fits","_SEQ.fits")
            if len(out_ext) >1:
                print "OUT_EXT=",out_ext
                log.debug("[reduceSeq] *** Creating final output file *WARPING* single output frames....***")
                #option 1: create a MEF with the results attached, but not warped
                #mef=misc.mef.MEF(outs)
                #mef.createMEF(self.out_file)
                #option 2(current): SWARP result images to register the N-extension into one wide-single extension
                log.debug("*** Coadding/Warping overlapped files....")
                
                swarp = astromatic.SWARP()
                swarp.config['CONFIG_FILE'] = self.config_dict['config_files']['swarp_conf'] 
                swarp.ext_config['COPY_KEYWORDS'] = 'OBJECT,INSTRUME,TELESCOPE,IMAGETYP,FILTER,FILTER1,FILTER2,SCALE,MJD-OBS'
                swarp.ext_config['IMAGEOUT_NAME'] = seq_result_outfile
                swarp.ext_config['WEIGHTOUT_NAME'] = self.out_file.replace(".fits",".weight.fits")
                swarp.ext_config['WEIGHT_TYPE'] = 'MAP_WEIGHT'
                swarp.ext_config['WEIGHT_SUFFIX'] = '.weight.fits'
                
                try:
                    swarp.run(out_ext, updateconfig=False, clean=False)
                except SWARPException, e:
                    log.error("Error while running SWARP")
                    raise e
                except Exception, e:
                    log.error("Unknow error while running SWARP: %s",str(e))
                    raise e
                
                files_created.append(seq_result_outfile)
                log.info("*** Obs. Sequence reduced. File %s created.  ***", 
                         seq_result_outfile)
                
            elif len(out_ext)==1:
                shutil.move(out_ext[0], seq_result_outfile)
                
                files_created.append(seq_result_outfile)
                log.info("*** Obs. Sequence reduced. File %s created.  ***", 
                         seq_result_outfile)
            else:
                log.error("[reduceSeq] No output files generated by the current Obj.Seq. \
                data reduction ....review your logs files")
            
        else:
            log.error("[reduceSeq] Cannot identify the type of the sequence to reduce ...")
            raise Exception("[reduceSeq] Cannot identify the type of the sequence to reduce ...")    
                    
        # Insert all created files into the DB                
        for file in files_created:
            if file!=None: self.db.insert(file)
        
        # not sure if is better to do here ??? 
        #self.purgeOutput()
            
        return files_created
 
    def reduceSet_deprecated(self, red_mode="quick"):
        """
        The main method for full DataSet reduction.
        
        Main steps:
        
         1. Build master calibration files
        
         2. Find out Objects/pointings/filter
        
         3. For each object/pointing/filter do :
         
            3.1 c=get_calibration_files (depend on filter/text)
            3.2 o = get_objects/pointings()
            3.3 exts=split(c,o)
            3.4 for each extension e do:
                out+=reduce_obj(o_e, c_e)
            3.5 create_joined_mef(out)
            3.6 Add objects to catalog (?)
            
         
        @return: Return a list of N files produced as result of the data reduction of
        the N sequecend found.
        
        @attention: Must be 
        """
       
        log.info("Starting Reduction of data set ...")
        
        # Clean old files in out_dir's
        self.cleanUpFiles([self.out_dir]) #"/data/out2","/data/out3","/data/out4"])
        
        # Init DB
        if self.db==None: self.__initDB()
        
        if red_mode=="quick": # we suppose it's Quick-Look mode, thus not use calibration files
            log.debug("Quick reduction mode, no calibration files will be built")
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
            
        
        #### Get all the SCIENCE sequences ###
        sequences = self.getObjectSequences()
        
        i = 0 # index for object sequences 
        out_ext = [] # it will store the partial extension reduced output filenames
        seq_result_outfile = "" # will store the full-frame out filename of the current reduced sequence  
        seq_outfile_list = [] # will store the full-frame result for each sequence-data-reduction 
        
        # For each object-sequence, split and reduce de sequence
        for obj_seq in sequences:  #sequences is a list of list, and obj_seq is a list
            log.debug("===> Reduction of obj_seq %d",i)
            if len(obj_seq)<4:
                log.info("Found and skipping  SHORT Obs. object sequence. Only %d frames found. Required >4 frames",len(obj_seq))
                continue # continue with next sequence !!
                #raise Exception("Found a short Obs. object sequence. Only %d frames found. Required >4 frames",len(obj_seq))
                #return []
            else:
                #avoid call getCalibFor() when red_mode="quick"
                if red_mode == "quick":
                    dark,flat,bpm = None,None,None
                else:
                    dark, flat, bpm = self.getCalibFor(obj_seq)
                    # return 3 filenames of master calibration frames (dark, flat, bpm)
                obj_ext, next = self.split(obj_seq) # it must return a list of list (one per each extension)
                dark_ext, cext = self.split([dark])
                flat_ext, cext = self.split([flat])
                bpm_ext, cext = self.split([bpm])
                parallel = self.config_dict['general']['parallel']
                
                if parallel==True:
                    log.info("Entering parallel data reduction ...")
                    try:
                        # Map the parallel process
                        n_cpus = self.config_dict['general']['ncpus']
                        results = pprocess.Map(limit=n_cpus, reuse=1) 
                        # IF reuse=0, it block the application !! I don't know why ?? 
                        # though in pprocess examples it works! 
                        calc = results.manage(pprocess.MakeReusable(self.reduceSingleObj))
                        for n in range(next):
                            log.info("===> (PARALLEL) Reducting extension %d", n+1)
                            ## At the moment, we have the first calibration file for each extension; what rule could we follow ?
                            if dark_ext==[]: mdark = None
                            else: mdark=dark_ext[n][0] # At the moment, we take the first calibration file found for each extension
                            if flat_ext==[]: mflat = None
                            else: mflat=flat_ext[n][0] # At the moment, we take the first calibration file found for each extension
                            if bpm_ext==[]: mbpm = None
                            else: mbpm = bpm_ext[n][0] # At the moment, we take the first calibration file found for each extension
                            
                            l_out_dir = self.out_dir + "/Q%02d" % (n+1)
                            if not os.path.isdir(l_out_dir):
                                try:
                                    os.mkdir(l_out_dir)
                                except OSError:
                                    log.error("Cannot create output directory %s",l_out_dir)
                            else: self.cleanUpFiles([l_out_dir])
                            
                            # async call to procedure
                            extension_outfilename = l_out_dir + "/" + os.path.basename(self.out_file.replace(".fits",".Q%02d.fits"% (n+1)))
                            calc( obj_ext[n], mdark, mflat, mbpm, red_mode, l_out_dir, extension_outfilename)
                        
                        # Here is where we WAIT (BLOCKING) for the results (iteration is a blocking call)
                        for result in results:
                            #print "RESULT=",result
                            out_ext.append(result)

                        log.critical("DONE PARALLEL REDUCTION ")
                            
                        #pool=multiprocessing.Pool(processes=2)
                        # !! ERROR !!, because multiprocessing.pool.map() does not support class methods
                        #pool.map(self.reduceSingleObj,[[obj_ext[0], mdark, mflat, mbpm, red_mode, self.out_dir+"/out_Q01.fits"],\
                        #                     [obj_ext[1], mdark, mflat, mbpm, red_mode, self.out_dir+"/out_Q02.fits"]])
                            
                    except Exception,e:
                        log.error("Error while parallel data reduction ! --> %s",str(e))
                        raise e
                    
                else:
                    for n in range(next):
                        log.info("Entering sequencial data reduction ...")    
                        ################## only for debug purposes
                        #if n!=2: continue
                        ################## end of debug 
                        log.info("===> (SERIAL) Reducting extension %d", n+1)
                        ## At the moment, we have the first calibration file for each extension; what rule could we follow ?
                        if dark_ext==[]: mdark=None
                        else: mdark=dark_ext[n][0]  # At the moment, we have the first calibration file for each extension
                        if flat_ext==[]: mflat=None
                        else: mflat=flat_ext[n][0]  # At the moment, we have the first calibration file for each extension
                        if bpm_ext==[]: mbpm=None
                        else: mbpm=bpm_ext[n][0]    # At the moment, we have the first calibration file for each extension
                        try:
                            out_ext.append(self.reduceSingleObj(obj_ext[n], mdark, mflat, mbpm, red_mode, out_dir=self.out_dir,\
                                                          output_file=self.out_dir+"/out_Q%02d.fits"%(n+1)))
                        except Exception,e:
                            log.error("Some error while reduction of extension %d of object sequence %d", n+1, i)
                            raise e
                            #continue
        
            # if all reduction were fine, now join/stich back the extensions in a widther frame
            seq_result_outfile = self.out_file.replace(".fits","_SEQ%02d.fits" %(i))
            if len(out_ext) >1 :
                log.debug("*** Creating final output file *WARPING* single output frames....***")
                #option 1: create a MEF with the results attached, but not warped
                #mef=misc.mef.MEF(outs)
                #mef.createMEF(self.out_file)
                #option 2(current): SWARP resulted images to register the N-extension into one wide-single extension
                log.debug("*** Coadding/Warping overlapped files....")
                swarp = astromatic.SWARP()
                swarp.config['CONFIG_FILE'] = self.config_dict['config_files']['swarp_conf'] 
                swarp.ext_config['COPY_KEYWORDS'] = 'OBJECT,INSTRUME,TELESCOPE,IMAGETYP,FILTER,FILTER1,FILTER2,SCALE,MJD-OBS'
                swarp.ext_config['IMAGEOUT_NAME'] = seq_result_outfile
                swarp.ext_config['WEIGHTOUT_NAME'] = self.out_file.replace(".fits",".weight.fits")
                swarp.ext_config['WEIGHT_TYPE'] = 'MAP_WEIGHT'
                swarp.ext_config['WEIGHT_SUFFIX'] = '.weight.fits'
                swarp.run(out_ext, updateconfig=False, clean=False)
                
                seq_outfile_list.append(seq_result_outfile)
                log.info("*** Obs. Sequence reduced. File %s created.  ***", 
                         seq_result_outfile)
                
            elif len(out_ext)==1:
                shutil.move(out_ext[0], seq_result_outfile)
                
                seq_outfile_list.append(seq_result_outfile)
                log.info("*** Obs. Sequence reduced. File %s created.  ***", 
                         seq_result_outfile)
            else:
                log.error("No output files generated by the current Obj.Seq. \
                data reduction ....review your logs files")
                
            i+=1
            
                    
            
        #log.critical(" TBD : for each object sequence, we need a different final output file !!!!")
        log.info("*******************************************************************")
        log.info("*** All Obs.Sequences (%d) has been reduced. Congratulations !!****", i)
        log.info("*******************************************************************")
        
        for f in seq_outfile_list: self.db.insert(f)
        self.db.ListDataSet()
        
        return seq_outfile_list
        
    def reduceSingleObj(self, obj_frames, master_dark, master_flat, master_bpm, 
                  red_mode, out_dir, output_file):
        
        """ 
        Main reduction procedure. 
        Given a set of object(science) frames and (optionally) master calibration files,
        run the data reduction of the observing object sequence, producing an reduced
        ouput frame if no error; otherwise return None or raise exception.
            
        NOTE: Currently this method only accepts single FITS files (not MEF), it means
              the splitting must be done previusly to call this method.
       
        Parameters
        ----------
        
        obj_frames : list
            list of files to be processed
        
        master_dark : str
            Master dark filename to be used in the processing
            
        master_flat : str
            Master flat filename to be used in the processing
            
        
        
        Returns
        -------
        
        output_file : str 
            If 'red_mode' = 'quick' or 'science', then the function return the 
            coadded frame obtained from the reduction, both 'quick' and 'science' 
            reduction mode.
            If 'red_mode' = 'lemon', the method return a file listing the files
            obtained of the pre-processing done (dark, flat and sky subtraction) 
        
        """
        
        log.info("##################################")
        log.info("#### Starting Object Data Reduction #####")
        log.info("#### MODE = %s  ", self.red_mode)
        log.info("#### OUT_DIR = %s ",out_dir)
        log.info("#### OUT_FILE = %s ", output_file)
        log.info(" ----------------------------------")
        c_filter = datahandler.ClFits(obj_frames[0]).getFilter()
        log.info("#### FILTER = %s", c_filter)
        log.info("#### MASTER_DARK = %s ", master_dark)
        log.info("#### MASTER_FLAT = %s ", master_flat)
        log.info("#### MASTER_BPM = %s ", master_bpm)
        log.info("##################################")
        #print "OBJS =",obj_frames
        
        # set the reduction mode
        if red_mode != None: self.red_mode = red_mode
        #else, keep the red_mode value from the object
        
        
        # Clean old files 
        #self.cleanUpFiles()
        
        # Change cwd to self.out_dir
        old_cwd = os.getcwd()
        os.chdir(out_dir) 
        
        # Copy/link source files (file or directory) to reduce to the working directory
        # and Initialize self.m_LAST_FILES
        if not os.path.dirname(obj_frames[0])==out_dir:
            misc.fileUtils.linkSourceFiles(obj_frames, out_dir)
            self.m_LAST_FILES = [out_dir+"/"+os.path.basename(file_i) for file_i in obj_frames]
        else:
            self.m_LAST_FILES = obj_frames
            
        #print "\nSOURCES TO BE REDUCED:"
        #print   "====================="
        #print self.m_LAST_FILES
        #print   "====================="
        
        ########################################################################
        # 0 - Some checks (filter, ....) 
        ########################################################################
        # TODO : it could/should be done in reduceSeq, to avoid the spliting ...??
        log.info("**** Data Validation ****")
        if self.check_data:
            if (self.checkData(chk_shape=True, chk_filter=True, chk_type=False, chk_expt=True, 
                               chk_itime=True, chk_ncoadd=True, chk_readmode=True)==True):
                log.debug("Data checking was OK !")
            else:
                raise Exception("Mismatch in data checking !")
        
        ########################################################################
        # 00 - Sort out data by MJD (self.m_LAST_FILES)
        ########################################################################
        # in principle, it is suppossed they are already ascending sorter, but...
        try:
            self.m_LAST_FILES = self.sortOutData() 
        except:
            raise
        
        ########################################################################
        # 000 - Find out dither mode
        ########################################################################
        try:
            self.obs_mode = self.getObsMode()  # overwrite initial given observing mode
        except:
            raise
        
        log.info( "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" )
        log.info( "OBSERVING SEQUENCE DETECTED. OBS_MODE= %s", self.obs_mode)
        log.info( "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" )
        
        # keep a copy of the original file names
        self.m_rawFiles = self.m_LAST_FILES        
        
        ########################################################################
        # 1 - Apply dark, flat to ALL files 
        ########################################################################
        if self.apply_dark_flat==1 and (master_dark!=None or master_flat!=None):
            log.info("**** Applying dark and Flat ****")
            res = reduce.ApplyDarkFlat(self.m_LAST_FILES, 
                                       master_dark, 
                                       master_flat, 
                                       out_dir)
            self.m_LAST_FILES = res.apply()
        
        ########################################################################
        # 2 - Compute Super Sky Flat-Field --> GainMap
        #      The 'local_master_flat' is ONLY used to build the GainMap, and
        #      can come from the 'master_flat' provided by the user or from the 
        #      just created 'super_flat' with sci images 
        ########################################################################
        if master_flat==None:
            try:
                # - Find out what kind of observing mode we have (dither, ext_dither, ...)
                log.info('**** Computing Super-Sky Flat-Field (local_master_flat) ****')
                local_master_flat = out_dir+"/superFlat.fits"
                if self.obs_mode=="dither":
                    log.debug("---> dither sequece <----")
                    misc.utils.listToFile(self.m_LAST_FILES, out_dir+"/files.list")
                    superflat = reduce.SuperSkyFlat(out_dir+"/files.list", 
                                                    local_master_flat, bpm=None, 
                                                    norm=True, 
                                                    temp_dir=self.temp_dir)
                    superflat.create()
                elif (self.obs_mode=="dither_on_off" or
                      self.obs_mode=="dither_off_on" or
                      self.obs_mode=="other"):
                    log.debug("----> EXTENDED SOURCE !!! <----")
                    sky_list = self.getSkyFrames()
                    misc.utils.listToFile(sky_list, out_dir+"/files.list")
                    superflat = reduce.SuperSkyFlat(out_dir+"/files.list", 
                                                    local_master_flat, bpm=None, 
                                                    norm=True, 
                                                    temp_dir=self.temp_dir)
                    superflat.create()                            
                else:
                    log.error("Dither mode not supported")
                    raise Exception("Error, dither mode not supported")
            except Exception,e:
                raise e    
        else:
            local_master_flat = master_flat 
            log.info("Using the given (dome or twlight) master flat")
                 
        ########################################################################
        # 3 - Compute Gain map and apply BPM
        ########################################################################
        log.info("**** Computing gain-map from ****")
        gainmap = out_dir+'/gain_'+self.m_filter+'.fits'
        # get gainmap parameters
        if self.config_dict:
            mingain = self.config_dict['gainmap']['mingain']
            maxgain = self.config_dict['gainmap']['maxgain']
            nxblock = self.config_dict['gainmap']['nxblock']
            nyblock = self.config_dict['gainmap']['nyblock']
            nsigma = self.config_dict['gainmap']['nsigma']
        else:
            mingain = 0.5
            maxgain = 1.5
            nxblock = 16
            nyblock = 16
            nsigma = 5
            
        # When gainmap is created (from dome or sky flats), it must be normalized
        # wrt mode of chip 1 to get gain differences, set bad pixels, 
        # outlier set =0 (e.g. pixels deviating >5 sigma from local median,
        # pixels deviating >30%(?),...
        # do_normalization=False because it is suppossed that FF is already normalized
        g = reduce.calGainMap.GainMap(local_master_flat, gainmap, bpm=master_bpm, 
                                    do_normalization=False, # because it is suppossed that FF is already normalized 
                                    mingain=mingain, 
                                    maxgain=maxgain, nxblock=nxblock,
                                    nyblock=nyblock, nsigma=nsigma)
        g.create() 
           
        ########################################################################
        # Add external Bad Pixel Map to gainmap (maybe from master DARKS,FLATS ?)
        ########################################################################
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
                        
        ########################################################################
        # 4 - First Sky subtraction (IRDR) - sliding window technique
        ########################################################################
        log.info("**** 1st Sky subtraction (without object mask) ****")
        misc.utils.listToFile(self.m_LAST_FILES, out_dir+"/skylist1.list")
        # return only filtered images; in exteded-sources, sky frames  are not included 
        self.m_LAST_FILES = self.skyFilter(out_dir+"/skylist1.list",
                                           gainmap, 'nomask', self.obs_mode)       

        ########################################################################
        # 4.1 - Divide by the master flat after sky subtraction ! 
        # some people think it is better do it now (M.J.Irwing, CASU)
        # But, dark is not applied/required, as it was implicity done when sky subtr.
        # However, it does not produce good results, it looks like if we undo the 
        # sky subtraction, because sky subtraction is quite related with flatfielding
        # For further details, check TN and Wei-Hao (SIMPLE) mails. 
        # So, it is implemented here only for academic purposes !
        ########################################################################
        if self.apply_dark_flat==2 and master_flat!=None:
            log.info("**** Applying Flat AFTER sky subtraction ****")
            res = reduce.ApplyDarkFlat(self.m_LAST_FILES, 
                                       None,  
                                       master_flat, 
                                       out_dir)
            self.m_LAST_FILES = res.apply()

        ########################################################################
        # 4.2 - LEMON connection - End here for LEMON-1 processing    
        ########################################################################
        """if self.red_mode=='lemon':
            misc.utils.listToFile(self.m_LAST_FILES, out_dir+"/files_skysub.list")
            log.info("1st Skysubtraction done !")
            #return out_dir+"/files_skysub.list"
        """                   
        ########################################################################
        # 5 - Quality assessment (FWHM, background, sky transparency, 
        # ellipticity, PSF quality)  
        ########################################################################
                            
        log.info("**** Data Quality Assessment **** (TBD)")                   

        ########################################################################
        # -- una prueba con astrowarp : no va mal, a simple vista da resultados 
        # parecidos, y en CPU tambien =, por tanto, opcion a considerar !!---
        # 6b - Computer dither offsets and coadd
        ########################################################################
        prueba = False  
        if prueba:
            if self.obs_mode!='dither' or self.red_mode=="quick":
                log.info("**** Doing Astrometric calibration and  coaddition result frame ****")
                #misc.utils.listToFile(self.m_LAST_FILES, out_dir+"/files_skysub.list")
                aw = reduce.astrowarp.AstroWarp(self.m_LAST_FILES, catalog="GSC-2.3", 
                coadded_file=output_file, config_dict=self.config_dict)
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
        
        ## -- fin prueba !!
        
        ########################################################################
        # 6 - Compute dither offsets from the first sky subtracted/filtered 
        # images using cross-correlation
        ########################################################################
        misc.utils.listToFile(self.m_LAST_FILES, out_dir+"/files_skysub.list")
        try:
            offset_mat = self.getPointingOffsets(out_dir+"/files_skysub.list", 
                                                 out_dir+'/offsets1.pap')                
        except Exception,e:
            log.error("Erron while getting pointing offsets. Cannot continue with data reduction...")
            raise e
        
        ########################################################################
        # 7 - First pass coaddition using offsets
        ########################################################################
        log.info("**** Initial coaddition of sky subtracted frames ****")
        fo = open(out_dir+'/offsets1.pap',"r")
        fs = open(out_dir+'/stack1.pap','w+')
        for line in fo:
            n_line = line.replace(".fits.objs", ".fits") 
            fs.write( n_line )
        fo.close()
        fs.close()    
        self.coaddStackImages(out_dir+'/stack1.pap', gainmap, 
                              out_dir+'/coadd1.fits','average')
    
        ########################################################################
        # End of first cycle: SINGLE REDUCTION (quick mode or extended object !) 
        ########################################################################
        if self.obs_mode!='dither' or self.red_mode=="quick":
            log.info("**** Doing Astrometric calibration of coadded result frame ****")
           
            reduce.astrowarp.doAstrometry(out_dir+'/coadd1.fits', output_file, 
                                           self.config_dict['astrometry']['catalog'], 
                                           config_dict=self.config_dict, 
                                           do_votable=True)
             
            log.info("Generated output file ==>%s", output_file)
            log.info("#########################################")
            log.info("##### End of QUICK data reduction ######")
            log.info("#########################################")
            return output_file 
        
        ########################################################################
        # 8 - Create master object mask
        ########################################################################

        log.info("************************")
        log.info(" START SECOND PASS      ")
        log.info("************************")

        log.info("**** Master object mask creation ****")
        obj_mask  = self.__createMasterObjMask(out_dir+'/coadd1.fits', 
                                               out_dir+'/masterObjMask.fits') 

        ########################################################################
        # 8.5 - Re-compute the gainmap taking into account the object mask
        ########################################################################
        # T O D O 
        
        ########################################################################
        # 9 - Second Sky subtraction (IRDR) using then OBJECT MASK
        ########################################################################
        log.info("**** Sky subtraction with 2nd object mask ****")
        # 9.1 Compound masked sky file list as input to IRDR::skyfilter()
        fs = open(out_dir+"/skylist2.pap","w+")
        i = 0
        j = 0
        for file in self.m_rawFiles:
            if self.apply_dark_flat==1 and master_flat!=None and master_dark!=None:
                line = file.replace(".fits","_D_F.fits") + " " + obj_mask + " "\
                + str(offset_mat[j][0]) + " " + str(offset_mat[j][1])
            elif self.apply_dark_flat==1 and master_flat!=None:
                line = file.replace(".fits","_F.fits") + " " + obj_mask + " "\
                + str(offset_mat[j][0]) + " " + str(offset_mat[j][1])
            elif self.apply_dark_flat==1 and master_dark!=None:
                line = file.replace(".fits","_D.fits") + " " + obj_mask + " "\
                + str(offset_mat[j][0]) + " " + str(offset_mat[j][1])
            else:
                line = file + " " + obj_mask + " " + str(offset_mat[j][0]) + \
                " " + str(offset_mat[j][1])
            fs.write(line+"\n")
            if (self.obs_mode=='dither_on_off' or 
                self.obs_mode=='dither_off_on') and i%2:
                j=j+1
            elif self.obs_mode=='dither':
                j=j+1
            i=i+1
        fs.close()
        self.m_LAST_FILES = self.skyFilter(out_dir+"/skylist2.pap", gainmap, 
                                           'mask', self.obs_mode)      
    
        ########################################################################
        # 9.1 - Remove crosstalk - (only if bright stars are present)    
        ########################################################################
        if self.config_dict['general']['remove_crosstalk']:
            log.info("**** Removing crosstalk ****")
            try:
                res = map ( reduce.dxtalk.remove_crosstalk, self.m_LAST_FILES, 
                            [None]*len(self.m_LAST_FILES), [True]*len(self.m_LAST_FILES))
                self.m_LAST_FILES = res
            except Exception,e:
                raise

        ########################################################################
        # 9.2 - Remove Cosmic Rays -    
        ########################################################################
        if self.config_dict['general']['remove_cosmic_ray']:
            try:
                res = map ( reduce.remove_cosmics.remove_cr, self.m_LAST_FILES, 
                            [None]*len(self.m_LAST_FILES), [True]*len(self.m_LAST_FILES),
                            [False]*len(self.m_LAST_FILES))
                self.m_LAST_FILES = res
            except Exception,e:
                raise e
        ########################################################################
        # 9.3 - LEMON connection - End here for LEMON processing    
        ########################################################################
        
        if self.red_mode=='lemon':
            misc.utils.listToFile(self.m_LAST_FILES, out_dir+"/files_skysub2.list")
            log.info("End of sequence LEMON-reduction. # %s # files created. ",
                     len(self.m_LAST_FILES))
            return out_dir+"/files_skysub2.list"
    
        ########################################################################
        # 9.4 - Divide by the master flat after sky subtraction ! (see notes above)
        # (the same task as above 4.2) --> HAS NO SENSE !!! only for a test ??? or a.l.a. O2k 
        ########################################################################
        if self.apply_dark_flat==2 and master_flat!=None:
            log.info("**** Applying Flat AFTER sky subtraction ****")
            res = reduce.ApplyDarkFlat(self.m_LAST_FILES, 
                                       None,  
                                       master_flat, 
                                       out_dir)
            self.m_LAST_FILES = res.apply()
            
        ########################################################################
        # 10a - Compute field distortion and final stack:
        #       1-Remove field distortion from individual images (SCAMP+SWARP)
        #       2-Coaddition of corrected field distortion images (SWARP)
        #       3-Final Astrometric calibration (SCAMP) of the coadded image
        ########################################################################
        _astrowarp = False
        if _astrowarp:
            print "astrowarp--->LAST_FILES=",self.m_LAST_FILES
            
            log.info("**** Astrometric calibration and stack of individual \
            frames to field distortion correction ****")
            aw = reduce.astrowarp.AstroWarp(self.m_LAST_FILES, 
                                            coadded_file=output_file, 
                                            config_dict=self.config_dict,
                                            do_votable=True)
            try:
                aw.run()
            except Exception,e:
                log.error("Some error while running Astrowarp....")
                raise e
            log.info("Sucessful end of Pipeline (I hope!)")
            return output_file
        
            
    
        ########################################################################
        # 10b.1 - Create second coadded image of the dithered stack using new sky 
        # subtracted frames (using the same offsets)
        ########################################################################
        log.info("**** Coadding image free distorion frames ****")
        self.coaddStackImages(out_dir+'/stack1.pap', gainmap, out_dir+'/coadd2.fits')
        reduce.imtrim.imgTrim(out_dir+'/coadd2.fits')
        
        ########################################################################
        # 10b.2 - Make Astrometry
        ########################################################################
        log.info("**** Computing Astrometric calibration of coadded (2nd) result frame ****")
        reduce.astrowarp.doAstrometry(out_dir+'/coadd2.fits', output_file, 
                                      self.config_dict['astrometry']['catalog'] , 
                                      config_dict=self.config_dict)
        
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


## {{{ http://code.activestate.com/recipes/327142/ (r1)
def printDict(aDict, br='\n', html=0,
            keyAlign='l',   sortKey=0,
            keyPrefix='',   keySuffix='',
            valuePrefix='', valueSuffix='',
            leftMargin=0,   indent=1 ):
    '''
return a string representive of aDict in the following format:
    {
     key1: value1,
     key2: value2,
     ...
     }

Spaces will be added to the keys to make them have same width.

sortKey: set to 1 if want keys sorted;
keyAlign: either 'l' or 'r', for left, right align, respectively.
keyPrefix, keySuffix, valuePrefix, valueSuffix: The prefix and
   suffix to wrap the keys or values. Good for formatting them
   for html document(for example, keyPrefix='<b>', keySuffix='</b>'). 
   Note: The keys will be padded with spaces to have them
         equally-wide. The pre- and suffix will be added OUTSIDE
         the entire width.
html: if set to 1, all spaces will be replaced with '&nbsp;', and
      the entire output will be wrapped with '<code>' and '</code>'.
br: determine the carriage return. If html, it is suggested to set
    br to '<br>'. If you want the html source code eazy to read,
    set br to '<br>\n'

version: 04b52
author : Runsun Pan
require: odict() # an ordered dict, if you want the keys sorted.
         Dave Benjamin 
         http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/161403
    '''
   
    if aDict:

        #------------------------------ sort key
        if sortKey:
            dic = aDict.copy()
            keys = dic.keys()
            keys.sort()
            ###aDict = odict()
            for k in keys:
                aDict[k] = dic[k]
            
        #------------------- wrap keys with ' ' (quotes) if str
        tmp = ['{']
        ks = [type(x)==str and "'%s'"%x or x for x in aDict.keys()]

        #------------------- wrap values with ' ' (quotes) if str
        vs = [type(x)==str and "'%s'"%x or x for x in aDict.values()] 

        maxKeyLen = max([len(str(x)) for x in ks])

        for i in range(len(ks)):

            #-------------------------- Adjust key width
            k = {1            : str(ks[i]).ljust(maxKeyLen),
                 keyAlign=='r': str(ks[i]).rjust(maxKeyLen) }[1]
            
            v = vs[i]        
            tmp.append(' '* indent+ '%s%s%s:%s%s%s,' %(
                        keyPrefix, k, keySuffix,
                        valuePrefix,v,valueSuffix))

        tmp[-1] = tmp[-1][:-1] # remove the ',' in the last item
        tmp.append('}')

        if leftMargin:
            tmp = [ ' '*leftMargin + x for x in tmp ]
          
        if html:
            return '<code>%s</code>' %br.join(tmp).replace(' ','&nbsp;')
        else:
            return br.join(tmp)     
    else:
        return '{}'

'''
Example:

>>> a={'C': 2, 'B': 1, 'E': 4, (3, 5): 0}

>>> print prnDict(a)
{
 'C'   :2,
 'B'   :1,
 'E'   :4,
 (3, 5):0
}

>>> print prnDict(a, sortKey=1)
{
 'B'   :1,
 'C'   :2,
 'E'   :4,
 (3, 5):0
}

>>> print prnDict(a, keyPrefix="<b>", keySuffix="</b>")
{
 <b>'C'   </b>:2,
 <b>'B'   </b>:1,
 <b>'E'   </b>:4,
 <b>(3, 5)</b>:0
}

>>> print prnDict(a, html=1)
<code>{
&nbsp;'C'&nbsp;&nbsp;&nbsp;:2,
&nbsp;'B'&nbsp;&nbsp;&nbsp;:1,
&nbsp;'E'&nbsp;&nbsp;&nbsp;:4,
&nbsp;(3,&nbsp;5):0
}</code>

>>> b={'car': [6, 6, 12], 'about': [15, 9, 6], 'bookKeeper': [9, 9, 15]}

>>> print prnDict(b, sortKey=1)
{
 'about'     :[15, 9, 6],
 'bookKeeper':[9, 9, 15],
 'car'       :[6, 6, 12]
}

>>> print prnDict(b, keyAlign="r")
{
        'car':[6, 6, 12],
      'about':[15, 9, 6],
 'bookKeeper':[9, 9, 15]
}
'''
## end of http://code.activestate.com/recipes/327142/ }}}


        
#################
## Important Note 
#################
#13-Jun-2012
#===========
# It it very important to create the pool of workers here, just after the 
# class definition and just before any call to PyRAF in order to avoid
# the bug found by V.Terron that seems to be in PyRAF 

# It seems that PyRAF is messing with the multiprocessing module. The reason
# why I think so is because two consecutive calls to difiphot.apphot.qphot work
# fine, as well as using it as the function passed to map_async. But calling
# qphot and then map_async, with the instantiation of the pool of workers done
# in-between, makes the scripts fail with the most arcane of errors, such as
# (yes, this is a real example):
#
# "Exception Exception OSErrorOSError: : ((1010, , ''NNoo cchhiilldd
# pprroocceesssseess'')) in in <bound method Subprocess.__del__ of <Subprocess
# '/iraf/iraf/noao/bin.linux/x_apphot.e -c', at 3c8c2d8>><bound method
# Subprocess.__del__ of <Subprocess '/iraf/iraf/noao/bin.linux/x_apphot.e -c',
# at 3c8c2d8>> ignored"
#
# However, everything goes as expected is the pool of workers is created before
# qphot is ever used, directly or indirectly, in the script. It seems the first
# execution of PyRAF's qphot is affecting how the multiprocessing module
# works. That's why we need to create the pool before PyRAF runs, so that we
# use the original, 'unmodified' code.

# 
#multiprocessing.freeze_support()
#pool = multiprocessing.Pool(2)#processes=multiprocessing.cpu_count())
#If processes is None then the number returned by cpu_count() is used.         

#15-Jun-2012
#===========
#Finally, it was not possible to create as a global variable the pool, becasue
#it there is problems with sharing a global pool (or whatever variable) with
#the Process created in the QL module for dispatching the on-line reduction.
#So, what we do now is create a new pool everytime we need to do any iraf
#depending task.
 

        
        
        
        
        
        
        
        
        
        
        
    
        
