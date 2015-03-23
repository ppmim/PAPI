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
################################################################################
    
#From system
import os
import os.path
import fileinput
import shutil
import tempfile
import dircache
import multiprocessing
import itertools
import math

# IRAF packages
import pyraf
from pyraf import iraf
#from iraf import noao

# Math module for efficient array processing
import numpy
import astropy.io.fits as fits
from astropy import wcs

# Log
import misc.paLog
from misc.paLog import log    
from misc.version import __version__


#PAPI packages 
import datahandler
import reduce
import reduce.checkQuality
import misc.fileUtils
import misc.utils
from reduce.makeobjmask import *
import misc.imtrim
import reduce.remove_cosmics
import reduce.astrowarp
import reduce.solveAstrometry
import misc.mef 
import astromatic
from astromatic.swarp import *
import datahandler.dataset
import misc.collapse
import correctNonLinearity



# If your parallel tasks are going to use the same instance of PyRAF (and thus 
# the same process cache), as in the case of running the entire parallel program 
# inside a single Python script via, say, the multiprocessing module, you will 
# want to turn off process caching. This turns off the ability for PyRAF to 
# use/re-use the same continually running, connected, IRAF task subprocess to do 
# the work of each call to the task. With no process caching allowed, each new 
# call you make to the IRAF task will start a new, separate IRAF executable, 
# which will live only as long as it is doing your work, which is what you want. 
# To turn off process caching, include this line in your Python script:
#
#           iraf.prcacheOff()
#   
#
# if you imported iraf from pyraf as above.
pyraf.iraf.prcacheOff()

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
        
        out_file: str
            Filename of the final result file (reduced seq.) produced
            If no outfile name is given, the result of each sequence reduced
            will be saved with a filename as: 'PANIC.[DATE-OBS].fits',
            where DATE-OBS is the keyword value of the first file in the sequence
        obs_mode : str
            'dither' - files belong to a dither pattern observation
            'other' - 
        
        dark : str
            if given, the filename of the master dark to be used for the reduction
        
        flat : str
            if given, the filename of the master flat to be used for the reduction
        
        bpm : str
            the Bad Pixel Mask to be used to the data reduction (to fix, to 
            to consider in gainmaps or nothing to do).

        red_mode: str
            reduction mode to run the pipeline

        group_by: str
            grouping mode of files (OT metadata based or Filter based)

        check_data: str
            whether to check the data matches Readout mode, ExpTime, NCOADDS, ...

        config_dict: dictionary
            dictionary containing the configuration parameters

        external_db_files : str or list
            File list used as an external calibration database.
            Then, if during the reduction of a ReductionSet(RS) no calibration 
            (dark, flat) are found in the current RS, then PAPI will look for 
            them into this directory.
            If the directory does not exists, or no calibration are found, then
            no calibrations will be used for the data reduction.
            Note that the calibrations into the current RS have always higher 
            priority than the ones in the external calibration DB.   
        
        temp_dir: str
            Pathname of the temporal directory
        
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
            raise ReductionSetException ("Empty file list, no files to reduce ...")
        else:
            for f in rs_filelist:
                if not os.path.exists(f) or not os.path.splitext(f)[1] == ".fits":
                    log.error("File %s does not exists or not has '.fits' extension", f)
                    raise ReductionSetException ("File %s does not exist or \
                    has not '.fits' extension"%f)

        
        # List containing the science data filenames to reduce
        self.rs_filelist = rs_filelist 

        # Main directories:

        # PAPI_HOME
        try:
            self.papi_home = os.environ['PAPI_HOME']
            if self.papi_home[-1]!='/':
                self.papi_home+='/'
        except Exception,e:
            log.error("Error, variable PAPI_HOME not defined.")
            raise e

        # Temporal directory for temporal files created during the data 
        # reduction process
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
        if out_file==None: 
            # if no filename, we take the DATE-OBS of the first file of the SEQ
            with fits.open(rs_filelist[0], ignore_missing_end=True) as myhdulist:
                if 'DATE-OBS' in myhdulist[0].header:
                    self.out_file = self.out_dir + "/PANIC." + \
                            myhdulist[0].header['DATE-OBS'] +".fits"
                else:
                    output_fd, self.out_file = tempfile.mkstemp(suffix='.fits', 
                                                        prefix='red_set_', 
                                                        dir=self.out_dir)
                    os.close(output_fd)
                    os.unlink(self.out_file) # we only need the name
        else:    
            self.out_file = out_file  # final reduced data file (out)
            
        
        self.obs_mode = obs_mode  # observing mode (dither, dither_on_off, dither_off_on....)
        self.master_dark = dark # master dark to use (input)
        self.master_flat = flat # master flat to use (input)
        self.apply_dark_flat = self.config_dict['general']['apply_dark_flat'] # 0=no, 1=before, 2=after
        self.coadd_mode = self.config_dict['general']['coadd_mode']
        
        #
        # Bad Pixel Maks
        # 
        self.bpm_mode = self.config_dict['bpm']['mode']
        if self.bpm_mode != 'none':
            self.master_bpm  = self.config_dict['bpm']['bpm_file']
        else:
            self.master_bpm  = bpm  
        
        # 
        # Non-Linearity correction
        # 
        self.non_linearity_apply = self.config_dict['nonlinearity']['apply']
        if self.non_linearity_apply:
            myhdulist = fits.open(rs_filelist[-1])
            if myhdulist[0].header['READMODE'] == 'line.interlaced.read':
                self.non_linearity_model = self.config_dict['nonlinearity']['model_lir']
            elif myhdulist[0].header['READMODE'] == 'fast-reset-read.read':
                self.non_linearity_model = self.config_dict['nonlinearity']['model_rrrmpia']
            else:
                log.warning("Non-Linearity model does not match. Correction de-activated.")
                self.non_linearity_apply = False
                self.non_linearity_model = None

        #
        # Reduction mode (lemon=for LEMON pipeline, quick=for QL, science=for 
        # science, lab=laboratory) 
        #
        self.red_mode = red_mode
        
        # Flag to decide how classification will be done, by OT (OB_ID, OB_PAT,
        # FILTER and TEXP kws) or by FILTER (FILTER kw)
        self.group_by = group_by.lower() 

        # Flag to indicate if data checking need to be done (see checkData() 
        # method)
        self.check_data = check_data 
                
        if self.config_dict:
            self.m_irdr_path = self.papi_home + self.config_dict['config_files']['irdr_bin']
            # update config values from config_dict
            # half width of sky filter window in frames
            self.HWIDTH = self.config_dict['skysub']['hwidth'] 
            # Minimum number of sky frames required in the sliding window for 
            # the sky subtraction
            self.MIN_SKY_FRAMES = self.config_dict['skysub']['min_frames']
            # Maximum exposition time difference (days) between two consecutive 
            # frames  
            self.MAX_MJD_DIFF = self.config_dict['general']['max_mjd_diff'] /86400.0 
            # Minimum overlap correlation fraction between offset translated 
            # images (from irdr::offset.c)
            self.MIN_CORR_FRAC = self.config_dict['offsets']['min_corr_frac'] 
        else:
            print "Program should not enter here  !!!"
            # Some "default" config values (see below how they are updated from 
            # the config_dict)
            # Environment variables
            """
            self.m_irdr_path = os.environ['IRDR_BIN']
            self.MAX_MJD_DIFF = (1/86400.0)*10*60 #6.95e-3  # Maximum seconds (10min=600secs aprox) of temporal distant allowed between two consecutive frames 
            self.HWIDTH = 2 #half width of sky filter window in frames
            self.MIN_SKY_FRAMES = 5  # minimun number of sky frames required in the sliding window for the sky subtraction
            self.MIN_CORR_FRAC = 0.1 # Minimun overlap correlation fraction between offset translated images (from irdr::offset.c)
            """
        
        # Variables holding current reduction status
        
        # (m_LAST_FILES) Contain the files as result of the last processing 
        # step (science processed frames).Properly initialized in 
        # reduceSingleObj()
        self.m_LAST_FILES = []   

        # Raw files (originals in the working directory) of the current 
        # sequence being reduced.
        self.m_rawFiles = []     
        # Filter of the current data set (m_LAST_FILES)
        self.m_filter = ""
        # Type (dark, flat, object, ...) of the current data set; should be 
        # always object !       
        self.m_type = ""
        # Exposition Time of the current data set files         
        self.m_expt = 0.0
        # Number of coadds of the current data set files        
        self.m_ncoadd = 0  
        # Integration Time of the current data set files      
        self.m_itime  = 0.0 
        # readout mode (must be the same for all data files)     
        self.m_readmode = ""     
        
        
        ## Local DataBase (in memory)
        self.db = None
        
        ## External DataBase (in memory) - optional
        # It is an optional external database provided when creating a ReductionSet 
        # instance and mainly can have calibration files required for the reduction 
        # of the current data set. It will be used during the on-line QL 
        # data reduction,e.g., tw_flats that require a master dark than can be 
        # in other RS. Mainly used in Quick-Look !!    
        self.ext_db = None
        
        # (optional) file list to build the external DB. We proceed this way,
        # because if we give a DB connection to the ReductionSet class instead
        # of a list of files, we can have problems because SQLite3 does not
        # support access from multiple  threads, and the RS.reduceSet() can
        # be executed from other thread than it was created.
        # --> It can be a <list of files> or a <directory name> having the files
        # Further info: http://stackoverflow.com/questions/393554/python-sqlite3-and-concurrency  
        #               http://www.sqlite.org/cvstrac/wiki?p=MultiThreading
        self.ext_db_files = []
        if external_db_files == None:
            if self.config_dict !=None:
                cal_dir = self.config_dict['general']['ext_calibration_db']
                if os.path.isdir(cal_dir):
                    for ifile in dircache.listdir(cal_dir):
                        if ifile.endswith(".fits") or ifile.endswith(".fit"):
                            self.ext_db_files.append((cal_dir + "/" + ifile).replace('//','/'))
        else:
            self.ext_db_files = external_db_files


        # Print config_dictonary values in log file for debugging
        log.info("[ReductionSet] CONFIGURATION VALUES OF RS ------------------")
        log.info("PAPI Version: %s" %__version__)
        log.info(printDict(self.config_dict))
        log.info("------------------------------------------------------------")
  
    def __initDB(self):
        """
        Initialize the Data Bases (local and external), loading the full frame 
        list.
        
        Notes
        ----- 
        If we initialize the DB in the constructor __init__(), we will have 
        problems if we use the DB from any method of ReductionSet and then from 
        other thread, i.e.,sqlite3 do not support access from multiple threads.
        For more info, see :
        
        http://stackoverflow.com/questions/393554/python-sqlite3-and-concurrency
        http://www.sqlite.org/cvstrac/wiki?p=MultiThreading
        
        """
        
        # Local DataBase (in memory)
        log.info("Initializing Local DataBase")
        instrument = self.config_dict['general']['instrument'].lower()
        try:
            self.db = datahandler.dataset.DataSet(self.rs_filelist, instrument)
            self.db.createDB()
            self.db.load()
        except Exception,e:
            log.error("Error while LOCAL data base initialization: \n %s"%str(e))
            raise Exception("Error while LOCAL data base initialization")
            
        
        self.m_LAST_FILES = self.rs_filelist # Later properly initialized in reduceSingleObj()
        
        if self.master_dark!=None: self.db.insert(self.master_dark)
        if self.master_flat!=None: self.db.insert(self.master_flat)
        if self.master_bpm !=None: self.db.insert(self.master_bpm)
        
        # self.db.ListDataSet()
        
        # External DataBase (in memory)
        if len(self.ext_db_files) > 0:
            log.info("Initializing External DataBase")
            try:
                self.ext_db = datahandler.dataset.DataSet(self.ext_db_files,
                                                          instrument)
                self.ext_db.createDB()
                self.ext_db.load()
                log.info("Calibration files found in External DB:")
                self.ext_db.ListDataSet()
            except Exception,e:
                log.error("Error while EXTERNAL data base initialization: \n %s"%str(e))
                raise Exception("Error while EXTERNAL data base initialization")

            # self.ext_db.ListDataSet()
        else:
            self.ext_db = None
                
        
    def checkData(self, chk_shape=True, chk_filter=True, chk_type=True, 
                  chk_expt=True, chk_itime=True, chk_ncoadd=True, chk_cont=True,
                  chk_readmode=True, chk_instrument=True, file_list=None):
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

        chk_instrument: bool
            Flag meaning if INSTRUME keyword must be checked to be equal.
            In addition, INSTRUME keyword must match with value specified
            in the configuration file.

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
        instrument_0 = f.getInstrument().lower()
        
        log.info("Values to check: FILTER=%s TYPE=%s EXPT=%s ITIME=%s NCOADD=%s READMODE=%s SHAPE=%s INST=%s",
                 filter_0, type_0, expt_0, itime_0, ncoadd_0, readmode_0, shape_0, instrument_0)
        
        self.m_filter = filter_0
        self.m_type = type_0
        self.m_expt = expt_0
        self.m_itime = itime_0
        self.m_ncoadd = ncoadd_0
        self.m_readmode = readmode_0
        self.m_instrument = instrument_0
        
        if (chk_instrument and 
            self.config_dict['general']['instrument'].lower()!=instrument_0):
            log.error("INSTRUMENT value mismatch -- %s "%instrument_0)
            return (False, "chk_instrument") 
        
        mismatch_filter = False
        mismatch_type = False
        mismatch_expt = False
        mismatch_itime = False
        mismatch_ncoadd = False
        mismatch_cont = False
        mismatch_readmode = False
        mismatch_shape = False
        mismatch_instrument = False
        
        prev_MJD = -1
        for file in files_to_check:
            fi = datahandler.ClFits( file )
            if chk_instrument and not mismatch_instrument: 
                if fi.getInstrument() != instrument_0:
                    log.debug("File %s does not match INSTRUMENT", file)
                    mismatch_instrument = True
                    info_mismatch = "chk_instrument"
                    break
            if chk_shape and not mismatch_shape: 
                if fi.shape != shape_0:
                    log.debug("File %s does not match SHAPE", file)
                    mismatch_shape = True
                    info_mismatch = "chk_shape"
                    break
            if chk_filter and not mismatch_filter: 
                if fi.getFilter() != filter_0:
                    log.debug("File %s does not match FILTER", file)
                    mismatch_filter = True
                    info_mismatch = "chk_filter"
                    break
            if chk_type and not mismatch_type: 
                if fi.getType() != type_0:
                    log.debug("File %s does not match TYPE", file)
                    mismatch_type = True
                    info_mismatch = "chk_type"
                    break
            if chk_expt and not mismatch_expt: 
                if fi.expTime() != expt_0:
                #if prev_MJD!=-1 and ((fi.expTime()+self.MAX_MJD_DIFF)<expt_0 or
                #    (fi.expTime()-self.MAX_MJD_DIFF)>expt_0):   # more relaxed situation
                    log.debug("File %s does not match EXPTIME", file)
                    mismatch_expt = True
                    info_mismatch = "chk_expt"
                    break
            if chk_itime and not mismatch_itime: 
                if fi.getItime() != itime_0:
                    log.debug("File %s does not match ITIME", file)
                    mismatch_itime = True
                    info_mismatch = "chk_itime"
                    break
            if chk_ncoadd and not mismatch_ncoadd: 
                if fi.getNcoadds() != ncoadd_0:
                    log.debug("File %s does not match NCOADD", file)
                    mismatch_ncoadd = True
                    info_mismatch = "chk_ncoadd"
                    break
            if chk_readmode and not mismatch_readmode: 
                if fi.getReadMode() != readmode_0:
                    log.debug("File %s does not match READMODE", file)
                    mismatch_readmode = True
                    info_mismatch = "chk_readmode"
                    break
            if chk_cont and not mismatch_cont:
                if prev_MJD!=-1 and (fi.getMJD()-prev_MJD)>self.MAX_MJD_DIFF:
                    log.error("Maximmun time distant between two consecutives frames exceeded. File= %s",file)
                    mismatch_cont = True
                    info_mismatch = "chk_cont"
                    break
                else:
                    prev_MJD=fi.getMJD()
                    
        if mismatch_shape or mismatch_filter or mismatch_type or mismatch_expt or  \
            mismatch_itime or mismatch_ncoadd or mismatch_cont:
            log.error("Data checking found a mismatch....check your data files....")
            #raise Exception("Error while checking data (filter, type, ExpT, Itime, NCOADDs, MJD)")
            return (False, info_mismatch)             
        else:    
            log.debug("All files match same file filter")
            return (True, None)
            
    def checkFilter(self):
        """
        Return true is all files in file have the same filter type, otherwise False.
        """
        f = datahandler.ClFits(self.m_LAST_FILES[0])
        filter_0 = f.getFilter()
        self.m_filter = filter_0
        for file in self.m_LAST_FILES:
            fi = datahandler.ClFits( file )
            if fi.getFilter() != filter_0:
                log.debug("File %s does not match file filter", file)
                return False
            
        log.debug("All files match same file filter")
        return True
        
        
    def checkType(self, type_to_check=None):
        """
        Return true is all files in file have the same type(science, dark, flat,
        ...), false otherwise
        
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
        
        new_frame_list = [] # a list of N list, where N=number of extension of the MEF 
        nExt = 0
        if frame_list==None or len(frame_list)==0 or frame_list[0]==None:
            log.debug("Nothing to split ! %s"%frame_list)
            return [],0
        
        first_img = datahandler.ClFits(frame_list[0])
        # First, we need to check if we have MEF files
        # 1) Not a MEF 
        if not first_img.isMEF():
            if first_img.isPANICFullFrame():
                try:
                    mef = misc.mef.MEF(frame_list)
                    (nExt, sp_frame_list) = mef.splitGEIRSToSimple(".Q%02d.fits", 
                                                                   out_dir=self.temp_dir)
                except Exception,e:
                    log.error("Some error while splitting PANIC data set. %s",str(e))
                    raise e   
            else:
                # No split is required
                log.debug("No split is required")
                nExt = 1
                new_frame_list.append(frame_list)
        # 2) Suppose we have MEF files ...
        else:
            if first_img.getInstrument() == "hawki":
                kws_to_cp = ['DATE','OBJECT','DATE-OBS','RA','DEC','RADECSYS','UTC','LST',
		          'UT','ST','AIRMASS','IMAGETYP','EXPTIME','TELESCOP','INSTRUME','MJD-OBS',
		          'FILTER', 'FILTER1', 'FILTER2', "HIERARCH ESO TPL ID", "HIERARCH ESO TPL EXPNO", 
                  'HIERARCH ESO TPL NEXP','NCOADDS','HIERARCH ESO DET NDIT', 'NDIT',
                  'HIERARCH ESO INS FILT1 NAME', 'HIERARCH ESO INS FILT2 NAME',
                  'BSCALE', 'BZERO',
                ]
                instr = 'hawki'  
            else: # PANIC
                kws_to_cp = ['DATE','OBJECT','DATE-OBS','RA','DEC','EQUINOX','LST',
                   'UT','AIRMASS','IMAGETYP','TELESCOP','INSTRUME','MJD-OBS',
                   'BSCALE', 'BZERO',
                   'CTIME','ITIME','NCOADDS','EXPTIME','T_FOCUS','READMODE',
                   'FILTER', 'OBS_TOOL', 'PROG_ID', 'OB_ID', 
                   'OB_NAME', 'OB_PAT', 'PAT_NAME','PAT_EXPN', 'PAT_NEXP',
                   'CASSPOS','PIXSCALE', 'LAMP', 'DET_ID',
                   'PAPITYPE', 'PAPIVERS','OBSERVER', 'ORIGIN'
                ]
                instr = 'panic'
            try:
                mef = misc.mef.MEF(frame_list)
                (nExt, sp_frame_list) = mef.doSplit(".Q%02d.fits", 
                                                    out_dir=self.temp_dir, 
                                                    copy_keyword=kws_to_cp,
                                                    instrument=instr)
            except Exception,e:
                log.debug("Some error while splitting data set. %s",str(e))
                raise e
            
        # now, generate the new output filenames        
        # In principle, it is not needed; we could use [sp_frame_list]
        for n in range(1,nExt+1):
            new_frame_list.append([self.temp_dir + "/" + 
                                   os.path.basename(file.replace(".fits", ".Q%02d.fits"%n)) 
                                   for file in frame_list])
            """
            for f in new_file_names:
                #if re.search(".*(\.Q01)(.fits)$", f):
                    sources.append(f)
            """
        #print "SP_FRAME_LIST=",sp_frame_list
        #print "NEW_FRAME_LIST=",new_frame_list

        return new_frame_list, nExt
    
    def getObsMode(self):
        """
        Return the type of dither sequece followed in the currect  
        'm_LAST_FILES' list. It could be:
            - dither (T-T-T-T-T- ....)
            - dither_on_off (T-S-T-S-T-S-T-S-T-....)
            - dither_off_on (S-T-S-T-S-T-S-T-S-....)
            - other  (non defined sequence,unknown) (e.g,  T-S-T-T-S-T-T-S-T-....)
            
        NOTE: To find out the dither/observing sequence then OBJECT keyword 
              will be checked (see ClFits class)
        """
                   
        mode = 'other' # default
                   
        #Step 1: we suppose list file is sorted out  by MJD
        #Step 2: get the data type of the first file and check sequence starting from this file
        fits_0 = datahandler.ClFits(self.m_LAST_FILES[0])
        fits_1 = datahandler.ClFits(self.m_LAST_FILES[1])
        i = 0
        if fits_0.isSky() and fits_1.isObject():
            mode = 'dither_off_on'
            # Then, we are going to suppose the sequence S-T-S-T-S- .... (dither_off_on)
            for file in self.m_LAST_FILES:
                if not i%2: # even (par)
                    fits = datahandler.ClFits(file)
                    if not fits.isSky():
                        return 'other'
                elif i%2: # odd
                    fits = datahandler.ClFits(file)
                    if not fits.isObject():
                        return 'other'
                i = i + 1         
        elif fits_0.isObject() and fits_1.isSky():
            # Then, we are going to suppose the sequence T-S-T-S-T- .... (dither_on_off)
            mode = 'dither_on_off'
            for file in self.m_LAST_FILES:
                if not i%2: # even (par)
                    fits = datahandler.ClFits(file)
                    if not fits.isObject():
                        return 'other'
                elif i%2: # odd
                    fits = datahandler.ClFits(file)
                    if not fits.isSky():
                        return 'other'
                i = i + 1
        elif fits_0.isObject() and fits_1.isObject():
            # check if all are objects ...
            # and if are, then, we are going to suppose the sequence T-T-T-T-T- .... (dither)
            mode = 'dither'
            for file in self.m_LAST_FILES:
                fits = datahandler.ClFits(file)
                if not fits.isObject():
                    return 'other'
        else:
            # it might be a lab test or whatever ...
            mode = 'other'
                        
        return mode
                        
    def getSkyFrames(self, list=None):
        """
        Given a list of files(data set), return the files identified as 
        'sky' frames in the m_LAST_FILES
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
        Given a list of frames belonging to a observing sequence for 
        a given object (star, galaxy, whatever),return the most recently created 
        calibration files (master dark,flat,bpm) in order to reduce the sequence.
        The search of the calibration files is done, firstly in the local DB, 
        but if no results, then in the external DB if it was provided.
        
        This routine is also used (copy) on the QL (mainGUI.py)

        Returns
        -------  
        A triplet with the calibration files (dark, flat, bpm) found; If more 
        than one master were found, the most recently created (according to MJD) 
        is returned.
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
        
        # Init the DB (ext_db is always initialited whether db!=None)  
        if self.db==None: 
            self.__initDB()

        # DARK - Does require equal EXPTIME Master Dark ???
        # First, look for a DARK_MODEL
        master_dark = self.db.GetFilesT('MASTER_DARK_MODEL', -1) 
        if len(master_dark)==0 and self.ext_db!=None:
            master_dark = self.ext_db.GetFilesT('MASTER_DARK_MODEL', -1)
        
        # Secondly (hopefully), try to find a MASTER_DARK with equal expTime
        if len(master_dark)==0:
            log.info("Now, trying to find a MASTER_DARK")
            master_dark = self.db.GetFilesT('MASTER_DARK', expTime)
        if len(master_dark)==0 and self.ext_db!=None:
            log.info("Last chance to find a MASTER_DARK in ext_DB")
            master_dark = self.ext_db.GetFilesT('MASTER_DARK', expTime)
        
             
        # FLATS - Do NOT require equal EXPTIME, but FILTER
        master_flat = self.db.GetFilesT('MASTER_DOME_FLAT', -1, filter)
        if master_flat==[]:
            master_flat = self.db.GetFilesT('MASTER_TW_FLAT', -1, filter)
        if len(master_flat)==0 and self.ext_db!=None:
            master_flat = self.ext_db.GetFilesT('MASTER_DOME_FLAT', -1, filter)
            if len(master_flat)==0:
                master_flat=self.ext_db.GetFilesT('MASTER_TW_FLAT', -1, filter)

        # BPM                
        master_bpm = self.db.GetFilesT('MASTER_BPM')
        if len(master_bpm)==0 and self.ext_db!=None:
            master_bpm = self.ext_db.GetFilesT('MASTER_BPM')

        log.debug("Master Darks found %s", master_dark)
        log.debug("Master Flats found %s", master_flat)
        log.debug("Master BPMs  found %s", master_bpm)
        
        
        # Return the most recently created (according to MJD order)
        if len(master_dark)>0: 
            r_dark = master_dark[-1]
            log.debug("First DARK candidate: %s"%r_dark)            
            r_dark = self.getBestShapedFrame(master_dark, sci_obj_list[0])
            log.debug("Second DARK candidate: %s"%r_dark)            
        else: 
            r_dark = None
        
        if len(master_flat)>0:
            r_flat = master_flat[-1]
            log.debug("First FLAT candidate: %s"%r_flat)            
            r_flat = self.getBestShapedFrame(master_flat, sci_obj_list[0])
            log.debug("Second FLAT candidate: %s"%r_flat)            
        else: 
            r_flat = None
        
        if len(master_bpm)>0: 
            r_bpm = master_bpm[-1]
            log.debug("First BPM candidate: %s"%r_bpm)            
            r_bpm = self.getBestShapedFrame(master_bpm, sci_obj_list[0])
            log.debug("Second BPM candidate: %s"%r_bpm)            
        else: 
            r_bpm = None

        return r_dark, r_flat, r_bpm
        
    def getBestShapedFrame(self, framelist, src_frame):
        """
        Given a list of frames (calibrations) sorted by MJD, return the frame that
        has the same number of extension (MEF) than src_frame and with the 
        same shape (dimensions).

        This routine is also used (copy) on the QL (mainGUI.py)

        Returns
        -------
        Returns the calibration file that match the src_frame, otherwise, None
        is returned.
        
        """

        candidate = None
        with fits.open(src_frame) as src_fd:
            for iframe in framelist:
                with fits.open(iframe) as ifd:
                    # not MEF
                    if len(ifd)==len(src_fd) and len(ifd)==1:
                        if ifd[0].data.shape==src_fd[0].shape:
                            return iframe
                        elif datahandler.ClFits(iframe).isMasterDarkModel():
                            # Exception,  DarkModel will have always 2 layers 
                            # per extension.
                            return iframe
                        else: 
                            continue
                    # MEF 
                    elif len(ifd)==len(src_fd) and len(ifd)>1:
                        # We should check each extension, but it would be strange
                        # to have extension with different shapes.
                        if ifd[1].data.shape==src_fd[1].shape:
                            return iframe
                        elif datahandler.ClFits(iframe).isMasterDarkModel():
                            # Exception, DarkModel will have always 2 layers 
                            # per extension
                            return iframe
                        else:
                            # dimesions do not fit
                            continue
                    # diferent number of extensions
                    else:
                        continue

            # If this point is reached, it means that any suitable file was found.
            return None


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
        grouping criteria given by self.group_by and orderted by MJD.
        
        Parameters
        ----------
        show: bool
            if True, print out in std output the found sequences.
        
        Returns
        -------
        A list of lists of sequence files and their types (DARK, TW_FLAT, 
        DOME_FLAT, SCIENCE) ; see ClFits class for further details).
        
        """
        
        log.debug("[getSequences] Looking for Data Sequences into the DataSet")
        seqs = []
        seq_types = []
        
        if self.group_by=='ot':
            seqs, seq_types = self.getOTSequences(show)
        elif self.group_by=='filter':
            if self.db==None: self.__initDB()
            seqs, seq_types = self.db.GetSequences(group_by='filter',
               max_mjd_diff = self.config_dict['general']['max_mjd_diff']/86400.0,
               max_ra_dec_diff = self.config_dict['general']['max_ra_dec_offset'],
               max_nfiles = self.config_dict['general']['max_num_files'])
        else:
            log.error("[getSequences] Wrong data grouping criteria")
            raise Exception("[getSequences] Found a not valid data grouping criteria %s"%(self.group_by))

        # Print out the groups. It can be used to copy&paste the sequence files
        # because is a "clean" list of full path-names of files, instead of
        # files printed in "debug" output 
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
                    #print file + " type= %s"%self.db.GetFileInfo(file)[2]
                    print file 
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
        Look for sequences (calib, science) in the current data set and print
        the results.
        The sequences must follow the PANIC Observing Tool schema, and basically 
        are detected by the pair {PAT_EXPN, PAT_NEXP} :
        
        OBS_TOOL, OB_ID, OB_PAT, PAT_EXPN, PAT_NEXP
         
        
        The algorith followed is the next:
        
            1. Sort out the data files by MJD
            2. Do
                2.1 Look for the first file with OBS_TOOL=True and PAT_EXPN=1
                    2.1.1 Group all next files up to PAT_EXPN=PAT_NEXP
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
            # return a list of SCIENCE (or SKY) frames grouped by FILTER
            (seq_list, seq_par) = self.db.GetSequences(group_by='filter',
                                                       max_mjd_diff=self.config_dict['general']['max_mjd_diff']/86400.0,
                                                       max_ra_dec_diff=self.config_dict['general']['max_ra_dec_offset'],
                                                       max_nfiles=self.config_dict['general']['max_num_files'])
                
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
        """
        Check all the current files in the ReductionSet list to find out if 
        they are all calibatrion frames"""
        
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
        For 'each' (really not each, depend on dither pattern, e.g., 
        extended sources) input image, a sky frame is computed by combining a 
        certain number of the closest images, then this sky frame is subtracted
        to the image and the result is divided by the master flat; 
        
        
        Parameters
        ----------
        list_file: str
            a text file containing the suited structure in 
            function of the observing_mode (the list shall be sorted by obs-date)

        gain_file: str
            gain map file used as bad pixel mask
        
        mask: str
            [mask|nomask] flag used to indicate if a object mask was 
            especified into the 'list_file'
        
        obs_mode: str
            [dither|other] dither or other(e.g., nodding) pattern 
            to process
          
        skymodel: str
            [median|min] sky model used for sky subtraction. Only 
            required if obs_mode=dither. (median=coarse fields, min=crowded fields)
           
        Returns
        -------
        The function generate a set of sky subtrated images (*.skysub.fits) and
        Return ONLY filtered images; when extended-source-mode ,sky 
        frames are not included in the returned file list. 
        The out-file-list is previously ordered by obs-data.             

        Notes
        -----                 
        - This function is a wrapper for skyfilter.c (IRDR). 
        - The detection of SKY/TARGET frames is done in 'skyfilter_general' based
        on next header keywords: OBJECT=SKY or IMAGETYP=SKY.

        """               
        
        # Skyfilter parameters
        halfnsky = self.HWIDTH  # value from config file [skysub.hwidth]
        destripe = 'none'
        out_files = []
        
        # get the skymodel
        if skymodel == None:
            skymodel = self.config_dict['skysub']['skymodel']
            
        if obs_mode == 'dither':
            # It comprises any dithering pattern for point-like objects
            skyfilter_cmd = self.m_irdr_path + '/skyfilter '+ list_file + \
            '  ' + gain_file + ' ' + str(halfnsky) + ' ' + mask + '  ' +  \
            destripe + '  ' + skymodel
        elif obs_mode == 'other':
            # It comprises nodding pattern for extended objects 
            # The detection of SKY/TARGET frames is done in 'skyfilter_general' based
            # on next header keywords: OBJECT=SKY or IMAGETYP=SKY.
            skyfilter_cmd = self.m_irdr_path + '/skyfilter_general ' + list_file \
            + '  ' + gain_file +' '+ str(halfnsky)+' '+ mask + '  ' + destripe
        elif obs_mode == 'dither_on_off': 
            # actually not used
            skyfilter_cmd = self.m_irdr_path + '/skyfilteronoff ' + list_file + \
            '  ' + gain_file + ' ' + str(halfnsky) + ' ' + mask + '  ' + destripe 
        elif obs_mode == 'dither_off_on':
            # actually not used 
            skyfilter_cmd = self.m_irdr_path + '/skyfilteroffon ' + list_file +\
            '  ' + gain_file + ' ' + str(halfnsky) + ' ' + mask + '  ' + destripe
        else:
            log.error("Observing mode not supported")
            raise
                  
        
        print "SKY_FILTER_CMD = ", skyfilter_cmd
        print "CMD_ARGS = ", skyfilter_cmd.split()
        args = skyfilter_cmd.split()
        
        """
        args = [self.m_irdr_path + '/skyfilter',
                list_file,
                gain_file,
                str(halfnsky), 
                mask,
                destripe
                ,skymodel]
        """
        
        output_lines = []
        try:
            output_lines  = subprocess.check_output(args, stderr = subprocess.STDOUT, 
                                                    shell=False, bufsize=0)
            log.debug("irdr::skyfilter command executed ...")
            sys.stdout.flush()
            sys.stderr.flush()
        except subprocess.CalledProcessError, e:
            log.critical("Error running skyfilter: %s"%output_lines)
            log.debug("Exception: %s"%str(e))
            raise Exception("Error running skyfilter")
        
        
        # look for sky subtracted images created by irdr::skyfilter            
        files = [line.split(" ")[0].replace("\n","") 
                    for line in fileinput.input(list_file)] # it takes into account the two kind of possible inputs files to skyfilter
        for file in files:
            # it takes into acount dither_on_off and other extended obs. patterns
            if os.path.exists(file + ".skysub"): 
                shutil.move(file.replace(".fits", ".fits.skysub"), 
                                file.replace(".fits", ".skysub.fits"))
                out_files.append(file.replace(".fits", ".skysub.fits"))
            
        # Sort-out data files by obs-data (actually, not required when obs_mode='dither')
        out_files = self.sortOutData(out_files) 
            
        return out_files
                                  
    
    def subtractNearSky(self, near_list=None, file_pos=0, out_filename=None):
        """
        Compute and subtract the nearest sky to the image in position 'fn' 
        in the given frame list (near_list). 
                     
        This function make use of skyfilter_single.c (IRDR)              
        
        INPUT
            file_pos : file position (1-N) in sci file list (the list is supposed 
            to be sorted by obs-date).
            (0 means all will be filtered)
        
        OUTPUT
            The function generate a sky subtrated image (*.skysub.fits)
        
        VERSION
            1.1, 20101103 by jmiguel@iaa.es
            1.2, 20110297 added out_dir parameter to skyfilter_single
    
        TODO: extended objects !!!!
        """
        
        log.debug("Start subtractNearSky")
        
        # default values
        if near_list==None:
            near_list = self.rs_filelist
        if out_filename==None: 
            out_filename = self.out_file
        
        
        # Some previous checks
        if file_pos < 0 or file_pos > len(near_list):
            log.error("Wrong frame number selected in near-sky subtraction")
            return None
        
        if len(near_list) < (self.HWIDTH + 1):
            log.error("Wrong number of sky frames provided. Min number of sky frame is %d", 
                      self.MIN_SKY_FRAMES)
            return None

        # 0.0 Check and collapse if required (cube images)
        near_list = misc.collapse.collapse(near_list, out_dir=self.temp_dir)
        
        # 0.1 Get the gain map
        if not self.master_flat or not os.path.exists( self.master_flat ):
            #raise Exception("Error, gain map file <%s> not found"%gain)
            # TODO: --> DONE try to compute GainMap using the given images !!!
            log.debug("---> creating gain map <----")
            output_fd, l_gainMap = tempfile.mkstemp(suffix='.fits', 
                                                    dir=self.out_dir)
            os.close(output_fd)
            os.unlink(l_gainMap) # we only need the name
            output_fd, files_list = tempfile.mkstemp(suffix='.list', 
                                                     dir=self.out_dir)
            os.close(output_fd)
            os.unlink(files_list) # we only need the name
            try:
                misc.utils.listToFile(near_list, files_list)
                # Note: we must normalize wrt chip 1, so norm=True
                superflat = reduce.SuperSkyFlat(files_list, l_gainMap, 
                                                bpm=None, norm=True, 
                                                temp_dir=self.temp_dir)
                superflat.create()
            except Exception,e:
                log.error("Error while creating gain map : %s", str(e))
                raise
        else: l_gainMap = self.master_flat

        
        # 0.2 Check if GainMap need to be split
        gain_ext, g_next = self.split([l_gainMap])
        if g_next == 1: gain_ext*=4
        
        # 1. Split MEF file (into the self.out_dir)
        obj_ext, next = self.split(near_list)
        # it must return a list of list (one per each extension) 
        out_ext = []
        
        # 2. Process each extension
        for n in range(next):
            log.debug("===> Processing extension %d", n+1)
            # Create the temp list file of nearest (ar,dec,mjd) from current 
            # selected science file
            listfile = self.out_dir + "/nearfiles.list"
            misc.utils.listToFile(obj_ext[n], listfile)
            print "NEAR_FILES=", obj_ext[n]
            print "GAIN_EXT_N=", gain_ext[n][0]
            # Call external app skyfilter (irdr)
            hwidth = self.HWIDTH
            
            cmd = self.m_irdr_path + "/skyfilter_single %s %s %d nomask none %d %s"\
                        %(listfile, gain_ext[n][0], hwidth, file_pos, self.out_dir)
            print "CMD=",cmd
            e = misc.utils.runCmd( cmd )
            if e==1: # success
                fname = self.out_dir + "/" + os.path.basename(obj_ext[n][file_pos-1].replace(".fits", (".fits.skysub")))
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
    
    def getWCSPointingOffsets(self, images_in,
                              p_offsets_file="/tmp/offsets.pap"):
      """
      Derive pointing offsets of each image taking as reference the first one and
      using WCS astrometric calibration of the images. Note that is **very** 
      important that the images are astrometrically calibrated.
      
        Parameters
        ----------
                
        images_in: str
            list of filename of images astrometrically calibrated
        p_offsets_file: str
            filename of file where offset will be saved (it can be used
            later by dithercubemean). The format of the file will be:
            
            /path/to/filename00   offsetX00   offsetY00
            /path/to/filename01   offsetX01   offsetY01
            ...
            /path/to/filename0N   offsetX0N   offsetY0N
            
        Returns
        -------    
        offsets: narray          
                two dimensional array (Nx2) with offsets (in pixles),
                where N=number of files.
                
        Notes:
            It assumed that the input images have a good enough astrometric
            calibration (hopefully obtained with Astrometry.net), and North
            is up and East is left, and no rotation angle exists.
            
        
      """
      
      log.info("Starting getWCSPointingOffsets....")

      
      # Very important the pixel scale in order to find out good offsets values !!
      pix_scale = self.config_dict['general']['pix_scale']
      # Init variables
      i = 0 
      offsets_mat = None
      ra0 = -1
      dec0 = -1
      offsets = numpy.zeros([len(images_in), 2] , dtype=numpy.float32)
      
      # Reference image
      ref_image = images_in[0]
      try:
            ref = datahandler.ClFits(ref_image)
            # If present, pix_scale in header is prefered
            pix_scale = ref.pixScale
            ra0 = ref.ra    
            dec0 = ref.dec
            log.debug("Ref. image: %s RA0= %s DEC0= %s PIXSCALE= %f"%(ref_image, ra0, dec0, pix_scale))
      except Exception,e:
          log.error("Cannot get reference RA_0,Dec_0 coordinates: %s"%str(e))
          raise e
        
      offset_txt_file = open(p_offsets_file, "w")
      for my_image in images_in:
            try:
                ref = datahandler.ClFits(my_image)
                ra = ref.ra
                dec = ref.dec
                log.debug("Image: %s RA[%d]= %s DEC[%d]= %s"%(my_image, i, ra, i, dec))
                
                # Assummed that North is up and East is left
                offsets[i][0] = ((ra - ra0)*3600 * math.cos(dec/57.29578)) / float(pix_scale)
                offsets[i][1] = ((dec0 - dec)*3600) / float(pix_scale)
                
                log.debug("offset_ra  = %s"%offsets[i][0])
                log.debug("offset_dec = %s"%offsets[i][1])
                
                offset_txt_file.write(my_image + "   " + "%.6f   %0.6f\n"%(offsets[i][0], offsets[i][1]))
                i+=1
            except Exception,e:
                log.error("Error computing the offsets for image %s. \n %s"%(my_image, str(e)))
                raise e
        
      offset_txt_file.close()
      
      # Write out offsets to file
      # numpy.savetxt(p_offsets_file, offsets, fmt='%.6f')
      log.debug("(WCS) Image Offsets (pixels): ")
      numpy.set_printoptions(suppress=True)
      log.debug(offsets)
      
      return offsets
    
    def getPointingOffsets (self, images_in=None, 
                            p_offsets_file='/tmp/offsets.pap'):
        """
        Derive pointing offsets between each image using SExtractor OBJECTS 
        (makeObjMask) and offsets (IRDR).
        
        Note: (from Infrared Imaging Data Reduction Software and Techniques, C.N.Sabbey)
        The approximate dither offsets stored inthe FITS header WCS information are 
        refined using cross-correlation analysis(offsets.c). 
        The non-zero (object) pixels of the reference frame object mask
        (SExtractor OBJECTS image) are stored in a pixel list (x, y, brightness), and
        this list is cross-correlated against the object mask images of the following frames
        in the dither set. The SExtractor OBJECTS image conveniently removes the
        background (important for cross-correlation methods) and identifies the object
        pixels more reliably than a simple thresholding algorithm (e.g., especially in
        images with a non-flat background, large noise, and cosmic rays). Using an
        object list in the cross-correlation focuses on the pixels that contribute to the
        cross-correlation signal and is faster than cross-correlating two images.
        The cross-correlation technique uses coordinate, magnitude, and shape information,
        and was found to be more reliable than matching object coordinate
        lists (the improvement was noticed in extreme cases, like Galactic center images
        and nearly empty fields with an extended galaxy). A subpixel offset measurement
        accuracy of about 0.1 pixels is obtained by fitting a parabola to the peak
        of the cross-correlation image. In terms of speed, this cross-correlation method
        was found to be around 10 times faster (for typical survey data and a relatively large
        search box of 100 pixels) than IRAF STSDAS crosscor.
        
        Although the success of irdr:offsets rate is aprox. 100%, failure is indicated 
        by an offset measurement corresponding exactly to the border of the search 
        area, or a small fraction of object pixels overlapping in the aligned data 
        images.
        
        Parameters
        ----------
                
        images_in: str
            filename of list file (if = None, then use all .skysub.fits files 
            in the out directory --> TB removed !) 
        p_offsets_file: str
            filename of file where offset will be saved (it can be used
            later by dithercubemean)
        
        Returns
        -------    
        offsets: narray          
                two dimensional array with offsets
                
        Notes:
            It assumed that North is up and East is left.
        """
            
        log.info("Starting getPointingOffsets....")
        offsets_mat = None
           
        # STEP 1: Create SExtractor OBJECTS images
        suffix = '_' + self.m_filter + '.skysub.fits'
        #output_list_file=self.out_dir+"/gpo_objs.pap"
        output_fd, output_list_file = tempfile.mkstemp(suffix='.pap', 
                                                       dir=self.out_dir)
        os.close(output_fd)
        os.unlink(output_list_file) # we only need the name
        
        log.debug("Creating OBJECTS images (SExtractor)....")
        
        if self.config_dict:
            mask_minarea = self.config_dict['offsets']['mask_minarea']
            mask_maxarea = self.config_dict['offsets']['mask_maxarea']
            mask_thresh = self.config_dict['offsets']['mask_thresh']
            satur_level = self.config_dict['offsets']['satur_level']
            single_p = self.config_dict['offsets']['single_point']
        else:
            mask_minarea = 5
            mask_maxarea = 200
            mask_thresh = 1.5
            satur_level = 300000
            single_p = True
            
        # In order to set a real value for satur_level, we have to check the
        # number of coadds of the images (NCOADDS or NDIT keywords).
        filelist = [line.replace( "\n", "")
                for line in fileinput.input(images_in)]
        try:
            pf = datahandler.ClFits(filelist[0])
            satur_level = int(satur_level) * int(pf.getNcoadds())
            log.critical("SAT_LEVEL=%s"%satur_level)
        except:
            log.warning("Error read NCOADDS value. Taken default value (=1)")
            satur_level = satur_level
            log.critical("SAT_LEVEL(def)=%s"%satur_level)

        if images_in==None:
            # we use the images ending with suffing in the output directory
            my_input =  self.out_dir + '*' + suffix
        elif os.path.isfile(images_in): 
            # we use the given list of images
            my_input = images_in
        else:
            log.error("Option not recognized !!!")
            raise Exception("Wrong input frames given")

        # Make mask    
        try:
            makeObjMask( my_input, mask_minarea, mask_maxarea, mask_thresh, 
                         satur_level, output_list_file, single_point=single_p)
        except Exception,e:
            log.error("Error making object mask")
            raise e
        
        # STEP 2: Compute dither offsets (in pixles) using cross-correlation technique ==> offsets
        #>mosaic objfiles.nip $off_err > offsets1.nip
        search_box = 50 # half_width of search box in arcsec (default 10)
        offsets_cmd = self.m_irdr_path+'/offsets '+ output_list_file + '  ' + str(search_box) + ' >' + p_offsets_file
        if misc.utils.runCmd( offsets_cmd )==0:
            log.critical("Some error while computing dither offsets")
            raise Exception("Some error while computing dither offsets")
        else:
            try:
                offsets_mat = numpy.loadtxt(p_offsets_file, usecols = (1,2,3)) # columns => (xoffset, yoffset, match fraction) in PIXELS
                # check if correlation overlap fraction is good enough for all offsets computed
                if (offsets_mat[:,2]<self.MIN_CORR_FRAC).sum()>0:
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
            input: file listing the file to coadd and the offsets between
                   each one.
                      
            gain: gain map file to use for the coaddition (it take into account 
                  the BPM)
            
            type_comb : type of combination to use (currently, only average 
                        available):
                        - average: calculate robust mean of stack (using weights)
                        - sum: arithmetic sum of the stack (without weights)
            
        OUTPUTS:
            output : coadded image (and the weight map .weight.fits)
            
        TODO: allow other kind of combination (median, ...see iraf.imcombine)
            
        """
                                                  
        log.info("Start coaddStackImages ...")                                          
        # STEP 1: Define parameters                                          
        input_file = input
        if input_file == None:
            log.error("Bad input file provided !")
            return
        
        if gain == None:
            gain_file = self.out_dir + "/gain_" + self.m_filter + ".fits"
        else:
            gain_file = gain
        
        if output == None:
            output_file = self.out_dir + "/coadd_" + self.m_filter + ".fits"     
        else:
            output_file = output
            
        weight_file = output_file.replace(".fits",".weight.fits")
        
        # STEP 2: Run the coadd                                           
        if type_comb=='average': # (use IRDR::dithercubemean)
            prog = self.m_irdr_path + "/dithercubemean "
            cmd  = prog + " " + input_file + " " + gain_file + " " + output_file + " " + weight_file 
        elif type_comb=='sum':
            # actually, weight_file is not used
            prog = self.m_irdr_path + "/dithercubemean "
            cmd  = prog + " " + input_file + " " + gain_file + " " + output_file + " " + weight_file + " sum " 
        #elif type_comb=='median': # (use IRAF::imcombine)
        else: 
            return (None,None)
        
        e = misc.utils.runCmd( cmd )
        if e==0:
            log.debug("Some error while running command %s", cmd)
            return (None,None)
        else:
            log.debug("Successful ending of coaddStackImages")
            return (output, weight_file)
        
                                              
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
            mask_maxarea = self.config_dict['skysub']['mask_maxarea']
            mask_thresh = self.config_dict['skysub']['mask_thresh']
            satur_level = self.config_dict['skysub']['satur_level']
            single_p = False # we cannot use single point, skyfilter need the whole object!
            dilate = self.config_dict['general']['dilate']
        else:
            mask_minarea = 5
            mask_maxarea = 0 # unlimited
            mask_thresh = 1.5
            satur_level = 55000
            single_p = False
            dilate = 0.5
        
        # In order to set a real value for satur_level, we have to check the
        # number of coadds of the images (NCOADDS or NDIT keywords).
        try:
            pf = datahandler.ClFits(input_file)
            satur_level = int(satur_level) * int(pf.getNcoadds())
        except:
            log.warning("Error read NCOADDS value. Taken default value (=1)")
            satur_level = satur_level

        # Call module makeObjMask
        makeObjMask(input_file, mask_minarea, mask_maxarea, mask_thresh, satur_level,
                   outputfile=self.out_dir + "/objmask_file.txt", single_point=single_p)
        
        if os.path.exists(input_file + ".objs"): 
            shutil.move(input_file + ".objs", output_master_obj_mask)
            log.debug("New Object mask created : %s", output_master_obj_mask)

        # STEP 2: dilate mask
        # mult. scale factor to expand object regions; default is 0.5 (ie, make 50%% larger)
        # For wide-fields with distortion, a dilate > 0.5 is recommended.
        if dilate > 0:
            log.info("Dilating image ....")
            prog = self.m_irdr_path + "/dilate "
            cmd  = prog + " " + output_master_obj_mask + " " + str(dilate)
            # dilate will overwrite the master object mask

            e = misc.utils.runCmd( cmd )
            if e == 0:
                log.debug("Some error while running command %s", cmd)
            else:
                log.debug("Successful ending of createMasterObjMask")
                

        return output_master_obj_mask
                                        
    
    def cleanUpFiles(self, list_dirs):
        """
        Clean up fomer files from the working directory, probably from the 
        last data reduction.
        It is called just before starting a data sequence reduction.
        """
        
        for out_dir in list_dirs:
            # We cannot remove .fits neither .skysub.fits becasue could have 
            # files of a previous reduction of a sequence, for example, when
            # we want to reduce a several sequences of one night, etc.
            misc.fileUtils.removefiles(out_dir + "/c_*", out_dir + "/dc_*",
                                       out_dir + "/*.nip", out_dir + "/*.pap" )
            misc.fileUtils.removefiles(out_dir + "/coadd*", out_dir + "/*.objs",
                                       out_dir + "/uparm*")
            misc.fileUtils.removefiles(out_dir + "/*.head", out_dir + "/*.list",
                                       out_dir + "/*.xml", out_dir + "/*.ldac",
                                       out_dir + "/*.png" )

    def purgeOutput(self):
        """
        Purge the output directory in order to remove all the intermediate files.
        It is called just after finishing a data sequence reduction.
        """
        
        log.info("Purging the output dir ...")
        
        out_dir = self.out_dir
        tmp_dir = self.temp_dir
                 
        misc.fileUtils.removefiles(out_dir+"/*.ldac",out_dir+"/py-sex*",
                                   out_dir+"/*.objs")
        misc.fileUtils.removefiles(out_dir+"/coadd1*", out_dir+"/*_D.fits",
                                       out_dir+"/*_F.fits", out_dir+"/*_D_F.fits" )
        misc.fileUtils.removefiles(out_dir+"/gain*.fits", out_dir+"/masterObjMask.fits",
                                       out_dir+"/*.pap", out_dir+"/*.list", 
                                       out_dir+"/superFlat.fits")
        misc.fileUtils.removefiles(out_dir+"/*.head", out_dir+"/*.txt",
                                       out_dir+"/*.xml")#, out_dir+"/*.png")
       
        misc.fileUtils.removefiles(out_dir + "/*.Q0?.fits")
 
        misc.fileUtils.removefiles(out_dir + "/*.Q0?.fits")
        
        # Temporal directory is emptied
        misc.fileUtils.removefiles(tmp_dir + "/*")
        
        # Remove extension directories
        for i in range(4):
            if os.path.exists(self.out_dir+"/Q%02d"%(i+1)):
                shutil.rmtree(self.out_dir+"/Q%02d"%(i+1), True)
        
    ############# Calibration Stuff ############################################
    def buildCalibrations(self):
        """
        Build the whole master calibration files from the currect calibrations 
        files found in the data set (darks, flats)
        """
        
        log.debug("Start builing the whole calibration files ...")
        # If not initialized, Init DB
        if self.db == None: self.__initDB()
        master_files = []
        
        try:
            master_files += self.reduceSet(self.red_mode, seqs_to_reduce=None, 
                                           types_to_reduce=['DARK','DOME_FLAT',
                                                           'SKY_FLAT'])
        except Exception,e:
            log.error("Some error while builing master calibration files...: %s", str(e))
            raise e
        
        finally:
            log.debug("Calibration files created : %s", master_files)
            return master_files
        
    def buildCalibrations_orig(self):
        """
        ANY MORE USED !!!
        
        Build the whole master calibration files from the currect calibrations 
        files found in the data set (darks, flats).
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
        ANY MORE USED !!!
        
        Look for master flats (sky, twlight,dome or all) files in the data set, 
        group them by FILTER and create the master gainmap, as many as found groups
        
        Return the list of gainmap created
        
        TODO: take into account the possibility to found several master flats 
        and then combine them to build a gain map; at the moment only the first 
        one found is used.
        """
        
        log.debug("Building GainMap for %s Flats", type)
        l_gainmaps = []
        # 1. Look for master flat frames
        full_flat_list = []
        if type=="all":
            full_flat_list = self.db.GetFilesT(type="MASTER_SKY_FLAT", 
                                               texp=-1, filter="ANY")
            full_flat_list+= self.db.GetFilesT(type="MASTER_TW_FLAT", 
                                               texp=-1, filter="ANY")
            full_flat_list+= self.db.GetFilesT(type="MASTER_DOME_FLAT",
                                               texp=-1, filter="ANY")
        elif type=="sky":
            full_flat_list = self.db.GetFilesT(type="MASTER_SKY_FLAT",
                                               texp=-1, filter="ANY")
        elif type=="twlight":
            full_flat_list = self.db.GetFilesT(type="MASTER_TW_FLAT", 
                                               texp=-1, filter="ANY")
        elif type=="dome":
            full_flat_list = self.db.GetFilesT(type="MASTER_DOME_FLAT",
                                               texp=-1, filter="ANY")
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
                #TODO: here we should check if we have more that one master flat, 
                # and if have, then combine them ... 
                # generate a random filename for the master, to ensure we do not
                # overwrite any file
                output_fd, outfile = tempfile.mkstemp(suffix='.fits', 
                                                    dir=self.out_dir)
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
                
                task = reduce.calGainMap.GainMap(group[0], outfile, bpm=None, 
                                               do_normalization=True,
                                               mingain=mingain, maxgain=maxgain, 
                                               nxblock=nxblock,nyblock=nyblock, 
                                               nsigma=nsigma)
                
                out = None
                out = task.create()
                l_gainmaps.append(out) # out must be equal to outfile
            except Exception,e:
                log.error("Some error while creating gainmap: %s",str(e))
                raise e
            if k<len(sorted_list):
                # reset the new group
                group = []
                last_filter = sorted_list[k][1]

        # insert products (gainmaps) into DB
        for f in l_gainmaps: self.db.insert(f)
        self.db.ListDataSet()  
        return l_gainmaps
    

    def buildMasterDarks(self):
        """
        ANY MORE USED !!!
        
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
                output_fd, outfile = tempfile.mkstemp(suffix='.fits', 
                                                        dir=self.out_dir)
                os.close(output_fd)
                os.unlink(outfile) # we only need the name
                task = reduce.calDark.MasterDark (group, self.temp_dir, 
                                                  outfile, texp_scale=False)
                out = task.createMaster()
                l_mdarks.append(out) # out must be equal to outfile
            except Exception,e:
                log.error("Some error while creating master dark: %s",str(e))
                log.error("Proceding to next dark group ...")
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
        ANY MORE USED !!!
        
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
                # Generate a random filename for the master, to ensure we do 
                # not overwrite any file
                output_fd, outfile = tempfile.mkstemp(suffix='.fits', 
                                                    dir=self.out_dir)
                os.close(output_fd)
                os.unlink(outfile) # we only need the name
                #outfile = self.out_dir+"/master_domeflat_%s.fits"%last_filter # added as suffix (FILTER)
                task = reduce.calDomeFlat.MasterDomeFlat(group, self.temp_dir, 
                                                            outfile, None)
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
        ANY MORE USED !!!
        
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
                    # Generate a random filename for the master, to ensure we do 
                    # not overwrite any file
                    output_fd, outfile = tempfile.mkstemp(suffix='.fits', 
                                                        dir=self.out_dir)
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
        ANY MORE USED !!!
        
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
                    # Generate a random filename for the master super flat, to 
                    # ensure we do not overwrite any file
                    output_fd, output_path = tempfile.mkstemp(suffix='.fits', 
                                                            dir=self.out_dir)
                    os.close(output_fd)
                    os.unlink(output_path) # we only need the name
                    #outfile = self.out_dir+"/master_superflat_%s.fits"%last_filter # added as suffix (FILTER)
                    superflat = reduce.SuperSkyFlat(seq, output_path, bpm=None, 
                        norm=False, temp_dir=self.temp_dir)
                    out=superflat.create()
                    l_mflats.append(out)
                except Exception,e:
                    log.error("Some error while creating master SuperFlat: %s",
                        str(e))
                    log.error("but, proceding with next group ...")
                    #raise e
                
        # insert products (master SuperFlats) into DB
        for f in l_mflats: self.db.insert(f)
        self.db.ListDataSet()  
        return l_mflats # a list of master super flats created
    
    def reduceSet(self, red_mode=None, seqs_to_reduce=None, 
                    types_to_reduce=['all']):
        """
        This is the main method for full DataSet reduction supposed it was 
        obtained with the PANIC OT or grouped as 'filter'. 
        
        Main steps:
        
         1. Get all OT/filter sequences
         
         2. For seq in Sequence
    
            ReduceSeq(seq)
            
         3. Insert the results into the local DB 
        
        Parameters
        ----------
         
        red_mode: str
            reduction mode (lemon, quick, science); default mode is 'quick'.
        
        seqs_to_reduce: list
            list of sequence number [0,N-1] to be reduced;
            default (None), all sequences found will be reduced.
        types_to_reduces: str
            Types of sequences to reduce (all, DARK, DOME_FLAT, TW_FLAT, SCIENCE). 
           
        Returns
        -------
            The result files of sequences successfully reduced.
        
        """
        
        log.debug("[reduceSet] Dataset reduction process...")

        reduced_sequences = 0
        files_created = []
        failed_sequences = 0
        failed_sequeces_files = []
        
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
            sequences, seq_types = self.reorder_sequences(sequences, seq_types)
        
        if len(sequences)==0:
            raise Exception("No well-defined sequence to process was found")
        
        k = 0
        for seq,type in zip(sequences, seq_types):
            if k in seqs_to_reduce and ('all' in types_to_reduce or type in types_to_reduce):
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
                    failed_sequeces_files.append(seq)
                    log.error("[reduceSet] Error, cannot reduce sequence : \n %s \n %s"%(str(seq),str(e)))
                    if len(sequences)==1:
                        raise e
                    elif len(sequences)>1:
                        log.debug("[reduceSet] Procceding to next sequence...")
                    
            k = k + 1
    
        # print out the results
        #failed_sequences = len(seqs_to_reduce)-reduced_sequences
        if failed_sequences==len(sequences):
            raise Exception("All sequences to reduce failed !")
        
        log.debug("[reduceSet] All sequences processed.")
        log.debug("[reduceSet] Files generated # %d #: ***"%len(files_created))
        for r_file in files_created: log.debug("\t    - %s"%r_file)
        log.debug("\t    Sequences failed  # %d #: ***"%failed_sequences)
        for seq_failed in failed_sequeces_files:
            log.debug("Seq. failed: \t    - %s\n"%seq_failed) 

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
        
        req_types_order = ['DARK', 'DOME_FLAT', 'TW_FLAT', 'SKY_FLAT', 
                           'FOCUS', 'SCIENCE']
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
    
        Parameters
        ----------
        sequence: list
            list of files of the sequence to be reduced
        
        type: str
            type of sequence (see ClFits.type) (currently not used !) 
        
        Returns
        -------
            List of filenames created by the reduction proccess

        Notes
        -----
        First step is Non-linearity correction because it must be done to all
        frames.
        Calibration files will not be splited for the building of the master 
        calibration file, but science files will. 
        In principle, not too much time could be saved if we split the calibration
        files during the building of master calibrations.
        However, the master calibration files will be split at the stage of
        the processing of scince files. 
        
        """
        
        
        log.debug("[reduceSeq] ** Starting reduction **")

        # print sequence in log file
        for i_file in sequence:
            log.debug(i_file)
        
        #
        # First of all, let see whether Non-linearity correction must be done
        # Note1: we do not need to 'split the extensions', they are processed 
        # one by one (serial) by NonLinearityCorrection(). However, due to
        # NonLinearityCorrection() need MEF files as input, it does the 
        # conversion to MEF if needed. 
        # Note2: NonLinearityCorrection() performs the NLC in parallel,
        # instead of processing the sequence file by file.
        if self.non_linearity_apply==True:
            master_nl = self.non_linearity_model # (lir or rrrmpia)
            try:
                log.info("**** Applying Non-Linearity correction ****")
                log.info("NLC Model: %s"%master_nl)
                nl_task = correctNonLinearity.NonLinearityCorrection(master_nl, 
                            sequence, out_dir=self.temp_dir, suffix='_LC')
                corr_sequence = nl_task.runMultiNLC()

            except Exception,e:
                log.error("Error while applying NL model: %s"%str(e))
                raise e
            
            sequence = corr_sequence

        
        # start creating the pool of process
        n_cpus = self.config_dict['general']['ncpus']
        #n_cpus = multiprocessing.cpu_count()

        # Finally, it was not possible to create 'pool' as a global variable because
        # there are problems with sharing a global pool (or whatever variable) with
        # the Process created in the QL module for dispatching the on-line reduction.
        # So, what we do now is create a new pool everytime we need to do any iraf
        # depending task. In general, we create the pool for whatever task in 
        # this function, but in serial reduction of SCI seqs, the pool will not
        # be used.
        pool = multiprocessing.Pool(processes=n_cpus)
        
        files_created = []
        
        # Take the first file of the sequence in order to find out the type of 
        # the sequence.
        # TODO: I should use 'type' parameter of this method instead to read the
        # type of the first file. I did not found any reason to not use 'type' !
        cfits = datahandler.ClFits(sequence[0])
        
        
        if cfits.isDark():
            log.debug("[reduceSeq] A Dark sequence is going to be reduced: \n%s"%str(sequence))
            try:
                # Generate (and create the file) a random filename for the master, 
                # to ensure we do not overwrite any file
                output_fd, outfile = tempfile.mkstemp(suffix='.fits', 
                                                      prefix='mDark_', 
                                                      dir=self.out_dir)
                os.close(output_fd)
                os.unlink(outfile) # we only need the name
                
                # Check and collapse if required (cube images)
                sequence = misc.collapse.collapse(sequence, out_dir=self.temp_dir)
                
                # Check for EXPT in order to know how to create the master dark 
                # (dark model or fixed EXPT)     
                # Orthodox master dark -- same EXPTIME & NCOADDS
                r = self.checkData(chk_shape=True, chk_filter=True, chk_type=True, 
                                   chk_expt=True, chk_itime=True, 
                                   chk_ncoadd=True, chk_cont=True, 
                                   chk_readmode=True, chk_instrument=True,
                                   file_list=sequence)
                                   
                if (r[0]==True):
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
                    
                    
                elif (r[0]==False and r[1]=="chk_expt"):
                    log.info("Found a dark series of frames with different EXPTIME: Dark model will be created")
                    use_dark_model = True
                    if use_dark_model==True:# and self.red_mode !="quick":
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
                else:
                    log.error("An error in data checking was found. Review your data.")
                    raise Exception("An error in data checking was found. Review your data.")
                        
                if out!=None: files_created.append(out) # out must be equal to outfile
            except Exception,e:
                log.error("[reduceSeq] Some error while creating master DARK: %s",str(e))
                raise e
        elif cfits.isDomeFlatON() or cfits.isDomeFlatOFF():
            log.debug("[reduceSeq] A DomeFlat sequence is going to be reduced: \n%s"%str(sequence))
            try:
                # Generate a random filename for the master, to ensure we do not
                # overwrite any file
                output_fd, outfile = tempfile.mkstemp(suffix='.fits', 
                                                      prefix='mDFlat_', 
                                                      dir=self.out_dir)
                os.close(output_fd)
                os.unlink(outfile) # we only need the name

                # Check and collapse if required (cube images)
                sequence = misc.collapse.collapse(sequence, out_dir=self.temp_dir)

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
        elif cfits.isTwFlat() or cfits.getType(False)=="DOME_FLAT":
            # Added support to proccess DOME_FLATs series (not ON/OFF dome_flat)
            log.debug("[reduceSeq] A TwFlat sequence is going to be reduced: \n%s"%str(sequence))
            try:
                # Look for the required MasterDark (any ExpTime);first in the 
                # local DB (current RS), and if anyone found, then in the 
                # external DB. 
                # Local (Initially, EXPTIME is not a constraint for MASTER_DARK_MODEL);
                master_dark_model = self.db.GetFilesT('MASTER_DARK_MODEL') 
                
                # External (ExpTime is not a constraint)
                if len(master_dark_model) == 0 and self.ext_db != None:
                    log.debug("No MasterDarkModel in current local DB. Trying in external (historic) DB...")
                    master_dark_model = self.ext_db.GetFilesT('MASTER_DARK_MODEL') 
                
                # If no dark_model, look for MASTER_DARKs with any EXPTIME
                if len(master_dark_model) == 0:
                    master_dark_model = None
                    # Look for recently created MASTER_DARKs (PAPI outputs are added to 'db')
                    master_darks = self.db.GetFilesT('MASTER_DARK', -1)
                    # Look for historic/old created MASTER_DARKs (set in config file) 
                    if self.ext_db != None:
                        master_darks += self.ext_db.GetFilesT('MASTER_DARK', -1)
                    if self.master_dark !=None: master_darks += [self.master_dark]
                else:
                    # could there be > 1 master darks, then use the last (mjd sorted)
                    master_dark_model = master_dark_model[-1]
                    master_darks = []

                # If some kind of master darks were found, then go to MasterTwilightFlat
                if len(master_darks) > 0 or master_dark_model != None:
                    log.debug("MASTER_DARKs = %s"%master_darks)
                    log.debug("MASTER_DARK_MODEL = %s"%master_dark_model)
                    
                    # Get filter name for filename
                    filter_name = cfits.getFilter()
                        
                    # generate a random filename for the masterTw, to ensure we do not overwrite any file
                    output_fd, outfile = tempfile.mkstemp(suffix='.fits', 
                                                          prefix='mTwFlat_' + filter_name + '_', 
                                                          dir=self.out_dir)
                        
                    os.close(output_fd)
                    os.unlink(outfile) # we only need the name

                    # Check and collapse if required (cube images)
                    sequence = misc.collapse.collapse(sequence, out_dir=self.temp_dir)

                    m_smooth = self.config_dict['twflats']['median_smooth']
                    
                    task = reduce.calTwFlat.MasterTwilightFlat(sequence, 
                                                               master_dark_model,
                                                               master_darks,
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
                    msg = "No MASTER_DARK or MASTER_DARK_MODEL found. Cannot build Master TwFlat" 
                    log.error(msg)
                    raise Exception(msg)
            except Exception,e:
                log.error("[reduceSeq] Some error while creating master TwFlat: %s",str(e))
                raise e
        elif cfits.isFocusSerie():
            #
            # NOTE: Focus series are not pre-reduced, ie., neither Dark nor Flat
            # Field is applied.
            #
            
            # Because has been found that the automatic focus evaluation based on SExtractor
            # does not work pretty good, it is deactivated for the momment.
            log.warning("[reduceSeq] Proccessing of Focus Serie deactivated")
            return files_created
            
            log.warning("[reduceSeq] Focus Serie is going to be reduced:\n%s"%str(sequence))
            try:
                # 
                # Generate a random filename for the pdf, to ensure we do not
                # overwrite any file.
                output_fd, outfile = tempfile.mkstemp(suffix='.pdf', 
                                                          prefix='focusSer_', 
                                                          dir=self.out_dir)
                os.close(output_fd)
                os.unlink(outfile) # we only need the name

                # Check and collapse if required (cube images)
                sequence = misc.collapse.collapse(sequence, out_dir=self.temp_dir)
                #
                misc.utils.listToFile(sequence, self.temp_dir+"/focus.list")
                pix_scale = self.config_dict['general']['pix_scale']
                satur_level =  self.config_dict['skysub']['satur_level']
                detector = self.config_dict['general']['detector']
                task = reduce.eval_focus_serie.FocusSerie(self.temp_dir+"/focus.list", 
                                                                   str(outfile),
                                                                   pix_scale, 
                                                                   satur_level,
                                                                   show=False,
                                                                   window=detector)
                red_parameters = ()
                result = pool.apply_async(task.eval_serie, red_parameters)
                result.wait()
                best_focus, out = result.get()
                
                pool.close()
                pool.join()

                if out!=None: 
                    files_created.append(out) # out must be equal to outfile
            except Exception,e:
                log.error("[reduceSeq] Error while processing Focus Series: %s",str(e))
                raise e
        elif cfits.isScience():
            l_out_dir = ''
            results = None
            out_ext = []
            log.info("[reduceSeq] Reduction of SCIENCE Sequence: \n%s"%str(sequence))
            if len(sequence) < self.config_dict['general']['min_frames']:
                log.info("[reduceSeq] Found a too SHORT Obs. object sequence.\n\
                 Only %d frames found. Required >%d frames"%(len(sequence),
                                                             self.config_dict['general']['min_frames']))
                raise Exception("Found a short Obs. object sequence. \n\
                Only %d frames found. Required >%d frames" %(len(sequence),
                                                            self.config_dict['general']['min_frames']))
            else:
                # Check and collapse if required (cube images)
                sequence = misc.collapse.collapse(sequence, out_dir=self.temp_dir)
                
                #
                # Get calibration files.
                # If no calibrations are found, the reduction continues
                # without dark subtraction and/or flat-fielding and/or BPM.
                dark, flat, bpm = None, None, None
                if self.red_mode == 'quick' or self.red_mode == 'quick-lemon':
                    # Quick-Mode: optionally calibrations are used.
                    if self.apply_dark_flat==1 or self.apply_dark_flat==2: 
                        dark, flat, bpm = self.getCalibFor(sequence)
                else:
                    # Science-Mode: always calibration are required !
                    # but if not found, it continues without them.
                    dark, flat, bpm = self.getCalibFor(sequence)
                    # Return 3 filenames of master calibration frames (dark, flat, bpm), 


                # Check and split files if required
                obj_ext, next = self.split(sequence) # it must return a list of list (one per each extension)
                dark_ext, cext = self.split([dark])
                flat_ext, cext = self.split([flat])
                bpm_ext, cext = self.split([bpm])
                parallel = self.config_dict['general']['parallel']
                
                # Select the detector to process (Q1, Q2, Q3, Q4, All)
                detector = self.config_dict['general']['detector']
                q = -1 # all
                if next==4:
                    if detector=='Q1': q = 0   # SG1
                    elif detector=='Q2': q = 1 # SG2
                    elif detector=='Q3': q = 2 # SG3
                    elif detector=='Q4': q = 4 # SG4
                    elif detector=='Q123': q = -4 # all except SG4
                    else: q = -1 # all detectors
                    if q >=0 :
                    #if q!=-1:
                        # A single detector will be reduced.
                        obj_ext = [obj_ext[q]]
                        if len(dark_ext) == next: dark_ext = [dark_ext[q]]
                        if len(flat_ext) == next: flat_ext = [flat_ext[q]]
                        if len(bpm_ext) == next: bpm_ext = [bpm_ext[q]]
                        # Reset to 1 the number of extensions
                        next = 1
                    else:
                        # Nothing to do, **all** detectors will be processed
                        pass

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
                            if next==1 and q>=0:
                            #if next==1 and q!=-1: # single detector processing
                                q_ext = q + 1
                            else:
                                q_ext = n + 1
                            log.info("[reduceSeq] ===> (PARALLEL) Reducting detector %d", q_ext)
                            
                            if q==-4 and q_ext==4:
                                # skip detector SG4
                                continue
                            
                            #
                            # For the moment, we have the first calibration file 
                            # for each extension; what rule could we follow ?
                            #
                            if dark_ext==[]: mdark = None
                            else: mdark = dark_ext[n][0] 
                            if flat_ext==[]: mflat = None
                            else: mflat = flat_ext[n][0]
                            if bpm_ext==[]: mbpm = None
                            else: mbpm = bpm_ext[n][0]
                            
                            l_out_dir = self.out_dir + "/Q%02d" % q_ext
                            if not os.path.isdir(l_out_dir):
                                try:
                                    os.mkdir(l_out_dir)
                                except OSError:
                                    log.error("[reduceSeq] Cannot create output directory %s",l_out_dir)
                            else: self.cleanUpFiles([l_out_dir])
                            
                            # async call to procedure
                            #extension_outfilename = l_out_dir + "/" + os.path.basename(self.out_file.replace(".fits",".Q%02d.fits"% (n+1)))
                            extension_outfilename = l_out_dir + "/" + "PANIC_SEQ_Q%02d.fits"% q_ext
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
                            #pool = multiprocessing.Pool(2) #  por que el 2 ???
                            results += [pool.map_async(self.calc,
                                                      [red_parameters])]         
                        # Here is where we WAIT (BLOCKING) for the results 
                        # (result.get() is a blocking call).
                        # If the remote call raised an exception then that 
                        # exception will be reraised by get().

                        for result in results:
                            result.wait()
                            # the 0 index is *ONLY* required if map_async is used !!!
                            out_ext.append(result.get()[0]) 

                        ##for result in results:
                        ##    out_ext.append(result)
                            
                        # Prevents any more tasks from being submitted to the pool. 
                        # Once all the tasks have been completed the worker 
                        # processes will exit.
                        pool.close()
                        # Wait for the worker processes to exit. One must call 
                        #close() or terminate() before using join().
                        pool.join()
                        log.critical("[reduceSeq] DONE PARALLEL REDUCTION ")
                    except Exception,e:
                        log.error("[reduceSeq] Error while parallel data reduction ! --> %s",str(e))
                        raise e
                    
                else:
                    ######## Serial #########
                    for n in range(next):
                        if next == 1 and q >=0 : # single detector processing
                            q_ext = q + 1
                        else:
                            q_ext = n + 1
                        
                        if q == -4 and q_ext == 4:
                            # skip detector SG4
                            continue

                        log.info("[reduceSeq] Entering SERIAL science data reduction ...")    
                        log.info("[reduceSeq] ===> (SERIAL) Reducting extension %d", q_ext)
                        
                        #
                        # For the moment, we take the first calibration file 
                        # for each extension; what rule could we follow ?
                        #                        
                        if dark_ext==[]: mdark = None
                        else: mdark = dark_ext[n][0]
                        
                        if flat_ext==[]: mflat = None
                        else: mflat = flat_ext[n][0]
                        
                        if bpm_ext==[]: mbpm = None
                        else: mbpm = bpm_ext[n][0]
                        
                        l_out_dir = self.out_dir + "/Q%02d" % q_ext
                        if not os.path.isdir(l_out_dir):
                            try:
                                os.mkdir(l_out_dir)
                            except OSError:
                                log.error("[reduceSeq] Cannot create output directory %s",l_out_dir)
                        else: self.cleanUpFiles([l_out_dir])
                        extension_outfilename = l_out_dir + "/" + "PANIC_SEQ_Q%02d.fits"% q_ext
                        
                        try:
                            out_ext.append(self.reduceSingleObj(obj_ext[n],
                                            mdark, mflat,
                                            bpm, self.red_mode,
                                            out_dir=self.out_dir,
                                            output_file = extension_outfilename))
                        except Exception,e:
                            log.error("[reduceSeq] Error while serial data reduction of extension %d of object sequence", q_ext)
                            raise e
        
            # LEMON mode
            # If red_mode is 'lemon', then no warping of frames is required.
            #
            if self.red_mode=='lemon' or self.red_mode=='quick-lemon':
                files_created = out_ext
                log.info("*** Obs. Sequence LEMON-reduced. ***") 
                return files_created
            
            # If all reduction were fine, now join/stick back the extensions in 
            # a wider frame,i.e., build the mosaic.
            #seq_result_outfile = self.out_file.replace(".fits","_SEQ.fits")
            #
            
            #
            # Get the output filename for the reduced sequence
            #            
            try:
                with fits.open(obj_ext[0][0], ignore_missing_end=True) as myhdulist:
                    # Add the PAPI version
                    myhdulist[0].header.set('PAPIVERS', 
                                            __version__, 'PANIC Pipeline version')
                    if 'DATE-OBS' in myhdulist[0].header:
                        seq_result_outfile = self.out_dir + "/PANIC." + myhdulist[0].header['DATE-OBS'] +".fits"
                    else:
                        output_fd, seq_result_outfile = tempfile.mkstemp(suffix='.fits', 
                                                            prefix='PANIC_SEQ_', 
                                                            dir=self.out_dir)
                        os.close(output_fd)
                        os.unlink(seq_result_outfile) # we only need the name
            except Exception,e:
                log.error("Error: %s"%str(e))
                raise e
            
            if len(out_ext)>1:
                log.debug("[reduceSeq] *** Creating final output file *WARPING* single output frames....***")
                #option 1: create a MEF with the results attached, but not warped
                #mef=misc.mef.MEF(outs)
                #mef.createMEF(seq_result_outfile)
                #option 2(current): SWARP result images to register the N-extension into one wide-single extension
                log.debug("*** Coadding/Warping overlapped files....")
                
                swarp = astromatic.SWARP()
                swarp.config['CONFIG_FILE'] = self.papi_home + self.config_dict['config_files']['swarp_conf'] 
                # Note: copy_keywords must be without spaces between keys and 'coma'
                swarp.ext_config['COPY_KEYWORDS'] = "OBJECT,INSTRUME,TELESCOPE,FILTER"\
                "IMAGETYP,FILTER,FILTER1,FILTER2,SCALE,MJD-OBS,HISTORY,NCOADDS,"\
                "NDIT,PAPIVERS"
                swarp.ext_config['IMAGEOUT_NAME'] = seq_result_outfile
                swarp.ext_config['WEIGHTOUT_NAME'] = seq_result_outfile.replace(".fits",".weight.fits")
                # TBC
                #swarp.ext_config['WEIGHT_TYPE'] = 'MAP_WEIGHT'
                #swarp.ext_config['WEIGHT_SUFFIX'] = '.weight.fits'
                swarp.ext_config['RESAMPLE'] = 'Y'
                
                
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
                if os.path.isfile(out_ext[0].replace(".fits", ".weight.fits")):
                    shutil.move(out_ext[0].replace(".fits", ".weight.fits"), 
                                seq_result_outfile.replace(".fits", ".weight.fits"))
                files_created.append(seq_result_outfile)
                log.info("*** Obs. Sequence reduced. File %s created.  ***", 
                         seq_result_outfile)
            else:
                log.error("[reduceSeq] No output files generated by the current Obj.Seq. \
                data reduction ....review your logs files")
            
            # Add the list of original raw-files to the result image header
            if len(out_ext)>0:
                new_frame = fits.open(seq_result_outfile, 'update')
                raw_frames = [os.path.basename(f) for f in sequence]
                new_frame[0].header.add_history("RAW_FRAMES= %s"%str(raw_frames))
                new_frame.close()
                
        else:
            log.error("[reduceSeq] Cannot identify the type of the sequence to reduce ...")
            raise Exception("[reduceSeq] Cannot identify the type of the sequence to reduce ...")    
        
        # Trim/crop the sticked created (science) images and insert them into 
        # the DB
        for file in files_created:
            # if source image was a science one, then products should be also 
            # science images.
            if cfits.isScience():
                # input image (and weight map) is overwritten
                misc.imtrim.imgTrim(file)
            if file!=None and os.path.splitext(file)[1]=='.fits':
                log.debug("Inserting result in DB: %s",file)
                self.db.insert(file)
        
        # not sure if is better to do here ??? 
        #self.purgeOutput()
            
        # Only one file it is expected to be returned
        return files_created
 
        
    def reduceSingleObj_NOT_USED_(self, obj_frames, master_dark, master_flat, master_bpm, 
                  red_mode, out_dir, output_file):
        
        """ 
        Main reduction procedure for a dither sequence. 
        Given a set of object(science) frames and (optionally) master calibration 
        files, run the data reduction of the observing object sequence, producing 
        an reduced ouput frame if no error; otherwise return None or raise exception.
            
        NOTE: Currently this method only accepts single FITS files (not MEF), 
        it means the splitting must be done previusly to call this method.
       
        Parameters
        ----------
        
        obj_frames : list
            list of files to be processed
        
        master_dark : str
            Master dark filename to be used in the processing
            
        master_flat : str
            Master flat filename to be used in the processing
            
        master_bpm: str
            Master BPM filename to be used in the processing (fixing or adding 
            to gainmap).

        red_mode: str
            Reduction mode 

        out_dir: str
            Output directory

        output_file: str
            Output filename to be create

        
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
            if (self.checkData(chk_shape=True, chk_filter=True, chk_type=False, 
                               chk_expt=True, chk_itime=True, chk_ncoadd=False, 
                               chk_readmode=True, chk_instrument=True)[0]==True):
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
            # overwrite initial given observing mode
            self.obs_mode = self.getObsMode()  
        except:
            raise
        
        log.info( "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" )
        log.info( "OBSERVING SEQUENCE DETECTED. OBS_MODE= %s", self.obs_mode)
        log.info( "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" )
        
        # keep a copy of the original file names
        self.m_rawFiles = self.m_LAST_FILES        

        ########################################################################
        # 0 - Bad Pixels; two options:
        # To Fix: replace with a bi-linear interpolation from nearby pixels.
        # To add to gainmap:  to set bad pixels to bkg level 
        # Both options are incompatible.
        ########################################################################
        if self.config_dict['bpm']['mode'].lower()=='none':
            master_bpm_4gain = None
            master_bpm_4fix = None
        elif self.config_dict['bpm']['mode'].lower()=='fix':
            master_bpm_4gain = None
            master_bpm_4fix = master_bpm
        elif self.config_dict['bpm']['mode'].lower()=='grab':
            master_bpm_4gain = master_bpm
            master_bpm_4fix = None
        else:
            master_bpm_4gain = None
            master_bpm_4fix = None

        ########################################################################
        # 1 - Apply dark, flat to ALL files 
        ########################################################################
        if self.apply_dark_flat==1 and \
            (master_dark!=None or master_flat!=None or master_bpm_4fix!=None):
            log.info("**** Applying Dark, Flat and BPM ****")
            res = reduce.ApplyDarkFlat(self.m_LAST_FILES, 
                                       master_dark, 
                                       master_flat,
                                       master_bpm_4fix, 
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
        gainmap = out_dir + '/gain_' + self.m_filter + '.fits'
        # get gainmap parameters
        if self.config_dict:
            mingain = self.config_dict['gainmap']['mingain']
            maxgain = self.config_dict['gainmap']['maxgain']
            nxblock = self.config_dict['gainmap']['nxblock']
            nyblock = self.config_dict['gainmap']['nyblock']
            nsigma = self.config_dict['gainmap']['nsigma']
        else:
            mingain = 0.2
            maxgain = 1.8
            nxblock = 16
            nyblock = 16
            nsigma = 5
            
        # When gainmap is created (from dome or sky flats), it must be normalized
        # wrt mode of chip 1 to get gain differences, set bad pixels, 
        # outlier set =0 (e.g. pixels deviating >5 sigma from local median,
        # pixels deviating >30%(?),...
        # do_normalization=False because it is suppossed that FF is already 
        # normalized.
        g = reduce.calGainMap.GainMap(local_master_flat, gainmap, 
                                    bpm=master_bpm_4gain, 
                                    do_normalization=False,  
                                    mingain=mingain, maxgain=maxgain,
                                    nxblock=nxblock, nyblock=nyblock,
                                    nsigma=nsigma)
        g.create() 
           
        ########################################################################
        # Add/combine external Bad Pixel Map to gainmap (maybe from master DARKS,FLATS ?)
        # BPMask ==> Bad pixeles >0, Good pixels = 0
        # GainMap ===> Bad pixels <=0, Good pixels = 1
        # Later, in skyfilter procedure bad pixels are set to bkg level.
        ########################################################################
        if master_bpm_4gain !=None:
            log.info("Combinning external BPM (%s) and Gainmap (%s)"%(master_bpm_4gain, gainmap))
            if not os.path.exists( master_bpm_4gain ):
                log.error("No external Bad Pixel Mask found. Cannot find file : %s"%master_bpm_4gain)
            else:
                # Convert badpix (>0) to 0 value, and goodpix (=0) to >1.0
                bpm_data = fits.getdata(master_bpm_4gain, header=False)
                badpix_p = numpy.where(bpm_data>0)
                gain_data, gh = fits.getdata(gainmap, header=True)
                gain_data[badpix_p] = 0
                gh.set('HISTORY','Combined with BPM:%s'%master_bpm_4gain)
                fits.writeto(gainmap, gain_data, header=gh, clobber=True)
        
        ########################################################################
        # 3b - Lab mode: it means only D,FF and FWHM estimation is done for 
        #      each raw frame.  
        ########################################################################
        if self.red_mode=='lab':
            log.info("**** Computing FWHM of each pre-reduced image ****")
            # First, get input parameter for FWHM estimation
            satur_level = self.config_dict['astrometry']['satur_level']
            pix_scale = self.config_dict['general']['pix_scale']
            fwhm_t = []
            std_t = []
            for r_file in self.m_LAST_FILES:
                try:
                    pf = datahandler.ClFits(r_file)
                    satur_level = int(satur_level) * int(pf.getNcoadds())
                except:
                    log.warning("Cannot get NCOADDS. Taken default (=1)")
                    satur_level = satur_level
                
                # Note that we could computer FWHM for a well-defined area
                # using edge_x and edge_y parameters.    
                cq = reduce.checkQuality.CheckQuality(r_file, isomin=10.0,
                    ellipmax=0.3, edge_x=100, edge_y=100, pixsize=pix_scale, 
                    gain=1, sat_level=satur_level)
                try:
                    (fwhm, std, k, k) = cq.estimateFWHM()
                    if fwhm>0 and fwhm<20:
                        log.info("File %s - FWHM = %s (pixels) std= %s"%(r_file, fwhm, std))
                        fwhm_t.append(fwhm)
                        std_t.append(std)
                    elif fwhm<0 or fwhm>20:
                        log.error("Wrong estimation of FWHM of file %s"%r_file)
                        log.error("Please, review your data.")           
                    else:
                        log.error("ERROR: Cannot estimate FWHM of file %s"%r_file)           
                except Exception,e:
                    log.error("ERROR: something wrong while computing FWHM")
                    raise e
            
            log.info("**********************************")
            log.info("*** Mean-FWHM = %s"%numpy.mean(fwhm_t))
            log.info("*** Mean-STD = %s"%numpy.mean(std_t))
            log.info("*** END of Lab data reduction ****")
            return self.m_LAST_FILES[0]
            
                        
        ########################################################################
        # 4 - First Sky subtraction (IRDR) - sliding window technique
        #     (and set bad pixels to bkg lvl)
        ########################################################################
        log.info("**** 1st Sky subtraction (without object mask) ****")
        misc.utils.listToFile(self.m_LAST_FILES, out_dir+"/skylist1.list")
        # Return only filtered images; in exteded-sources, sky frames  are not 
        # included.
        #  
        self.m_LAST_FILES = self.skyFilter(out_dir+"/skylist1.list",
                                           gainmap, 'nomask', self.obs_mode)
        
        ########################################################################
        # 4.1 - Remove crosstalk - (only if bright stars are present)
        #       (if red_mode!=quick, then crosstalk will be done later)
        ########################################################################
        if ((self.red_mode=='quick' or self.red_mode=='quick-lemon') and 
            self.config_dict['general']['remove_crosstalk']):
            log.info("**** Removing crosstalk ****")
            try:
                res = map(reduce.dxtalk.remove_crosstalk, self.m_LAST_FILES, 
                            [None]*len(self.m_LAST_FILES), 
                            [True]*len(self.m_LAST_FILES))
                self.m_LAST_FILES = res
            except Exception,e:
                raise e

        ########################################################################
        # 4.2 - Divide by the master flat after sky subtraction ! 
        # some people think it is better do it now (M.J.Irwing, CASU)
        # But, dark is not applied/required, as it was implicity done when sky 
        # subtraction.
        # However, it does not produce good results, it looks like if we undo the 
        # sky subtraction, because sky subtraction is quite related with 
        # flatfielding.
        # For further details, check TN and Wei-Hao (SIMPLE) mails. 
        # So, it is implemented here only for academic purposes !
        ########################################################################
        if self.apply_dark_flat==2 and master_flat!=None:
            log.info("**** Applying Flat AFTER sky subtraction ****")
            res = reduce.ApplyDarkFlat(self.m_LAST_FILES, 
                                       None,  
                                       master_flat,
                                       None, 
                                       out_dir)
            self.m_LAST_FILES = res.apply()

        ########################################################################
        # 4.3 - LEMON connection - End here for Quick-LEMON-1 processing    
        ########################################################################
        if self.red_mode=='quick-lemon':
            output_fd, papi_output = \
                tempfile.mkstemp(prefix = out_dir,
                                 suffix = '.papi_output', text = True)
            os.close(output_fd)
            os.unlink(papi_output)
            misc.utils.listToFile(self.m_LAST_FILES, papi_output)
            log.info("1st Skysubtraction done !")
            return papi_output
                           
        ########################################################################
        # 5 - Quality assessment (FWHM, background, sky transparency, 
        # ellipticity, PSF quality)  
        ########################################################################
                            
        log.info("**** Data Quality Assessment **** (TBD)")                   

        ########################################################################
        # -- una prueba con astrowarp : no va bien; a simple vista da resultados 
        # parecidos, y en CPU tambien =, por tanto, opcion a considerar, pero
        # no funciona en todos los casos, pues en cuanto falla SCAMP no se
        # puede hacer el stack !
        # No obstante, con Astrometry.net igual si se puede conseguir algo estable ...
        # 6b - Computer dither offsets and coadd
        ########################################################################
        #self.coadd_mode = 'swarp'
        #self.coadd_mode = 'dithercubemean'
        if self.coadd_mode=='swarp':
            if self.obs_mode!='dither' or self.red_mode=="quick":
                log.info("**** Doing Astrometric calibration and coaddition result frame ****")
                #misc.utils.listToFile(self.m_LAST_FILES, out_dir+"/files_skysub.list")
                aw = reduce.astrowarp.AstroWarp(self.m_LAST_FILES, catalog="USNO-B1", 
                              coadded_file=output_file, config_dict=self.config_dict)
                try:
                    aw.run(engine=self.config_dict['astrometry']['engine'])
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
        # 6new - Compute dither offsets from the first sky subtracted/filtered 
        # images using astrometric calibration (Astrometric.Net).
        ########################################################################
        misc.utils.listToFile(self.m_LAST_FILES, out_dir+"/files_skysub.list")
        new_files = []
        if True:
            for my_file in self.m_LAST_FILES:
                # Run astrometric calibration
                try:
                    solved = reduce.solveAstrometry.solveField(my_file, 
                                  out_dir, # self.temp_dir produces collision
                                  self.config_dict['general']['pix_scale'])
                except Exception,e:
                    raise Exception("[solveAstrometry] Cannot solve Astrometry for file: %s \n%s"%(my_file, str(e)))
                else:
                    # Rename the file
                    out_filename = my_file.replace(".fits", ".ast.fits")
                    new_files.append(out_filename)
                    shutil.move(solved, out_filename)
            
            self.m_LAST_FILES = new_files
            try:
                offset_mat = self.getWCSPointingOffsets(self.m_LAST_FILES, 
                                                    out_dir + '/offsets1.pap')                
            except Exception,e:
                log.error("Error while getting WCS pointing offsets. Cannot continue with data reduction...")
                raise e
        
        ########################################################################
        # 6old - Compute dither offsets from the first sky subtracted/filtered 
        # images using cross-correlation (SExtractor + offsets)
        ########################################################################
        if False:
            misc.utils.listToFile(self.m_LAST_FILES, out_dir+"/files_skysub.list")
            try:
                offset_mat = self.getPointingOffsets(out_dir+"/files_skysub.list", 
                                                    out_dir+'/offsets1.pap')                
            except Exception,e:
                log.error("Error while getting pointing offsets. Cannot continue with data reduction...")
                raise e
        
        ########################################################################
        # 7 - First pass coaddition using offsets
        ########################################################################
        log.info("**** Initial coaddition of sky subtracted frames ****")
        fo = open(out_dir+'/offsets1.pap',"r")
        fs = open(out_dir+'/stack1.pap','w+')
        for line in fo:
            n_line = line.replace(".fits.objs", ".fits") 
            fs.write(n_line)
        fo.close()
        fs.close()
        self.coaddStackImages(out_dir+'/stack1.pap', gainmap, 
                              out_dir+'/coadd1.fits','average')
         
        #        
        # Remove (trim or crop) any frame around the image due to the coadd process        
        # Finally, it was decided to be done at the end of reduceSeq(), because
        # astrometric calibration (SWARP) is also adding a new border that 
        # need to be removed as well.
        #misc.imtrim.imgTrim(out_dir+'/coadd1.fits')
                             
        ########################################################################
        # End of first cycle: SINGLE REDUCTION (quick mode or extended object !) 
        ########################################################################
        if self.obs_mode!='dither' or self.red_mode=="quick":
            log.info("**** Doing Astrometric calibration of coadded result frame ****")
           
            if self.config_dict['astrometry']['engine']=='SCAMP':
                try:
                    reduce.astrowarp.doAstrometry(out_dir+'/coadd1.fits', output_file, 
                                           self.config_dict['astrometry']['catalog'], 
                                           config_dict=self.config_dict, 
                                           do_votable=False)
                except Exception,e:
                    raise Exception("[astrowarp] Cannot solve Astrometry %s"%str(e))
            else:
                try:
                    solved = reduce.solveAstrometry.solveField(out_dir+'/coadd1.fits', 
                                                    out_dir, # self.temp_dir produces collision
                                                    self.config_dict['general']['pix_scale'])

                except Exception,e:
                    raise Exception("[solveAstrometry] Cannot solve Astrometry %s"%str(e))
                else:
                    # Rename the file
                    shutil.move(solved, output_file)


            log.info("Generated output file ==>%s", output_file)

            if self.config_dict['general']['estimate_fwhm']:
                log.info("**** FWHM estimation of coadded_1 result frame ****")
                satur_level = self.config_dict['astrometry']['satur_level']
                pix_scale = self.config_dict['general']['pix_scale']
                
                try:
                    pf = datahandler.ClFits(output_file)
                    satur_level = int(satur_level) * int(pf.getNcoadds())
                except:
                    log.warning("Cannot get NCOADDS. Taken default (=1)")
                    satur_level = satur_level
                    
                cq = reduce.checkQuality.CheckQuality(output_file, isomin=10.0,
                    ellipmax=0.3, edge_x=200, edge_y=200, pixsize=pix_scale, 
                    gain=4.15, sat_level=satur_level)
                try:
                    (fwhm, std, k, k) = cq.estimateFWHM()
                    if fwhm>0 and fwhm<20:
                        log.info("File %s - FWHM = %s (pixels) std= %s"%(output_file, fwhm, std))
                    elif fwhm<0 or fwhm>20:
                        log.error("Wrong estimation of FWHM of file %s"%output_file)
                        log.error("Please, review your data.")           
                    else:
                        log.error("ERROR: Cannot estimate FWHM of file %s"%output_file)           
                except Exception,e:
                    log.error("ERROR: something wrong while computing FWHM")
                    raise e
                            
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
        # TODO 
        
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
                j = j+1
            elif self.obs_mode=='dither':
                j = j+1
            i = i+1
        
        fs.close()
        self.m_LAST_FILES = self.skyFilter(out_dir+"/skylist2.pap", gainmap, 
                                           'mask', self.obs_mode)      
    
        ########################################################################
        # 9.1 - Remove crosstalk - (only if bright stars are present)    
        ########################################################################
        if self.config_dict['general']['remove_crosstalk']:
            log.info("**** Removing crosstalk ****")
            try:
                res = map(reduce.dxtalk.remove_crosstalk, self.m_LAST_FILES, 
                         [None]*len(self.m_LAST_FILES), [True]*len(self.m_LAST_FILES))
                self.m_LAST_FILES = res
            except Exception,e:
                raise e

        ########################################################################
        # 9.2 - Remove Cosmic Rays -    
        ########################################################################
        if self.config_dict['general']['remove_cosmic_ray']:
            try:
                res = map(reduce.remove_cosmics.remove_cr, self.m_LAST_FILES, 
                            [None]*len(self.m_LAST_FILES), [True]*len(self.m_LAST_FILES),
                            [False]*len(self.m_LAST_FILES))
                self.m_LAST_FILES = res
            except Exception,e:
                raise e
        ########################################################################
        # 9.3 - LEMON connection - End here for LEMON processing    
        ########################################################################
        
        if self.red_mode=='lemon':
            output_fd, papi_output = \
                tempfile.mkstemp(prefix = out_dir,
                                 suffix = '.papi_output', text = True)
            os.close(output_fd)
            os.unlink(papi_output)
            misc.utils.listToFile(self.m_LAST_FILES, papi_output)
            log.info("End of sequence LEMON-reduction. # %s # files created. ",
                     len(self.m_LAST_FILES))
            return papi_output
    
        #######################################################################
        # 9.4 - Divide by the master flat after sky subtraction ! (see notes above)
        # (the same task as above 4.2) --> HAS NO SENSE !!! only for a test ??? or a.l.a. O2k 
        #######################################################################
        if self.apply_dark_flat==2 and master_flat!=None:
            log.info("**** Applying Flat AFTER sky subtraction ****")
            res = reduce.ApplyDarkFlat(self.m_LAST_FILES, 
                                       None,  
                                       master_flat,
                                       None, 
                                       out_dir)
            self.m_LAST_FILES = res.apply()
            
        #######################################################################
        # 10a - Compute field distortion and make final stack:
        #       1-Remove field distortion from individual images (SCAMP+SWARP)
        #       2-Coaddition of corrected field distortion images (SWARP)
        #       3-Final Astrometric calibration (SCAMP) of the coadded image
        #######################################################################
        if self.coadd_mode=='swarp':
            print "astrowarp--->LAST_FILES=",self.m_LAST_FILES
            
            log.info("**** Astrometric calibration and stack of individual \
            frames to field distortion correction ****")
            aw = reduce.astrowarp.AstroWarp(self.m_LAST_FILES, 
                                            coadded_file=output_file, 
                                            config_dict=self.config_dict,
                                            do_votable=False)
            try:
                aw.run(engine=self.config_dict['astrometry']['engine'])
            except Exception,e:
                log.error("Some error while running Astrowarp....")
                raise e
            log.info("Successful end of Pipeline (I hope!)")
            return output_file
        
    
        #######################################################################
        # 10b.1 - Create second coadded image of the dithered stack using new  
        # sky subtracted frames (using the same offsets)
        #######################################################################
        log.info("**** Coadding image free distorion frames ****")
        self.coaddStackImages(out_dir + '/stack1.pap', gainmap, 
                              out_dir + '/coadd2.fits')
        
        #        
        # Remove (trim or crop) any frame around the image due to the coadd process        
        # Finally, it was decided to be done at the end of reduceSeq(), because
        # astrometric calibration (SWARP) is also adding a new border that 
        # need to be removed as well.        
        #misc.imtrim.imgTrim(out_dir+'/coadd2.fits')
        
        #######################################################################
        # 10b.2 - Make final Astrometry
        #######################################################################
        log.info("**** Computing Astrometric calibration of coadded (2nd) result frame ****")
        
        if self.config_dict['astrometry']['engine']=='SCAMP':
            try:
                reduce.astrowarp.doAstrometry(out_dir+'/coadd2.fits', output_file, 
                                      self.config_dict['astrometry']['catalog'], 
                                      config_dict=self.config_dict,
                                      do_votable=False,
                                      resample=True, # means remove field distorion
                                      subtract_back=True)
            except Exception,e:
                raise Exception("[reductionset] Cannot solve Astrometry %s"%str(e))
        else:
            try:
                solved = reduce.solveAstrometry.solveField(out_dir+'/coadd2.fits', 
                                                    self.temp_dir,
                                                    self.config_dict['general']['pix_scale'])
            except Exception,e:
                raise Exception("[reductionset] Cannot solve Astrometry %s"%str(e))
            else:
                # Rename the file
                shutil.move(solved, output_file)


        log.info("Generated output file ==>%s", output_file)
        
        os.chdir(old_cwd)
        
        log.info("##################################")
        log.info("##### End of data reduction ######")
        log.info("##################################")
        
        return output_file 

    def reduceSingleObj(self, obj_frames, master_dark, master_flat, master_bpm, 
                  red_mode, out_dir, output_file):
        
        """ 
        Main reduction procedure for a dither sequence. 
        Given a set of object(science) frames and (optionally) master calibration 
        files, run the data reduction of the observing object sequence, producing 
        an reduced ouput frame if no error; otherwise return None or raise exception.
            
        NOTE: Currently this method only accepts single FITS files (not MEF), 
        it means the splitting must be done previusly to call this method.
       
        Parameters
        ----------
        
        obj_frames : list
            list of files to be processed
        
        master_dark : str
            Master dark filename to be used in the processing
            
        master_flat : str
            Master flat filename to be used in the processing
            
        master_bpm: str
            Master BPM filename to be used in the processing (fixing or adding 
            to gainmap).

        red_mode: str
            Reduction mode 

        out_dir: str
            Output directory

        output_file: str
            Output filename to be create

        
        Returns
        -------
        
        output_file : str 
            If 'red_mode' = 'quick' or 'science', then the function return the 
            coadded frame obtained from the reduction, both 'quick' and 'science' 
            reduction mode.
            If 'red_mode' = 'lemon', the method return a file listing the files
            obtained of the pre-processing done (dark, flat and sky subtraction) 
        
        """
        
        c_filter = datahandler.ClFits(obj_frames[0]).getFilter()
        log.info("#########################################")
        log.info("#### Starting Object Data Reduction #####")
        log.info("#### MODE = %s  ", self.red_mode)
        log.info("#### OUT_DIR = %s ",out_dir)
        log.info("#### OUT_FILE = %s ", output_file)
        log.info(" ----------------------------------")
        log.info("#### FILTER = %s", c_filter)
        log.info("#### MASTER_DARK = %s ", master_dark)
        log.info("#### MASTER_FLAT = %s ", master_flat)
        log.info("#### MASTER_BPM = %s ", master_bpm)
        log.info("#########################################")
        
        # if given, set the reduction mode, else use default mode.
        if red_mode != None: self.red_mode = red_mode
        
        
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
            
        
        ########################################################################
        # 0 - Some checks (filter, ....) 
        ########################################################################
        # TODO : it could/should be done in reduceSeq, to avoid the spliting ...??
        log.info("**** Data Validation ****")
        if self.check_data:
            if (self.checkData(chk_shape=True, chk_filter=True, chk_type=False, 
                               chk_expt=True, chk_itime=True, chk_ncoadd=False, 
                               chk_readmode=True, chk_instrument=True)[0]==True):
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
            # overwrite initial given observing mode
            self.obs_mode = self.getObsMode()  
        except:
            raise
        
        log.info( "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" )
        log.info( "OBSERVING SEQUENCE DETECTED. OBS_MODE= %s", self.obs_mode)
        log.info( "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" )
        
        # keep a copy of the original file names
        self.m_rawFiles = self.m_LAST_FILES        

        ########################################################################
        # 0 - Bad Pixels; two options:
        # To Fix: replace with a bi-linear interpolation from nearby pixels.
        # To add to gainmap:  to set bad pixels to bkg level 
        # Both options are incompatible.
        ########################################################################
        if self.config_dict['bpm']['mode'].lower()=='none':
            master_bpm_4gain = None
            master_bpm_4fix = None
        elif self.config_dict['bpm']['mode'].lower()=='fix':
            master_bpm_4gain = None
            master_bpm_4fix = master_bpm
        elif self.config_dict['bpm']['mode'].lower()=='grab':
            master_bpm_4gain = master_bpm
            master_bpm_4fix = None
        else:
            master_bpm_4gain = None
            master_bpm_4fix = None

        ########################################################################
        # 1 - Apply dark, flat to ALL files 
        ########################################################################
        if self.apply_dark_flat==1 and \
            (master_dark!=None or master_flat!=None or master_bpm_4fix!=None):
            log.info("**** Applying Dark, Flat and BPM ****")
            res = reduce.ApplyDarkFlat(self.m_LAST_FILES, 
                                       master_dark, 
                                       master_flat,
                                       master_bpm_4fix, 
                                       out_dir)
            self.m_LAST_FILES = res.apply()
        
        ########################################################################
        # 2 - Compute Super Sky Flat-Field --> GainMap
        #      The 'local_master_flat' is ONLY used to build the GainMap, and
        #      can come from the 'master_flat' provided by the user or from the 
        #      just created 'super_flat' with sci images 
        ########################################################################
        if master_flat == None:
            try:
                # - Find out what kind of observing mode we have (dither, ext_dither, ...)
                log.info('**** Computing Super-Sky Flat-Field (local_master_flat) ****')
                local_master_flat = out_dir+"/superFlat.fits"
                if self.obs_mode == "dither":
                    log.debug("---> dither sequece <----")
                    misc.utils.listToFile(self.m_LAST_FILES, out_dir+"/files.list")
                    superflat = reduce.SuperSkyFlat(out_dir+"/files.list", 
                                                    local_master_flat, bpm=None, 
                                                    norm=True, 
                                                    temp_dir=self.temp_dir)
                    superflat.create()
                elif (self.obs_mode == "dither_on_off" or
                      self.obs_mode == "dither_off_on" or
                      self.obs_mode == "other"):
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
        gainmap = out_dir + '/gain_' + self.m_filter + '.fits'
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
        # do_normalization=False because it is suppossed that FF is already 
        # normalized.
        g = reduce.calGainMap.GainMap(local_master_flat, gainmap, 
                                    bpm=master_bpm_4gain, 
                                    do_normalization=False,
                                    mingain=mingain, maxgain=maxgain,
                                    nxblock=nxblock, nyblock=nyblock,
                                    nsigma=nsigma)
        g.create() 
           
        ########################################################################
        # Add/combine external Bad Pixel Map to gainmap (maybe from master DARKS,FLATS ?)
        # BPMask ==> Bad pixeles >0, Good pixels = 0
        # GainMap ===> Bad pixels <=0, Good pixels = 1
        # Later, in skyfilter procedure bad pixels are set to bkg level.
        ########################################################################
        if master_bpm_4gain != None:
            log.info("Combinning external BPM (%s) and Gainmap (%s)"%(master_bpm_4gain, gainmap))
            if not os.path.exists( master_bpm_4gain ):
                log.error("No external Bad Pixel Mask found. Cannot find file : %s"%master_bpm_4gain)
            else:
                # Convert badpix (>0) to 0 value, and goodpix (=0) to >1.0
                bpm_data = fits.getdata(master_bpm_4gain, header=False)
                badpix_p = numpy.where(bpm_data>0)
                gain_data, gh = fits.getdata(gainmap, header=True)
                gain_data[badpix_p] = 0
                gh.set('HISTORY','Combined with BPM:%s'%master_bpm_4gain)
                fits.writeto(gainmap, gain_data, header=gh, clobber=True)
        
        ########################################################################
        # 3b - Lab mode: it means only D,FF and FWHM estimation is done for 
        #      each raw frame.
        ########################################################################
        if self.red_mode == 'lab':
            log.info("**** Computing FWHM of each pre-reduced image ****")
            # First, get input parameter for FWHM estimation
            satur_level = self.config_dict['astrometry']['satur_level']
            pix_scale = self.config_dict['general']['pix_scale']
            fwhm_t = []
            std_t = []
            for r_file in self.m_LAST_FILES:
                try:
                    pf = datahandler.ClFits(r_file)
                    satur_level = int(satur_level) * int(pf.getNcoadds())
                except:
                    log.warning("Cannot get NCOADDS. Taken default (=1)")
                    satur_level = satur_level
                
                # Note that we could compute FWHM for a well-defined area
                # using edge_x and edge_y parameters.    
                cq = reduce.checkQuality.CheckQuality(r_file, isomin=10.0,
                    ellipmax=0.3, edge_x=100, edge_y=100, pixsize=pix_scale, 
                    gain=1, sat_level=satur_level)
                try:
                    (fwhm, std, k, k) = cq.estimateFWHM()
                    if fwhm>0 and fwhm<20:
                        log.info("File %s - FWHM = %s (pixels) std= %s"%(r_file, fwhm, std))
                        fwhm_t.append(fwhm)
                        std_t.append(std)
                    elif fwhm<0 or fwhm>20:
                        log.error("Wrong estimation of FWHM of file %s"%r_file)
                        log.error("Please, review your data.")           
                    else:
                        log.error("ERROR: Cannot estimate FWHM of file %s"%r_file)           
                except Exception,e:
                    log.error("ERROR: something wrong while computing FWHM")
                    raise e
            
            log.info("**********************************")
            log.info("*** Mean-FWHM = %s"%numpy.mean(fwhm_t))
            log.info("*** Mean-STD = %s"%numpy.mean(std_t))
            log.info("*** END of Lab data reduction ****")
            return self.m_LAST_FILES[0]
            
                        
        ########################################################################
        # 4 - First Sky subtraction (IRDR) - sliding window technique
        #     (and set bad pixels to bkg lvl)
        ########################################################################
        log.info("**** 1st Sky subtraction (without object mask) ****")
        misc.utils.listToFile(self.m_LAST_FILES, out_dir + "/skylist1.list")
        # Return only filtered images; in extended-sources, sky frames  are not 
        # included.
        #  
        self.m_LAST_FILES = self.skyFilter(out_dir + "/skylist1.list",
                                           gainmap, 'nomask', self.obs_mode)
        
        ########################################################################
        # 4.1 - Remove crosstalk - (only if bright stars are present)
        #       (if red_mode!=quick, then crosstalk will be done later)
        ########################################################################
        # Always remove crosstalk if enabled in config.
        if self.config_dict['general']['remove_crosstalk']:
        #if ((self.red_mode=='quick' or self.red_mode=='quick-lemon') and 
        #    self.config_dict['general']['remove_crosstalk']):
            log.info("**** Removing crosstalk ****")
            try:
                res = map(reduce.dxtalk.remove_crosstalk, self.m_LAST_FILES, 
                            [None]*len(self.m_LAST_FILES), 
                            [True]*len(self.m_LAST_FILES))
                self.m_LAST_FILES = res
            except Exception,e:
                raise e      

        ########################################################################
        # 4.2 - Divide by the master flat after sky subtraction ! 
        # some people think it is better do it now (M.J.Irwing, CASU)
        # But, dark is not applied/required, as it was implicity done when sky 
        # subtraction.
        # However, it does not produce good results, it looks like if we undo the 
        # sky subtraction, because sky subtraction is quite related with 
        # flatfielding.
        # For further details, check TN and Wei-Hao (SIMPLE) mails. 
        # So, it is implemented here only for academic purposes !
        ########################################################################
        if self.apply_dark_flat==2 and master_flat!=None:
            log.info("**** Applying Flat AFTER sky subtraction ****")
            res = reduce.ApplyDarkFlat(self.m_LAST_FILES, 
                                       None,  
                                       master_flat,
                                       None, 
                                       out_dir)
            self.m_LAST_FILES = res.apply()
        
        ########################################################################
        # Preliminary Astrometric calibration of sky-subtracted frames.
        # ######################################################################
        if self.config_dict['offsets']['method'] == 'wcs':
            log.info("**** Preliminary Astrometric calibration ****")
            new_files = []
            for my_file in self.m_LAST_FILES:
                # Run astrometric calibration
                try:
                    solved = reduce.solveAstrometry.solveField(my_file, 
                                out_dir, # self.temp_dir produces collision
                                self.config_dict['general']['pix_scale'])
                except Exception,e:
                    raise Exception("[solveAstrometry] Cannot solve Astrometry for file: %s \n%s"%(my_file, str(e)))
                else:
                    # Rename the file
                    out_filename = my_file.replace(".fits", ".ast.fits")
                    new_files.append(out_filename)
                    shutil.move(solved, out_filename)
            
            self.m_LAST_FILES = new_files
        
        ########################################################################
        # 4.3 - LEMON connection - End here for Quick-LEMON-1 processing    
        ########################################################################
        if self.red_mode == 'quick-lemon':
            output_fd, papi_output = \
                tempfile.mkstemp(prefix = out_dir,
                                 suffix = '.papi_output', text = True)
            os.close(output_fd)
            os.unlink(papi_output)
            misc.utils.listToFile(self.m_LAST_FILES, papi_output)
            log.info("Quick-LEMON output generated !")
            return papi_output
         
        ########################################################################
        # 5 - Compute dither offsets:
        #        - using astrometric calibration (wcs)
        #        - using cross-correlation (irdr.offsets)
        ########################################################################
        self.offsets_method = self.config_dict['offsets']['method']
        if self.offsets_method == 'wcs':
            log.info("Computing dither offsets using astrometric calibration")
            try:
                offset_mat = self.getWCSPointingOffsets(self.m_LAST_FILES, 
                                                        out_dir + '/offsets1.pap')                
            except Exception,e:
                log.error("Error while computing WCS pointing offsets. Cannot continue with data reduction...")
                raise e
        else:
            log.info("Computing dither offsets using cross-correlation")
            misc.utils.listToFile(self.m_LAST_FILES, out_dir + "/files_skysub.list")
            try:
                offset_mat = self.getPointingOffsets(out_dir + "/files_skysub.list", 
                                                    out_dir + '/offsets1.pap')
            except Exception,e:
                log.error("Error while getting pointing offsets. Cannot continue with data reduction...")
                raise e
        
        ########################################################################
        # End of first cycle: SINGLE REDUCTION (quick mode or extended object !)    
        # 6 - Preliminary coaddition and registering of the stack. Two options:
        #     (1) We use SCAMP + SWARP
        #     (2) We use irdr:dithercubemean() + SCAMP/Astrometry.net
        ########################################################################
        # In case of 'science' reduction (ie. 2-skysubtraction passes), 
        ## 'swarp' cannnot be used to build the 1st-coadd due to field distorion.
        # We need to keep the distortion in order to have a good objectMask registering
        # for the 2nd-skysubtraction. In addition, for wide fields, the matching
        # is not perfect, but can be improved dilating (0.5) the objectMask.
        if self.coadd_mode == 'swarp' and (self.red_mode == "quick" or 
                                           self.obs_mode != 'dither'):
            log.info("**** Doing 1st Stack Coaddition (swarp)****")
            # astrometric calibration + coadd
            aw = reduce.astrowarp.AstroWarp(self.m_LAST_FILES, catalog="USNO-B1", 
                        coadded_file=output_file, config_dict=self.config_dict)
            try:
                aw.run(engine=self.config_dict['astrometry']['engine'])
            except Exception,e:
                log.error("Some error while running Astrowarp: %s"%str(e))
                raise e
        else:
            log.info("**** Doing 1st Stack Coaddition (dithercubemean) ****")
            # create input file for dithercubemean
            fo = open(out_dir + '/offsets1.pap', "r")
            fs = open(out_dir + '/stack1.pap', 'w+')
            for line in fo:
                n_line = line.replace(".fits.objs", ".fits") 
                fs.write(n_line)
            fo.close()
            fs.close()
            
            # Dithercubemean
            self.coaddStackImages(out_dir + '/stack1.pap', gainmap, 
                                  out_dir + '/coadd1.fits', 'sum')
            
            # Astrometry of coadded stack
            log.info("**** Doing Astrometric calibration of 1st coadd stack ****")
            if self.config_dict['astrometry']['engine']=='SCAMP':
                try:
                    reduce.astrowarp.doAstrometry(out_dir + '/coadd1.fits', output_file, 
                                          self.config_dict['astrometry']['catalog'], 
                                          do_votable=False)
                except Exception,e:
                    raise Exception("[astrowarp] Cannot solve Astrometry %s"%str(e))
            else:
                try:
                    solved = reduce.solveAstrometry.solveField(out_dir + '/coadd1.fits', 
                                                    out_dir, # self.temp_dir produces collision
                                                    self.config_dict['general']['pix_scale'])

                except Exception,e:
                    raise Exception("[solveAstrometry] Cannot solve Astrometry %s"%str(e))
                else:
                    # Rename the file
                    shutil.move(solved, output_file)
                    
        # ###########################
        # Data quality estimation
        # ###########################
        if self.config_dict['general']['estimate_fwhm']:
            log.info("**** FWHM estimation of coadded_1 result frame ****")
            satur_level = self.config_dict['astrometry']['satur_level']
            pix_scale = self.config_dict['general']['pix_scale']
            
            try:
                pf = datahandler.ClFits(output_file)
                satur_level = int(satur_level) * int(pf.getNcoadds())
            except:
                log.warning("Cannot get NCOADDS. Taken default (=1)")
                satur_level = satur_level
                
            cq = reduce.checkQuality.CheckQuality(output_file, isomin=10.0,
                ellipmax=0.3, edge_x=200, edge_y=200, pixsize=pix_scale, 
                gain=4.15, sat_level=satur_level)
            try:
                (fwhm, std, k, k) = cq.estimateFWHM()
                if fwhm>0 and fwhm<20:
                    log.info("File %s - FWHM = %s (pixels) std= %s"%(output_file, fwhm, std))
                elif fwhm<0 or fwhm>20:
                    log.error("Wrong estimation of FWHM of file %s"%output_file)
                    log.error("Please, review your data.")           
                else:
                    log.error("ERROR: Cannot estimate FWHM of file %s"%output_file)           
            except Exception,e:
                log.error("ERROR: something wrong while computing FWHM")
                raise e

        ########################################################################
        # End of QUICK reduction
        ########################################################################
        if self.obs_mode!='dither' or self.red_mode=="quick":
            log.info("Generated output file ==>%s", output_file)
            log.info("#########################################")
            log.info("##### End of QUICK  data reduction ######")
            log.info("#########################################")
            return output_file
       
        ########################################################################
        # 8 - Create master object mask of the 1st stack coadd.
        ########################################################################
        log.info("************************")
        log.info(" START SECOND PASS      ")
        log.info("************************")

        # Build masterObjMask of the first coadded stack
        log.info("**** Master object mask creation ****")
        obj_mask = self.__createMasterObjMask(output_file, 
                                               out_dir + '/masterObjMask.fits')
        # Prueba
        # Create a object mask for each sky-subtracted file
        obj_mask_skysub = []
        for im in self.m_LAST_FILES:
            res = self.__createMasterObjMask(im,
                                            out_dir + "/" + os.path.basename(im).replace(".fits",".obj.fits"))
            obj_mask_skysub.append(res)
            
        ########################################################################
        # Re-compute dither offsets taking into account the object mask
        ########################################################################
        try:
            # All the files (sky-subtrated and first coadd) are already 
            # astrometrically calibrated.
            temp_list = [obj_mask] + self.m_LAST_FILES
            offset_mat = self.getWCSPointingOffsets(temp_list, 
                                                        out_dir + '/offsets1.pap')
        except Exception,e:
            log.error("Error while computing WCS pointing offsets. "
                "Cannot continue with data reduction...")
            raise e
        
        ########################################################################
        # 8.5 - Re-compute the gainmap taking into account the object mask
        ########################################################################
        # TODO ? In principle, it is not necessary.
        
        ########################################################################
        # 9 - Second Sky subtraction (IRDR) using then OBJECT MASK
        ########################################################################
        log.info("**** Sky subtraction with 2nd object mask ****")
        # 9.1 Compound masked sky file list as input to IRDR::skyfilter()
        fs = open(out_dir + "/skylist2.pap", "w+")
        i = 0
        # (j=1) Due to the offsets were computed with the obj_mask, the first 
        # offset (=0,0) is skipped.
        j = 1 
        ########## PRUEBA ################
        offset_mat.fill(0)
        ### Fin prueba ###
        for file in self.m_rawFiles:
            # In case of whatever T-S-T-S-... sequence, only T frames should be used;
            # however, the second pass of skyfilter (with object mask)
            # has no sense for this type of sequences.
            # if datahandler.ClFits(file).isSky():
            #    continue
            if self.apply_dark_flat == 1 and master_flat != None and master_dark != None:
                line = file.replace(".fits","_D_F.fits") + " " + obj_mask_skysub[i] + " "\
                + str(offset_mat[j][0]) + " " + str(offset_mat[j][1])
            elif self.apply_dark_flat == 1 and master_flat != None:
                line = file.replace(".fits","_F.fits") + " " + obj_mask_skysub[i] + " "\
                + str(offset_mat[j][0]) + " " + str(offset_mat[j][1])
            elif self.apply_dark_flat == 1 and master_dark != None:
                line = file.replace(".fits","_D.fits") + " " + obj_mask_skysub[i] + " "\
                + str(offset_mat[j][0]) + " " + str(offset_mat[j][1])
            else:
                line = file + " " + obj_mask_skysub[i] + " " + str(offset_mat[j][0]) + \
                " " + str(offset_mat[j][1])
            
            fs.write(line + "\n")
            if (self.obs_mode == 'dither_on_off' or 
                self.obs_mode == 'dither_off_on') and i%2:
                j = j + 1
            elif self.obs_mode == 'dither':
                j = j + 1
            i = i + 1
            
        # Close file
        fs.flush()
        os.fsync(fs.fileno())
        fs.close()
        self.m_LAST_FILES = self.skyFilter(out_dir + "/skylist2.pap", gainmap, 
                                           'mask', self.obs_mode)
        
        ########################################################################
        # 9.1 - Remove crosstalk - (only if bright stars are present)    
        ########################################################################
        if self.config_dict['general']['remove_crosstalk']:
            log.info("**** Removing crosstalk ****")
            try:
                res = map(reduce.dxtalk.remove_crosstalk, self.m_LAST_FILES, 
                         [None]*len(self.m_LAST_FILES), [True]*len(self.m_LAST_FILES))
                self.m_LAST_FILES = res
            except Exception,e:
                raise e

        ########################################################################
        # 9.2 - Remove Cosmic Rays -    
        ########################################################################
        if self.config_dict['general']['remove_cosmic_ray']:
            try:
                res = map(reduce.remove_cosmics.remove_cr, self.m_LAST_FILES, 
                            [None]*len(self.m_LAST_FILES), [True]*len(self.m_LAST_FILES),
                            [False]*len(self.m_LAST_FILES))
                self.m_LAST_FILES = res
            except Exception,e:
                raise e
	
    
        #######################################################################
        # 9.3 - Divide by the master flat after sky subtraction ! (see notes above)
        # (the same task as above 4.2) --> HAS NO SENSE !!! only for a test ??? or a.l.a. O2k 
        #######################################################################
        if self.apply_dark_flat==2 and master_flat!=None:
            log.info("**** Applying Flat AFTER sky subtraction ****")
            res = reduce.ApplyDarkFlat(self.m_LAST_FILES, 
                                       None,  
                                       master_flat,
                                       None, 
                                       out_dir)
            self.m_LAST_FILES = res.apply()

	########################################################################
        # Preliminary Astrometric calibration of sky-subtracted frames.
        # ######################################################################
        log.info("**** Preliminary Astrometric calibatrion (2nd sky) ****")
        new_files = []
        for my_file in self.m_LAST_FILES:
            # Run astrometric calibration
            try:
                solved = reduce.solveAstrometry.solveField(my_file, 
                               out_dir, # self.temp_dir produces collision
                               self.config_dict['general']['pix_scale'])
            except Exception,e:
                raise Exception("[solveAstrometry] Cannot solve Astrometry for file: %s \n%s"%(my_file, str(e)))
            else:
                # Rename the file
                out_filename = my_file.replace(".fits", ".ast.fits")
                new_files.append(out_filename)
                shutil.move(solved, out_filename)
                log.debug("New file calibrated: %s"%out_filename)
        
        self.m_LAST_FILES = new_files
        
        ########################################################################
        # 9.4 - LEMON connection - End here for LEMON processing    
        ########################################################################
        
        if self.red_mode=='lemon':
            output_fd, papi_output = \
                tempfile.mkstemp(prefix = out_dir,
                                 suffix = '.papi_output', text = True)
            os.close(output_fd)
            os.unlink(papi_output)
            misc.utils.listToFile(self.m_LAST_FILES, papi_output)
            log.info("End of sequence LEMON-reduction. # %s # files created. ",
                     len(self.m_LAST_FILES))
            return papi_output
        
        #######################################################################
        # 10 - Compute field distortion and make final stack:
        #       1-Remove field distortion from individual images (SCAMP+SWARP)
        #       2-Coaddition of corrected field distortion images (SWARP)
        #       3-Final Astrometric calibration (SCAMP) of the coadded image
        #######################################################################
        if self.coadd_mode == 'swarp':
            log.info("**** Doing Final Stack Coaddition (swarp)****")
            # Build stack using SWARP
            aw = reduce.astrowarp.AstroWarp(self.m_LAST_FILES, catalog="USNO-B1", 
                         coadded_file=output_file, config_dict=self.config_dict,
                         resample=True, subtract_back=True)
            try:
                aw.run(engine=self.config_dict['astrometry']['engine'])
            except Exception,e:
                log.error("Some error while running Astrowarp....")
                raise e
        else:
            log.info("**** Doing Final Stack Coaddition (dithercubemean) ****")
            
            # create input file for dithercubemean
            fo = open(out_dir + '/offsets1.pap', "r")
            fs = open(out_dir + '/stack2.pap', 'w+')
            for line in fo:
                if ".skysub." in line:
                    fs.write(line)
            fo.close()
            fs.close()
            
            # Dithercubemean
            self.coaddStackImages(out_dir + '/stack2.pap', gainmap, 
                                    out_dir + '/coadd2.fits', 'average')
                
            # Final Astrometric calibration
            log.info("**** Doing Astrometric calibration of Final coadded stack ****")
            if self.config_dict['astrometry']['engine']=='SCAMP':
                try:
                    reduce.astrowarp.doAstrometry(out_dir + '/coadd2.fits', output_file, 
                                     self.config_dict['astrometry']['catalog'], 
                                     config_dict=self.config_dict, 
                                     do_votable=False,
                                     resample=True, # it means remove field distorion
                                     subtract_back=True)
                except Exception,e:
                    raise Exception("[astrowarp] Cannot solve Astrometry %s"%str(e))
            else:
                try:
                    solved = reduce.solveAstrometry.solveField(out_dir+'/coadd2.fits', 
                                       out_dir, # self.temp_dir produces collision
                                       self.config_dict['general']['pix_scale'])

                except Exception,e:
                    raise Exception("[solveAstrometry] Cannot solve Astrometry %s"%str(e))
                else:
		  # Rename the file
                  shutil.move(solved, output_file)
        

        log.info("Generated output file ==>%s", output_file)
        
        os.chdir(old_cwd)
        
        log.info("##################################")
        log.info("##### End of data reduction ######")
        log.info("##################################")
        
        return output_file
      
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
 
