#! /usr/bin/env python

# Copyright (c) 2008-2015 IAA-CSIC  - All rights reserved. 
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
# mainGUI (implement main GUI functions for PANIC QL-pipeline)
#
# mainGUI.py
#
################################################################################


    ####################################################################
    #                                                                  #
    #                   This is where it all starts.                   #
    #         The body of this program is at the end of this file      #
    #                                                                  #
    ####################################################################

# system modules
import locale
locale.setlocale(locale.LC_ALL, '')
locale.setlocale(locale.LC_NUMERIC, 'C')
print "LC_NUMERIC =", locale.getlocale(locale.LC_NUMERIC) 

import matplotlib
# Next is needed in order to avoid a crash/deadlock when running 
# pyraf graphs and matplotlib.pyplot graphs
# For 'TkAgg' backend (default) produces the crash.
matplotlib.use('QT4Agg')

import sys
import os
import os.path
import fnmatch
import shutil
import time
import math
import tempfile
from optparse import OptionParser
import datetime
import dircache

# Tell PyRAF to skip all graphics initialization and run in terminal-only mode (=1).
# Otherwise (=0) we will get annoying warning messages (such as "could not open
# XWindow display" or "No graphics display available for this session") when
# working at a remote terminal or at a terminal without any X Windows support.
# Any tasks which attempt to display graphics will fail, of course, but we are
# not going to make use of any of them, anyway.

# What is check is if 'PYRAF_NO_DISPLAY' in os.environ:, so no definition !!
#os.environ['PYRAF_NO_DISPLAY'] = '0'

# PANIC modules
import reduce
import reduce.calTwFlat
import reduce.calBPM_2
import reduce.calBPM
import reduce.applyDarkFlat
import reduce.checkQuality
import reduce.astrowarp
import reduce.reductionset as RS
import misc.fileUtils
import misc.utils as utils
import misc.mef
import datahandler
import misc.display as display
import astromatic
import photo.photometry

#Log
import misc.paLog
from misc.paLog import log
# Fits
import astropy.io.fits as fits
from astropy import coordinates as coord
from astropy import units as u

# IRAF packages

# When PyRAF is imported, it creates, unless it already exists, a pyraf/
# directory for cache in the current working directory. It also complains that
# "Warning: no login.cl found" if this IRAF file cannot be found either. 
# To avoid these two annoying messages, and do not clutter the filesystem with pyraf/
# directories, if $HOME/iraf/login.cl exists, it is used and pyraf/ directory
# is created there. 
from pyraf import iraf
from iraf import noao
from iraf import mscred


# Multiprocessing
from multiprocessing import Process, Queue

from PyQt4 import QtCore, QtGui
from PyQt4.QtCore import QString
from PyQt4.QtCore import QTimer
from PyQt4.QtGui import QFileDialog
from PyQt4.QtGui import QTreeWidgetItem
from PyQt4.QtGui import QTreeWidgetItemIterator
from PyQt4.QtCore import Qt
from PyQt4.QtGui import QApplication, QCursor
from PyQt4.QtGui import *


import commissioning.runStarfocus as focus
import commissioning.getImageOffsets as off



from misc.version import __version__
#-------------------------------------------------------------------------------
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
    
    #Pickle methods properly, including class methods.
    
    func_name = method.im_func.__name__
    obj = method.im_self
    cls = method.im_class
    return _unpickle_method, (func_name, obj, cls)

def _unpickle_method(func_name, obj, cls):
    
    #Unpickle methods properly, including class methods.
    
    for cls in cls.mro():
        try:
            func = cls.__dict__[func_name]
        except KeyError:
            #pass
            raise
        else:
            break
    return func.__get__(obj, cls)        

import copy_reg 
import types 

copy_reg.pickle(types.MethodType,  
    _pickle_method,  
    _unpickle_method)  

#-------------------------------------------------------------------------------


from PyQt4 import QtCore, QtGui, uic

#compile on-the-fly the .ui file
form_class, base_class = uic.loadUiType('panicQL.ui')


class MainGUI(QtGui.QMainWindow, form_class):
    def __init__(self, source_dir="/tmp/data", output_dir="/tmp/out", 
                 temp_dir="/tmp", config_opts=None, *args):
        super(MainGUI, self).__init__(*args)

        #panicQL.__init__(self)
        
        self.setupUi(self)
        
        # Short-keys
        #self.fileShow.setShortcut(QKeySequence("Ctrl+S"))

        self.config_opts = config_opts
              
        ## Init member variables
        self.timer_dc = None
        # File used to run the iraf.obsutil.starfocus task
        self.focus_tmp_file = os.path.expanduser("~") + "/iraf/focus_seq.txt"
        
        # Init main directories
        if source_dir!= None:
            self.m_default_data_dir = source_dir
        else:      
            self.m_default_data_dir = os.environ['PANIC_DATA']
            
        self._fitsGeirsWritten_ = os.path.expanduser("~") + "/tmp/fitsGeirsWritten"
        if not os.path.exists(self._fitsGeirsWritten_):
            self._fitsGeirsWritten_ = "Error: Cannot find ~/tmp/fitsGeirsWritten"


        self.m_sourcedir =  self.m_default_data_dir 
        self.m_outputdir = output_dir 
        self.m_tempdir = temp_dir
        self._ini_cwd = os.getcwd()
        
        # PAPI_HOME
        try:
            self.papi_home = os.environ['PAPI_HOME']
            if self.papi_home[-1]!='/':
                self.papi_home+='/'
        except Exception,e:
            log.error("Error, variable PAPI_HOME not defined.")
            raise e


        # Create LoggingConsole
        self.logConsole = LoggingConsole(self.textEdit_log, self.textEdit_log_2)
        
        # Check we have R/W access to temp and out directories
        if not os.access(self.m_tempdir,os.R_OK|os.W_OK):
            log.error("Directory %s has not R/W access",self.m_tempdir)
            self.logConsole.error(str(QString("[WARNING] Directory %1 has not R/W access")
                                      .arg(self.m_tempdir)))
        if not os.access(self.m_outputdir,os.R_OK|os.W_OK):
            log.error("Directory %s has not R/W access",self.m_outputdir)
            self.logConsole.error(str(QString("[WARNING] Directory %1 has not R/W access")
                                      .arg(self.m_outputdir)))
            
    
        self.m_frameList_dark = ''
        self.m_frameList_dflat = ''
        self.m_frameList_sflat = ''
        self.m_masterDark = ''
        self.m_masterFlat = ''
        self.m_masterMask = ''
        self.m_masterNLC = ''
        self.m_popup_l_sel = []

        # Stuff to detect end of an observation sequence to know if data reduction could start
        self.curr_sequence = []  # list having the files of the current sequence received
        self.isOTrunning = True
        self.last_filter = ''  # filter name (J,H,Ks, ...) of the last dataframe received
        self.last_ob_id = -1   # Obseving Block ID (unique) of the last dataframe received
        # last image type (DARK, DFLAT, TWFLAT, SCIENCE, ...) of last received image
        # Note: in that variable, we do NOT distinguish DFLAT_ON and D_FLAT_OFF
        self.last_img_type = None 
        self.last_ra = -1
        self.last_dec = -1
        self.last_filename = None  # last filename of the FITS received
        self.MAX_POINT_DIST = 1000 # minimun distance (arcsec) to consider a telescope pointing to a new target
        self.MAX_READ_ERRORS = 5 # maximun number of tries while reading a new detetected FITS files 
        
        # GUI properties
        self.m_listView_first_item_selected = None #QListViewItem
        self.m_listView_item_selected = '' 
        self.m_show_imgs = False
        self.m_processing = False
        self.proc_started = False # QL processing activated (True) or deactivated (False)
        self._proc = None # variable to handle QProcess tasks
        self._process = None # pointer to the current running multiprocessing.Process
        self.read_error_files = {} # a dictionary to track the error while reading/detecting FITS files
        
        
        # Load external calibrations files
        self.ext_calib_dir = self.config_opts['general']['ext_calibration_db']

        # Default run mode
        if self.config_opts['quicklook']['run_mode']=="None": 
            self.comboBox_QL_Mode.setCurrentIndex(0)
        elif self.config_opts['quicklook']['run_mode']=="Lazy": 
            self.comboBox_QL_Mode.setCurrentIndex(1)
        elif self.config_opts['quicklook']['run_mode']=="PreReduction": 
            self.comboBox_QL_Mode.setCurrentIndex(2)
        else: 
            self.comboBox_QL_Mode.setCurrentIndex(0)
        
        # Set default values in any case (OT or Filter)
        self.lineEdit_ra_dec_near_offset.setText(str(self.config_opts['general']['max_ra_dec_offset']))
        self.lineEdit_time_near_offset.setText(str(self.config_opts['general']['max_mjd_diff']))
        self.lineEdit_max_num_files.setText(str(self.config_opts['general']['max_num_files']))
            
        self.group_by = self.config_opts['general']['group_by'].lower()
        if self.group_by=='ot':
            self.checkBox_data_grouping_OT.setCheckState(Qt.Checked)
            self.lineEdit_ra_dec_near_offset.setEnabled(False)
            self.lineEdit_time_near_offset.setEnabled(False)
            self.lineEdit_max_num_files.setEnabled(False)
        else:
            self.checkBox_data_grouping_OT.setCheckState(Qt.Unchecked)
            self.lineEdit_ra_dec_near_offset.setEnabled(True)
            self.lineEdit_time_near_offset.setEnabled(True)
            self.lineEdit_max_num_files.setEnabled(True)    
        
        self.logConsole.info("Welcome to the PANIC QuickLook tool version %s"%__version__)
        self.logConsole.info("Instrument: %s"%self.config_opts['general']['instrument'])
        self.logConsole.info("Group by: %s"%self.group_by)
        

        self.__initializeGUI()
        self.createActions()

        # Init in memory Database
        # -----------------------
        # DataBase for input files (in memory)
        self.inputsDB = None

        # DataBase for output files (in memory)
        self.outputsDB = None
        # Init DBs
        self.__initDBs()
        
        # Data Collectors initialization
        # ------------------------------
        self.file_pattern = str(self.lineEdit_filename_filter.text())
        
        if os.path.basename(self.m_sourcedir) == 'save_CA2.2m.log': 
            s_mode = 'geirs-file'
        elif os.path.basename(self.m_sourcedir) == 'fitsGeirsWritten': 
            s_mode = 'geirs-file2'
        else: 
            s_mode = 'dir'
        
        self.dc = None
        self.dc = datahandler.DataCollector(s_mode, self.m_sourcedir, 
                                            self.file_pattern , 
                                            self.new_file_func)
        # Data collector for output files
        self.dc_outdir = None # Initialized in checkOutDir_slot()
        #datahandler.DataCollector("dir", self.m_outputdir, self.file_pattern, self.new_file_func_out)
        
        ## Task management using Threads (ExecTaskThread)
        ## ----------------------------------------------
        self._task = None                       # Pointer to task/thread that is running 
        #self._task_event = threading.Event()    # Event to synchronize ExecTaskThread and ExecTaskThread
        self._task_info = None
        self._task_info_list = []               # Queue-list where task status is saved
        
        self._task_timer = QTimer( self )
        self.connect( self._task_timer, QtCore.SIGNAL("timeout()"), self.checkLastTask )
        self._task_timer.start(1000)    # 1 second continuous timer
        
        
        ## Task management using Processing Queue management
        ## -------------------------------------------------
        # Create queues
        #freeze_support()
        self._task_queue = Queue()
        self._done_queue = Queue()
        
        # Timer for DoneQueue (tasks already done)
        self._queue_timer_done = QTimer( self )
        self.connect( self._queue_timer_done, QtCore.SIGNAL("timeout()"), 
                      self.checkDoneQueue )
        self._queue_timer_done.start(1000)    # 1 second continuous timer
        
        # Timer for TaskQueue (pending tasks)
        self._queue_timer_todo = QTimer( self )
        self.connect( self._queue_timer_todo, QtCore.SIGNAL("timeout()"), 
                      self.taskRunner )
        self._queue_timer_todo.start(1000)    # 1 second continuous timer
        
        
    def __initializeGUI(self):
        """This method initializes some values in the GUI and in the members variables"""

        ## Some signal-connections
        _fromUtf8 = QtCore.QString.fromUtf8
        QtCore.QObject.connect(self.listView_dataS, 
                               QtCore.SIGNAL(_fromUtf8("itemClicked(QTreeWidgetItem*,int)")), 
                               self.selected_file_slot)
        #QtCore.QObject.connect(self.listView_dataS, 
        #                       QtCore.SIGNAL(_fromUtf8("itemSelectionChanged()")), 
        #                       self.selected_file_slot)
        QtCore.QObject.connect(self.listView_dataS, 
                               QtCore.SIGNAL(_fromUtf8("itemDoubleClicked(QTreeWidgetItem*,int)")), 
                               self.display_slot)
        
        ## Init Panel widgets values
        self.lineEdit_sourceD.setText(self.m_sourcedir)
        self.lineEdit_outputD.setText(self.m_outputdir)
        self.lineEdit_tempD.setText(self.m_tempdir)
        
        ## Init calibration files
        self.lineEdit_masterDark.setText(QString(self.m_masterDark))
        self.lineEdit_masterFlat.setText(QString(self.m_masterFlat))
        self.lineEdit_masterMask.setText(QString(self.m_masterMask))
        self.lineEdit_masterNLC.setText(QString(self.m_masterNLC))
        
        
    def __initDBs(self):
        """
        This method initializes the in memory DBs used for input and output files
        """
        
        # Read the instrument for DB creation
        instrument = self.config_opts['general']['instrument'].lower()
        
        # Create Inputs-DB
        try:
            self.inputsDB = datahandler.dataset.DataSet(None, instrument)
            self.inputsDB.createDB()
        except Exception,e:
            log.error("Error while INPUT data base initialization: \n %s"%str(e))
            raise Exception("Error while INPUT data base initialization")
        
        # Create Outputs-DB
        try:
            self.outputsDB = datahandler.dataset.DataSet(None, instrument)
            self.outputsDB.createDB()
            # Insert/load the external calibration files
            self.load_extenal_calibs(self.ext_calib_dir)
        except Exception,e:
            log.error("Error during OUTPUT data base initialization: \n %s"%str(e))
            raise Exception("Error during OUTPUT data base initialization")    
        
    def load_extenal_calibs(self, calib_dir):
        """
        Load all the calibration files found in the specified path.
        """
        
        log.info("Loading External Calibrations into DB from: %s "%calib_dir)
        
        if os.path.isdir(calib_dir):
            for ifile in dircache.listdir(calib_dir):
                if ifile.endswith(".fits") or ifile.endswith(".fit"):
                    log.debug("Inserting file %s"%(calib_dir + "/" + ifile))
                    self.outputsDB.insert(calib_dir + "/" + ifile)
                    
    def new_file_func_out(self, filename):
        """
        Callback used when a new file is detected in output directory.
        """
                     
        self.new_file_func(filename, fromOutput=True)
            
    def new_file_func(self, filename, fromOutput=False):
        """ 
        Function executed when a new file is detected into the data source 
        dir or into the out_dir.
        """
        
        log.debug("New file (in source or out_dir) notified: %s",filename)
        ####################################################
        ### Check if __last__ file is from source directory
        ### Actually, the file was already detected, but this call is only used 
        ### to update the ListView and avoid updating the ListView each time
        ### a file is detected.
        ### Note that 'filename' was already inserted in the DB; this
        ### extra call is only used to update the ListView.
        if filename.endswith("__last__"):
            log.debug("Updating the ListView")
            self.slot_classFilter()
            return 


        ####################################################
        ### Check if deleted file is from source directory
        if filename.endswith("__deleted__"):
            log.debug("File %s disappeared from source directory. Deleted from DB"%filename.replace("__deleted__",""))
            try:
                if fromOutput: self.outputsDB.delete(filename.replace("__deleted__",""))
                else: self.inputsDB.delete(filename.replace("__deleted__",""))
            except:
                log.error("Some error while deleting file %s"%filename.replace("__deleted__",""))
            self.slot_classFilter()
            return
        
        ######################
        ## Insert into DB
        ######################
        inserted = False
        try:
            if fromOutput: inserted = self.outputsDB.insert(filename)
            else: inserted = self.inputsDB.insert(filename)
        except Exception,e:
            log.error("Error while inserting file %s"%filename)
            self.logConsole.warning("Error inserting file [%s]"%(str(e)))
            raise e

        # If file insertion into DB failed, at least nothing else to do...        
        if not inserted:
            self.logConsole.warning("Error inserting file [%s]"%filename) 
            return
            
        #########################
        # An alternative method to update the ListView; it allows keep the view 
        # filtered.
        #self.slot_classFilter() 
        #########################
        
        ## Update Last frame widget
        self.lineEdit_last_file.setText(str(os.path.basename(filename)))
        self.last_filename = filename

        if fromOutput:
            #nothing else to do ...
            return 

        # Check if end of observing sequence (science or calibration), then 
        # start processing.
        end_seq = False
        seq = []
        seqType = ''
        (end_seq, seq, seqType) = self.checkEndObsSequence(filename)
        if end_seq:
            log.debug("Detected end of observing sequence: [%s]"%(seqType))
            self.logConsole.warning("Detected end of observing sequence [%s]"%(seqType))       
        
        
        # Display new file detected
        if self.getDisplayMode()==1 or self.getDisplayMode()==3:
            display.showFrame(filename)
        
        ##################################################################
        ##  If selected, file or sequence processing will start ....
        #   TODO : !! Not completelly finished !!!
        ##################################################################
        if self.proc_started:
            if self.comboBox_QL_Mode.currentText()=="None":
                return
            elif self.comboBox_QL_Mode.currentText().contains("Pre-reduction"):
                if end_seq: self.processFiles(seq)
            elif self.comboBox_QL_Mode.currentText().contains("Lazy"):
                if end_seq and seqType!="SCIENCE":
                    # Build master calibrations
                    self.processFiles(seq)
                else:
                    self.processLazy(filename)
                return
    
    def processLazy(self, filename):
        """
        Do some operations to the last file detected. It depends on:
        
            - checkBox_subDark_FF_BPM
            - checkBox_subLastFrame
            - checkBox_subSky
            
            TODO
            - checkBox_data_grouping
            - comboBox_AstromCatalog
            - comboBox_pre_skyWindow
            
        """
        
        log.debug("[processLazy] Starting to process the file %s",filename)
        
        (date, ut_time, type, filter, texp, detector_id, 
         run_id, ra, dec, object, mjd) = self.inputsDB.GetFileInfo(filename)
        
        
        # ONLY SCIENCE frames will be pre-processed 
        if type!="SCIENCE":
            return

        self.logConsole.info("[processLazy] Processing file %s "%filename)
        
        # ###########################################################################################
        # According to what options have been selected by the user, we do a processing or other...
        if self.checkBox_subDark_FF_BPM.isChecked():
            try:
                # Look for (last received) calibration files
                mDark, mFlat, mBPM = self.getCalibFor([filename])
                
                # Both master_dark and master_flat are optional
                if mDark or mFlat:
                    # Put into the queue the task to be done
                    func_to_run = reduce.ApplyDarkFlat([filename], mDark, mFlat, 
                                                     mBPM, self.m_outputdir,
                                                     bpm_action='grab') # fix is a heavy process for QL
                    params = ()
                    log.debug("Inserting in queue the task ....")
                    self._task_queue.put([(func_to_run.apply, params)])
                    
                else:
                    self.logConsole.error("[processLazy] Cannot find the appropriate master calibration for file %s"%filename)
                    #QMessageBox.critical(self, 
                    #                     "Error", 
                    #                     "Error, cannot find the master calibration files")
            except Exception, e:
                QMessageBox.critical(self, "Error", "Error while processing file.  %s"%str(e))
                #self.m_processing = False
                #QApplication.restoreOverrideCursor()
                raise e    
        # ##########################################################################################
        elif self.checkBox_subLastFrame.isChecked():
            #Change to working directory
            os.chdir(self.m_tempdir)  # -- required ???
            #Create working thread that process the file
            try:
                ltemp = self.inputsDB.GetFilesT('SCIENCE') # (mjd sorted)
                if len(ltemp)>1:
                    last_file = ltemp[-2] # actually, the last in the list is the current one (filename=ltemp[-1])
                    # Change cursor
                    #QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
                    #self.m_processing = False    # Pause autochecking coming files - ANY MORE REQUIRED ?, now using a mutex in thread !!!!
                    #self._task = self.mathOp
                    #thread = reduce.ExecTaskThread(self._task, 
                    #                               self._task_info_list,  
                    #                               [filename, last_file],
                    #                               '-',None)
                    #thread.start()
                    
                    # Put into the queue the task to be done
                    func_to_run = mathOp
                    _suffix = "_" + os.path.basename(last_file)
                    out_filename = self.m_outputdir + "/" + os.path.basename(filename).replace(".fits", _suffix) 
                    params = ([filename, last_file], '-', 
                                         out_filename, self.m_tempdir)
                    self._task_queue.put([(func_to_run, params)])
                    
                else:
                    log.debug("Cannot find the previous file to subtract by")
            except Exception,e:
                QMessageBox.critical(self, "Error", "Error while processing file.  %s"%str(e))
                #self.m_processing = False
                #QApplication.restoreOverrideCursor()
                raise e
        # ##########################################################################################
        elif self.checkBox_subSky.isChecked():
            try:
                self.subtract_nearSky_slot(True)
            except Exception,e:
                self.m_processing = False # ANY MORE REQUIRED ?
                raise e    
        
            
    def getCalibFor(self, sci_obj_list):
        """
        Given a list of frames belonging to a observing sequence for 
        a given object (star, galaxy, whatever),return the most recently created 
        calibration files (master dark,flat,bpm) in order to reduce the sequence.
        The search of the calibration files is done, firstly in the local DB, but
        if no results, then in the external DB if it was provided.
          
        Returns
        -------
          
        Returns 3 calibration files (dark, flat, bpm); If more than one master
        were found, the most recently created (according to MJD) is returned.
        If some master were not found, None is returned.
    
        Notes
        -----
        This function is also equally implemented in ReductionSet 
        """
        
        log.debug("Looking for calibration files into DB")
        
        master_dark = [] # we'll get a list of master dark candidates
        master_flat = [] # we'll get a list of master flat candidates
        master_bpm = []
        
        obj_frame = datahandler.ClFits(sci_obj_list[0])
        # We take as sample, the first frame in the list, but all frames must
        # have the same features (expT,filter,ncoadd, readout-mode, ...)
        expTime = obj_frame.expTime()
        filter = obj_frame.getFilter()
        
        #self.inputsDB.ListDataSet()
        #self.outputsDB.ListDataSet()
        
        # DARK - Do NOT require equal EXPTIME Master Dark ???
        # First, look for a DARK_MODEL, then MASTER_DARK
        master_dark = self.inputsDB.GetFilesT('MASTER_DARK_MODEL', -1) 
        if len(master_dark)==0 and self.outputsDB!=None:
            master_dark = self.outputsDB.GetFilesT('MASTER_DARK_MODEL', -1)
            if len(master_dark)==0:
                master_dark = self.outputsDB.GetFilesT('MASTER_DARK', expTime)
        
        # FLATS - Do NOT require equal EXPTIME, but FILTER
        master_flat = self.inputsDB.GetFilesT('MASTER_DOME_FLAT', -1, filter)
        if len(master_flat)==0:
            master_flat = self.inputsDB.GetFilesT('MASTER_TW_FLAT', -1, filter)
        if len(master_flat)==0 and self.outputsDB!=None:
            master_flat = self.outputsDB.GetFilesT('MASTER_DOME_FLAT', -1, filter)
            if len(master_flat)==0:
                master_flat = self.outputsDB.GetFilesT('MASTER_TW_FLAT', -1, filter)

        # BPM: it is read from config file
        #if self.config_opts['bpm']['mode']!='none':
        master_bpm.append(self.config_opts['bpm']['bpm_file'])
        
        """
        master_bpm = self.inputsDB.GetFilesT('MASTER_BPM')
        if len(master_bpm)==0 and self.outputsDB!=None:
            master_bpm = self.outputsDB.GetFilesT('MASTER_BPM')
        """

        log.debug("Master Darks found %s", master_dark)
        log.debug("Master Flats found %s", master_flat)
        log.debug("Master BPMs  found %s", master_bpm)
        
        # Return the most recently created (according to MJD order)
        if len(master_dark)>0: 
            r_dark = master_dark[-1]
            log.debug("First DARK candidate: %s"%r_dark)            
            r_dark = self.getBestShapedFrame(master_dark, sci_obj_list[0])
            log.debug("Final DARK candidate: %s"%r_dark)            
        else: 
            r_dark = None
        
        if len(master_flat)>0:
            r_flat = master_flat[-1]
            log.debug("First FLAT candidate: %s"%r_flat)            
            r_flat = self.getBestShapedFrame(master_flat, sci_obj_list[0])
            log.debug("Final FLAT candidate: %s"%r_flat)            
        else: 
            r_flat = None
            
        if len(master_bpm)>0:
            r_bpm = master_bpm[-1]
            log.debug("First BPM candidate: %s"%r_bpm)            
            r_bpm = self.getBestShapedFrame(master_bpm, sci_obj_list[0])
            log.debug("Final BPM candidate: %s"%r_bpm)            
        else: 
            r_bpm = None
        
        return r_dark, r_flat, r_bpm
        
    
    def getBestShapedFrame(self, framelist, src_frame):
        """
        Given a list of frames (calibrations) sorted by MJD, return the frame that
        has the same number of extension (MEF) than src_frame and with the 
        same shape (dimensions).

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

    #####################################################
    ### SLOTS ###########################################
    #####################################################
    
    def pushB_sel_updateCalibs_slot(self):
        """
        """
        
        self._update_master_calibrations()
        
    def setDataSourceDir_slot(self):
        """
        Slot function called when the "Input Dir" button is clicked
        """

        source = QFileDialog.getExistingDirectory( self, 
                                                   "Choose existing Directory for source files", 
                                                   self.m_default_data_dir, 
                                                   QFileDialog.ShowDirsOnly
                                                   | QFileDialog.DontResolveSymlinks)
        # source can only be a directory
        
        ###source = QFileDialog.getOpenFileNames( "Source data log (*.log)", 
        ###                                       self.m_default_data_dir, self, 
        ###                                       "Source Dialog","select source")
        ## NOTE: 'source' can be a file or a directory
        
        if (not source):
            return
        else:
            dir = str(source) # required when QFileDialog.getExistingDirectory
            ###dir = str(source[0]) # required when QFileDialog.getOpenFileNames
            if dir==self.m_outputdir:
                self.logConsole.error("Error, Input and Output directories cannot be the same.")
                return
            if (self.m_sourcedir != dir and self.m_outputdir!=dir):
                self.logConsole.info("+Source : " + dir)
                self.lineEdit_sourceD.setText(dir)
                self.m_sourcedir = str(dir)
                self.checkBox_geirsFile.setChecked(False)
                self.checkBox_currentNight.setChecked(False)
                ## Create DataCollector for a path     
                self.file_pattern = str(self.lineEdit_filename_filter.text())
                if self.dc!=None: del self.dc
                if os.path.isfile(dir) and os.path.basename(dir)=="save_CA2.2m.log":
                    self.dc = datahandler.DataCollector("geirs-file", str(dir), 
                                                        self.file_pattern , 
                                                        self.new_file_func)  
                elif os.path.isfile(dir) and os.path.basename(dir)=="fitsfiles.corrected":
                    self.dc = datahandler.DataCollector("geirs-file2", str(dir), 
                                                        self.file_pattern , 
                                                        self.new_file_func)  
                elif os.path.isdir(dir):
                    self.dc = datahandler.DataCollector("dir", str(dir), 
                                                        self.file_pattern , 
                                                        self.new_file_func)
                ## Activate the autochecking of new files
                self.checkBox_autocheck.setChecked(True)
                ##Create QTimer for the data collector
                if self.timer_dc !=None and self.timer_dc.isActive():
                    # Already created in a former user action 
                    # Stop it to avoid collision
                    self.timer_dc.stop()
                
                # Check for new files just now !
                QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
                self.dc.check()
                # Restore the pointer-icon                    
                QApplication.restoreOverrideCursor()                   
                
                # Activate DataCollector timer                 
                if self.timer_dc!=None:
                    self.timer_dc.setSingleShot(False)
                    self.timer_dc.start(1500)
                else:
                    self.timer_dc = QTimer( self )
                    self.connect( self.timer_dc, 
                                 QtCore.SIGNAL("timeout()"), 
                                 self.checkFunc )
                    self.timer_dc.setSingleShot(False)
                    self.timer_dc.start(1500) ## 1 seconds continuous timer
                    
            else:
                #The same dir, nothing to do
                pass
    
    def setGEIRS_Input_slot(self):
        """Called when check-button for Input GEIRS-file is clicked"""
        
        ## Activate or deactivate the autochecking of new files
        if self.checkBox_geirsFile.isChecked():
            self.lineEdit_sourceD.setText(self._fitsGeirsWritten_)
            self.m_sourcedir = self._fitsGeirsWritten_
            if self.dc!=None: del self.dc
            self.dc = datahandler.DataCollector("dir", self.m_sourcedir, 
                                                       self.file_pattern , 
                                                       self.new_file_func)
        else:
            self.lineEdit_sourceD.setText("")
            self.m_sourcedir = ""

    def setCurrentNight_Input_slot(self):
        """Called when check-button for Input checkBox_currentNight is clicked"""
        
        ## Activate or deactivate the autochecking of new files
        if self.checkBox_currentNight.isChecked():
            if datetime.datetime.utcnow().hour > 0 and datetime.datetime.utcnow().hour < 8:
                currentDate = (datetime.datetime.utcnow()-datetime.timedelta(days=1)).isoformat().split('T')[0]
            else:
                currentDate = datetime.datetime.today().isoformat().split('T')[0]
            self.m_sourcedir = "/data1/PANIC/" + currentDate
            self.lineEdit_sourceD.setText(self.m_sourcedir)
            if self.dc != None: del self.dc
            self.dc = datahandler.DataCollector("dir", self.m_sourcedir, 
                                                       self.file_pattern , 
                                                       self.new_file_func)
        else:
            self.lineEdit_sourceD.setText("")
            self.m_sourcedir = ""
            
    def autocheck_slot(self):
        """Called when check-button for Input dir is clicked"""
        
        ## Activate or deactivate the autochecking of new files
        if self.checkBox_autocheck.isChecked():
            if self.timer_dc !=None and not self.timer_dc.isActive():
                # Already created in a former user action 
                self.timer_dc.setSingleShot(False)
                self.timer_dc.start(1500)
            else:
                ##Create QTimer for the data collector
                self.timer_dc = QTimer( self )
                self.connect( self.timer_dc, QtCore.SIGNAL("timeout()"), 
                             self.checkFunc )
                self.timer_dc.setSingleShot(False)
                self.timer_dc.start(1500) ## 1,5 seconds continuous timer
        else:
            self.checkBox_outDir_autocheck.setChecked(False)
            #Stop the DataCollector timer
            self.timer_dc.stop()
            print "*** Stopped DataCollector ***"
            
    def checkOutDir_slot(self):
        """Called when check-button for Output dir is clicked"""
        
        ## Activate or deactivate the autochecking of new files
        if self.checkBox_outDir_autocheck.isChecked():
            if self.dc_outdir==None:
                self.dc_outdir = datahandler.DataCollector("dir", 
                                                           self.m_outputdir, 
                                                           self.file_pattern, 
                                                           self.new_file_func_out)
            # Source dir check is required 
            if not self.checkBox_autocheck.isChecked():
                self.checkBox_autocheck.setChecked(True)
                self.autocheck_slot()
                
    def checkDoneQueue(self):
        """
        Check Queue of done tasks launched with Process in taskRunner.
        Function called periodically (every 1 sec).
        """
        if not self._done_queue.empty():
            log.debug("Something new in the DoneQueue !")
            try:
                r = self._done_queue.get()
                self.logConsole.debug("[checkDoneQueue] Task finished")
                log.info("Got from Queue: %s"%str(r))
                
                if r!=None:
                    if type(r)==type(Exception()):
                        self.logConsole.error(str(r))
                        self.logConsole.error("No processing results obtained")
                    elif type(r)==type(list()):
                        if len(r)==0 :
                            self.logConsole.info(str(QString("No value returned")))
                            #QMessageBox.information(self, "Info", "No value returned")
                        else:
                            str_list = ""
                            #print "FILES CREATED=",self._task_info._return
                            #display.showFrame(r) #_return is a file list
                            for i_file in r:
                                if i_file.endswith(".fits"):
                                    if self.getDisplayMode()>=2: 
                                        display.showFrame(i_file)
                                elif i_file.endswith(".pdf"):
                                    # Todo
                                    log.debug("PDF display not yet implemented.")
                                    continue
                                #display.showFrame(file)
                                str_list+=" +File: " + str(i_file)+"\n"
                                #!!! keep up-date the out DB for future calibrations !!!
                                # Because some science sequences could need the
                                # master calibration created by a former reduction,
                                # and only if apply_master_dark flat is activated,
                                # the last produced file is inserted into the output DB
                                # However, in order to avoid twice insert into
                                # outputDB (although I think it should not be a
                                # problem), if the checkBox for the outputs is 
                                # activated on the GUI, the DB insertion will be 
                                # done there (I hope), and not here !
                                if not self.checkBox_outDir_autocheck.isChecked():   
                                    log.debug("Updating DB...")
                                    self.outputsDB.insert(i_file)
                                #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                            self.logConsole.debug(str(QString("%1 file/s created : \n %2")
                                                      .arg(len(r))
                                                      .arg(str(str_list))))
                    elif type(r)==type("") and os.path.isfile(r):
                        self.logConsole.debug(str(QString("New file %1 created ")
                                                  .arg(r)))
                        if r.endswith(".fits"):
                            if self.getDisplayMode()>=2: display.showFrame(r)
                        # Keep updated the out-DB for future calibrations
                        # See comments above
                        if not self.checkBox_outDir_autocheck.isChecked():
                            log.debug("Updating DB...")
                            self.outputsDB.insert(r)
                        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    else:
                        # Cannot identify the type of the results...for sure
                        # something was wrong...
                        self.logConsole.error("No processing results obtained")
                else:
                    self.logConsole.warning("Nothing returned; processing maybe failed !")
            except Exception,e:
                raise Exception("Error while checking_task_info_list: %s"%str(e))
            finally:
                # Anyway, restore cursor
                QApplication.restoreOverrideCursor()
                self.m_processing = False
                # Return to the previus working directory
                os.chdir(self._ini_cwd)
                # Good moment to update the master calibration files in widgets   
                self._update_master_calibrations()
                
                    
    def checkLastTask(self):
        """
        Receiver function signaled by QTimer 'self._task_timer' object.
        Funtion called periodically (every 1 sec) to check last task results
        launched with the ExecTaskThread function.
        """
        
        if len(self._task_info_list)>0:
            self.logConsole.debug("Task finished")
            try:
                self._task_info = self._task_info_list.pop()
                if self._task_info._exit_status == 0: # EXIT_SUCCESS, all was OK
                    self.logConsole.info("Process successfully finished  !")
                    if self._task_info._return!=None:
                        if type(self._task_info._return)==type(list()):
                            if len(self._task_info._return)==0 :
                                self.logConsole.info(str(QString("No value returned")))
                                #QMessageBox.information(self, "Info", "No value returned")
                            else:
                                str_list = ""
                                #print "FILES CREATED=",self._task_info._return
                                if self.getDisplayMode()>=2: 
                                    display.showFrame(self._task_info._return) #_return is a file list
                                for file in self._task_info._return:
                                    #display.showFrame(file)
                                    str_list+=str(file)+"\n"
                                    #!!! keep up-date the out DB for future calibrations !!!
                                    # Because some science sequences could need the
                                    # master calibration created by a former reduction,
                                    # and only if apply_master_dark flat is activated,
                                    # the last produced file is inserted into the output DB
                                    # However, in order to avoid twice insert into
                                    # outputDB (although I think it should not be a
                                    # problem), if the checkBox for the outputs is 
                                    # activated on the GUI, the DB insertion will be 
                                    # done there (I hope), and not here !
                                    if not self.checkBox_outDir_autocheck.isChecked():
                                        log.debug("Updating DB...")
                                        if (os.path.isfile(file) and
                                            file.endswith(".fits") or 
                                            file.endswith(".fit")):
                                            self.outputsDB.insert(file)
                                    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                self.logConsole.debug(str(QString("%1 files created: \n %2")
                                                          .arg(len(self._task_info._return))
                                                          .arg(str(str_list))))
                                QMessageBox.information(self,"Info", 
                                                        QString("%1 files created: \n %2")
                                                        .arg(len(self._task_info._return))
                                                        .arg(str(str_list)))
                        elif type(self._task_info._return)==type("") and \
                            os.path.isfile(self._task_info._return):
                            self.logConsole.debug(str(QString("New file %1 created.")
                                                      .arg(self._task_info._return)))
                            if self._task_info._return.endswith(".fits"):
                                if self.getDisplayMode()>=2:
                                    display.showFrame(self._task_info._return)
                            # Keep updated the out-DB for future calibrations
                            # See comments above
                            if not self.checkBox_outDir_autocheck.isChecked():
                                log.debug("Updating DB...")
                                ret_file = self._task_info._return
                                if (os.path.isfile(ret_file) and
                                    ret_file.endswith(".fits") or 
                                    ret_file.endswith(".fit")):
                                    self.outputsDB.insert(self._task_info._return)
                            #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        else:
                            # Cannot identify the type of the results...for sure
                            # something was wrong...
                            # Anycase, we show it !
                            self.logConsole.info(str(QString("Value returned : %1")
                                                      .arg(str(self._task_info._return))))
                            #self.logConsole.error("No processing results obtained !")
                    else:
                        self.logConsole.info("Nothing returned !")
                else:
                    self.logConsole.error(str(QString("Sequence processing failed \n %1")
                                              .arg(str(self._task_info._exc))))
                    QMessageBox.critical(self, 
                                         "Error", 
                                         "Error while running task. "+ str(self._task_info._exc))
            except Exception,e:
                raise Exception("Error while checking_task_info_list: %s"%str(e))
            finally:
                #Anyway, restore cursor
                QApplication.restoreOverrideCursor()
                self.m_processing = False
                # Return to the previus working directory
                os.chdir(self._ini_cwd)
                #Good moment to update the master calibration files    
                self._update_master_calibrations()

            
        
    def _update_master_calibrations(self):
        """
        Query the **outputsDB** to update the master calibrations files with the 
        last calibration files received, and then update:
            
            self.m_masterDark
            self.m_masterFlat
            self.m_masterMask
            self.m_masterNLC
            
        and display them on the "Calibrations" view Panel.

        However, for the QL reduction (Lazy or Pre-reduction) we query for newer
        and properly Calibrations that fit the SCI properties (TEXT, FILTER, etc).

        So, this "_update_master_calibrations" only has sense to know the last 
        calibrations found in the DB. 
        """
        
        log.debug("Updating master calibration files received")
        
        master_dark = [] # we'll get a list of master dark candidates
        master_flat = [] # we'll get a list of master flat candidates
        master_bpm = None 
        master_nlc = None
        filter = 'ANY'
        
        # DARK
        # Do NOT require equal EXPTIME Master Dark - only get the last one
        master_dark = self.outputsDB.GetFilesT('MASTER_DARK', -1) 
        
        # FLATS (for any filter)
        master_flat = self.outputsDB.GetFilesT('MASTER_DOME_FLAT', -1, filter)
        if master_flat==[]:
            master_flat = self.outputsDB.GetFilesT('MASTER_TW_FLAT', -1, filter)
        
        # BPM                
        # master_bpm = self.outputsDB.GetFilesT('MASTER_BPM')
        #if self.config_opts['bpm']['mode']!='none':
        master_bpm  = self.config_opts['bpm']['bpm_file']

        # NLC                
        # master_nlc = self.outputsDB.GetFilesT('MASTER_NLC')
        #if self.config_opts['nonlinearity']['apply']!=False:
        master_nlc  = self.config_opts['nonlinearity']['model_lir']
            
        """
        log.debug("Master Darks found %s", master_dark)
        log.debug("Master Flats found %s", master_flat)
        log.debug("Master BPMs  found %s", master_bpm)
        """
        # Return the most recently created (according to MJD order)
        if len(master_dark)>0: 
            self.m_masterDark = master_dark[-1]
            self.lineEdit_masterDark.setText(QString(self.m_masterDark))
            texp = self.outputsDB.GetFileInfo(self.m_masterDark)[4]
            self.lineEdit_masterDark_texp.setText(QString(texp))
        #else: self.m_masterDark = None
        if len(master_flat)>0: 
            self.m_masterFlat = master_flat[-1]
            self.lineEdit_masterFlat.setText(QString(self.m_masterFlat))
            filter = self.outputsDB.GetFileInfo(self.m_masterFlat)[3]
            self.lineEdit_masterFlat_Filter.setText(QString(filter))
        #else: self.m_masterFlat = None
        if master_bpm: 
            self.m_masterMask = master_bpm
            self.lineEdit_masterMask.setText(QString(self.m_masterMask))
        if master_nlc: 
            self.m_masterNLC = master_nlc
            self.lineEdit_masterNLC.setText(QString(self.m_masterNLC))
            self.lineEdit_readout_mode.setText(QString("unknown"))
        #else: self.m_masterMask = None
                

        
         
    def checkEndObsSequence(self, filename):
        """
        Check if the given filename is the end of an observing sequence (calib or science),
        if it is, the method returns True and a list having the files which belong to the list
        and it means the reduction could start.
        Otherwise, False will be returned (and the unfinished currenct list) and
        no reduction can still be done.
        
        Notes
        -----
        It even works for calibration frame sequences (dark, flats, ...)

        Returns
        -------
        A triplet as follow:
            - True if EOS was found, otherwise False
            - if True, the list of files of the sequence
            - if True, the type of sequence (DARK, FLAT, SCIENCE, etc ...)
            
        """
        
        endSeq = False
        retSeq = []
        typeSeq = ''
        
        # Read the FITS file
        fits = datahandler.ClFits(filename, check_integrity=False)
        #only for debug !!
        log.info("Current_FILTER= %s, Last_FILTER=%s, Current_OB_ID=%s, Last_OB_ID=%s",
                 fits.getFilter(), self.last_filter, fits.getOBId(), self.last_ob_id )
        #############################################
        # Based on the meta-data provided by the OT
        # Option based on number of expositions (NOEXP) in the pattern and the 
        # exposition number (EXPNO) of the current frame; It even works for 
        # calibration frames sequences (dark, flats, ...)
        # Added checking of type mismatch of sequence files
        log.info("EXPNO= %s, NOEXPO= %s", fits.getExpNo(), fits.getNoExp())
        if fits.isFromOT():
            log.debug("Checking OT keywords...")
            self.curr_sequence.append(filename)
            # General case for OT observations
            if fits.getExpNo()==1:
                #Start of sequence
                self.curr_sequence = [filename]
                endSeq, retSeq, typeSeq = False, self.curr_sequence, fits.getType()
                #due to start of a new sequnece, we need to update 'last_' values here
                #because some values (last_img_type) are checked below
                self.last_ra = fits.ra
                self.last_dec = fits.dec
                self.last_filter = fits.getFilter()
                self.last_ob_id = fits.getOBId()
                self.last_img_type = fits.getType()
                log.debug("Start OS-0 %s detected : %s"%(typeSeq, str(retSeq)))
                self.logConsole.warning("Start of observing sequence [%s]"%(typeSeq))
                
            if fits.getExpNo()==fits.getNoExp() and fits.getExpNo()!=-1 \
            and (fits.getType()==self.last_img_type or self.last_img_type==None \
                or fits.isSky() or self.last_img_type=='SKY'):
                # End of sequence
                retSeq = self.curr_sequence[:] # very important !! Lists are mutable objects !
                self.curr_sequence = []
                endSeq, retSeq, typeSeq = True, retSeq, fits.getType()
                log.debug("EOS-0 %s detected : %s"%(typeSeq, str(retSeq)))
            else:
                if self.last_img_type!=None and \
                fits.getType()!=self.last_img_type \
                and not fits.isDomeFlat() and not (fits.isSky() or self.last_img_type=='SKY'):
                    log.debug("last_type=%s  new_type=%s"%(self.last_img_type, fits.getType()))
                    # skip/remove the current/last file, type mismatch !
                    self.curr_sequence = self.curr_sequence[:-1]
                    log.error("Detected file type mismatch. File %s skipped in sequence !"%(filename))
                # Mid of sequence, continue adding file
                endSeq, retSeq, typeSeq = False, self.curr_sequence, fits.getType()
        # ############################################
        # We suppose data is obtained using GEIRS+MIDAS_scripts observations
        # POINT_NO(=OB_ID), DITH_NO, EXPO_NO keyword will be checked for sequece detection
        # TODO: TBC !!! I don't know if it will work with calibration sequences.
        else:
            log.debug("Checking GEIRS keywords...")
            if self.last_ob_id==-1: # first time
                # Start of sequence 
                self.curr_sequence.append(filename)
                log.debug("Adding file to sequence: %s",str(self.curr_sequence))
                endSeq, retSeq, typeSeq = False, self.curr_sequence, fits.getType()
            elif fits.getOBId()!=self.last_ob_id \
                or fits.getFilter()!= self.last_filter \
                or fits.getType()!=self.last_img_type:
                retSeq = self.curr_sequence[:] 
                # End of sequence, and then reset the sequence list
                self.curr_sequence = [filename]
                endSeq, retSeq, typeSeq = True, retSeq, fits.getType()
                self.last_ob_id = -1 # re-init
                log.debug("EOS-1 %s detected : %s"%(typeSeq, str(retSeq)))
            else:
                # Check RA,Dec coordinates for distance
                ra_point_distance = self.last_ra-fits.ra
                dec_point_distance = self.last_dec-fits.dec
                dist = math.sqrt((ra_point_distance*ra_point_distance)+(dec_point_distance*dec_point_distance))
                if dist>self.MAX_POINT_DIST:
                    retSeq = self.curr_sequence[:] # very important !! Lists are mutable objects ! 
                    # End of sequence, then reset the sequence list
                    self.curr_sequence = [filename]
                    endSeq,retSeq,typeSeq = True,retSeq,fits.getType()
                    self.last_ob_id = -1 # re-init
                    log.debug("EOS-2 %s detected : %s"%(typeSeq, str(self.retSeq)))
                else:
                    # Mid of sequence, continue adding file
                    self.curr_sequence.append(filename)
                    log.debug("Adding file to sequence: %s",str(self.curr_sequence))
                    endSeq, retSeq, typeSeq = False, self.curr_sequence, fits.getType()
        
        # And finally, before return, update 'last'_values
        self.last_ra = fits.ra
        self.last_dec = fits.dec
        self.last_filter = fits.getFilter()
        self.last_ob_id = fits.getOBId()
        self.last_img_type = fits.getType()
        
        return endSeq,retSeq,typeSeq
                
        """
        # next option for sequence detection is based on FILTER and OB_ID
        if self.isOTrunning:
            if self.last_ob_id==-1: # first time
                self.last_filter=fits.getFilter()
                self.last_ob_id=fits.getOBId()
                self.curr_sequence.append(filename)
                return False,self.curr_sequence
                
            if self.last_filter!=fits.getFilter() or self.last_ob_id!=fits.getOBId():
                self.last_filter=fits.getFilter()
                self.last_ob_id=fits.getOBId()
                seq=self.curr_sequence
                #reset the sequence list
                self.curr_sequence=[filename]
                return True,seq
            else:
                self.curr_sequence.append(filename)
                return False,self.curr_sequence
        """
     
    def checkFunc(self):
        """
        Receiver function signaled by QTimer 'self.timer_dc' object.
        Funtion called periodically to check for new files.
        """

        if True:
            self.dc.check()
            if self.checkBox_outDir_autocheck.isChecked():
                self.dc_outdir.check()


    def getDisplayMode(self):
        """
        Read the 'comboBox_show_imgs' and return the option selected:

        """

        if self.comboBox_show_imgs.currentText().contains("None"):
            return 0 
        elif self.comboBox_show_imgs.currentText().contains("Only new files"):
            return 1
        elif self.comboBox_show_imgs.currentText().contains("Only results"):
            return 2
        elif self.comboBox_show_imgs.currentText().contains("All"):
            return 3

    def data_grouping_slot(self):
        """
        Slot called when the 'checkBox_data_grouping_OT' is clicked.
        If checked, then data grouping will be done using OB_ID, OB_PAT, 
        FILTER keywords (or whatever we decide at the moment).
        If not checked, tha data grouping will be done using RA,Dec values.
        """
        
        if self.checkBox_data_grouping_OT.isChecked():
            self.group_by = 'ot'
            self.lineEdit_ra_dec_near_offset.setEnabled(False)
            self.lineEdit_time_near_offset.setEnabled(False)
            self.lineEdit_max_num_files.setEnabled(False)
        else:
            self.group_by = 'filter'
            self.lineEdit_ra_dec_near_offset.setEnabled(True)
            self.lineEdit_time_near_offset.setEnabled(True)
            self.lineEdit_max_num_files.setEnabled(True)      
                  
    def setOutputDir_slot(self):
        """Select Ouput Directory for processing products"""
       
        dir = QFileDialog.getExistingDirectory( self, 
                                                "Choose existing Directory for output files", 
                                                self.m_outputdir, 
                                                QFileDialog.ShowDirsOnly
                                                | QFileDialog.DontResolveSymlinks) 
        
        
        # CAREFUL: <sourcedir> must be DISTINT  to <outputdir>, infinite loop
        if dir and self.m_outputdir!=str(dir) and self.m_sourcedir!=str(dir):
            self.lineEdit_outputD.setText(dir)
            self.m_outputdir=str(dir)
            self.logConsole.info("+Output dir : " + self.m_outputdir)
            
            ##Create DataCollector for a path     
            self.file_pattern = str(self.lineEdit_filename_filter.text())
            if os.path.isdir(self.m_outputdir):
                self.dc_outdir = datahandler.DataCollector("dir", 
                                                           self.m_outputdir, 
                                                           self.file_pattern , 
                                                           self.new_file_func_out)
            
            ## Activate the autochecking of new files
            self.checkBox_outDir_autocheck.setChecked(True)
                    
        else:
            #The same dir, nothing to do
            pass
        
    def setTempDir_slot(self):
        """
        Select Temp Directory for temporal files while processing
        """
        
        dir = QFileDialog.getExistingDirectory( self, 
                                                "Choose existing Directory for tmp files", 
                                                self.m_tempdir, 
                                                QFileDialog.ShowDirsOnly
                                                | QFileDialog.DontResolveSymlinks)
         
        # CAREFUL: <sourcedir> must be DISTINT  to <tempdir>, infinite loop
        if dir and self.m_sourcedir!=dir:
            self.lineEdit_tempD.setText(dir)
            self.m_tempdir = str(dir)

    def setDarks_slot(self):
        """
        Select master dark to use in data reduction
        """
        
        filelist = QFileDialog.getOpenFileNames( self,
                                                "Select master dark file",
                                                self.m_default_data_dir,
                                                "FITS files (*.fit*)")
        
        
        print "FL=", str(filelist.first())
        a = []
        for fs in filelist:
            a.append(str(fs))
            
        self.m_frameList_dark = a
        
        if not filelist.isEmpty():
            for filename in filelist:
                self.listBox_darks.insertItem(0, str(filename))

    def setDFlats_slot(self):
        """
        Select Master dome flat to be used
        """
                            
        filenames = QFileDialog.getOpenFileNames( self,
                                                "Select Master dome flat to use",
                                                self.m_default_data_dir,
                                                "FITS files (*.fit*)")
        if not filenames.isEmpty():
            for file in filenames:
                self.listBox_domeF.insertItem(0, str(file))  

    def clear_mainlist_slot(self):
        """
        Remove all files from ListView, DataCollector and DB (source and OUTS !).
        """
        
        self.checkBox_geirsFile.setChecked(False)
        self.checkBox_currentNight.setChecked(False)
        self.checkBox_autocheck.setChecked(False)
        self.checkBox_outDir_autocheck.setChecked(False)
        
        # Deactivate timers
        self.autocheck_slot()
        
        self.listView_dataS.clear()
        if self.dc!=None: 
            self.dc.Clear()
        if self.dc_outdir!=None: 
            self.dc_outdir.Clear()

        self.inputsDB.clearDB()
        self.outputsDB.clearDB()
        
        # But keep the external calibration files on outputsDB
        self.load_extenal_calibs(self.ext_calib_dir)
        
    def add_slot(self):
        """
        Add a new file to the main list panel, but not to the DataCollector 
        list (dc), so it might be already inside. 
        Note: see the adding of "__last__" suffix.
        """     
        filenames = QFileDialog.getOpenFileNames( self,
                                                "Select one or more files to add",
                                                self.m_default_data_dir,
                                                "FITS files (*.fit*)")

        if self.comboBox_classFilter.currentText()=="OUTS":
            ifFromOut = True
        else:
            ifFromOut = False

        if not filenames.isEmpty():
            for ifile in filenames:
                self.new_file_func(str(ifile), ifFromOut)
        
            # Due to the trick to avoid the overload of the QL when a full
            # directory is loaded (see datacollector::findNewFiles) we 
            # used an extra call with the suffix "__last__" to update the View
            self.new_file_func(str(ifile+"__last__"), ifFromOut)

    def del_slot(self):
        """     
        Delete the current selected file from the main list view panel, 
        but we do not remote from the DataCollector neither file system.
        """
        
        self.m_popup_l_sel = []
        listViewItem = QTreeWidgetItemIterator(self.listView_dataS,
                                               QTreeWidgetItemIterator.Selected)
        while listViewItem.value():
            fileName = str(listViewItem.value().text(0))
            indx = self.listView_dataS.indexOfTopLevelItem(listViewItem.value())
            self.listView_dataS.takeTopLevelItem(indx)
            #self.dc.Clear(fileName)
            if self.comboBox_classFilter.currentText()=="OUTS":
                self.outputsDB.delete( fileName )
            else:
                self.inputsDB.delete( fileName )
            #listViewItem +=1 # because of remove, we do not need an increment
          
        #self.inputsDB.ListDataSetNames()
    
    def display_slot(self):
        
        if self.m_listView_item_selected:
            display.showFrame(self.m_listView_item_selected)
        else:
            print "Error, no file selected !"
            
    def filename_filter_slot(self):
        """ Modify filename filter for the data collector"""

        new_filter, ok = QInputDialog.getText(self, "New filter",
                                              "Enter filename filter:")
        if ok and not new_filter.isEmpty():
            self.lineEdit_filename_filter.setText( str(new_filter) )
            self.dc.SetFileFilter( str(new_filter) )
            self.dc_outdir.SetFileFilter( str(new_filter) )

    def slot_classFilter(self):
        """ Filter files on main ListView"""
        
        if self.comboBox_classFilter.currentText()=="GROUP":
            self.listView_dataS.clear()
            sequences = []
            seq_types = []
            
            #Look for sequences
            ra_dec_near_offset = self.lineEdit_ra_dec_near_offset.text().toInt()[0] #arcsec
            time_near_offset = self.lineEdit_time_near_offset.text().toInt()[0]/86400.0 #day units
            max_number_files = self.lineEdit_max_num_files.text().toInt()[0]
            
            sequences, seq_types = self.inputsDB.GetSequences(self.group_by,
                                                    max_mjd_diff=time_near_offset,
                                                    max_ra_dec_diff=ra_dec_near_offset, 
                                                    max_nfiles=max_number_files) 
            
            #
            # Look for un-groupped files and build a group/sequence with them
            #
            temp = set([])
            for lista in sequences:
                temp = temp.union(set(lista))
            un_groupped = set(self.inputsDB.GetFiles()) - temp
            """    
            for file in self.inputsDB.GetFiles():
                for lista in sequences:
                    if file in lista:
                        temp.add(file)
                #temp.union( set([file for lista in sequences if file in lista]) )
            un_groupped = set(self.inputsDB.GetFiles()) - temp
            """
            if len(un_groupped)>0:
                sequences.append(list(un_groupped))
                seq_types.append("UNKNOWN")
            
            k = 0
            for seq in sequences:
                elem = QTreeWidgetItem( self.listView_dataS )
                #elem.setText(0, "TYPE="+str(seq_types[k]))
                #elem.setText(0, "OB_ID="+str(seq[0])+" ** OB_PAT="+str(seq[1])+" ** FILTER="+str(seq[2]) + " ** #imgs="+str(len(fileList[k])) ) # OB_ID + OB_PAT + FILTER
                for file in seq:
                    (date, ut_time, type, filter, texp, detector_id, run_id, 
                     ra, dec, object, mjd) = self.inputsDB.GetFileInfo(file)
		    #nCoadd = fits.getval(file, "NCOADDS", ext=0)
                    if file==seq[0]:
			nCoadd = fits.getval(file, "NCOADDS", ext=0)
			nObject = fits.getval(file,"OBJECT", ext=0) 
                        #the first time, fill the "tittle" of the group 
                        elem.setText(0, "TYPE="+str(seq_types[k]) + 
                                     "  ** FILTER=" + str(filter) + 
                                     "  ** TEXP=" + str(texp) + 
                                     "  ** OBJ=" + str(nObject) + 
                                     " ** #imgs=" + str(len(seq)))
                        
                    e_child = QTreeWidgetItem(elem)
                    e_child.setText (0, str(file))
                    e_child.setText (1, str(type))
                    e_child.setText (2, str(filter))
                    e_child.setText (3, str(texp))
                    e_child.setText (4, str(date)+"::"+str(ut_time))
                    e_child.setText (5, str(object))
                    if dec<0: sign = -1;
                    else: sign = 1;
                    c = coord.ICRS(ra=ra, dec=dec ,unit=(u.degree, u.degree))
                    str_ra = "%02d:%02d:%04.1f"%(c.ra.hms[0], c.ra.hms[1], c.ra.hms[2])
                    str_dec = "%02d:%02d:%02.0f"%(c.dec.dms[0], c.dec.dms[1]*sign, c.dec.dms[2]*sign)
                    e_child.setText (6, str(str_ra))
                    e_child.setText (7, str(str_dec))
                k+=1
        else:
            #############################################    
            ## No grouping, only filtering by given type
            #############################################
            self.listView_dataS.clear()
            db = None
            if str(self.comboBox_classFilter.currentText())=="INPUTS":
                fileList = self.inputsDB.GetFilesT("ANY")
                db = self.inputsDB
            elif str(self.comboBox_classFilter.currentText())=="ALL":
                fileList = self.inputsDB.GetFilesT("ANY")
                db = self.inputsDB
                #db.ListDataSet()
            elif str(self.comboBox_classFilter.currentText())=="OUTS":
                fileList = self.outputsDB.GetFilesT("ANY")
                db = self.outputsDB
                #db.ListDataSet()
            elif str(self.comboBox_classFilter.currentText())=="SKY_FLAT":
                fileList = self.inputsDB.GetFilesT("SKY_FLAT")
                db = self.inputsDB
            elif str(self.comboBox_classFilter.currentText())=="FOCUS":
                fileList = self.inputsDB.GetFilesT("FOCUS")
                db = self.inputsDB
            elif str(self.comboBox_classFilter.currentText())=="MASTERS":
                fileList = self.outputsDB.GetFilesT("MASTER_DARK")
                fileList += self.outputsDB.GetFilesT("MASTER_DARK_MODEL")
                fileList += self.outputsDB.GetFilesT("MASTER_TW_FLAT")
                fileList += self.outputsDB.GetFilesT("MASTER_DOME_FLAT")
                fileList += self.outputsDB.GetFilesT("MASTER_SKY_FLAT")
                db = self.outputsDB
            else:
                type = str(self.comboBox_classFilter.currentText())
                fileList = self.inputsDB.GetFilesT(type)
                db = self.inputsDB

            elem = None
            for file in fileList:
                elem = QTreeWidgetItem( self.listView_dataS )
                (date, ut_time, type, filter, texp, detector_id, run_id, ra, 
                 dec, object, mjd)= db.GetFileInfo(file)
                elem.setText (0, str(file))
                elem.setText (1, str(type))
                elem.setText (2, str(filter))
                elem.setText (3, str(texp))
                elem.setText (4, str(date)+"::"+str(ut_time))
                elem.setText (5, str(object))
                c = coord.ICRS(ra=ra, dec=dec ,unit=(u.degree, u.degree))
                str_ra = "%02d:%02d:%04.1f"%(c.ra.hms[0], c.ra.hms[1], c.ra.hms[2])
                str_dec = "%02d:%02d:%02.0f"%(c.dec.dms[0], c.dec.dms[1], c.dec.dms[2])
                elem.setText (6, str(str_ra))
                elem.setText (7, str(str_dec))
            
            # In addition, if "ALL" is selected, we show the OUTS as well
            if str(self.comboBox_classFilter.currentText())=="ALL":
                elem = None
                fileList = self.outputsDB.GetFilesT("ANY")
                db = self.outputsDB
                for file in fileList:
                    elem = QTreeWidgetItem( self.listView_dataS )
                    (date, ut_time, type, filter, texp, detector_id, run_id, ra, 
                     dec, object, mjd)=db.GetFileInfo(file)
                    elem.setText (0, str(file))
                    elem.setText (1, str(type))
                    elem.setText (2, str(filter))
                    elem.setText (3, str(texp))
                    elem.setText (4, str(date)+"::"+str(ut_time))
                    elem.setText (5, str(object))
                    c = coord.ICRS(ra=ra, dec=dec ,unit=(u.degree, u.degree))
                    str_ra =  "%02d:%02d:%04.1f"%(c.ra.hms[0], c.ra.hms[1], c.ra.hms[2])
                    str_dec = "%02d:%02d:%02.0f"%(c.dec.dms[0], c.dec.dms[1], c.dec.dms[2])
                    elem.setText (6, str(str_ra))
                    elem.setText (7, str(str_dec))
            
            if elem:
                self.listView_dataS.setCurrentItem(elem)
                
                      
#########################################################################
###### Pop-Up ###########################################################
#########################################################################
    def createActions(self):
        """
        Create actions used in the Pop-Up menu (contextMenuEvent).
        Shurtcut: Note that, for letters, the case used in the specification 
        string does not matter.
        """
        
        #### Main Pop-up actions
        self.dispAct = QtGui.QAction("&Display image", self,
            shortcut="Ctrl+D",
            statusTip="Display current selected image", 
            triggered=self.display_slot)
        
        self.copyAct = QtGui.QAction("&Copy files to clipboard", self,
            shortcut="Ctrl+C",
            statusTip="Copy current selected files to clipboard", 
            triggered=self.copy_sel_files_slot)
        
        self.ditherAct = QtGui.QAction("&Show Dither pattern", self,
            shortcut="Ctrl+p",
            statusTip="Show plot with dithter pattern of selected files.", 
            triggered=self.show_dither_pattern_slot)
        
        self.toTextFileAct = QtGui.QAction("&Copy files to text file", self,
            shortcut="Shift+T",
            statusTip="Copy current selected files text file.", 
            triggered=self.copy_files2TextFile_slot)

        self.mDarkAct = QtGui.QAction("&Build Master Dark", self,
            shortcut="Ctrl+M",
            statusTip="Build master dark with selected files", 
            triggered=self.createMasterDark_slot)
        
        self.mDFlatAct = QtGui.QAction("&Build Master Dome Flat", self,
            shortcut="Ctrl+F",
            statusTip="Build master dome flat with selected files", 
            triggered=self.createMasterDFlat_slot)
        
        self.mDTwFlatAct = QtGui.QAction("&Build Master Twlight Flat", self,
            shortcut="Ctrl+T",
            statusTip="Build master twlight flat with selected files", 
            triggered=self.createMasterTwFlat_slot)
        
        self.mGainMapAct = QtGui.QAction("&Build Gain Map", self,
            shortcut="Ctrl+G",
            statusTip="Build master gain map with selected files", 
            triggered=self.createGainMap_slot)
        
        self.mBPMAct = QtGui.QAction("&Build BPM", self,
            shortcut="Ctrl+B",
            statusTip="Build Bad Pixel Map with selected files", 
            triggered=self.createBPM_slot)

        self.mBPM_appAct = QtGui.QAction("&Apply and Show BPM", self,
            shortcut="Shift+B",
            statusTip="Apply BPM and show the masked pixels in a temp file.", 
            triggered=self.applyBPM_slot)

        self.mDF_appAct = QtGui.QAction("&Apply Dark & FlatField & BPM", self,
            shortcut="Ctrl+A",
            statusTip="Apply Dark, Flat and BPM", 
            triggered=self.applyDarkFlat)

        self.mFocusEval = QtGui.QAction("&Focus evaluation", self,
            shortcut="Ctrl+F",
            statusTip="Run a telescope focus evaluation of a focus serie.", 
            triggered=self.focus_eval)

        self.subOwnSkyAct = QtGui.QAction("Subtract own-sky", self,
            shortcut="Ctrl+S",
            statusTip="Subtract own-sky", 
            triggered=self.subtract_ownSky_slot)
        
        self.subNearSkyAct = QtGui.QAction("Subtract near-Sky", self,
            shortcut="Ctrl+N",
            statusTip="Subtract near sky", 
            triggered=self.subtract_nearSky_slot)

        self.quickRedAct = QtGui.QAction("Quick-Reduction", self,
            shortcut="Ctrl+Q",
            statusTip="Run quick data reduction of selected files", 
            triggered=self.do_quick_reduction_slot)
        
        self.astroAct = QtGui.QAction("Astrometric Calib.", self,
            shortcut="Shift+A",
            statusTip="Run astrometric calibration of selected file", 
            triggered=self.do_raw_astrometry)
        
        self.photoAct = QtGui.QAction("Photometric Calib.", self,
            shortcut="Shift+P",
            statusTip="Run photometric calibration of selected file", 
            triggered=self.do_raw_photometry)
        
        self.statsAct = QtGui.QAction("Statistics", self,
            shortcut="Shift+S",
            statusTip="Get statistincs of selected file", 
            triggered=self.show_stats_slot)
        
        self.fwhmAct = QtGui.QAction("FWHM mean estimation", self,
            shortcut="Ctrl+W",
            statusTip="Get FWHM selected file", 
            triggered=self.fwhm_estimation_slot)
        
        self.bckgroundAct = QtGui.QAction("Background estimation", self,
            shortcut="Shift+B",
            statusTip="Get image background of selected file", 
            triggered=self.background_estimation_slot)
            
        # Sub-menu Math actions
        self.subAct = QtGui.QAction("&Subtract Images", self,
            shortcut="Ctrl+-",
            statusTip="Subtract selected files", 
            triggered=self.subtractFrames_slot)
        
        self.combAct = QtGui.QAction("Combine Images (median)", self,
            shortcut="Ctrl+*",
            statusTip="Median combine selected files", 
            triggered=self.combFrames_slot)
        
        self.sumAct = QtGui.QAction("Sum Images", self,
            shortcut="Ctrl++",
            statusTip="Sum selected files", 
            triggered=self.sumFrames_slot)
        
        self.divAct = QtGui.QAction("Divide Images", self,
            shortcut="Ctrl+/",
            statusTip="Divide selected files", 
            triggered=self.divideFrames_slot)

        # Sub-menu FITS actions       
        self.mef2singleAct = QtGui.QAction("MEF2Single", self,
            #shortcut="Ctrl+M",
            statusTip="Convert a MEF to a single file", 
            triggered=self.MEF2Single_slot)

        self.single2mefAct = QtGui.QAction("Single2MEF", self,
            #shortcut="Ctrl+",
            statusTip="Convert a single file to MEF file", 
            triggered=self.Single2MEF_slot)


        self.splitMEFAct = QtGui.QAction("Split MEF", self,
            #shortcut="Ctrl+s",
            statusTip="Split MEF file into 4-single files.", 
            triggered=self.splitMEF_slot)
        
        self.splitSingleAct = QtGui.QAction("Split Single", self,
            #shortcut="Ctrl+j",
            statusTip="Split single file into 4-single files.", 
            triggered=self.splitSingle_slot)
        
        self.collapseCubeAct = QtGui.QAction("Collapse Cube", self,
            #shortcut="Ctrl+j",
            statusTip="Collapse FITS with a cube of N layers.", 
            triggered=self.collapseCube_slot)
        
        # Group-menu actions
        self.redSeqAct = QtGui.QAction("Reduce Sequence", self,
            shortcut="Ctrl+R",
            statusTip="Reduce the selected sequence", 
            triggered=self.reduceSequence_slot)
        
        self.copyGAct = QtGui.QAction("&Copy files to clipboard", self,
            shortcut="Shift+C",
            statusTip="Copy current selected files to clipboard", 
            triggered=self.copy_files_slot)
        
        
    def listView_popup_slot(self, listItem, mouse_pointer):
        """
        deprecated with PyQt4!
        """
        pass

    def contextMenuEvent(self, event):
        """
        Show pop-up menu
        """
        
        #### Get items selected in the ListView
        self.m_popup_l_sel = []
        listViewItem = QTreeWidgetItemIterator(self.listView_dataS, 
                                               QTreeWidgetItemIterator.Selected)
        father = None
        while listViewItem.value():
            self.m_popup_l_sel.append(str(listViewItem.value().text(0)))
            # if a parent of the group, no popup to show 
            if (self.comboBox_classFilter.currentText()=="GROUP" and \
                listViewItem.value().childCount()!=0): # is a parent of a group
                father = listViewItem.value()
                self.m_listView_first_item_selected = father
                break
            listViewItem +=1

        if len(self.m_popup_l_sel)<=0:
            return
            
        
        # we have selected a group father
        if father: # we have selected a group father
            #### Create the 'Group' Popup menu
            popUpMenu = QMenu()
            popUpMenu.addAction(self.redSeqAct)
            popUpMenu.addAction(self.copyGAct)
               
            group_files = []
            for ic in range(father.childCount()):
                group_files.append(father.child(ic).text(0))
                
        else:
            #### Create the 'Files' Popup menu 
            popUpMenu = QtGui.QMenu("QL popup menu", self)
            popUpMenu.addAction(self.dispAct)
            popUpMenu.addAction(self.copyAct)
            popUpMenu.addAction(self.toTextFileAct)
            popUpMenu.addAction(self.ditherAct)
            popUpMenu.addSeparator()
            popUpMenu.addAction(self.mDarkAct)
            popUpMenu.addAction(self.mDFlatAct)
            popUpMenu.addAction(self.mDTwFlatAct)
            popUpMenu.addAction(self.mGainMapAct)
            popUpMenu.addAction(self.mBPMAct)
            popUpMenu.addSeparator()
            popUpMenu.addAction(self.mDF_appAct)
            popUpMenu.addAction(self.mBPM_appAct)
            popUpMenu.addAction(self.mFocusEval)
            popUpMenu.addAction(self.subOwnSkyAct)
            popUpMenu.addAction(self.subNearSkyAct)
            popUpMenu.addAction(self.quickRedAct)
            popUpMenu.addSeparator()
            popUpMenu.addAction(self.astroAct)
            popUpMenu.addAction(self.photoAct)
            popUpMenu.addSeparator()
            popUpMenu.addAction(self.statsAct)
            popUpMenu.addAction(self.fwhmAct)
            popUpMenu.addAction(self.bckgroundAct)
            
            newAction = popUpMenu.addAction("Math")
            subMenu_arith = QtGui.QMenu("Popup Submenu", self)
            subMenu_arith.addAction(self.sumAct)
            subMenu_arith.addAction(self.subAct)
            subMenu_arith.addAction(self.divAct)
            subMenu_arith.addAction(self.combAct)
            newAction.setMenu(subMenu_arith)
            
            newAction = popUpMenu.addAction("FITS")
            subMenu_fits = QtGui.QMenu("Popup Submenu", self)
            subMenu_fits.addAction(self.mef2singleAct)
            subMenu_fits.addAction(self.single2mefAct)
            subMenu_fits.addAction(self.splitMEFAct)
            subMenu_fits.addAction(self.splitSingleAct)
            subMenu_fits.addAction(self.collapseCubeAct)
            
            newAction.setMenu(subMenu_fits)
            

            ## Disable some menu items depeding of the number of item selected 
            ## in the list view.
            ## Non mentioned actions are always Enabled (e.g.: self.copyGAct,
            ##  self.splitAct, self.joinAct)
            ## Then, only is needed to set non enabled action.
            
            # Restore default values, because an former call could have changed
            # the state.
            self.dispAct.setEnabled(True)
            self.ditherAct.setEnabled(True)
            self.mDarkAct.setEnabled(True)
            self.mDFlatAct.setEnabled(True)
            self.mDTwFlatAct.setEnabled(True)
            self.mGainMapAct.setEnabled(True)
            self.mBPMAct.setEnabled(True)
            self.mDF_appAct.setEnabled(True)
            self.mBPM_appAct.setEnabled(True)
            self.mFocusEval.setEnabled(True)
            self.subOwnSkyAct.setEnabled(True)
            self.subNearSkyAct.setEnabled(True)
            self.quickRedAct.setEnabled(True)
            self.astroAct.setEnabled(True)
            self.photoAct.setEnabled(True)
            self.fwhmAct.setEnabled(True)
            self.bckgroundAct.setEnabled(True)
            self.subAct.setEnabled(True)
            self.sumAct.setEnabled(True)
            self.divAct.setEnabled(True)
            self.combAct.setEnabled(True)
            
            if len(self.m_popup_l_sel)==1:
                #print "#SEL_1=",len(self.m_popup_l_sel)
                #self.dispAct.setEnabled(False)
                self.mDarkAct.setEnabled(False)
                self.mDFlatAct.setEnabled(False)
                self.mDTwFlatAct.setEnabled(False)
                self.mGainMapAct.setEnabled(False)
                self.mBPMAct.setEnabled(False)
                self.mFocusEval.setEnabled(False)
                #self.subNearSkyAct.setEnabled(False)
                self.subAct.setEnabled(False)
                self.sumAct.setEnabled(False)
                self.divAct.setEnabled(False)
                self.combAct.setEnabled(False)
                
            if len(self.m_popup_l_sel)>1:
                #print "#SEL_2",len(self.m_popup_l_sel)
                self.dispAct.setEnabled(False)
                self.subOwnSkyAct.setEnabled(False)
                self.subNearSkyAct.setEnabled(False)
                self.astroAct.setEnabled(False)
                self.photoAct.setEnabled(False)
                #self.fwhmAct.setEnabled(False)
                self.bckgroundAct.setEnabled(False)
                self.mBPM_appAct.setEnabled(False)

            if len(self.m_popup_l_sel)!=2:
                self.subAct.setEnabled(False)
                # It allowed to sum/combine 2 or more images
                #self.sumAct.setEnabled(False) 
                self.divAct.setEnabled(False)
                
            if len(self.m_popup_l_sel)<5:
                pass
                
        ## Finally, execute the popup
        popUpMenu.exec_(event.globalPos())
        
    def reduceSequence_slot(self):
        """
        Run the data reduction of the current grouped-sequence selected in the
        ListView when the GROUP view is selected.
        """
        
        group_files = []
        #child = self.m_listView_first_item_selected.firstChild()
        #while child:
        #    group_files.append(str(child.text(0)))
        #    child = child.nextSibling()
        father = self.m_listView_first_item_selected
        for ic in range(father.childCount()):
                group_files.append(str(father.child(ic).text(0)))
        
        try:
            self.processFiles(group_files, interactive=True)
        except Exception,e:
            QMessageBox.critical(self, "Error", 
                                 "Error while group data reduction: \n%s"%str(e))
            raise e

    def copy_files_slot(self):
        """
        Copy the file list fullnames belonging to the current group to the 
        clipboard in order to allow copy&paste
        """
        
        text = ""
        father = self.m_listView_first_item_selected
        for ic in range(father.childCount()):
            text = text + str(father.child(ic).text(0)) + "\n"
        
        clipboard = QApplication.clipboard()
        clipboard.setText(text)
    
    def copy_sel_files_slot(self):
        """
        Copy the selected file list fullnames to the clipboard in order to 
        allow copy&paste.
        """
        
        text = ""
        for filename in self.m_popup_l_sel:
            text = text + str(filename) + "\n"
        
        clipboard = QApplication.clipboard()
        clipboard.setText(text)    
    
    def copy_files2TextFile_slot(self):
        """
        Copy selected files into a text file. In addition, the name of the text 
        file is copied to the clipboard.
        """

        textFile = QFileDialog.getOpenFileName( self,
                                            "Select output text file",
                                            self.m_tempdir,   
                                            "(*.txt*")
        textFile = str(textFile)
        if str(textFile):
            self.genFileList(self.m_popup_l_sel, textFile)
            clipboard = QApplication.clipboard()
            clipboard.setText(textFile)
            QMessageBox.information(self, "Info", "Text file %s created and copied to clipboard"%(textFile))

    def show_dither_pattern_slot(self):
        """
        Show plot with the dither pattern followed by selected files.
        """
        
        # The sort out of the files by MJD is done 
        # in getWCSPointingOffsets().
        try:    
            offsets = off.getWCSPointingOffsets(self.m_popup_l_sel, "/dev/null")
        except Exception,e:
            msg = "Error, cannot find out the image offsets. %s"%str(e)
            log.error(msg)
            QMessageBox.information(self, "Info", msg)
        else:
            # Draw plot with dither pathern
            off.draw_offsets(offsets, self.config_opts['general']['pix_scale'], 1.0)

    def MEF2Single_slot(self):
        """
        Convert a MEF file with 4-extensions to a single file of 4kx4k.
        """
        
        QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
        
        for file in self.m_popup_l_sel:
            try:
                mef = misc.mef.MEF([file])
                res = mef.doJoin(".join.fits", output_dir=self.m_outputdir)[1]
            except Exception,e:
                log.debug("Cannot convert MEF to Single file %s. Maybe it's not a MEF file", str(e))
                QMessageBox.critical(self, "Error", "Cannot convert MEF to Single file : %s \n Maybe it's not a MEF file"%(file))
            else:
                line = "File generated: %s"%res
                self.logConsole.info(QString(str(line)))
        
        QApplication.restoreOverrideCursor()
                
    def Single2MEF_slot(self):
        """
        Convert a GEIRS single file (4kx4k) to a MEF file with 4-extensions,
        2kx2k each.
        """
        
        QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
        
        for file in self.m_popup_l_sel:
            try:
                mef = misc.mef.MEF([file])
                res = mef.convertGEIRSToMEF(out_dir=self.m_outputdir)[1]
            except Exception,e:
                log.debug("Cannot convert Single to MEF file %s. Maybe it's not a single file", str(e))
                QMessageBox.critical(self, "Error", "Cannot convert Single to MEF file : %s \n Maybe it's not a single file"%(file))
            else:
                line = "File generated: %s"%res
                self.logConsole.info(QString(str(line)))
        
        QApplication.restoreOverrideCursor()
        
    def splitMEF_slot(self):
        """
        Split each MEF selected file from the list view into NEXT separate 
        single FITS files, where NEXT is number of extensions.
        As result, NEXT files should be created
        """
        QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
        
        for file in self.m_popup_l_sel:
            try:
                mef = misc.mef.MEF([file])
                res = mef.doSplit(".Q%02d.fits", out_dir=self.m_outputdir)[1]
            except Exception,e:
                log.debug("Cannot split file %s. Maybe it's not a MEF file", str(e))
                QMessageBox.critical(self, "Error", "Cannot split file : %s \n Maybe it's not a MEF file"%(file))
                #self.logConsole.info(QString(str(line)))
            else:
                line = "File generated: %s"%res
                self.logConsole.info(QString(str(line)))
                
        QApplication.restoreOverrideCursor()
    
    def splitSingle_slot(self):
        """
        Split the Single selected files (4kx4k) from the list view into 4 separate 
        single FITS files or 2kx2k. The 4-detector-frame can be a cube of data.
        As result, 4 files should be created.
        """
        
        QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
        
        for file in self.m_popup_l_sel:
            try:
                mef = misc.mef.MEF([file])
                res = mef.splitGEIRSToSimple(out_dir=self.m_outputdir)[1]
            except Exception,e:
                msg = "Cannot split file %s. Maybe it's not a 4kx4k Single file"%(file, str(e))
                log.debug(msg)
                QMessageBox.critical(self, "Error", msg)
                self.logConsole.info(QString(str(msg)))
            else:
                line = "File generated: %s"%res
                self.logConsole.info(QString(str(line)))
        
        QApplication.restoreOverrideCursor()
    
    def collapseCube_slot(self):
        """
        Collapse a FITS cube of N layers or planes. No recentering or layers 
        is done. FITS can be MEF or SEF.
        """
        
        QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
        
        for file in self.m_popup_l_sel:
            try:
                new_file = misc.collapse.collapse([file], out_dir=self.m_outputdir)[0]
            except Exception,e:
                msg = "Cannot split file %s. Maybe it's not a Cube file"%(file, str(e))
                log.debug(msg)
                QMessageBox.critical(self, "Error", msg)
                self.logConsole.info(QString(str(msg)))
            else:
                line = "File generated: %s"%new_file
                self.logConsole.info(QString(str(line)))
        
        QApplication.restoreOverrideCursor()
    
    def selected_file_slot(self, listItem):
        """
        To know which item is selected 
        """
                
        if listItem:
            if (self.comboBox_classFilter.currentText()=="GROUP" and 
                listItem.childCount()!=0): # it'a a parent
                self.m_listView_item_selected = None # for consistency
                return
            else: self.m_listView_item_selected = str(listItem.text(0))
    
    def imexam_slot(self):
        """
        Imexam the currect filename selected - not used because it fails !!
        """
        
        try:
            #self.m_processing = False    # Pause autochecking coming files - ANY MORE REQUIRED ?, now using a mutex in thread !!!!
            thread = reduce.ExecTaskThread(self.run_imexam, 
                                           self._task_info_list, 
                                           self.m_popup_l_sel[0] )
            thread.start()
        except Exception,e:
            QMessageBox.critical(self, "Error", 
                                 "Error while Imexam with image %s"%(self.m_popup_l_sel[0]))
            raise e
        
    def run_imexam(self, filename):
        """ 
        Run Imexam the currect filename. First we start the DS9 display  
        
        DOES NOT WORK !!! --> try with multiprocessing.Process !
        
        Console error
        ~~~~~~~~~~~~~
        "
        DS9 not running
        display frame (1:) (3): error in background error handler:
        out of stack space (infinite loop?)
        while executing
        "::tcl::Bgerror {out of stack space (infinite loop?)} {-code 1 -level 0 -errorcode NONE -errorinfo {out of stack space (infinite loop?)
        while execu..."
        DS9 not running
        "
        
        """
        
        display.startDisplay()
        try:
            # It's a sync call, i mean, a blocking call that doesn't return 
            # until is finished
            iraf.imexam(filename) 
        except Exception,e:
            log.error("Error while Imexam with image %s",filename)
            raise e
            
        
    ######### End Pup-Up ########################################################    
    
    def fileExit(self):
        #First, stop the DataCollector timer
        if self.checkBox_autocheck.isChecked():
            self.timer_dc.stop()
        
        os.system("killall ds9")
        sys.exit(0)
        
    def fileOpen(self):
        self.setDataSourceDir_slot()
        return
    
    def listbox_dataSource_popup(self):
        """
        Not finished ! To be completed
        """
        #QMessageBox.critical(self,"Error", "Could not save the current document")
        popUpMenu = QMenu()
        popUpMenu.addAction("Create Master Dark")
        popUpMenu.addAction("Create Master Dome-Flat")
        popUpMenu.addAction("Create Master Sky-Flat")
        #actions["act1"].addTo(fileMenu)
        popUpMenu.exec_(QCursor.pos())
        
        return

    def iraf_console_slot(self, isForFocus=False):
        """Open a IRAF console session"""
        
        # To avoid that the 'self.focus_tmp_file' file could be left due
        # some error/abort while running the iraf.startfocus procedure,
        # the file is removed if the call to this routine is not for
        # the focus evaluation with iraf.starfocus.
        if not isForFocus and os.path.exists(self.focus_tmp_file):
	    os.unlink(self.focus_tmp_file)
	    
        os.system("cd $HOME/iraf;/usr/local/bin/xgterm -title IRAF -cr red -ms blue -sb -sl 1000 -geometry 100x30 -bg grey -fg black -e cl &")
        
    def start_ds9_slot(self):
        """
        Start DS9 display, and if exixts, open the currect selected file in 
        the list view.
        """
        
        if self.m_listView_item_selected!="":
            display.showFrame(self.m_listView_item_selected)
        else:
            display.startDisplay()
            
        self.logConsole.info("DS9 launched !")
        
    def start_aladin_slot(self):
        """Start Aladin tool"""

        os.system('echo "load %s ;sync; get vizier(2mass)" |/usr/local/Aladin/Aladin -nobanner &' %(self.m_listView_item_selected))
        self.logConsole.info("Aladin launched !")
        
        # utils.runCmd does not allow launch in background !!

    def pushB_sel_masterDark_slot(self):
        """
        Called to select manually a master dark.
        """
        source = QFileDialog.getOpenFileName( self,
                                             "Select master DARK file", 
                                              self.m_default_data_dir, 
                                              "(*.fit*")

        if str(source):
            fits = datahandler.ClFits(str(source))
            if fits.getType()!='MASTER_DARK' or \
                fits.getType()!='MASTER_DARK_MODEL':
                res = QMessageBox.information(self, "Info", 
                                            QString("Selected frame does not look an MASTER DARK.\n Continue anyway?"), 
                                            QMessageBox.Ok, QMessageBox.Cancel)
                if res==QMessageBox.Cancel:
                    return
                  
            self.lineEdit_masterDark.setText(source)
            self.lineEdit_masterDark_texp.setText(str(fits.expTime()))
            self.m_masterDark = str(source)

    def pushB_sel_masterFlat_slot(self):
        """
        Called to select manually a master flat.
        """
        source = QFileDialog.getOpenFileName( self,
                                              "Select master FLAT file",
                                              self.m_default_data_dir, 
                                              "(*.fit*)")

        if str(source):
            fits = datahandler.ClFits(str(source))
            if fits.getType()!='MASTER_SKY_FLAT' or \
                fits.getType()!='MASTER_DOME_FLAT' or \
                fits.getType()!='MASTER_TW_FLAT':
                res = QMessageBox.information(self, "Info", 
                                            QString("Selected frame does not look an MASTER FLAT.\n Continue anyway?"), 
                                            QMessageBox.Ok, QMessageBox.Cancel)
                if res==QMessageBox.Cancel:
                    return
            self.lineEdit_masterFlat.setText(source)
            self.lineEdit_masterFlat_Filter.setText(fits.getFilter())
            self.m_masterFlat = str(source)
               

    def pushB_sel_masterMask_slot(self):
        """
        Called to select manually a master BPM.
        """
        source = QFileDialog.getOpenFileName( self,
                                            "Select Master BPM",
                                            self.m_default_data_dir, 
                                            "(*.fit*)")

        if str(source):
            self.lineEdit_masterMask.setText(source)
            self.m_masterMask = str(source)                       

    def pushB_sel_masterNLC_slot(self):
        """
        Called to select manually a master NLC.
        """
        source = QFileDialog.getOpenFileName( self,
                                            "Select Master NLC", 
                                            self.m_default_data_dir, 
                                            "(*.fit*)")

        if str(source):
            self.lineEdit_masterNLC.setText(source)
            self.m_masterNLC = str(source)
            
    def pushB_subtract_last2_slot(self):
        """
        Slot called when button 'Subtract-Last2' is clicked.
        
        The two newest frames received (MJD) will be subtracted
        """
        
        #(type , filter) = self.inputsDB.GetFileInfo(self.last_filename)[2:4]
        #if type != 'SCIENCE':
        #only SCIENCE frames are processed
        #pass
        #ltemp = self.inputsDB.GetFilesT('SCIENCE') # (mjd sorted)
        
        ltemp = self.inputsDB.GetFilesT('ANY') # (mjd sorter)
        if len(ltemp)>1:
            #get the  last 2 files (included the current one)
            last2_files = ltemp[-2:]
            if (self.inputsDB.GetFileInfo(last2_files[0])[3] == 
                self.inputsDB.GetFileInfo(last2_files[1])[3]):
                try:
                    #Put into the queue the task to be done
                    func_to_run = mathOp
                    params = (last2_files, "-", "/tmp/sub.fits", self.m_tempdir)
                    self._task_queue.put([(func_to_run, params)])
                except:
                    QMessageBox.critical(self, "Error", "Error while subtracting files")
                    raise
                
            else:
                self.logConsole.error("Last two frames have different FILTER. Cannot subtract each other")
                log.warning("Last two frames have different FILTER. Cannot subtract each other")
                return
        else:
            self.logConsole.info("Not enough files for subtraction")
            log.debug("Not enough files for subtraction")
            return
        
        
    
    def genFileList( self, file_list, outFilename ):
        """ 
        Generate a file 'outFilename' listing the files passed as a 
        python-list in the file_list parameter
        """
        fo = open(outFilename,"w")
        for my_file in file_list:
            fo.write(my_file+"\n")
        fo.close()
                             
    ############################################################################
    ############# PROCESSING STAFF #############################################
    ############################################################################
        
    def subtractFrames_slot(self):
        """
        This method is called to subtract two images selected from the File 
        List View.
        """


        if (len(self.m_popup_l_sel)!=2):
            QMessageBox.critical(self, "Error", "You need to select  2 files")
        else:
            outFilename = QFileDialog.getSaveFileName(self,
                                                      "Choose a filename to save under",
                                                      self.m_outputdir + "/sub.fits", 
                                                      "*.fits")
            if not outFilename.isEmpty():
                try:
                    #Change cursor
                    QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
                    # Pause autochecking coming files - ANY MORE REQUIRED ?, 
                    # now using a mutex in thread !!!!
                    self.m_processing = False    
                    thread = reduce.ExecTaskThread(mathOp, 
                                                   self._task_info_list, 
                                                   self.m_popup_l_sel,'-', 
                                                   str(outFilename))
                    thread.start()
                except:
                    QMessageBox.critical(self, "Error", 
                                         "Error while subtracting files")
                    raise
        
    def combFrames_slot(self):
        """
        This methot is called to **combine** the images selected from the 
        File List View.
        
        Note: The type of combinning operation to the pixels is a median after
        a sigma reject algorith.

        """

        if (len(self.m_popup_l_sel)<2):
            QMessageBox.critical(self, "Error", "You need to select  at least 2 files")
        else:
            outFilename = QFileDialog.getSaveFileName(self,
                                                      "Choose a filename to save under",
                                                      self.m_outputdir+"/sum.fits", 
                                                      "*.fits")
            if not outFilename.isEmpty():
                try:
                    #Change cursor
                    QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
                    self.m_processing = False    
                    # Pause autochecking coming files - ANY MORE REQUIRED ?, 
                    # now using a mutex in thread !!!!
                    thread = reduce.ExecTaskThread(mathOp, 
                                                   self._task_info_list, 
                                                   self.m_popup_l_sel,
                                                   'combine', str(outFilename))
                    thread.start()
                except:
                    QMessageBox.critical(self, "Error", 
                                         "Error while adding files")
                    raise
                
    def sumFrames_slot(self):
        """
        This methot is called to **sum** the images selected from the 
        File List View.

        """

        if (len(self.m_popup_l_sel)<2):
            QMessageBox.critical(self, "Error", "You need to select  at least 2 files")
        else:
            outFilename = QFileDialog.getSaveFileName(self,
                                                      "Choose a filename to save under",
                                                      self.m_outputdir+"/sum.fits", 
                                                      "*.fits")
            if not outFilename.isEmpty():
                try:
                    #Change cursor
                    QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
                    self.m_processing = False    
                    # Pause autochecking coming files - ANY MORE REQUIRED ?, 
                    # now using a mutex in thread !!!!
                    thread = reduce.ExecTaskThread(mathOp, 
                                                   self._task_info_list, 
                                                   self.m_popup_l_sel,
                                                   '+', str(outFilename))
                    thread.start()
                except:
                    QMessageBox.critical(self, "Error", 
                                         "Error while adding files")
                    raise
                
    def divideFrames_slot(self):
        """
        This methot is called to divide two images selected from the File List View
        """

        if (len(self.m_popup_l_sel)!=2):
            QMessageBox.critical(self, "Error", "You need to select  2 files")
        else:
            outFilename = QFileDialog.getSaveFileName(self,
                                                      "Choose a filename to save under",
                                                      self.m_outputdir + "/div.fits", 
                                                      "*.fits")
            if not outFilename.isEmpty():
                try:
                    #Change cursor
                    QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
                    # Pause autochecking coming files - ANY MORE REQUIRED ?, 
                    # now using a mutex in thread !!!!
                    self.m_processing = False
                    thread=reduce.ExecTaskThread(mathOp, 
                                                 self._task_info_list, 
                                                 self.m_popup_l_sel, 
                                                 '/',  str(outFilename))
                    thread.start()
                except:
                    QMessageBox.critical(self, "Error", 
                                         "Error while dividing files")
                    raise

    def createMasterDark_slot(self):
        """
        Called to create a master dark frame.
        """

        if len(self.m_popup_l_sel)<3:
            QMessageBox.information(self,"Info","Not enough frames !")
            return

        #listItems = self.listBox_darks.findItems(QString("*"), 
        #                                    Qt.MatchWrap | Qt.MatchWildcard)
        # Build list
        #l_list = [ str(item.text()) for item in listItems]

        #if len(l_list)<3:
        #    QMessageBox.information(self,"Info","Not enough frames !")
        #    return
        
        #print "LIST =", self.m_popup_l_sel
        outfileName = QFileDialog.getSaveFileName(self,
                                                  "Choose a filename to save under",
                                                  self.m_outputdir + "/master_dark.fits", 
                                                  "*.fits")
        # if not outfileName.isEmpty():
        #     try:
        #         texp_scale = False
        #         QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
        #         self._task = reduce.calDark.MasterDark(self.m_popup_l_sel, 
        #                                                self.m_tempdir, 
        #                                                str(outfileName), 
        #                                                texp_scale)
        #         thread = reduce.ExecTaskThread(self._task.createMaster, 
        #                                        self._task_info_list)
        #         thread.start()
        #     except Exception, e:
        #         QApplication.restoreOverrideCursor() 
        #         QMessageBox.critical(self, "Error", 
        #                              "Error while creating master Dark. \n"+str(e))
        #         raise e
        # else:
        #     pass

    
        if not outfileName.isEmpty():
            try:
                self.processFiles(self.m_popup_l_sel, group_by=self.group_by, 
                                  outfilename=str(outfileName))
            except Exception,e:
                QMessageBox.critical(self, "Error", 
                                 "Error while creating master Dark: \n%s"%str(e))
                raise e
        else:
            pass     

    def createMasterDFlat_slot(self):
        
        """
        Create a master dome flat field from selected frames
        """
        
        #listItems = self.listBox_domeF.findItems(QString("*"), 
        #                                         Qt.MatchWrap | Qt.MatchWildcard)
        # Build list
        #l_list = [ str(item.text()) for item in listItems]

        
        if len(self.m_popup_l_sel)>2:
            outfileName = QFileDialog.getSaveFileName(self,
                                                      "Choose a filename to save under",
                                                      self.m_outputdir + "/master_dflat.fits", 
                                                      "*.fits")
            if not outfileName.isEmpty():
                try:
                    QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
                    self._task = reduce.calDomeFlat.MasterDomeFlat(
                                        self.m_popup_l_sel, 
                                        self.m_tempdir, 
                                        str(outfileName))
                    thread=reduce.ExecTaskThread(self._task.createMaster, 
                                                 self._task_info_list)
                    thread.start()
                except:
                    QApplication.restoreOverrideCursor()
                    QMessageBox.critical(self, "Error", "Error while creating master Dome Flat")
                    raise
        else:
            QMessageBox.information(self,"Info","Error, not enough frames (>2) !")

    def createMasterTwFlat_slot(self):
      
        """
        Create a master twilight flat field from selected frames
        """
        #listItems = self.listBox_SkyF.findItems(QString("*"), 
        #                                         Qt.MatchWrap | Qt.MatchWildcard)
        # Build list
        #l_list = [ str(item.text()) for item in listItems]
        
        if len(self.m_popup_l_sel)>2:
            outfileName = QFileDialog.getSaveFileName(self,
                           "Choose a filename to save under",
                           self.m_outputdir+"/master_twflat.fits", 
                           "*.fits")
            
            if not outfileName.isEmpty():
                # Look for MasterDarks and MasterDarkModel
                dark_model = []
                darks = []
                dark_model = self.outputsDB.GetFilesT('MASTER_DARK_MODEL')
                # If no files, returns an empy list
                darks = self.outputsDB.GetFilesT('MASTER_DARK', -1)
                
                print "DARKS=", darks
                print "DM=",dark_model
                try:
                    QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
                    
                    self._task = reduce.calTwFlat.MasterTwilightFlat(
                        flat_files=self.m_popup_l_sel,
                        master_dark_model=dark_model,
                        master_dark_list=darks,
                        output_filename=str(outfileName))
                    
                    thread = reduce.ExecTaskThread(self._task.createMaster, 
                                                 self._task_info_list)
                    thread.start()
                except Exception,e:
                    QApplication.restoreOverrideCursor()
                    msg = "Error creating master Twilight Flat file: %s"%str(e)
                    log.error(msg)
                    self.logConsole.error(msg)
                    QMessageBox.critical(self, "Error",  msg)
        else:
            QMessageBox.information(self,"Info","Error, not enough frames (>2) !")
    
    def createGainMap_slot(self):

        """ 
        Create a gain-map using the own science files selected on the main list view
        
        TODO: Check the selected images are SCIENCE frames with same filter !!! 
        and expTime, .....
         
        """
        
        if len(self.m_popup_l_sel)<=1:
            QMessageBox.information(self,"Info","Not enough frames !")
            return

        outfileName = QFileDialog.getSaveFileName(self,
                                                  "Choose a filename to save under",
                                                  self.m_outputdir+"/gainmap.fits", 
                                                   "*.fits")
        if not outfileName.isEmpty():
            flat_type = self.inputsDB.GetFileInfo(self.m_popup_l_sel[0])[2]
            if flat_type.count("SCIENCE"):
                try:
                    QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
                    self._task = reduce.calGainMap.SkyGainMap(self.m_popup_l_sel, 
                                                              str(outfileName), 
                                                              None, self.m_tempdir)
                    thread = reduce.ExecTaskThread(self._task.create, self._task_info_list)
                    thread.start()
                except Exception, e:
                    QApplication.restoreOverrideCursor()
                    QMessageBox.critical(self, "Error", "Error while creating Gain Map. "+str(e))
                    raise e
            elif flat_type.count("DOME_FLAT"):
                try:
                    QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
                    self._task = reduce.calGainMap.DomeGainMap(self.m_popup_l_sel, 
                                                              str(outfileName), 
                                                              None)
                    thread = reduce.ExecTaskThread(self._task.create, self._task_info_list)
                    thread.start()
                except Exception, e:
                    QApplication.restoreOverrideCursor()
                    QMessageBox.critical(self, "Error", "Error while creating Gain Map. "+str(e))
                    raise e
            elif flat_type.count("SKY_FLAT"):
                try:
                    QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
                    self._task = reduce.calGainMap.TwlightGainMap(self.m_popup_l_sel,
                                                              self.m_masterDark, 
                                                              str(outfileName), 
                                                              None, self.m_tempdir)
                    thread = reduce.ExecTaskThread(self._task.create, self._task_info_list)
                    thread.start()
                except Exception, e:
                    QApplication.restoreOverrideCursor()
                    QMessageBox.critical(self, "Error", "Error while creating Gain Map. "+str(e))
                    raise e
            else:
                QMessageBox.information(self,"Info","Cannot build GainMap, image type not supported.")
                log.error("Cannot build GainMap, image type <%s> not supported."%str(flat_type))
        else:
            pass                     
                
    def subtract_ownSky_slot(self):
        """ Subtract OWN image sky background using SExtrator tool
            NOTE: this operation support MEF files, but the background image
            created by SExtractor might not have all the keywords required (AR, DEC, INSTRUMENT, ...)
        """
        
        if not self.m_listView_item_selected:
            return
        fits = datahandler.ClFits(self.m_listView_item_selected)
        
        if fits.getType()=='SCIENCE':
            sex_config = self.papi_home + self.config_opts['config_files']['sextractor_conf']
            sex_param = self.papi_home + self.config_opts['config_files']['sextractor_param']
            sex_conv = self.papi_home + self.config_opts['config_files']['sextractor_conv']
            minarea =  self.config_opts['skysub']['mask_minarea']
            threshold = self.config_opts['skysub']['mask_thresh']
            input_file = self.m_listView_item_selected
            out_file = self.m_outputdir+"/"+os.path.basename(input_file.replace(".fits",".skysub.fits"))
        
            #cmd="sex %s -c %s -FITS_UNSIGNED Y -DETECT_MINAREA %s  -DETECT_THRESH %s  -CHECKIMAGE_TYPE -BACKGROUND -CHECKIMAGE_NAME %s" % (input_file, sex_config, str(minarea), str(threshold), out_file)  
            #cmd="sex %s -c %s -DETECT_MINAREA %s -DETECT_THRESH %s -CHECKIMAGE_TYPE -BACKGROUND -CHECKIMAGE_NAME %s -PARAMETERS_NAME %s -FILTER_NAME %s"% (input_file, sex_config, str(minarea), str(threshold), out_file, sex_param, sex_conv )  
            
            """
            ***
                Be sure that SExtractor working directory has R/W in order to create
                temporal files created by SExtractor !!!
            ***
            """
            self.m_processing = True
            
            # Next implementation might not work properly because sex.run() 
            # doesn't return the output file created !
            #
            input_file = self.m_listView_item_selected
            out_file = self.m_outputdir+"/"+os.path.basename(input_file.replace(".fits",".skysub.fits"))
            sex = astromatic.SExtractor()
            #sex.config['CONFIG_FILE']="/disk-a/caha/panic/DEVELOP/PIPELINE/PANIC/trunk/config_files/sex.conf"
            #sex.ext_config['CHECKIMAGE_TYPE'] = "OBJECTS"
            sex.config['CATALOG_TYPE'] = "NONE"
            sex.config['CHECKIMAGE_TYPE'] = "-BACKGROUND"
            sex.config['CHECKIMAGE_NAME'] = out_file
            sex.config['DETECT_THRESH'] = self.config_opts['skysub']['mask_thresh']
            sex.config['DETECT_MINAREA'] = self.config_opts['skysub']['mask_minarea']
            try:
                QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
                updateconfig = True
                clean = False
                sex.run(input_file, updateconfig, clean)
            except Exception,e:
                QMessageBox.critical(self, "Error", "Error while running Sextractor"+str(e))
                raise e
            else:
                if self.getDisplayMode()>=2:
                    display.showFrame(out_file)
            finally:
                QApplication.restoreOverrideCursor()
                self.m_processing = False
        else:
            log.error("Sorry, selected file does not look a science file  ")
            QMessageBox.information(self, "Info", "Sorry, selected file does not look a science file ")                             
                         
    def subtract_nearSky_slot(self, last_files=False):
        """ 
        Subtract nearest sky using skyfiler from IRDR package:

        Parameters
        ----------
            - If last_files==False, it looks for the nearest files to the 
            currently selected file.
            - If last_files==True, the last N-1 files received are used to 
            compute the sky to subtract.  
        
        Returns
        -------
        Nothing or Exception in case of error.

        """
      
        
        # #####################################
        # Look for files in search-radius mode
        # #####################################
        if not last_files:
            # Compute search-radius of nearest frames
            ra_dec_near_offset = self.lineEdit_ra_dec_near_offset.text().toInt()[0]/3600.0
            time_near_offset = self.lineEdit_time_near_offset.text().toInt()[0]/86400.0
            #print "RA_DEC_OFFSET", ra_dec_near_offset
            #print "TIME_NEAR_OFFSET", time_near_offset
          
            # Look for NEAREST science files
            if self.m_listView_item_selected:
                fits = datahandler.ClFits(self.m_listView_item_selected)
                if fits.getType()=='SCIENCE':
                    near_list = self.inputsDB.GetFiles('ANY', 'SCIENCE', -1, 
                                                       fits.filter, fits.mjd, 
                                                       fits.ra, fits.dec,  
                                                       ra_dec_near_offset*2, 
                                                       time_near_offset, runId=0)
                    # For the moment, the minimun number of nearest is >0
                    if len(near_list)==0:
                        QMessageBox.information(self, "Info", "Not enough science frames found")  
                        return
                else:
                    QMessageBox.information(self, "Info", "Selected frame is not a science frame") 
                    return
                       
                #Create master list of nearest frames (ar,dec,mjd) to the current selected science file
                p = 1
                file_n = -1
                my_list = ""
                for file in near_list:
                    my_list = my_list + file + "\n"
                    if (file==self.m_listView_item_selected): 
                        file_n = p
                    else: 
                        p+=1
                
                # Ask user validation
                res = QMessageBox.information(self, "Info", 
                                              QString("Selected near frames are:\n %1")
                                              .arg(my_list), 
                                              QMessageBox.Ok, QMessageBox.Cancel)
                if res==QMessageBox.Cancel:
                    return
                     
        # ##################################    
        # Take the lastest N-files received
        # ##################################
        else:
            (type , filter) = self.inputsDB.GetFileInfo(self.last_filename)[2:4]
            if type != 'SCIENCE':
                #only SCIENCE frames are processed
                pass
            ltemp = self.inputsDB.GetFilesT('SCIENCE',-1, filter) # (mjd sorted)
            if len(ltemp)>4:
                #get the last 5 files (included the current one)
                near_list = ltemp[-5:] # actually, the last in the list is the current one (filename=ltemp[-1])
                file_n = len(near_list) # range used by subtractNearSky is [1-N]
            else:
                self.logConsole.info("Not enough files for sky subtraction (>4)")
                log.debug("Not enough files for sky subtraction (>5)")
                return

        # ##################################    
        # Now, compute the sky-subtraction
        # ##################################
                            
        #Change to working directory
        os.chdir(self.m_tempdir)
        #Change cursor
        QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
        #Create working thread that compute sky-frame
        try:
            self._task = RS.ReductionSet( [str(item) for item in near_list], 
                                          self.m_outputdir,
                                          out_file=self.m_outputdir + "/skysub.fits",
                                          obs_mode="dither", dark=None, 
                                          flat=self.m_masterFlat,
                                          bpm=None, red_mode="quick",
                                          group_by=self.group_by, 
                                          check_data=True, 
                                          config_dict=self.config_opts)
            
            thread = reduce.ExecTaskThread(self._task.subtractNearSky, 
                                           self._task_info_list, 
                                           self._task.rs_filelist, 
                                           file_n # range used by subtractNearSky is [1-N]
                                           )
            thread.start()
        except:
            #Anyway, restore cursor
            # Although it should be restored in checkLastTask, could happend an exception while creating the class RS,
            # thus the ExecTaskThread can't restore the cursor
            QApplication.restoreOverrideCursor() 
            QMessageBox.critical(self, "Error", "Error while subtracting near sky")
            raise 
              
    def createBPM_slot(self):
        """ 
        Create a Bad Pixel Mask from a set of selected files (flats)
        It is run interactively.

        NOTE: not sure if it works ??
        """
                    
        if len(self.m_popup_l_sel)>3:
            outfileName = QFileDialog.getSaveFileName(self,
                                                      "Choose a filename so save under",
                                                      "/tmp/BadPixelMask", 
                                                      "BPM (*.pl)", )
            if not outfileName.isEmpty():
                show = True
                try:
                    self.genFileList(self.m_popup_l_sel, "/tmp/bpm.list")
                    dark = "/tmp/master_dark.fits"
                    lsig = 10
                    hsig = 10
                    QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
                    self._task = reduce.calBPM_2.BadPixelMask( "/tmp/bpm.list", 
                                                               dark, 
                                                               str(outfileName), 
                                                               lsig, hsig, 
                                                               self.m_tempdir)
                    thread = reduce.ExecTaskThread(self._task.create, 
                                                   self._task_info_list)
                    thread.start()
                except:
                    #Restore cursor
                    QApplication.restoreOverrideCursor()
                    QMessageBox.critical(self, "Error", "Not suitable frames to compute BPM.\n You need flat_off and flat_on frames")
                    raise
        else:
            QMessageBox.critical(self, "Error","Error, not suitable frames selected (flat fields)")

    def applyDarkFlat(self):
        """
        Apply to the selected files the master Dark and master Flat found in the
        database. The master Dark and FlatField are searched individualy for each
        selected file.

        The original files are not modified, but new output files with _D_F.fits
        suffix are created.
        """

        # Ask for master calibration 
        msgBox = QMessageBox()
        msgBox.setText("        Apply Dark, Flat-Field and BPM")
        msgBox.setInformativeText("Do you want to <AutoSearch> the calibration files or use <Defaults> ones ?")
        button_autosearch = msgBox.addButton("AutoSearch", QMessageBox.ActionRole)
        button_defaults = msgBox.addButton("Defaults/Select", QMessageBox.ActionRole)
        button_cancel = msgBox.addButton("Cancel", QMessageBox.ActionRole)
        msgBox.setDefaultButton(button_autosearch)
        
        msgBox.exec_()
        
        calibrations = dict(dark='', flat='', bpm='', nlc='')
        auto_search = False
        if msgBox.clickedButton() == button_defaults or msgBox.clickedButton()==button_autosearch:
            if msgBox.clickedButton() == button_defaults and self.checkBox_use_defCalibs.isChecked():
                calibrations['dark'] = self.m_masterDark
                calibrations['flat'] = self.m_masterFlat
                calibrations['bpm'] = self.m_masterMask
                calibrations['nlc'] = self.m_masterNLC
                force_apply = False
                
            elif msgBox.clickedButton()== button_defaults and not self.checkBox_use_defCalibs.isChecked():
                # Select master DARK
                # In this case, a 'simple' file (whatever) could be used as a MASTER_DARK/MASTER_FLAT
                # Nooo, I do not think so !
                force_apply = False 
                filenames = QFileDialog.getOpenFileNames(self,
                                                "Select master DARK to use",
                                                self.m_outputdir,
                                                "FITS files (*.fit*)")
                if not filenames.isEmpty():
                    calibrations['dark'] = str(filenames[0])

                # Select master FLAT
                filenames = QFileDialog.getOpenFileNames(self,
                                                "Select master FLAT to use",
                                                self.m_outputdir,
                                                "FITS files (*.fit*)")
                
                if not filenames.isEmpty():
                    calibrations['flat'] = (str(filenames[0]))
                
                # Select master BPM: by the momment, it is read from config file.
                if self.config_opts['bpm']['mode']!='none':
                    calibrations['bpm'] = self.config_opts['bpm']['bpm_file']
      
            elif msgBox.clickedButton()== button_autosearch:
                # We will look for default cailbration for each file
                # Default calibration files must be a MASTER_DARK/MASTER_FLAT
                force_apply = False
                auto_search = True
        else:
            # Cancel button pressed
            return

        # Select BPM action
        msgBox = QMessageBox()
        msgBox.setText("        BPM Action to do")
        msgBox.setInformativeText("Do you want to <Set NaNs (=0)> or <Fix> BPM ?")
        button_nans = msgBox.addButton("NaNs", QMessageBox.ActionRole)
        button_fix = msgBox.addButton("Fix Bab Pixels", QMessageBox.ActionRole)
        button_cancel = msgBox.addButton("Cancel", QMessageBox.ActionRole)
        msgBox.setDefaultButton(button_nans)
        msgBox.exec_()
        if msgBox.clickedButton()==button_nans:
            bpm_mode = 'grab' # set NaNs
        elif msgBox.clickedButton()==button_fix:
            bpm_mode = 'fix'
        else:
            bpm_mode = 'none'
            
        # Now, start dark subtraction and Flat-Fielding...
        if len(self.m_popup_l_sel)>0:
            for filename in self.m_popup_l_sel:
                try:
                    mDark = None
                    mFlat = None
                    mBPM = None
                    
                    #bpm_mode = self.config_opts['bpm']['mode']
                    #if self.comboBox_bpm_action.currentText()=="None":
                    #    bpm_mode = 'none'
                    #elif self.comboBox_bpm_action.currentText()=="Set NaNs":
                    #    bpm_mode = 'grab'
                    #elif self.comboBox_bpm_action.currentText()=="Fix":
                    #    bpm_mode = 'fix'
                    
                    if auto_search:
                        # Look for (last received) calibration files
                        mDark, mFlat, mBPM = self.getCalibFor([filename])
                    else:
                        if calibrations['dark']!='':
                            mDark = calibrations['dark']
                        if calibrations['flat']!='':
                            mFlat = calibrations['flat']
                        if bpm_mode!='none':
                            mBPM = self.config_opts['bpm']['bpm_file']

                    # Show message with calibrations found
                    # Ask for confirmation
                    msg = "DARK= %s \n FLAT= %s \n BPM= %s"%(mDark, mFlat, mBPM)
                    resp = QMessageBox.information(self, "Info", 
                        QString("Calibrations found are: \n %1").arg(str(msg)),
                                     QMessageBox.Ok, QMessageBox.Cancel)
                    if resp==QMessageBox.Cancel:
                        return
                    
                    log.debug("Source file: %s"%filename)
                    log.debug("Calibrations to use - DARK: %s   FLAT: %s  BPM: %s"%(mDark, mFlat, mBPM))
                    self.logConsole.info("Calibrations to use:")
                    self.logConsole.info("+ Dark :%s"%mDark)
                    self.logConsole.info("+ Flat :%s"%mFlat)
                    self.logConsole.info("+ BPM (mode=%s) :%s"%(bpm_mode,mBPM))
                    
                    # Both master_dark and master_flat and BPM are optional
                    if mDark or mFlat or mBPM:
                        # Put into the queue the task to be done
                        func_to_run = reduce.ApplyDarkFlat([filename], 
                                                         mDark, mFlat, mBPM, 
                                                         self.m_outputdir,
                                                         bpm_action=bpm_mode, # fix is a heavy process for QL
                                                         force_apply=force_apply)
                        params = ()
                        log.debug("Inserting in queue the task ....")
                        self._task_queue.put([(func_to_run.apply, params)])
                        
                    else:
                        self.logConsole.error("[ApplyDarkFlat] Cannot find the appropriate master calibrations for file %s"%filename)
                        #QMessageBox.critical(self, 
                        #                     "Error", 
                        #                     "Error, cannot find the master calibration files")
                except Exception, e:
                    QMessageBox.critical(self, "Error", "Error while processing file.  %s"%str(e))
                    #self.m_processing = False
                    #QApplication.restoreOverrideCursor()
                    raise e 

    def applyBPM_slot(self):
        """
        Apply to the selected file the BPM defined in the config file and :
        
         1) 'grab': show in Red the masked pixels on DS9 (for this, user 
             should select red from menu Edit->Preferences->Blank/Inf/NaN color).
         2) 'fix': Bad Pixels are replaced with a bi-linear interpolation 
             from nearby pixels.
        
        Note that, the original selected file is not modified, all BPM is applied
        to a new output file.

        This routine can be useful to view on QL how the BPM affect our data.
        """

        if len(self.m_popup_l_sel)>0:
            
            msgBox = QMessageBox()
            msgBox.setText("Method selection:")
            msgBox.setInformativeText("Do you want to <Grab (set to NaN)>  or <Fix> the Bad Pixels?")
            button_grab = msgBox.addButton("Grab (set to NaN)", QMessageBox.ActionRole)
            button_fix = msgBox.addButton("Fix (bi-linear interpol)", QMessageBox.ActionRole)
            msgBox.setDefaultButton(button_grab)
            msgBox.exec_()

            if msgBox.clickedButton()== button_fix: bpm_mode = 'fix'
            elif msgBox.clickedButton()== button_grab: bpm_mode = 'grab'
            else: return
            
            #init_outdir = self.m_outputdir + "/" + \
            #        os.path.basename(self.m_popup_l_sel[0]).replace(".fits","_BPM.fits")
            #outfileName = QFileDialog.getSaveFileName(self,
            #                                          "Choose a filename so save under",
            #                                          init_outdir, 
            #
            #                                      "fits (*.fits)", )
            
            #if not outfileName.isEmpty():
            # Because in principle it is a quick proceduce, we do not use
            # the processing queue. We run it on the current thread.
            QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
            try:
                master_bpm =  self.config_opts['bpm']['bpm_file']
                if os.path.isfile(master_bpm):
                    task = reduce.ApplyDarkFlat(self.m_popup_l_sel, 
                                            None, None, 
                                            master_bpm, 
                                            self.m_outputdir,
                                            bpm_action=bpm_mode)
                    result = task.apply()
                    #result = reduce.calBPM.applyBPM(
                    #                self.m_popup_l_sel[0],
                    #                master_bpm,
                    #                str(outfileName),
                    #                False)
                    
                    msg = "New file created: %s"%result[0]
                    self.logConsole.info(msg)
                    
                    if self.getDisplayMode()>=2:
                      display.showFrame(result[0])
                else:
                    msg = "Master BPM '%s' does not exist."%master_bpm
                    QMessageBox.critical(self, "Error", msg)
            except Exception,e:
                msg = "Error applying BPM: '%s'"%(str(e))
                self.logConsole.info(msg)
                QMessageBox.critical(self, "Error", msg)
                raise e
            finally:
                QApplication.restoreOverrideCursor()

    def focus_eval(self):
        """
        Ask what method is prefered, iraf.starfocus or SExtractor based.
        """

        msgBox = QMessageBox()
        msgBox.setText("Method selection:")
        msgBox.setInformativeText("Do you want to use <iraf.starfocus> or <SExtractor>?")
        button_iraf = msgBox.addButton("IRAF", QMessageBox.ActionRole)
        button_sex = msgBox.addButton("SExtractor", QMessageBox.ActionRole)
        msgBox.setDefaultButton(button_iraf)
        msgBox.exec_()

        if msgBox.clickedButton()== button_iraf: self.focus_eval_iraf()
        elif msgBox.clickedButton()==button_sex: self.focus_eval_sex()
        else: pass

        
    def focus_eval_iraf(self):
        """
        Run the focus evaluation procedure of a set of files of focus serie.
        It is run **interactively**.
        """
        
        if len(self.m_popup_l_sel)>3:
            with fits.open(self.m_popup_l_sel[0]) as myfits:
                if len(myfits)>1:
                    # Convert MEF-files to join-files (iraf.starfocus only support that kind of files)
                    try:
                        mef = misc.mef.MEF(self.m_popup_l_sel)
                        my_files = mef.doJoin(".join.fits", output_dir=self.m_outputdir)[1]
                    except Exception,e:
                        log.debug("Cannot convert MEF to Single file %s. Maybe it's not a MEF file", str(e))
                        QMessageBox.critical(self, "Error", "Cannot convert MEF to Single file : %s \n Maybe it's not a MEF file"%(file))
                    else:
                        line = "Files generated: \n%s\n"%my_files
                        self.logConsole.info(QString(str(line)))
                else:
                    # We have non-MEF files
                    my_files = list(self.m_popup_l_sel)
            
            text_file = self.focus_tmp_file  # ~/iraf/focus_seq.txt
            iraf_logfile = self.m_tempdir + "/starfocus.log"
            # copy the list and swap first and center position
            my_files[0], my_files[len(my_files)/2] = my_files[len(my_files)/2], my_files[0]    
            # copy list to file; if file exists, overwrite
            self.genFileList(my_files, text_file)
            try:
                ## New approach
                #os.system("/home/panic/DEVELOP/papi/commissioning/runStarfocus.py -s %s &"%text_file)
                #focus.runFocusEvaluation(text_file, "", iraf_logfile)
                # if not started, launch DS9
                display.startDisplay()
                # Launch IRAF; note that next call is asynchronous, that is,
                # it returns inmediately after launch IRAF.
                self.iraf_console_slot(True)
            except Exception,e:
                os.unlink(text_file)
                log.error("Error, cannot run iraf.obsutil.starfocus(): %s"%str(e))
                QMessageBox.critical(self, "Error", str(e))
            else:
                # text_file is removed by papi_ql_user.cl script
                pass
        else:
            QMessageBox.critical(self, "Error","Error, not enough number (>3) of frames selected")

    def focus_eval_sex(self):
        """
        Run the focus evaluation procedure of a set of files of focus serie.
        It is run **non interactively**.
        
        NOTE: it only works for MEF files !!!

        """
        
        if len(self.m_popup_l_sel)>3:
            with fits.open(self.m_popup_l_sel[0]) as myfits:
                if len(myfits)!=5:
                    msg = "This routine only works for MEF files. Run FITS->Single2MEF."
                    self.logConsole.error(msg)
                    QMessageBox.critical(self, "Error", msg)
                    return
                
            init_outdir = self.m_tempdir + "/focus_eval.pdf"
            outfileName = QFileDialog.getSaveFileName(self,
                                                      "Choose a filename so save under",
                                                      init_outdir, 
                                                      "pdf (*.pdf)", )
            if not outfileName.isEmpty():
                # Ask for detector to use 
                msgBox = QMessageBox()
                msgBox.setText("Detector selection:")
                msgBox.setInformativeText("Do you want to use <ALL> detectors or <SELECT> one ?")
                button_all = msgBox.addButton("ALL", QMessageBox.ActionRole)
                button_sg1 = msgBox.addButton("SG1/Q1", QMessageBox.ActionRole)
                button_sg2 = msgBox.addButton("SG2/Q2", QMessageBox.ActionRole)
                button_sg3 = msgBox.addButton("SG3/Q3", QMessageBox.ActionRole)
                button_sg4 = msgBox.addButton("SG4/Q4", QMessageBox.ActionRole)
                msgBox.setDefaultButton(button_sg2)
                msgBox.exec_()
        
                if msgBox.clickedButton()== button_all: detector = 'all'
                elif msgBox.clickedButton()==button_sg1: detector = 'Q1'
                elif msgBox.clickedButton()==button_sg2: detector = 'Q2'
                elif msgBox.clickedButton()==button_sg3: detector = 'Q3'
                elif msgBox.clickedButton()==button_sg4: detector = 'Q4'
                else: detector = 'all'
                    
                #
                show = True
                try:
                    self.genFileList(self.m_popup_l_sel, "/tmp/focus.list")
                    QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
                    pix_scale = self.config_opts['general']['pix_scale']
                    satur_level =  self.config_opts['skysub']['satur_level']
                    #detector = self.config_opts['general']['detector']
                    self._task = reduce.eval_focus_serie.FocusSerie("/tmp/focus.list", 
                                                               str(outfileName),
                                                               pix_scale, 
                                                               satur_level,
                                                               show=True,
                                                               window=detector)

                    best_focus, outfile = self._task.eval_serie()
                    self.logConsole.info(str(QString("Best Focus (mm)= %1 --> File = %2")
                                            .arg(best_focus)
                                            .arg(os.path.basename(outfile))))
                    QApplication.restoreOverrideCursor()
                    
                    #thread = reduce.ExecTaskThread(self._task.eval_serie, 
                    #                                  self._task_info_list)
                    #thread.start()
                except Exception,e:
                    # Restore the cursor
                    QApplication.restoreOverrideCursor()
                    msg = "Cannot evaluate focus series: %s"%(str(e))
                    self.logConsole.error(msg)
                    QMessageBox.critical(self, "Error", msg)
                    raise e
        else:
            QMessageBox.critical(self, "Error","Error, not enough number of frames selected")

    def show_stats_slot(self):
        """Show image statistics in the log console of the files selected"""
        
        length_fn = len(self.m_popup_l_sel[0])
        msg = "FILE" + (length_fn+25)*" " + "MEAN     MODE       STDDEV      MIN        MAX"
        self.logConsole.info(msg)
        for file in self.m_popup_l_sel:
            values = (iraf.mscstat (images=file,
            fields="image,mean,mode,stddev,min,max",format='no',Stdout=1))
            for line in values:
                #line=os.path.basename(values[0])
                self.logConsole.info(str(line))
        
    def background_estimation_slot(self):
        """ 
        Give an background estimation of the current selected image. 
        """
        
        pix_scale = self.config_opts['general']['pix_scale']
        cq = reduce.checkQuality.CheckQuality(self.m_popup_l_sel[0],
                                              pixsize=pix_scale)
        try:     
            img = cq.estimateBackground(self.m_outputdir+"/bckg.fits")
            
            values = (iraf.mscstat(images=img,
            fields="image,mean,mode,stddev,min,max",format='yes',Stdout=1))
            #file,mean,mode,stddev,min,max=values[0].split()
            self.logConsole.info(str(QString("Background estimation :")))
            for line in values:
                self.logConsole.info(str(line))
                #self.logConsole.info(QString("Background estimation MEAN= %1   MODE=%2    STDDEV=%3    MIN=%4         MAX=%5").arg(mean).arg(mode).arg(stddev).arg(min).arg(max))
            if self.getDisplayMode()>=2:
                display.showFrame(img)
        except Exception,e:
            self.logConsole.error("ERROR: something wrong while computing background")
            raise e
          
    def fwhm_estimation_slot (self):
        """ Give an FWHM estimation of the current selected image """
        
        # Ask for detector to be used 
        msgBox = QMessageBox()
        msgBox.setText("Detector selection:")
        msgBox.setInformativeText("Do you want to use <ALL> detectors or <SELECT> one ?")
        button_all = msgBox.addButton("ALL", QMessageBox.ActionRole)
        button_sg1 = msgBox.addButton("SG1/Q1", QMessageBox.ActionRole)
        button_sg2 = msgBox.addButton("SG2/Q2", QMessageBox.ActionRole)
        button_sg3 = msgBox.addButton("SG3/Q3", QMessageBox.ActionRole)
        button_sg4 = msgBox.addButton("SG4/Q4", QMessageBox.ActionRole)
        msgBox.setDefaultButton(button_sg2)
        msgBox.exec_()
        
        if msgBox.clickedButton()== button_all: detector = 'all'
        elif msgBox.clickedButton()==button_sg1: detector = 'Q1'
        elif msgBox.clickedButton()==button_sg2: detector = 'Q2'
        elif msgBox.clickedButton()==button_sg3: detector = 'Q3'
        elif msgBox.clickedButton()==button_sg4: detector = 'Q4'
        else: detector = 'all'
        #
        
        for ifile in self.m_popup_l_sel:
            file_info = self.inputsDB.GetFileInfo(ifile)
            if file_info!=None and (file_info[2]!='SCIENCE' and file_info[2]!='FOCUS'):
                QMessageBox.warning(self, "Warning", "Selected file is not SCIENCE type.")
            elif file_info==None:
                if self.outputsDB.GetFileInfo(ifile)[2]!='SCIENCE':
                    QMessageBox.warning(self, "Warning", "Selected file is not SCIENCE type.")
            
            pix_scale = self.config_opts['general']['pix_scale']
            satur_level = self.config_opts['skysub']['satur_level']
            iso_min_size = self.config_opts['skysub']['mask_minarea']
            cq = reduce.checkQuality.CheckQuality(ifile,
                                                  sat_level=satur_level, 
                                                  isomin=iso_min_size,
                                                  ellipmax=0.9, # basically, no limit !
                                                  pixsize=pix_scale,
                                                  window=detector)
            try:
                fwhm,std,k,k = cq.estimateFWHM()
                if fwhm>0:
                    self.logConsole.info(str(QString("%1  FWHM = %2 (pixels) std= %3")
                                             .arg(os.path.basename(ifile))
                                             .arg(fwhm)
                                             .arg(std)))
                else:
                    self.logConsole.error("ERROR: Cannot estimate FWHM of selected image.")           
            except Exception,e:
                self.logConsole.error("ERROR: Cannot estimate FWHM of selected image: %s"%str(e))
                #raise e
    
    def do_quick_reduction_slot(self):
        """ 
        Run a quick-reduction (subtract sky, shift and align) from the set of 
        frames, selected by user or automatic search based of nearest 
        (ra,dec, mjd) ones.
        """
        
        file_n = 0 # really, not used for the moment
        file_list = []
        
        ra_dec_near_offset = self.lineEdit_ra_dec_near_offset.text().toInt()[0]/3600.0
        time_near_offset = self.lineEdit_time_near_offset.text().toInt()[0]/86400.0
        print "RA_DEC_OFFSET", ra_dec_near_offset
        print "TIME_NEAR_OFFSET", time_near_offset
        
        # CASE 1: Automatic search for nearest frames (ra, dec, mjd) 
        if len(self.m_popup_l_sel)==1:
            QMessageBox.information(self, "Info", "Only one file was selected, automatic file grouping will be done.")
            fits = datahandler.ClFits(self.m_listView_item_selected)
            if fits.getType()=='SCIENCE':
                near_list = self.inputsDB.GetFiles('ANY', 'SCIENCE', -1, 
                                                   fits.filter, 
                                                   fits.mjd, fits.ra, fits.dec, 
                                                   ra_dec_near_offset*2, 
                                                   time_near_offset, runId=0)
                #print "NEAR_LIST=", near_list
                # For the moment, the minimun number of nearest is >0
                if len(near_list)==0:
                    QMessageBox.information(self, "Info", "Not enough science frames found")  
                    return
            else:
                QMessageBox.information(self, "Info", "Selected frame is not a science frame") 
                return
            
            # Create nearest file list (ar,dec,mjd) from current selected science file
            view_list = ""
            i = 1
            for file in near_list:
                if file==self.m_listView_item_selected: 
                    file_n=i
                else: 
                    i=i+1
                if (datahandler.ClFits(file).getType()!='SCIENCE'):
                    QMessageBox.critical(self, "Error", 
                                         QString("File %1 is not a science frame")
                                         .arg(file))
                else:
                    file_list.append(file)
                    view_list +=  file + "\n"
            
            resp = QMessageBox.information(self, "Info", 
                                         QString("Selected near frames are:\n %1")
                                         .arg(view_list),
                                         QMessageBox.Ok, QMessageBox.Cancel)
            if resp==QMessageBox.Cancel:
                return
                    
        # CASE 2: Stack frames selected by user in the list_view
        elif len(self.m_popup_l_sel)>1 and len(self.m_popup_l_sel)<5:
            QMessageBox.information(self,"Info","Error, not enough frames selected to reduce (>4) !")
            return    
        elif len(self.m_popup_l_sel)>=5:
            # Create file list from current selected science files
            for file in self.m_popup_l_sel :
                if (datahandler.ClFits(file).getType()!='SCIENCE'):
                    QMessageBox.critical(self, "Error", QString("File %1 is not a science frame")
                                         .arg(file))
                    return
                else: file_list.append(file)
            file_n=-1 # actually, not used
            
        # Select the name of the output result file
        outfileName = QFileDialog.getSaveFileName(self,
                                                  "Choose a filename to save under",
                                                  self.m_outputdir+"/red_result.fits", 
                                                  "*.fits")
        if outfileName.isEmpty(): return # nothig to do !
           
        #Create working thread that compute sky-frame
        if len(file_list)>1:
            try:
                #because choosen files might not be a complete sequence we
                #set group_by to 'filter'
                self.processFiles(file_list, group_by='filter', 
                                  outfilename=str(outfileName))
            except Exception, e:
                QMessageBox.critical(self, "Error", "Error while building stack")
                raise e
        
        
    def do_raw_astrometry(self):
        """
        Compute an astrometric solution for the selected file in the 
        main list view panel.
        Basically, it is SCAMP is called and compared to 2MASS catalog 
        (or GSC2.3, USNO-B) to compute astrometric calibration. 
        A new WCS is added to the header.
        
        Returns
        -------
        None, but produce an astrometrically calibrated image.
          
        Todo
        ---- 
        We should check the image is pre-reduced or at least, with 
        enough objects.
        """
        
        if len(self.m_popup_l_sel)==1:
            fits = datahandler.ClFits(self.m_listView_item_selected)
            if fits.isMEF():
                  QMessageBox.information(self,"Info", QString("Sorry, but it only "
                    "works for single detertor images (2kx2k). "
                    "Run FITS->Split MEF."))
                  return
            if fits.isPANICFullFrame():
                  QMessageBox.information(self,"Info", QString("Sorry, but it only "
                    "works for single detertor images (2kx2k). "
                    "Run FITS->Split Single."))
                  return
            if True: #fits.getType()=='SCIENCE': 
                # Run astrometry parameters
                out_file = self.m_outputdir+"/"+os.path.basename(self.m_listView_item_selected.replace(".fits",".wcs.fits"))
                # Catalog
                if self.comboBox_AstromCatalog.currentText().contains("2MASS"): 
                    catalog = "2MASS"
                elif self.comboBox_AstromCatalog.currentText().contains("USNO-B1"): 
                    catalog = "USNO-B1"
                elif self.comboBox_AstromCatalog.currentText().contains("GSC-2.2"): 
                    catalog = "GSC-2.2"
                else: 
                    catalog = "2MASS"
                
                # Change to working directory
                os.chdir(self.m_tempdir)
                # Change cursor
                QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
                # Create working thread that compute sky-frame
                try:
                    self.logConsole.info("Starting Astrometric calibration (%s)"%catalog)
                    if self.comboBox_AstromEngine.currentText()=="SCAMP":
                        self.logConsole.info(" +Engine: SCAMP")
                        thread = reduce.ExecTaskThread(reduce.astrowarp.doAstrometry, 
                                                       self._task_info_list,
                                                       #args of the doAstrometry function
                                                       self.m_listView_item_selected, 
                                                       out_file, catalog, 
                                                       self.config_opts)
                    else: # Astrometry.net
                        self.logConsole.info(" +Engine: Astrometry.net")
                        thread = reduce.ExecTaskThread(reduce.solveAstrometry.solveField,
                                                    self._task_info_list,
                                                    # args of the solveAstrometry
                                                    self.m_listView_item_selected, 
                                                    self.m_outputdir,
                                                    self.config_opts['general']['pix_scale'])
                    thread.start()
                except:
                    QMessageBox.critical(self, "Error", "Error while computing astrometry.")
                    raise
            else:
                QMessageBox.information(self,"Info", QString("Sorry, but you need a science frame."))
        
    def do_raw_photometry(self):
        """
        Compute an rough photometric solution for the selected file in the main 
        list view panel.
        Basically, It is compared to 2MASS catalog (or other) to compute photometric
        calibration, and obtain a zero-point estimation.

        Returns
        -------
        None, but produce an photometric fitting curve and a Zero point
        estimation.
        
        Todo
        ----  
        Check the image is pre-reduced or at least, with enough objects 
        """
        
        if len(self.m_popup_l_sel)==1:
            fits = datahandler.ClFits(self.m_listView_item_selected)
            if fits.getType()=='SCIENCE': 
                ## Run photometry parameters
                out_file = self.m_outputdir+"/" + \
                    os.path.basename(self.m_listView_item_selected.replace(".fits","_photo.pdf"))
                # Catalog
                if self.comboBox_AstromCatalog.currentText().contains("2MASS"): 
                    catalog = "2MASS"
                elif self.comboBox_AstromCatalog.currentText().contains("USNO-B1"): 
                    catalog = "USNO-B1"
                elif self.comboBox_AstromCatalog.currentText().contains("GSC-2.2"): 
                    catalog = "GSC-2.2"
                else: 
                    catalog = "2MASS"
                
                #Change to working directory
                os.chdir(self.m_tempdir)
                #Change cursor
                QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
                #Create working thread that compute sky-frame
                try:
                    snr = 10.0
                    zero_point = 0.0
                    self.logConsole.info("Starting Photometric calibration (%s)"%catalog)
                    thread = reduce.ExecTaskThread(photo.photometry.doPhotometry,
                                                   self._task_info_list,
                                                   #args of the doPhotometry function
                                                   self.m_listView_item_selected,
                                                   self.config_opts['general']['pix_scale'],
                                                   catalog,
                                                   out_file,
                                                   snr,
                                                   zero_point)
                    thread.start()
                except:
                    QMessageBox.critical(self, "Error", "Error while computing photometry.")
                    raise
            else:
                QMessageBox.information(self,"Info", QString("Sorry, but you need a science frame."))
    
      
    def createCalibs_slot(self):
        
        """
        Given the current dataset files, compute all the master calibration 
        files found in the DB dataset.
        """
        
        fileList = self.inputsDB.GetFilesT("ANY")
        
        if len(fileList) > 1:
            # Ask for confirmation
            resp = QMessageBox.information(self, "Info", 
                                         QString("All calibrations will be built into OutDir using files of InputDir"),
                                         QMessageBox.Ok, QMessageBox.Cancel)
            if resp==QMessageBox.Cancel:
                return
            
            # Create the RS            
            try:
                calib_db_files = None
                params = [(fileList, self.m_outputdir, 
                        None,
                        "dither", None, 
                        None, None, "quick",
                        self.group_by, True,
                        self.config_opts,
                        calib_db_files, None)]
            
                log.debug("Let's create Calibrations ...")
                func_to_run = RS.ReductionSet(*(params[0])).buildCalibrations
                log.debug("ReductionSet created !")
                self._task_queue.put([(func_to_run,())])
                log.debug("New task queued")

            except Exception,e:
                QMessageBox.critical(self, "Error", "Error while building master calibrations files")
                print e
                raise Exception("Error creating calibrations")
        else:
            QMessageBox.information(self,"Info", QString("No files found"))


    def helpIndex(self):
        """Show help on the default browser"""
        
        import webbrowser
        new = 2 # open in a new tab, if possible

        # open a public URL, in this case, the webbrowser docs
        url = "http://www.iaa.es/~jmiguel/PANIC/PAPI/html/index.html"
        try:
            webbrowser.open(url,new=new)
        except Exception,e:
            log.erro("Cannot open URL: %s"%str(e))

        # open an HTML file on my own (Windows) computer
        #url = "file://X:/MiscDev/language_links.html"
        #webbrowser.open(url,new=new)
    
    def testSlot(self):
        """
        Not implemeted
        """
        
        msgBox = QMessageBox()
        msgBox.setText("Starting to process files.")
        msgBox.setInformativeText("Which files do you wish to process ?")
        button1 = msgBox.addButton("Process all", QMessageBox.ActionRole)
        button2 = msgBox.addButton("Process new", QMessageBox.ActionRole)
        button3 = msgBox.addButton("Cancel", QMessageBox.ActionRole)
        msgBox.setDefaultButton(button2)
        
        msgBox.exec_()
        
        if msgBox.clickedButton()== button1:
            print "Button1"
        elif msgBox.clickedButton()== button2:
            print "Button2"
        elif msgBox.clickedButton()== button3:
            print "Button3"
        else:
            print "[Default] button3"
        
        pass
    
    def pushB_start_stop_slot(self):
        """
        Method called when "Start/Stop-processing button is clicked
        """
        
        # Change the button color and label
        if self.proc_started == True:
            self.proc_started = False
            self.pushButton_start_proc.setText("Start Processing")
            #Set Green background
            self.pushButton_start_proc.setAutoFillBackground(True)
            self.pushButton_start_proc.setStyleSheet("background-color: rgb(0, 128, 0); color: rgb(0, 0, 0)")
            self.stopProcessing()
        else:    
            # Ask if we with to process only the next files coming or also the
            # current ones in the source_dir and then the new ones
            msgBox = QMessageBox()
            msgBox.setText("Starting to process files.")
            msgBox.setInformativeText("Which files do you wish to process ?")
            button_all = msgBox.addButton("Process all", QMessageBox.ActionRole)
            button_new = msgBox.addButton("Process new", QMessageBox.ActionRole)
            button_cancel = msgBox.addButton("Cancel", QMessageBox.ActionRole)
            msgBox.setDefaultButton(button_new)
            
            msgBox.exec_()
            
            if msgBox.clickedButton()== button_all or msgBox.clickedButton()==button_new:
                self.proc_started = True
                self.pushButton_start_proc.setText("Stop Processing")
                #Set Red background
                self.pushButton_start_proc.setAutoFillBackground(True)
                self.pushButton_start_proc.setStyleSheet("background-color: rgb(255, 0, 0); color: rgb(0, 0, 0)")
                if msgBox.clickedButton()== button_all:
                    self.processFiles()
                
 
    def stopProcessing(self):
        """Stop the processing throwing away all the tasks in the queue"""
        while not self._task_queue.empty():
            self._task_queue.get_nowait()
            log.debug("Queue empted")
        
        if self._process!=None and self._process.is_alive():
            self._process.terminate()
        #else:
        #    print "da igual, lo paro"
        #    self._process.terminate()
        
        self.m_processing = False
        QApplication.restoreOverrideCursor()
        log.info("Processing stopped !")
        self.logConsole.debug("Processing stopped !")
            
        
    def worker1(self, input, output):
        #print "arg_1=",input.get()
        for func, args in iter(input.get, 'STOP'):
            #result = calculate(func, args)
            print "ARGS=", args
            output.put(RS.ReductionSet(args).func())
            
    def processFiles(self, files=None, group_by=None, outfilename=None,
                     interactive=False):
        """
        Process the files provided; if any files were given, all the files 
        in the current Source List View (but not the output files) will be
        processed. The processing task will be inserted into the _task_queue
        to be processed as soon as the current (if any) processing has finished
        (m_processing=False).  
        
        Parameters
        ----------
        files: list
            files list to be processed; if None, all the files in the
            current List View will be used, doing a previous data grouping.
        
        group_by: str
            if 'ot' the data grouping will be done using OT keywords,
            of if equal to 'filter', data grouping will be done using
            AR, Dec and Filter keywords. 
        
        outfile: str
            Filename of the final output file (aligned & stacked) produced.
            If no outfile name is given, the result of each sequence reduced
            will be saved with a filename as: 'PANIC.[DATE-OBS].fits',
            where DATE-OBS is the keyword value of the first file in the sequence.
        
        Notes
        -----
        When the RS object is created, we provide the outputDB as the 
        external DB for the RS; this way fomer master calibration files can be
        used for the current data reduction (e.g., TwFlats who need a master dark)
        """
        
        if files==None:
            # All the files in ListView, even the generated output 
            # files (but they will not be processed)
            files = self.inputsDB.GetFiles() 

        if len(files)==0:
            return
        
        log.debug("Starting to process files...")
        self.logConsole.info("++ Starting to process next [%d] files :"%len(files))
        for file in files:
            self.logConsole.info("     - " + file)
            #self.logConsole.info(QString("    - %1").arg(file))
        self.logConsole.info("... processing sequence ...")
            
        # Create working thread that process the files
        try:
    
            ###self._task = RS.ReductionSet(files, self.m_outputdir, out_file=outfilename,
            ###                                obs_mode="dither", dark=None, 
            ###                                flat=None, bpm=None, red_mode="quick",
            ###                                group_by=self.group_by, check_data=True,
            ###                                config_dict=self.config_opts,
            ###                                external_db_files=self.outputsDB.GetFiles())
            # provide the outputDB files as the external calibration files for the RS 
            
            # New approach using multiprocessing
            #rs_filelist, out_dir=None, out_file=None, obs_mode="dither", 
            #     dark=None, flat=None, bpm=None, red_mode="quick", 
            #     group_by="ot", check_data=True, config_dict=None, 
            #     external_db_files=None, temp_dir = None,
            
            if group_by==None: l_group_by = self.group_by
            else: l_group_by = group_by
            
            # Here, it is decided if last calibration files will be used
            # GUI preferences have higher priority than values selected in
            # config file. 
            if self.checkBox_pre_subDark_FF.checkState()==Qt.Checked:
                calib_db_files = self.outputsDB.GetFiles()
                self.config_opts['general']['apply_dark_flat'] = 1
                log.debug("ext-calibretion DB loaded")
            else:
                calib_db_files = None
                self.config_opts['general']['apply_dark_flat'] = 0

            #
            # Load config values from Setup Tab
            #
            if self.checkBox_pre_appNLC.isChecked():
                self.config_opts['nonlinearity']['apply'] = True
            else:
                self.config_opts['nonlinearity']['apply'] = False
            
            if self.checkBox_pre_appBPM.isChecked():
                self.config_opts['bpm']['mode'] = 'grab' # mark Bad Pixels to NaN
            else:
                self.config_opts['bpm']['mode'] = 'none'

            #
            if self.comboBox_AstromEngine.currentText()=="SCAMP":
                self.config_opts['astrometry']['engine'] = "SCAMP"
            else:
                self.config_opts['astrometry']['engine'] = "AstrometryNet"
            #
            if self.comboBox_pre_skyWindow.currentText()=="2-frames": nhw = 1
            elif self.comboBox_pre_skyWindow.currentText()=="4-frames": nhw = 2
            elif self.comboBox_pre_skyWindow.currentText()=="6-frames": nhw = 3
            elif self.comboBox_pre_skyWindow.currentText()=="8-frames": nhw = 4
            else: nhw = 2
            self.config_opts['skysub']['hwidth'] = nhw
            
            # Select detector (map SGi to Qi)
            if self.comboBox_detector.currentText()=="All": detector = 'all'
            elif self.comboBox_detector.currentText()=="SG1": detector = 'Q1'
            elif self.comboBox_detector.currentText()=="SG2": detector = 'Q2'
            elif self.comboBox_detector.currentText()=="SG3": detector = 'Q3'
            elif self.comboBox_detector.currentText()=="SG4": detector = 'Q4'
            elif self.comboBox_detector.currentText()=="SG123": detector = 'Q123'
            else: detector = 'all'
            self.config_opts['general']['detector'] = detector

            

            #
            params = [(files, self.m_outputdir, outfilename,
                      "dither", None, 
                      None, None, "quick",
                      l_group_by, True,
                      self.config_opts,
                      calib_db_files, None)]
            
            log.debug("Let's create a ReductionSet ...")

            # Show the calibration to be used
            if self.config_opts['general']['apply_dark_flat']!=0:
                # Look for calibs
                md, mf, mb = self.getCalibFor(files)
                
                # Show message with calibrations found
                # and ask for confirmation.
                if interactive:
                    msg = "DARK= %s \n FLAT= %s \n BPM= %s"%(md, mf, mb)
                    resp = QMessageBox.information(self, "Info", 
                            QString("Calibrations found are: \n %1").arg(str(msg)),
                                        QMessageBox.Ok, QMessageBox.Cancel)
                    if resp==QMessageBox.Cancel:
                        return
                
                self.logConsole.info("Calibrations to use:")
                self.logConsole.info("+ Dark=%s"%md)
                self.logConsole.info("+ Flat=%s"%mf)
                self.logConsole.info("+ BPM =%s"%mb)
            else:
                self.logConsole.info("No calibrations used.")

            func_to_run = RS.ReductionSet(*(params[0])).reduceSet
            log.debug("ReductionSet created !")
            self._task_queue.put([(func_to_run,())])
            log.debug("New task queued") 
            #self._task_queue.put(params) #default function supposed!
            
            ##Process(target=self.worker, 
            ##        args=(self._task_queue, self._done_queue)).start()
                    
        except Exception,e:
            QMessageBox.critical(self, "Error", "Error while processing Obs. Sequence: \n%s"%str(e))
            raise e # Para que seguir elevando la excepcion ?
        
    def taskRunner(self):
        """
        Procedure that continuisly in checking the queue of pending tasks to be 
        done. The results are obtained later at checkDoneQueue().

        The tasks are processed sequentially, that is, a new proccesing start only
        when the previous one has finished. 
        """
 
        # Update the number of tasks (not necessarialy equals to number of
        # sequences) in to queue to be processed !
        # Current processing will also be taken into account.
        self.lineEdit_queue_size.setText(str(self._task_queue.qsize() + 
                                            int(self.m_processing)))
        
        if not self._task_queue.empty() and not self.m_processing:
            log.debug("Something new in the TaskQueue !")
            try:
                self.logConsole.debug("Starting to process queued task")
                #Change cursor
                QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
                #'mutex' variable 
                self.m_processing = True
                self._process = Process(target=self.worker, 
                    args=(self._task_queue, self._done_queue))
                self._process.start()
            except Exception,e:
                #NOTE: I think this point will never be reached !!!
                log.error("Error in task Runner: %s"%(str(e)))
                self.logConsole.debug("Error in taskRunner")
                self.m_processing = False
                QApplication.restoreOverrideCursor()
            finally:
                log.debug("End of taskRunner")
                 
    def worker_original(self, input, output):
        """
        NOT USED - Callback function used by Process task
        """
        
        args = input.get()
        print "ARGS=",args
        
        try:
            output.put(RS.ReductionSet(*(args[0])).reduceSet())
            log.info("[workerd] task done !")
        except Exception,e:
            log.error("[worker] Error while processing task")
            output.put(None) # the DoneQueue will detect it
        finally:
            self.m_processing = False
            log.debug("Worker finished its task !")
    
    def worker(self, input, output):
        """
        Callback function used by Process task.
        It is run in a separate/parallel process, similar to a thread.
        
        Parameters
        ----------
        
        input : Queue
            Input queue of tasks to be done
        output : Queue
            Output queue of tasks already done
            
        """
        
        # Get the next task from the queue
        func, args = input.get()[0]
        #print "FUNC=",func
        #print "ARGS=",args
        
        log.debug("[worker] Worker start to work ...")
        try:
            output.put(func(*args))
            log.info("[worker] task done !")
        except Exception,e:
            log.error("[worker] Error while processing task: %s"%str(e))
            # Because excetions cannot be catched in taskRunner due to a 
            # multiprocessing.Process exception is inserted in output queue 
            # and then recognized
            output.put(e) # the DoneQueue timer will detect it
        finally:
            self.m_processing = False
            log.debug("Worker finished its task !")
    
    # Menus stuff functions 
    def editCopy(self):
        print "panicQL.editCopy(): Not implemented yet"

    def fileSave(self):
        print "panicQL.fileSave(): Not implemented yet"

    def fileSaveAs(self):
        print "panicQL.fileSaveAs(): Not implemented yet"

    def filePrint(self):
        print "panicQL.filePrint(): Not implemented yet"

    def editUndo(self):
        print "panicQL.editUndo(): Not implemented yet"

    def editRedo(self):
        print "panicQL.editRedo(): Not implemented yet"

    def editCut(self):
        print "panicQL.editCut(): Not implemented yet"

    def editPaste(self):
        print "panicQL.editPaste(): Not implemented yet"

    def editFind(self):
        print "panicQL.editFind(): Not implemented yet"

    def helpContents(self):
        print "panicQL.helpContents(): Not implemented yet"

    def helpAbout(self):
        QMessageBox.about(self,
                          "PANIC Quick-Look Tool",
"""
PQL version: %s\nCopyright (c) 2009-2014 IAA-CSIC  - All rights reserved.\n
Author: Jose M. Ibanez. (jmiguel@iaa.es)
Instituto de Astrofisica de Andalucia, IAA-CSIC

This software is part of PAPI (PANIC Pipeline)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""%__version__)

    def setSFlat_slot(self):
        print "panicQL.setSFlat_slot(): Not implemented yet"

################################################################################


class LoggingConsole (object):
    
    def __init__ (self, textEdit1=None, textEdit2=None, *a, **k):
        
        super (LoggingConsole, self).__init__ (*a, **k)

        self.textEdit_w1 = textEdit1
        self.textEdit_w2 = textEdit2

    def append(self, message='', tag=None):
        
        if self.textEdit_w1 is not None: 
            self.textEdit_w1.append(self.logMsg(message, tag))
        if self.textEdit_w2 is not None:
            self.textEdit_w2.append(self.logMsg(message, tag))
    
    def info(self, message):
        self.append(message, "INFO")
    
    def warning(self, message):
        self.append(message, "WARNING")           

    def error(self, message):
        self.append(message, "ERROR")

    def debug(self, message):
        self.append(message, "DEBUG")
    
    
    def logMsg(self, msg, tag=None):
        """
        Format a msg in order to appear as a log message
        """
        
        if tag=="INFO":
            if self.textEdit_w1 is not None:
                self.textEdit_w1.setTextColor(QColor("blue"))
                self.textEdit_w2.setTextColor(QColor("blue"))
        elif tag=="WARNING":
            if self.textEdit_w1 is not None:
                self.textEdit_w1.setTextColor(QColor("darkRed"))
                self.textEdit_w2.setTextColor(QColor("darkRed"))
        elif tag=="ERROR":
            if self.textEdit_w1 is not None:
                self.textEdit_w1.setTextColor(QColor("red"))
                self.textEdit_w2.setTextColor(QColor("red"))
        elif tag=="DEBUG":
            if self.textEdit_w1 is not None:
                self.textEdit_w1.setTextColor(QColor("darkGreen"))
                self.textEdit_w2.setTextColor(QColor("darkGreen"))
    
            
        s_time = time.strftime("[%Y%m%d-%H:%M:%S] ", time.gmtime())
        #s_time = "["+datetime.datetime.utcnow().isoformat()+"] "
        
        return s_time + " " + msg


################################################################################
# Some functions 
################################################################################
# Because mathOp is aused with the queue of process, it cannot belong to MainGUI
# class or an error in multiprocessing:
#    'objects should only be shared between processes through inheritance'
#            
def mathOp(files, operator, outputFile=None, tempDir=None):
    """
    This method will do the math operation (+,-,/, combine) specified with the 
    input files.

    Note: The type of combinning operation to the pixels is a median after
    a sigma reject algorithm.
    """
    
    log.debug("Start mathOp")
    
    if tempDir==None:
        t_dir = "/tmp"
    else:
        t_dir = tempDir
    
    if outputFile==None:
        output_fd, outputFile = tempfile.mkstemp(suffix='.fits', dir=t_dir)
        os.close(output_fd)

    if operator!='+' and operator!='-' and operator!='/' and operator!='combine':
        log.error("Math operation not supported")
        return None

    try:
        # Remove an old output file (might it happen ?)
        misc.fileUtils.removefiles(outputFile)
        ## MATH operation '+'
        if (operator=='combine' and len(files)>1):
            log.debug("Files to combine = [%s]", files )
            misc.utils.listToFile(files, t_dir+"/files.tmp") 
            # Very important to not scale the frames, because it could 
            # produce wrong combined images due to outliers (bad pixels)
            iraf.mscred.combine(input=("@"+(t_dir+"/files.tmp").replace('//','/')),
                     output=outputFile,
                     combine='median',
                     ccdtype='',
                     reject='sigclip',
                     lsigma=3,
                     hsigma=3,
                     subset='no',
                     scale='none'
                     #masktype='none'
                     #verbose='yes'
                     #scale='exposure',
                     #expname='EXPTIME'
                     #ParList = _getparlistname ('flatcombine')
                     )
        elif (operator=='+' and len(files)>1):
            log.debug("Files to sum = [%s]", files )
            misc.utils.listToFile(files, t_dir+"/files.tmp") 
            # Very important to not scale the frames, because it could 
            # produce wrong combined images due to outliers (bad pixels)
            iraf.mscred.combine(input=("@"+(t_dir+"/files.tmp").replace('//','/')),
                     output=outputFile,
                     combine='sum',
                     ccdtype='',
                     reject='none',
                     lsigma=3,
                     hsigma=3,
                     subset='no',
                     scale='none',
                     weight = 'none'
                     #masktype='none'
                     #verbose='yes'
                     #scale='exposure',
                     #expname='EXPTIME'
                     #ParList = _getparlistname ('flatcombine')
                     )
        ## MATH operation '-,/,*' and just 2 files
        elif len(files)==2:
            iraf.mscarith(operand1=files[0],
                      operand2=files[1],
                      op=operator,
                      result=outputFile,
                      verbose='yes'
                      )
        else:
            log.error("Operation not allowed")
            return None
    except Exception,e:
        log.error("[mathOp] An erron happened while math operation with FITS files")
        raise e
    
    log.debug("mathOp result : %s"%outputFile)
    
    return outputFile
