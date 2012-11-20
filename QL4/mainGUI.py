#! /usr/bin/env python

# Copyright (c) 2008-2012 IAA-CSIC  - All rights reserved. 
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
# Last update 20/Nov/2012
#
################################################################################


    ####################################################################
    #                                                                  #
    #                   This is where it all starts.                   #
    #         The body of this program is at the end of this file      #
    #                                                                  #
    ####################################################################




# system modules
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

# PANIC modules
import reduce
import reduce.calTwFlat
import reduce.calBPM_2
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
#import runQtProcess as QtProc # not need it anymore !!

#Log
import misc.paLog
from misc.paLog import log



# IRAF packages
from pyraf import iraf
from iraf import noao
from iraf import imred
from iraf import ccdred
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

__version__ = "1.0.2"

#-------------------------------------------------------------------------------
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
            
        self.config_opts = config_opts
              
        ## Init member variables
        self.timer_dc = None
        # Init main directories
        if source_dir!= None:
            self.m_default_data_dir = source_dir
        else:      
            self.m_default_data_dir = os.environ['PANIC_DATA']
            
        self.m_sourcedir =  self.m_default_data_dir 
        self.m_outputdir = output_dir 
        self.m_tempdir = temp_dir
        self._ini_cwd = os.getcwd()
        
        # Create LoggingConsole
        self.logConsole = LoggingConsole(self.textEdit_log, self.textEdit_log_2)
        
        # check we have R/W access to temp and out directories
        if not os.access(self.m_tempdir,os.R_OK|os.W_OK):
            log.error("Directory %s has not R/W access",self.m_tempdir)
            self.logConsole.error(str(QString("[WARNING] Directory %1 has not R/W access").arg(self.m_tempdir)))
        if not os.access(self.m_outputdir,os.R_OK|os.W_OK):
            log.error("Directory %s has not R/W access",self.m_outputdir)
            self.logConsole.error(str(QString("[WARNING] Directory %1 has not R/W access").arg(self.m_outputdir)))
            
    
        self.m_frameList_dark = ''
        self.m_frameList_dflat = ''
        self.m_frameList_sflat = ''
        self.m_masterDark = '/tmp/master_dark.fits'
        self.m_masterFlat = '/tmp/master_flat.fits'
        self.m_masterMask = '/tmp/master_mask.fits'


        # Stuff to detect end of an observation sequence to know if data reduction could start
        self.curr_sequence = []  # list having the files of the current sequence received
        self.isOTrunning = True
        self.last_filter = ''  # filter name (J,H,Ks, ...) of the last dataframe received
        self.last_ob_id = -1   # Obseving Block ID (unique) of the last dataframe received
        self.last_img_type = None # image type (DARK, FLAT, SCIENCE, ...) of last received image
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
        self.read_error_files = {} # a dictionary to track the error while reading/detecting FITS files
        
        
        # Default run mode
        if  self.config_opts['quicklook']['run_mode']=="None": 
            self.comboBox_QL_Mode.setCurrentIndex(0)
        elif self.config_opts['quicklook']['run_mode']=="Lazy": 
            self.comboBox_QL_Mode.setCurrentIndex(1)
        elif self.config_opts['quicklook']['run_mode']=="PreReduction": 
            self.comboBox_QL_Mode.setCurrentIndex(2)
        else: 
            self.comboBox_QL_Mode.setCurrentIndex(0)
        
        self.group_by = self.config_opts['general']['group_by'].lower()
            
        
        self.logConsole.info("Welcome to the PANIC QuickLook tool version %s"%__version__)
        

        self.__initializeGUI()
        self.createActions()
        
        ## Init in memory Database
        ## -----------------------
        #datahandler.dataset.initDB()
        #DataBase for input files (in memory)
        self.inputsDB = None
        #DataBase for output files (in memory)
        self.outputsDB = None
        self.__initDBs()
        
        ## Data Collectors initialization
        ## ------------------------------
        self.file_pattern = str(self.lineEdit_filename_filter.text())
        # Data collector for input files
        #self.dc = datahandler.DataCollector("geirs-file2", self.m_sourcedir, 
        #                                    self.file_pattern , self.new_file_func)
        
        if os.path.basename(self.m_sourcedir)=='save_CA2.2m.log': s_mode = 'geirs-file'
        elif os.path.basename(self.m_sourcedir)=='fitsfiles.corrected': s_mode = 'geirs-file2'
        else: s_mode = 'dir'
        
        self.dc = datahandler.DataCollector(s_mode, self.m_sourcedir, 
                                            self.file_pattern , self.new_file_func)
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
        self.connect( self._queue_timer_done, QtCore.SIGNAL("timeout()"), self.checkDoneQueue )
        self._queue_timer_done.start(1000)    # 1 second continuous timer
        
        # Timer for TaskQueue (pending tasks)
        self._queue_timer_todo = QTimer( self )
        self.connect( self._queue_timer_todo, QtCore.SIGNAL("timeout()"), self.TaskRunner )
        self._queue_timer_todo.start(1000)    # 1 second continuous timer
        
        ##Start display (DS9)
        #display.startDisplay()
        #time.sleep(1) # wait until display is up
        
    def __initializeGUI(self):
        """This method initializes some values in the GUI and in the members variables"""

        ## Some signal-connections
        _fromUtf8 = QtCore.QString.fromUtf8
        QtCore.QObject.connect(self.listView_dataS, 
                               QtCore.SIGNAL(_fromUtf8("itemClicked(QTreeWidgetItem*,int)")), 
                               self.selected_file_slot)
        
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
        
    def __initDBs(self):
        """
        This method initializes the in memory DBs used for input and output files
        """
        
        # Inputs DB
        try:
            self.inputsDB = datahandler.dataset.DataSet(None)
            self.inputsDB.createDB()
            #self.inputsDB.load()
        except Exception,e:
            log.error("Error while INPUT data base initialization: \n %s"%str(e))
            raise Exception("Error while INPUT data base initialization")
        
        # Outputs DB
        try:
            self.outputsDB = datahandler.dataset.DataSet(None)
            self.outputsDB.createDB()
            #self.outputsDB.load()
        except Exception,e:
            log.error("Error while OUTPUT data base initialization: \n %s"%str(e))
            raise Exception("Error while OUTPUT data base initialization")    
        
    def new_file_func_out(self, filename):
        """
        Callback used when a new file is detected in output dir
        """
                     
        self.new_file_func(filename, fromOutput=True)
            
    def new_file_func(self, filename, fromOutput=False):
        """ Function executed when a new file is detected into the data source dir or into the out_dir"""
        
        log.debug("New file (in source or out_dir) notified: %s",filename)
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
        #datahandler.dataset.initDB()
        if fromOutput: self.outputsDB.insert(filename)
        else: self.inputsDB.insert(filename)
        ## Query DB
        #(date, ut_time, type, filter, texp, detector_id, run_id, ra, dec, object, mjd)=self.inputsDB.GetFileInfo(filename)
        #fileinfo=self.inputsDB.GetFileInfo(str(dir)+"/"+filename)
        #print "FILEINFO= ", fileinfo
        
        #########################
        # an alternative method to update the ListView; it allows keep the view filter
        self.slot_classFilter() 
        #########################
        ## Update the ListView table
        """elem = QListViewItem( self.listView_dataS )
        elem.setText (0, str(filename))
        elem.setText (1, str(type))
        elem.setText (2, str(filter))
        elem.setText (3, str(texp))
        elem.setText (4, str(date)+"::"+str(ut_time))
        elem.setText (5, str(object))
        elem.setText (6, str(ra))
        elem.setText (7, str(dec))
        """
        ## Update Last frame widget
        self.lineEdit_last_file.setText(str(os.path.basename(filename)))
        self.last_filename = filename

        if fromOutput:
            #nothing else to do ...
            return 

        ## Check if end of observing sequence (science or calibration), then start processing
        end_seq = False
        seq = []
        seqType = ''
        (end_seq, seq, seqType) = self.checkEndObsSequence(filename)
        if end_seq:
            log.debug("Detected end of observing sequence: [%s]"%(seqType))
            self.logConsole.warning("Detected end of observing sequence [%s]"%(seqType))       
        
        
        ##################################################################
        ##  If selected, file or sequence processing will start ....
        #   TODO : !! Not completelly well designed !!!
        ##################################################################
        if self.proc_started:
            if self.comboBox_QL_Mode.currentText()=="None":
                return
            elif self.comboBox_QL_Mode.currentText().contains("Pre-reduction") and end_seq:
                #self.processSeq(seq)
                self.processFiles(seq)
                return
            elif self.comboBox_QL_Mode.currentText().contains("Lazy"):
                self.processLazy(filename)
                return
    
    def processLazy(self, filename):
        """
        Do some  operations to the last file detected. It depend on:
        
            - checkBox_show_imgs
            - checkBox_subDark
            - checkBox_appFlat
            - checkBox_subLastFrame
            
            ToDo:
            - checkBox_appBPM
            - checkBox_subSky
            - checkBox_appSAstrom
            
            Other:
            - checkBox_data_grouping
            - checkBox_doRegrig
            - comboBox_AstromCatalog
            - comboBox_skyWindow
            
        """
        
        log.debug("[processLazy] Starting to process the file %s",filename)
        
        (date, ut_time, type, filter, texp, detector_id, 
         run_id, ra, dec, object, mjd) = self.inputsDB.GetFileInfo(filename)
        
        # ONLY SCIENCE frames will be processed 
        if type!="SCIENCE":
            return

        self.logConsole.info(" [processLazy] Processing file %s "%filename)
        
        # ###########################################################################################
        # According to what options have been selected by the user, we do a processing or other...
        if self.checkBox_subDark.isChecked() or self.checkBox_appFlat.isChecked():
            try:
                log.debug("En-queue the task ....")
                #Look for (last received) calibration files
                mDark, mFlat, mBPM = self.getCalibFor([filename])
                # Both master_dark and master_flat are optional
                if mDark or mFlat:
                    #Change cursor
                    #QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
                    #self.m_processing = False    # Pause autochecking coming files - ANY MORE REQUIRED ?, now using a mutex in thread !!!!
                    #self._task = reduce.ApplyDarkFlat([filename], mDark, mFlat, 
                    #                                  self.m_outputdir)
                    #thread = reduce.ExecTaskThread(self._task.apply, 
                    #                             self._task_info_list)
                    #thread.start()
                    
                    #Put into the queue the task to be done
                    func_to_run = reduce.ApplyDarkFlat([filename], mDark, mFlat, 
                                                     self.m_outputdir)
                    params = ()
                    self._task_queue.put([(func_to_run.apply, params)])
                    
                else:
                    QMessageBox.critical(self, "Error", "Error, cannot find the master calibration files")
            except Exception, e:
                QMessageBox.critical(self, "Error", "Error while processing file.  %s"%str(e))
                #self.m_processing = False
                #self.setCursor(Qt.arrowCursor)
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
                    #Change cursor
                    #QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
                    #self.m_processing = False    # Pause autochecking coming files - ANY MORE REQUIRED ?, now using a mutex in thread !!!!
                    #self._task = self.mathOp
                    #thread = reduce.ExecTaskThread(self._task, 
                    #                               self._task_info_list,  
                    #                               [filename, last_file],
                    #                               '-',None)
                    #thread.start()
                    
                    #Put into the queue the task to be done
                    func_to_run = mathOp 
                    params = ([filename, last_file],'-', None, self.m_tempdir)
                    self._task_queue.put([(func_to_run, params)])
                    
                else:
                    log.debug("Cannot find the previous file to subtract by")
            except Exception,e:
                QMessageBox.critical(self, "Error", "Error while processing file.  %s"%str(e))
                #self.m_processing = False
                #self.setCursor(Qt.arrowCursor)
                raise e
        # ##########################################################################################
        elif self.checkBox_subSky.isChecked():
            try:
                self.subtract_nearSky_slot(True)
            except Exception,e:
                self.m_processing = False # ANY MORE REQUIRED ?,
                raise e    
        # ########################################################################################
        elif self.checkBox_show_imgs.isChecked():
            display.showFrame(filename)
            
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
        master_bpm = [] # we'll get a list of master flat candidates
        
        obj_frame = datahandler.ClFits(sci_obj_list[0])
        # We take as sample, the first frame in the list, but all frames must
        # have the same features (expT,filter,ncoadd, readout-mode, ...)
        expTime = obj_frame.expTime()
        filter = obj_frame.getFilter()
        
        #self.inputsDB.ListDataSet()
        #self.outputsDB.ListDataSet()
        
        #DARK - Do NOT require equal EXPTIME Master Dark ???
        master_dark = self.inputsDB.GetFilesT('MASTER_DARK_MODEL', -1) 
        if len(master_dark)==0 and self.outputsDB!=None:
            master_dark = self.outputsDB.GetFilesT('MASTER_DARK', -1) 
        #FLATS - Do NOT require equal EXPTIME, but FILTER
        master_flat = self.inputsDB.GetFilesT('MASTER_DOME_FLAT', -1, filter)
        if master_flat==[]:
            master_flat = self.inputsDB.GetFilesT('MASTER_TW_FLAT', -1, filter)
        if len(master_flat)==0 and self.outputsDB!=None:
            master_flat = self.outputsDB.GetFilesT('MASTER_DOME_FLAT', -1, filter)
            if len(master_flat)==0:
                master_flat=self.outputsDB.GetFilesT('MASTER_TW_FLAT', -1, filter)

        #BPM                
        master_bpm = self.inputsDB.GetFilesT('MASTER_BPM')
        if len(master_bpm)==0 and self.outputsDB!=None:
            master_bpm = self.outputsDB.GetFilesT('MASTER_BPM')

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
        
    #####################################################
    ### SLOTS ###########################################
    #####################################################

    def findOS_slot(self):
        """ ONLY FOR A TEST - to be removed """
        
        parList,fileList = self.inputsDB.GetSeqFiles()
        k=0
        for seq in parList:
            # create the father
            elem = QTreeWidgetItem( self.listView_OS )
            elem.setText (0, "OB_ID="+str(seq[0])+" OB_PAT="+seq[1]+" FILTER="+seq[2] ) # OB_ID + OB_PAT + FILTER
            #elem.setText (1, str(seq[1])) # OB_PAT
            #elem.setText (2, str(seq[2])) # FILTER
            for file in fileList[k]:
                (date, ut_time, type, filter, texp, detector_id, run_id, ra, 
                 dec, object, mjd) = self.inputsDB.GetFileInfo(file)
                e_child = QTreeWidgetItem (elem)
                e_child.setText (0, str(file))
                e_child.setText (1, str(type))
                e_child.setText (2, str(filter))
                e_child.setText (3, str(texp))
                e_child.setText (4, str(date)+"::"+str(ut_time))
                e_child.setText (5, str(object))
                e_child.setText (6, str(ra))
                e_child.setText (7, str(dec))
            k+=1    
    
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
            self.logConsole.info("+Source : " + dir)
            if (self.m_sourcedir != dir and self.m_outputdir!=dir):
                self.lineEdit_sourceD.setText(dir)
                self.m_sourcedir = str(dir)
                ##Create DataCollector for a path     
                self.file_pattern = str(self.lineEdit_filename_filter.text())
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
                    pass
                else:
                    self.timer_dc = QTimer( self )
                    self.connect( self.timer_dc, QtCore.SIGNAL("timeout()"), self.checkFunc )
                    self.timer_dc.start(1500) ## 1 seconds single-shoot timer
                    
            else:
                #The same dir, nothing to do
                pass
    
    def setDataSourceDir_slot_B(self):
        """
        Deprecated
        """
        
        dir = ""
        dir = QFileDialog.getExistingDirectory( self, "Read existing Directory", 
                                                   self.m_default_data_dir, 
                                                   QFileDialog.ShowDirsOnly
                                                   | QFileDialog.DontResolveSymlinks)
        if (not dir):
            return

        self.lineEdit_sourceD.setText(dir)
        self.m_sourcedir=str(dir)
        filelist=os.listdir(str(dir))
        filelist=fnmatch.filter(filelist, "*.fit*")
       
        for filename in filelist:
            #print "FILENAME= " , filename
            self.inputsDB.insert(str(dir)+"/"+filename)
            (date, ut_time, type, filter, texp, detector_id, run_id, 
             mjd) = self.inputsDB.GetFileInfo(str(dir)+"/"+filename)
            #fileinfo=self.inputsDB.GetFileInfo(str(dir)+"/"+filename)
            #print "FILEINFO= ", fileinfo
            elem = QTreeWidgetItem(self.listView_dataS)
            elem.setText (0, str(filename))
            elem.setText (1, str(type))
            elem.setText (2, str(filter))
            elem.setText (3, str(texp))
            elem.setText (4, str(date)+"::"+str(ut_time))
            
            if self.m_show_imgs: 
                display.showFrame(str(dir)+"/"+filename) 
            
            #self.listView_dataS.insertItem(str(filename))
        self.inputsDB.ListDataSet()

    def autocheck_slot(self):
        """Called when check-button for Input dir is clicked"""
        
        ## Activate or deactivate the autochecking of new files
        if self.checkBox_autocheck.isChecked():
            if self.timer_dc !=None and not self.timer_dc.isActive():
                # Already created in a former user action 
                self.timer_dc.start(1500, False)
            else:
                ##Create QTimer for the data collector
                self.timer_dc = QTimer( self )
                self.connect( self.timer_dc, QtCore.SIGNAL("timeout()"), self.checkFunc )
                self.timer_dc.start(1500) ## 1 seconds single-shot timer
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
        Check Queue of done tasks launched with Process in TaskRunner.
        Funtion called periodically (every 1 sec).
        """
        if not self._done_queue.empty():
            log.debug("Something new in the DoneQueue !")
            try:
                r = self._done_queue.get()
                self.logConsole.debug("[checkDoneQueue] Task finished")
                log.info("Get from Queue: %s"%str(r))
                
                if r!=None:
                    if type(r)==type(list()):
                        if len(r)==0 :
                            self.logConsole.info(str(QString("No value returned")))
                            #QMessageBox.information(self, "Info", "No value returned")
                        else:
                            str_list = ""
                            #print "FILES CREATED=",self._task_info._return
                            display.showFrame(r) #_return is a file list
                            for file in r:
                                #display.showFrame(file)
                                str_list+=str(file)+"\n"
                                #!!! keep up-date the out DB for future calibrations !!!
                                # Because some science sequences could neen the
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
                                    self.outputsDB.insert(file)
                                #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                            self.logConsole.debug(str(QString("%1 files created: \n %1").arg(len(r)).arg(str(str_list))))
                            QMessageBox.information(self,"Info", QString("%1 files created: \n %1").arg(len(r)).arg(str(str_list)))
                    elif type(r)==type("") and os.path.isfile(r):
                        self.logConsole.debug(str(QString(">>New file %1 created ").arg(r)))
                        if r.endswith(".fits"):
                            display.showFrame(r)
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
                    self.logConsole.warning("Nothing returned; processing may be failed !")
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
                                display.showFrame(self._task_info._return) #_return is a file list
                                for file in self._task_info._return:
                                    #display.showFrame(file)
                                    str_list+=str(file)+"\n"
                                    #!!! keep up-date the out DB for future calibrations !!!
                                    # Because some science sequences could neen the
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
                                        self.outputsDB.insert(file)
                                    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                self.logConsole.debug(str(QString("%1 files created: \n %1").arg(len(self._task_info._return)).arg(str(str_list))))
                                QMessageBox.information(self,"Info", QString("%1 files created: \n %1").arg(len(self._task_info._return)).arg(str(str_list)))
                        elif type(self._task_info._return)==type("") and \
                                    os.path.isfile(self._task_info._return):
                            self.logConsole.debug(str(QString(">>New file %1 created ").arg(self._task_info._return)))
                            if self._task_info._return.endswith(".fits"):
                                display.showFrame(self._task_info._return)
                            # Keep updated the out-DB for future calibrations
                            # See comments above
                            if not self.checkBox_outDir_autocheck.isChecked():
                                log.debug("Updating DB...")
                                self.outputsDB.insert(self._task_info._return)
                            #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        else:
                            # Cannot identify the type of the results...for sure
                            # something was wrong...
                            self.logConsole.error("No processing results obtained")
                    else:
                        self.logConsole.info("Nothing returned !")
                else:
                    self.logConsole.error(str(QString("Sequence processing failed \n %1").arg(str(self._task_info._exc))))
                    QMessageBox.critical(self, "Error", "Error while running task. "+str(self._task_info._exc))
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
            
        and display them on the "Calibrations" view Panel
        """
        
        log.debug("Updating master calibration files received")
        
        master_dark = [] # we'll get a list of master dark candidates
        master_flat = [] # we'll get a list of master flat candidates
        master_bpm = [] # we'll get a list of master flat candidates
        
        filter = 'ANY'
        
        # DARK
        # Do NOT require equal EXPTIME Master Dark ???
        master_dark = self.outputsDB.GetFilesT('MASTER_DARK', -1) 
        # FLATS (for any filter)
        master_flat = self.outputsDB.GetFilesT('MASTER_DOME_FLAT', -1, filter)
        if master_flat==[]:
            master_flat = self.outputsDB.GetFilesT('MASTER_TW_FLAT', -1, filter)
        # BPM                
        master_bpm = self.outputsDB.GetFilesT('MASTER_BPM')
        
        """
        log.debug("Master Darks found %s", master_dark)
        log.debug("Master Flats found %s", master_flat)
        log.debug("Master BPMs  found %s", master_bpm)
        """
        # Return the most recently created (according to MJD order)
        if len(master_dark)>0: 
            self.m_masterDark = master_dark[-1]
            self.lineEdit_masterDark.setText(QString(self.m_masterDark))
        #else: self.m_masterDark = None
        if len(master_flat)>0: 
            self.m_masterFlat = master_flat[-1]
            self.lineEdit_masterFlat.setText(QString(self.m_masterFlat))
            filter = self.outputsDB.GetFileInfo(self.m_masterFlat)[3]
            self.lineEdit_masterFlat_Filter.setText(QString(filter))
        #else: self.m_masterFlat = None
        if len(master_bpm)>0: 
            self.m_masterMask = master_bpm[-1]
            self.lineEdit_masterMask.setText(QString(self.m_masterMask))
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
        fits = datahandler.ClFits(filename)
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
            and (fits.getType()==self.last_img_type or self.last_img_type==None):
                #End of sequence
                retSeq = self.curr_sequence[:] # very important !! Lists are mutable objects !
                self.curr_sequence = []
                endSeq, retSeq, typeSeq = True, retSeq, fits.getType()
                log.debug("EOS-0 %s detected : %s"%(typeSeq, str(retSeq)))
            else:
                if self.last_img_type!=None and fits.getType()!=self.last_img_type:
                    #skip/remove the current/last file, type mismatch !
                    self.curr_sequence = self.curr_sequence[:-1]
                    log.error("Detected file type mismatch. File %s skipped"%(filename))
                #Mid of sequence, continue adding file
                endSeq, retSeq, typeSeq = False, self.curr_sequence, fits.getType()
        # ############################################
        # We suppose data is obtained using GEIRS+MIDAS_scripts observations
        # POINT_NO(=OB_ID), DITH_NO, EXPO_NO keyword will be checked for sequece detection
        # TODO: TBC !!! I don't know if it will work with calibration sequences
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
                retSeq = self.curr_sequence[:]  # very important !! Lists are mutable objects ! or clist = copy.copy(list)
                #End of sequence, and then reset the sequence list
                self.curr_sequence = [filename]
                endSeq, retSeq, typeSeq = True, retSeq, fits.getType()
                self.last_ob_id = -1 # re-init
                log.debug("EOS-1 %s detected : %s"%(typeSeq, str(retSeq)))
            else:
                #Check RA,Dec coordinates for distance
                ra_point_distance = self.last_ra-fits.ra
                dec_point_distance = self.last_dec-fits.dec
                dist = math.sqrt((ra_point_distance*ra_point_distance)+(dec_point_distance*dec_point_distance))
                if dist>self.MAX_POINT_DIST:
                    retSeq = self.curr_sequence[:] # very important !! Lists are mutable objects ! 
                    #End of sequence, then reset the sequence list
                    self.curr_sequence = [filename]
                    endSeq,retSeq,typeSeq = True,retSeq,fits.getType()
                    self.last_ob_id = -1 # re-init
                    log.debug("EOS-2 %s detected : %s"%(typeSeq, str(self.retSeq)))
                else:
                    #Mid of sequence, continue adding file
                    self.curr_sequence.append(filename)
                    log.debug("Adding file to sequence: %s",str(self.curr_sequence))
                    endSeq, retSeq, typeSeq = False, self.curr_sequence, fits.getType()
        
        #and finally, before return, update 'last'_values
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
        Funtion called periodically to check for new files (only if no frame 
        is being processed)
        """
        if True:
            self.dc.check()
            if self.checkBox_outDir_autocheck.isChecked():
                self.dc_outdir.check()
                    
    def show_images_slot(self):
        if self.checkBox_show_imgs.isChecked():
            self.m_show_imgs=True
        else:
            self.m_show_imgs=False
            
    def data_grouping_slot(self):
        """
        Slot called when the 'checkBox_data_grouping' is clicked.
        If checked, then data grouping will be done using OB_ID, OB_PAT, 
        FILTER keywords (or whatever we decide at the moment).
        If not checked, tha data grouping will be done using RA,Dec values.
        """
        
        if self.checkBox_data_grouping.isChecked():
            self.lineEdit_ra_dec_near_offset.setEnabled(False)
            self.lineEdit_time_near_offset.setEnabled(False)
        else:
            self.lineEdit_ra_dec_near_offset.setEnabled(True)
            self.lineEdit_time_near_offset.setEnabled(True)      
                  
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
            self.m_tempdir=str(dir)

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
        Remove all files from ListView, DataCollector and DB 
        """
         
        self.listView_dataS.clear()
        self.dc.Clear()
        self.inputsDB.clearDB()


    def add_slot(self):
        """
        Add a new file to the main list panel, but not to the DataCollector 
        list (dc), so it might be already inside 
        """     
        filenames = QFileDialog.getOpenFileNames( self,
                                                "Select one or more files to add",
                                                self.m_default_data_dir,
                                                "FITS files (*.fit*)")

        if not filenames.isEmpty():
            for file in filenames:
                self.new_file_func(str(file))
        
    def del_slot(self):
        """ 
        Delete the current selected file from the main list view panel, 
        but we do not remote from the DataCollector neither file system
        """
        
        self.m_popup_l_sel = []
        listViewItem = QTreeWidgetItemIterator(self.listView_dataS,
                                               QTreeWidgetItemIterator.Selected)
        while listViewItem.value():
            fileName = str(listViewItem.value().text(0))
            indx = self.listView_dataS.indexOfTopLevelItem(listViewItem.value())
            self.listView_dataS.takeTopLevelItem(indx)
            #self.dc.Clear(fileName)
            self.inputsDB.delete( fileName )
            #listViewItem +=1 # because of remove, we do not need an increment
          
        #self.inputsDB.ListDataSetNames()
    
    def display_slot(self):
        
        if self.m_listView_item_selected:
            print "Display image !!"
            display.showFrame(self.m_listView_item_selected)
        else:
            print "Error, no file selected !"
            
    def filename_filter_slot(self):
        """ Modify filename filter for the data collector"""

        new_filter, ok = QInputDialog.getText("New filter","Enter filename filter:")
        if ok and not new_filter.isEmpty():
            self.lineEdit_filename_filter.setText( str(new_filter) )
            self.dc.SetFileFilter( str(new_filter) )

    def slot_classFilter(self):
        """ Filter files on main ListView"""
        
        if self.comboBox_classFilter.currentText()=="GROUP":
            self.listView_dataS.clear()
            sequences = []
            seq_types = []
            
            #Look for sequences
            sequences, seq_types = self.inputsDB.GetSequences(self.group_by) 
            
            #look for un-groupped files
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
                    if file==seq[0]:
                        #the first time, fill the "tittle" of the group 
                        elem.setText(0, "TYPE="+str(seq_types[k]) + 
                                     "  ** FILTER="+str(filter) + 
                                     "  ** TEXP="+str(texp) + 
                                     " ** #imgs="+str(len(seq)))
                        
                    e_child = QTreeWidgetItem(elem)
                    e_child.setText (0, str(file))
                    e_child.setText (1, str(type))
                    e_child.setText (2, str(filter))
                    e_child.setText (3, str(texp))
                    e_child.setText (4, str(date)+"::"+str(ut_time))
                    e_child.setText (5, str(object))
                    e_child.setText (6, str(ra))
                    e_child.setText (7, str(dec))
                k+=1
        else:
            #############################################    
            ## No grouping, only filtering by given type
            #############################################
            self.listView_dataS.clear()
            db = None
            if str(self.comboBox_classFilter.currentText())=="ALL":
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
            else:
                type = str(self.comboBox_classFilter.currentText())
                fileList = self.inputsDB.GetFilesT(type)
                db = self.inputsDB

            elem = None
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
                elem.setText (6, str(ra))
                elem.setText (7, str(dec))
            
            if elem:
                self.listView_dataS.setCurrentItem(elem)
                
            # filtering
            """
            self.listView_dataS.clearSelection()
            it = QListViewItemIterator (self.listView_dataS)
            listViewItem = it.current()
            while listViewItem:
                if listViewItem.text(1).contains(self.comboBox_classFilter.currentText()) or self.comboBox_classFilter.currentText()=="ALL":
                    listViewItem.setVisible(True)
                    listViewItem.setSelectable(True)
                else:
                    listViewItem.setVisible(False)
                    listViewItem.setSelectable(False)
                it+=1
                listViewItem = it.current() 
            """           
#########################################################################
###### Pop-Up ###########################################################
#########################################################################
    def createActions(self):
        """
        Create actions used in the Pop-Up menu (contextMenuEvent)
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
            shortcut="Ctrl+P",
            statusTip="Build Bad Pixel Map with selected files", 
            triggered=self.createBPM_slot)

        self.subOwnSkyAct = QtGui.QAction("Subtract self-sky", self,
            shortcut="Ctrl+S",
            statusTip="Subtract self-sky", 
            triggered=self.subtract_ownSky_slot)
        
        self.subNearSkyAct = QtGui.QAction("Subtract near-Sky", self,
            shortcut="Ctrl+N",
            statusTip="Subtract near sky", 
            triggered=self.subtract_nearSky_slot)
        
        self.quickRedAct = QtGui.QAction("Quick-Reduction", self,
            shortcut="Ctrl+Q",
            statusTip="Run quick data reduction of selected files", 
            triggered=self.do_quick_reduction_slot)
        
        self.stackFrameAct = QtGui.QAction("Stack Frames", self,
            shortcut="Ctrl+F",
            statusTip="Stack selected files", 
            triggered=self.createStackedFrame_slot)
        
        self.astroAct = QtGui.QAction("Astrometric Calib.", self,
            shortcut="Ctrl+a",
            statusTip="Run astrometric calibration of selected file", 
            triggered=self.do_raw_astrometry)
        
        self.photoAct = QtGui.QAction("Photometric Calib.", self,
            shortcut="Ctrl+p",
            statusTip="Run photometric calibration of selected file", 
            triggered=self.do_raw_photometry)
        
        self.statsAct = QtGui.QAction("Statistics", self,
            shortcut="Ctrl+s",
            statusTip="Get statistincs of selected file", 
            triggered=self.show_stats_slot)
        
        self.fwhmAct = QtGui.QAction("FWHM", self,
            shortcut="Ctrl+W",
            statusTip="Get FWHM selected file", 
            triggered=self.fwhm_estimation_slot)
        
        self.bckgroundAct = QtGui.QAction("Background estimation", self,
            shortcut="Ctrl+b",
            statusTip="Get image background of selected file", 
            triggered=self.background_estimation_slot)
            
        # Sub-menu Math actions
        self.subAct = QtGui.QAction("&Subtract Images", self,
            shortcut="Ctrl+-",
            statusTip="Subtract selected files", 
            triggered=self.subtractFrames_slot)
        
        self.sumAct = QtGui.QAction("Sum Images", self,
            shortcut="Ctrl++",
            statusTip="Sum selected files", 
            triggered=self.sumFrames_slot)
        
        self.divAct = QtGui.QAction("Div Images", self,
            shortcut="Ctrl+/",
            statusTip="Divide selected files", 
            triggered=self.divideFrames_slot)

        # Sub-menu fits actions       
        self.splitAct = QtGui.QAction("Split MEF", self,
            shortcut="Ctrl+m",
            statusTip="Split currect selected MEF file", 
            triggered=self.splitMEF_slot)
        
        self.joinAct = QtGui.QAction("Join MEF", self,
            shortcut="Ctrl+j",
            statusTip="Join currect selected files into a MEF", 
            triggered=self.joinMEF_slot)
        
        # Group-menu actions
        self.redSeqAct = QtGui.QAction("Reduce Sequence", self,
            shortcut="Ctrl+r",
            statusTip="Reduce the selected sequence", 
            triggered=self.reduceSequence_slot)
        
        self.copyGAct = QtGui.QAction("&Copy files to clipboard", self,
            shortcut="Ctrl+c",
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
            
        #print "LIST=", self.m_popup_l_sel
        
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
            popUpMenu = QMenu()
            popUpMenu.addAction(self.dispAct)
            popUpMenu.addAction(self.copyAct)
            popUpMenu.addSeparator()
            popUpMenu.addAction(self.mDarkAct)
            popUpMenu.addAction(self.mDFlatAct)
            popUpMenu.addAction(self.mDTwFlatAct)
            popUpMenu.addAction(self.mGainMapAct)
            popUpMenu.addAction(self.mBPMAct)
            popUpMenu.addSeparator()
            popUpMenu.addAction(self.subOwnSkyAct)
            popUpMenu.addAction(self.subNearSkyAct)
            popUpMenu.addAction(self.quickRedAct)
            popUpMenu.addAction(self.stackFrameAct)
            popUpMenu.addSeparator()
            popUpMenu.addAction(self.astroAct)
            popUpMenu.addAction(self.photoAct)
            popUpMenu.addSeparator()
            popUpMenu.addAction(self.statsAct)
            popUpMenu.addAction(self.fwhmAct)
            popUpMenu.addAction(self.bckgroundAct)
            
            newAction = popUpMenu.addAction("Math")
            subMenu_arith = QtGui.QMenu("Popup Submenu", self)
            subMenu_arith.addAction(self.subAct)
            subMenu_arith.addAction(self.sumAct)
            subMenu_arith.addAction(self.divAct)
            newAction.setMenu(subMenu_arith)
            
            newAction = popUpMenu.addAction("FITS")
            subMenu_fits = QtGui.QMenu("Popup Submenu", self)
            subMenu_fits.addAction(self.splitAct)
            subMenu_fits.addAction(self.joinAct)
            newAction.setMenu(subMenu_fits)
            

            ## Disable some menu items depeding of the number of item selected 
            ## in the list view
            if len(self.m_popup_l_sel)==1:
                self.copyAct.setEnabled(False)
                self.mDarkAct.setEnabled(False)
                self.mDFlatAct.setEnabled(False)
                self.mDTwFlatAct.setEnabled(False)
                self.mGainMapAct.setEnabled(False)
                self.mBPMAct.setEnabled(False)
                self.subAct.setEnabled(False)
                self.sumAct.setEnabled(False)
                self.divAct.setEnabled(False)
            elif len(self.m_popup_l_sel)>1:
                self.dispAct.setEnabled(False)
                self.mBPMAct.setEnabled(False)
                self.subOwnSkyAct.setEnabled(False)
            elif len(self.m_popup_l_sel)>2:
                self.subAct.setEnabled(False)
                self.sumAct.setEnabled(False)

            if len(self.m_popup_l_sel)<5:
                self.subNearSkyAct.setEnabled(False)
            
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
            return self.processFiles(group_files)
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
        allow copy&paste
        """
        
        text = ""
        for filename in self.m_popup_l_sel:
            text = text + str(filename) + "\n"
        
        clipboard = QApplication.clipboard()
        clipboard.setText(text)    
            
    def splitMEF_slot(self):
        """
        Split each MEF selected file from the list view into NEXT separate 
        single FITS file, where NEXT is number of extensions.
        As result, NEXT files should be created
        """
        for file in self.m_popup_l_sel:
            try:
                mef = misc.mef.MEF([file])
                mef.doSplit(".%02d.fits")
            except Exception,e:
                log.debug("Cannot split file %s. Maybe it's not a MEF file", str(e))
                QMessageBox.critical(self, "Error", "Cannot split file : %s \n Maybe it's not a MEF file"%(file))
                #self.logConsole.info(QString(str(line)))
    
    def joinMEF_slot(self):
        """
        Join/stitch a MEF file to stitch them together
        As result, one single FITS file with the stitched image should be created
        """
        for file in self.m_popup_l_sel:
            try:
                mef = misc.mef.MEF([file])
                mef.doJoin(".join.fits")
            except Exception,e:
                log.debug("Cannot stitch file %s. Maybe it's not a MEF file",str(e))
                QMessageBox.critical(self, "Error", "Cannot stitch file: %s. \n Maybe it's not a MEF file"%(file))
    
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
        Imexam the currect filename selected - not used because it fails 
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

    def iraf_console_slot(self):
        """Open a IRAF console session"""
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
            
        #os.system("/usr/local/bin/ds9 %s &" %((self.m_listView_item_selected)))
        self.logConsole.info("DS9 launched !")
        
    def start_aladin_slot(self):
        """Start Aladin tool"""

        os.system('echo "load %s ;sync; get vizier(2mass)" |/usr/local/bin/aladin -nobanner &' %(self.m_listView_item_selected))
        self.logConsole.info("Aladin launched !")
        
        # utils.runCmd does not allow launch in background !!

    def pushB_sel_masteDark_slot(self):
        """
        Called to create a master dark - to be deprecated
        """
        source = QFileDialog.getOpenFileName( self,
                                              self.m_default_data_dir, 
                                              "Select master dark file",  
                                              "(*.fit*")

        if str(source):
            self.lineEdit_masterDark.setText(source)
            self.m_masterDark = str(source)

    def pushB_sel_masterFlat_slot(self):
        """
        Called to create a master flat - to be deprecated
        """
        source = QFileDialog.getOpenFileName( self,
                                              self.m_default_data_dir, 
                                              "Select master flat file",  
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
        source=QFileDialog.getOpenFileName( self,
                                            self.m_default_data_dir, 
                                            "Select Master BPM", 
                                            "(*.pl*, *.fit*)")

        if str(source):
            self.lineEdit_masterMask.setText(source)
            self.m_masterMask = str(source)                       

    def pushB_subtract_last2_slot(self):
        """
        Slot called when button 'Subtract-Last2' is clicked.
        
        The two newest SCIENCE frames received (MJD) will be subtracted
        """
        
        #(type , filter) = self.inputsDB.GetFileInfo(self.last_filename)[2:4]
        #if type != 'SCIENCE':
        #only SCIENCE frames are processed
        #pass
        ltemp = self.inputsDB.GetFilesT('SCIENCE') # (mjd sorted)
        if len(ltemp)>1:
            #get the  last 2 files (included the current one)
            last2_files = ltemp[-2:]
            if (self.inputsDB.GetFileInfo(last2_files[0])[3] == 
                self.inputsDB.GetFileInfo(last2_files[1])[3]):
                try:
                    #Change cursor
                    #QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
                    #self.m_processing = False    # Pause autochecking coming files - ANY MORE REQUIRED ?, now using a mutex in thread !!!!
                    #thread = reduce.ExecTaskThread(self.mathOp, 
                    #                               self._task_info_list, 
                    #                               last2_files,'-', 
                    #                               "/tmp/sub.fits")
                    #thread.start()
                    
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
        """This method is called to subtract two images selected from the File List View"""


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
                    thread = reduce.ExecTaskThread(self.mathOp, 
                                                   self._task_info_list, 
                                                   self.m_popup_l_sel,'-', 
                                                   str(outFilename))
                    thread.start()
                except:
                    QMessageBox.critical(self, "Error", 
                                         "Error while subtracting files")
                    raise
        
    def sumFrames_slot(self):
        """This methot is called to sum two images selected from the File List View"""

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
                    thread = reduce.ExecTaskThread(self.mathOp, 
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
                    thread=reduce.ExecTaskThread(self.mathOp, 
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
        Called when Create button pressed in tab 'Calibrations'
        """

        listItems = self.listBox_darks.findItems(QString("*"), 
                                                 Qt.MatchWrap | Qt.MatchWildcard)
        # Build list
        l_list = [ str(item.text()) for item in listItems]

        if len(l_list)<3:
            QMessageBox.information(self,"Info","Not enough frames !")
            return

        outfileName = QFileDialog.getSaveFileName(self,
                                                  "Choose a filename to save under",
                                                  self.m_outputdir + "/master_dark.fits", 
                                                  "*.fits")
        if not outfileName.isEmpty():
            try:
                texp_scale = False
                QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
                self._task = reduce.calDark.MasterDark(l_list, 
                                                       self.m_tempdir, 
                                                       str(outfileName), 
                                                       texp_scale)
                thread = reduce.ExecTaskThread(self._task.createMaster, 
                                               self._task_info_list)
                thread.start()
            except Exception, e:
                self.setCursor(Qt.arrowCursor)
                QMessageBox.critical(self, "Error", 
                                     "Error while creating master Dark. \n"+str(e))
                raise e
        else:
            pass

    
    def do_quick_reduction_slot(self):
        """ 
        Run a quick reduction of the current user selected science files 
        in the list view panel.
        """
            
        if len(self.m_popup_l_sel)<5:
            QMessageBox.information(self,"Info","Error, not enough frames selected to reduce (>=5) !")
            return
              
        #Create working thread that compute sky-frame
        try:
            self.processFiles(self.m_popup_l_sel)
        except Exception,e:
            QMessageBox.critical(self, "Error", "Error while Quick data reduction: \n%s"%str(e))
            raise e

    def createMasterDFlat_slot(self):
        
        """
        Create a master dome flat field from selected frames
        """
        
        listItems = self.listBox_domeF.findItems(QString("*"), 
                                                 Qt.MatchWrap | Qt.MatchWildcard)
        # Build list
        l_list = [ str(item.text()) for item in listItems]

        if len(l_list)>2:
            outfileName = QFileDialog.getSaveFileName(self,
                                                      "Choose a filename to save under",
                                                      self.m_outputdir + "/master_dflat.fits", 
                                                      "*.fits")
            if not outfileName.isEmpty():
                try:
                    QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
                    self._task = reduce.calDomeFlat.MasterDomeFlat(l_list, 
                                                                   self.m_tempdir, 
                                                                str(outfileName))
                    thread=reduce.ExecTaskThread(self._task.createMaster, 
                                                 self._task_info_list)
                    thread.start()
                except:
                    self.setCursor(Qt.arrowCursor)
                    QMessageBox.critical(self, "Error", "Error while creating master Dome Flat")
                    raise
        else:
            QMessageBox.information(self,"Info","Error, not enough frames (>2) !")

    def createMasterTwFlat_slot(self):
      
        """
        Create a master twilight flat field from selected frames
        """
        listItems = self.listBox_SkyF.findItems(QString("*"), 
                                                 Qt.MatchWrap | Qt.MatchWildcard)
        # Build list
        l_list = [ str(item.text()) for item in listItems]
        
        if len(l_list)>2:
            outfileName = QFileDialog.getSaveFileName(self,
                                                      "Choose a filename to save under",
                                                      self.m_outputdir+"/master_twflat.fits", 
                                                      "*.fits") 
            if not outfileName.isEmpty():
                try:
                    QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
                    self._task = reduce.calTwFlat.MasterTwilightFlat (l_list, 
                                                                      self.m_masterDark, 
                                                                      str(outfileName))
                    thread = reduce.ExecTaskThread(self._task.createMaster, 
                                                 self._task_info_list)
                    thread.start()
                except:
                    self.setCursor(Qt.arrowCursor)
                    log.error("Error creating master Twilight Flat file")
                    raise
    
    def createGainMap_slot(self):

        """ 
        Create a gain-map using the own science files selected on the main list view
        
        TODO: Check the selected images are SCIENCE frames with same filter !!! and expTime, .....
         
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
                    self.setCursor(Qt.arrowCursor)
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
                    self.setCursor(Qt.arrowCursor)
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
                    self.setCursor(Qt.arrowCursor)
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
            sex_config = self.config_opts['config_files']['sextractor_conf']
            sex_param = self.config_opts['config_files']['sextractor_param']
            sex_conv = self.config_opts['config_files']['sextractor_conv']
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
            
            # Next implementation does not work because sex.run() doesn't return
            # the output file created !
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
                display.showFrame(out_file)
            finally:
                self.setCursor(Qt.arrowCursor)
        else:
            log.error("Sorry, selected file does not look a science file  ")
            QMessageBox.information(self, "Info", "Sorry, selected file does not look a science file ")                             
                         
    def subtract_nearSky_slot(self, last_files=False):
        """ Subtract nearest sky using skyfiler from IRDR package"""
      
        if self.checkBox_data_grouping.isChecked():
            log.error("Sorry, data grouping with POINT_NO, DITH_NO, EXPO_NO not yet implemented ")
            QMessageBox.information(self, "Info", "Sorry, data grouping with POINT_NO, DITH_NO, EXPO_NO not yet implemented ")
            return
        
        # #####################################
        # Look for files in search-radius mode
        # #####################################
        if not last_files:
            # Compute search-radius of nearest frames
            ra_dec_near_offset = self.lineEdit_ra_dec_near_offset.text().toInt()[0]/3600.0
            time_near_offset = self.lineEdit_time_near_offset.text().toInt()[0]/86400.0
            print "RA_DEC_OFFSET", ra_dec_near_offset
            print "TIME_NEAR_OFFSET", time_near_offset
          
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
                        file_n=p
                    else: 
                        p+=1
                
                # Ask user validation
                res = QMessageBox.information(self, "Info", 
                                              QString("Selected near frames are:\n %1").arg(my_list), 
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
                file_n = len(near_list)
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
                                          out_file=self.m_outputdir+"/skysub.fits",
                                          obs_mode="dither", dark=None, 
                                          flat=self.m_masterFlat,
                                          bpm=None, red_mode="quick",
                                          group_by=self.group_by, 
                                          check_data=True, 
                                          config_dict=self.config_opts)
            
            thread = reduce.ExecTaskThread(self._task.subtractNearSky, 
                                           self._task_info_list, 
                                           self._task.rs_filelist, file_n)
            thread.start()
        except:
            #Anyway, restore cursor
            # Although it should be restored in checkLastTask, could happend an exception while creating the class RS,
            # thus the ExecTaskThread can't restore the cursor
            QApplication.restoreOverrideCursor() 
            QMessageBox.critical(self, "Error", "Error while subtracting near sky")
            raise 
              
    def createBPM_slot(self):
        """ Create a Bad Pixel Mask from a set of selected files (flats)
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


    def show_stats_slot(self):
        """Show image statistics in the log console of the files selected"""
        
        self.logConsole.info("FILE                            MEAN         MODE       STDDEV       MIN       MAX  ")
        for item in self.m_popup_l_sel:
            values = (iraf.mscstat (images=item,
            fields="image,mean,mode,stddev,min,max",format='no',Stdout=1))
            for line in values:
                #line=os.path.basename(values[0])
                self.logConsole.info(str(QString(str(line))))
        
    def background_estimation_slot(self):
        """ Give an background estimation of the current selected image """
        
        cq = reduce.checkQuality.CheckQuality(self.m_popup_l_sel[0])
        try:     
            img=cq.estimateBackground(self.m_outputdir+"/bckg.fits")
            
            values = (iraf.mscstat (images=img,
            fields="image,mean,mode,stddev,min,max",format='yes',Stdout=1))
            #file,mean,mode,stddev,min,max=values[0].split()
            self.logConsole.info(str(QString("Background estimation :")))
            for line in values:
                self.logConsole.info(str(line))
                #self.logConsole.info(QString("Background estimation MEAN= %1   MODE=%2    STDDEV=%3    MIN=%4         MAX=%5").arg(mean).arg(mode).arg(stddev).arg(min).arg(max))
            
            display.showFrame(img)
        except Exception,e:
            self.logConsole.error("ERROR: something wrong while computing background")
            raise e
          
    def fwhm_estimation_slot (self):
        """ Give an FWHM estimation of the current selected image """
        
        cq = reduce.checkQuality.CheckQuality(self.m_popup_l_sel[0])
        try:
            fwhm,std=cq.estimateFWHM()
            if fwhm>0:
                self.logConsole.info(str(QString("FWHM = %1 (pixels) std= %2").arg(fwhm).arg(std)))
            else:
                self.logConsole.error("ERROR: Cannot estimage FWHM with selected image")           
        
        except Exception,e:
            self.logConsole.error("ERROR: something wrong while computing background")
            raise e
        
    def createStackedFrame_slot(self):
        """ 
        Compute a stacked frame (subtract sky, shift and align) from a set of 
        nearest (ra,dec, mjd) frames, selected by user or automatic search.
        Actually, it is basically a quick-reduction !!! 
        (see do_quick_reduction_slot), with the only difference that here we 
        can look for the files to be reduced in automatic way or use the ones 
        selected  by the user (as in do_quick_reduction_slot)
        """
        
        file_n = 0 # really, not used for the moment
        file_list=[]
        
        ra_dec_near_offset = self.lineEdit_ra_dec_near_offset.text().toInt()[0]/3600.0
        time_near_offset = self.lineEdit_time_near_offset.text().toInt()[0]/86400.0
        print "RA_DEC_OFFSET", ra_dec_near_offset
        print "TIME_NEAR_OFFSET", time_near_offset
        
        # CASE 1: Automatic search for nearest frames (ra, dec, mjd) 
        if len(self.m_popup_l_sel)==1:
            if self.checkBox_data_grouping.isChecked():
                log.error("Sorry, data grouping with POINT_NO, DITH_NO, EXPO_NO not yet implemented ")
                QMessageBox.information(self, "Info", "Sorry, data grouping with POINT_NO, DITH_NO, EXPO_NO not yet implemented ")
                return
            QMessageBox.information(self, "Info", "Only one file was selected, automatic file grouping will be done.")
            fits = datahandler.ClFits(self.m_listView_item_selected)
            if fits.getType()=='SCIENCE':
                near_list = self.inputsDB.GetFiles('ANY', 'SCIENCE', -1, fits.filter, fits.mjd, fits.ra, fits.dec,  ra_dec_near_offset*2, time_near_offset, runId=0)
                print "NEAR_LIST=", near_list
                # For the moment, the minimun number of nearest is >0
                if len(near_list)==0:
                    QMessageBox.information(self, "Info", "Not enough science frames found")  
                    return
            else:
                QMessageBox.information(self, "Info", "Selected frame is not a science frame") 
                return
            
            # Create nearest file list (ar,dec,mjd) from current selected science file
            view_list = ""
            i=1
            for file in near_list:
                if file==self.m_listView_item_selected: 
                    file_n=i
                else: 
                    i=i+1
                if (datahandler.ClFits(file).getType()!='SCIENCE'):
                    QMessageBox.critical(self, "Error", 
                                         QString("File %1 is not a science frame").arg(file))
                else:
                    file_list.append(file)
                    view_list +=  file + "\n"
            
            resp=QMessageBox.information(self, "Info", 
                                         QString("Selected near frames are:\n %1").arg(view_list),
                                         QMessageBox.Ok, QMessageBox.Cancel)
            if resp==QMessageBox.Cancel:
                return
                    
        # CASE 2: Stack frames selected by user in the list_view
        elif len(self.m_popup_l_sel)>1 and len(self.m_popup_l_sel)<4:
            QMessageBox.information(self,"Info","Error, not enough frames selected to reduce (>4) !")
            return    
        elif len(self.m_popup_l_sel)>=4:
            # Create file list from current selected science files
            for file in self.m_popup_l_sel :
                if (datahandler.ClFits(file).getType()!='SCIENCE'):
                    QMessageBox.critical(self, "Error", QString("File %1 is not a science frame").arg(file))
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
                self.processFiles(file_list, group_by='filter') 
            except Exception, e:
                QMessageBox.critical(self, "Error", "Error while building stack")
                raise e
        
        
    def createSuperMosaic_slot(self):
        """
        TODO : here is only a scheme of the routine !!!!
        Create an stitched wide-mosaic of next frames
        """
        QMessageBox.information(self, "Info", "Sorry, funtion not yet implemented !")
        
        # Next code is only a raw template 
        # Check list lenght
        if len(self.m_popup_l_sel)<1:
            QMessageBox.information(self, "Info", "Not enough science pre-reduced frames selected")
            return
        
        if len(self.m_popup_l_sel)==1:
            QMessageBox.information(self, "Info", "Only one file was selected, automatic file grouping will be done.")
            return
            # NEED TO BE DONE
        # Stack frame user selected  
        elif len(self.m_popup_l_sel)>1:
            # Create file list from current selected science files
            listfile = self.m_tempdir+"/mosaic_files.list"      
            file_lst = open( listfile, "w" )
            for file in self.m_popup_l_sel :
                if (datahandler.ClFits(file).getType()!='SCIENCE_REDUCED'):
                    QMessageBox.critical(self, "Error", QString("File %1 is not a science reduced frame").arg(file))
                    file_lst.close()
                    return
                else:
                    file_lst.write( file +"\n")
            file_lst.close()
        
        # Call external script (SWARP)
        os.chdir(self.m_tempdir)
        #Change to working directory
        os.chdir(self.m_tempdir)
        #Restore cursor
        QApplication.restoreOverrideCursor()

        # Prepare SWARP call
        swarp = astromatic.SWARP()
        swarp.config['CONFIG_FILE'] = self.config_dict['config_files']['swarp_conf']
        basename, extension = os.path.splitext(self.m_popup_l_sel[0])
        swarp.ext_config['HEADER_SUFFIX'] = extension + ".head"  # very important !
        if not os.path.isfile(self.m_popup_l_sel[0]+".head"):
            raise Exception ("Cannot find required .head file")
            
        swarp.ext_config['COPY_KEYWORDS'] = 'OBJECT,INSTRUME,TELESCOPE,IMAGETYP,FILTER,FILTER1,FILTER2,SCALE,MJD-OBS,HISTORY'
        swarp.ext_config['IMAGEOUT_NAME'] = os.path.dirname(self.m_popup_l_sel[0]) + "/coadd_tmp.fits"
        if os.path.isfile(basename + ".weight" + extension):
            swarp.ext_config['WEIGHT_TYPE'] = 'MAP_WEIGHT'
            swarp.ext_config['WEIGHT_SUFFIX'] = '.weight' + extension
            swarp.ext_config['WEIGHTOUT_NAME'] = os.path.dirname(self.m_popup_l_sel[0]) + "/coadd_tmp.weight.fits"
        
        """
        #IMPORTANT: Rename the .head files in order to be found by SWARP
        #but it's much efficient modify the HEADER_SUFFIX value (see above)
        for aFile in self.input_files:
            basename, extension = os.path.splitext(aFile)
            shutil.move(aFile+".head", basename +".head")  # very important !!
        """    
        
        #Create working thread that compute sky-frame
        try:
            updateconfig = False
            clean = False
            thread = reduce.ExecTaskThread(swarp.run, 
                                           self._task_info_list,
                                           #args of the doAstrometry function
                                           self.m_popup_l_sel, updateconfig, clean) 
            thread.start()
        except:
            QMessageBox.critical(self, "Error", "Error while running SWARP")
            raise
                
      
      
    def do_raw_astrometry(self):
        """
        @summary: Compute an astrometric solution for the selected file in the 
        main list view panel.
        Basically, it is SCAMP is called and compared to 2MASS catalog 
        (or GSC2.3, USNO-B) to compute astrometric calibration. 
        A new WCS is added to the header.
        
        @return: None, but produce an astrometrically calibrated image.
          
        @todo: we should check the image is pre-reduced or at least, with 
        enough objects 
        """
        
        if len(self.m_popup_l_sel)==1:
            fits = datahandler.ClFits(self.m_listView_item_selected)
            if fits.getType()=='SCIENCE': 
                ## Run astrometry parameters
                out_file = self.m_outputdir+"/"+os.path.basename(self.m_listView_item_selected.replace(".fits",".wcs.fits"))
                # Catalog
                if self.comboBox_AstromCatalog.currentText().contains("2MASS"): catalog="2MASS"
                elif self.comboBox_AstromCatalog.currentText().contains("USNO-B1"): catalog="USNO-B1"
                elif self.comboBox_AstromCatalog.currentText().contains("GSC-2.2"): catalog="GSC-2.2"
                else: catalog="2MASS"
                
                #Change to working directory
                os.chdir(self.m_tempdir)
                #Change cursor
                QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
                #Create working thread that compute sky-frame
                try:
                    thread = reduce.ExecTaskThread(reduce.astrowarp.doAstrometry, 
                                                   self._task_info_list,
                                                   #args of the doAstrometry function
                                                   self.m_listView_item_selected, 
                                                   out_file, catalog, 
                                                   self.config_opts)
                    thread.start()
                except:
                    QMessageBox.critical(self, "Error", "Error while subtracting near sky")
                    raise
            else:
                QMessageBox.information(self,"Info", QString("Sorry, but you need a reduced science frame."))
        
    def do_raw_photometry(self):
        """
        @summary: Compute an rough photometric solution for the selected file 
        in the main list view panel.
        Basically, It is compared to 2MASS catalog (or other) to compute photometric
        calibration, and obtain a zero-point estimation.

        @return: None, but produce an photometric fitting curve and a Zero point
        estimation.
          
        @todo: we should check the image is pre-reduced or at least, with 
        enough objects 
        """
        
        if len(self.m_popup_l_sel)==1:
            fits = datahandler.ClFits(self.m_listView_item_selected)
            if fits.getType()=='SCIENCE': 
                ## Run astrometry parameters
                out_file = self.m_outputdir+"/" + \
                    os.path.basename(self.m_listView_item_selected.replace(".fits","_photo.pdf"))
                # Catalog
                if self.comboBox_AstromCatalog.currentText().contains("2MASS"): catalog="2MASS"
                elif self.comboBox_AstromCatalog.currentText().contains("USNO-B1"): catalog="USNO-B1"
                elif self.comboBox_AstromCatalog.currentText().contains("GSC-2.2"): catalog="GSC-2.2"
                else: catalog="2MASS"
                
                #Change to working directory
                os.chdir(self.m_tempdir)
                #Change cursor
                QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
                #Create working thread that compute sky-frame
                try:
                    thread = reduce.ExecTaskThread(photo.photometry.doPhotometry,
                                                   self._task_info_list,
                                                   #args of the doPhotometry function
                                                   self.m_listView_item_selected,
                                                   catalog,
                                                   out_file)
                    thread.start()
                except:
                    QMessageBox.critical(self, "Error", "Error while subtracting near sky")
                    raise
            else:
                QMessageBox.information(self,"Info", QString("Sorry, but you need a reduced science frame."))
    
    def createCalibs_slot(self):
        
        """
        Given the current data set files, compute the whole master calibration files
        """
        
        fileList = self.inputsDB.GetFilesT("ANY")
        
        if len(fileList)>1:
            #Change cursor
            QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
            
            try:
                self._task = RS.ReductionSet( fileList, self.m_outputdir, out_file=self.m_outputdir+"/red_result.fits", \
                                            obs_mode="dither", dark=None, flat=None, bpm=None, red_mode="quick",\
                                            group_by=self.group_by, check_data=True, config_dict=self.config_opts )
                
                thread=reduce.ExecTaskThread(self._task.buildCalibrations, self._task_info_list)
                thread.start()
            except:
                #Anyway, restore cursor
                # Although it should be restored in checkLastTask, could happend an exception while creating the class RS,
                # thus the ExecTaskThread can't restore the cursor
                QApplication.restoreOverrideCursor() 
                QMessageBox.critical(self, "Error", "Error while building  master calibrations files")
                raise
    
    def testSlot(self):
        """
        Not implemeted
        """
        pass
    
    def pushB_start_stop_slot(self):
        """
        Method called when "Start/Stop-processing button is clicked
        """
        
        # Change the button color and label
        if self.proc_started == True:
            self.proc_started = False
            self.pushButton_start_proc.setText("Start Processing")
            self.pushButton_start_proc.setPaletteBackgroundColor(QColor(205,64,64))
        else:    
            self.proc_started = True
            self.pushButton_start_proc.setText("Stop Processing")
            self.pushButton_start_proc.setPaletteBackgroundColor(QColor(34,139,34))
            
            # Ask if we with to process only the next files coming or also the
            # current ones in the source_dir and then the new ones
            res = QMessageBox.information(self, "Info", 
                                QString("Process also the current files in the \
source directory (If no, only the new ones coming will be processed) ?"), 
                                QMessageBox.Ok, QMessageBox.Cancel)
            
            if res==QMessageBox.Cancel:
                return
            else:
                self.processFiles()
            
 
        """
        it = QListViewItemIterator (self.listView_OS)
        listViewItem = it.current()
        while listViewItem: 
            if listViewItem.isSelected():
                print "1st CHILD=", listViewItem.firstChild()
                
            it+=1
            listViewItem = it.current()
        """
        """
        return 
        
        name="Pepe"
        msg = QString ("Me llamo %1").arg( name )
        QMessageBox.critical(self,"Error",msg)

        self.listView_1.addColumn( "C_X" );
        elem = QListViewItem( self.listView_1 )
        elem.setText(0,"HOLA")
        elem.setText(1,"ADIOS")
        
        item =QStyleSheetItem( self.textEdit_log.styleSheet(), "mytag" )
        item.setColor( QColor("red") )
        item.setFontWeight( QFont.Bold )
        item.setFontUnderline( True )
        self.textEdit_log.append("HOLA me llamo <mytag> Jose Miguel </mytag>")
        """
    
    
    def worker1(self, input, output):
        #print "arg_1=",input.get()
        for func, args in iter(input.get, 'STOP'):
            #result = calculate(func, args)
            print "ARGS=", args
            output.put(RS.ReductionSet(args).func())
            
    def processFiles(self, files=None, group_by=None):
        """
        @summary: Process the files provided; if any files were given, all the files 
        in the current Source List View (but not the output files) will be
        processed. The processing task will be inserted into the _task_queue
        to be processed as soon as the current (if any) processing has finished
        (m_processing=False).  
        
        @param files: files list to be processed; if None, all the files in the
        current List View will be used !
        
        @note: When the RS object is created, we provide the outputDB as the 
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
            
        #Create working thread that process the files
        try:
            # generate a random filename for the master, to ensure we do not overwrite any file
            output_fd, outfilename = tempfile.mkstemp(suffix='.fits', prefix='redObj_', dir=self.m_outputdir)
            os.close(output_fd)
            os.unlink(outfilename) # we only need the name
            ###self._task = RS.ReductionSet(files, self.m_outputdir, out_file=outfilename,
            ###                                obs_mode="dither", dark=None, 
            ###                                flat=None, bpm=None, red_mode="quick",
            ###                                group_by=self.group_by, check_data=True,
            ###                                config_dict=self.config_opts,
            ###                                external_db_files=self.outputsDB.GetFiles())
            # provide the outputDB files as the external calibration files for the RS 
            
            log.debug("ReductionSet created !")
            # New approach using multiprocessing
            #rs_filelist, out_dir=None, out_file=None, obs_mode="dither", 
            #     dark=None, flat=None, bpm=None, red_mode="quick", 
            #     group_by="ot", check_data=True, config_dict=None, 
            #     external_db_files=None, temp_dir = None,
            
            if group_by==None: l_group_by = self.group_by
            else: l_group_by = group_by
            
            params = [(files, self.m_outputdir, outfilename,
                      "dither", None, 
                      None, None, "quick",
                      l_group_by, True,
                      self.config_opts,
                      self.outputsDB.GetFiles(), None)]
            
            func_to_run = RS.ReductionSet(*(params[0])).reduceSet
            self._task_queue.put([(func_to_run,())])
            #self._task_queue.put(params) #default function supposed!
            
            ##Process(target=self.worker, 
            ##        args=(self._task_queue, self._done_queue)).start()
                    
        except Exception,e:
            QMessageBox.critical(self, "Error", "Error while processing Obs. Sequence: \n%s"%str(e))
            raise e # Para que seguir elevando la excepcion ?
        
    def TaskRunner(self):
        """
        Procedure that continuisly in checking the task of pending task to be done
        """
        if not self._task_queue.empty() and not self.m_processing:
            log.debug("Something new in the TaskQueue !")
            try:
                self.logConsole.debug("Starting to process queued task")
                #Change cursor
                QApplication.setOverrideCursor(QCursor(Qt.WaitCursor))
                #'mutex' variable 
                self.m_processing = True
                Process(target=self.worker, 
                    args=(self._task_queue, self._done_queue)).start()
            except Exception,e:
                log.error("Error in task Runner: %s"%(str(e)))
                self.logConsole.debug("Error in TaskRunner")
                self.m_processing = False
                self.setCursor(Qt.arrowCursor)
            finally:
                log.debug("End of TaskRunner")
                 
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
        Callback function used by Process task
        """
        
        func, args = input.get()[0]
        print "FUNC=",func
        print "ARGS=",args
        
        log.debug("[worker] Worker start to work ...")
        try:
            output.put(func(*args))
            log.info("[worker] task done !")
        except Exception,e:
            log.error("[worker] Error while processing task\n %s"%str(e))
            output.put(None) # the DoneQueue will detect it
        finally:
            self.m_processing = False
            log.debug("Worker finished its task !")
    
    def mathOp(self,files, operator, outputFile=None, tempDir=None):
        """
        This method will do the math operation (+,-,/) specified with the 
        input files.
        """
        
        log.debug("Start mathOp")
        
        if tempDir==None:
            print "TEMP DIR is None !!"
            t_dir = "/tmp"
        else:
            t_dir = tempDir
        
        if outputFile==None:
            output_fd, outputFile = tempfile.mkstemp(suffix='.fits', dir=t_dir)
            os.close(output_fd)
    
        if operator!='+' and operator!='-' and operator!='/':
            log.error("Math operation not supported")
            return None
    
        try:
            # Remove an old output file (might it happen ?)
            misc.fileUtils.removefiles(outputFile)
            ## MATH operation '+'
            if (operator=='+' and len(files)>2):
                log.debug("Frame list to combine = [%s]", files )
                misc.utils.listToFile(files, t_dir+"/files.tmp") 
                iraf.mscred.combine(input=("@"+(t_dir+"/files.tmp").replace('//','/')),
                         output=outputFile,
                         combine='average',
                         ccdtype='',
                         reject='sigclip',
                         lsigma=3,
                         hsigma=3,
                         subset='no',
                         scale='mode'
                         #masktype='none'
                         #verbose='yes'
                         #scale='exposure',
                         #expname='EXPTIME'
                         #ParList = _getparlistname ('flatcombine')
                         )
            ## MATH operation '-,/,*'
            else:
                iraf.mscarith(operand1=files[0],
                          operand2=files[1],
                          op=operator,
                          result=outputFile,
                          verbose='yes'
                          )
        except Exception,e:
            log.error("[mathOp] An erron happened while math operation with FITS files")
            raise e
        
        log.debug("mathOp result : %s"%outputFile)
        
        return outputFile

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

    def helpIndex(self):
        print "panicQL.helpIndex(): Not implemented yet"

    def helpContents(self):
        print "panicQL.helpContents(): Not implemented yet"

    def helpAbout(self):
        print "panicQL.helpAbout(): Not implemented yet"

    def setSFlat_slot(self):
        print "panicQL.setSFlat_slot(): Not implemented yet"

################################################################################


class LoggingConsole (object):
    
    def __init__ (self, textEdit1=None, textEdit2=None, *a, **k):
        
        super (LoggingConsole, self).__init__ (*a, **k)

        self.textEdit_w1 = textEdit1
        self.textEdit_w2 = textEdit2
        """
        ## Create log tags
        if self.textEdit_w1 is not None:
            # Error
            item = QStyleSheetItem( self.textEdit_w1.styleSheet(), "error_tag" )
            item.setColor( QColor("red") )
            item.setFontWeight( QFont.Bold )
            item.setFontUnderline( False )
            # Warning
            item = QStyleSheetItem( self.textEdit_w1.styleSheet(), "warning_tag" )
            item.setColor( QColor("blue") )
            item.setFontWeight( QFont.Bold )
            item.setFontUnderline( False )
            # Info
            item = QStyleSheetItem( self.textEdit_w1.styleSheet(), "info_tag" )
            item.setColor( QColor("black") )
            item.setFontWeight( QFont.Bold )
            item.setFontUnderline( False )
            # Debug
            item = QStyleSheetItem( self.textEdit_w1.styleSheet(), "debug_tag" )
            item.setColor( QColor("#1C731C") ) # dark green
            item.setFontWeight( QFont.Bold )
            item.setFontUnderline( False )
        
        if self.textEdit_w2 is not None:
            # Error
            item = QStyleSheetItem( self.textEdit_w2.styleSheet(), "error_tag" )
            item.setColor( QColor("red") )
            item.setFontWeight( QFont.Bold )
            item.setFontUnderline( False )
            # Warning
            item = QStyleSheetItem( self.textEdit_w2.styleSheet(), "warning_tag" )
            item.setColor( QColor("blue") )
            item.setFontWeight( QFont.Bold )
            item.setFontUnderline( False )
            # Info
            item = QStyleSheetItem( self.textEdit_w2.styleSheet(), "info_tag" )
            item.setColor( QColor("black") )
            item.setFontWeight( QFont.Bold )
            item.setFontUnderline( False )
            # Debug
            item = QStyleSheetItem( self.textEdit_w2.styleSheet(), "debug_tag" )
            item.setColor( QColor("#1C731C") ) # dark green
            item.setFontWeight( QFont.Bold )
            item.setFontUnderline( False )
        """

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
#Because mathOp is used with the queue of process, it cannot belong to MainGUI
#class or an error in multiprocessing:
#  'objects should only be shared between processes through inheritance'
#            
def mathOp(files, operator, outputFile=None, tempDir=None):
    """
    This method will do the math operation (+,-,/) specified with the 
    input files.
    """
    
    log.debug("Start mathOp")
    
    if tempDir==None:
        print "TEMP DIR is None !!"
        t_dir = "/tmp"
    else:
        t_dir = tempDir
    
    if outputFile==None:
        output_fd, outputFile = tempfile.mkstemp(suffix='.fits', dir=t_dir)
        os.close(output_fd)

    if operator!='+' and operator!='-' and operator!='/':
        log.error("Math operation not supported")
        return None

    try:
        # Remove an old output file (might it happen ?)
        misc.fileUtils.removefiles(outputFile)
        ## MATH operation '+'
        if (operator=='+' and len(files)>2):
            log.debug("Frame list to combine = [%s]", files )
            misc.utils.listToFile(files, t_dir+"/files.tmp") 
            iraf.mscred.combine(input=("@"+(t_dir+"/files.tmp").replace('//','/')),
                     output=outputFile,
                     combine='average',
                     ccdtype='',
                     reject='sigclip',
                     lsigma=3,
                     hsigma=3,
                     subset='no',
                     scale='mode'
                     #masktype='none'
                     #verbose='yes'
                     #scale='exposure',
                     #expname='EXPTIME'
                     #ParList = _getparlistname ('flatcombine')
                     )
        ## MATH operation '-,/,*'
        else:
            iraf.mscarith(operand1=files[0],
                      operand2=files[1],
                      op=operator,
                      result=outputFile,
                      verbose='yes'
                      )
    except Exception,e:
        log.error("[mathOp] An erron happened while math operation with FITS files")
        raise e
    
    log.debug("mathOp result : %s"%outputFile)
    
    return outputFile