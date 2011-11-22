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
# mainGUI (implement main GUI functions for PANIC QL-pipeline)
#
# mainGUI.py
#
# Last update 13/Sep/2011
#
################################################################################


    ####################################################################
    #                                                                  #
    #                   This is where it all starts.                   #
    #         The body of this program is at the end of this file      #
    #                                                                  #
    ####################################################################


# system modules
from qt import *
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
from panicQL import *
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
from runQtProcess import *
#Log
import misc.paLog
from misc.paLog import log


# Interaction with FITS files
import pyfits

# IRAF packages
from pyraf import iraf
from iraf import noao
from iraf import imred
from iraf import ccdred
from iraf import mscred

# Math module for efficient array processing
import numpy

  
class MainGUI(panicQL):

            
    def __init__(self, source_dir="/tmp/data", output_dir='/tmp/out', 
                 temp_dir='/tmp', config_opts=None):               
        
        panicQL.__init__(self)
        
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
        self.last_img_type = '' # image type (DARK, FLAT, SCIENCE, ...) of last received image
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
        if  self.config_opts['quicklook']['run_mode']=="None": self.comboBox_QL_Mode.setCurrentItem(0)
        elif self.config_opts['quicklook']['run_mode']=="Lazy": self.comboBox_QL_Mode.setCurrentItem(1)
        elif self.config_opts['quicklook']['run_mode']=="PreReduction": self.comboBox_QL_Mode.setCurrentItem(2)
        else: self.comboBox_QL_Mode.setCurrentItem(0)
        
        self.group_by = self.config_opts['general']['group_by'].lower()
            
        
        self.logConsole.info("Wellcome to the PANIC QuickLook tool v1.0" )
        

        self.__initializeGUI()
        
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
        self.dc = datahandler.DataCollector("geirs-file2", self.m_sourcedir, 
                                            self.file_pattern , self.new_file_func)
        # Data collector for output files
        self.dc_outdir = None # Initialized in checkOutDir_slot()
        #datahandler.DataCollector("dir", self.m_outputdir, self.file_pattern, self.new_file_func_out)
        
        ## Task management
        ## ---------------
        self._task = None                       # Pointer to task/thread that is running 
        #self._task_event = threading.Event()    # Event to syncronize ExecTaskThread and ExecTaskThread
        self._task_info = None
        self._task_info_list = []               # Queue-list where task status is saved
        
        self._task_timer = QTimer( self )
        self.connect( self._task_timer, SIGNAL("timeout()"), self.checkLastTask )
        self._task_timer.start( 1000, False )    # 1 second continuous timer
        
        ##Start display (DS9)
        #display.startDisplay()
        #time.sleep(1) # wait until display is up
        
    def __initializeGUI(self):
        """This method initializes some values in the GUI and in the members variables"""

        ## Init Panel widgets values
        self.lineEdit_sourceD.setText(self.m_sourcedir)
        self.lineEdit_outputD.setText(self.m_outputdir)
        self.lineEdit_tempD.setText(self.m_tempdir)
        ## Init calibration files
        self.lineEdit_masterDark.setText(QString(self.m_masterDark))
        self.lineEdit_masterFlat.setText(QString(self.m_masterFlat))
        self.lineEdit_masterMask.setText(QString(self.m_masterMask))
        
        
        ## PRUEBAS !!!! #####
        #self.listView_config.setItem(1,1, QComboBox())
        elem = QCheckListItem( self.listView_config, "check1" )
        elem.setRenameEnabled(0,1)
        elem.setText (0, "test")
        #self.listView_config.setRenameEnabled(1,0)
        self.listView_config.setSorting(-1)
        ### FIN DE PRUEBAS #####
    
    def __initDBs(self):
        """This method initializes the in memory DBs used for input and output files"""
        
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
        """Callback used when a new file is detected in output dir"""
                     
        self.new_file_func(filename, fromOutput=True)
        
            
    def new_file_func(self, filename, fromOutput=False):
        """ Function executed when a new file is detected into the data source dir or into the out_dir"""
        
        
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
        
        #######################################
        ### otherwise, it must be a new file !!
        log.debug( "New file detected --> %s. Going to verification...", filename)
        # (Try) to check if the FITS file writing has finished -- 
        try:
            temp = pyfits.open(filename,"readonly", ignore_missing_end=True) # since some problems with O2k files 
            #temp.verify(option='exception')  # it is a too severe checking !!!
            temp.close()
            misc.utils.fits_simple_verify(filename)
            
            if filename in self.read_error_files:
                log.debug("New try to read %s was successful ! ", filename)
                #self.read_error_files.remove(filename)
                del self.read_error_files[filename] 
            temp.close()
        except Exception,e:
            log.error("Error while opening file %s. Maybe writing is not finished: %s", filename, str(e))
            if filename not in self.read_error_files:
                log.debug("Removing file from DC.dirlist %s", filename)
                if fromOutput: self.dc_outdir.remove(filename) # it means the file could be detected again by the DataCollector
                else: self.dc.remove(filename)
                #self.read_error_files.append(filename) # to avoid more than two tries of opening the file
                self.read_error_files[filename]=1 # init the error counter
            else:
                if self.read_error_files[filename]>self.MAX_READ_ERRORS:
                    #definitely, we discard the file
                    log.error("Definitely discarted file %s", filename)
                    self.logConsole.error(str(QString("File %1 corrupted and discarted").arg(filename)))
                    del self.read_error_files[filename] # could be interesting to have a list of discarted files ??
                    #self.read_error_files.remove(filename)
                else:
                    log.debug("File %s still has some read error ...will retry to read it again ...", filename)
                    self.read_error_files[filename]+=1
                    if fromOutput: self.dc_outdir.remove(filename) # it means the file could be detected again by the DataCollector
                    else: self.dc.remove(filename) 
            return
           
        self.logConsole.append("New file (verified) detected in source: " + filename)
        
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
            
 
 
        """    
        elif self.comboBox_QL_Mode.currentText().contains("Lazy"):
            if self.m_show_imgs:
                display.showFrame(filename)
        #QL-Mode: Pre-reduction
        elif self.comboBox_QL_Mode.currentText().contains("Pre-reduction"):
            ## Process
            self.processSeq(seq)
            return
        #QL-Mode: UserDef_1
        elif self.comboBox_QL_Mode.currentText().contains("UserDef_1"):
            return ## TO BE DONE
        #QL-Mode: UserDef_2
        elif self.comboBox_QL_Mode.currentText().contains("UserDef_2"):
            return ## TO BE DONE
        #QL-Mode: UserDef_3
        elif self.comboBox_QL_Mode.currentText().contains("UserDef_3"):
            return ## TO BE DONE
        """
    
    def processSeq(self, obsSequence):
        
        """ DEPRECATED !!!
        
            Process the observing sequence received according with the QL pipeliene recipes
            The sequence could be a calib or science sequence.
            
            @param obsSequence: a list files belonging to the observing sequence
            
            @deprecated: can be replaced by processFiles
             
        """
        
        log.debug("Starting to process Observation Sequence...")
        self.logConsole.info("++ Starting to process Observation Sequence :")
        
        for file in obsSequence:
            self.logConsole.info("     - " + file)
            #self.logConsole.info(QString("    - %1").arg(file))
        self.logConsole.info("... processing sequence ...")
            
        #Change to working directory
        os.chdir(self.m_tempdir)
        
        #Change cursor
        self.setCursor(Qt.waitCursor)
        
        # Pause autochecking coming files - ANY MORE REQUIRED ?, now using a 
        # mutex in thread !!!!
        self.m_processing = False    
        
        #Create working thread that process the obsSequence
        try:
            # generate a random filename for the master, to ensure we do not overwrite any file
            output_fd, outfilename = tempfile.mkstemp(suffix='.fits', prefix='redObj_', dir=self.m_outputdir)
            os.close(output_fd)
            os.unlink(outfilename) # we only need the name
            self._task = RS.ReductionSet( obsSequence, self.m_outputdir, out_file=outfilename, \
                                            obs_mode="dither", dark=None, flat=None, bpm=None, red_mode="quick", \
                                            group_by=self.group_by, check_data=True, config_dict=self.config_opts)
            
            if self._task.isaCalibSet():
                log.debug("It's a calib sequence what is going to be reduced !")
                self.logConsole.info("It is a CALIBRATION sequence")
                #thread = reduce.ExecTaskThread(self._task.buildCalibrations, self._task_info_list)
                thread = reduce.ExecTaskThread(self._task.reduceSet, self._task_info_list, "quick")
            else:
                log.debug("It's a science sequence what is going to be reduced !")
                self.logConsole.info("It is a SCIENCE sequence")
                thread = reduce.ExecTaskThread(self._task.reduceSet, self._task_info_list, "quick")
            thread.start()
        except Exception,e:
            #Anyway, restore cursor
            # Although it should be restored in checkLastTask, could happend an exception while creating the class RS,
            # thus the ExecTaskThread can't restore the cursor
            self.setCursor(Qt.arrowCursor) 
            QMessageBox.critical(self, "Error", "Error while processing Obs. Sequence: \n%s"%str(e))
            self.m_processing = False
            raise e

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
        
        log.debug("Starting to process the file %s",filename)
        
        (date, ut_time, type, filter, texp, detector_id, 
         run_id, ra, dec, object, mjd) = self.inputsDB.GetFileInfo(filename)
        
        # ONLY SCIENCE frames will be processed 
        if type!="SCIENCE":
            return

        self.logConsole.info(" ... processing file %s..."%filename)
        
        # ###########################################################################################
        # According to what options have been selected by the user, we do a processing or other...
        if self.checkBox_subDark.isChecked() or self.checkBox_appFlat.isChecked():
            #Change to working directory
            os.chdir(self.m_tempdir)  # -- required ???
            #Create working thread that process the file
            try:
                # Master dark with ANY EXPT, and scaled later
                # could be > 1 master darks, but used the last(mjd sorted)
                master_dark = self.inputsDB.GetFilesT('MASTER_DARK', -1) 
                # could be > 1 master flat (tw,dome), but used the last(mjd sorted)
                master_flat = self.inputsDB.GetFilesT('MASTER_TW_FLAT',-1, filter) 
                # Both master_dark and master_flat are optional
                if len(master_dark)>0 or len(master_flat)>0:
                    if len(master_dark)>0: mDark = master_dark[-1] # most recently
                    else: mDark = None
                    if len(master_flat)>0: mFlat = master_flat[-1] # most recently
                    else: mFlat = None
                    #Change cursor
                    self.setCursor(Qt.waitCursor)
                    self.m_processing = False    # Pause autochecking coming files - ANY MORE REQUIRED ?, now using a mutex in thread !!!!
                    self._task = reduce.ApplyDarkFlat([filename], mDark, mFlat, 
                                                      self.m_outputdir)
                    thread = reduce.ExecTaskThread(self._task.apply, 
                                                 self._task_info_list)
                    thread.start()
                else:
                    QMessageBox.critical(self, "Error", "Error, cannot find the master calibration files")
            except:
                QMessageBox.critical(self, "Error", "Error while processing file.  %s"%str(e))
                self.m_processing = False
                self.setCursor(Qt.arrowCursor)
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
                    self.setCursor(Qt.waitCursor)
                    self.m_processing = False    # Pause autochecking coming files - ANY MORE REQUIRED ?, now using a mutex in thread !!!!
                    self._task = self.mathOp
                    thread = reduce.ExecTaskThread(self._task, 
                                                   self._task_info_list,  
                                                   [filename, last_file],
                                                   '-',None)
                    thread.start()
                    
                else:
                    log.debug("Cannot find the previous file to subtract by")
            except:
                QMessageBox.critical(self, "Error", "Error while processing file.  %s"%str(e))
                self.m_processing = False
                self.setCursor(Qt.arrowCursor)
                raise e
        # ##########################################################################################
        elif self.checkBox_subSky.isChecked():
            print "QUE PASSSSSSSS !!!!"
            try:
                self.subtract_nearSky_slot(True)    
            except Exception,e:
                self.m_processing = False # ANY MORE REQUIRED ?,
                raise e    
        # ########################################################################################
        elif self.checkBox_show_imgs.isChecked():
            display.showFrame(filename)
            
        
    #####################################################
    ### SLOTS ###########################################
    #####################################################

    def findOS_slot(self):
        """ ONLY FOR A TEST - to be removed """
        
        parList,fileList = self.inputsDB.GetSeqFiles()
        k=0
        for seq in parList:
            # create the father
            elem = QListViewItem( self.listView_OS )
            elem.setText (0, "OB_ID="+str(seq[0])+" OB_PAT="+seq[1]+" FILTER="+seq[2] ) # OB_ID + OB_PAT + FILTER
            #elem.setText (1, str(seq[1])) # OB_PAT
            #elem.setText (2, str(seq[2])) # FILTER
            for file in fileList[k]:
                (date, ut_time, type, filter, texp, detector_id, run_id, ra, 
                 dec, object, mjd) = self.inputsDB.GetFileInfo(file)
                e_child = QListViewItem(elem)
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
        
        source = QFileDialog.getExistingDirectory( self.m_default_data_dir, 
                                                 self,"get existing directory", 
                                                 "Choose a directory",True )
        
        #source=QFileDialog.getOpenFileNames( "Source log (*.log)", self.m_default_data_dir, self, "Source Dialog","select source")
        ## NOTE: 'dir' can be a file or a directory
        
        if (not source):
            return
        else:
            #dir=str(source[0])
            dir = str(source)
            self.logConsole.info("+Source : " + dir)
            if (self.m_sourcedir != dir and self.m_outputdir!=dir):
                self.lineEdit_sourceD.setText(dir)
                self.m_sourcedir = str(dir)
                ##Create DataCollector for a path     
                self.file_pattern = str(self.lineEdit_filename_filter.text())
                if os.path.isfile(dir):
                    self.dc = datahandler.DataCollector("geirs-file", str(dir), self.file_pattern , self.new_file_func)  
                elif os.path.isdir(dir):
                    self.dc = datahandler.DataCollector("dir", str(dir), self.file_pattern , self.new_file_func)
                    #self.dc = datahandler.DataCollector("geirs-file2", str(dir), self.file_pattern , self.new_file_func)	
                ## Activate the autochecking of new files
                self.checkBox_autocheck.setChecked(True)
                ##Create QTimer for the data collector
                if self.timer_dc !=None and self.timer_dc.isActive():
                    # Already created in a former user action 
                    pass
                else:
                    self.timer_dc = QTimer( self )
                    self.connect( self.timer_dc, SIGNAL("timeout()"), self.checkFunc )
                    self.timer_dc.start( 1500, False ) ## 1 seconds single-shoot timer
                    
            else:
                #The same dir, nothing to do
                pass
    
    def setDataSourceDir_slot_B(self):
        dir=""
        dir=QFileDialog.getExistingDirectory( self.m_default_data_dir, self,
                                             "get existing directory", "Choose a directory",True )
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
            elem = QListViewItem( self.listView_dataS )
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
                self.connect( self.timer_dc, SIGNAL("timeout()"), self.checkFunc )
                self.timer_dc.start( 1500, False ) ## 1 seconds single-shot timer
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
            
    def checkLastTask(self):
        """
        Receiver function signaled by QTimer 'self._task_timer' object.
        Funtion called periodically (every 1 sec) to check last task results.
        """
        
        if len(self._task_info_list)>0:
            self.logConsole.debug("Task finished")
            try:
                self._task_info = self._task_info_list.pop()
                if self._task_info._exit_status == 0: # EXIT_SUCCESS, all was OK
                    self.logConsole.info("Process successful finished  !")
                    if self._task_info._return!=None:
                        if type(self._task_info._return)==type(list()):
                            if len(self._task_info._return)==0:
                                self.logConsole.info(str(QString("Any value returned")))
                                QMessageBox.information(self, "Info", "Any value returned")
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
                            self.logConsole.error("Any processing results obtained")
                else:
                    self.logConsole.error(str(QString("Sequence processing failed \n %1").arg(str(self._task_info._exc))))
                    QMessageBox.critical(self, "Error", "Error while running task. "+str(self._task_info._exc))
            except Exception,e:
                raise Exception("Error while checking_task_info_list: %s"%str(e))
            finally:
                #Anyway, restore cursor
                self.setCursor(Qt.arrowCursor)
                self.m_processing = False
                # Return to the previus working directory
                os.chdir(self._ini_cwd)
                #Good moment to update the master calibration files    
                self._update_master_calibrations()

            
        
    def _update_master_calibrations(self):
        """
        Query the **outputsDB** to update the master calibrations files with the last 
        calibration files received, and then update: 
            
            self.m_masterDark
            self.m_masterFlat
            self.m_masterMask
            
        and display them on the "Calibrations" view Panel
        """
        
        log.debug("Updating master calibration files")
        
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
        
        Note: It even works for calibration frame sequences (dark, flats, ...)

        @return: a triplet as follow:
            - True if EOS was found, otherwise False
            - if True, the list of files of the sequence
            - if True, the type of sequence (DARK, FLAT, SCIENCE, etc ...)
            
        """
        
        endSeq=False
        retSeq=[]
        typeSeq=''
        
        # Read the FITS file
        fits = datahandler.ClFits(filename)
        #only for debug !!
        log.info("Current_FILTER= %s, Last_FILTER=%s, Current_OB_ID=%s, Last_OB_ID=%s",
                 fits.getFilter(), self.last_filter, fits.getOBId(), self.last_ob_id )
        #############################################
        # Based on the meta-data provided by the OT
        # Option based on number of expositions (NOEXP) in the pattern and the exposition 
        # number (EXPNO) of the current frame; It even works for calibration frames sequences (dark, flats, ...)
        log.info("EXPNO= %s, NOEXPO= %s", fits.getExpNo(), fits.getNoExp())
        if fits.isFromOT():
            log.debug("Checking OT keywords...")
            # General case for OT observations
            if fits.getExpNo()==fits.getNoExp() and fits.getExpNo()!=-1:
                self.curr_sequence.append(filename)
                retSeq = self.curr_sequence[:] # very important !! Lists are mutable objects !
                self.curr_sequence = []
                endSeq, retSeq, typeSeq = True, retSeq, fits.getType()
                log.debug("EOS-0 %s detected : %s"%(typeSeq, str(retSeq)))
            else:
                self.curr_sequence.append(filename)
                endSeq, retSeq, typeSeq = False, self.curr_sequence, fits.getType()
        # ############################################
        # We suppose data is obtained using GEIRS+MIDAS_scripts observations
        # POINT_NO, DITH_NO, EXPO_NO keyword will be checked for sequece detection
        # TODO: TBC !!! I don't know if it will work with calibration sequences
        else:
            log.debug("Checking GEIRS keywords...")
            if self.last_ob_id==-1: # first time 
                self.curr_sequence.append(filename)
                log.debug("Adding file to sequence: %s",str(self.curr_sequence))
                endSeq, retSeq, typeSeq = False, self.curr_sequence, fits.getType()
            elif fits.getOBId()!=self.last_ob_id \
                or fits.getFilter()!= self.last_filter \
                or fits.getType()!=self.last_img_type:
                retSeq = self.curr_sequence[:]  # very important !! Lists are mutable objects ! or clist = copy.copy(list)
                #reset the sequence list
                self.curr_sequence = [filename]
                endSeq, retSeq, typeSeq = True, retSeq, fits.getType()
                log.debug("EOS-1 %s detected : %s"%(typeSeq, str(retSeq)))
            else:
                ra_point_distance = self.last_ra-fits.ra
                dec_point_distance = self.last_dec-fits.dec
                dist = math.sqrt((ra_point_distance*ra_point_distance)+(dec_point_distance*dec_point_distance))
                if dist>self.MAX_POINT_DIST:
                    retSeq = self.curr_sequence[:] # very important !! Lists are mutable objects ! 
                    #reset the sequence list
                    self.curr_sequence = [filename]
                    endSeq,retSeq,typeSeq = True,retSeq,fits.getType()
                    log.debug("EOS-2 %s detected : %s"%(typeSeq, str(self.retSeq)))
                else:
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
        Funtion called periodically to check for new files (only if no frame is being processed)
        """
        #print "---------------------->m_processing=", self.m_processing
        if not self.m_processing:
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
        If checked, then data grouping will be done using OB_ID, OB_PAT, FILTER keywords
        (or whatever we decide at the moment).
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
       
        dir=QFileDialog.getExistingDirectory( self.m_outputdir, self,
                                             "get existing directory", "Choose a directory",True )
        
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
        """Select Temp Directory for temporal files while processing"""
        dir = QFileDialog.getExistingDirectory( self.m_default_data_dir, self,
                                             "get existing directory", 
                                             "Choose a directory",True )
        # CAREFUL: <sourcedir> must be DISTINT  to <tempdir>, infinite loop
        if dir and self.m_sourcedir!=dir:
            self.lineEdit_tempD.setText(dir)
            self.m_tempdir=str(dir)

    def setDarks_slot(self):
        #dir=QFileDialog.getExistingDirectory( self.m_default_data_dir, self,
        #                                     "get existing directory", "Choose a directory",True )
        filelist=QFileDialog.getOpenFileNames( "FITS files (*.fit*)", self.m_default_data_dir, self, "FileDialog","select dark files")
        
        print "FL=", str(filelist.first())
        a=[]
        for fs in filelist:
            a.append(str(fs))
            
        self.m_frameList_dark=a
        
        if not filelist.isEmpty():
            for filename in filelist:
                self.listBox_darks.insertItem(str(filename))

    def setDFlats_slot(self):
        
        filenames=QFileDialog.getOpenFileNames( "FITS files (*.fit*)", 
                                                self.m_default_data_dir, 
                                                self, "FileDialog", 
                                                "select files")

        if not filenames.isEmpty():
            for file in filenames:
                self.listBox_domeF.insertItem(str(file))  

    def clear_mainlist_slot(self):
      
        """ Remove all files from ListView, DataCollector and DB """
         
        self.listView_dataS.clear()
        self.dc.Clear()
        self.inputsDB.clearDB()


    def add_slot(self):
      
        """
        Add a new file to the main list panel, but not to the DataCollector 
        list (dc), so it might be already inside 
        """     
        filenames = QFileDialog.getOpenFileNames( "FITS files (*.fit*)", 
                                                  self.m_default_data_dir, 
                                                  self, "FileDialog",
                                                  "select files")

        if not filenames.isEmpty():
            for file in filenames:
                self.new_file_func(str(file))
        
    def del_slot(self):
        
        """ 
        Delete the current selected file from the main list view panel, 
        but we do not remote from the DataCollector neither file system
        """
        
        self.m_popup_l_sel = []
        it=QListViewItemIterator(self.listView_dataS)
        listViewItem = it.current()
        while listViewItem: 
            if listViewItem.isSelected():
                fileName=str(listViewItem.text(0))
                self.listView_dataS.takeItem(listViewItem)
                #self.dc.Clear(fileName)
                self.inputsDB.delete( fileName )
                ##self.m_popup_l_sel.append(str(listViewItem.text(0)))
                print "SELECTED TO REMOVE=", fileName
            else:
                it+=1
            listViewItem = it.current()
          
        #self.inputsDB.ListDataSetNames()
    
    def display_slot(self):
      
        if self.m_listView_item_selected:
            display.showFrame(self.m_listView_item_selected)

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
            """
            files = self.inputsDB.GetFiles()
            # Now, we instantiate a ReductionSet, only to get the sequences, but
            # not for data reduction purposes. 
            try:
                rs = RS.ReductionSet( files, self.m_outputdir, out_file=None, \
                                            obs_mode="dither", dark=None, 
                                            flat=None, bpm=None, red_mode="quick", \
                                            group_by=self.config_opts['general']['group_by'], 
                                            check_data=True, config_dict=self.config_opts)
                
                sequences, seq_types = rs.getSequences()
                del rs # we do not need anymore
            except Exception, e:
                log.error("Error while creating Reduction set")
                raise e
            """
            #look for sequences
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
                elem = QListViewItem( self.listView_dataS )
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
                        
                    e_child = QListViewItem(elem)
                    e_child.setText (0, str(file))
                    e_child.setText (1, str(type))
                    e_child.setText (2, str(filter))
                    e_child.setText (3, str(texp))
                    e_child.setText (4, str(date)+"::"+str(ut_time))
                    e_child.setText (5, str(object))
                    e_child.setText (6, str(ra))
                    e_child.setText (7, str(dec))
                k+=1
        elif self.comboBox_classFilter.currentText()=="GROUP_OLD":
            self.listView_dataS.clear()
            parList,fileList = self.inputsDB.GetSeqFiles()
            k = 0
            for seq in parList:
                # create the father
                elem = QListViewItem( self.listView_dataS )
                # OB_ID + OB_PAT + FILTER
                elem.setText(0, "OB_ID="+str(seq[0]) + " ** OB_PAT=" + 
                             str(seq[1])+" ** FILTER=" + str(seq[2]) + 
                             " ** #imgs="+str(len(fileList[k])) ) 
                #elem.setText (1, str(seq[1])) # OB_PAT
                #elem.setText (2, str(seq[2])) # FILTER
                for file in fileList[k]:
                    (date, ut_time, type, filter, texp, detector_id, 
                     run_id, ra, dec, object, mjd) = self.inputsDB.GetFileInfo(file)
                    e_child = QListViewItem(elem)
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
            db=None
            if str(self.comboBox_classFilter.currentText())=="ALL":
                fileList = self.inputsDB.GetFilesT("ANY")
                db=self.inputsDB
                #db.ListDataSet()
            elif  str(self.comboBox_classFilter.currentText())=="OUTS":
                fileList = self.outputsDB.GetFilesT("ANY")
                db=self.outputsDB
                #db.ListDataSet()
            elif str(self.comboBox_classFilter.currentText())=="SKY_FLAT":
                fileList = self.inputsDB.GetFilesT("SKY_FLAT")
                db=self.inputsDB
            else:
                type=str(self.comboBox_classFilter.currentText())
                fileList = self.inputsDB.GetFilesT(type)
                db=self.inputsDB

            elem=None
            for file in fileList:
                elem = QListViewItem( self.listView_dataS )
                (date, ut_time, type, filter, texp, detector_id, run_id, ra, dec, object, mjd)=db.GetFileInfo(file)
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
            it=QListViewItemIterator (self.listView_dataS)
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

    def listView_popup_slot(self, listItem, mouse_pointer):
        
        #### Get items selected in the ListView
        self.m_popup_l_sel = []
        it = QListViewItemIterator (self.listView_dataS)
        listViewItem = it.current()
        father = None
        while listViewItem: 
            if listViewItem.isSelected():
                self.m_popup_l_sel.append(str(listViewItem.text(0)))
                # if a parent of the group, no popup to show 
                if (self.comboBox_classFilter.currentText()=="GROUP" and \
                    listViewItem.firstChild()!=None): # is a parent of a group
                    father = listViewItem
                    self.m_listView_first_item_selected = father
                    break
            it+=1
            listViewItem = it.current()
        
        if len(self.m_popup_l_sel)<=0:
            return
        #print "LIST=", self.m_popup_l_sel
        
        # we have selected a group father
        if father: # we have selected a group father
            #### Create the Group Popup menu
            popUpMenu = QPopupMenu()
            popUpMenu.insertItem("Reduce Obs. Sequece", self.reduceSequence_slot, 0, 1 )
            popUpMenu.insertItem("Copy files to clipboard", self.copy_files_slot, 0, 2 )
               
            group_files = []
            child = father.firstChild()
            while child:
                group_files.append(str(child.text(0)))
                child = child.nextSibling()
            #print "CHILDS=",group_files
        else:
            #### Create the Files Popup menu 
            popUpMenu = QPopupMenu()
            popUpMenu.insertItem("Display image", self.display_slot, 0, 1 )
            popUpMenu.insertItem("Copy files to clipboard", self.copy_sel_files_slot, 0, 20 )
            popUpMenu.insertSeparator()
            popUpMenu.insertItem("Create Master Dark",  self.createMasterDark_slot, 0, 2 )
            popUpMenu.insertItem("Create Master Dome-Flat", self.createMasterDFlat_slot, 0, 3)
            popUpMenu.insertItem("Create Master Twilight-Flat", self.createMasterTwFlat_slot, 0, 4 )
            popUpMenu.insertItem("Create Gain Map (SuperFlat)", self.createGainMap_slot, 0, 5 )
            popUpMenu.insertItem("Create Bad Pixel Mask", self.createBPM_slot, 0, 6 )
            popUpMenu.insertSeparator()
            popUpMenu.insertSeparator()
            popUpMenu.insertItem("Subtract (own) sky", self.subtract_ownSky_slot, 0, 7 )
            popUpMenu.insertItem("Subtract nearest sky", self.subtract_nearSky_slot, 0, 8 )
            popUpMenu.insertItem("Quick single Pre-Reduction", self.do_quick_reduction_slot, 0, 9 )
            popUpMenu.insertItem("Stack-Shift and Align", self.createStackedFrame_slot, 0, 10 )
            #popUpMenu.insertItem("Build super-Mosaic", self.createSuperMosaic_slot, 0, 11)
            popUpMenu.setItemEnabled(12, False)
            popUpMenu.insertSeparator()
            popUpMenu.insertSeparator()
            popUpMenu.insertItem("Rough Astrometry", self.do_raw_astrometry, 0, 13)
            popUpMenu.insertItem("Rough Photometry", self.do_raw_photometry, 0, 14)
            popUpMenu.insertSeparator()
            popUpMenu.insertSeparator()
            popUpMenu.insertItem("Show Stats", self.show_stats_slot, 0, 15 )
            popUpMenu.insertItem("FWHM estimation", self.fwhm_estimation_slot, 0, 16 )
            popUpMenu.insertItem("Background estimation", self.background_estimation_slot, 0, 17 )
            popUpMenu.insertItem("Imexam", self.imexam_slot, 0, 18 )
            popUpMenu.insertSeparator()
            # Math sub-menu
            subPopUpMenu = QPopupMenu()
            subPopUpMenu.insertItem("Substract Frames", self.subtractFrames_slot, 0, 1 )
            subPopUpMenu.insertItem("Sum Frames", self.sumFrames_slot, 0, 2 )
            subPopUpMenu.insertItem("Divide Frames", self.divideFrames_slot, 0, 3 )
            popUpMenu.insertItem("Math", subPopUpMenu, 0, 19 )
            # Fits sub-menu
            subPopUpMenu2 = QPopupMenu()
            subPopUpMenu2.insertItem("Split MEF file", self.splitMEF_slot, 0, 1 )
            subPopUpMenu2.insertItem("Join MEF file", self.joinMEF_slot, 0, 2 )
            #subPopUpMenu2.insertItem("Slice cube", self.sliceCube_slot, 0, 3 )
            #subPopUpMenu2.insertItem("Coadd cube", self.coaddCube_slot, 0, 4 )
            popUpMenu.insertItem("Fits", subPopUpMenu2, 0, 20 )
            
    
            ## Disable some menu items depeding of the number of item selected in the list view
            if len(self.m_popup_l_sel)==1:
                popUpMenu.setItemEnabled(2, False)
                popUpMenu.setItemEnabled(3, False)
                popUpMenu.setItemEnabled(4, False)
                popUpMenu.setItemEnabled(5, False)
                popUpMenu.setItemEnabled(6, False)
                #popUpMenu.setItemEnabled(10, False)
                subPopUpMenu.setItemEnabled(1,False)
                subPopUpMenu.setItemEnabled(2,False)
                #subPopUpMenu.setItemEnabled(3,False)
            elif len(self.m_popup_l_sel)>1:
                popUpMenu.setItemEnabled(1, False)
                popUpMenu.setItemEnabled(7, False)
                popUpMenu.setItemEnabled(8, False)
                #popUpMenu.setItemEnabled(11, False)
            elif len(self.m_popup_l_sel)>2:
                subPopUpMenu.setItemEnabled(1,False)
                subPopUpMenu.setItemEnabled(3,False)
            
            if len(self.m_popup_l_sel)<5:
                popUpMenu.setItemEnabled(9, False)
            
        ## Finally, execute the popup
        popUpMenu.exec_loop(QCursor.pos())   
    
    def reduceSequence_slot(self):
        """
        Run the data reduction of the current grouped-sequence selected in the
        ListView when the GROUP view is selected.
        """
        
        group_files = []
        child = self.m_listView_first_item_selected.firstChild()
        while child:
            group_files.append(str(child.text(0)))
            child = child.nextSibling()
        #print "CHILDS=",group_files
        
        #Change to working directory
        os.chdir(self.m_tempdir)
        #Change cursor
        self.setCursor(Qt.waitCursor)
        #Create working thread that compute sky-frame
        try:
            self._task = RS.ReductionSet (group_files, self.m_outputdir, 
                                        out_file=None,
                                        obs_mode="dither", dark=None, flat=None, 
                                        bpm=None, red_mode="quick", group_by=self.group_by, 
                                        check_data=True, config_dict=self.config_opts,
                                        external_db_files=self.outputsDB.GetFiles())
            
            thread = reduce.ExecTaskThread(self._task.reduceSet, self._task_info_list, "quick")
            thread.start()
        except Exception,e:
            #Anyway, restore cursor
            # Although it should be restored in checkLastTask, could happend an exception while creating the class RS,
            # thus the ExecTaskThread can't restore the cursor
            self.setCursor(Qt.arrowCursor) 
            QMessageBox.critical(self, "Error", "Error while group data reduction: \n%s"%str(e))
            raise e
    
    def copy_files_slot(self):
        """
        Copy the file list fullnames belonging to the current group to the 
        clipboard in order to allow copy&paste
        """
        
        text=""
        child = self.m_listView_first_item_selected.firstChild()
        while child:
            text=text+str(child.text(0))+"\n"
            child=child.nextSibling()    
        
        clipboard = QApplication.clipboard()
        clipboard.setText(text)

    def copy_sel_files_slot(self):
        """
        Copy the selected file list fullnames to the clipboard in order to 
        allow copy&paste
        """
        
        text=""
        for filename in self.m_popup_l_sel:
            text=text+str(filename)+"\n"
        
        clipboard = QApplication.clipboard()
        clipboard.setText(text)    
            
    def splitMEF_slot(self):
        """Split each MEF selected file from the list view into NEXT separate single FITS file, where NEXT is number of extensions.
           As result, NEXT files should be created
        """
        for file in self.m_popup_l_sel:
            try:
                mef = misc.mef.MEF([file])
                mef.doSplit(".%02d.fits")
            except Exception,e:
                log.debug("Cannot split file %s. Maybe it's not a MEF file",str(e))
                QMessageBox.critical(self, "Error", "Cannot split file : %s \n Maybe it's not a MEF file"%(file))
                #self.logConsole.info(QString(str(line)))
    
    def joinMEF_slot(self):
        """Join/stitch a MEF file to stitch them together
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
        """To know which item is selected """
                
        if listItem:
            if self.comboBox_classFilter.currentText()=="GROUP" and listItem.firstChild()!=None: # it'a a parent
                self.m_listView_item_selected=None # for consistency
                return
            else: self.m_listView_item_selected=str(listItem.text(0))
    
    def imexam_slot(self):
        """Imexam the currect filename selected """
        
        try:
            #self.m_processing = False    # Pause autochecking coming files - ANY MORE REQUIRED ?, now using a mutex in thread !!!!
            thread = reduce.ExecTaskThread(self.run_imexam, self._task_info_list, self.m_popup_l_sel[0] )
            thread.start()
        except Exception,e:
            QMessageBox.critical(self, "Error", "Error while Imexam with image %s"%(self.m_popup_l_sel[0]))
            raise e
        
    def run_imexam(self, filename):
        """ 
        Run Imexam the currect filename. First we start the DS9 display  
        
        DOES NOT WORK !!!
        
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
        #QMessageBox.critical(self,"Error", "Could not save the current document")
        popUpMenu = QPopupMenu()
        popUpMenu.insertItem("Create Master Dark")
        popUpMenu.insertItem("Create Master Dome-Flat")
        popUpMenu.insertItem("Create Master Sky-Flat")
        #actions["act1"].addTo(fileMenu)
        popUpMenu.exec_loop(QCursor.pos())
        return

    def iraf_console_slot(self):
        """Open a IRAF console session"""
        os.system("cd $HOME/iraf;/usr/local/bin/xgterm -title IRAF -cr red -ms blue -sb -sl 1000 -geometry 100x30 -bg grey -fg black -e cl &")
        
    def start_ds9_slot(self):
        """Start DS9 display, and if exixts, open the currect selected file in the list view"""
        
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
        #if utils.runCmd("/usr/local/bin/aladin &")==0: # some error
            #log.error("Error starting Aladin application")
            #raise Exception("Error executing command Aladin")
            #QMessageBox.warning(self,"Error","Error starting Alading application")
        #else:    
            #self.logConsole.info("Aladin launched !!!")

    def pushB_sel_masteDark_slot(self):
        source=QFileDialog.getOpenFileName( self.m_default_data_dir, "Master file (*.fit*)",  self, "Source Dialog","select master")

        if str(source):
            self.lineEdit_masterDark.setText(source)
            self.m_masterDark=str(source)

    def pushB_sel_masterFlat_slot(self):
        """
        @summary: to be done
        """
        source=QFileDialog.getOpenFileName( self.m_default_data_dir, "Master file (*.fit*)",  self, "Source Dialog","select master")

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
        source=QFileDialog.getOpenFileName( self.m_default_data_dir, "Master file (*.fit*)",  self, "Source Dialog","select master")

        if str(source):
            self.lineEdit_masterMask.setText(source)
            self.m_masterMask=str(source)                       
    
    def genFileList( self, file_list, outFilename ):
      
      """ Generate a file 'outFilename' listing the files passed as a python-list in the file_list parameter"""
      fo = open(outFilename,"w")
      for file in file_list:
        fo.write(file+"\n")
      fo.close()
                             
    ############################################################################
    ############# PROCESSING STAFF #############################################
    ############################################################################
        
    def subtractFrames_slot(self):
        """This method is called to subtract two images selected from the File List View"""

        if (len(self.m_popup_l_sel)!=2):
            QMessageBox.critical(self, "Error", "You need to select  2 files")
        else:
            outFilename = QFileDialog.getSaveFileName(self.m_outputdir+"/sub.fits", "*.fits", self, "Save File dialog")
            if not outFilename.isEmpty():
                try:
                    #Change cursor
                    self.setCursor(Qt.waitCursor)
                    self.m_processing = False    # Pause autochecking coming files - ANY MORE REQUIRED ?, now using a mutex in thread !!!!
                    thread = reduce.ExecTaskThread(self.mathOp, 
                                                   self._task_info_list, 
                                                   self.m_popup_l_sel,'-', 
                                                   str(outFilename))
                    thread.start()
                except:
                    QMessageBox.critical(self, "Error", "Error while subtracting files")
                    raise

    def mathOp(self, files, operator, outputFile=None):
        """This method will do the math operation (+,-,/) specified with the input files"""
        
        log.debug("Start mathOp")
        
        if outputFile==None:
            output_fd, outputFile = tempfile.mkstemp(suffix='.fits', dir=self.m_outputdir)
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
                misc.utils.listToFile(files, self.m_outputdir+"/files.tmp") 
                iraf.mscred.combine(input=("@"+(self.m_outputdir+"/files.tmp").replace('//','/')),
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
        
    def sumFrames_slot(self):
        """This methot is called to sum two images selected from the File List View"""

        if (len(self.m_popup_l_sel)<2):
            QMessageBox.critical(self, "Error", "You need to select  at least 2 files")
        else:
            outFilename = QFileDialog.getSaveFileName(self.m_outputdir+"/sum.fits", "*.fits", self, "Save File dialog")
            if not outFilename.isEmpty():
                try:
                    #Change cursor
                    self.setCursor(Qt.waitCursor)
                    self.m_processing = False    # Pause autochecking coming files - ANY MORE REQUIRED ?, now using a mutex in thread !!!!
                    thread = reduce.ExecTaskThread(self.mathOp, self._task_info_list, self.m_popup_l_sel,'+', str(outFilename))
                    thread.start()
                except:
                    QMessageBox.critical(self, "Error", "Error while adding files")
                    raise
                

    def divideFrames_slot(self):
        """This methot is called to divide two images selected from the File List View"""

        if (len(self.m_popup_l_sel)!=2):
            QMessageBox.critical(self, "Error", "You need to select  2 files")
        else:
            outFilename = QFileDialog.getSaveFileName(self.m_outputdir+"/div.fits", "*.fits", self, "Save File dialog")
            if not outFilename.isEmpty():
                try:
                    #Change cursor
                    self.setCursor(Qt.waitCursor)
                    self.m_processing = False    # Pause autochecking coming files - ANY MORE REQUIRED ?, now using a mutex in thread !!!!
                    thread=reduce.ExecTaskThread(self.mathOp, self._task_info_list, self.m_popup_l_sel,'/', str(outFilename))
                    thread.start()
                except:
                    QMessageBox.critical(self, "Error", "Error while dividing files")
                    raise

    def createMasterDark_slot(self):

        if len(self.m_popup_l_sel)<=1:
            QMessageBox.information(self,"Info","Not enought frames !")
            return

        outfileName = QFileDialog.getSaveFileName(self.m_outputdir+"/master_dark.fits", "*.fits", self, "Save File dialog")
        if not outfileName.isEmpty():
            try:
                texp_scale = False
                self.setCursor(Qt.waitCursor)
                self._task = reduce.calDark.MasterDark(self.m_popup_l_sel, self.m_tempdir, str(outfileName), texp_scale)
                thread = reduce.ExecTaskThread(self._task.createMaster, self._task_info_list)
                thread.start()
            except Exception, e:
                self.setCursor(Qt.arrowCursor)
                QMessageBox.critical(self, "Error", "Error while creating master Dark. \n"+str(e))
                raise e
        else:
            pass

    
    def do_quick_reduction_slot(self):
        """ 
        Run a quick reduction of the current user selected science files 
        in the list view panel.
        """
            
        if len(self.m_popup_l_sel)<5:
            QMessageBox.information(self,"Info","Error, not enought frames selected to reduce (>=5) !")
            return
              
        #Change to working directory
        os.chdir(self.m_tempdir)
        #Change cursor
        self.setCursor(Qt.waitCursor)
        #Create working thread that compute sky-frame
        try:
            self._task = RS.ReductionSet( self.m_popup_l_sel, self.m_outputdir, 
                                          out_file=None,
                                          obs_mode="dither", dark=None, 
                                          flat=None, bpm=None, red_mode="quick",
                                          group_by=self.group_by, check_data=True, 
                                          config_dict=self.config_opts)
                                            
            thread = reduce.ExecTaskThread(self._task.reduceSet, self._task_info_list, "quick")
            thread.start()
        except Exception,e:
            #Anyway, restore cursor
            # Although it should be restored in checkLastTask, could happend an exception while creating the class RS,
            # thus the ExecTaskThread can't restore the cursor
            self.setCursor(Qt.arrowCursor) 
            QMessageBox.critical(self, "Error", "Error while Quick data reduction: \n%s"%str(e))
            raise e

    def createMasterDFlat_slot(self):
        
        """Create a master dome flat field from selected frames"""
        
        if len(self.m_popup_l_sel)>2:
            outfileName = QFileDialog.getSaveFileName(self.m_outputdir+"/master_dflat.fits", "*.fits", self, "Save File dialog")
            if not outfileName.isEmpty():
                try:
                    self.setCursor(Qt.waitCursor)
                    self._task = reduce.calDomeFlat.MasterDomeFlat (self.m_popup_l_sel, self.m_tempdir , str(outfileName))
                    thread=reduce.ExecTaskThread(self._task.createMaster, self._task_info_list)
                    thread.start()
                except:
                    self.setCursor(Qt.arrowCursor)
                    QMessageBox.critical(self, "Error", "Error while creating master Dome Flat")
                    raise
        else:
            QMessageBox.information(self,"Info","Error, not enought frames (>2) !")

    def createMasterTwFlat_slot(self):
      
        """Create a master twilight flat field from selected frames"""
        
        if len(self.m_popup_l_sel)>2:
            outfileName = QFileDialog.getSaveFileName(self.m_outputdir+"/master_twflat.fits", 
                                                      "*.fits", 
                                                      self, 
                                                      "Save File dialog")
            if not outfileName.isEmpty():
                try:
                    self.setCursor(Qt.waitCursor)
                    self._task = reduce.calTwFlat.MasterTwilightFlat (self.m_popup_l_sel, 
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
            QMessageBox.information(self,"Info","Not enought frames !")
            return

        outfileName = QFileDialog.getSaveFileName (self.m_outputdir+"/gainmap.fits", 
                                                   "*.fits", self, 
                                                   "Save File dialog")
        if not outfileName.isEmpty():
            flat_type = self.inputsDB.GetFileInfo(self.m_popup_l_sel[0])[2]
            if flat_type.count("SCIENCE"):
                try:
                    self.setCursor(Qt.waitCursor)
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
                    self.setCursor(Qt.waitCursor)
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
                    self.setCursor(Qt.waitCursor)
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
            #out_file = "/tmp/skysub.fits"
        
            #cmd="sex %s -c %s -FITS_UNSIGNED Y -DETECT_MINAREA %s  -DETECT_THRESH %s  -CHECKIMAGE_TYPE -BACKGROUND -CHECKIMAGE_NAME %s" % (input_file, sex_config, str(minarea), str(threshold), out_file)  
            cmd="sex %s -c %s -DETECT_MINAREA %s -DETECT_THRESH %s -CHECKIMAGE_TYPE -BACKGROUND -CHECKIMAGE_NAME %s -PARAMETERS_NAME %s -FILTER_NAME %s"% (input_file, sex_config, str(minarea), str(threshold), out_file, sex_param, sex_conv )  
            
            """
            ***
                Be sure that SExtractor working directory has R/W in order to create
                temporal files created by SExtractor !!!
            ***
            """
            #Change to working directory
            os.chdir(self.m_tempdir)
            #Change cursor
            self.setCursor(Qt.waitCursor) # restored checkLastTask
            #Call external script (papi)
            self.m_processing = True
            self._proc=RunQtProcess(cmd, self.logConsole, self._task_info_list, out_file)      
            self._proc.startCommand()

            """
            # Next implementation does not work because sex.run() doesn't retuen
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
                self.setCursor(Qt.waitCursor)
                updateconfig=True
                clean=False
                thread=reduce.ExecTaskThread(sex.run, self._task_info_list, \
                                             # args for Sextractor:run() method
                                             input_file, updateconfig, clean)
                thread.start()    
            except Exception,e:
                self.setCursor(Qt.arrowCursor)
                QMessageBox.critical(self, "Error", "Error while sky subtraction "+str(e))
                raise e
            """
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
                        QMessageBox.information(self, "Info", "Not enought science frames found")  
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
                    if (file==self.m_listView_item_selected): file_n=p
                    else: p+=1
                
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
                self.logConsole.info("Not enought files for sky subtraction (>4)")
                log.debug("Not enought files for sky subtraction (>5)")
                return

        # ##################################    
        # Now, compute the sky-subtraction
        # ##################################
                            
        #Change to working directory
        os.chdir(self.m_tempdir)
        #Change cursor
        self.setCursor(Qt.waitCursor)
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
            
            thread = reduce.ExecTaskThread(self._task.subtractNearSky, self._task_info_list, self._task.rs_filelist, file_n)
            thread.start()
        except:
            #Anyway, restore cursor
            # Although it should be restored in checkLastTask, could happend an exception while creating the class RS,
            # thus the ExecTaskThread can't restore the cursor
            self.setCursor(Qt.arrowCursor) 
            QMessageBox.critical(self, "Error", "Error while subtracting near sky")
            raise 
              
    def createBPM_slot(self):
        """ Create a Bad Pixel Mask from a set of selected files (flats)
        """
        if len(self.m_popup_l_sel)>3:
            outfileName = QFileDialog.getSaveFileName("/tmp/BadPixelMask", "*.pl", self, "Save File dialog")
            if not outfileName.isEmpty():
                show = True
                try:
                    self.genFileList(self.m_popup_l_sel, "/tmp/bpm.list")
                    dark = "/tmp/master_dark.fits"
                    lsig = 10
                    hsig = 10
                    self.setCursor(Qt.waitCursor)
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
                    self.setCursor(Qt.arrowCursor)
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
            Compute a stacked frame (subtract sky, shift and align) from a set of nearest (ra,dec, mjd) 
            frames, selected by user or automatic search.
            Actually, it is basically a quick-reduction !!! (see do_quick_reduction_slot), with the only difference
            that here we can look for the files to be reduced in automatic way or use the ones selected  by the user (as in do_quick_reduction_slot)
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
                    QMessageBox.information(self, "Info", "Not enought science frames found")  
                    return
            else:
                QMessageBox.information(self, "Info", "Selected frame is not a science frame") 
                return
            
            # Create nearest file list (ar,dec,mjd) from current selected science file
            view_list = ""
            i=1
            for file in near_list:
                if file==self.m_listView_item_selected: file_n=i
                else: i=i+1
                if (datahandler.ClFits(file).getType()!='SCIENCE'):
                    QMessageBox.critical(self, "Error", QString("File %1 is not a science frame").arg(file))
                else:
                    file_list.append(file)
                    view_list +=  file + "\n"
            
            resp=QMessageBox.information(self, "Info", QString("Selected near frames are:\n %1").arg(view_list), \
                                            QMessageBox.Ok, QMessageBox.Cancel)
            if resp==QMessageBox.Cancel:
                return
                    
        # CASE 2: Stack frames selected by user in the list_view
        elif len(self.m_popup_l_sel)>1 and len(self.m_popup_l_sel)<4:
            QMessageBox.information(self,"Info","Error, not enought frames selected to reduce (>4) !")
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
        outfileName = QFileDialog.getSaveFileName(self.m_outputdir+"/red_result.fits", "*.fits", self, "Save File dialog")
        if outfileName.isEmpty(): return # nothig to do !
           
        #Change to working directory
        os.chdir(self.m_tempdir)
        #Change cursor
        self.setCursor(Qt.waitCursor)
        #Create working thread that compute sky-frame
        if len(file_list)>1:
            try:
                self._task = RS.ReductionSet( file_list, self.m_outputdir, out_file=self.m_outputdir+"/red_result.fits", \
                                            obs_mode="dither", dark=None, flat=None, bpm=None, red_mode="quick", \
                                            group_by=self.group_by, check_data=True, config_dict=self.config_opts)
                
                thread = reduce.ExecTaskThread(self._task.reduceSet, self._task_info_list, "quick")
                thread.start()
            except:
                #Anyway, restore cursor
                # Although it should be restored in checkLastTask, could happend an exception while creating the class RS,
                # thus the ExecTaskThread can't restore the cursor
                self.setCursor(Qt.arrowCursor) 
                QMessageBox.critical(self, "Error", "Error while building stack")
                raise
        
        
    def createSuperMosaic_slot(self):
        """
        TBD : here is only a scheme of the routine !!!!
        Create an stitched wide-mosaic of next frames
        """
        QMessageBox.information(self, "Info", "Sorry, funtion not yet implemented !")
        
        # Next code is only a raw template 
        # Check list lenght
        if len(self.m_popup_l_sel)<1:
            QMessageBox.information(self, "Info", "Not enought science pre-reduced frames selected")
            return
        
        if len(self.m_popup_l_sel)==1:
            QMessageBox.information(self, "Info", "Only one file was selected, automatic file grouping will be done.")
            return
            # NEED TO BE DONE
        # Stack frame user selected  
        elif len(self.m_popup_l_sel)>1:
            # Create file list from current selected science files
            listfile=self.m_tempdir+"/mosaic_files.list"      
            file_lst= open( listfile, "w" )
            for file in self.m_popup_l_sel :
                if (datahandler.ClFits(file).getType()!='SCIENCE_REDUCED'):
                    QMessageBox.critical(self, "Error", QString("File %1 is not a science reduced frame").arg(file))
                    file_lst.close()
                    return
                else:
                    file_lst.write( file +"\n")
            file_lst.close()
        
        # Call external script (SWARP)
        cwd=os.getcwd()
        os.chdir(self.m_tempdir)
        cmd=self.m_tempdir+"/build_mosaic %s" %(listfile)
        #Change to working directory
        os.chdir(self.m_tempdir)
        #Change cursor
        self.setCursor(Qt.waitCursor) # restored in checkLastTask
        # Call external script (papi)
        self._proc = RunQtProcess(cmd, self.logConsole, self._task_info_list, self.m_outputdir+"/mosaic.fits" )      
        self._proc.startCommand()
      
      
    def do_raw_astrometry(self):
        """
        @summary: Compute an astrometric solution for the selected file in the 
        main list view panel.
        Basically, it is SCAMP is called and compared to 2MASS catalog 
        (or GSC2.3, USNO-B) to compute astrometric calibration. 
        A new WCS is added to the header.
        
        @return: None, but produce an astrometrically calibrated image.
          
        @todo: we should check the image is pre-reduced or at least, with 
        enought objects 
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
                self.setCursor(Qt.waitCursor)
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
        enought objects 
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
                self.setCursor(Qt.waitCursor)
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
            self.setCursor(Qt.waitCursor)
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
                self.setCursor(Qt.arrowCursor) 
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
                                source directory (If no, only the new ones coming \
                                will be processed) ?"), 
                                QMessageBox.Ok, QMessageBox.Cancel)
            if res==QMessageBox.Cancel:
                return
            else:
                self.processFiles()
            
 
        """
        it=QListViewItemIterator (self.listView_OS)
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
        
    def processFiles(self, files=None):
        """
        @summary: Process the files provided; if any files were given, all the files 
        in the current Source List View (but not the output files) will be
        processed.
        
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
        self.logConsole.info("++ Starting to process next files :")
        for file in files:
            self.logConsole.info("     - " + file)
            #self.logConsole.info(QString("    - %1").arg(file))
        self.logConsole.info("... processing sequence ...")
            
        #Change to working directory
        os.chdir(self.m_tempdir)
        #Change cursor
        self.setCursor(Qt.waitCursor)
        self.m_processing = False    # Pause autochecking for coming files - ANY MORE REQUIRED ?, now using a mutex in thread !!!!
        #Create working thread that process the files
        try:
            # generate a random filename for the master, to ensure we do not overwrite any file
            output_fd, outfilename = tempfile.mkstemp(suffix='.fits', prefix='redObj_', dir=self.m_outputdir)
            os.close(output_fd)
            os.unlink(outfilename) # we only need the name
            self._task = RS.ReductionSet( files, self.m_outputdir, out_file=outfilename,
                                            obs_mode="dither", dark=None, 
                                            flat=None, bpm=None, red_mode="quick",
                                            group_by=self.group_by, check_data=True,
                                            config_dict=self.config_opts,
                                            external_db_files=self.outputsDB.GetFiles())
            # provide the outputDB files as the external calibration files for the RS 
            
            log.debug("ReductionSet created !")
            thread = reduce.ExecTaskThread(self._task.reduceSet, self._task_info_list)
            thread.start()
        except Exception,e:
            #Anyway, restore cursor
            # Although it should be restored in checkLastTask, could happend an exception while creating the class RS,
            # thus the ExecTaskThread can't restore the cursor
            self.setCursor(Qt.arrowCursor) 
            QMessageBox.critical(self, "Error", "Error while processing Obs. Sequence: \n%s"%str(e))
            self.m_processing = False
            raise e
        
      
################################################################################


class LoggingConsole (object):
    
    def __init__ (self, textEdit1=None, textEdit2=None, *a, **k):
        
        super (LoggingConsole, self).__init__ (*a, **k)

        self.textEdit_w1 = textEdit1
        self.textEdit_w2 = textEdit2
        
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
        

    def append(self, message='', tag=None):
        
        if self.textEdit_w1 is not None: 
            self.textEdit_w1.append(self.logMsg(message,tag))
        if self.textEdit_w2 is not None:
            self.textEdit_w2.append(self.logMsg(message,tag))
    
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
            prefix = "<info_tag>"
            suffix = "</info_tag>"
        elif tag=="WARNING":
            prefix = "<warning_tag>"
            suffix = "</warning_tag>"
        elif tag=="ERROR":
            prefix = "<error_tag>"
            suffix = "</error_tag>"
        elif tag=="DEBUG":
            prefix = "<debug_tag>"
            suffix = "</debug_tag>"
        else:
            prefix = ""
            suffix = ""    
            
            
        #s_time = time.strftime("[%Y%m%d-%H:%M:%S] ", time.gmtime())
        s_time = "["+datetime.datetime.utcnow().isoformat()+"] "
        prefix += s_time 
        
        return prefix + msg + suffix
