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
# runGUI (run main GUI for PANIC pipeline)
#
# runGUI.py
#
# Last update 22/Sep/2009
#
################################################################################


    ####################################################################
    #                                                                  #
    #                   This is where it all starts.                   #
    #         The body of this program is at the end of this file      #
    #                                                                  #
    ####################################################################


from qt import *



from panicQL import *
import sys
import os
import subprocess
import threading
import os.path
import fnmatch
import shutil
import time


import reduce
import reduce.calTwFlat
import reduce.calBPM_2
import reduce.checkQuality
import misc.fileUtils
import misc.utils as utils
import datahandler
import misc.display as display
from runQtProcess import *


# Interact with FITS files
import pyfits

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
  
  
class MainGUI(panicQL):

            
    def __init__(self):               
        
        log.debug("Here start Quick-Look tool !!!")
        panicQL.__init__(self)
        ## Init member variables
        self.m_default_data_dir = os.environ['PANIC_DATA']
        self.m_frameList_dark = ''
        self.m_frameList_dflat = ''
        self.m_frameList_sflat = ''
        self.m_masterDark = '/tmp/master_dark.fits'
        self.m_masterFlat = '/tmp/master_flat.fits'
        self.m_masterMask = '/tmp/master_mask.fits'


        self.m_sourcedir = ''
        self.m_outputdir = '/tmp/out/'
        self.m_tempdir = '/tmp/'
        #self.m_panic_dir = '/disk-a/caha/panic/DEVELOP/PIPELINE/PANIC/'       # PANIC SW directory
        self.m_papi_dir  = os.environ['PAPI_HOME']  # PAnic PIpeline directory
        self._ini_cwd = os.getcwd()
        self._last_cmd=None             # Last external shell command executed or started
        self._last_cmd_output=None      # Output (file) that last shell command should generate

        self.m_listView_first_item_selected=''
        self.m_listView_item_selected='' 
        self.m_show_imgs = False
        self.m_proc_imgs = False
        self.m_processing = False
        self._proc = None # variable to handle QProcess tasks
        
        ## Create log tags
        # Error
        item =QStyleSheetItem( self.textEdit_log.styleSheet(), "error_tag" )
        item.setColor( QColor("red") )
        item.setFontWeight( QFont.Bold )
        item.setFontUnderline( True )
        # Warning
        item =QStyleSheetItem( self.textEdit_log.styleSheet(), "warning_tag" )
        item.setColor( QColor("blue") )
        item.setFontWeight( QFont.Bold )
        item.setFontUnderline( True )
        # Info
        item =QStyleSheetItem( self.textEdit_log.styleSheet(), "info_tag" )
        item.setColor( QColor("black") )
        item.setFontWeight( QFont.Bold )
        item.setFontUnderline( True )
        
        
        self.textEdit_log.append("Wellcome to the <info_tag> PANIC QuickLook tool (v1.0) ! </info_tag>")
        

        self.initialize()
        
        ## Init in memory Database
        datahandler.dataset.initDB()
        
        ## Data Collectors init
        self.file_pattern = str(self.lineEdit_filename_filter.text())
        self.dc=datahandler.DataCollector("dir", self.m_sourcedir, self.file_pattern , self.new_file_func)
        self.dc_outdir=None
        #datahandler.DataCollector("dir", self.m_outputdir, self.file_pattern, self.new_file_func_out)
        
        # Task management
        self._task = None                       # Pointer to task/thread that is running 
        #self._task_event = threading.Event()    # Event to syncronize ExecTaskThread and ExecTaskThread
        self._task_info=None
        self._task_info_list = []               # Queue-list where task status is saved
        
        self._task_timer = QTimer( self )
        self.connect( self._task_timer, SIGNAL("timeout()"), self.checkLastTask )
        self._task_timer.start( 1000, False )    # 1 second continuous timer
        
        ##Start display (DS9)
        #display.startDisplay()
        #time.sleep(1) # wait until display is up
        
    def initialize(self):
        """This method will initialize some values in the GUI and in the members variables"""

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
    
    def new_file_func(self, filename, process=True):
        """ Function executed when a new file is detected into the data source dir or into the out_dir"""
        
        log.debug( "New file detected--> %s", filename)
        self.textEdit_log.append("New file detected in source-->  " + filename)
        ## Insert into DB
        #datahandler.dataset.initDB()
        datahandler.dataset.filesDB.insert(filename)
        ## Query DB
        (date, ut_time, type, filter, texp, detector_id, run_id, ra, dec, object)=datahandler.dataset.filesDB.GetFileInfo(filename)
        #fileinfo=datahandler.dataset.filesDB.GetFileInfo(str(dir)+"/"+filename)c
        #print "FILEINFO= ", fileinfo
        ## Show into ListView table
        elem = QListViewItem( self.listView_dataS )
        elem.setText (0, str(filename))
        elem.setText (1, str(type))
        elem.setText (2, str(filter))
        elem.setText (3, str(texp))
        elem.setText (4, str(date)+"::"+str(ut_time))
        elem.setText (5, str(object))
        elem.setText (6, str(ra))
        elem.setText (7, str(dec))
        
        ## Last frame
        self.lineEdit_last_file.setText(str(os.path.basename(filename)))
        
        if process:
            ###############
            #QL-Mode
            ###############
            #QL-Mode: None
            if self.comboBox_QL_Mode.currentText()=="None":
                return
            #QL-Mode: Lazy
            elif self.comboBox_QL_Mode.currentText().contains("Lazy"):
                if self.m_show_imgs:
                    display.showFrame(filename)
                ## Process
                if self.m_proc_imgs:
                    self.m_processing = True    # Pause autochecking files
                    self.process(filename)
                    self.m_processing = False   # Resume autochecking files
            
            #QL-Mode: Pre-reduction
            elif self.comboBox_QL_Mode.currentText().contains("Pre-reduction"):
                return ## TO BE DONE
            #QL-Mode: UserDef_1
            elif self.comboBox_QL_Mode.currentText().contains("UserDef_1"):
                return ## TO BE DONE
            #QL-Mode: UserDef_2
            elif self.comboBox_QL_Mode.currentText().contains("UserDef_2"):
                return ## TO BE DONE
            #QL-Mode: UserDef_3
            elif self.comboBox_QL_Mode.currentText().contains("UserDef_3"):
                return ## TO BE DONE
        
    def new_file_func_out(self, filename):
        """Callback used when a new file is detected in output dir"""
                     
        self.new_file_func(filename, process=False)
        
    def process(self, filename):
        ## Process the new image received with the QL pipeliene recipes
        print "------------------------------------------>m_processing(process)=", self.m_processing
        self.textEdit_log.append("++Processing file : " + filename)
        self.QL2(filename, self.textEdit_log)


    #####################################################
    ### SLOTS ###########################################
    #####################################################

    def setDataSourceDir_slot(self):
        #dir=""
        source=QFileDialog.getExistingDirectory( self.m_default_data_dir, self,"get existing directory", "Choose a directory",True )
        
        #source=QFileDialog.getOpenFileNames( "Source log (*.log)", self.m_default_data_dir, self, "Source Dialog","select source")
        ## NOTE: 'dir' can be a file or a directory
        
        if (not source):
            return
        else:
            #dir=str(source[0])
            dir=str(source)
            self.textEdit_log.append("+Source : " + dir)
            if (self.m_sourcedir != dir):
                self.lineEdit_sourceD.setText(dir)
                self.m_sourcedir = str(dir)
                ##Create DataCollector for a path     
                self.file_pattern = str(self.lineEdit_filename_filter.text())
                #print "FIL_PATTERN = ", fil_pat
                if os.path.isfile(dir):
                    self.dc = datahandler.DataCollector("geirs-file", str(dir), self.file_pattern , self.new_file_func)  
                elif os.path.isdir(dir):
                    self.dc = datahandler.DataCollector("dir", str(dir), self.file_pattern , self.new_file_func)
                
                ## Activate the autochecking of new files
                self.checkBox_autocheck.setChecked(True)
                ##Create QTimer for the data collector
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
            datahandler.dataset.filesDB.insert(str(dir)+"/"+filename)
            (date, ut_time, type, filter, texp, detector_id, run_id)=datahandler.dataset.filesDB.GetFileInfo(str(dir)+"/"+filename)
            #fileinfo=datahandler.dataset.filesDB.GetFileInfo(str(dir)+"/"+filename)
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
        datahandler.dataset.filesDB.ListDataSet()

    def autocheck_slot(self):
        """Called when check-button for Input dir is clicked"""
        
        ## Activate or deactivate the autochecking of new files
        if self.checkBox_autocheck.isChecked():
            ##Create QTimer for the data collector
            self.timer_dc = QTimer( self )
            self.connect( self.timer_dc, SIGNAL("timeout()"), self.checkFunc )
            self.timer_dc.start( 1500, False ) ## 1 seconds single-shot timer
        else:
            self.checkBox_outDir_autocheck.setChecked(False)
            #Stop the DataCollector timer
            self.timer_dc.stop()
            print "stoped!!!!"
            
    def checkOutDir_slot(self):
        """Called when check-button for Output dir is clicked"""
        
        ## Activate or deactivate the autochecking of new files
        if self.checkBox_outDir_autocheck.isChecked():
            if self.dc_outdir==None:
                self.dc_outdir=datahandler.DataCollector("dir", self.m_outputdir, self.file_pattern, self.new_file_func_out)
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
            try:
                self._task_info=self._task_info_list.pop()
                if self._task_info._exit_status == 0: # EXIT_SUCCESS, all was OK
                    if self._task_info._return!=None:
                        QMessageBox.information(self,"Info", QString("File %1 created").arg(self._task_info._return))
                        display.showFrame(self._task_info._return)
                else:
                    QMessageBox.critical(self, "Error", "Error while running task.  "+str(self._task_info._exc))
                #Anyway, restore cursor
                self.setCursor(Qt.arrowCursor)
                # Return to the previus working directory
                os.chdir(self._ini_cwd)
            except:
                raise Exception("Error while checking _task_info_list")
     
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
            
    def process_slot(self):
        if self.checkBox_process.isChecked():
            self.m_proc_imgs=True
        else:
            self.m_proc_imgs=False
                    
    def show_images_slot(self):
        if self.checkBox_show_imgs.isChecked():
            self.m_show_imgs=True
        else:
            self.m_show_imgs=False
            
    def data_grouping_slot(self):
      
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
        
        if dir and self.m_outputdir!=str(dir):
            self.lineEdit_outputD.setText(dir)
            self.m_outputdir=str(dir)
            self.textEdit_log.append("+Output dir : " + self.m_outputdir)
            
            ##Create DataCollector for a path     
            self.file_pattern = str(self.lineEdit_filename_filter.text())
            if os.path.isdir(self.m_outputdir):
                self.dc_outdir = datahandler.DataCollector("dir", self.m_outputdir, self.file_pattern , self.new_file_func_out)
            
            ## Activate the autochecking of new files
            self.checkBox_outDir_autocheck.setChecked(True)
                    
        else:
            #The same dir, nothing to do
            pass
        
    def setTempDir_slot(self):
        dir=QFileDialog.getExistingDirectory( self.m_default_data_dir, self,
                                             "get existing directory", "Choose a directory",True )
        if dir:
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
        
        filenames=QFileDialog.getOpenFileNames( "FITS files (*.fit*)", self.m_default_data_dir, self, "FileDialog","select files")

        if not filenames.isEmpty():
            for file in filenames:
                self.listBox_domeF.insertItem(str(file))  

    def clear_mainlist_slot(self):
      
        """ Remove all files from ListView, DataCollector and DB """
         
        self.listView_dataS.clear()
        self.dc.Clear()
        datahandler.dataset.filesDB.clearDB()


    def add_slot(self):
      
        """Add a new file to the main list panel, but not to the DataCollector list (dc), so it might be already inside """     
        filenames=QFileDialog.getOpenFileNames( "FITS files (*.fit*)", self.m_default_data_dir, self, "FileDialog","select files")

        if not filenames.isEmpty():
            for file in filenames:
                self.new_file_func(str(file))
        
    def del_slot(self):
        
        """ Delete the current selected file from the main list view panel, but we do not remote from the DataCollector neither file system"""
        
        self.m_popup_l_sel = []
        it=QListViewItemIterator (self.listView_dataS)
        listViewItem = it.current()
        while listViewItem: 
            if listViewItem.isSelected():
                fileName=str(listViewItem.text(0))
                self.listView_dataS.takeItem(listViewItem)
                #self.dc.Clear(fileName)
                datahandler.dataset.filesDB.delete( fileName )
                ##self.m_popup_l_sel.append(str(listViewItem.text(0)))
                print "SELECTED TO REMOVE=", fileName
            else:
                it+=1
            listViewItem = it.current()
          
        #datahandler.dataset.filesDB.ListDataSetNames()
    
    def display_slot(self):
      
        if self.m_listView_item_selected:
            display.showFrame(self.m_listView_item_selected)

    def filename_filter_slot(self):
        """ Modify filename filter for the data collector"""

        new_filter, ok = QInputDialog.getText("New filter","Enter filename filter:")
        if ok and not new_filter.isEmpty():
            self.lineEdit_filename_filter.setText( str(new_filter) )
            self.dc.SetFileFilter( str(new_filter) )

    def slot_classFilter(self, filter_string):
        """ Filter files on main ListView"""
        
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
                       
#########################################################################
###### Pop-Up ###########################################################
#########################################################################

    def listView_popup_slot(self, listItem, mouse_pointer):
        
        #### Get items selected in the ListView
        self.m_popup_l_sel = []
        it=QListViewItemIterator (self.listView_dataS)
        listViewItem = it.current()
        while listViewItem: 
            if listViewItem.isSelected():
                self.m_popup_l_sel.append(str(listViewItem.text(0)))
            it+=1
            listViewItem = it.current()
        
        if len(self.m_popup_l_sel)<=0:
            return
        
        #print "LIST=", self.m_popup_l_sel
        
        # Get the first one selected
        self.m_listView_first_item_selected = self.m_popup_l_sel[0]
        #### Create the Popup menu 
        popUpMenu = QPopupMenu()
        id=popUpMenu.insertItem("Display image", self.display_slot, 0, 1 )
        popUpMenu.insertSeparator()
        popUpMenu.insertItem("Create Master Dark",  self.createMasterDark_slot, 0, 2 )
        popUpMenu.insertItem("Create Master Dome-Flat", self.createMasterDFlat_slot, 0, 3)
        popUpMenu.insertItem("Create Master Twilight-Flat", self.createMasterTwFlat_slot, 0, 4 )
        popUpMenu.insertItem("Create SuperSky-Flat", self.createSkyFlat_slot_2, 0, 5 )
        popUpMenu.insertItem("Create Bad Pixel Mask", self.createBPM_slot, 0, 6 )
        popUpMenu.insertSeparator()
        popUpMenu.insertSeparator()
        popUpMenu.insertItem("Subtract (own) sky", self.subtract_sky_slot, 0, 7 )
        popUpMenu.insertItem("Subtract nearest sky", self.subtract_nearSky_slot, 0, 8 )
        popUpMenu.insertItem("Quick single Pre-Reduction", self.do_quick_reduction_slot, 0, 9 )
        popUpMenu.insertItem("Stack-Shift and Align", self.createStackedFrame_slot, 0, 10 )
        popUpMenu.insertItem("Build super-Mosaic", self.createSuperMosaic_slot, 0, 11)
        popUpMenu.setItemEnabled(12, False)
        popUpMenu.insertSeparator()
        popUpMenu.insertItem("Raw astrometry", self.do_raw_astrometry, 0, 13)
        popUpMenu.insertSeparator()
        popUpMenu.insertItem("Show Stats", self.show_stats_slot, 0, 14 )
        popUpMenu.insertItem("FWHM estimation", self.fwhm_estimation_slot, 0, 15 )
        popUpMenu.insertItem("Background estimation", self.background_estimation_slot, 0, 16 )
        popUpMenu.insertItem("Test", self.testSlot, 0, 17 )
        popUpMenu.insertSeparator()
        subPopUpMenu = QPopupMenu()
        subPopUpMenu.insertItem("Substract Frames", self.subtractFrames_slot, 0, 1 )
        subPopUpMenu.insertItem("Sum Frames", self.sumFrames_slot, 0, 2 )
        subPopUpMenu.insertItem("Divide Frames", self.divideFrames_slot, 0, 3 )
        popUpMenu.insertItem("Math", subPopUpMenu, 0, 20 )
        

        ## Disable some menu items depeding of the number of item selected in the list view
        if len(self.m_popup_l_sel)==1:
            popUpMenu.setItemEnabled(2, False)
            popUpMenu.setItemEnabled(3, False)
            popUpMenu.setItemEnabled(4, False)
            popUpMenu.setItemEnabled(5, False)
            popUpMenu.setItemEnabled(6, False)
            popUpMenu.setItemEnabled(10, False)
            subPopUpMenu.setItemEnabled(1,False)
            subPopUpMenu.setItemEnabled(2,False)
            subPopUpMenu.setItemEnabled(3,False)
        elif len(self.m_popup_l_sel)>1:
            popUpMenu.setItemEnabled(1, False)
            popUpMenu.setItemEnabled(8, False)
            popUpMenu.setItemEnabled(11, False)
        elif len(self.m_popup_l_sel)>2:
            subPopUpMenu.setItemEnabled(1,False)
            subPopUpMenu.setItemEnabled(3,False)
            
            
        ## Finally, execute the popup
        popUpMenu.exec_loop(QCursor.pos())   



    def selected_file_slot(self, listItem):
        """To know which item is selected """
                
        if listItem:
            self.m_listView_item_selected=str(listItem.text(0))
            
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
        self.textEdit_log.append("<info_tag> DS9 launched !!! </info_tag>")
        
    def start_aladin_slot(self):
        """Start Aladin tool"""
        
        
        os.system('echo "load %s ;sync; get vizier(2mass)" |/usr/local/bin/aladin -nobanner &' %(self.m_listView_item_selected))
        self.textEdit_log.append("<info_tag> Aladin launched !!! </info_tag>")
        
        # utils.runCmd does not allow launch in background !!
        #if utils.runCmd("/usr/local/bin/aladin &")==0: # some error
            #log.error("Error starting Aladin application")
            #raise Exception("Error executing command Aladin")
            #QMessageBox.warning(self,"Error","Error starting Alading application")
        #else:    
            #self.textEdit_log.append("<info_tag> Aladin launched !!! </info_tag>")

    def pushB_sel_masteDark_slot(self):
        source=QFileDialog.getOpenFileName( self.m_default_data_dir, "Master file (*.fit*)",  self, "Source Dialog","select master")

        if str(source):
            self.lineEdit_masterDark.setText(source)
            self.m_masterDark=str(source)

    def pushB_sel_masterFlat_slot(self):
        source=QFileDialog.getOpenFileName( self.m_default_data_dir, "Master file (*.fit*)",  self, "Source Dialog","select master")

        if str(source):
            fits = datahandler.ClFits(str(source))
            if fits.getType()!='MASTER_SKY_FLAT' or fits.getType()!='MASTER_DOME_FLAT' or fits.getType()!='MASTER_TW_FLAT':
                res=QMessageBox.information(self, "Info", QString("Selected frame does not look an MASTER FLAT.\n Continue anyway?"), QMessageBox.Ok, QMessageBox.Cancel)
                if res==QMessageBox.Cancel:
                    return
            self.lineEdit_masterFlat.setText(source)
            self.lineEdit_masterFlat_Filter.setText(fits.getFilter())
            self.m_masterFlat=str(source)
               

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
                             
#####################################################################################################
############# PROCESSING STAFF ######################################################################
#####################################################################################################
        
    def subtractFrames_slot(self):
        """This method is called to subtract two images selected from the File List View"""

        if (len(self.m_popup_l_sel)!=2):
            QMessageBox.critical(self, "Error", "You need to select  2 files")
        else:
            rb = reduce.ReductionBlock (self.m_popup_l_sel)
            outFilename = QFileDialog.getSaveFileName(self.m_outputdir+"/sub.fits", "*.fits", self, "Save File dialog")
            if not outFilename.isEmpty():
                try:
                    rb.mathOp('-', str(outFilename))
                except:
                    QMessageBox.critical(self, "Error", "Error while subtracting files")
                    raise
                else:
                    display.showFrame(str(outFilename))
                    self.textEdit_log.append("New file created: " + str(outFilename))


    def sumFrames_slot(self):
        """This methot is called to sum two images selected from the File List View"""

        if (len(self.m_popup_l_sel)<2):
            QMessageBox.critical(self, "Error", "You need to select  at least 2 files")
        else:
            rb = reduce.ReductionBlock (self.m_popup_l_sel)
            outFilename = QFileDialog.getSaveFileName(self.m_outputdir+"/sum.fits", "*.fits", self, "Save File dialog")
            if not outFilename.isEmpty():
                try:
                    rb.mathOp('+', str(outFilename))
                except:
                    QMessageBox.critical(self, "Error", "Error while summing files")
                    raise
                else:
                    display.showFrame(str(outFilename))
                    self.textEdit_log.append("New file created: " + str(outFilename))

    def divideFrames_slot(self):
        """This methot is called to divide two images selected from the File List View"""

        if (len(self.m_popup_l_sel)!=2):
            QMessageBox.critical(self, "Error", "You need to select  2 files")
        else:
            rb = reduce.ReductionBlock (self.m_popup_l_sel)
            outFilename = QFileDialog.getSaveFileName(self.m_outputdir+"/div.fits", "*.fits", self, "Save File dialog")
            if not outFilename.isEmpty():
                try:
                    rb.mathOp('/', str(outFilename))
                except:
                    QMessageBox.critical(self, "Error", "Error while dividing files")
                    raise
                else:
                    display.showFrame(str(outFilename))
                    self.textEdit_log.append("New file created: " + str(outFilename))

    def createMasterDark_slot(self):

        if len(self.m_popup_l_sel)<=1:
            QMessageBox.information(self,"Info","Not enought frames !")
            return

        outfileName = QFileDialog.getSaveFileName(self.m_outputdir+"/master_dark.fits", "*.fits", self, "Save File dialog")
        if not outfileName.isEmpty():
            try:
                texp_scale=False
                self.setCursor(Qt.waitCursor)
                self._task=reduce.calDark.MasterDark (self.m_popup_l_sel, "/tmp", str(outfileName), texp_scale)
                thread=reduce.ExecTaskThread(self._task.createMaster, self._task_info_list)
                thread.start()
            except Exception, e:
                self.setCursor(Qt.arrowCursor)
                QMessageBox.critical(self, "Error", "Error while creating master Dark. "+str(e))
                raise e
        else:
            pass

    def do_quick_reduction_slot(self):
        # Run quick-reduction mode with the
        self.QL2(self.m_listView_first_item_selected, self.textEdit_log)

    def createMasterDFlat_slot(self):
        
        """Create a master dome flat field from selected frames"""
        
        if len(self.m_popup_l_sel)>2:
            outfileName = QFileDialog.getSaveFileName(self.m_outputdir+"/master_dflat.fits", "*.fits", self, "Save File dialog")
            if not outfileName.isEmpty():
                try:
                    self.setCursor(Qt.waitCursor)
                    self._task = reduce.calDomeFlat.MasterDomeFlat (self.m_popup_l_sel, "/tmp", str(outfileName))
                    thread=reduce.ExecTaskThread(self._task.createMaster, self._task_info_list)
                    thread.start()
                except:
                    QMessageBox.critical(self, "Error", "Error while creating master Dome Flat")
                    raise
        else:
            QMessageBox.information(self,"Info","Error, not enought frames (>2) !")

    def createMasterTwFlat_slot(self):
      
        """Create a master twilight flat field from selected frames"""
        
        if len(self.m_popup_l_sel)>2:
            outfileName = QFileDialog.getSaveFileName(self.m_outputdir+"/master_twflat.fits", "*.fits", self, "Save File dialog")
            if not outfileName.isEmpty():
                try:
                    self.setCursor(Qt.waitCursor)
                    self._task = reduce.calTwFlat.MasterTwilightFlat (self.m_popup_l_sel, self.m_masterDark, str(outfileName))
                    thread=reduce.ExecTaskThread(self._task.createMaster, self._task_info_list)
                    thread.start()
                except:
                    log.error("Error creating master Twilight Flat file")
                    raise
        
    def createSkyFlat_slot_1(self):
        """ DOES NOT WORK FINE , TO_REVIEW!!!!"""
        
        rb = reduce.ReductionBlock (self.m_popup_l_sel)
        if len(self.m_popup_l_sel)>2:
            outfileName = QFileDialog.getSaveFileName(self.m_outputdir+"/super_sky_flat.fits", "*.fits", self, "Save File dialog")
            if not fileName.isEmpty():
                try:
                    rb.createSuperFlat(self.m_masterDark, str(outfileName))
                except:
                    QMessageBox.critical(self, "Error", "Error while creating super Flat")
                    raise
                else:
                    QMessageBox.information(self,"Info",QString("Super Flat  %1 created").arg(str(outfileName)))
                    display.showFrame(str(outfileName))
            else:
                pass
        else:
            QMessageBox.information(self,"Info","Error, not enought frames (>2) !")
    
    def createSkyFlat_slot_2(self):

        """ 
        Create a sky flat using the own science files selected on the main list view
        
        TODO: Check the selected images are SCIENCE frames with same filter !!! and expTime, .....
        
        """
        
        # Create file list with current selected science files
        filename=self.m_tempdir+"/files.list"      
        file_lst= open( filename, "w" )
        for file in self.m_popup_l_sel :
            file_lst.write( file +"\n")
        file_lst.close()
        # Call external apps
        cmd=self.m_papi_dir + "/mksuperflat.py" + " " + filename + " /tmp/skyFlat.fits"
        e=utils.runCmd(cmd)
        
        if e==1: # No error
            cmd2=self.m_papi_dir +"/irdr/bin/gainmap /tmp/skyFlat.fits /tmp/gain.fits 5 16 16 0.7 1.4"
            e=utils.runCmd(cmd2)  
            if e==1: # No error
                misc.fileUtils.removefiles("/tmp/skyFlat.fits")
                self.textEdit_log.append("<info_tag> Sky Super-flat /tmp/gain.fits created successful!!! </info_tag>")
                QMessageBox.information(self,"Info", QString("Super-Flat created (%1) !").arg("/tmp/gain.fits"))
                display.showFrame("/tmp/gain.fits")
                
            else:
                log.error("Error creating master Super Flat file")
                QMessageBox.critical(self, "Error", "Error while creating Gain Map")
        else:
            log.error("Error creating master Super Flat file")
            QMessageBox.critical(self, "Error", "Error while creating Super Flat")            
                
    def subtract_sky_slot(self):
        """ Subtract own image sky using SExtrator tool"""
        
        
        fits = datahandler.ClFits(self.m_listView_item_selected)
        if fits.getType()=='SCIENCE':
            sex_config = os.environ['IRDR_BASEDIR']+"/src/config/default.sex"
            minarea =  5
            threshold = 2.0
            input_file = self.m_listView_item_selected
            out_file = self.m_listView_item_selected.replace(".fits",".skysub.fits")
            cmd="sex %s -c %s -FITS_UNSIGNED Y -DETECT_MINAREA %s  -DETECT_THRESH %s  -CHECKIMAGE_TYPE -BACKGROUND -CHECKIMAGE_NAME %s" % (input_file, sex_config, str(minarea), str(threshold), out_file)  
            
            #Change to working directory
            os.chdir(self.m_papi_dir)
            #Change cursor
            self.setCursor(Qt.waitCursor)
            # Call external script (papi)
            self._proc=RunQtProcess(cmd, self.textEdit_log, self._task_info_list, out_file)      
            self._proc.startCommand()
        
        else:
            log.error("Sorry, selected file does not look a science file  ")
            QMessageBox.information(self, "Info", "Sorry, selected file does not look a science file ")                             
                         
    def subtract_nearSky_slot(self):
        """ Subtract nearest sky using skyfiler from IRDR package"""
      
        if self.checkBox_data_grouping.isChecked():
            log.error("Sorry, data grouping with POINT_NO, DITH_NO, EXPO_NO not yet implemented ")
            QMessageBox.information(self, "Info", "Sorry, data grouping with POINT_NO, DITH_NO, EXPO_NO not yet implemented ")
            return
        
        ra_dec_near_offset = self.lineEdit_ra_dec_near_offset.text().toInt()[0]/3600.0
        time_near_offset = self.lineEdit_time_near_offset.text().toInt()[0]/86400.0
        print "RA_DEC_OFFSET", ra_dec_near_offset
        print "TIME_NEAR_OFFSET", time_near_offset
      
        
        # Look for near science files
        if self.m_listView_item_selected:
            fits = datahandler.ClFits(self.m_listView_item_selected)
            if fits.getType()=='SCIENCE':
                near_list = datahandler.dataset.filesDB.GetFiles('ANY', 'SCIENCE', -1, fits.filter, fits.mjd, fits.ra, fits.dec,  ra_dec_near_offset*2, time_near_offset, runId=0)
                print "NEAR_LIST=", near_list
                # For the moment, the minimun number of nearest is >0
                if len(near_list)==0:
                    QMessageBox.information(self, "Info", "Not enought science frames found")  
                    return
            else:
                QMessageBox.information(self, "Info", "Selected frame is not a science frame") 
                return
            
            # Create file list nearest (ar,dec,mjd) from current selected science file
            filename=self.m_tempdir+"/files.list"      
            file_lst= open( filename, "w" )
            i=1
            my_list=""
            for f in near_list:
                file=str(f[0])
                my_list=my_list+file+"\n"
                if (file==self.m_listView_item_selected):
                    filen=i
                else:
                    i=i+1
                if (datahandler.ClFits(file).getType()!='SCIENCE'):
                    QMessageBox.critical(self, "Error", QString("File %1 is not a science frame").arg((file)))
                else:
                    file_lst.write( file +"\n")
            file_lst.close()
            res=QMessageBox.information(self, "Info", QString("Selected near frames are:\n %1").arg(my_list), QMessageBox.Ok, QMessageBox.Cancel)
            if res==QMessageBox.Cancel:
                return 
      
            # Call external apps skyfilter
            #gain=self.m_papi_dir+"/gain.fits"
            gain=self.m_masterFlat
            if not os.path.exists( gain ):
                QMessageBox.information(self, "Info", QString("No gain map %1 found. Please, select one on 'Calibrations' tab").arg(gain))
                return
                    
            #filen=near_list.index(self.m_listView_item_selected)
            hwidth=2
            cmd=self.m_papi_dir+"/irdr/bin/skyfilter_single %s %s %d nomask none %d" %(filename, gain, hwidth,filen)
            #Change to working directory
            os.chdir(self.m_papi_dir)
            #Change cursor
            self.setCursor(Qt.waitCursor)
            # Call external script (papi)
            self._proc=RunQtProcess(cmd, self.textEdit_log, self._task_info_list, self.m_listView_item_selected+".skysub")      
            self._proc.startCommand()
            
    
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
                    self._task = reduce.calBPM_2.BadPixelMask( "/tmp/bpm.list", dark, str(outfileName), lsig, hsig)
                    thread=reduce.ExecTaskThread(self._task.create, self._task_info_list)
                    thread.start()
                except:
                    QMessageBox.critical(self, "Error", "Not suitable frames to compute BPM.\n You need flat_off and flat_on frames")
                    raise
        else:
            QMessageBox.critical(self, "Error","Error, not suitable frames selected (flat fields)")


    def show_stats_slot(self):
        """Show image statistics in the log console of the files selected"""
        
        self.textEdit_log.append("<info_tag>FILE                     MEAN         MODE       STDDEV       MIN       MAX  </info_tag>")
        for item in self.m_popup_l_sel:
            values = (iraf.imstat (images=item,
            fields="image,mean,mode,stddev,min,max",format='no',Stdout=1))
            line=os.path.basename(values[0])
            self.textEdit_log.append(QString(str(line)))
        
    def background_estimation_slot(self):
        """ Give an background estimation of the current selected image"""
        
        cq = reduce.checkQuality.CheckQuality(self.m_popup_l_sel[0])
        try:     
            img=cq.estimateBackground(self.m_outputdir+"bckg.fits")
            
            values = (iraf.imstat (images=img,
            fields="image,mean,mode,stddev,min,max",format='no',Stdout=1))
            file,mean,mode,stddev,min,max=values[0].split()
            
            self.textEdit_log.append(QString("<info_tag> Background estimation MEAN= %1    MODE=%2    STDDEV=%3    MIN=%4         MAX=%5</info_tag>").arg(mean).arg(mode)
            .arg(stddev).arg(min).arg(max))
            
            display.showFrame(img)
        except:
            self.textEdit_log.append("<error_tag> ERROR: something wrong while computing background </error_tag>")
            raise
          
    def fwhm_estimation_slot (self):
        """ Give an FWHM estimation of the current selected image"""
        
        cq = reduce.checkQuality.CheckQuality(self.m_popup_l_sel[0])
        fwhm,std=cq.estimateFWHM()
        if fwhm>0:
            self.textEdit_log.append(QString("<info_tag> FWHM = %1 (pixels) std= %2 </info_tag>").arg(fwhm).arg(std))
        else:
            self.textEdit_log.append("<error_tag> ERROR: Cannot estimage FWHM with selected image  </error_tag>")           
            
    def createStackedFrame_slot(self):
        """ Compute a stacked frame (shift and aligned) from a set of nearest (ra,dec, mjd) frames, selected by user or automatic search"""
        
        filen=0 # really, not used for the moment
        
        # Check list lenght
        if len(self.m_popup_l_sel)<1:
            QMessageBox.information(self, "Info", "Not enought science frames selected")
            return
        
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
                near_list = datahandler.dataset.filesDB.GetFiles('ANY', 'SCIENCE', -1, fits.filter, fits.mjd, fits.ra, fits.dec,  ra_dec_near_offset*2, time_near_offset, runId=0)
                print "NEAR_LIST=", near_list
                # For the moment, the minimun number of nearest is >0
                if len(near_list)==0:
                    QMessageBox.information(self, "Info", "Not enought science frames found")  
                    return
            else:
                QMessageBox.information(self, "Info", "Selected frame is not a science frame") 
                return
            
           # Create nearest file list (ar,dec,mjd) from current selected science file
            filename=self.m_tempdir+"/files.list"      
            file_lst= open( filename, "w" )
            i=1
            my_list=""
            for f in near_list:
                file=str(f[0])
                my_list=my_list+file+"\n"
                if (file==self.m_listView_item_selected):
                    filen=i
                else:
                    i=i+1
                if (datahandler.ClFits(file).getType()!='SCIENCE'):
                    QMessageBox.critical(self, "Error", QString("File %1 is not a science frame").arg(file))
                else:
                    file_lst.write( file +"\n")
            file_lst.close() 
            resp=QMessageBox.information(self, "Info", QString("Selected near frames are:\n %1").arg(my_list), QMessageBox.Ok, QMessageBox.Cancel)
            if resp==QMessageBox.Cancel:
                return
                    
        # CASE 2: Stack frames  selected by user in the list_view
        elif len(self.m_popup_l_sel)>1:
            # Create file list from current selected science files
            filename=self.m_tempdir+"/files.list"      
            file_lst= open( filename, "w" )
            for file in self.m_popup_l_sel :
                if (datahandler.ClFits(file).getType()!='SCIENCE'):
                    QMessageBox.critical(self, "Error", QString("File %1 is not a science frame").arg(file))
                    file_lst.close()
                    return
                else:
                    file_lst.write( file +"\n")
            file_lst.close()
            filen=0
            
        #Change to working directory
        cwd=os.getcwd()
        os.chdir(self.m_papi_dir)
        cmd=self.m_papi_dir+"/papi_v1 %s single dither %s" %(filename,self.m_outputdir)
        #Change cursor
        self.setCursor(Qt.waitCursor)
        
        #err=utils.runCmd ( cmd )
        # Call external script (papi)
        self._proc=RunQtProcess(cmd, self.textEdit_log, self._task_info_list, self.m_outputdir+"/single.fits")      
        self._proc.startCommand() # asyncronous call than launch a subprocess to run the 'cmd' command
        
        #if err==0:
            #print "Some error while Stack generation (papi_v1) ..."
            #print "OUT > ", err
            #QMessageBox.critical(self,"Error", QString("Some error while sky subtraction : %1").arg(err))
        #else:
            ##Copy to output dir the stacked frame (output from PAPI)
            #shutil.move(self.m_papi_dir+"/single.fits", self.m_outputdir)
            #print "PASO move!!!"  
            #display.showFrame(self.m_outputdir+"/single.fits")
        
        # Return to the previus working directory
        #os.chdir(cwd)
        
    def createSuperMosaic_slot(self):
      
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
            filename=self.m_tempdir+"/mosaic_files.list"      
            file_lst= open( filename, "w" )
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
        os.chdir(self.m_papi_dir)
        cmd=self.m_papi_dir+"/build_mosaic %s" %(filename)
        #Change to working directory
        os.chdir(self.m_papi_dir)
        #Change cursor
        self.setCursor(Qt.waitCursor)
        # Call external script (papi)
        self._proc=RunQtProcess(cmd, self.textEdit_log, self._task_info_list, out_filename)      
        self._proc.startCommand()
      
      
    def do_raw_astrometry(self):
        """Compute an astrometric solution for the selected file in the main list view panel"""
        
        if len(self.m_popup_l_sel)==1:
            fits = datahandler.ClFits(self.m_listView_item_selected)
            if fits.getType()=='SCIENCE':
                # Call external script (papi)
                cmd=self.m_papi_dir+"/astrometry_scamp.pl 2mass regrid %s" %(self.m_listView_item_selected)
                
                #Change to working directory
                os.chdir(self.m_outputdir)
                #Change cursor
                self.setCursor(Qt.waitCursor)
                # Call external script (papi)
                self._proc=RunQtProcess(cmd, self.textEdit_log, self._task_info_list, None)      
                self._proc.startCommand()
                
            else:
                QMessageBox.information(self,"Info", QString("Sorry, but you need a reduced science frame."))
        
    

    def testSlot(self):

        rb = reduce.ReductionBlock (self.m_popup_l_sel)
        if rb.createBadPixelMask(self.m_outputdir+"/badPixelMask"):
            QMessageBox.information(self,"Info",QString("BPM file %1 created").arg(self.m_outputdir+"/badPixelMask.pl"))
        else:
            QMessageBox.information(self,"Info","Error, building Bad Pixel Mask (BPM) !")
        
        self.textEdit_log.append("<info_tag> Test finished !!! </info_tag>")
        
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

        
    
    #######################################################################################
    # Quick-Look 2: Preprocess frame. If the frame is  MEF, it will be stripped and
    # processed in parallel. Later, the result frame will be displayed  into DS9 display.
    #######################################################################################
    def QL2 ( self, frame_to_reduce, clog ):

        n_ext = 0
        mef_filenames = []
        result_frames = []

        if ( frame_to_reduce.endswith(".fits") ):
            clog.append('Start QuickLook-2 processing for frame: %s' %frame_to_reduce )
            f=datahandler.ClFits(frame_to_reduce)
            clog.append("Frame Type=%s" %f.type )

            # A SCIENCE frame
            if (f.type=="SCIENCE"):
                # A MEF SCIENCE frame
                if (f.mef==True):
                    # Split the Multi-Extension FITS file
                    n_ext = misc.fileUtils.splitMEF( frame_to_reduce, mef_filenames  )
                    par = True # for tests
                    #############################
                    #Execute PARALLEL reduction
                    #############################
                    if par:
                        #Synchronous call to parallel reduction
                        result_frames=self.do_parallel_reduc_MEF ( mef_filenames )

                    ##############################
                    #Execute SEQUENTIAL reduction
                    ##############################
                    else:
                        #Synchronous call to parallel reduction
                        result_frames=self.do_seq_reduc_MEF ( mef_filenames )
                    ###############################
                    #Finaly, DISPLAY the frames
                    ###############################
                    display.showSingleFrames( result_frames )
                #A Single Frame
                else:
                    # Single reduction
                    r=reduce.SimpleReduce()
                    #source_frame = '/disk-a/caha/panic/DATA/data_mat/QL1/orion0021_x4_1.fits'
                    #master_dark  = '/disk-a/caha/panic/DATA/data_mat/out/master_dark.fits'
                    #master_flat  = '/disk-a/caha/panic/DATA/data_mat/out/master_normflat.fits'
                    #out_dir      =  '/disk-a/caha/panic/DATA/data_mat/out/'
                    
                    appMask = False
                    try:
                        frame_out = r.run( frame_to_reduce, self.m_masterDark, self.m_masterFlat, self.m_outputdir+'/result.fits', self.m_tempdir, appPixMask=False )
                        #Finaly, display the frame
                        display.showFrame (frame_out)
                    except:
                        log.error("Error while single reduction")
                        QMessageBox.critical(self, "Error", "Error while single reduction of selected frame")
                        
                    #frame_out = frame_to_reduce.replace(".fits","_D_F_S.fits")
                    #threadsmod.ReduceThread(i, filenames[i-1], dark_frame, flat_frame, out_frame)
                    #prep = tasks.subtractSky1( frame, False)
                    #prep  = tasks.test2( frame, frame , True)
            # Not is a SCIENCE frame, it is a calibration, then show it directly
            else:
                clog.append('------>Frame is not a science frame, then it only will be displayed' )
                # Multi-Extension FITS file
                if f.mef==True:
                    #ds9 -zscale 'foo.fits[1]' 'foo.fits[2]' 'foo.fits[3]' foo.fits[4]' -tile
                    display.showFrame(frame_to_reduce)
                # Single FITS file
                else:
                    display.showFrame(frame_to_reduce)

      
    ################################################################################
    def do_parallel_reduc_MEF( self, filenames ):
      """Do a parallel reduction of a MEF file launching one thread for each extension"""
    
      log.debug("Start do_parallel_reduc_MEF")

      start = time.time()
      threads = []
      source_frames = ['/disk-a/caha/panic/DATA/data_mat/orion0021.fits','/disk-a/caha/panic/DATA/data_mat/orion0022.fits', '/disk-a/caha/panic/DATA/data_mat/orion0023.fits','/disk-a/caha/panic/DATA/data_mat/orion0024.fits']

      #self.m_masterDark   = '/disk-a/caha/panic/DATA/data_mat/out/master_dark.fits'
      #self.m_flat_frame   = '/disk-a/caha/panic/DATA/data_mat/out/master_normflat.fits'
      #self.out_frame    = '/disk-a/caha/panic/DATA/data_mat/out/prueba1.fits'
      appMask      = True
      result_frames=[]

      for i in range(1, len(filenames)+1):
        result_frames.append( filenames[i-1].replace(".fits", "_out.fit"))
        threads.append( reduce.ReduceThread(i, filenames[i-1], self.m_masterDark, self.m_masterFlat, result_frames[i-1], self.m_outputdir) )
        threads[i-1].start()

      for t in threads:
          t.join()

      print "Finalizaron todas las HEBRAS en %f secs !!!" %(time.time()-start)
      print "RESULT_FRAMES=", result_frames
      return result_frames


    ################################################################################
    def do_seq_reduc_MEF(self, filenames ):

      start = time.time()

      threads = []

      #dark_frame   = '/disk-a/caha/panic/DATA/data_mat/out/master_dark.fits'
      #flat_frame   = '/disk-a/caha/panic/DATA/data_mat/out/master_normflat.fits'
      #out_frame    = '/disk-a/caha/panic/DATA/data_mat/out/prueba1.fits'
      #appMask      = True
      result_frames = []
      
      for i in range(1, len(filenames)+1):
        result_frames.append( filenames[i-1].replace(".fits", "_out.fit"))
        t=reduce.ReduceThread(i, filenames[i-1], self.m_masterDark, self.m_masterFlat, result_frames[i-1], self.m_outputdir)
        t.start()
        t.join()

      print "Finalizaron todas las HEBRAS (secuencialmente) en %f secs !!!" %(time.time()-start)
      print "RESULT_FRAMES=", result_frames

      return result_frames

    ################################################################################

    def do_parallel_reduc_MEF_2( filenames ):
    # Doing Asynchronous callback with Pyro!! NOT USED !!! only for a test  !!

      import Pyro.naming, Pyro.core
      from Pyro.errors import NamingError

      start = time.time()

      objs = []

      source_frames = ['/disk-a/caha/panic/DATA/data_mat/orion0021.fits','/disk-a/caha/panic/DATA/data_mat/orion0022.fits', '/disk-a/caha/panic/DATA/data_mat/orion0023.fits','/disk-a/caha/panic/DATA/data_mat/orion0024.fits']
      dark_frame   = '/disk-a/caha/panic/DATA/data_mat/out/master_dark.fits'
      flat_frame   = '/disk-a/caha/panic/DATA/data_mat/out/master_normflat.fits'
      out_frame    = '/disk-a/caha/panic/DATA/data_mat/out/prueba1.fits'
      appMask      = False

      # locate the NS
      locator = Pyro.naming.NameServerLocator()
      #print 'Searching Name Server...',
      ns = locator.getNS()

      for i in range(1, len(source_frames)+1):

        # resolve the Pyro object
        #print 'finding object'
        try:
          name='sreduce_%d' %i
          URI=ns.resolve(name)
        #print 'URI:',URI
        except NamingError,x:
          print 'Couldn\'t find object, nameserver says:',x
          self.error=1
          raise

        # create a proxy for the Pyro object, and return that
        obj = Pyro.core.getProxyForURI(URI)
        obj._setOneway('run') # Make asynchronous callback 'run'
        error = obj.run ( source_frames[i-1], dark_frame, flat_frame, out_frame, False)
        objs.append(obj)


      #for i in range(1, len(filenames)+1):
      #  while objs[i-1].run_status!=1:
      #    pass


      print "Finalizaron todas las HEBRAS en %f secs !!!" %(time.time()-start)




################################################################################
#########   MAIN   #############################################################  
################################################################################

if __name__ == "__main__":
    app = QApplication(sys.argv)
    f = MainGUI()
    f.show()
    app.setMainWidget(f)
    app.exec_loop()
