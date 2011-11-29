#!/usr/bin/env python

################################################################################
#
# DataCollector (PANIC DRS component)
#
# DataCollector.py
#
# Created     : 30/Oct/2008     jmiguel@iaa.es
# Last update : 12/Nov/2008     jmiguel@iaa.es
# 
################################################################################


import sys
import os
import fnmatch
import string
#import utils
import time
#import fileUtils
#import misc.dataset
#from threading import Thread
import copy
import fileinput
import glob

#PAPI modules
import datahandler
		
#Logging        
from misc.paLog import log       

        
class DataCollector (object):
	
    """
        \brief
        Class that implement the data receiver FITS data files comming from GEIRS  
    """
    
    def __init__(self, mode, source, filename_filter, p_callback_func):
    
        self.mode = mode
        self.filename_filter = filename_filter
        self.source = source
        self.callback_func = p_callback_func
		    
        self.dirlist = [] 
	    # Define the two lists containing the filenames of unprocessed and reduced
	    # files.
        self.newfiles     = []
        self.reducedfiles = []
        # next variable will include the files that was not able to read, in order
        # to avoid re-read again
        self.bad_files_found = []

        #Some flags
        self.stop = False
    
    def check(self):
	    
	    #print 'Start checking ...'
	    self.autoCheckFiles()
	    #print '...END checking'
    
    def Clear(self):
        self.dirlist = [] 
        self.newfiles = []
        self.reducedfiles = []
        self.bad_files_found = []
        
    def remove(self, pathname):
        """
        Remove a file from the current 'self.dirlist', so it could be again detected
        """
        
        self.dirlist.remove(pathname)
        #self.newfiles.remove(pathname)
        
    def SetFileFilter(self, new_filter):
        self.filename_filter = new_filter
        
    def autoCheckFiles(self):
        """
        Automatic checking of new files
        """
        # Is the source is exists (directory or file)
        if os.path.exists(self.source):	    
            # If this is the first time this routine is called, prepare the
            # directory list. 
            # ----> Nothing special to do ????  <------
            # Check for new files in the input directory
            self.findNewFiles()
        else:
            print "Source %s does NOT exists :"%self.source
                
    def __listFiles(self, dirpath):
        """
        Get directory listing sorted by creation date
        """
        #a = [os.path.join(dirpath, s) for s in os.listdir(dirpath)
        #    if os.path.isfile(os.path.join(dirpath, s))]
        a = [os.path.join(dirpath, s) for s in glob.glob(dirpath+"/"+self.filename_filter) \
              if os.path.isfile(os.path.join(dirpath, s))]
         
        a.sort(key=lambda s: os.path.getmtime(s))
        return a
       
    def __listFilesMJD(self, dirpath):
        """
    	-- NOT USED -- !!!
    	
    	Sort out input data files by MJD
    	
    	@NOTE: be careful, it could be a heavy routine    
        """
        
        dir_files = [os.path.join(dirpath, s) for s in glob.glob(dirpath+"/"+self.filename_filter) \
                     if os.path.isfile(os.path.join(dirpath, s))]
    	
        dataset = []
        
        for file in dir_files:
            try:
                fits = datahandler.ClFits(file)
            except Exception,e:
                print "[__listFilesMJD] Error reading file %s , skipped..."%(file)
                print str(e)
                self.bad_files_found.append(file)      
            else:
                dataset.append((file, fits.getMJD()))
        
        dataset = sorted(dataset, key=lambda data_file: data_file[1])          
        sorted_files = []
        for tuple in dataset:
            sorted_files.append(tuple[0])
    
    	return sorted_files
    
    def __sortFilesMJD(self, i_files):
        """
        Sort out input data files by MJD
        
        @NOTE: be careful, it could be a heavy routine    
        """
        
        dataset = []
        
        for file in i_files:
        	# filter out files already detected as bad files
        	if file not in self.bad_files_found:
	            try:
	                fits = datahandler.ClFits(file)
	            except Exception,e:
	                print "[__sortFilesMJD] Error reading file %s , skipped..."%(file)
	                print str(e)
	                self.bad_files_found.append(file)      
	            else:
	                dataset.append((file, fits.getMJD()))
        
        dataset = sorted(dataset, key=lambda data_file: data_file[1])          
        sorted_files = []
        for tuple in dataset:
            sorted_files.append(tuple[0])
    
        return sorted_files

    def read_GEIRS_fitsLog(self, type=1, start_datetime=None, end_datetime=None):
        """
        @summary: Function to read at self.source the GEIRS log files listing 
        the FITS files created. There are two options: 
            - ~GEIRS/log/save_CA2.2m.log
            - ~/tmp/fitsfiles.corrected
        
        @param type: 1 if source is a 'save_CA2.2m.log' file or 2 if it is a 
        fitsfiles.corrected.
        
        @param start_date: datetime object for the first file to consider, 
        where start to look for.
        @param end_date: datetime object for the last file to consider, 
        where end to look for.
        
        @return: A list with all the files read from the log file
                
        """  
    
        if start_datetime==None:
            l_start_datetime = datetime.datetime.now()-datetime.timedelta(days=1)    
        if end_datetime==None:
            l_end_datetime = datetime.datetime.now()
        if l_end_datetime < l_start_datetime:
             print "[DC] Error, end_datetime < start_datetime !"
             return []
         
         
        # Read the file contents from a generated GEIRS file
        if type == 1:
            for line in fileinput.input(self.source):
                sline = string.split(line)
                if sline[0]!="#":
                    contents.append(sline[6])
                    #print "FILE = ", sline[6]
        # To read ~/tmp/fitsfiles.corrected
        elif type == 2:
            # Read the file contents from a generated GEIRS file
            for line in fileinput.input(self.source):
                sline = string.split(line)
                if sline[0][0]!="#" and sline[1]!="ERROR":
                    contents.append(sline[1])
                    #print "FILE = ", sline[6]
                
        
        
    def findNewFiles(self):
    
        """
	    Find files that are not yet listed in the dirlist
	    ,and match these with the filename filter pattern
        """
	    
        # Retrieve the current value of the filename pattern from the widget
        pattern = self.filename_filter
        contents = []
        
        if self.mode=="dir":
            # Read the directory contents
            #contents = [os.path.join(self.source, file) for file in os.listdir(self.source)]
            contents = self.__listFiles(self.source)
        elif self.mode=="file":
            # Read the file contents
            contents = [line for line in fileinput.input(self.source)]
        # To read ~/GEIRS/log/save_CA2.2m.log
        elif self.mode=="geirs-file":
            # Read the file contents from a generated GEIRS file
            for line in fileinput.input(self.source):
                sline = string.split(line)
                if sline[0]!="#":
                    contents.append(sline[6])
                    #print "FILE = ", sline[6]
        # To read ~/tmp/fitsfiles.corrected
        elif self.mode=="geirs-file2":
            # Read the file contents from a generated GEIRS file
            for line in fileinput.input(self.source):
                sline = string.split(line)
                if sline[0][0]!="#" and sline[1]!="ERROR":
                    contents.append(sline[1])
                    #print "FILE = ", sline[6]
                
	    
        # Check the obtained list of files agains the existing directory list
	    # Remove files that already existed in the directory list
    
	    # ORDER of if-statements is important!
    
        iterlist = copy.copy(self.dirlist) # lists are mutables !
        for file in iterlist:
    
            if file not in contents:
                # Hmm... a strange situation. Apparently a file listed in self.dirlist
                # DISappeared from the directory. Adjust the lists accordingly
                print '[DC] File %s disappeared from directory - updating lists' % file
                self.dirlist.remove(file)
                self.callback_func(file+"__deleted__")
				
                ## new---
                #fc=fits.FitsFile(file)
                #fc.recognize()
                #fits.detected_files.remove( file )
                # Show into the frameLog
                #frameLog.put_warning (os.path.basename(file) + " --- " + "Removed !" )
                                ## wen---
                
                # Do NOT swap the following two statements!
                
                # Is this file already in the list of processed files?
                if file in self.reducedfiles: self.reducedfiles.remove(file)
                
                # Is this file already in the list of unprocessed files?
                if file in self.newfiles: self.newfiles.remove(file)
	            
				
            if file in contents:
                # Normal situation, all files in self.dirlist are also in the current
                # directory listing. Remove these files one-by-one from the list, so that
                # the remaining files are those that are the new files in the input
                # directory (this time).
                contents.remove(file)
			
				
        # ## 2011-09-12
        # Before adding to dirlist and process the new files, we sort out by MJD
        # ## 
        contents = self.__sortFilesMJD(contents)
        
        # Now loop over the remaining files
        for file in contents:

            # And append these to 'self.dirlist' for future reference
            self.dirlist.append(file)
            
            # Only look at the filename (disregard from directory path)
            basename = os.path.basename(file)
            
            # Make sure that (1) the 'file' is not a directory entry, (2) that it
            # matches the filename filter pattern, and (3) that it not yet in
            # the unprocessed or processed list.
            if   ( (not os.path.isdir(file))
			        and (fnmatch.fnmatch(basename, pattern))
			        and (file not in self.reducedfiles)    
			        and (file not in self.newfiles) 
			        ):
              
                #print '[DC] Appending file to the queue'
			    
                ## new---
                
	                #fc=fits.FitsFile(file)
	                ##fc.recognize()
                ##fits.detected_files.append( file )
                ## Append to the frameLog
	                #fileList.put (fc.filename , fc.type )
	                #misc.dataset.filesDB.insert(file)
	                #misc.dataset.filesDB.ListDataSet()
                
                ## wen---
                
                self.callback_func(file)
                # Only then, send message to the receiver client
                #self.newfiles.append(file) # removed line--> jmiguel - 2010-11-12
                print '[DC] Found new file ....%s' %file
            else:
                print "[DC] Warning, %s not a compliant file !!" %file
	
    
if __name__ == "__main__":
    def imprime_fichero(un_fichero):
        print un_fichero
    print 'DataCollector sample started ...'	 
    dr = DataCollector("geirs-file", "/home/panicmgr/GEIRS/log/save_CA2.2m.log", "*.fits", imprime_fichero)
    dr.check()
    print 'Datacollector sample finished successful'
    
	
