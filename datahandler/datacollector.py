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
		
#Logging        
from misc.paLog import log       

        
class DataCollector ():
	
    """
        \brief
        Class that implement the data receiver FITS data files comming from GEIRS  
    """
    source=""
    mode="" # file or directory mode
    filename_filter="*.fits"
    checkID=""
    
    
    def __init__(self, mode, source, filename_filter, p_callback_func):
    
	    self.mode=mode
	    self.filename_filter = filename_filter
	    self.source = source
	    self.callback_func = p_callback_func
		    
	    self.dirlist = [] 
	    # Define the two lists containing the filenames of unprocessed and reduced
	    # files.
	    self.newfiles     = []
	    self.reducedfiles = []
	    self.stop = False
    
    def check(self):
	    
	    #print 'Start checking ...'
	    self.autoCheckFiles()
	    #print '...END checking'
    
    def Clear(self):
        self.dirlist = [] 
        self.newfiles     = []
        self.reducedfiles = []
        
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
        # Is the autocheck mode is a dir and exists
        if os.path.exists(self.source):	    
            # If this is the first time this routine is called, prepare the
            # directory list.
            if (self.checkID is None):
                if self.mode=="dir":
                    # Store the current contents of the input directory for future reference
                    #self.dirlist =  [os.path.join(self.source, file) for file in os.listdir(self.source)]
                    self.dirlist = self.__listFiles(self.source)
                elif self.mode=="file":
                    self.dirlist = [line for line in fileinput.input(self.source)]
                print 'Autochecking started'		
            # Check for new files in the input directory
            self.findNewFiles()
        else:
            # Reset the ID of the waiting loop
            self.checkID = None
            #messageLog.put('Autochecking stopped', 5)
            print 'Autochecking stopped'
                
    def __listFiles(self, dirpath):
        """
        Get directory listing sorted by creation date
        """
        a = [os.path.join(dirpath, s) for s in os.listdir(dirpath)
            if os.path.isfile(os.path.join(dirpath, s))]
        a.sort(key=lambda s: os.path.getmtime(s))
        return a
           
    def findNewFiles(self):
    
        """
	    Find files that are not yet listed in the dirlist
	    ,and match these with the filename filter pattern
        """
	    
        # Retrieve the current value of the filename pattern from the widget
        pattern = self.filename_filter
        contents=[]
        
        if self.mode=="dir":
            # Read the directory contents
            #contents = [os.path.join(self.source, file) for file in os.listdir(self.source)]
            contents = self.__listFiles(self.source)
        elif self.mode=="file":
            # Read the file contents
            contents = [line for line in fileinput.input(self.source)]
        elif self.mode=="geirs-file":
            # Read the file contents from a generated GEIRS file
            for line in fileinput.input(self.source):
                sline = string.split(line)
                if sline[0]!="#":
                    contents.append(sline[6])
                    #print "FILE = ", sline[6]
        
                
	    
        # Check the obtained list of files agains the existing directory list
	    # Remove files that already existed in the directory list
    
	    # ORDER of if-statements is important!
    
        iterlist = copy.copy(self.dirlist)
        for file in iterlist:
    
            if file not in contents:
                # Hmm... a strange situation. Apparently a file listed in self.dirlist
                # DISappeared from the directory. Adjust the lists accordingly
                print 'File %s disappeared from directory - updating lists' % file
                self.dirlist.remove(file)
                
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
                if file in self.newfiles:     self.newfiles.remove(file)
	            
    
            if file in contents:
                # Normal situation, all files in self.dirlist are also in the current
                # directory listing. Remove these files one-by-one from the list, so that
                # the remaining files are those that are the new files in the input
                # directory (this time).
                contents.remove(file)

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
              
                print '[DC] Appending file to the queue'
			    
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
                print "Warning, %s not a compliant file !!" %file
	
    
if __name__ == "__main__":
    print 'DataCollector sample started ...'	 
    dr = DataCollector("geirs-file", "/disk-a/caha/panic/GEIRS/log/save_CA2.2m.log", "*.fits", 0)
    dr.check()
    print 'Datacollector sample finished successful'
    
	
