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
import datetime as dt
import logging

#PAPI modules
import datahandler
		
#Logging        
from misc.paLog import log       

        
class DataCollector (object):
	
    """
    Class that implement the data receiver FITS data files comming from GEIRS  
    """
    
    def __init__(self, mode, source, filename_filter, p_callback_func):
        """
        Initialize object
        
        Parameters
        ----------
        mode : str
            'dir' - Read/inspect the files files from the given directory in 'source'
            'file' - Read/inspect the files from the given file given in 'source'
            'geirs-file' - Read the file contents from a generated GEIRS file of 
            type 1, i.e, ~GEIRS/log/save_CA2.2m.log
            'geirs-file2' - Read the file contents from a generated GEIRS file of 
            type 1, i.e., ~/tmp/fitsfiles.corrected
            
        source: str
            directory name or filename to be read for inputs
            
        filename_filter : str
            filter for files to look for (i.e: *.fits)
            
        p_callback_func : str
            function name to be executed each time a new file is read/detected.
            This function usually has as first parameter the filename just read.
            
        """
        
        self.mode = mode
        self.filename_filter = filename_filter
        self.source = source
        self.callback_func = p_callback_func
        
        self.dirlist = [] 
        # Define the two lists containing the filenames of unprocessed and reduced
        # files.
        self.newfiles     = []
        self.reducedfiles = []
        # next variable includes the files that was not able to read, 
        # in order to not try to read them again
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
        self.bad_files_found = [] # we'll try to read again
        
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
        Get directory list sorted by creation date.
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
        
        """            
        # Now, try to re-read the bad_files_found
        # In fact the next code is no necessary because it is done in clfits:recognize()
        sleep_time = 0.0
        while len(self.bad_files_found)>0 and sleep_time<5.0:
                sleep_time +=0.5
                print "[__sortFilesMJD] Start to re-read bad_files_found (sleep_time=%f)"%sleep_time
                for file in self.bad_files_found:
                    time.sleep(sleep_time)
                    try:
                        fits = datahandler.ClFits(file)
                        self.bad_files_found.remove(file)
                        print "[__sortFilesMJD] file successfully re-read !"
                    except Exception,e:
                        print "[__sortFilesMJD] Error re-reading file %s , skipped..."%(file)
                        print str(e)
                    else:
                        dataset.append((file, fits.getMJD()))
        #Files that failed to re-read are discarded for ever
        self.discarded_files += self.bad_files_found
        self.bad_files_found = []
        """
        
        dataset = sorted(dataset, key=lambda data_file: data_file[1])          
        sorted_files = []
        for l_tuple in dataset:
            sorted_files.append(l_tuple[0])
    
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
    
        contents = []
        
        if start_datetime==None:
            # by default, we limit the files one day old
            l_start_datetime = dt.datetime.now()-dt.timedelta(days=1)    
        if end_datetime==None:
            l_end_datetime = dt.datetime.now()
        if l_end_datetime < l_start_datetime:
            print "[DC] Error, end_datetime < start_datetime !"
            return []
        
        #print "start_date = %s"%l_start_datetime
        #print "end_date = %s"%l_end_datetime
         
        # Read the file contents from a generated GEIRS file
        if type == 1:
            for line in fileinput.input(self.source):
                sline = string.split(line)
                if sline[0]!="#":
                    try:
                        line_date = dt.datetime.strptime(sline[0]+" "+sline[1],
                                                     "%Y-%m-%d %H:%M:%S")
                    except ValueError,e:
                        print "Error, cannot read datetime stamp in log file line :",line
                        continue
                    if line_date > l_start_datetime and line_date < l_end_datetime:
                        contents.append(sline[6])
                        print "FILE = ", sline[6]
                    else:
                        pass
                        #print "File too old ...."
                        #print "file datetime = %s"%line_date
        # To read ~/tmp/fitsfiles.corrected
        elif type == 2:
            # Read the file contents from a generated GEIRS file
            for line in fileinput.input(self.source):
                sline = string.split(line)
                if sline[0][0]!="#" and sline[1]!="ERROR":
                    try:
                        line_date = dt.datetime.strptime(sline[0],"%Y-%m-%d_%Hh%Mm%S")
                    except ValueError,e:
                        print "Error, cannot read datetime stamp in log file line",line
                        continue
                    if line_date > l_start_datetime and line_date < l_end_datetime:
                        contents.append(sline[1])
                        #print "FILE = ", sline[1]
                    else:
                        pass
                        #print "File too old ...."
                        #print "file datetime = %s"%line_date
        
        return contents        
        
        
    def findNewFiles(self):
        """
        Find files that are not yet listed in the dirlist
        ,and match these with the filename filter pattern
        """

        # Retrieve the current value of the filename pattern from the widget
        pattern = self.filename_filter
        contents = []
        
        try:
            # To Read the directory contents
            if self.mode=="dir":
                #contents = [os.path.join(self.source, file) for file in os.listdir(self.source)]
                contents = self.__listFiles(self.source)
            # To Read a simple text file
            elif self.mode=="file":
                contents = [line for line in fileinput.input(self.source)]
            # To read ~/GEIRS/log/save_CA2.2m.log
            elif self.mode=="geirs-file":
                contents = self.read_GEIRS_fitsLog(type=1)
                """
	            # Read the file contents from a generated GEIRS file
	            for line in fileinput.input(self.source):
	                sline = string.split(line)
	                if sline[0]!="#":
	                    contents.append(sline[6])
	                    #print "FILE = ", sline[6]
	            """
            # To read ~/tmp/fitsfiles.corrected
            elif self.mode=="geirs-file2":
                contents = self.read_GEIRS_fitsLog(type=2)
                """
	            # Read the file contents from a generated GEIRS file
	            for line in fileinput.input(self.source):
	                sline = string.split(line)
	                if sline[0][0]!="#" and sline[1]!="ERROR":
	                    contents.append(sline[1])
	                    #print "FILE = ", sline[6]
	            """    
         
        except Exception, e:
            print "Some error while reading source  %s "%self.source
            return

        # Check the obtained list of files agains the existing directory list
        # Remove files that already existed in the directory list
    
        # ORDER of if-statements is important!
    
        iterlist = copy.copy(self.dirlist) # lists are mutable !
        for file in iterlist:
            if file not in contents:
                # Hmm... a strange situation. Apparently a file listed in self.dirlist
                # DISappeared from the directory. Adjust the lists accordingly
                print '[DC] File %s disappeared from directory - updating lists' % file
                self.dirlist.remove(file)
                self.callback_func(file+"__deleted__")
                
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
        # Only when mode=dir, because it is supposed in 'file's-modes are already sorted
        if self.mode=="dir":
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
            if ( (not os.path.isdir(file))
			        and (fnmatch.fnmatch(basename, pattern))
			        and (file not in self.reducedfiles)    
			        and (file not in self.newfiles)
                    and (file not in self.bad_files_found)
			        ):
              
                # Try-to-read the file (only for integrity check)
                try:
                    # Now, try to read the file
                    fits = datahandler.ClFits(file)
                    del fits
                except Exception,e:
                    print "[findNewFiles] Error reading file %s , skipped..."%(file)
                    print str(e)
                    self.remove(file)
                    self.bad_files_found.append(file)
                else:
                    self.callback_func(file)
                    # Only then, send message to the receiver client
                    #self.newfiles.append(file) # removed line--> jmiguel - 2010-11-12
                    #print '[DC] Found new file ....%s' %file
            else:
                print "[DC] Warning, %s not a compliant file or does not exist !!" %file

        # if it is the last file of a full-directory list, then notify QL to
        # update the ListView.
        # It is done to avoid overloading the QL update of the ListView each time
        # when a long directory is loaded in one go. Howerver, last file is not 
        # inserted twice.
        if len(contents)>0:
            self.callback_func(contents[-1]+"__last__")


    
if __name__ == "__main__":
    def imprime_fichero(un_fichero):
        print un_fichero
    print 'DataCollector sample started ...'	 
    dr = DataCollector("geirs-file", "/home/panicmgr/GEIRS/log/save_CA2.2m.log", "*.fits", imprime_fichero)
    dr.check()
    print 'Datacollector sample finished successful'
    

