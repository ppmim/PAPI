#! /usr/bin/env python
#encoding:UTF-8

# Copyright (c) 2008-2019 Jose M. Ibanez All rights reserved.
# Institute of Astrophysics of Andalusia, IAA-CSIC
#
# This file is part of PAPI
#
# PAPI is free software: you can redistribute it and/or modify
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
# DataCollector (PANIC DRS component)
#
# datacollector.py
#
# Created     : 30/Oct/2008     jmiguel@iaa.es
# 
################################################################################

import os
import fnmatch
import string
#import misc.dataset
import copy
import fileinput
import glob
import datetime as dt

# PAPI modules
import datahandler
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
            type 2, i.e., ~/tmp/fitsGeirsWritten (old fitsfiles.corrected)
            
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
        
        # Define the two lists containing the filenames of unprocessed and
        # reduced files.
        self.newfiles = []
        self.reducedfiles = []
        
        # Next variable includes the files that was not able to read, 
        # in order to not try to read them again
        self.bad_files_found = []

        # Files that gave and error and need to be re-read next time, upto 
        # '_n_retries_' times
        self.pend_to_read = {}
        self._n_retries_ = 20

        # Some flags
        self.stop = False
    
    def check(self):
        self.autoCheckFiles()

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
            print("Source %s does NOT exists :" % self.source)
                
    def __listFiles(self, dirpath):
        """
        Get directory list sorted by creation date.
        
        Note: only .fits or .fit files are considered, not matter 
        what kind of filename_filter is set.
        """
        
        #a = [os.path.join(dirpath, s) for s in os.listdir(dirpath)
        #    if os.path.isfile(os.path.join(dirpath, s))]
        a = [os.path.join(dirpath, s) for s in glob.glob(dirpath + "/" + self.filename_filter) \
              if os.path.isfile(os.path.join(dirpath, s)) and ( 
                 os.path.join(dirpath, s).endswith('.fits') or os.path.join(dirpath, s).endswith('.fit')) ]
         
        a.sort(key=lambda s: os.path.getmtime(s))
        return a
       
    def __sortFilesMJD(self, i_files):
        """
        Sort out input data files by MJD
        
        NOTE 1: this routine takes into account whether a file is still being saved
        and then try to read it again later upto '_n_retries_' times.

        NOTE 2: be careful, it could be a heavy routine    
        """
        
        dataset = []
        for file in i_files:
            # filter out files already detected as bad files
            if file not in self.bad_files_found:
                try:
                    fits = datahandler.ClFits(file, check_integrity=True)
                except IOError as e:
                    if file in self.pend_to_read:
                        if  self.pend_to_read[file] < self._n_retries_:
                            self.pend_to_read[file] = self.pend_to_read[file] + 1
                        else:
                            # definitely, file is discarted
                            self.bad_files_found.append(file)
                            del self.pend_to_read[file]
                            print("[__sortFilesMJD] Definitely file %s , is discarted"%(file))
                    else:
                        self.pend_to_read[file] = 1
                except Exception as e:
                    print("[__sortFilesMJD] Error reading file %s , skipped..." %(file))
                    print(str(e))
                    self.bad_files_found.append(file)      
                else:
                    dataset.append((file, fits.getMJD()))
        
        
        dataset = sorted(dataset, key=lambda data_file: data_file[1])          
        sorted_files = []
        for l_tuple in dataset:
            sorted_files.append(l_tuple[0])
    
        return sorted_files

    def read_GEIRS_fitsLog(self, type=1, start_datetime=None, end_datetime=None):
        """
        @summary: Function to read at self.source the GEIRS log files listing 
        the FITS files created. There are two options: 
            - ~/GEIRS/log/save_CA2.2m.log
            - ~/tmp/fitsGeirsWritten
        
        @param type: 1 if source is a 'save_CA2.2m.log' file or 2 if it is a 
        'fitsGeirsWritten'.
        
        @param start_date: datetime object for the first file to consider, 
        where start to look for.
        @param end_date: datetime object for the last file to consider, 
        where end to look for.
        
        @return: A list with all the files read from the log file
                
        """  
    
        contents = []
        
        if start_datetime == None:
            # by default, we limit the files one day old
            l_start_datetime = dt.datetime.now() - dt.timedelta(days=1)    
        if end_datetime == None:
            l_end_datetime = dt.datetime.now()
        if l_end_datetime < l_start_datetime:
            print("[DC] Error, end_datetime < start_datetime !")
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
                    except ValueError as e:
                        print("Error, cannot read datetime stamp in log file line :", line)
                        print(str(e))
                        continue
                    if line_date > l_start_datetime and line_date < l_end_datetime:
                        contents.append(sline[6])
                        print("FILE = ", sline[6])
                    else:
                        pass
                        # print "File too old ...."
                        # print "file datetime = %s"%line_date
        # To read ~/tmp/fitsGeirsWritten
        elif type == 2:
            # Read the file contents from a generated GEIRS file
            for line in fileinput.input(self.source):
                sline = string.split(line)
                if sline[0][0]!="#" and sline[1]!="ERROR":
                    try:
                        # The lines in the log file have two fields: timestamp+filename
                        # timestap = `date --iso-8601=seconds` --> 2014-10-23T12:32:38+0200
                        # Due to datetime.strptime does not work with %z directive,
                        # the UTC offset is skipped. 
                        line_date = dt.datetime.strptime(sline[0].split("+")[0],"%Y-%m-%dT%H:%M:%S")
                    except ValueError as e:
                        print("Error, cannot read datetime stamp in log file line", line)
                        print(sline[0])
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
            if self.mode == "dir":
                #contents = [os.path.join(self.source, file) for file in os.listdir(self.source)]
                contents = self.__listFiles(self.source)
            # To Read a simple text file
            elif self.mode == "file":
                contents = [line for line in fileinput.input(self.source)]
            # To read ~/GEIRS/log/save_CA2.2m.log
            elif self.mode == "geirs-file":
                contents = self.read_GEIRS_fitsLog(type=1)
                """
	            # Read the file contents from a generated GEIRS file
	            for line in fileinput.input(self.source):
	                sline = string.split(line)
	                if sline[0]!="#":
	                    contents.append(sline[6])
	                    #print "FILE = ", sline[6]
	            """
            # To read ~/tmp/fitsGeirsWritten 
            elif self.mode == "geirs-file2":
                contents = self.read_GEIRS_fitsLog(type=2)
                """
	            # Read the file contents from a generated GEIRS file
	            for line in fileinput.input(self.source):
	                sline = string.split(line)
	                if sline[0][0]!="#" and sline[1]!="ERROR":
	                    contents.append(sline[1])
	                    #print "FILE = ", sline[6]
	            """    
         
        except Exception as e:
            print("Some error while reading source  %s " % self.source)
            return

        # Check the obtained list of files agains the existing directory list
        # Remove files that already existed in the directory list
    
        # ORDER of if-statements is important!
    
        iterlist = copy.copy(self.dirlist) # lists are mutable !
        for file in iterlist:
            if file not in contents:
                # Hmm... a strange situation. Apparently a file listed in self.dirlist
                # DISappeared from the directory. Adjust the lists accordingly
                print('[DC] File %s disappeared from directory - updating lists' % file)
                self.dirlist.remove(file)
                self.callback_func(file + "__deleted__")
                
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
        if self.mode == "dir":
            # check of data integrity is also done at sortFilesMJD
            contents = self.__sortFilesMJD(contents)
            
        # Now loop over the remaining files (the new files arrived !)
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
                    if self.mode != "dir": # for dir, it was already done above in __sortFilesMJD()
                        fits = datahandler.ClFits(file, check_integrity=True)
                        del fits
                except IOError as e:
                    # file is still being saved ! should not happen when reading GEIRS logs
                    print("[findNewFiles-1] Error reading file %s , skipped..." %(file))
                    print(str(e))
                    self.remove(file)
                    self.bad_files_found.append(file)    
                except Exception as e:
                    print("[findNewFiles-2] Error reading file %s , skipped..."%(file))
                    print(str(e))
                    self.remove(file)
                    self.bad_files_found.append(file)
                else:
                    #print "New File to be inserted : %s"%file
                    self.callback_func(file)
                    # Only then, send message to the receiver client
                    #self.newfiles.append(file) # removed line--> jmiguel - 2010-11-12
                    #print '[DC] Found new file ....%s' %file
            else:
                print("[DC] Warning, %s not a compliant file or does not exist !!" %file)

        # if it is the last file of a full-directory list, then notify QL to
        # update the ListView.
        # It is done to avoid overloading the QL update of the ListView each time
        # when a long directory is loaded in one go. Howerver, last file is not 
        # inserted twice.
        if len(contents) > 0:
            self.callback_func(contents[-1] + "__last__")


if __name__ == "__main__":
    def imprime_fichero(f):
        print(f)
    print('DataCollector sample started ...')
    dr = DataCollector("geirs-file", "/home/panicmgr/GEIRS/log/save_CA2.2m.log",
                       "*.fits", imprime_fichero)
    dr.check()
    print('Datacollector sample finished successful')
    

