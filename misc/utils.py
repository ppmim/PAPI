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
# PANICtool
#
# utils.py
#
# Last update 24/09/2008
#
################################################################################

"""
   Utils module containing some useful class and funtions used in PANICtool
"""

# Import necessary modules

# System modules
import os
import time
import subprocess
import fileinput


# PAPI modules
from misc.paLog import log


class clock:

    def __init__(self):
        self.start=0

    def tic(self):
        self.start=time.time()

    def tac(self):
        lap=time.time()-self.start
        return "Elapsed time(s): %f" %lap
    
    
################################################################################
# Some useful functions for string handling (IRAF related)
################################################################################

def stringToList( a_string_list ):
    """ This function converts string list from IRAF format to python
        list format, using as separator ','
    """
    ret=[]
    # Split into a list
    ret = a_string_list.split(',')
    #Then, strip blanks from begining and end
    ret = [elem.strip() for elem in ret]
    if ret.count('')>0:
        ret.remove('')
    
    return ret
    
def listToString( a_list ):
    """ This function convert a list of string into a single string separating
        each list element with a ',' (IRAF file input format)
    """
    a_string=''
    n_tot = len(a_list)
    n=0
    for elem in a_list:
        n=n+1
        if n!=n_tot:
            a_string += elem + ' , '
        else:
            a_string += elem
    
    return a_string
        
def listToFile ( a_list, output_file=None ):
    """ This function write a list of filenames into a file """
    
    if output_file==None:
        temp_fd, temp_path = tempfile.mkstemp()
    else:        
        temp_fd=open(output_file,"w+") # truncate the file if exists
    
    for filename in a_list:
        temp_fd.write(filename.replace('//','/')+"\n")
    
    temp_fd.close()
    
    if output_file!=None: return output_file
    else: return temp_path
              
def fileToList ( input_file, out_filenames_list=None):
    """ This function dump the filenames from a input_file into a list"""
                  
    filelist=[line.replace( "\n", "") for line in fileinput.input(input_file)]
    if out_filenames_list!=None:
        out_filenames_list=filelist              
    
    return filelist
     
def runCmd( str_cmd, p_shell=True ):
    """ 
        DESCRIPTION
                A wrapper to run system commands  
        INPUTS
                str_cmd      - Command string to be executed in the shell
                p_shell      - if True (default), command will be executed through the shell, 
                                and all cout/cerr messages will be available
                             - if False, exception is the only way to find out problems during the call
          
        OUTPUT
                Return 0 if some errors 
                Retuen 1 if all was OK
                
        TODO 
                - allow to launch commands in background 
                - best checking of error when shell=True        
    """
           
    log.debug("Running command : %s \n", str_cmd)
    try:
        p = subprocess.Popen(str_cmd, bufsize=0, shell=True, stdin=subprocess.PIPE, 
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE, 
                             close_fds=True)
    except:
        log.error("Some error while running command...")
        raise 
    
    #Warning
    #We use communicate() rather than .stdin.write, .stdout.read or .stderr.read 
    #to avoid deadlocks due to any of the other OS pipe buffers filling up and 
    #blocking the child process.(Python Ref.doc)

    (stdoutdata, stderrdata) = p.communicate()
    err = stdoutdata + " " + stderrdata

    if len(err)>1:
        print "[runCmd]: STDOUT + STDERR = ", err
    
    
    # IMPORTANT: Next checking (only available when shell=True) not always detect all kind of errors !!
    if (err.lower().count('error ') or err.lower().count('error:') 
      or err.count('Segmentation fault') or err.count("command not found")
      or err.count('No source found')
      or err.count("No such file or directory")
      ):
        log.error("An error happened while running command --> %s \n" %err)
        return 0 # ERROR
    else:
        return 1 # NO ERROR


