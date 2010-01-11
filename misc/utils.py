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

# PAPI modules
import misc.paLog
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
        
def runCmd( str_cmd, p_shell=True ):
    """ 
        DESCRIPTION
                A wrapper to run system commands  
        INPUTS
                str_cmd      - Command string to be executed in the shell
                p_shell      - if True (default), command will be executed through the shell, and all cout/cerr messages will be available
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
        p = subprocess.Popen(str_cmd, shell=p_shell, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
    except:
        log.error("Some error while running command...")
        raise 
       
    err=p.stderr.read()
    out=p.stdout.read()

    # IMPORTANT: Next checking (only available when shell=True) not always detect all kind of errors !!
    if (err.count('error') or err.count('ERROR') or out.count('error') or err.count('Segmentation fault') or err.count("command not found")
      or err.count("No such file or directory")):
        log.error( "An error happened while running command --> %s \n", err)
        #sys.exit(3)
        return 0 # ERROR
    
    else:
        return 1 # NO ERROR
        
        
                   