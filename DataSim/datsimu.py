#!/usr/bin/env python

################################################################################
#
# DatSimu
#
# datsimu.py
#
# Last update 08/Oct/2009
#
################################################################################


    ####################################################################
    #                                                                  #
    #                   This is where it all starts.		       #
    #         The body of this program is at the end of this file      #
    #								       #
    ####################################################################


"""
   Main routines for the automatic data simulator (ADS).
"""

_version     = "1.0.0"
_date        = "19-03-2008"
_author      = "Jose M. Ibanez (jmiguel@iaa.es)"


_minversion_numpy   = "1.0.1"
_minversion_pyfits  = "1.1"
_minversion_pyraf   = "1.4"
_minversion_biggles = "1.6.4"


################################################################################

# Import necessary modules

#import config
#import messageLog
import os
import sys
import time
import dircache
import pyfits
import shutil
import fileinput
from optparse import OptionParser

#PAPI modules
import datahandler
import misc.mef

#################################################################################        

def run(source, dest_path, in_data_type="all", delay=1.0, test=False, mef=False):

    print "Start the Data Simulator..."

        
    # read source files
    list_s = []    
    if os.path.isfile(source):
        for file in fileinput.input(source):
            file = file.replace( "\n", "")
            print "[file] parsing file ->",file
            if file.endswith(".fits") or file.endswith(".fit"):
                list_s.append(file)
                print "MY_FILE",file
            else: print "Something is wrong......"
    elif os.path.isdir(source):
        for file in dircache.listdir(source):
            print "[dir] parsing file ->",file
            if file.endswith(".fits") or file.endswith(".fit"):
                list_s.append(source+"/"+file)           
                print "File added to list..."
 
    print "Sorting out files ...",list_s
    # sort out files
    list_s = sortOutData(list_s)
    #print "LIST",list_s
    
    print "Copying files ...."
    # procced to copy to destiny
    for frame in list_s:
        toCopy = False
        if ( frame.endswith(".fits") or frame.endswith(".fit") ):
            data = pyfits.open(frame, ignore_missing_end=True)
            read_type = data[0].header['OBJECT']
            print "FILE= %s ,TYPE = %s" %(frame, read_type)
            if  ( in_data_type == "" or in_data_type == "all" ):
                toCopy = True
            elif  (read_type.count(in_data_type)>0): #re.compile("dark",re.IGNORECASE).search(read_type, 1)):
                toCopy = True
        
            if (toCopy == True):
                if ( not test ):
                    # If MEF is activated 
                    if ( mef ):
                        filelist = frame*4 
                        fname = frame.replace(".fits","_x4.fits")
                        
                        mef = misc.mef.MEF(filelist)
                        if mef.createMEF(dest_path + "/" + fname):
                            print 'Created  file %s' %(dest_path + "/" + fname)
                            time.sleep(float(delay))
                    else:
                        try:
                            shutil.copy(frame, dest_path )
                            print 'Copied %s file to %s' %(frame, dest_path)
                            time.sleep(float(delay))
                        except Exception,e:
                            print "I/O error: file not copied"
                # Only a test
                else:
                    print 'Test to copy  %s file to %s' %(frame,dest_path)
        else:
	    print "File not ending in [.fit|.fits]" 
          
    print "END the Data Simulator"
    
def sortOutData(list):
    """
    Sort out input data files by MJD
    """
    
    dataset = []
        
    for file in list:
        try:
            fits = datahandler.ClFits(file)
        except:
            print "Error reading file %s , skipped..."%(file)      
        else:
            dataset.append((file, fits.getMJD()))
        
    dataset = sorted(dataset, key=lambda data_file: data_file[1])          
    sorted_files = []
    for tuple in dataset:
        sorted_files.append(tuple[0])
    
    return sorted_files
    


################################################################################
# main
if __name__ == "__main__":
    
    # Get and check command-line options
    usage = "usage: %prog [options] arg1 arg2"
    parser = OptionParser(usage)
    
    parser.add_option ("-s", "--source",
                  action = "store", dest = "source",
                  help = "Source (pathname of file-list) of data frames")
    
    parser.add_option ("-d", "--destiny",
                  action = "store", dest = "dest_path",
                  help = "Destiny pathname of data frames")
    
    parser.add_option ("-t", "--type",
                  action = "store", dest = "type", type="str",
                  help = "Type of files to copy(all,dark,flat,science) [default=%default]",
                  default="all")
    
    parser.add_option ("-D", "--Delay",
                  action = "store", dest = "delay", default=1.0, type=int,
                  help = "Delay between file copy [default=%default]")
    
    parser.add_option ("-T", "--test",
                  action = "store_true", dest = "test", default=False,
                  help = "Does only a test, not doing the file copying [default=%default]")
    
    parser.add_option("-M", "--mef",
                  action = "store_true", dest = "mef", default = False,
                  help = "Build MEF file prior files are copied [default=%default]")
    

    (options, args) = parser.parse_args()
    
    
    if not options.source or len(args)!=0: 
    # args is the leftover positional arguments after all options have been processed
        parser.print_help()
        parser.error("wrong number of arguments " )
        
    run(options.source, options.dest_path, options.type, options.delay, 
        options.test, options.mef)

################################################################################
