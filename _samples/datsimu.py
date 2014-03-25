#!/usr/bin/env python

################################################################################
#
# DatSimu
#
# DataSimu.py
#
# Last update 19/March/2008
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
import threading
import dircache
import getopt
import pyfits
import shutil


from threading import Thread


#################################################################################        

def run(args):

    print "Start the Data Simulator..."

    # Init the variables with default values
    source_path = ""
    dest_path = ""
    in_data_type = "all"
    delay = 1.0
    test  = False
    mef   = False
    
    # Get and check command-line options
    try:
        opts, args = getopt.getopt(args, "-s -d -t -D -m", ['source=','dest=','type=','delay=',"test", "mef"])
    except getopt.GetoptError:
        # print help information and exit:
        usage()
        sys.exit(2)

    #if (opts==[]):
    #    usage()
    #    sys.exit(2)
    
    print opts
    
    for option, parameter in opts:
        if option in ("-s", "--source"):
            source_path = parameter
            print "SOURCE directory =", source_path
        if option in ("-d", "--dest"):
            dest_path = parameter
            print "DEST directory =", dest_path
        if option in ("-t", "--type"):
            in_data_type = parameter
            print "in_data_type=", in_data_type
        if option in ("-D", "--delay"):
            delay = parameter
            print "delay=", delay
        if option in ("--test"):
            test=True
            print "Simulation  = True"
        if option  in ("--mef"):
            mef=True
            print "MEF = True"
            
    if  source_path=="" or dest_path=="":
        usage()
        sys.exit(2)
        
    list_s = dircache.listdir(source_path)
    for frame in list_s:
        toCopy=False
        if ( frame.endswith(".fits") ):
            data = pyfits.open(source_path + "/" + frame)
            read_type = data[0].header['OBJECT']
            print "FILE= %s ,TYPE = %s" %(frame, read_type)
            if  ( in_data_type == "" or in_data_type == "all" ):
                toCopy=True
            elif  (read_type.count(in_data_type)>0): #re.compile("dark",re.IGNORECASE).search(read_type, 1)):
                toCopy=True
        
            if  (toCopy == True):
                if ( not test ):
                    # If MEF is activated 
                    if ( mef ):
                        filelist[0]=source_path + "/" + frame
                        filelist[1]=filelist[0]
                        filelist[2]=filelist[0]
                        filelist[3]=filelist[0]
                        fname = frame.replace(".fits","_x4.fits")
                        if createMEF(filelist, dest_path + "/" + fname, 4 ):
                            print 'Created  file %s' %(dest_path + "/" + fname)
                            time.sleep(float(delay))
                    else:
                        try:
                            shutil.copyfile(source_path + "/" + frame, dest_path + "/" + frame)
                            print 'Copied %s file to %s' %(frame,dest_path)
                            time.sleep(float(delay))
                        except ValueError:
                            print "I/O error: file not copied"
                # Only a test
                else:
                    print 'Test to copy  %s file to %s' %(frame,dest_path)
        
          
    print "END the Data Simulator"

def usage ():
     print "Unknown command line parameter. Required parameters are : "
     print "-s / --source=		Source of data frames"
     print "-d / --dest=                Destiny of data frames"
     print "Optional parameter are : "
     print "-t / --type=                Type of data (default 'all')"
     print "-D / --Delay=               Delay in seconds between each copied image (default '1')"
     print "--test                      Test to copy files, not copying files"
     print "--mef                       Create Multi-(4)Extension FITS from a single file" 


#################################################################################
"""
Create a Multi Extension FITS from a set of 4 frames given
"""
def createMEF( filelist, filename, next=4):

    #STEP 1: Check all required frames are given
    if ( len(filelist) != next ):
        print "Error, Not enough number frames given. Required frames=", next
        return 0

    #STEP 2: Start the MEF creation
    #outfile = pyfits.open(filename, "append")
    first_file=True
    i=1
    for file in filelist:
        try:
            infile=pyfits.open(file)
        except IOError:
            #raise Exception("ERROR, Cannot open file")
            print "ERROR, Cannot open file " + file
            return 0

        print "File given ", file
        # Append the Primary Header    
        if first_file:
            phdr = pyfits.PrimaryHDU(header=infile[0].header,data=None) 
            phdr.header.update(key='EXTEND', value=pyfits.TRUE, after='NAXIS', comment="FITS dataset may contain extensions")
            phdr.header.update(key='NEXTEND', value=next, after='EXTEND',comment="Number of standard extensions")
            phdr.header.update(key='FILENAME', value=filename)
            #hdulist = pyfits.HDUList([pyfits.PrimaryHDU(header=phdr, data=None)])
            hdulist = pyfits.HDUList()
            hdulist.append(phdr)
            first_file=False

        #Append the ith Extension
        _name="IMAGE %d" %i
        hdu = pyfits.ImageHDU(data=infile[0].data, header=infile[0].header , name=_name) 
        hdulist.append(hdu) 
        infile.close()
        del infile, hdu
        i = i + 1

    # Write out the HDUlist (clobber=True, then overwrite output file if  exists)
    hdulist.verify('silentfix')
    hdulist.writeto(filename, output_verify='ignore', clobber=True)
    hdulist.close(output_verify='ignore')
    del hdulist
    print "End of createMEF"

    return 1

#################################################################################  

class myThread(Thread):
    
    def __init__(self):
        threading.Thread.__init__(self)
        
    def run(self):
            
        #The task to do
        
        print 'CHILD: run_darkcombine  finished OK'
        
################################################################################      
if __name__=="__main__":

    filelist=['/disk-a/caha/panic/DATA/data_mat/QL1/orion0028.fits','/disk-a/caha/panic/DATA/data_mat/QL1/orion0029.fits','/disk-a/caha/panic/DATA/data_mat/QL1/orion0030.fits','/disk-a/caha/panic/DATA/data_mat/QL1/orion0031.fits']
    
    #createMEF(filelist, "/disk-a/caha/panic/DATA/testMEF.fits", 4)

    run(sys.argv[1:])

#    aThread1 = myThread()
#    aThread1.start()
    
#    aThread2 = myThread()
#    aThread2.start()
    
#    while aThread1.isAlive():
#        print 'Thread1...parent\'s heartbeat...'
#        sleep(2)
        
#    while aThread1.isAlive():
#        print 'Thread2...parent\'s heartbeat...'
#        sleep(2)
    
    
#    print 'PARENT: Thread finished, but did apsum?'
        
################################################################################
