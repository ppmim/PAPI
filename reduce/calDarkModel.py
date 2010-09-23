#!/usr/bin/env python
################################################################################
#
# PANICtool
#
# calDarkModel.py
#
# Created    : 17/06/2009    jmiguel@iaa.es
# Last update: 22/06/2009    jmiguel@iaa.es
#              03/03/2010    jmiguel@iaa.es Added READMODE checking 
#
################################################################################
# TODO: increase running speed !!!
################################################################################
# Import necessary modules

import getopt
import sys
import os
import logging
import fileinput
import time

import misc.fileUtils
import misc.utils as utils
import datahandler

# Interact with FITS files
import pyfits
import numpy

# Logging
from misc.paLog import log

class MasterDarkModel:
    """
    \brief Class used to build and manage a master calibration dark model
    
    \par Class:
        MasterDarkModel
    \par Purpose:
        Create a master Dark Model from a list for dark files
        
    \par Description:
         An input series dark exposures with a range of exposure times is given. A linear fit is done at each pixel position of data number versus exposure time. A each pixel position in the output map represents the slope of the fit done at that position and is thus the dark current expressed in units of data numbers per second.   
    \par Language:
        PyRaf
    \param input_data
        A list of dark files
    \param bpm
        Input bad pixel mask or NULL
    \retval 0
        If no error, a fits file (nx*ny) with 2 planes (extensions)
        plane 0 = bias
        plane 1 = dark current in DN/sec
        
        DARKCURRENT The median dark current in data numbers per second found from the median value of the output dark current map.
    \author
        JMIbannez, IAA-CSIC
        
    """
    def __init__(self, input_files, output_dir, output_filename="/tmp/mdarkmodel.fits", bpm=None):
        self.__input_files=input_files
        self.__output_file_dir=output_dir
        self.__output_filename=output_filename  # full filename (path+filename)
        self.__bpm=bpm
    
    def createDarkModel(self):
      
        """
        \brief Create a master DARK model from the dark file list
        """   
        log.debug("Start createDarkModel")
        start_time = time.time()
        t=utils.clock()
        t.tic()
        
        # Get the user-defined list of dark frames
        framelist = self.__input_files
        
        # STEP 0: Determine the number of darks frames to combine
        try:    
            nframes = len(framelist)
        except IndexError:
            log.error("No DARK frames defined")
            raise
        
        if not os.path.exists(os.path.dirname(self.__output_filename)):
            raise NameError, 'Wrong output path'
        if not self.__output_filename:
            log.error("Combined DARK frame not defined")
            raise "Wrong output filename"
    
        # Change to the source directory
        base, infile   = os.path.split(self.__output_filename) 
        
        darks=numpy.zeros(nframes,dtype=numpy.int)
         
        # STEP 1: Check TYPE(dark),READMODE and read the EXPTIME of each frame
        #print "FRAMELIST= %s" %framelist
        i=0
        f_readmode=-1
        for iframe in framelist:
            fits=datahandler.ClFits(iframe)
            log.debug("Frame %s EXPTIME= %f TYPE= %s " %(iframe, fits.exptime, fits.type)) 
            # Check TYPE (dark)
            if  not fits.isDark():
                log.warning("Warning: Task 'createDarkModel' found a non dark frame. Skipping %s", iframe)
                darks[i]=0
            else:
                # Check READMODE
                if ( f_readmode!=-1 and (f_readmode!= fits.getReadMode() )):
                    log.error("Error: Task 'createMasterDark' finished. Found a DARK frame with different  READMODE")
                    darks[i]=0  
                    #continue
                    raise Exception("Found a DARK frame with different  READMODE") 
                else: 
                    f_readmode=fits.getReadMode()
                    darks[i]=1
                
            i=i+1
        log.debug('All frames checked')   
        
        naxis1=fits.naxis1
        naxis2=fits.naxis2            
        ndarks=(darks==1).sum()
        
        if ndarks<2:
            log.error('Dark frameset doesnt have enough frames. At least 2 dark frames are needed')
            return False
        
        #Initialize some storage arrays
        temp=numpy.zeros([ndarks, naxis1, naxis2], dtype=numpy.float)
        out=numpy.zeros([2, naxis1, naxis2], dtype=numpy.float)
        times=numpy.zeros(ndarks, dtype=numpy.float)
        
        #loop the images
        counter=0
        for i in range(0,nframes):
            if darks[i]==1:
                file=pyfits.open(framelist[i])
                f=datahandler.ClFits ( framelist[i] )
                temp[counter, :,:]=file[0].data
                times[counter]=float(f.expTime())
                print "EXPTIME=", times[counter]
                counter=counter+1
                file.close()
                
        print "Now polyfit..."
        #now collapse and fit the data
        slopes=numpy.zeros(naxis1*naxis2, dtype=numpy.float)
        bias=numpy.zeros(naxis1*naxis2, dtype=numpy.float)
        for i in range(0, naxis1):
            for j in range(0, naxis2):
                #result=numpy.polyfit(times, temp[:, i,j], deg=1) # result==> Polynomial coefficients, highest power first.
                b=(temp[counter-1,i,j]-temp[0,i,j])/(times[counter-1]-times[0])
                a=temp[0,i,j]-b*times[0]
                result=numpy.array([b,a])
                out[:, i,j]=result
                slopes[j+i*naxis2]=result[0]
                bias[j+i*naxis2]=result[1]
                #print "ROUND=%s %s result=%s" %(i,j, result)
        #out=numpy.where(out==0, numpy.polyfit(times, temp, deg=1), 0)   
        
        #Get the median value of the dark current                 
        median_dark_current=numpy.mean(slopes)    #numpy.median(out[0,:,:])
        median_bias=numpy.mean(bias)
        print "MEDIAN_DARK_CURRENT=", median_dark_current
        print "MEDIAN BIAS=", median_bias    
        
        misc.fileUtils.removefiles( self.__output_filename )               
        # Write result in a FITS
        hdu = pyfits.PrimaryHDU()
        hdu.scale('float32') # importat to set first data type
        hdu.data=out     
        hdulist = pyfits.HDUList([hdu])
        hdu.header.update('OBJECT','MASTER_DARK_CURRENT')
        hdu.header.add_history('Dark model based on %s' % framelist)
        hdulist.writeto(self.__output_filename)
        hdulist.close(output_verify='ignore')
        #--
        #hdu.data=out[0,:,:]
        #f_dark=pyfits.HDUList([hdu])
        #f_dark.writeto("/tmp/darkc.fits")
        #hdu.data=out[1,:,:]
        #f_bias=pyfits.HDUList([hdu])
        #f_bias.writeto("/tmp/bias.fits")
        #--
        log.debug('Saved DARK Model to %s' , self.__output_filename)
        log.debug("createDarkModel' finished %s", t.tac() )
        
        return True
        
################################################################################
# Functions       
def usage ():
    print "Unknown command line parameter. Required parameters are : "
    print "-s / --source=      Source file list of data frames"
    print "-o / --out=         Output dark model filename "

################################################################################
# main
if __name__ == "__main__":
    print 'Start MasterDarkModel....'
    # Get and check command-line options
    args = sys.argv[1:]
    source_file_list = ""
    output_filename = ""
    try:
        opts, args = getopt.getopt(sys.argv[1:], "s:o:", ["source=","out="])
    except getopt.GetoptError:
        # print help information and exit:
        usage()
        sys.exit(1)

    
    for option, parameter in opts:
        print "OPTION=", opts
        if option in ("-s", "--source"):
            source_file_list = parameter
            print "Source file list =", source_file_list
            if not os.path.exists(os.path.dirname(source_file_list)):
                print 'Error, file list does not exists'
                sys.exit(1)
        if option in ("-o", "--out"):
            output_filename = parameter
            print "Output file =", output_filename
            
    if  source_file_list=="" or output_filename=="":
        usage()
        sys.exit(3)
    
    filelist=[line.replace( "\n", "") for line in fileinput.input(source_file_list)]
    #filelist=['/disk-a/caha/panic/DATA/ALHAMBRA_1/A0408060036.fits', '/disk-a/caha/panic/DATA/ALHAMBRA_1/A0408060037.fits']
    print "Files:",filelist
    mDark = MasterDarkModel(filelist,"/tmp",output_filename)
    mDark.createDarkModel()
    
        
