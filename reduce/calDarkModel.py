#!/usr/bin/env python

# Copyright (c) 2009-2012 IAA-CSIC  - All rights reserved. 
# Author: Jose M. Ibanez. 
# Instituto de Astrofisica de Andalucia, IAA-CSIC
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

from optparse import OptionParser
import sys
import os
import fileinput
import time

import misc.fileUtils
import misc.utils as utils
import datahandler
import misc.robust as robust

# Interact with FITS files
import pyfits
import numpy

#import scipy.stats.stats


# Logging
from misc.paLog import log

class MasterDarkModel(object):
    """
    Class used to build and manage a master calibration dark model
        
    As input a series of dark exposures with a range of exposure times is given. 
    A linear fit is done at each pixel position of data number versus 
    exposure time. A each pixel position in the output map represents the 
    slope of the fit done at that position and is thus the dark current 
    expressed in units of data numbers per second.   
    
    Parameters
    ----------
    input_data: list
        A list of dark files
    temp_dir: str
    	Directory for temporal files
    output_filename: str
    	Filename for the master dark obtained
    bpm: str
        Input bad pixel mask or NULL (optional)
    
    Returns
    -------
        If no error, a fits file (nx*ny) with 2 planes (extensions)
        plane 0 = dark current in DN/sec
        plane 1 = bias
        
        DARKCURRENT The median dark current in data numbers per second found 
        from the median value of the output dark current map.
    
    TODO
    ---- 
        - Data model for MEF files (PANIC)
        
    """
    def __init__(self, input_files, temp_dir='/tmp/', 
                 output_filename="/tmp/mdarkmodel.fits", 
                 bpm=None, show_stats=True):
        
        self.__input_files = input_files
        self.__output_filename = output_filename  # full filename (path+filename)
        self.__bpm = bpm
        self.show_stats = show_stats
        
    
    def createDarkModel(self):
      
        """
        Create a master DARK model from the dark file list
        """   
        log.debug("Start createDarkModel")
        start_time = time.time()
        t = utils.clock()
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
            log.error("Output DARK frame not defined")
            raise Exception("Wrong output filename")
    
        # Change to the source directory
        base, infile = os.path.split(self.__output_filename) 
        
        darks = numpy.zeros(nframes,dtype=numpy.int)
         
        # STEP 1: Check TYPE(dark),READMODE and read the EXPTIME of each frame
        #print "FRAMELIST= %s" %framelist
        i = 0
        f_readmode = -1
        for iframe in framelist:
            fits = datahandler.ClFits(iframe)
            log.debug("Frame %s EXPTIME= %f TYPE= %s " %(iframe, fits.exptime, fits.type)) 
            # Check TYPE (dark)
            if  not fits.isDark():
                log.warning("Warning: Task 'createDarkModel' found a non dark frame. Skipping %s", iframe)
                darks[i] = 0
            else:
                # Check READMODE
                if ( f_readmode!=-1 and (f_readmode!= fits.getReadMode() )):
                    log.error("Error: Task 'createMasterDark' finished. Found a DARK frame with different  READMODE")
                    darks[i] = 0  
                    #continue
                    raise Exception("Found a DARK frame with different READMODE") 
                else: 
                    f_readmode = fits.getReadMode()
                    darks[i] = 1
                
            i = i+1
        log.debug('All frames checked')   
        
        naxis1 = fits.naxis1
        naxis2 = fits.naxis2            
        ndarks = (darks==1).sum()
        
        if ndarks<2:
            log.error('Dark frameset doesnt have enough frames. At least 2 dark frames are needed')
            raise Exception("Dark sequence is too short. Al least 2 dark frames are needed !")
        
        #Initialize some storage arrays
        temp = numpy.zeros([ndarks, naxis1, naxis2], dtype=numpy.float)
        out = numpy.zeros([2, naxis1, naxis2], dtype=numpy.float)
        times = numpy.zeros(ndarks, dtype=numpy.float)
        
        #loop the images
        counter = 0
        for i in range(0, nframes):
            if darks[i]==1:
                file = pyfits.open(framelist[i])
                f = datahandler.ClFits ( framelist[i] )
                temp[counter, :,:] = file[0].data
                times[counter] = float(f.expTime())
                _mean = numpy.mean(file[0].data)
                _robust_mean = robust.mean(file[0].data)
                _median = numpy.median(file[0].data)
                _mode = 3*_median - 2*_mean
                #m_mode = my_mode(file[0].data)[0][0]
                #scipy_mode = scipy.stats.stats.mode ( file[0].data ) # extremely low !!
                if self.show_stats:
                    log.info("Dark frame TEXP=%s , ITIME=%s ,MEAN_VALUE=%s , MEDIAN=%s ROBUST_MEAN=%s"%(f.expTime(), f.getItime(), _mean, _median, _robust_mean))
                counter = counter+1
                file.close()
                
        log.debug("Now fitting the dark model...")
        #now collapse and fit the data
        slopes = numpy.zeros(naxis1*naxis2, dtype=numpy.float)
        bias = numpy.zeros(naxis1*naxis2, dtype=numpy.float)
        for i in range(0, naxis1):
            for j in range(0, naxis2):
                #result=numpy.polyfit(times, temp[:, i,j], deg=1) # result==> Polynomial coefficients, highest power first.
                b = (temp[counter-1,i,j]-temp[0,i,j])/(times[counter-1]-times[0]) # slope
                a = temp[0,i,j]-b*times[0] # intercept (bias)
                result = numpy.array([b,a])
                out[:, i,j] = result
                slopes[j+i*naxis2] = result[0]
                bias[j+i*naxis2] = result[1]
                #print "ROUND=%s %s result=%s" %(i,j, result)
        #out=numpy.where(out==0, numpy.polyfit(times, temp, deg=1), 0)   
        
        #Get the median value of the dark current                 
        median_dark_current = numpy.median(slopes)    #numpy.median(out[0,:,:])
        median_bias = numpy.median(bias) 
        
        log.info("MEDIAN_DARK_CURRENT = %s"%median_dark_current)
        log.info("MEDIAN BIAS = %s"% median_bias)    
        
        misc.fileUtils.removefiles( self.__output_filename )               

        # Write result in a FITS
        hdu = pyfits.PrimaryHDU()
        hdu.scale('float32') # important to set first data type
        hdu.data = out

        # copy some keywords 
        hdr0 = pyfits.getheader( framelist[numpy.where(darks==1)[0][0]]  )
        hdu.header.update('INSTRUME', hdr0['INSTRUME'])
        hdu.header.update('TELESCOP', hdr0['TELESCOP'])
        hdu.header.update('CAMERA', hdr0['CAMERA'])
        hdu.header.update('MJD-OBS', hdr0['MJD-OBS'])
        hdu.header.update('DATE-OBS', hdr0['DATE-OBS'])
        hdu.header.update('DATE', hdr0['DATE'])
        hdu.header.update('UT', hdr0['UT'])
        hdu.header.update('LST', hdr0['LST'])
        hdu.header.update('ORIGIN', hdr0['ORIGIN'])
        hdu.header.update('OBSERVER', hdr0['OBSERVER'])
        #
        hdu.header.update('PAPITYPE','MASTER_DARK_MODEL')
        hdu.header.add_history('Dark model based on %s' % framelist)
        hdulist = pyfits.HDUList([hdu])
        # write FITS
        try:
            hdulist.writeto(self.__output_filename)
            hdulist.close(output_verify='ignore')
        except Exception,e:
            log.error("Error writing dark model %s"%self.__output_filename)
            raise e
        
        log.debug('Saved DARK Model to %s' , self.__output_filename)
        log.debug("createDarkModel' finished %s", t.tac() )

        
        return self.__output_filename
        
################################################################################
def my_mode(data):
    """
    An easy (efficient and precise ??) way to find out the mode stats of an array
    (not used)
    """
    
    counts = {}
    for x in data.flatten():
        counts[x] = counts.get(x,0) + 1
    
    maxcount = max(counts.values())
    modelist = []
    
    for x in counts:
        if counts[x] == maxcount:
            modelist.append(x)
    
    return modelist,maxcount

# main
if __name__ == "__main__":
    usage = "usage: %prog [options]"
    desc = """
This module receives a series of FITS images (darks) with increasing exposure 
time and creates the master dark model and computes several statistics.
"""
    
    parser = OptionParser(usage, description=desc)
                  
    parser.add_option("-s", "--source",
                  action="store", dest="source_file_list",
                  help="Source file listing the filenames of dark frames.")
    
    parser.add_option("-o", "--output",
                  action="store", dest="output_filename", 
                  help="final coadded output image")
    
    parser.add_option("-S", "--show_stats",
                  action="store_true", dest="show_stats", default=False,
                  help="Show frame stats [default False]")    
    
    (options, args) = parser.parse_args()
    
    if len(sys.argv[1:])<1:
       parser.print_help()
       sys.exit(0)

    if not options.source_file_list or not options.output_filename:
        parser.print_help()
        parser.error("Incorrect number of arguments " )
    
    
    if os.path.isdir(options.source_file_list) or os.path.isdir(options.output_filename):
        parser.print_help()
        parser.error("Source and output must be a file, not a directory")
        
    filelist = [line.replace( "\n", "") 
                for line in fileinput.input(options.source_file_list)]
    
    try:
        mDark = MasterDarkModel(filelist, "/tmp", 
                                options.output_filename,
                                options.show_stats)
        mDark.createDarkModel()
    except Exception,e:
        log.error("Error computing dark model: %s"%str(e))
        sys.exit(0)
    
    
        
        
