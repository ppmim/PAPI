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
#              18/02/2014    jmiguel@iaa.es Speeded up with numpy.polynomial.polynomial.polyfit 
#                            and added support for MEFs.
#                            Changed structure of planes, p0=bias, p1=dark_curr
#
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
import astropy.io.fits as fits
import numpy

#import scipy.stats.stats


# Logging
from misc.paLog import log
from misc.version import __version__

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
        plane 0 = bias 
        plane 1 = dark current in DN/sec
        
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
        
        if not os.path.exists( os.path.abspath(os.path.join(self.__output_filename, os.pardir))):
            raise NameError, 'Wrong output path'
        if not self.__output_filename:
            log.error("Output DARK frame not defined")
            raise Exception("Wrong output filename")
    
        # Change to the source directory
        base, infile = os.path.split(self.__output_filename) 
        
        darks = numpy.zeros(nframes, dtype=numpy.int)
         
        # STEP 1: Check TYPE(dark), READMODE and read the EXPTIME of each frame
        # print "FRAMELIST= %s" %framelist
        i = 0
        f_readmode = -1
        f_n_extensions = -1
        for iframe in framelist:
            myfits = datahandler.ClFits(iframe)
            log.debug("Frame %s EXPTIME= %f TYPE= %s " %(iframe, myfits.exptime, myfits.type)) 
            # Check TYPE (dark)
            if not myfits.isDark():
                log.warning("Warning: Task 'createDarkModel' found a non dark frame. Skipping %s", iframe)
                darks[i] = 0
            else:
                # Check READMODE
                if (f_readmode !=-1 and (f_readmode != myfits.getReadMode() )):
                    log.error("Error: Task 'createMasterDark' finished. Found a DARK frame with different  READMODE")
                    darks[i] = 0  
                    #continue
                    raise Exception("Found a DARK frame with different READMODE") 
                else: 
                    f_readmode = myfits.getReadMode()
                    f_n_extensions = myfits.getNExt()
                    #log.debug("NEXT= %s"%(f_n_extensions))
                    darks[i] = 1
                
            i = i + 1
        log.debug('All frames checked')   
        
        naxis1 = myfits.naxis1
        naxis2 = myfits.naxis2            
        ndarks = (darks==1).sum()
        
        if ndarks < 2:
            log.error('Dark frameset doesnt have enough frames. At least 2 dark frames are needed')
            raise Exception("Dark sequence is too short. Al least 2 dark frames are needed !")
        
        # Initialize some storage arrays
        times = numpy.zeros(ndarks, dtype=numpy.float32)
        temp = numpy.zeros([ndarks, f_n_extensions, naxis1, naxis2], dtype=numpy.float32)
        out = numpy.zeros([2, f_n_extensions, naxis1, naxis2], dtype=numpy.float32)
        

        # Loop the images
        counter = 0 # counter for the number of darks
        for i in range(0, nframes):
            file = fits.open(framelist[i])
            f = datahandler.ClFits(framelist[i])
            if darks[i] == 1: # is a good dark
                for i_ext in range(0, f_n_extensions):
                    if f_n_extensions == 1:
                        temp[counter, 0, :,:] = file[0].data
                    else:
                        log.debug("Found MEF file")
                        temp[counter, i_ext, :,:] = file[i_ext + 1].data
                _mean = numpy.mean(temp[counter])
                _robust_mean = robust.r_nanmean(temp[counter].reshape(naxis1 * naxis2 * f_n_extensions))
                _median = robust.r_nanmedian(temp[counter])
                _mode = 3 * _median - 2 * _robust_mean
                if self.show_stats:
                    log.info("Dark frame TEXP=%s , ITIME=%s ,MEAN(not used)=%s , ROBUST_MEDIAN=%s ROBUST_MEAN=%s" % (f.expTime(), f.getItime(), _mean, _median, _robust_mean))
                times[counter] = float(f.expTime())
                counter = counter + 1
                file.close()

        log.debug("Now fitting the dark model...")
        # now collapse and fit the data
        # polyfit returns polynomial coefficients ordered from low to high.
        # It means, 0-coeff => bias, 1-coeff => dark_current 
        fit = numpy.polynomial.polynomial.polyfit(times, 
                            temp.reshape(len(times), naxis1 * naxis2 * f_n_extensions), deg=1)

        # Get the median value of the dark current                 
        median_dark_current = robust.r_nanmedian(fit[1])
        median_bias = robust.r_nanmedian(fit[0])

        log.info("MEDIAN_DARK_CURRENT = %s" % median_dark_current)
        log.info("MEDIAN BIAS = %s" % median_bias)    
        
        misc.fileUtils.removefiles( self.__output_filename )               

        # Write result in a FITS
        hdulist = fits.HDUList()
        hdr0 = fits.getheader(framelist[numpy.where(darks==1)[0][0]])
        prihdu = fits.PrimaryHDU (data = None, header = None)
        try:
            prihdu.header.set('INSTRUME', hdr0['INSTRUME'])
            prihdu.header.set('TELESCOP', hdr0['TELESCOP'])
            prihdu.header.set('CAMERA', hdr0['CAMERA'])
            prihdu.header.set('MJD-OBS', hdr0['MJD-OBS'])
            prihdu.header.set('DATE-OBS', hdr0['DATE-OBS'])
            prihdu.header.set('DATE', hdr0['DATE'])
            prihdu.header.set('UT', hdr0['UT'])
            prihdu.header.set('LST', hdr0['LST'])
            prihdu.header.set('ORIGIN', hdr0['ORIGIN'])
            prihdu.header.set('OBSERVER', hdr0['OBSERVER'])
        except Exception,e:
            log.warning("%s"%str(e))

        prihdu.header.set('PAPITYPE','MASTER_DARK_MODEL')
        prihdu.header.set('PAPIVERS', __version__, 'PANIC Pipeline version')
        
        prihdu.header.add_history('Dark model based on %s' % framelist)
        prihdu.header.add_history('Plane 0: bias ; Plane 1: dark current')
        
        if f_n_extensions>1:
            prihdu.header.set('EXTEND', True, after = 'NAXIS')
            prihdu.header.set('NEXTEND', f_n_extensions)
            prihdu.header.set('FILENAME', self.__output_filename)
            hdulist.append(prihdu)
            for i_ext in range(0, f_n_extensions):
                hdu = fits.PrimaryHDU()
                hdu.scale('float32') # important to set first data type
                hdu.data = fit.reshape(2, f_n_extensions, naxis1, naxis2)[:, i_ext, :, :]
                hdulist.append(hdu)
                del hdu
        else:
            prihdu.scale('float32') # important to set first data type
            prihdu.data = fit.reshape(2, 1, naxis1, naxis2)[:, 0, :, :]
            hdulist.append(prihdu)
         
        
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
    
    
        
        
