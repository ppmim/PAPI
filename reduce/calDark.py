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
# calDark.py
#
# Created    : 07/11/2008    jmiguel@iaa.es
# Last update: 29/09/2009    jmiguel@iaa.es - Scale by EXPTIME
#              11/12/2009    jmiguel@iaa.es - Rename output filename to include EXPTIME and NCOADDS and use ClFits class
#              14/12/2009    jmiguel@iaa.es - Skip non DARK frames and cotinue working with the good ones (good_frames)
#              02/03/2010    jmiguel@iaa.es - added READEMODE checking
#              14/09/2010    jmiguel@iaa.es - added support to MEF files, calling mscred.darkcombine subrutine instead of imred.darkcombine
#              07/03/2011    jmiguel@iaa.es - Added Stats output and normalization (divide master dark by the TEXP to get a master dark in ADU/s units)
#
# TODO
#  - checking of ITIME ( and not only EXPTIME, NCOADDS )
################################################################################

################################################################################
# Import necessary modules

import getopt
import sys
import os
import logging
import fileinput
import time
import shutil
from optparse import OptionParser

import misc.fileUtils
import misc.utils as utils
import datahandler

# Pyraf modules
from pyraf import iraf
from iraf import noao
from iraf import imred
from iraf import ccdred
from iraf import mscred

# Interact with FITS files
import pyfits

# Import Pyro core
#import Pyro.core
#import Pyro.naming

# Logging
from misc.paLog import log

class MasterDark:
    """
    \brief Class used to build and manage a master calibration dark 
    
    \par Class:
        MasterDark
    \par Purpose:
        Create a master Dark from a list for dark files (single or MEF files)
    \par Description:
            
    \par Language:
        PyRaf
    \param data
        A list of dark files
    \param bpm
        Input bad pixel mask or NULL
    \retval 0
        If no error
    \author
        JMIbannez, IAA-CSIC
        
    """
    def __init__(self, file_list, temp_dir, output_filename="/tmp/mdark.fits", 
                 texp_scale=False, bpm=None, normalize=False):
        """Constructor"""
        
        self.__file_list = file_list
        self.__output_filename = output_filename  # full filename (path+filename)
        self.__temp_dir = temp_dir #temporal dir used for temporal/intermediate files
        self.__bpm = bpm
        self.m_min_ndarks = 3
        self.m_texp_scale = texp_scale
        self.m_normalize = normalize
        

    def createMaster(self):
      
        """
        \brief Create a master DARK from the dark file list
        """   
        log.debug("Start createMaster")
        t = utils.clock()
        t.tic()
        
        
        # Get the user-defined list of dark frames
        framelist=self.__file_list
        
        # STEP 0: Determine the number of darks frames to combine
        try:    
            nframes = len(framelist)
        except IndexError:
            log.error("No DARK frames defined")
            raise
        
        if nframes<self.m_min_ndarks:
            log.error("Not enought number of dark frames (>%s) to compute master dark: %s",self.m_min_ndarks, framelist)
            raise Exception("Not enought number of dark frames (>=%s) to compute master dark" %(self.m_min_ndarks))
            
            
        if not os.path.exists(os.path.dirname(self.__output_filename)):
            raise NameError, 'Wrong output path'
        if not self.__output_filename:
            log.error("Combined DARK frame not defined")
            raise "Wrong output filename"
    
        # Change to the source directory
        base, infile   = os.path.split(self.__output_filename)
        iraf.chdir(base)
    
        
        # STEP 1: Check the EXPTIME, TYPE(dark) of each frame
        f_expt=-1.0
        f_type=''
        f_ncoadds=-1
        f_readmode=-1
        good_frames=[]
        for iframe in framelist:
            f=datahandler.ClFits ( iframe )
            log.debug("Frame %s EXPTIME= %f TYPE= %s NCOADDS= %s REAMODE= %s" %(iframe, f.expTime(), f.getType(), f.getNcoadds(), f.getReadMode() )) 
            if not f.isDark():
                log.error("Error: Task 'createMasterDark' finished. Frame %s is not 'DARK'",iframe)
                raise Exception("Found a non DARK frame") 
                #continue
            else:        
                # Check EXPTIME, TYPE(dark) and READMODE
                if ( not self.m_texp_scale and f_expt!=-1 and 
                     (int(f.expTime()) != int(f_expt) or  
                      f.getType()!=f_type or 
                      f.getNcoadds()!=f_ncoadds or 
                      f.getReadMode()!=f_readmode)  ):
                    log.error("Error: Task 'createMasterDark' finished. Found a DARK frame (%s)with different EXPTIME, NCOADDS or READMODE",iframe)
                    #continue
                    raise Exception("Found a DARK frame with different EXPTIME or NCOADDS or READMODE") 
                else: 
                    f_expt  = f.expTime()
                    f_ncoadds= f.getNcoadds()
                    f_type  = f.getType()
                    f_readmode = f.getReadMode()
                    good_frames.append(iframe.replace("//","/"))
                                        
        log.debug('Right, dark frames with same type are: %s', good_frames)   
    
        if self.m_texp_scale:
            scale_str='exposure'
        else:
            scale_str='none'
        
        
        # Cleanup : Remove old masterdark
        misc.fileUtils.removefiles(self.__output_filename)
        tmp1 = self.__temp_dir + "/dark_tmp.fits"
        misc.fileUtils.removefiles(tmp1)
        
        #Add TEXP and NCOADD to master filename
        if f_ncoadds==-1: f_ncoadds=1
        self.__output_filename = self.__output_filename.replace(".fits","_%d_%d.fits"%(f_expt, f_ncoadds))
        
        
        # Call the noao.imred.ccdred task through PyRAF
        misc.utils.listToFile(good_frames, self.__temp_dir+"/files.list") 
        iraf.mscred.darkcombine(input = "@"+(self.__temp_dir+"/files.list").replace('//','/'),
                        output = tmp1.replace('//','/'),
                        combine = 'average',
                        ccdtype = '',
                        process = 'no',
                        reject = 'minmax',
                        nlow = '0',
                        nhigh = '1',
                        nkeep = '1',
                        scale = scale_str,
                        #expname = 'EXPTIME'
                        #ParList = _getparlistname('darkcombine')
                        )
        
    
    	if self.m_normalize:
    	    log.debug("Normalizing master dark to 1 sec")
    	    # divide master dark by the TEXP to get a master dark in ADU/s units
    	    texp=datahandler.ClFits(tmp1).expTime()
    	    iraf.mscred.mscarith(operand1 = tmp1,
    				operand2 = texp,
    				op = '/',
    				result =self.__output_filename,
    				verbose = 'no'
    				)
    	else:
    	    shutil.move(tmp1, self.__output_filename)
    	    
        darkframe = pyfits.open(self.__output_filename,'update')
        #Add a new keyword-->PAPITYPE
        darkframe[0].header.update('PAPITYPE','MASTER_DARK','TYPE of PANIC Pipeline generated file')
	darkframe[0].header.update('IMAGETYP','MASTER_DARK','TYPE of PANIC Pipeline generated file')
	if 'PAT_NEXP' in darkframe[0].header:
		darkframe[0].header.update('PAT_NEXP',1,'Number of position into the current dither pattern')
        if self.m_normalize:
            darkframe[0].header.update('EXPTIME',1.0)
            darkframe[0].header.update('ITIME', 1.0)
            darkframe[0].header.update('NCOADDS',1)
            
        darkframe.close(output_verify='ignore') # This ignore any FITS standard violation and allow write/update the FITS file    
        
        log.debug('Saved master DARK to %s' , self.__output_filename)
        log.debug("createMasterDark' finished %s", t.tac() )
        
    	# Get some stats from master dark (mean/median/rms)
        print "Stats:"
        print "----- "
        values = (iraf.mscstat (images=self.__output_filename,
                                fields="image,mean,mode,stddev,min,max",
                                format='yes',Stdout=1))
        for line in values:
    	       print line
            
        return self.__output_filename
        

################################################################################
# main
if __name__ == "__main__":
    print 'Start MasterDark....'
    # Get and check command-line options
    
    usage = "usage: %prog [options] arg1 arg2 ..."
    parser = OptionParser(usage)
    
                  
    parser.add_option("-s", "--source",
                  action="store", dest="source_file_list",
                  help="Source file list of data frames. It can be a file or directory name.")
    
    parser.add_option("-o", "--output",
                  action="store", dest="output_filename", help="final coadded output image")
    
    parser.add_option("-n", "--normalize",
                  action="store_true", dest="normalize", default=False,
                  help="normalize master dark to 1 sec [default False]")
    
    parser.add_option("-e", "--scale",
                  action="store_true", dest="texp_scale", default=False,
                  help="scale raw frames by TEXP [default False]")
    
    parser.add_option("-v", "--verbose",
                  action="store_true", dest="verbose", default=True,
                  help="verbose mode [default]")
    
    (options, args) = parser.parse_args()
    
    
    if not options.source_file_list or not options.output_filename:
	parser.print_help()
        parser.error("incorrect number of arguments " )
    
    
    filelist=[line.replace( "\n", "") for line in fileinput.input(options.source_file_list)]
    #filelist=['/disk-a/caha/panic/DATA/ALHAMBRA_1/A0408060036.fits', '/disk-a/caha/panic/DATA/ALHAMBRA_1/A0408060037.fits']
    #print "Files:",filelist
    
    try:
	mDark = MasterDark(filelist,"/tmp", options.output_filename, options.texp_scale, None, options.normalize)
	mDark.createMaster()
    except Exception,e:
	log.error("Task failed. Some error was found")
	raise e
	
    
    
        
