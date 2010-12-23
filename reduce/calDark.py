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
    def __init__(self, file_list, output_dir, output_filename="/tmp/mdark.fits", texp_scale=False, bpm=None):
        self.__file_list=file_list
        self.__output_file_dir=output_dir
        self.__output_filename=output_filename  # full filename (path+filename)
        self.__bpm=bpm
        self.m_min_ndarks = 3
        self.m_texp_scale = texp_scale
    
    def createMaster(self):
      
        """
        \brief Create a master DARK from the dark file list
        """   
        log.debug("Start createMaster")
        start_time = time.time()
        t=utils.clock()
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
                log.error("Error: Task 'createMasterDark' finished. Frame type is not 'DARK'.")
                raise Exception("Found a non DARK frame") 
                #continue
            else:        
                # Check EXPTIME, TYPE(dark) and READMODE
                if ( not self.m_texp_scale and f_expt!=-1 and (int(f.expTime()) != int(f_expt) or  f.getType()!=f_type or f.getNcoadds()!=f_ncoadds or f.getReadMode()!=f_readmode)  ):
                    log.error("Error: Task 'createMasterDark' finished. Found a DARK frame with different EXPTIME, NCOADDS or READMODE")
                    #continue
                    raise Exception("Found a DARK frame with different EXPTIME or NCOADDS or READMODE") 
                else: 
                    f_expt  = f.expTime()
                    f_ncoadds= f.getNcoadds()
                    f_type  = f.getType()
                    f_readmode = f.getReadMode()
                    good_frames.append(iframe.replace("//","/"))
                                        
        log.debug('Right, dark frames with same type are: %s', good_frames)   
    
        # Cleanup : Remove old masterdark
        misc.fileUtils.removefiles(self.__output_filename)
        if self.m_texp_scale:
            scale_str='exposure'
        else:
            scale_str='none'
        
        if f_ncoadds==-1: f_ncoadds=1
        self.__output_filename=self.__output_filename.replace(".fits","_%d_%d.fits"%(f_expt,f_ncoadds))
        
        # Call the noao.imred.ccdred task through PyRAF
        
        misc.utils.listToFile(good_frames, self.__output_file_dir+"/files.list") 
        iraf.mscred.darkcombine(input="@"+(self.__output_file_dir+"/files.list").replace('//','/'),
                        output=self.__output_filename,
                        combine='average',
                        ccdtype='',
                        process='no',
                        reject='minmax',
                        nlow='0',
                        nhigh='1',
			nkeep='1',
                        scale=scale_str,
                        #expname='EXPTIME'
                        #ParList = _getparlistname('darkcombine')
                        )
                        
        #outdata[0].header.add_history('Averaged %i frames to obtain combined DARK' % nframes)
        
        darkframe = pyfits.open(self.__output_filename,'update')
        #darkframe[0].header.add_history('Combined images by averaging (%s files) ' % good_frames)
        #Add a new keyword-->PAPITYPE
        darkframe[0].header.update('PAPITYPE','MASTER_DARK','TYPE of PANIC Pipeline generated file')
        #darkframe[0].header.update('OBJECT','MASTER_DARK')
        darkframe.close(output_verify='ignore') # This ignore any FITS standar violation and allow write/update the FITS file    
    
        log.debug('Saved master DARK to %s' , self.__output_filename)
        log.debug("createMasterDark' finished %s", t.tac() )
        
        return self.__output_filename
        
################################################################################
# Functions       
def usage ():
    print "The required parameters are : "
    print "-s / --source=     Source file list of data frames"
    print "-o / --out=        Output master filename "
    print "\n"
    print "Optional parameters :"
    print "-t                 Scale by TEXP "

################################################################################
# main
if __name__ == "__main__":
    print 'Start MasterDark....'
    # Get and check command-line options
    args = sys.argv[1:]
    source_file_list = ""
    output_filename = ""
    texp_scale = False
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "s:o:t", ["source=","out=","t"])
    except getopt.GetoptError:
        # print help information and exit:
        usage()
        sys.exit(1)

    
    for option, parameter in opts:
        if option in ("-s", "--source"):
            source_file_list = parameter
            if not os.path.exists(os.path.dirname(source_file_list)):
                print 'Error, file list does not exists'
                sys.exit(1)
        if option in ("-o", "--out"):
            output_filename = parameter
            print "Output file =", output_filename
        if option in ("-t"):
            texp_scale = True
            print "TEXP scale =", texp_scale
            
    if  source_file_list=="" or output_filename=="":
        usage()
        sys.exit(3)
    
    filelist=[line.replace( "\n", "") for line in fileinput.input(source_file_list)]
    #filelist=['/disk-a/caha/panic/DATA/ALHAMBRA_1/A0408060036.fits', '/disk-a/caha/panic/DATA/ALHAMBRA_1/A0408060037.fits']
    #print "Files:",filelist
    mDark = MasterDark(filelist,"/tmp", output_filename, texp_scale)
    mDark.createMaster()
    
        
