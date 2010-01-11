#
# PANICtool
#
# coadd.py
#
# Created    : 02/12/2009    jmiguel@iaa.es
# Last update: 02/12/2009    jmiguel@iaa.es - 
#
# TODO
#  - 
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

# Pyraf modules
from pyraf import iraf
from iraf import noao
from iraf import imred
from iraf import ccdred

# Interact with FITS files
import pyfits

# Import Pyro core
import Pyro.core
import Pyro.naming

# Logging
from misc.paLog import log

class Coadd:
    """
    \brief Class used to coadd a set of dithered data
    
    \par Class:
        Coadd
    \par Purpose:
        Coadd a set of dithered data files using SWARP
    \par Description:
            
    \par Language:
        PyRaf
    \param data
        A list of calibrated files (dark subtracted, flatfield, sky subtracted and re-gridded )
    \retval 0
        If no error
    \author
        JMIbannez, IAA-CSIC (jmiguel@iaa.es)
        
    """
    def __init__(self, file_list, output_dir, output_filename="/tmp/coadd.fits"):
        self.__file_list=file_list
        self.__output_file_dir=output_dir
        self.__output_filename=output_filename  # full filename (path+filename)
    
    def coadd(self):
      
        """
        \brief Coadd data set 
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
            log.error("Not enought number of dark frames (>%s) to compute master dark",self.m_min_ndarks)
            raise Exception("Not enought number of dark frames")
            
            
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
        for iframe in framelist:
            f=pyfits.open(iframe)
            log.debug("Frame %s EXPTIME= %f TYPE= %s " %(iframe, f[0].header['EXPTIME'],f[0].header['OBJECT'])) 
            if self.m_texp_scale==True:
                #Only check TYPE (dark)
                f_type  =f[0].header['OBJECT']
                if not f_type.count('dark'):
                    log.error("Error: Task 'createMasterDark' finished. Frame type is not 'DARK'.")
                    f.close()
                    raise SystemExit
            else:        
                # Or Check EXPTIME, TYPE (dark)
                if ( f_expt!=-1 and (int(f[0].header['EXPTIME']) != int(f_expt) or  f[0].header['OBJECT']!=f_type )  ):
                    log.error("Error: Task 'createMasterDark' finished. Found a DARK frame with different EXPTIME")
                    f.close()
                    raise Exception("Found a DARK frame with different EXPTIME") 
                else:
                    f_expt  =f[0].header['EXPTIME']
                    f_type  =f[0].header['OBJECT']
                    if not f_type.count('dark'):
                        log.error("Error: Task 'createMasterDark' finished. Frame type is not 'DARK'.")
                        f.close()
                        raise Exception("Found a frame not being a DARK frame")
                
        log.debug('Right, all frames are same type')   
    
        m_framelist=""
        for iframe in framelist:
            m_framelist+=iframe+ ' , '
            
        # Cleanup : Remove an old masterdark
        misc.fileUtils.removefiles(self.__output_filename)
        if self.m_texp_scale:
            scale_str='exposure'
        else:
            scale_str='none'
        # Call the noao.imred.ccdred task through PyRAF
        iraf.darkcombine(input=m_framelist,
                        output=self.__output_filename,
                        combine='average',
                        ccdtype='none',
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
        darkframe[0].header.add_history('Combined images by averaging (%s files) ' % m_framelist)
        #Add a new keyword-->PIP_TYPE
        darkframe[0].header.update('hierarch PAPI.TYPE','MASTER_DARK','TYPE of PANIC Pipeline generated file')
        darkframe[0].header.update('OBJECT','MASTER_DARK')
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