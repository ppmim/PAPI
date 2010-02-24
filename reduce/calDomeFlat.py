#!/usr/bin/env python
################################################################################
#
# PANICtool
#
# calDomeFlat.py
#
# Created    : 14/11/2008    jmiguel@iaa.es
# Last update: 29/09/2009    jmiguel@iaa.es
#              11/12/2009    jmiguel@iaa.es - Include the use of ClFits class, and add filter name to the output filename
#              14/12/2009    jmiguel@iaa.es - Skip non DOME flats and cotinue working with the good ones
#              12/02/2010    jmiguel@iaa.es - check min number of dome flats 
# TODO:
#    - include BPM creation
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

# Interact with FITS files
import pyfits

# Import Pyro core
#import Pyro.core
#import Pyro.naming

# Logging
from misc.paLog import log

class MasterDomeFlat:
    """
    \brief Class used to build and manage a master calibration dome flat
    
    \par Class:
        MasterDomeFlat
    \par Purpose:
        Create a master Dome flat field
    \par Description:
        
    \par Language:
        PyRaf
    \param data
        A list of dome on/off flat fields
    \param bpm
        Input bad pixel mask or NULL
    \retval 0
        If no error
    \author
        JMIbannez, IAA-CSIC
  
    """
    
    def __init__(self, input_files, output_dir, output_filename="/tmp/mdflat.fits", config_par=None):
        self.__input_files = input_files
        self.__output_file_dir = output_dir
        self.__output_filename = output_filename  # full filename (path+filename)
        self.__config_par = config_par
        self.MIN_FLATS = 4
        
    
    def createMaster(self):
      
        """
        \brief Create a master Dome FLAT from the dome flat file list
        """   
        log.debug("Start createMasterDomeFlat")
        
        start_time = time.time()
        t=utils.clock()
        t.tic()
    
        # Cleanup old files
        
        misc.fileUtils.removefiles(self.__output_filename)
        domelist_lampon  = []
        domelist_lampoff = []
        
        # Get the user-defined list of dark frames
        framelist = self.__input_files
        
        
        # Determine the number of Flats frames to combine
        try:
            nframes = len(framelist[0])
        except IndExError:
            log.error("No FLAT frames defined")
            raise 'No FLAT frames defined'
        
        
        if not os.path.exists(os.path.dirname(self.__output_filename)):
            log.error('Directory of combined FLAT frame does not exist')
            raise 'Directory of combined FLAT frame does not exist'
        if not self.__output_filename :
            log.error('Combined FLAT frame not defined')
            raise 'Combined FLAT frame not defined'
    
        
        # Change to the source directory
        base, infile   = os.path.split(self.__output_filename)
        iraf.chdir(base)
    
        m_framelist=""
        for iframe in framelist:
            m_framelist+=iframe+ ' , '
    
        # STEP 1: Check the EXPTIME , TYPE(dome) and FILTER of each Flat frame
        f_expt=-1
        f_type=''
        f_filter=''
        for iframe in framelist:
            f=datahandler.ClFits ( iframe )
            print "Flat frame %s EXPTIME= %f TYPE= %s FILTER= %s" %(iframe, f.expTime(),f.getType(), f.getFilter())
            # Check EXPTIME
            if ( f_expt!=-1 and ( int(f.expTime()) != int(f_expt) or f.getFilter()!=f_filter )) :
                log.error("Error: Task 'createMasterDomeFlat' finished. Found a FLAT frame with \n different FILTER or EXPTIME. Skipped.")
                continue
            else: 
                f_expt=f.expTime()
                if f.isDomeFlat():
                    f_filter=f.getFilter()
                else:
                    log.error("Error,  frame %s does not look a Dome Flat field" %(iframe))
                    raise Exception("Error, frame %s does not look a Dome Flat field" %(iframe))
        
            # Separate lamp ON/OFF dome flats  
            if f.isDomeFlatON():
                domelist_lampon.append(iframe)
            elif f.isDomeFlatOFF():
                domelist_lampoff.append(iframe)
            else:
                log.error("Error: Task 'createMasterDomeFlat' finished. Found a FLAT frame with different Flat-Field type (should be domeflat on/off).Skipped")
        
        
        log.info('Right, all flat frames separated as:')
        log.info('DOME FLATS LAMP OFF (#%d) %s: ' , len(domelist_lampon), domelist_lampon )
        log.info('DOME FLATS LAMP ON  (#%d) %s: ' , len(domelist_lampoff), domelist_lampoff )
        log.info('Filter=%s , TEXP=%f ' , f_filter, f_expt)
        
        if len(domelist_lampon) < self.MIN_FLATS:
            log.error("Error, not enought lamp_on flats. At least %s are requered" %(self.MIN_FLATS))
            raise Exception("Error, not enought lamp_on flats. At least %s are requered" %(self.MIN_FLATS))
        
        if len(domelist_lampoff) < self.MIN_FLATS:
            log.error("Error, not enought lamp_on flats. At least %s are requered"%(self.MIN_FLATS))
            raise Exception("Error, not enought lamp_off flats. At least %s are requered" %(self.MIN_FLATS))
    
        #Clobber existing output images
        iraf.clobber='yes'
        
        # NOTE : We do not subtract any MASTER_DARK, it is not required for DOME FLATS (it is done implicitly)
    
        # STEP 2: Make the combine of Flat LAMP-ON frames
        # - Build the frame list for IRAF
        log.debug("Combining Flat LAMP-ON frames...")
        m_framelist=''
        for iframe in domelist_lampon:
            m_framelist+=iframe+ ' , '
        flat_lampon=self.__output_file_dir + "/flat_lampON.fits"
        misc.fileUtils.removefiles(flat_lampon)
        # - Call IRAF task
        iraf.flatcombine(input=m_framelist,
                        output=flat_lampon,
                        combine='median',
                        ccdtype='none',
                        process='no',
                        reject='sigclip',
                        subset='yes',
                        scale='mode'
                        #verbose='yes'
                        #scale='exposure',
                        #expname='EXPTIME'
                        #ParList = _getparlistname ('flatcombine')
                        )
    
        # STEP 3: Make the combine of Flat LAMP-OFF frames
        # - Build the frame list for IRAF
        log.debug("Combining Flat LAMP-OFF frames...")    
        m_framelist=''
        for iframe in domelist_lampoff:
            m_framelist+=iframe+ ' , '
        flat_lampoff=self.__output_file_dir + "/flat_lampOFF.fits"
        misc.fileUtils.removefiles(flat_lampoff)
        # - Call IRAF task
        iraf.flatcombine(input=m_framelist,
                        output=flat_lampoff,
                        combine='median',
                        ccdtype='none',
                        process='no',
                        reject='sigclip',
                        subset='yes',
                        scale='mode',
                        #verbose='yes'
                        #scale='exposure',
                        #expname='EXPTIME'
                        #ParList = _getparlistname ('flatcombine')
                        )
    
        # STEP 4 : Subtract lampON-lampOFF (implicit dark subtraction)
        flat_diff=self.__output_file_dir+"/flat_lampON_OFF.fits"
        log.debug("Subtractic Flat ON-OFF frames...") 
        # Remove an old masternormflat
        misc.fileUtils.removefiles(flat_diff)
        iraf.imarith(operand1 = flat_lampon,
                    operand2 = flat_lampoff,
                    op = '-',
                    result =flat_diff,
                    verbose = 'yes'
                    )
            
        # STEP 5: Normalize the flat-field
        # Compute the mean of the image
        log.debug("Normalizing master flat frame...")
        mean = float(iraf.imstat (
            images=flat_diff,
            fields='mean',Stdout=1)[1])
        
        # Cleanup: Remove temporary files
        misc.fileUtils.removefiles(flat_lampoff, flat_lampon)
        # Compute normalized flat
        self.__output_filename=self.__output_filename.replace(".fits","_%s.fits"%(f_filter))
        misc.fileUtils.removefiles(self.__output_filename)
        iraf.imarith(operand1=flat_diff,
                    operand2=mean,
                    op='/',
                    result=self.__output_filename,
                    )
    
        
        # Change back to the original working directory
        iraf.chdir()
        
        flatframe = pyfits.open(self.__output_filename,'update')
        flatframe[0].header.add_history('Computed normalized master dome flat (lamp_on-lamp_off)' )
        flatframe[0].header.add_history('lamp_on  files: %s' %domelist_lampon )
        flatframe[0].header.add_history('lamp_off files: %s' %domelist_lampoff )
        #Add a new keyword-->PIP_TYPE
        flatframe[0].header.update('PATYPE','MASTER_DOME_FLAT','TYPE of PANIC Pipeline generated file')
        flatframe[0].header.update('OBJECT','MASTER_DOME_FLAT')
        flatframe.close(output_verify='ignore') # This ignore any FITS standar violation and allow write/update the FITS file
        
        log.debug(t.tac() )
        log.debug('Saved master FLAT to %s' ,  self.__output_filename )
    
        return self.__output_filename
        
################################################################################
# Functions       
def usage ():
    print "Create 'master dome flat' procedure:"
    print "Unknown command line parameter. Required parameters are : "
    print "-s / --source=      Source file list of data frames"
    print "-o / --out=         Output master filename "

################################################################################
# main
if __name__ == "__main__":
    print 'Start MasterFlat....'
    # Get and check command-line options
    args = sys.argv[1:]
    source_file_list = ""
    output_filename = ""
    
    try:
        opts, args = getopt.getopt(args, "s:o:", ['source=','out='])
    except getopt.GetoptError:
        # print help information and exit:
        usage()
        sys.exit(1)

    
    for option, parameter in opts:
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
    #output_filename="/tmp/out/out.fits"
    #filelist=['/disk-a/caha/panic/DATA/ALHAMBRA_1/A0408060036.fits', '/disk-a/caha/panic/DATA/ALHAMBRA_1/A0408060037.fits']

    print "Files:",filelist
    mDFlat = MasterDomeFlat(filelist, "/tmp", output_filename)
    mDFlat.createMaster()
    
        
