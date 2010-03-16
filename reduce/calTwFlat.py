#!/usr/bin/env python
################################################################################
#
# PANICtool
#
# calTwilightFlat.py
#
# Created    : 19/05/2009    jmiguel@iaa.es
# Last update: 29/09/2009    jmiguel@iaa.es
#              14/12/2009    jmiguel@iaa.es  - Check NCOADDS; use ClFits class; Skip non TW flats and cotinue working with the good ones
#              03/03/2010    jmiguel@iaa.es Added READMODE checking 
#
# TODO:
#   - use of dark model to subtract right dark current
#   - take into account BPM !!!
#   - compute automatically the level of counts the twilight flats should have (lthr,hthr)
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

import numpy

# Interact with FITS files
import pyfits

# Import Pyro core
import Pyro.core
import Pyro.naming

# Logging
from misc.paLog import log

class MasterTwilightFlat:
    """
    \brief Class used to build and manage a master calibration twilight flat
    \par Class:
        MasterTwilightFlat
    \par Purpose:
        Create a master Twilight flat field
    \par Description:
        
        1. Check the  TYPE(twilight) and FILTER of each Flat frame
           If any frame on list missmatch the FILTER, then the master twflat will skip this frame and contiune with then next ones
           EXPTIME do not need be the same, so EXPTIME scaling with 'mode' will be done
           
           1.1: Check either over or under exposed frames
        
        2. We subtract a proper MASTER_DARK, it is required for TWILIGHT FLATS because they might have diff EXPTIMEs
        
        3. Make the combine (with sigclip rejection) of dark subtracted Flat frames scaling by 'mode'
        
        4. Normalize the tw-flat dividing by the mean value
        
    \par Language:
        PyRaf
    \param data
        A list file of twilight flat fields filenames
        Input bad pixel mask or NULL
    \param mdark
        Master current dark model to dark subtraction (required)
    \param output_filename  
        Master Tw Flat created 
    \param lthr
        Low threshold to identify good twilight flats (default 100)
    \param hthr
        High threshold to identify good twilight flats (default 100000)
    \param bpm
        Bad Pixel mask to use (optional)
    \retval median
        When all goes well
    \retval 0
        If no error
    \author
        JMIbannez, IAA-CSIC
  
    """
    def __init__(self, flat_files, dark_model, output_filename="/tmp/mtwflat.fits", lthr=100, hthr=100000, bpm=None):
        self.__input_files = flat_files
        self.__master_dark = dark_model
        self.__output_file_dir = os.path.dirname(output_filename)
        self.__output_filename = output_filename  # full filename (path+filename)
        self.__bpm = bpm
        
        self.m_MIN_N_GOOD=2
        self.m_lthr=lthr
        self.m_hthr=hthr
        self.m_min_flats=5
    
    def createMaster(self):
      
        """
        \brief Create a master Tw FLAT from the flat file list
        """   
        log.debug("Start createMasterTwilightFlat")
        
        start_time = time.time()
        t=utils.clock()
        t.tic()
    
        # Cleanup old files
        
        misc.fileUtils.removefiles(self.__output_filename)
        
        # Get the user-defined list of flat frames
        framelist = self.__input_files
        
        
        # Determine the number of Flats frames to combine
        try:
            nframes = len(framelist[0])
        except IndExError:
            log.error("No FLAT frames defined")
            raise ExError('No FLAT frames defined')
        
        if nframes<self.m_min_flats:
            log.error("Not enought number of flat frames (>%s) to compute master tw-flat",self.m_min_flats)
            return False
        
        if not os.path.exists(os.path.dirname(self.__output_filename)):
            log.error('Directory of combined FLAT frame does not exist')
            raise ExError, 'Directory of combined FLAT frame does not exist'
        if not self.__output_filename :
            log.error('Combined FLAT frame not defined')
            raise ExError('Combined FLAT frame not defined')
    
        
        # Change to the source directory
        base, infile   = os.path.split(self.__output_filename)
        iraf.chdir(base)
    
        
        # STEP 1: Check the  TYPE(twilight) and FILTER,READEMODE of each Flat frame
        # If any frame on list missmatch the FILTER, then the master twflat will be aborted
        # EXPTIME do not need be the same, so EXPTIME scaling will be done
        f_expt=-1
        f_type=''
        f_filter=''
        f_ncoadds=-1
        f_readmode=''
        good_frames=[]
        
        for iframe in framelist:
            f=datahandler.ClFits ( iframe )
            print "Flat frame %s EXPTIME= %f TYPE= %s FILTER= %s" %(iframe, f.expTime(),f.getType(), f.getFilter())
            if ( f_expt!=-1 and (f.getFilter()!=f_filter or f.getType()!=f_type or f.getReadMode()!=f_readmode)) :
                log.error("Task 'createMasterTwFlat' Found a FLAT frame with different FILTER or TYPE.")
                raise Exception("Error, frame %s has different FILTER or TYPE" %(iframe))
                #continue
            else: 
                f_expt=f.expTime()
                f_filter=f.getFilter()
                f_readmode=f.getReadMode()
                if f.isTwFlat():
                    f_type=f.getType()
                else:
                    log.error("Error, frame %s does not look a TwiLight Flat field" %(iframe))
                    raise Exception("Error, frame %s does not look a TwiLight Flat field" %(iframe))
            
            # STEP 1.1: Check either over or under exposed frames
            mean = float(iraf.imstat (
                            images=iframe+"[300:700,300:700]",
                            fields='mean',Stdout=1)[1])
            print "File %s filter[ %s ]  EXTP=%f TYPE=%s mean_window=%f" %(iframe, f_filter, f_expt, f_type, mean)
            if mean>self.m_lthr and mean<self.m_hthr:
                good_frames.append(iframe)
            else:
                log.error("Frame %s skipped, either over or under exposed" %(iframe))
            
            
        if len(good_frames)>self.m_MIN_N_GOOD:
            log.info('Found %d flat frames with same filter [%s] and type:\n', len(good_frames), f_filter)
            for e in good_frames:
                log.info("--->%s",e)            
        else:
            log.error("Error, not enought good frames, exiting....")
            raise Exception("Error, not enought good flat frames")
                
        #Clobber existing output images
        iraf.clobber='yes'
        
        # STEP 2: We subtract a proper MASTER_DARK, it is required for TWILIGHT FLATS because they might have diff EXPTIMEs
        # Prepare input list on IRAF string format
        m_framelist=utils.listToString(good_frames)
            
        # Prepare dark subtracted flat list
        m_darksublist=m_framelist.replace(".fits", "_D.fits")
         
        
        mdark = pyfits.open(self.__master_dark)
        t_dark=datahandler.ClFits ( self.__master_dark ).expTime()
        for iframe in good_frames:
            # Remove an old dark subtracted flat frames
            misc.fileUtils.removefiles(iframe.replace(".fits","_D.fits"))
            
            # Build master dark with proper EXPTIME and subtract (???? I don't know how good is this method !!!)
            f = pyfits.open(iframe)
            t_flat=datahandler.ClFits ( iframe ).expTime()
            #pr_mdark = (numpy.array(mdark[0].data, dtype=numpy.double)/float(mdark[0].header['EXPTIME']))*float(f[0].header['EXPTIME'])
            f[0].data = f[0].data - mdark[0].data*float(t_flat/t_dark)
            f[0].header.add_history('Dark subtracted %s (interpolated)' %self.__master_dark)
            #a=numpy.reshape(f[0].data, (2048*2048,))
            #print "MODE=", 3*numpy.median(a)-2*numpy.mean(a)
            #print "MEAN=" , numpy.mean(a)
            
            # Write output to outframe (data object actually still points to input data)
            try:
                f.writeto(iframe.replace(".fits","_D.fits"), output_verify='ignore')
            except IOError:
                raise ExError('Cannot write output to %s' % iframe.replace(".fits", "_D.fits"))
                     
            f.close()
                    
        
        # STEP 3: Make the combine of dark subtracted Flat frames scaling by 'mode'
        # - Build the frame list for IRAF
        log.debug("Combining Twilight Flat frames...")
        comb_flat_frame=(self.__output_file_dir+"/comb_tw_flats.fits").replace("//","/")
        print "COMB=", comb_flat_frame
        misc.fileUtils.removefiles(comb_flat_frame)
        # - Call IRAF task
        iraf.flatcombine(input=m_darksublist,
                        output=comb_flat_frame,
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
        
        # STEP 4: Normalize the flat-field
        # Compute the mean of the image
        log.debug("Normalizing master flat frame...")
        median = float(iraf.imstat (
            images=comb_flat_frame,
            fields='midpt',Stdout=1)[1])
        
        # Cleanup: Remove temporary files
        misc.fileUtils.removefiles(self.__output_filename)
        # Compute normalized flat
        iraf.imarith(operand1=comb_flat_frame,
                    operand2=median,
                    op='/',
                    pixtype='real',
                    result=self.__output_filename,
                    )
    
        
        # Change back to the original working directory
        iraf.chdir()
        
        flatframe = pyfits.open(self.__output_filename,'update')
        flatframe[0].header.add_history('Computed normalized master twilight flat' )
        flatframe[0].header.add_history('Twilight files: %s' %framelist )
        #Add a new keyword-->PIP_TYPE
        flatframe[0].header.update('PIP_TYPE','MASTER_TW_FLAT','TYPE of PANIC Pipeline generated file')
        flatframe[0].header.update('OBJECT','MASTER_TW_FLAT','Master Twilight flat field')
        flatframe.close(output_verify='ignore') # This ignore any FITS standar violation and allow write/update the FITS file
        
        log.debug(t.tac() )
        log.debug('Saved master TW_FLAT to %s' ,  self.__output_filename )
    
        return self.__output_filename
        
        
        
        
################################################################################
# Functions       
def usage ():
    print "Create 'master twilight flat' procedure:"
    print "Unknown command line parameter. Required parameters are : "
    print "-s / --source=      Source file list of data frames"
    print "-d / --dark=        Master dark current model to dark subtraction"
    print "-o / --out=         Output master filename "

################################################################################
# main
if __name__ == "__main__":
    print 'Start MasterTwFlat....'
    # Get and check command-line options
    args = sys.argv[1:]
    source_file_list = ""
    output_filename = ""
    dark_file =""
    try:
        opts, args = getopt.getopt(args, "s:d:o:", ['source=','dark=', 'out='])
    except getopt.GetoptError:
        # print help information and exit:
        usage()
        sys.exit(1)

    
    for option, parameter in opts:
        if option in ("-s", "--source"):
            source_file_list = parameter
            print "Source file list =", source_file_list
            if not os.path.isfile(source_file_list):
                print 'Error, file list does not exists :', source_file_list
                sys.exit(1)
        if option in ("-d", "--dark"):
            dark_file = parameter
            print "Dark file =", dark_file
        if option in ("-o", "--out"):
            output_filename = parameter
            print "Output file =", output_filename
            
    if  source_file_list=="" or output_filename=="" or dark_file=="":
        usage()
        sys.exit(3)
    
    filelist=[line.replace( "\n", "") for line in fileinput.input(source_file_list)]
    #filelist=['/disk-a/caha/panic/DATA/ALHAMBRA_1/A0408060036.fits', '/disk-a/caha/panic/DATA/ALHAMBRA_1/A0408060037.fits']
    print "Files:",filelist
    mTwFlat = MasterTwilightFlat(filelist,dark_file, output_filename)
    mTwFlat.createMaster()
    
        