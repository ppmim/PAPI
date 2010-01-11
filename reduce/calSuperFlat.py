#!/usr/bin/env python
################################################################################
#
# PANICtool
#
# calDomeFlat.py
#
# Created    : 14/11/2008    jmiguel@iaa.es
# Last update: 30/09/2009    jmiguel@iaa.es
#              16/12/2009    jmiguel@iaa.es - Use of ClFits; check NCOADDS; added --force and -n options
#                                           - modify dark subtration, now is done after combine (much efficient)
#
# TODO:
#  - work with dark models instead of a fixed master dark
#  - work  with diff EXPTIMEs
#  - create a BPM
#  - avoid frames with low sky background
#  - implement combine in memory using numpy (without IRAF.imcombine)
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
import Pyro.core
import Pyro.naming

# Logging
from misc.paLog import log

class MasterSuperFlat:
    """
    \brief Class used to build and manage a master calibration super flat from a list of night SCIENCE frames
    \par Class:
        MasterSuperFlat
    \par Purpose:
        Create a master normalized night Super Flat Field
    \par Description:
        
        1. Check the TYPE(science) and FILTER and EXPTIME of each frame
           If any frame on list missmatch the FILTER,TYPE,EXPTIME,NCOADDS then the master super-flat will be aborted
        
        2. Make the combine of Flat frames scaling by 'mode'
        
        3. If required, we subtract a proper MASTER_DARK using the provided master_dark file (optional)

        4. If required, normalize the super-flat dividing by the mean value (default True)
        
    \par Language:
        PyRaf
    \param data
        A list of science frames
    \param bpm
        Input bad pixel mask or NULL
    \param mdark
        Master dark to subtract (optional)
    \retval median
        When all goes well
    \retval 0
        If no error
    \author
        JMIbannez, IAA-CSIC
    """
    
    def __init__(self, input_files, output_filaname="/tmp/msflat.fits", master_dark=None, output_dir="/tmp", normal=True, force=False):
        """
        \brief Init the object
        \param input_file - data science files, non dark subtracted
        \param output_filename - output master super flat produced 
        \param master_dark - master dark to be subtracted to each science frame
        \param working_dir - dir used for computations and intermediate files
        """
        
        self.__input_files = input_files
        self.__output_file_dir = output_dir
        self.__master_dark = master_dark
        self.__output_filename = output_filename  # full filename (path+filename)
        
        self.m_min_flats = 5
        self.m_normaliz = normal
        self.m_force = force
        
    def createMaster(self):
      
        """
        \brief Create a master SUPER FLAT from a list for science files
        """   
        log.debug("Start createMasterSuperFlat")
        start_time = time.time()
        t=utils.clock()
        t.tic()
    
        # Cleanup old files
        misc.fileUtils.removefiles(self.__output_filename)
        
        # Get the user-defined list of dark frames
        framelist = self.__input_files
        
        
        # Determine the number of frames to combine
        try:
            nframes = len(framelist[0])
        except IndExError:
            raise ExError('No frames found')
        
        if nframes<self.m_min_flats:
            log.error("Not enought number of flat frames (>%s) to compute master super-flat",self.m_min_flats)
            return False
        
        if not os.path.exists(os.path.dirname(self.__output_filename)):
            raise ExError, 'Directory of combined FLAT frame does not exist'
        if not self.__output_filename :
            raise ExError('Combined FLAT frame not defined')
    
        # Change to the source directory
        base, infile   = os.path.split(self.__output_filename)
        iraf.chdir(base)
    
        m_framelist=""
        for iframe in framelist:
            m_framelist+=iframe+ ' , '
        #print "LIST of Flat frames: ", m_framelist    
    
        # STEP 1: Check the EXPTIME , TYPE(science) and FILTER of each input science frame
        f_expt=-1
        f_type=''
        f_filter=''
        f_ncoadds=-1
        for iframe in framelist:
            f=datahandler.ClFits ( iframe )
            print "Flat frame %s EXPTIME= %f TYPE= %s FILTER= %s" %(iframe, f.expTime(),f.getType(), f.getFilter())
        
            # Check EXPTIME and FILTER and TYPE (object)
            if ( f_expt!=-1 and ( int(f.expTime()) != int(f_expt) or f.getFilter()!=f_filter or f.getType()!=f_type or f.getNcoadds()!=f_ncoadds) ) :
                log.error("Error: Task 'createMasterSkyFlat' finished. Found a frame with \n different FILTER or EXPTIME or NCOADDS or TYPE")
                raise Execption("Error, found a frame with different TYPE, FILTER or EXPTIME or NCOADDS")
            else:
                f_expt=f.expTime()
                f_filter=f.getFilter()
                f_type = f.getType()
                f_ncoadds=f.getNcoadds()
        
        log.info('OK, all  frames same type:')
        log.info('Filter=%s , TEXP=%f TYPE=%s' , f_filter, f_expt, f_type)
        
    
        #Clobber existing output images
        iraf.clobber='yes'
        
        # STEP 2: Make the combine of dark subtracted Flat frames scaling by 'mode'
        # - Build the frame list for IRAF
        log.debug("Combining sky science frames...")
        comb_flat_frame=(self.__output_file_dir+"/comb_sp_flats.fits").replace("//","/")
        misc.fileUtils.removefiles(comb_flat_frame)
        i_darksublist=utils.listToString(framelist)
        # - Call IRAF task
        iraf.imcombine(input=i_darksublist,
                        output=comb_flat_frame,
                        combine='median',
                        offset='none',
                        reject='sigclip',
                        lsigma=2.5,
                        hsigma=2.5,
                        scale='median',
                        zero='none'
                        #verbose='yes'
                        #scale='exposure',
                        #expname='EXPTIME'
                        #ParList = _getparlistname ('flatcombine')
                        )
                        
        # STEP 3 : If required, we subtract a valid MASTER_DARK (same EXPTIME and NCOADDS) 
        if self.__master_dark!=None:
            #check master dark (here only check TYPE and EXPTIME and NCOADDS) 
            fdark=datahandler.ClFits ( master_dark )
            texp_dark=fdark.expTime()
            if not fdark.isMasterDark() or (not self.m_force & (fdark.getNcoadds()!=f_ncoadds or texp_dark!=f_expt)):
                log.error("Error, wrong master dark provided: %s. \nCheck TYPE or NCOADDS or EXPTIME", master_dark)
                raise Exception("Error, wrong master dark provided -->%s" %(master_dark))
            else:
                ds_superflat=comb_flat_frame.replace(".fits","_D.fits")
                # Remove any old dark subtracted flat frames
                misc.fileUtils.removefiles(ds_superflat)
                
                # Build master dark with proper EXPTIME and subtract (???? I don't know how well is this method !!!)
                f = pyfits.open(comb_flat_frame)
                mdark = pyfits.open(self.__master_dark)
                if f_expt!=texp_dark: #only for information
                    log.warning("Warning, EXPTIME does not match. Scaling dark before subtraction.")                  
                #pr_mdark = (numpy.array(mdark[0].data, dtype=numpy.double)/float(mdark[0].header['EXPTIME']))*float(f[0].header['EXPTIME'])
                f[0].data = f[0].data - mdark[0].data*float(f_expt/texp_dark)
                f[0].header.add_history('Dark subtracted %s' %self.__master_dark)
                   
                #a=numpy.reshape(f[0].data, (2048*2048,))
                #print "MODE=", 3*numpy.median(a)-2*numpy.mean(a)
                #print "MEAN=" , numpy.mean(a)
                
                # Write output to outframe (data object actually still points to input data)
                try:
                    f.writeto(ds_superflat, output_verify='ignore')
                    f.close()
                    last_out=ds_superflat
                except IOError:
                    raise ExError('Cannot write output to %s' % ds_superflat)
        else:
            last_out=comb_flat_frame    
        
        if self.m_normaliz:
            # STEP 4: Normalize the flat-field
            # Compute the mean of the image
            log.debug("Normalizing master flat frame...")
            median = float(iraf.imstat (
                images=("'"+last_out+"[100:1900,100:1900]'").replace('//','/'),
                fields='midpt',format='no',Stdout=1)[0])
            
            # Cleanup: Remove temporary files
            misc.fileUtils.removefiles(self.__output_filename)
            # Compute normalized flat
            iraf.imarith(operand1=last_out,
                    operand2=median,
                    op='/',
                    pixtype='real',
                    result=self.__output_filename,
                    )
        else:
            #only rename the last output
            shutil.move(last_out, self.__output_filename)
                
        # Change back to the original working directory
        iraf.chdir()
        
        flatframe = pyfits.open(self.__output_filename,'update')
        flatframe[0].header.add_history('Computed normalized master twilight flat' )
        flatframe[0].header.add_history('Twilight files: %s' %framelist )
        #Add a new keyword-->PIP_TYPE
        flatframe[0].header.update('PIP_TYPE','MASTER_TW_FLAT','TYPE of PANIC Pipeline generated file')
        flatframe[0].header.update('OBJECT','MASTER_TW_FLAT')
        flatframe.close(output_verify='ignore') # This ignore any FITS standar violation and allow write/update the FITS file
        
        log.debug(t.tac() )
        log.debug('Saved master SUPER_FLAT to %s' ,  self.__output_filename )
    
        return self.__output_filename
        
################################################################################
# Functions       
def usage ():
    print "Create 'master super flat' frame:"
    print "Required parameters are : "
    print "-s / --source=      Source file list of data frames"
    print "-o / --out=         Output master filename "
    print "-d / --dark=        Master dark frame (optional)"
    print "-n /                Do not normalize the output super flat (optional)"
    print "-f / --force        Force dark subtraction, even EXPTIME or NCOADDS does not match (optional, default=False)"

################################################################################
# main
if __name__ == "__main__":
    print 'Start MasterSuperFlat....\n'
    # Get and check command-line options
    args = sys.argv[1:]
    source_file_list = ""
    output_filename = ""
    master_dark=None # default
    normaliz=True
    force=False
    
    try:
        opts, args = getopt.getopt(args, "s:o:d:nf", ['source=', 'out=','dark='])
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
        if option in ("-d", "--dark"):
            master_dark = parameter
            print "Master Dark =", master_dark
        if option in ("-n"):
            normaliz=False    
        if option in ("-f", "--force"):
            force=True
                
    if  source_file_list=="" or output_filename=="":
        usage()
        sys.exit(3)
    
    filelist=[line.replace( "\n", "") for line in fileinput.input(source_file_list)]
    #filelist=['/disk-a/caha/panic/DATA/ALHAMBRA_1/A0408060036.fits', '/disk-a/caha/panic/DATA/ALHAMBRA_1/A0408060037.fits']
    mSFlat = MasterSuperFlat(filelist,output_filename, master_dark, "/tmp", normaliz, force)
    mSFlat.createMaster()
    
        