#!/usr/bin/env python

################################################################################
#
# SimpleReduce and ReductionBlock classes
#
# reducemod.py
#
# Last update 30/04/2008
#
################################################################################


"""
   Main routines for the simple automatic data reduction
"""

_version     = "1.0.0"
_date        = "30-04-2008"
_author      = "Jose M. Ibanez (jmiguel@iaa.es)"


_minversion_numpy   = "1.0.1"
_minversion_pyfits  = "1.1"
_minversion_pyraf   = "1.4"
_minversion_biggles = "1.6.4"


################################################################################
# Import necessary modules

import getopt
import sys
import os
import shutil
import datetime
import subprocess
import glob


import misc.fileUtils
import misc.utils as utils
import datahandler

# From reduce...
import calDark
import calDomeFlat
from makeobjmask import *
from mksuperflat import *
from imtrim import *

import numpy
#import scipy
#import scipy.ndimage
#import scipy.ndimage.filters
#import scipy.stats.mstats


# Pyraf modules
from pyraf import iraf
from iraf import noao
from iraf import imred
from iraf import ccdred
from iraf import images
from iraf import tv
#from iraf import mscred OJO !!!! puede hacer que no funcione el ccdred normal (combine, ...)

# Interact with FITS files
import pyfits

# Import Pyro core
#import Pyro.core
#import Pyro.naming

# Logging
import misc.paLog
from misc.paLog import log


################################################################################
class SimpleReduce:
 
    """
    \brief
        Do a simple data reduction to a single frame
        Subtract sky from an image using only its own background for sky computing
            
    \par Description:
         
         Algorith to do a simple and fast reduction of an IR image
         -------------------------- 
            
    \par Language:
        Python
    \param data
        
    \retval 0
        If no error, the reduced image
    \author
        JMIbannez, IAA-CSIC
    """
    
    def __init__(self):
        log.debug("SimpleReduce Object created !!")
        run_status=0
        pass
        
    def run (self, frame_in, master_dark, master_flat, result_frame="/tmp/result.fits", out_tmp_dir='/tmp/', appPixMask=False):
    
        run_status = 0
        frame_out = ''
        
        if  frame_in=="" or master_dark=="" or master_flat=="" or result_frame=="":
            raise "Error, wrong input frames entered "
        
        
        t=misc.utils.clock()
        t.tic()
    
        log.debug("[SimpleReduce]:Running simple reduce proceduce...")
        
        frameoutD = frame_in.replace(".fits","_D.fits")
        frameoutD = frameoutD.replace(os.path.dirname(frameoutD), out_tmp_dir)
        frameoutD = frameoutD.replace("//","/")
        
        fc=datahandler.ClFits(frame_in)
        ff=datahandler.ClFits(master_flat)
        if ( fc.filter!=ff.filter ):
            log.error("Error reducing frame %s. Filter missmatch with calibration frames" %frame_in)
            raise "Error reducing frame"
    
    
        print "TEMP_DIR ======================>", frameoutD
        #Change to the source directory
        iraf.chdir(os.path.dirname(out_tmp_dir))
    
        # Remove an old output file
        misc.fileUtils.removefiles(frameoutD)
    
        # Check if master dark exp_time is compliance with science file
        
        # STEP 1: Subtract masterdark
        iraf.imarith(operand1=frame_in,
                    operand2=master_dark,
                    op='-',
                    result=frameoutD,
                    verbose='yes'
                    )
        log.debug( "SUBTRACT DARK FINISHED ---> %s", t.tac())
        
        # STEP 2: Divide by normalized master FlatField (OK !!)
        ##frameoutF=frameoutD.replace("_D.fits","_D_F.fits")
        ##frameoutF = frameoutF.replace(os.path.dirname(frameoutF), out_tmp_dir)
        ##frameoutF = frameoutF.replace("//","/")
    
        #misc.fileUtils.removefiles(frameoutF)
        ##iraf.imarith(operand1=frameoutD,
        ##             operand2=master_flat,
        ##             op='/',
        ##             result=frameoutF,
        ##             verbose='yes'
        ##             )
        ##log.debug("DIVIDE BY FLAT FINISHED ---> %s",t.tac())
        
        # STEP 3: Compute sky background with only one frame (itself)
        frame_out = frameoutD.replace("_D.fits","_D_S.fits")
        frame_out = frame_out.replace(os.path.dirname(frame_out), out_tmp_dir)
        frame_out = frame_out.replace("//","/")
        log.debug("OUTPUT_FRAME ----> %s", frame_out)
        
        print "frameoutD=", frameoutD
        print "frame_out=", frame_out
        self.computeSkyBkg( frameoutD, frame_out )
        
        log.debug("Background computed ---> %s", t.tac())
        
        
        # STEP 3-2: Divide by normalized master FlatField (OK !!)
        frameoutF = frame_out.replace("_D_S.fits","_D_S_F.fits")
        frameoutF = frameoutF.replace(os.path.dirname(frameoutF), out_tmp_dir)
        frameoutF = frameoutF.replace("//","/")
    
        #misc.fileUtils.removefiles(frameoutF)
        iraf.imarith(operand1=frame_out,
                    operand2=master_flat,
                    op='/',
                    result=frameoutF,
                    verbose='yes'
                    )
        log.debug("DIVIDE BY FLAT FINISHED ---> %s",t.tac())
        
        # STEP 4: Apply pixel mask
        #if appPixMask:
        #    applyPixelMask( frame_out, "/disk-a/caha/panic/DATA/data_alh/output/BadPixMask-20080409192131.pl")
            
            
        # STEP 5: Cleanup & rename output frame
        print "Files to remove : ", frameoutD + " " +  frame_out
        misc.fileUtils.removefiles(frameoutD, frame_out)
        shutil.move(frameoutF, result_frame)
    
        log.debug("Simple reduce finished. %s" , t.tac() )
        run_status=1 # Sucessfull finished
        
        return result_frame

    ################################################################################
    def applyPixelMask(  frame_in, mask ):
        """
        Apply a previously created pixel mask to a frame
        """
    
        log.info ('Start applyPixelMask; applying mask %s to frame %s', mask, frame_in)
        t=misc.utils.clock()
        t.tic()
        
        base=os.path.dirname(frame_in)
        iraf.chdir(base)
        
        iraf.fixpix(images=frame_in,
                    masks=mask,
                    verbose='yes'
                    )
        
        log.debug( "applyPixelMask finished. %s", t.tac() )
    
    ################################################################################
    def computeSkyBkg ( self, image_in, image_out, output_type='-BACKGROUND', params=None, sconfig=None):
     
        """
        \brief Set SExtractor parameters and configuration and run.
        
        \param image_in: Detection (and measurement) image
        \param image_out:  Measurement image
        \param output_type: (optional) background |-background | objects | -objects (see SExtractor doc)
        \param params: (optional) extra output parameters (Parameters object)
        \param config: (optional) configuration options (Config object)
        
        NOTE: If either params or config is not specified, the defaults will be
        used
        """
        
        log.debug("Computing sky background...")
        
        sextractor_dir = self.m_terapix_path
        sextractor_conf_dir = self.m_papi_path+"/irdr/src/config"
        
        params= "-PARAMETERS_NAME " + sextractor_conf_dir + "/default.param " + \
                " -FILTER Y " + \
                " -FILTER_NAME " + sextractor_conf_dir + "/default.conv " + \
                " -FITS_UNSIGNED Y " + \
                " -CATALOG_NAME catalog.cat " + \
                " -DETECT_MINAREA  15 " + \
                " -DETECT_THRESH 5.0 " + \
                " -CHECKIMAGE_TYPE  " + output_type + \
                " -CHECKIMAGE_NAME " +  image_out
        
        try:        
            command = sextractor_dir + "/sex -c " + sextractor_conf_dir + "/default.sex " + " " + image_in + " " + params
            print "COMMAND=", command
            if utils.runCmd (command)==0: # there was some error
                log.debug("Some error on computeSkyBkg" )
                raise Exception("Some error on computeSkyBkg")
        except:
            raise Exception("Some error on computeSkyBkg")
        
        log.debug("end sky background computation")     
     
################################################################################
class ReductionBlock:
    """
    \brief A class that contain a set of files with the same type, filter, mode, texp,... to reduce (science or calib)
    """
    
    m_file_list = []         # original file list
    m_file_list_p = []       # last processed file list
    m_last_output_file = ''  # last output file generated
    m_cl_file_list = []      # original classified (CFits) file list
    m_history = []           # history of operations done to the file list (darked, flattened, avg_combined, med_combined, ...)
    m_output_dir = None      # it will be used for output files and as a working dir (temp, ...)
    m_images_in = None       # Cube for input data
    
    m_chipcode = 1           # Id of the chip (1,2,3,4)
    m_chip_offsets = [1.0, 1.0, 1.0, 1.0]    # matrix coefficients: relative QE sensitivity offsets (filter depending)
    m_filter = None          # Filter of the science data to be reduced (and so, the concerning calibration files:flats, chip_offsets, ...)
    m_naxis1 = 0             # Axis 1 dimension
    m_naxis2 = 0             # Axis 2 dimension
    m_memory = False         # By default, images are not stored in memory, only a pointer to its header
    
    
    def __init__(self, file_list, output_dir="/tmp/", memory=False):
        
        """\param file_list : a python list of pathnames
           \param output_dir :
           \param memory : if True, the images from the file list will be stored into memory 
        """
        
        log.info("Start Reduction block")
        
        if file_list==[]:
            log.error("Empty file list, cannot create Reduction Block")
            raise ExError,"Empty file list"
        
        # Init some variables
        self.m_file_list = file_list
        self.m_file_list_p = self.m_file_list # initially, are the same  
        self.m_nfiles = len(file_list)
        self.m_memory = memory
        self.m_output_dir = output_dir # it will be used for output files and as a working dir (temp, ...)
        self.m_terapix_path = os.environ['TERAPIX']
        self.m_papi_path = os.environ['PAPI_HOME']
        
        log.debug("File_list = %s", file_list)
        self.m_cl_file_list = []
        
        if self.m_memory:
            # Dump images into memory         
            i=0
            for file in self.m_file_list:
                try:
                    log.debug("RB-FILE = %s", file)
                    f = datahandler.ClFits(file)
                    if i==0: # only first time 
                        self.m_naxis1 = f.getNaxis1()
                        self.m_naxis2 = f.getNaxis2()
                        self.m_filter = f.getFilter()
                        #create cube matrix for all the data. Note how naxis are shapped into memory matrix !!        
                        self.m_images_in = numpy.zeros([self.m_nfiles, self.m_naxis2, self.m_naxis1], dtype=numpy.float32)
                    f.printClass()
                except:
                    log.error("Found not well formed FITS file: %s", file)
                    raise
                else:
                    self.m_cl_file_list.append(f)
                    t=pyfits.open(file)
                    self.m_images_in[i,:,:]=t[0].data
                    t.close()
                    i=i+1
           
            if (self.m_nfiles!=i):
                log.error("Some file was not read successfully. Review your files ...:")
                raise
        else:
            for file in self.m_file_list:
                try:
                    log.debug("RB-FILE = %s", file)
                    f = datahandler.ClFits(file)
                    self.m_cl_file_list.append(f)
                except:
                    log.error("Found not well formed FITS file: %s", file)
                    raise
                        
        
         
        # checkFilter
        #if not self.checkFilter():
        #    log.error("Frames dataset does not have all the same filter")
        #    raise ExError, "Filter matching error"
         
        # checkType
        #if not self.checkType():
        #    log.error("Frames dataset does not have all the same type")
        #    raise ExError, "Type matching error"
           
     
    def getList(self):
        return self.m_file_list
    
    def getProcList(self):
        return self.m_file_list_p
    
    def setNewList(self, new_file_list):
        self.m_file_list = new_file_list
        
    def checkFilter(self):
        """Return true is all files in file have the same filter type, false otherwise
        
        \return True or False
        """
        filter_0 = self.m_cl_file_list[0].filter
        for file in self.m_cl_file_list:
            if file.filter != filter_0:
                log.debug("File %s does not match file filter", file)
                return False
        
        log.debug("All files match same file filter")
        return True
        
    def checkExpTime(self):
        """Return true is all files in file have the same Exposition Time, false otherwise
        
        \return True or False
        """
        exptime_0 = self.m_cl_file_list[0].exptime
        for file in self.m_cl_file_list:
            if file.exptime != exptime_0:
                log.debug("File %s does not match file exptime", file)
                return False
        
        log.debug("All files match same file exptime")
        return True
      
    def checkType(self, type_to_check=None):
        """
        Return true is all files in file have the same type(science, dark, flat, ...), false otherwise
        
        \param type_to_check (\c string) type to check to
        \return True or False
        """
        if type_to_check==None:
            type_0 = self.m_cl_file_list[0].type
        else:
            type_0 = type_to_check
             
        for file in self.m_cl_file_list:
            if file.type != type_0:
                log.debug("File %s does not match file type %s", file.pathname, type_0)
                return False
        
        log.debug("All files match same file type")
        return True
                
        
    #def createMasterDark(self, outputfile):
        #"""
        #Create a master dark file from the reduction block files
        
        #\param outputfile (\c string) Name of the outputfile name for the master dark to be created
        #\return Return True master dark file was created, False otherwise
        #"""
        
        #if (self.checkExpTime() and self.checkType("DARK") and self.checkFilter() ):
            #try:
                #md = calDark.MasterDark (self.m_file_list_p, os.path.dirname(outputfile), outputfile)
                #md.createMaster()
                #self.m_last_output_file = outputfile
            #except:
                #log.error("Error creating master DARK file")
                #raise
            #else:
                #return
        #else:
            #log.error("Cannot create master dark. Files not same type")
            #raise ExError, "File mismatch"
        
    #def createMasterDFlat(self, outputfile):
        #"""
        #Create a master dome flat file from the reduction block files
        
        #\param outputfile (\c string) Name of the outputfile name for the master to be created
        #\return Return True master file was created, False otherwise
        #"""
        
        ## TODO : search for dome flats (on and off)
        #if True:
            #try:
                #md = calDomeFlat.MasterDomeFlat (self.m_file_list_p, os.path.dirname(outputfile), outputfile)
                #md.createMaster()
                #self.m_last_output_file = outputfile
            #except:
                #log.error("Error creating master Dome Flat file")
                #raise
            #else:
                #return True
        #else:
            #log.error("Cannot create master frame. Files not same type")
            #raise ExError, "Frames mismatch"
            #return False
                
    def createSuperFlat(self, master_dark, output_file):
        """
        Create a master super flat file from the reduction block files (science type)
        
        \param outputfile (\c string) Name of the outputfile name for the master to be created
        \return Return True master file was created, False otherwise
        """
        # Step 0: Check master dark is appropiate
        # @TODO
        
        log.info("Start createSuperFlat")
        # Remove an old files
        misc.fileUtils.removefiles(output_file)
        
        # Step 1: Check file group is valid
        if (self.checkExpTime() and self.checkType("SCIENCE") and self.checkFilter() ):
            try:
                out_1 = self.subtractDark(master_dark)
                rb_1 = ReductionBlock(out_1)   
                superFlat_1 = self.m_output_dir + 'superF1.fits'
                rb_1.combineFrames(superFlat_1)
                out_2 = rb_1.applyFlat(superFlat_1)              
                rb_2 = ReductionBlock(out_2)
                masks = rb_2.createObjMasks()
                rb_2.combineFrames(output_file)
                
                self.m_last_output_file = output_file
            except:
                log.error("Error creating SuperFlat file")
                raise
            else:
                log.info("Superflat %s created successfully !!!", output_file)
        else:
            log.error("Cannot create Super Flat. Files not same type")
            raise ExError, "Not right frames types"
                   
    ################################################################################
    def computeSkyBkg ( self, image_in, image_out, output_type='-BACKGROUND', params=None, sconfig=None):
     
        """
        \brief Set SExtractor parameters and configuration and run.
        
        \param image_in: Detection (and measurement) image
        \param image_out:  Measurement image
        \param output_type: (optional) background |-background | objects | -objects (see SExtractor doc)
        \param params: (optional) extra output parameters (Parameters object)
        \param config: (optional) configuration options (Config object)
        
        NOTE: If either params or config is not specified, the defaults will be
        used
        """
     
        log.debug("Computing sky background...")
        
        sextractor_dir = self.m_terapix_path
        sextractor_conf_dir = self.m_papi_path+"/irdr/src/config"
        
        params= "-PARAMETERS_NAME " + sextractor_conf_dir + "/default.param " + \
                " -FILTER Y " + \
                " -FILTER_NAME " + sextractor_conf_dir + "/default.conv " + \
                " -FITS_UNSIGNED Y " + \
                " -CATALOG_NAME catalog.cat " + \
                " -DETECT_MINAREA  15 " + \
                " -DETECT_THRESH 5.0 " + \
                " -CHECKIMAGE_TYPE  " + output_type + \
                " -CHECKIMAGE_NAME " +  image_out
        
        
        try:
            command = sextractor_dir + "/sex -c " + sextractor_conf_dir + "/default.sex " + " " + image_in + " " + params
        except:
            raise
        if utils.runCmd (command)==0: # there was some error
            log.debug("Some error on computeSkyBkg" )
            raise "Some error on computeSkyBkg"
        
        log.debug("end sky background computation")
        
                    
    def I_createSuperFlat(self, images_in=None, output_file='/tmp/gain.fits', bad_pixel_mask=None):
        """
            Create a super flat and gain map using the dithered images. (IRAF)
                         
            This function is a wrapper for mksuperflat.py that uses imcombine(IRAF) and gainmap(IRDR)              
        
            VERSION
                1.0, 20090909 by jmiguel@iaa.es
        
        """
        
        log.info('Start I_createSuperFlat...')
                                  
        if images_in==None:
            images = self.m_images_in
        else:
            images = images_in
            
                                  
        #STEP 1: Write out current images to FITS files
        output_dir='/tmp/'
        suffix='_sf_'+self.m_filter+'.fits'
        listfile=output_dir+'/sf_'+self.m_filter+'.pap'
        # Remove old files
        if os.path.exists(listfile): os.remove(listfile)
        for file in glob.glob(output_dir+"*_sf_"+self.m_filter+".fits"):
            os.remove(file)
            
        self.writeOutData(images, output_dir, suffix, listfile)
            
        # STEP 2: Create superflat using IRAF combine
        superflat_file = output_dir+'sf_'+self.m_filter+'.fits'
        log.info("SUPERFLAT=%s", superflat_file)
        makesuperflat( listfile, superflat_file )
                                        
        # STEP 3: Create gain map and set to 0 bad pixels using gainmap.c from IRDR
        #gainmap sflat.fits gain.fits $nsig $nxblock $nyblock $mingain $maxgain
        nxblock=16
        nyblock=16
        #The next values are to find out bad pixels 
        nsig=5
        mingain=0.7
        maxgain=1.3
        if output_file == None:
            gainfile = output_dir+'gain_'+self.m_filter+'.fits'
        else:
            gainfile = output_file
            
        gain_cmd=self.m_papi_path+'/irdr/bin/gainmap '+ superflat_file + '  ' + gainfile +' '+ str(nsig) + '  ' + str(nxblock) + '  '+ str(nyblock) + '  ' + str(mingain) + '  ' + str(maxgain) 
        if utils.runCmd( gain_cmd )==0:
            log.error("Some error while creating gainmap ....")
            return
                         
        # STEP 4 : Add external bad pixel mask (optional)
        # Check the badmask file exists and
        if bad_pixel_mask!=None and not os.path.exists( bpm_file ):
            log.debug('No external Bad Pixel Mask found. Cannot find file : "%s"' %bad_pixel_mask)
        elif bad_pixel_mask!=None:
            iraf.imarith(operand1=gainfile,
                  operand2=bad_pixel_mask,
                  op='*',
                  result='/tmp/gain_bpm.fits',
                  verbose='yes'
                  )
            #os.remove( gain_file )      
            os.move( '/tmp/gain_bpm.fits', gain_file )     
                   
    def MI_skyFilter (self, images_in=None, source='mem', gainfile=None, mask='nomask'):
        """
            For each input image, a sky frame is computed by combining a certain number of the closest images, 
            then this sky frame is subtracted to the image and the result is divided by the master flat; 
                         
            This function is a wrapper for skyfilter.c (IRDR)              
        
            VERSION
                1.0, 20090909 by jmiguel@iaa.es
        
            TODO: extended objects !!!!
        """
        
        log.info('Start MI_skyFilter...')
                                  
        if source=='mem':
            if images_in==None:
                images = self.m_images_in
            else:
                images = images_in
            
            #STEP 1: Write out current images to FITS files
            output_dir=self.m_output_dir
            suffix='_skf_'+self.m_filter+'.fits'
            listfile=output_dir+'skyfilter_'+self.m_filter+'.pap'
            # Remove old files
            if os.path.exists(listfile): os.remove(listfile)
            for file in glob.glob(output_dir+"*_skf_"+self.m_filter+".fits"):
                os.remove(file)
                
            self.writeOutData(images, output_dir, suffix, listfile)
        
        elif source=='file':
            # Check if the input is a file, then with master object mask and offsets to compute a better sky                                
            if os.path.isfile(images_in):
                listfile=images_in
            else:
                log.error("File %s does not exists !!", listfile)      
        
            
        if gainfile==None:
            m_gainfile=self.m_output_dir+'gain_'+self.m_filter+'.fits'
                                       
        # STEP 2: Call skyfilter
        halfnsky=4
        destripe='none'
            
        gain_cmd=self.m_papi_path+'/irdr/bin/skyfilter '+ listfile + '  ' + m_gainfile +' '+ str(halfnsky)+' '+ mask + '  ' + destripe 
        if utils.runCmd( gain_cmd )==1: # All was OK
            # Rename output files
            for file in glob.glob(self.m_output_dir+'/*.fits.skysub'):
                shutil.move(file, file.replace('.fits.skysub', '.skysub.fits'))                          
                
    def M_createSuperFlat(self , images_in=None, output_file=None, mask=None, fast=1):
        """
        ;
        ; NAME
        ;   M_createSuperFlat
        ;
        ; PURPOSE
        ;   Derive super sky flat from the image cube
        ;
        ; INPUT
        ;   images_in - image cube that contains the unflattened raw data.
        ;       The format should be nframes * nx * ny 
        ;       If None, then class member m_images_in is used
        ;
        ; OPTIONAL INPUT
        ;   mask - integer cube with dimensions identical to the data cube.
        ;       1 in the mask represent good data, 0 represent bad data.
        ;       Useful to prevent bias caused by stars in the fast mode.
        ;
        ;   output_file - if provided, a fits file with super flat is created
        ;
        ; KEYWORD
        ;   fast - set to 1 for fast flatfield (default).  In this way,
        ;       the cube will be averaged along its 0th dimenssion to
        ;       form a flat.  Set to 0 to use median instead of average.
        ;       This is slower but should be used if there is not a good
        ;       object mask.
        ;
        ; VERSION
        ;   1.0, 20090723 by jmiguel@iaa.es
        ;
        """
        
        log.info("Start M_createSuperFlat")
        
        # Get a copy of the input images, so any modification is done into original data
        if images_in==None:
            images = numpy.copy(self.m_images_in)
        else:
            images = numpy.copy(images_in)
            
        nimages = images.shape[0]
            
        # Step 1: Check file group is valid
        if (self.checkType("SCIENCE") and self.checkFilter() ):
            # Fast mode: cube average
            if fast==1:
                weights=numpy.ones(images.shape, dtype=numpy.int)
                A = ~numpy.isfinite(images)
                weights[A]=0
                master_mask = 1.0
                
                if mask!=None:
                    weights[mask==0]=0
                    images=numpy.where(mask==0, numpy.nan)
                     
                    master_mask = mask.sum(axis=0)*1.0 # only to float convert
                    master_mask = numpy.where(master_mask < nimages, numpy.nan, master_mask) # only bad pixel if is in each image
                    master_mask = master_mask/nimages # set bad pixel to 1
                
                for i in range(0, nimages):
                    images[i,:,:] = images[i,:,:]/numpy.median( images[i,:,:].reshape( images[i,:,:].shape[0]*images[i,:,:].shape[1] ) )*master_mask     
                        
                sflat = images.sum(axis=0)/weights.sum(axis=0)
            
            #Normal mode: cube median
            else:
                sflat = numpy.zeros([images.shape[1], images.shape[0]])
                master_mask = 1.0
                
                if mask!=None:
                    images=numpy.where(mask==0, numpy.nan)
                     
                    master_mask = mask.sum(axis=0)*1.0 # only to float convert
                    master_mask = numpy.where(master_mask < nimages, numpy.nan, master_mask) # only bad pixel if is in each image
                    master_mask = master_mask/nimages # set bad pixel to 1
                
                for i in range(0, nimages):
                    images[i,:,:] = images[i,:,:]/numpy.median( images[i,:,:].reshape( images[i,:,:].shape[0]*images[i,:,:].shape[1] ) )*master_mask                
                
                #for i in range(0, images.shape[1]):
                #    for j in range (0, images.shape[2]):
                #        sflat[i,j]= numpy.median( images[*,i,j] )
                        
                sflat = numpy.median( images, axis=0 )
                        
            images = 0.0  # to free memory ??
            sflat = numpy.where (sflat< 0.4 , numpy.nan, numpy.where(sflat>1.6, numpy.nan, sflat))       
             
            # Write out to file
            if output_file!=None:  
                misc.fileUtils.removefiles(output_file)
                hdu = pyfits.PrimaryHDU(sflat)
                #hdu.scale('float32') # importat to set first data type
                #hdu.data=sflat     
                hdulist = pyfits.HDUList([hdu])
                hdu.header.update('OBJECT','SUPER_FLAT')
                hdulist.writeto(output_file)
                hdulist.close(output_verify='ignore')
                    
                       
    def subtractDark(self, master_dark):
        """Subtract a master dark to each file of the reduction block
        """
        
        t=misc.utils.clock()
        t.tic()
        
        log.info("Start subtractDark")
        if not os.path.exists(master_dark):
            raise ExError, 'Master frame does not exist'
        
        framelist = ''
        framelist_out = ''
        framelist = utils.listToString(self.m_file_list_p)
        #for iframe in self.m_file_list_p:
        #    framelist += iframe + ' , '
            
        base, infile   = os.path.split(self.m_file_list_p[0])    
        framelist_out = framelist.replace(".fits","_D.fits")
        framelist_out = framelist_out.replace(base, self.m_output_dir).replace('//','/')  # IRAF interpreta '//' como ''
        log.debug("FrameList_Out = %s", framelist_out)
        # Remove an old files ( we need to conver to python list)
        ret_list = []
        ret_list = [i.replace(".fits","_D.fits").replace(base, self.m_output_dir).replace('//','/') for i in self.m_file_list_p]
        for f in ret_list:
            misc.fileUtils.removefiles(f)

        # STEP 1: Subtract masterdark
        iraf.imarith(operand1=framelist,
                  operand2=master_dark,
                  op='-',
                  result=framelist_out,
                  verbose='yes'
                  )
        # STEP 2: Check dark EXPTIME is complaint with files EXPTIME
        fdark = datahandler.ClFits(master_dark)
        if (not fdark.exptime == self.m_cl_file_list[0].exptime):
            log.debug("Adjusting EXPTIME factor ---> %s", float(fdark.exptime/self.m_cl_file_list[0].exptime))
            iraf.imarith(operand1=framelist_out,
                  operand2=float(fdark.exptime/self.m_cl_file_list[0].exptime), #ESTO ESTA MAAAALLLL !!!!!!!!!!!!!!!!!!!!!!
                  op='*',
                  result=framelist_out,
                  verbose='yes'
                  )
        log.debug("SUBTRACT DARK FINISHED ---> %s", t.tac())
        
        self.m_file_list_p = ret_list
        self.m_history.append("DARK")
        return ret_list
    
    def M_subtracDark(self, master_dark):
        """Subtract a master dark to each file of the reduction block
           The master dark is properly scaled before subtraction
           
           TODO: apply dark model !!!!!!!
        """  
        log.debug("Dark subtraction")
        
        
        mdark=pyfits.open(master_dark)
        for i in range(0,self.m_nfiles):
            self.m_images_in[i,:,:]=self.m_images_in[i,:,:]-mdark[0].data*float(self.m_cl_file_list[i].expTime()/mdark[0].header['EXPTIME'])
            self.m_cl_file_list[i].my_header.add_history('Dark Subtracted =%s' %master_dark)
      
    def M_appFlatField(self,  master_flat, images_in=None):
        """apply a master flat field to each file of the images_in data cube
        """
        log.debug("Applying flat-field")
          
        
        if images_in==None:
            images = self.m_images_in
        else:
            images = images_in
              
        mflat=pyfits.open(master_flat)
        for i in range(0,self.m_nfiles):
            images[i,:,:]=images[i,:,:]/mflat[0].data
            images[i,:,:]=numpy.where(numpy.isnan(images[i,:,:]), 0, images[i,:,:])
            images[i,:,:]=numpy.where(numpy.isneginf(images[i,:,:]), 0, images[i,:,:])
            images[i,:,:]=numpy.where(numpy.isposinf(images[i,:,:]), 0, images[i,:,:])
            if numpy.isneginf(images[i,:,:].any()):
                log.error("Neg INF   !!!!!!!!!!!!!")    
            if images_in==None:
                self.m_cl_file_list[i].addHistory('Flat field applied =%s' %master_flat)
        mflat.close()
                
            
    def computeSkyBackground_i(self, i_frame, n_frames, output_type='BACKGROUND', output_file=None ):
        """\brief Compute the sky background for the frame 'i_frame', using the nearest 2*n_frames
           \param i_frame index referenced frame to compute the related sky-background
           \param n_frames number of nearest 2*frames to use for the sky backgorund
           \param output_type: background | sub_background | objects
           \param ouput_sky sky background returned
        """
         
        
            
        # First, some checks
        # Check data type is SCIENCE or SKY (TODO)
        if not self.checkType('SCIENCE') and not self.checkType('DOME_FLAT_LAMP_OFF'):
            log.error("Error, Type mismatching")
            raise ExError,"Type mismatching"
            
         
        # Check 'i_frame' is in m_file_list
        if i_frame < 0 or i_frame>len(self.m_file_list):
            log.error("Not provided suitable parameters")
            raise ExError,"Bad parameters"
        
        # ******************************************************************* 
        # SINGLE : Compute the own background, not using anyother near frame
        # *******************************************************************
        if n_frames==0:
            self.computeSkyBkg ( self.m_file_list[i_frame], output_file, output_type)
            return    

        # ******************************************************************* 
        # MULTI 
        # *******************************************************************
        # adjust and check there are enought frames 
        m_n_frames = 0 
        if 2*n_frames>len(self.m_file_list)-1:
            m_n_frames = len(self.m_file_list)-1
        else: m_n_frames = n_frames
                
        log.debug("M_N_FRAMES=%d",m_n_frames)        
        # Sort the framelist
        self.sortByDateObs()
            
            
        # Compound sky list (PRE+POST frames)
        # *** Init some counts
        sky_list = []
        pre_pend  = m_n_frames - i_frame # pending pre-frames to add (from the post-i)
        if pre_pend<0: 
            print "PR_P=%d" %pre_pend
            pre_pend=0
        post_pend = m_n_frames - (len(self.m_file_list)-1- i_frame) # pending post-frames to add (from the pre-i)
        if post_pend<0:
            print "PS_P=%d" %post_pend 
            post_pend=0
        # En teoria, no puedo tener pre_pend!=0 y post_pend!=0
        
        log.debug("PRE_PEND=%d",pre_pend)
        log.debug("POST_PEND=%d",post_pend)
        
        # *** PRE-frames 
        if (i_frame-(m_n_frames + post_pend))>=0:
            sky_list = self.m_cl_file_list[(i_frame-(m_n_frames+post_pend)):i_frame]
        else: # not enought pre-frames
            sky_list = self.m_cl_file_list[0:i_frame]
        
        log.debug("SKY_LIST(pre)=%s",sky_list)
        
        # *** POST-frames 
        if (i_frame+(m_n_frames+pre_pend))<=len(self.m_file_list)-1:
            sky_list += self.m_cl_file_list[i_frame+1:i_frame+(m_n_frames+pre_pend+1)]
        else: # not enought post-frames
            sky_list += self.m_cl_file_list[i_frame+1:]
            
        log.debug("SKY_LIST(post)=%s",sky_list)
        sl = []
        for k in sky_list:
            sl.append(k.pathname) 
        # Now, create a temp reduction block
        tmp_rb = ReductionBlock(sl)
        if output_file==None:
            output_file = '/tmp/sky_i_%i_n_%i.fits' %(i_frame, m_n_frames)
        tmp_rb.combineFrames(output_file)
        
        # Check output_type 
        if output_type=='-BACKGROUND':
            l_output_file = output_file.replace('.fits','_S.fits')
            misc.fileUtils.removefiles(l_output_file)
            # Subtract background
            iraf.imarith(operand1=self.m_file_list[i_frame],
                  operand2=output_file,
                  op='-',
                  result=l_output_file,
                  verbose='yes'
                  )
        log.debug("BACKGROUND SUBTRACTED ---> %s", l_output_file)
        #Rename the output file according our input parameters
        shutil.copy(l_output_file, output_file)
        
        self.m_last_output_file  = l_output_file 
            
            
        
    def simpleReduce(self, i_frame, n_frames,  output_file, master_flat=None, bpmask_file=None, do_astrometry=False):
        """
            TODO
        """
        otemp_1 = '/tmp/temp1.fits'
        otemp_2 = '/tmp/temp2.fits'
        
        #Clean up
        misc.fileUtils.removefiles(otemp_1, otemp_2, output_file)
        
        #**************************************        
        #*** Reduction with a single frame
        #**************************************        
        if n_frames==0: 
            # STEP 1: Divide by master flat
            if master_flat!=None:
                iraf.imarith(operand1=self.m_file_list[i_frame],
                    operand2=master_flat,
                    op='/',
                    result=otemp_1,
                    verbose='yes'
                    )
                log.debug("Flat fielding FINISHED")
            else:
                otemp_1 = self.m_file_list[i_frame]
                
            # STEP 2:Compute sky background with N-nearest and subtract (dark subtraction implicit)
            tl=[]
            tl.append(otemp_1)
            tmp_rb = ReductionBlock(tl)
            tmp_rb.computeSkyBackground_i(0, 0, '-BACKGROUND', otemp_2)
                
            # STEP 3: Apply Bad Pixel Mask
            if bpmask_file!=None:
                tl=[]
                tl.append(otemp_2)     
                tmp_rb = ReductionBlock(tl)
                tmp_rb.applyBadPixelMask([0], bpmask_file)
        
        #**************************************        
        #**** Reduction with N-nearest frames
        #**************************************
        elif n_frames>0:
            # STEP 1:Compute sky background with N-nearest and subtract (dark subtraction implicit)
            self.computeSkyBackground_i(i_frame, n_frames, '-BACKGROUND', otemp_1 )
            # STEP 2: Divide by master flat
            if master_flat!=None:
                iraf.imarith(operand1=otemp_1,
                    operand2=master_flat,
                    op='/',
                    result=otemp_2,
                    verbose='yes'
                    )
                log.debug("Flat fielding FINISHED")
            else:
                otemp_2 = otemp_1
                
            # STEP 3: Apply Bad Pixel Mask
            if bpmask_file!=None:
                tl=[]
                tl.append(otemp_2)      
                tmp_rb = ReductionBlock(tl)
                tmp_rb.applyBadPixelMask([0], bpmask_file)
                
        # Last, rename output file
        shutil.move(otemp_2, output_file)
        self.m_last_output_file  = output_file
        
    ################################################################################## 
    def createObjMasks(self):
        """Find objects into the frames and create an object masks for each file; 
            write the mask file into the object frame.
        """
        log.debug("TODO")
        
    ##################################################################################
    def applyFlat(self, master_flat):
        """Divide all the frames into the Reduction Block by a given master (normalized) flat field
        """
                 
        t=misc.utils.clock()
        t.tic()
        
        if not os.path.exists(master_flat):
            raise ExError, 'Master frame does not exist'
        
        # Check filter
        df = datahandler.ClFits(master_flat)
        if self.m_cl_file_list[0].filter!= df.filter:
            raise ExError, 'Master frame does not match frame dataset filter'
        
        framelist = ''
        framelist_out = ''
        framelist = utils.listToString(self.m_file_list_p)
        #for iframe in self.m_file_list_p:
        #    framelist += iframe + ' , '
            
        base, infile = os.path.split(self.m_file_list_p[0])    
        framelist_out = framelist.replace(".fits","_DF.fits")
        framelist_out = framelist_out.replace(base, self.m_output_dir).replace('//','/')
        # Remove an old files ( we need to conver to python list)
        ret_list = []
        ret_list = [i.replace(".fits","_DF.fits").replace(base, self.m_output_dir).replace('//','/') for i in self.m_file_list_p]
        for f in ret_list:
            misc.fileUtils.removefiles(f)
        # STEP 1: Divide master flat
        iraf.imarith(operand1=framelist,
                  operand2=master_flat,
                  op='/',
                  result=framelist_out,
                  verbose='yes'
                  )
        log.debug("Flat fielding FINISHED ---> %s", t.tac())
        
        self.m_file_list_p = ret_list
        self.m_history.append("FLAT")
        return ret_list

    ################################################################################## 
    def combineFrames(self, output_file, type='median', mask=None, params=None):
        """Combine all the frames into the Reduction Block using 'type' operation (median, mean=average)
        """
        # TODO -->masks !!!
                 
        t=misc.utils.clock()
        t.tic()
        
        # Some checks
        # Type operation
        if type!='median'and type!='mean':
            raise ExError, 'Combine operation type not supported'
        
        # Filter and Type checking 
        if not self.checkFilter() or not self.checkExpTime():
            log.error("Frames dataset does not have all the same filter or ExpTime")
            raise ExError, "Filter matching error"
         
        # checkType
        if not self.checkType():
            log.error("Frames dataset does not have all the same type")
            raise ExError, "Type matching error"
        
        framelist = ''
        framelist = utils.listToString(self.m_file_list_p)
        log.debug("Frame list to combine = [%s]", framelist)
        #for iframe in self.m_file_list:
        #    framelist += iframe + ' , '
        
        # Set iraf working directory    
        base, infile   = os.path.split(self.m_file_list_p[0])    
        iraf.chdir(base)
        # Remove an old file
        misc.fileUtils.removefiles(output_file)
        # STEP 1: Combine data files
        iraf.combine(input=framelist,
                     output=output_file,
                     combine=type,
                     ccdtype='none',
                     reject='sigclip',
                     lsigma=3,
                     hsigma=3,
                     subset='yes',
                     scale='mode'
                     #masktype='none'
                     #verbose='yes'
                     #scale='exposure',
                     #expname='EXPTIME'
                     #ParList = _getparlistname ('flatcombine')
                     )

        log.debug("Frame combining  FINISHED ---> %s", t.tac())
        
        self.m_last_output_file = output_file
        self.m_history.append("COMBINE")
        return True
    
    ##################################################################################
    def createBadPixelMask(self, show=False, output_file=None):
        """Make a bad pixel mask (hot and cold pixels) from a set of images ( lamp_on and lamp_off frames and darks)
        """
        
        log.debug('createBadPixelMask started')
        t=misc.utils.clock()
        t.tic()
        
        darks_frames=[]
        flats_off_frames=[]
        flats_on_frames=[]
    
        
        #Change to the source directory
        base, infile   = os.path.split(self.m_file_list_p[0])    
        iraf.chdir(self.m_output_dir)
        
        #STEP 1: classify/split the frames in 3 sets (DOME_FLAT_LAMP_ON, DOME_FLAT_LAMP_OFF, DARKS)
        #        and create string list for IRAF tasks
        for iframe in self.m_cl_file_list:
            if iframe.getType()=='DOME_FLAT_LAMP_OFF':
                flats_off_frames.append(iframe.pathname)
            elif iframe.getType()=='DOME_FLAT_LAMP_ON':
                flats_on_frames.append(iframe.pathname)
            elif iframe.getType()=='DARK':
                darks_frames.append(iframe.pathname)
            else:
                # reject the frame
                log.error("DISCARTING: Frame %s is neither dome_flat nor a dark frame", iframe.pathname)
                
        #Check whether there are enought calib frames
        if (len(flats_off_frames)<1 or len(flats_on_frames)<1 or abs(len(flats_off_frames)-len(flats_off_frames))>10 
            or len(darks_frames)<3):
            log.error("There are not enought calib frames for create BPM !!")
            raise ExError, "Not enought calib frames"
        
        #STEP 2: Create the master dark
        try:
            rb_dark = ReductionBlock(darks_frames)
            master_dark = '/tmp/master_dark.fits'
            rb_dark.createMasterDark(master_dark)
        except:
            log.error("Error building master dark")
            raise
            
        
        #STEP 3: Subtrac the master dark to each dome flat
        try:       
            rb_flats_off = ReductionBlock(flats_off_frames)
            rb_flats_off.subtractDark(master_dark)
            rb_flats_on = ReductionBlock(flats_on_frames)
            rb_flats_on.subtractDark(master_dark)
        except:
            log.error("Error subtracting master dark")
            raise 
        
        #STEP 4: Combine dome dark subtracted flats (on/off)
        try:
            flat_off_comb = self.m_output_dir + 'flat_off_comb.fits'
            rb_flats_off.combineFrames(flat_off_comb, 'median')
            flat_on_comb = self.m_output_dir + 'flat_on_comb.fits'
            rb_flats_on.combineFrames(flat_on_comb, 'median')
        except:
            log.error("Error while combining flats frames")
            raise 
        
        #STEP 5: Compute flat_low/flat_high
        flat_ratio = self.m_output_dir + 'flat_ratio.fits'
        misc.fileUtils.removefiles(flat_ratio)
        log.debug( 'Computing FlatRatio')
        iraf.imarith(operand1=flat_off_comb,
                 operand2=flat_on_comb,
                 op='/',
                 result=flat_ratio,
                 )
        log.debug("Flat Ratio built")
                     
        #STEP 6: Create BPM
        log.debug( 'Now, computing bad pixel mask')
        dt = datetime.datetime.now()
        if output_file==None:
            outF = self.m_output_dir + 'BadPixMask'+dt.strftime("-%Y%m%d%H%M%S")
        else:
            outF = output_file
        misc.fileUtils.removefiles(outF + ".pl")
        iraf.ccdmask(image=flat_ratio,
                 mask=outF,
                 lsigma=4,
                 hsigma=4,
                 )
        
        # Clean up tmp files
        misc.fileUtils.removefiles(flat_off_comb, flat_on_comb, flat_ratio)
        # Change back to the original working directory
        iraf.chdir()
        
        if show:
            iraf.display(image=outF, frame=1)
                        
        
        log.info('Bad pixel mask created')
        log.debug("Time elapsed : [%s]" , t.tac() )

        return True

    ##################################################################################
    def applyBadPixelMask( self, i_frames, mask ):
        """
        \brief Apply a previously created pixel mask to frame elements indexed by 'i_frames' list in the list of frames
        \param i_frames list of indexes
        \param mask badpixel mask file to apply
        \return True or False
        """
        
        t=misc.utils.clock()
        t.tic()
     
        base=os.path.dirname(self.m_file_list[i_frames[0]])
        iraf.chdir(base)
        
        for i in i_frames:
            if i>=0 and i<len(self.m_file_list):
                log.debug ('Start applyBadPixelMask; applying mask %s to frame %s', mask, self.m_file_list[i])
                iraf.fixpix(images=self.m_file_list[i],
                     masks=mask,
                     verbose='yes'
                     )
            else:
                log.error("Index error, pixel mask to applied to frame")
                     
                     
        log.debug( "applyBadPixelMask finished. %s", t.tac() )

    ##################################################################################
    def mathOp( self, operator, outputFile=None ):
        """
        \brief Apply a math operation to the selected files
        \param operator (-, +, /)
        \param outputFile (by default /tmp/op.fits)
        \return True or False
        """
        t=misc.utils.clock()
        t.tic()

        log.debug("Start mathOp")
        if outputFile==None:
            outputFile = '/tmp/op.fits'

        if operator!='+' and operator!='-' and operator!='/':
            return

        # Remove an old output file
        misc.fileUtils.removefiles(outputFile)
        if (operator=='+' and len(self.m_file_list_p)>2):
            framelist = ''
            framelist = utils.listToString(self.m_file_list_p)
            log.debug("Frame list to combine = [%s]", framelist)
            iraf.combine(input=framelist,
                     output=outputFile,
                     combine='average',
                     ccdtype='none',
                     reject='sigclip',
                     lsigma=3,
                     hsigma=3,
                     subset='yes',
                     scale='mode'
                     #masktype='none'
                     #verbose='yes'
                     #scale='exposure',
                     #expname='EXPTIME'
                     #ParList = _getparlistname ('flatcombine')
                     )
        else:
            """iraf.mscarith(operand1=self.m_file_list_p[0],
                      operand2=self.m_file_list_p[1],
                      op=operator,
                      result=outputFile,
                      extname="",
                      verbose='yes'
                      )
            """
            iraf.imarith(operand1=self.m_file_list_p[0],
                      operand2=self.m_file_list_p[1],
                      op=operator,
                      result=outputFile,
                      verbose='yes'
                      )
        log.debug("Operation FINISHED ---> %s", t.tac())
        return

    ################################################################################## 
    def SNAP_like_reduction(self, master_flat=None):
       """ Make a reduction block procedure based on the SNAP pipeline"""

        ### TODO ###
        
        
    ############### Private functions ################################################ 
    def sortByDateObs(self):
        """ Sort the object frame list by DATE-OBS keyword"""
        
        self.m_cl_file_list.sort(self.compDateObs)
        log.debug("List sorted...")
        for f in self.m_cl_file_list:
            print "File: %s  DateObs: %s" %(f.filename,f.datetime_obs)
        
    def compDateObs(self, e1, e2):
        """ Compare routine used to sort the frame classified list"""
        
        if e1.datetime_obs < e2.datetime_obs:
            return -1
        elif e1.datetime_obs > e2.datetime_obs:
            return +1
        else:
            return 0
        
####################################################################################################################################
    def read_red_parameters(self):
        """ Read main reduction parameters """
        # TBD   
        
        log.debug("Reading reduction parameters ...")
        pass
    
    def checkData(self):
        """ Check data files are ok (filter, date, exists, size, tc ...) """
        # TBD
        log.debug("Checking data...")
        pass       
    
    def getMean(self):
        """only for a test"""
           
        #for i in range(0, self.m_nfiles):
        mean=numpy.mean(self.m_images_in[:,:,:])
        log.debug("MEAN value = %s", mean)
       
    def writeOutData(self, images_in=None, output_dir='/tmp/', suffix=None, file_list=None):
        """ Write out in memory processed data images into a new output FITS files with updated header
            read back from the original FITS files.
        
        INPUTS
                images     -  3 dimensional array with images to write out
                
                output_dir -  directory where FITS images will be created
                      
                suffix     -  optional suffix to use in the filename output
                      
                file_list  -  optional file name for a file listing the FITS files created. If None, no file list will be created
        OUTPUTS
        
        """
           
        #PENDIENTE COPIAR HEADER ORIGINAL actualizado !!!
        log.debug("Writing out data ...")
        
        if images_in==None:
            images = self.m_images_in
        else:
            images = images_in
            
        if file_list!=None:
            f_list = open(file_list, "w")
            
                
        for i in range(0, self.m_nfiles):
            if suffix!=None:
                output_filename=(output_dir+self.m_cl_file_list[i].filename).replace('.fits', suffix)
            else:
                output_filename=output_dir+self.m_cl_file_list[i].filename
            
            misc.fileUtils.removefiles(output_filename)
            hdu = pyfits.PrimaryHDU(images[i,:,:], self.m_cl_file_list[i].my_header)
            hdu.scale('float32') # importat to set first data type
            hdulist = pyfits.HDUList([hdu])
            hdu.header.update('OBJECT','SCIENCE')
            hdulist.writeto(output_filename)
            hdulist.close(output_verify='ignore')
            log.debug("New file created: %s ", output_filename)
            
            if file_list!=None:
                f_list.write(output_filename+"\n")
             
        if file_list!=None: f_list.close()
              
    def M_appChipOffset(self, images_in=None):
        """Apply chip offset relative QE (sensitivity) correction. Only for CHIPs != 1 """
                  
        log.debug("Applying chip-offsets")
                  
        if images_in==None:
            images = self.m_images_in
        else:
            images = images_in
                      
        if self.m_chipcode!=1:
            for i in range(0,self.m_images_in.shape[0]):
                images[i,:,:] = images[i,:,:]/self.m_chip_offsets[self.m_chipcode-1]
                #Replace bad values with 0
                images[i,:,:]=numpy.where(numpy.isnan(images[i,:,:]), 0, images[i,:,:])
                images[i,:,:]=numpy.where(numpy.isneginf(images[i,:,:]), 0, images[i,:,:])
                images[i,:,:]=numpy.where(numpy.isposinf(images[i,:,:]), 0, images[i,:,:])
    
    def M_subtracMedSky(self, images_in=None, fast=0, mask=None, output_file=None):
        """ Compute the median sky of the cube of images and then subtract it to each one.
            
            NOTE:I have not very clear what it does, because the 'hard' sky background is not removed.
            In the manual, it said that 'remove residual background pattern'
            
            ESTA RUTINA NO FUNCIONA BIEN !!!!! no se que "tipo" de cielo restaaaaaaaaaaa!!!!!
        """
          
        log.info("Start M_subtracMedSky")
        
        if images_in==None:
            images_ = self.m_images_in
        else:
            images_ = images_in
            
        # Get a copy of images_
        images = numpy.copy(images_)
        nimages = images.shape[0]
        
        # Step 1: Check file group is valid
        if (self.checkType("SCIENCE") and self.checkFilter() ):
            # Fast mode: cube average
            if fast==1:
                weights=numpy.ones(images.shape, dtype=numpy.int)
                A = ~numpy.isfinite(images)
                weights[A]=0
                master_mask = 1.0
                
                if mask!=None:
                    weights[mask==0]=0
                    images=numpy.where(mask==0, numpy.nan)
                    master_mask = mask.sum(axis=0)*1.0 # only to float convert
                    master_mask = numpy.where(master_mask < nimages, numpy.nan, master_mask) # only bad pixel if is in each image
                    master_mask = master_mask/nimages # set bad pixel to 1
                
                for i in range(0, nimages):
                    images[i,:,:] = images[i,:,:]-numpy.median( images[i,:,:].reshape( images[i,:,:].shape[0]*images[i,:,:].shape[1] ) )*master_mask     
                        
                sky_bg = images.sum(axis=0)/weights.sum(axis=0)
            
            #Normal mode: cube median (fast==0)
            else: 
                sky_bg = numpy.zeros([images.shape[1], images.shape[0]])
                master_mask = 1.0
                
                if mask!=None:
                    images=numpy.where(mask==0, numpy.nan)
                     
                    master_mask = mask.sum(axis=0)*1.0 # only to float convert
                    master_mask = numpy.where(master_mask < nimages, numpy.nan, master_mask) # only bad pixel if is in each image
                    master_mask = master_mask/nimages # set bad pixel to 1
                
                for i in range(0, nimages):
                    log.debug("------------------>MEDIAN (antes)= %s",  numpy.median( images[i,:,:],axis=None ))
                    images[i,:,:] = images[i,:,:]-numpy.median( images[i,:,:], axis=None ) *master_mask    
                    log.debug("------------------>MEDIAN (despues) = %s",  numpy.median( images[i,:,:],axis=None ))           
                
                sky_bg = numpy.median( images, axis=0 )
                        
            images = 0.0  # to free memory ??
            log.debug("SKY background is ---> %s", numpy.mean(sky_bg))
             
            # modify/update images subtracting sky background
            for i in range(0,nimages):
                images_[i,:,:] = images_[i,:,:]-sky_bg
                 
            # Write out to file sky_bg
            if output_file!=None:
                print "SHAPE= %s", sky_bg.shape   
                misc.fileUtils.removefiles(output_file)
                hdu = pyfits.PrimaryHDU(sky_bg)
                hdu.scale('float32') # importat to set first data type
                #hdu.data=sky_bg     
                hdulist = pyfits.HDUList([hdu])
                hdu.header.update('OBJECT','SKY_BG')
                hdulist.writeto(output_file)
                hdulist.close(output_verify='ignore')  
      
    def M_backgroundSubtraction(self, images_in):
        """DESCRIPTION
                Compute image background level (sky) and subtract it
                
           INPUTS
                image_in      two dimensional array
                      
           OUTPUTS
        """
            
        # TO BE DONE !!!!!!!
            
    
    def M_cosmicRemoval(self, images_in=None):
        """DESCRIPTION
                Fast cosmic ray removal using a sigma filter 
                
           INPUTS
                image_in      two dimensional array
                      
           OUTPUTS
        """
            
        # TO BE DONE !!!!!!!
            
        if images_in==None:
            images = self.m_images_in
        else:
            images = images_in
                
        n_sigma = 3.5
        box_width = 5        
        for i in range(0, self.m_nfiles):
            images[i,:,:] = scipy.ndimage.filters.median_filter( images[i,:,:], box_width )#scipy.ndimage.filters.gaussian_filter(images[i,:,:], 3) 
            #scipy.filters.median_filter( images[i,:,:], box_width )
            
            
            
                
    def M_getPointingOffsets (self, images_in=None, type='sex', p_min_area=5, p_mask_thresh=2):
        """DESCRIPTION
                Derive pointing offsets between each image using SExtractor OBJECTS (type=sex) or only header info (type=head)
                
           INPUTS
                images_in       cube of two dimensional arrays with images
                 
                type            define which type of method will be used to find out the offsets. It can be 'sex' or 'head'
                 
                p_min_area      minimun area for a detected object
                 
                p_mask_thresh   threshold for object detection
                      
           OUTPUTS
                offsets         two dimensional array with offsets
                
           NOTE
                Assumes that North is up and East is left
        """
        if images_in==None:
            images = self.m_images_in
        else:
            images = images_in
            
        if type=='sex':
            # STEP 1: Write out memory images to FITS files 
            output_dir=self.m_output_dir
            suffix='_gpo.fits'
            self.writeOutData(images, output_dir, suffix)
                
            # STEP 2: Create SExtractor OBJECTS images (call external app-> offsets)
            # First, check if already OBJECTS images are created
            if not os.path.exists( output_dir+"/mose_objs.pap" ):
                log.debug("Creaing OBJECTS images ....")
                min_area=p_min_area #SExtractor DETECT_MINAREA
                mask_thresh=p_mask_thresh #SExtractor DETECT_THRESH
                output_list_file=output_dir+"/gpo_objs.pap"
                makeObjMask( output_dir+'*'+suffix , min_area, mask_thresh, output_list_file)
            else:
                output_list_file=output_dir+"/mose_objs.pap"
                log.debug("BIENNNNNN no tengo que hacer nadaaaaaaaaaaaa!!!!")
                                            
                            
            # STEP 3: Compute dither offsets (in pixles) using cross-correlation technique ==> offsets
            #>mosaic objfiles.nip $off_err > offsets1.nip
            search_box=10 # half_width of search box in arcsec (default 10)
            offsets_file='/tmp/offsets.pap'
            offsets_cmd=self.m_papi_path+'/irdr/bin/offsets '+ output_list_file + '  ' + str(search_box) + ' >'+offsets_file
            if utils.runCmd( offsets_cmd )==0:
                log.error("Some error while computing dither offsets ....")
                return
            
            # STEP 4: Read back computed offsets from file into memory array-matrix with a row for each file offsets
            # Row 0 is corresponding to the reference frame  
            try:
	            m_offsets = numpy.loadtxt(offsets_file, usecols = (1,2,3)) # columns => (xoffset, yoffset, match fraction) in PIXELS  
            except IOError:
	            log.debug("Any offsets read. Check images ....")     
        
        elif type=='head':
            log.error("sorry, but not yet implemented ....")
        else:
            log.error("Option not recognized")
                        
    def I_getPointingOffsets (self, images_in=None,  p_min_area=5, p_mask_thresh=2, offsets_file='/tmp/offsets.pap'):
        """DESCRIPTION
                Derive pointing offsets between each image using SExtractor OBJECTS (makeObjMask) and offsets (IRDR)
                
           INPUTS
                images_in        filename of list file
                 
                p_min_area       minimun area for a detected object
                 
                p_mask_thresh    threshold for object detection
                      
           OUTPUTS
                offsets          two dimensional array with offsets
                
           NOTE
                Assumes that North is up and East is left
        """
            
        log.info("Starting I_getPointingOffsets....")
            
        # STEP 1: Create SExtractor OBJECTS images (call external app offsets.c)
        output_dir=self.m_output_dir
        suffix='_'+self.m_filter+'.skysub.fits'
        log.debug("Creaing OBJECTS images (SExtractor)....")
        min_area=p_min_area #SExtractor DETECT_MINAREA
        mask_thresh=p_mask_thresh #SExtractor DETECT_THRESH
        output_list_file=output_dir+"gpo_objs.pap"
        
        if images_in==None:
            makeObjMask( output_dir+'*'+suffix , min_area, mask_thresh, output_list_file)
        elif os.path.isfile(images_in):
            makeObjMask( images_in , min_area, mask_thresh, output_list_file)
        else:
            log.error("Option not recognized !!!")    
                            
        # STEP 2: Compute dither offsets (in pixles) using cross-correlation technique ==> offsets
        #>mosaic objfiles.nip $off_err > offsets1.nip
        search_box=10 # half_width of search box in arcsec (default 10)
        offsets_cmd=self.m_papi_path+'/irdr/bin/offsets '+ output_list_file + '  ' + str(search_box) + ' >'+offsets_file
        if utils.runCmd( offsets_cmd )==0:
            log.error ("Some error while computing dither offsets")
            return
        
        # STEP 3: Read back computed offsets from file into memory array-matrix with a row for each file offsets
        # Row 0 is corresponding to the reference frame  
        try:
            m_offsets = numpy.loadtxt(offsets_file, usecols = (1,2,3)) # columns => (xoffset, yoffset, match fraction) in PIXELS  
        except IOError:
            log.debug("Any offsets read. Check images ....")     
    
        log.debug("END of I_getPointingOffsets")
                 
                   
    def M_sky(self, image, skymode, skysigma):
        """DESCRIPTION
                Compute image sky level using the 'sky' task of the STSDAS Dither package to estimate sky value
                
           INPUTS
                image      two dimensional array
                      
           OUTPUTS
                skymode    Mode of the image (sky background)
                skysigma   Sigma of the sky mode        
        """
            
        # TO BE DONE !!!!!!!
        
        # set all iraf'y setup stuff
        iraf.flpr();iraf.flpr()
        iraf.stsdas()
        iraf.analysis()
        iraf.dither()
        iraf.unlearn(iraf.dither.sky)
        iraf.reset(imtype='fits')
        iraf.set(tmp='./')
        
        # set 'sky' task params
        iraf.dither.sky.input = simplefits
        iraf.dither.sky.masks = ''
        iraf.dither.sky.lower = skyest - 2.4*rms
        # iraf.dither.sky.upper = skyest + 1.5*rms  -- 12/Apr/2002,jpb
        iraf.dither.sky.upper = skyest + 1.6*rms
        iraf.dither.sky.dq = ''
        iraf.dither.sky.subsky = iraf.no         # do not subtract sky!
        # iraf.dither.sky.width = 1.9*rms  -- 12/Apr/2002,jpb
        iraf.dither.sky.width = 2.0*rms
        iraf.dither.sky.bunit = 'counts'
        iraf.dither.sky.expname = ''
        iraf.dither.sky.skyname = ''
        iraf.dither.sky.verbose = iraf.yes
        iraf.dither.sky.tempdir = './'
        iraf.dither.sky.skyvalue = 0.
        iraf.dither.sky.mode = 'h'
        
        # run sky using mode estimator
        iraf.dither.sky.stat = 'mode'
        iraf.dither.sky()
        sky_mode = iraf.dither.sky.skyvalue
        GlobalBlock.logfile.write(' IRAF dither.sky mode estimate: '+str(sky_mode))
        
        if irafskytoo == 2:
            # run sky using mean estimator
            iraf.dither.sky.stat = 'mean'
            iraf.dither.sky()
            sky_mean = iraf.dither.sky.skyvalue
            GlobalBlock.logfile.write(' IRAF dither.sky mean estimate: '+str(sky_mean))
            skyraf = min(sky_mean,sky_mode)
        else:
            skyraf = sky_mode
        
        return skyraf
  
      
    def M_maskObjects (self, masks, nsigma=0 ):
        """DESCRIPTION
                Find and create a (cube) object mask of the m_images_in
           
           INPUTS
                nsigma    number of sigmas for thresholding      
           OUTPUTS
                masks     Cube of objects mask 
              
        """
                                  
        """PENDIENTE !!!!!!!!!!! TO BE COMPLETED !!!!!!!! ====> SEE M_maskObjectsSE !!!!! """                      
         
        if nsigma==0:
            nsigma=3.0
                                  
        masks = numpy.ones(self.m_images_in.shape, dtype=numpy.int)
        mask  = numpy.ones(self.m_images_in.shape[1], self.m_images_in.shapes[2], dtype=numpy.float)
        
        for i in range(0,self.m_images_in.shape[0]):
            mask[:,:]=1.0
            image = self.m_images_in[i,:,:]
            
            A = numpy.isfinite(image)
            if numpy.where(A==True)[0].shape[0] > 30: # A.sum()>30
                self.M_sky(image, bg, rms)
                image_smooth = smooth(image,3)
                thres = nsigma*rms/3.0  # 3-pix square aperture, n sigma
                A = numpy.where(image_smooth>bg+thres)
                if A.sum() !=-1:
                    mask[A] = 0.0
                    
                mask = numpy.median(mask, axis=None)
                masks[i,:,:] = mask.astype(int)
                
         
    def M_maskObjectsSE (self, masks, images_in=None, p_min_area=10, p_mask_thresh=5.0, remove_objs=False ):
        """DESCRIPTION
                Find and create a (cube) object mask of the m_images_in based on SExtractor
           
           INPUTS
                images_in       images to look for objects
                
                nsigma          number of sigmas for thresholding
           
           OPTIONAL INPUTS
                p_min_area      SExtractor DETECT_MINAREA (minimun object area)
                           
                p_mask_thresh   SExtractor DETECT_THRESH
                
                remove_objs     Flat to indicate if remove or not tha .obj files created after masks creation
                
           OUTPUTS
                masks           Cube of objects mask 
              
        """
                                  
        log.info('Start M_maskObjectSE...')
                                  
        if images_in==None:
            images = self.m_images_in
        else:
            images = images_in
            
        masks = numpy.ones(images.shape, dtype=numpy.int)
        #mask  = numpy.ones(images.shape[1], self.m_images_in.shapes[2], dtype=numpy.float)                          
                                  
        #STEP 1: Write out current images to FITS files
        output_dir=self.m_output_dir
        suffix='_mose.fits'
        # Remove old files
        for file in glob.glob(output_dir+"/*_mose.fits"):
            os.remove(file)
        self.writeOutData(images, output_dir, suffix)
        
        #STEP 2: Create Object masks
        min_area=p_min_area #SExtractor DETECT_MINAREA
        mask_thresh=p_mask_thresh #SExtractor DETECT_THRESH
        output_list_file=output_dir+"/mose_objs.pap"
        makeObjMask( output_dir+'/*'+suffix , min_area, mask_thresh, output_list_file)            
                    
        #STEP 3: Read .fits.obj files and dump into mask matrix cube
        for i in range(0, self.m_nfiles):
            f=pyfits.open((output_dir+self.m_cl_file_list[i].filename).replace(".fits","_mose.fits.objs"))
            # Sextractor mask definition      ==> >0 in the pixel represent object pixel, and 0 represent sky pixel 
            # But now, our mask definition is ==> 1 in the mask represent sky data, 0 represent object or bad data
            masks[i,:,:]=numpy.where(f[0].data>0, 0, 1)
            f.close()
            
            if (0):# only for debug
                output='/tmp/mask_1.fits'
                misc.fileUtils.removefiles(output)
                hdu = pyfits.PrimaryHDU(masks[1,:,:])
                hdu.scale('int32') # importat to set first data type
                hdulist = pyfits.HDUList([hdu])
                hdu.header.update('OBJECT','MASK')
                hdulist.writeto(output)
                hdulist.close(output_verify='ignore')
                log.debug("New file created: %s ", output)    
        
            #STEP 4: Remove temporary files
            if remove_objs:
                misc.fileUtils.removefiles((output_dir+self.m_cl_file_list[i].filename).replace(".fits","_mose.fits.objs"))
                                              
    def M_coaddCubeImages(self, gain="/tmp/gain.fits", output="/tmp/coadd.fits"):
        """  Coadd the cube/stack of dithered FITS images using dithercubemean from IRDR"""
                                                  
        prog = self.m_papi_path+"/irdr/bin/dithercubemean"
        gain_file=gain
        output_file=output
        weight_file=output_file.replace(".fits",".weight.fits")
        
        # STEP 1: Dump the images to FITS files
        self.writeOutData(self.m_images_in, "/tmp/", None, None)
        
        # STEP 2: Create file list file with : filename + dither offsets values computed
        f_img = open(images_file)
        log.error("TO BE DONE !!!!!")
        
        # STEP 2: Run the coadd                                          
        cmd  = prog + " " + images_file + " " + gain_file + " " + output_file + " " + weight_file 
        e=utils.runCmd( cmd )
        if e==0:
            log.debug("Some error while running command %s", cmd)
        else:
            log.debug("Succesful ending of M_coaddCubeImages")
                                                      
    def I_coaddStackImages(self, input='/tmp/stack.pap', gain='/tmp/gain.fits', output='/tmp/coadd.fits'):
        """  Coadd the stack of dithered FITS images listed in 'input_file' using dithercubemean from IRDR"""
                                                  
        log.info("Start I_coaddStackImages ...")                                          
        # STEP 1: Define parameters                                          
        prog = self.m_papi_path+"/irdr/bin/dithercubemean "
        input_file=input
        if input_file==None:
            log.error("Bad input file provided !")
            return
        if gain==None:
            gain_file=self.m_output_dir+"/gain_"+self.m_filter+".fits"
        else:
            gain_file=gain
        if output==None:
            output_file=self.m_output_dir+"/coadd_"+self.m_filter+".fits"     
        else:
            output_file=output
            
        weight_file=output_file.replace(".fits",".weight.fits")
        
        # STEP 2: Run the coadd (dithercubemean)                                          
        cmd  = prog + " " + input_file + " " + gain_file + " " + output_file + " " + weight_file 
        e=utils.runCmd( cmd )
        if e==0:
            log.debug("Some error while running command %s", cmd)
        else:
            log.debug("Succesful ending of coaddStackImages")
    
    def I_createMasterObjMask( self, input_file, output_master_obj_mask ):
        """ Create a master object mask from a input file using SExtractor and then dilate the object file by a certain factor to remove also the  undetected object tails (dilate.c)"""
                                                             
                                                             
        log.info("Start I_createMasterObjMask....")
                                                             
        # STEP 1: create mask                                                     
        makeObjMask( input_file+"*", 5, 2.0)
        if os.path.exists(input_file+".objs"): shutil.move(input_file+".objs", output_master_obj_mask)
        # STEP 2: dilate mask
        prog = self.m_papi_path+"/irdr/bin/dilate "
        scale = 0.5 #mult. scale factor to expand object regions; default is 0.5 (ie, make 50%% larger)
        cmd  = prog + " " + input_file + " " + str(scale)  
        """e=utils.runCmd( cmd )
        if e==0:
            log.debug("Some error while running command %s", cmd)
        else:
            log.debug("Succesful ending of createMasterObjMask")
        """                                                             
                                                                 
    def M_showStats( self ):
        """ Show some stats of the images ..."""
                                                
        for i in range(0, self.m_nfiles):
            #mode=scipy.stats.mstats.mode(self.m_images_in[i,:,:],axis=None)[0].data
            mean=numpy.mean(self.m_images_in[i,:,:],axis=None)
            median=numpy.median(self.m_images_in[i,:,:],axis=None)
            min=numpy.min(self.m_images_in[i,:,:],axis=None)
            max=numpy.max(self.m_images_in[i,:,:],axis=None)
            log.debug("STATS of image %s: mean-> %s  |median-> %s |min-> %s | max-> %s ", self.m_cl_file_list[i].filename, mean, median, min, max) 
                                                        
                                                
                            
    def startReduction(self, type='double', mosaic='dither', dark=None, ext_flat=None):
        """ 
        \par Description:
            Main routine to start PANIC reduction of a frame set of science frames (point or extended objects)         
        \input type
            Specify what type of reduction will be done (single, double, full) 
        \param mosaic
            Specify the type of observation is going to be reduced (dithered sequence or object-sky sequence)
        \param dark
            =1 if we want dark subtraction
        \param ext_flat
            =1 if we provide an external (dome) flat, then we won't use super sky flat    
        \retval 0
            If no error
        \author
            JMIbannez, IAA-CSIC
        """
     
        log.info("Start Data Reduction !!!")
        
        # Clean old files
        for file in glob.glob(self.m_output_dir+"*.objs"): os.remove(file)
        for file in glob.glob(self.m_output_dir+"*.pap"): os.remove(file)
        
        # Read main parameters
        self.read_red_parameters()
        # Check data files are ok (filter, date, exists, etc ...) and then read data into 3D-matrix 'images_in'
        if self.checkFilter and self.checkType():
            # DARK subtraction
            if dark!=None:
                self.M_subtracDark(dark)
            
            # FLAT-FIELDING
            if ext_flat!=None:
                # External master flat (dome-flat)      
                self.M_appFlatField(ext_flat, self.m_images_in)
                self.M_appChipOffset(self.m_images_in)
                #c_images_in = numpy.copy(self.m_images_in)
                #self.M_subtracMedSky(c_images_in, 0, None,'/tmp/sky_bg_1.fits') # ojo, no tiene en cuenta las N-nearest
                #masks=None
                #self.M_maskObjectsSE(masks, c_images_in, 15, 7.5)
                #self.M_subtracMedSky(self.m_images_in, 0, masks, '/tmp/sky_bg_2.fits')
                self.M_getPointingOffsets(self.m_images_in, 'sex', 15, 7)
                self.M_coaddCubeImages()
                self.M_showStats()
            else:
                # SNAP steps
                
                # STEP 1: Create super flat (or gain map)
                self.I_createSuperFlat(self.m_images_in, None, None)
                
                # STEP 2: Subtract sky and divide by super flat (gain map)
                self.MI_skyFilter ('/tmp/sf_KS.pap', 'file', None, 'nomask')
                
                # STEP 3: Compute dither offsets from the first sky subtracted/filtered images using cross-correlation
                self.I_getPointingOffsets (None, 5, 2, '/tmp/offsets1.pap')
                
                # STEP 4: Create first coadded image of the dithered stack using preliminary offsets
                # 4.1 :create the file of sky subtracted images with the offsets computed
                fo=open('/tmp/offsets1.pap')
                fs=open('/tmp/stack1.pap','w')
                for line in fo:
                    n_line = line.replace(".skysub.fits.objs", ".skysub.fits") 
                    fs.write( n_line )
                fo.close()
                fs.close()    
                self.I_coaddStackImages('/tmp/stack1.pap', None, '/tmp/coadd1.fits')
                
                
                # STEP 5: Create master object mask
                self.I_createMasterObjMask('/tmp/coadd1.fits', '/tmp/masterObjMask.fits')
                
                # STEP 6: For each image, a new sky image is computed and subtracted by taking into account the new master object mask
                # 6.1: create the 
                f_temp=open('/tmp/offsets1.pap')
                f_out_skf=open('/tmp/skyfilter_2.pap','w')
                for line in f_temp:
                    n_line = line.split()[0].replace(".skysub.fits.objs", ".fits") + " "+ "/tmp/masterObjMask.fits" + " " + line.split()[1] + " " + line.split()[2] +"\n"
                    #print "LINE=",n_line
                    f_out_skf.write(n_line)
                f_out_skf.close()
                f_temp.close()         
                self.MI_skyFilter( '/tmp/skyfilter_2.pap', 'file' , None, 'mask')
                
                # STEP 7:Create second coadded image of the dithered stack using new sky subtracted frames
                self.I_coaddStackImages('/tmp/stack1.pap', None, '/tmp/coadd2.fits')
                imgTrim("/tmp/coadd2.fits")
                
                """ SIMPLE-WIRCAM steps
                #Super-Flat-Field
                super_flat='/tmp/super_flat.fits'
                self.M_createSuperFlat(self.m_images_in, super_flat, mask=None, fast=0)
                c_images_in = numpy.copy(self.m_images_in)
                self.M_appFlatField( super_flat, c_images_in)
                masks=None
                self.M_maskObjectsSE( masks, c_images_in, 15, 1.5)
                self.M_createSuperFlat(self.m_images_in, super_flat, masks, fast=0)
                self.M_appFlatField( super_flat, self.m_images_in)
                self.M_appChipOffset(self.m_images_in)
                """
            
            #self.writeOutData()
        
################################################################################
         
class ExAbort(Exception): 
  """Module-specific error, to be raised when an the reduction is to be aborted
     in a friendly way. This error should be caught by the task manager, who
     then stops the processing.
  """
  def __init__(self, *args):
    self.args = args


class ExError(Exception): 
  """Module-specific error, to be raised when an 'expected' exception occurs
     during the processing of a task. This error should then be caught by the
     task manager, avoiding nasty error output.
  """
  def __init__(self, *args):
    self.args = args

        
#################################################################################
### MAIN
#################################################################################        
if __name__ == "__main__":
    print 'main: reduce module ....'
    # Get and check command-line options
    args = sys.argv[1:]
    st ='/tmp/f1.fits , /tmp/f2.fits , /tmp/f3.fits , /tmp/f4.fits , '
    lis = ['/tmp/f5.fits','/tmp/f4.fits', '/tmp/f1.fits','/tmp/f2.fits', '/tmp/f3.fits', '/tmp/f6.fits', '/tmp/f7.fits']
    lis2 = ['/disk-a/caha/panic/DATA/prueba1/orion0020.fits', '/disk-a/caha/panic/DATA/prueba1/orion0021.fits', '/disk-a/caha/panic/DATA/prueba1/orion0022.fits']
    lis3 = ['/disk-a/caha/panic/DATA/ALHAMBRA_1/A0408160047.fits',
    '/disk-a/caha/panic/DATA/ALHAMBRA_1/A0408160048.fits',
    '/disk-a/caha/panic/DATA/ALHAMBRA_1/A0408160049.fits',
    '/disk-a/caha/panic/DATA/ALHAMBRA_1/A0408160050.fits', 
    '/disk-a/caha/panic/DATA/ALHAMBRA_1/A0408160051.fits', 
    '/disk-a/caha/panic/DATA/ALHAMBRA_1/A0408160052.fits', 
    '/disk-a/caha/panic/DATA/ALHAMBRA_1/A0408160053.fits', 
    '/disk-a/caha/panic/DATA/ALHAMBRA_1/A0408160054.fits',
    '/disk-a/caha/panic/DATA/ALHAMBRA_1/A0408160055.fits',
    '/disk-a/caha/panic/DATA/ALHAMBRA_1/A0408160056.fits',
    '/disk-a/caha/panic/DATA/ALHAMBRA_1/A0408160057.fits',
    '/disk-a/caha/panic/DATA/ALHAMBRA_1/A0408160058.fits',
    '/disk-a/caha/panic/DATA/ALHAMBRA_1/A0408160059.fits',
    '/disk-a/caha/panic/DATA/ALHAMBRA_1/A0408160060.fits']
    lis4 = ['/disk-a/caha/panic/DATA/ALHAMBRA_1/A0408160047.fits',
    '/disk-a/caha/panic/DATA/ALHAMBRA_1/A0408160048.fits',
    '/disk-a/caha/panic/DATA/ALHAMBRA_1/A0408160050.fits']
    lis5= ['/disk-a/caha/panic/DATA/SIMU_PANIC_2/SI/ferM_0001.fits','/disk-a/caha/panic/DATA/SIMU_PANIC_2/SI/ferM_0002.fits', '/disk-a/caha/panic/DATA/SIMU_PANIC_2/SI/ferM_0003.fits', '/disk-a/caha/panic/DATA/SIMU_PANIC_2/SI/ferM_0004.fits', '/disk-a/caha/panic/DATA/SIMU_PANIC_2/SI/ferM_0005.fits']
     
    #lis2 = ['/disk-a/caha/panic/DATA/prueba1/orion0020.fits']   
    #ret1 = stringToList(st)
    #ret2 = listToString(lis)
    
    # Tests for Reduction Block
    #rb = ReductionBlock (lis2)
    #rb.simpleReduce(6, 4,  '/tmp/p1.fits', '/tmp/master_dflat.fits', '/tmp/BadPixMask-20081127181927.pl') #, '/tmp/bpm.pl', do_astrometry=False)
    
    #rb.createBadPixelMask()
    #rb.sortByDateObs()
    #rb.computeSkyBackground_i(3, 0,'-BACKGROUND')
    """
    for i in range(len(lis)):
        for w in range(3):
            print "*******  I=%d W=%d  ***************" %(i,w+1)
            rb.computeSkyBackground_i(i,w+1)
       """     
    """rb.computeSkyBackground_i(0,3)
    rb.computeSkyBackground_i(1,3)
    rb.computeSkyBackground_i(2,3)
    rb.computeSkyBackground_i(3,3)
    rb.computeSkyBackground_i(4,3)
    rb.computeSkyBackground_i(5,3)
    """
    
    sr=SimpleReduce()
    sr.run("/disk-a/caha/panic/DEVELOP/PIPELINE/PANIC/PAPI/simu_1/p4/ferM_0032.fits", "/tmp/master_dark.fits", '/tmp/master_sflat.fits', "/tmp/result.fits", '/tmp/',False)
    
        
    """
    rb = ReductionBlock(lis5)
    rb.getMean()
    #rb.M_createSuperFlat('/tmp/super_flat1.fits', None, 0)
    #rb.M_subtracMedSky(0, None,'/tmp/sky_bg.fits')
    #rb.startReduction('double','dither','/tmp/master_dark.fits','/tmp/master_sflat.fits')
    rb.startReduction('double','dither',None, None)
    rb.getMean()
    """
    
    
    
    
    
