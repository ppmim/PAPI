#!/usr/bin/env python

################################################################################
#
# PAPI (PAnic PIpeline)
#
# papi.py
#
# Last update 04/March/2010
#
################################################################################


    ####################################################################
    #                                                                  #
    #                   This is where it all starts.                   #
    #         The body of this program is at the end of this file      #
    #                                                                  #
    ####################################################################

__date__ = "$Date: 2010-03-05 11:48:08 +0100 (Fri, 05 Mar 2010) $"
__author__ = "$Author: panic $"
__revision__ = "$Rev: 21 $"

    
#From system
import sys
import os
import os.path
import fnmatch
import time
from optparse import OptionParser
import fileinput
import glob
import shutil

# Interact with FITS files
import pyfits

# IRAF packages
from pyraf import iraf
from iraf import noao
from iraf import imred
from iraf import ccdred

# Math module for efficient array processing
import numpy

#Log
import misc.paLog
from misc.paLog import log    

#PAPI packages 
import datahandler
import reduce
import misc.fileUtils
import PAPI.linkSources as papi
import misc.utils
from reduce.makeobjmask import *
import reduce.imtrim


class ReductionSet:
    def __init__(self, list_file, out_dir, out_file, obs_mode, dark=None, flat=None, bpm=None):
        """ Init function """
        
        # Input values
        self.list_file = list_file # a file containing a list of data filenames 
        self.out_dir = out_dir     # directory where all output will be written
        self.out_file = out_file   # final reduced data file (out)
        self.obs_mode = obs_mode   # observing mode (dither, dither_on_off, dither_off_on....)
        self.master_dark = dark    # master dark to use (input)
        self.master_flat = flat    # master flat to use (input)
        self.bpm = bpm             # master Bad Pixel Mask to use (input)
        
        # Environment variables
        self.m_terapix_path = os.environ['TERAPIX']
        self.m_papi_path = os.environ['PAPI_HOME']
        
        
        self.m_LAST_FILES = []   # Contain the files as result of the last processing step
        self.m_rawFiles = []     # Raw files (originals in the working directory)
        self.m_filter = ""
        self.m_n_files = ""
        
        
        
    def checkData(self):
        """ This function makes some test with the data, previusly to start the reduction"""
                
        # check files currently in self.m_LAST_FILES
        if not self.checkFilter():
            raise Exception("Error while checking data")
        
        # check type (does not work with dither_on_off/dither_off_on)
        #if not self.checkType():
        #    raise Exception("Error while checking data")
          
        
    def checkFilter(self):
        """Return true is all files in file have the same filter type, false otherwise
        
        \return True or False
        """
        f = datahandler.ClFits(self.m_LAST_FILES[0])
        filter_0 = f.getFilter()
        self.m_filter = filter_0
        for file in self.m_LAST_FILES:
            fi=datahandler.ClFits( file )
            if fi.getFilter() != filter_0:
                log.debug("File %s does not match file filter", file)
                return False
            
        log.debug("All files match same file filter")
        return True
        
        
    def checkType(self, type_to_check=None):
        """
        Return true is all files in file have the same type(science, dark, flat, ...), false otherwise
        
        \param type_to_check (\c string) type to check to
        \return True or False
        """
        
        if type_to_check==None:
            f = datahandler.ClFits(self.m_LAST_FILES[0])
            type_0 = f.getType()
        else:
            type_0 = type_to_check
             
        for file in self.m_LAST_FILES:
            fi = datahandler.ClFits(file)
            if fi.getType() != type_0:
                log.debug("File %s does not match file type %s", file, type_0)
                return False
            
        log.debug("All files match same file type")
        return True
                 
    def sortOutData(self, list=None):
        """
        Sort out input data files by MJD
        """
        
        dataset=[]
        if list==None:
            m_list=self.m_LAST_FILES
        else:
            m_list=list
            
        for file in m_list:
            fits=datahandler.ClFits(file) 
            dataset.append((file, fits.getMJD()))
            
        dataset=sorted(dataset, key=lambda data_file: data_file[1])          
        sorted_files=[]
        for tuple in dataset:
            sorted_files.append(tuple[0])
        
        return sorted_files                                
                     
    def getObsMode(self):
        """
        Return the type of dither sequece followed in the currect  'm_LAST_FILES' list. It could be:
            - dither (T-T-T-T-T- ....)
            - dither_on_off (T-S-T-S-T-S-T-S-T-....)
            - dither_off_on (S-T-S-T-S-T-S-T-S-....)
            - other
        """
                   
        mode='other' # default
                   
        #Step 1: we suppose list file is sorted out  by MJD
        #Step 2: get the data type of the first file and check sequence starting from this file
        fits_0=datahandler.ClFits(self.m_LAST_FILES[0])
        fits_1=datahandler.ClFits(self.m_LAST_FILES[1])
        i=0
        if fits_0.isSky() and fits_1.isObject():
            mode='dither_off_on'
            # Then, we are going to suppose the sequence S-T-S-T-S- .... (dither_off_on)
            for file in self.m_LAST_FILES:
                if not i%2: #even
                    fits=datahandler.ClFits(file)
                    if not fits.isSky():
                        return 'other'
                elif i%2: #odd
                    fits=datahandler.ClFits(file)
                    if not fits.isObject():
                        return 'other'
                i=i+1         
        elif fits_0.isObject() and fits_1.isSky():
            # Then, we are going to suppose the sequence T-S-T-S-T- .... (dither_on_off)
            mode='dither_on_off'
            for file in self.m_LAST_FILES:
                if not i%2: #even
                    fits=datahandler.ClFits(file)
                    if not fits.isObject():
                        return 'other'
                elif i%2: #odd
                    fits=datahandler.ClFits(file)
                    if not fits.isSky():
                        return 'other'
                i=i+1                 
        elif fits_0.isObject() and fits_1.isObject():
            # Then, we are going to suppose the sequence T-T-T-T-T- .... (dither)
            mode='dither'
            for file in self.m_LAST_FILES:
                fits=datahandler.ClFits(file)
                if not fits.isObject():
                    return 'other'
        else:
            mode='other'
                        
        return mode
                        
    def getSkyFrames(self, list=None):
        """
        Return the files identified as 'sky' frames in the m_LAST_FILES
        """
       
        sky_list=[]
        if list==None:
            m_list=self.m_LAST_FILES
        else:
            m_list=list
            
        for file in m_list:
            fits=datahandler.ClFits(file)
            if fits.isSky():
                sky_list.append(file)
                
        return sky_list         
                         
    def skyFilter( self, list_file, gain_file, mask='nomask', obs_mode='dither' ):
        """
            For each input image, a sky frame is computed by combining a certain number of the closest images, 
            then this sky frame is subtracted to the image and the result is divided by the master flat; 
                         
            This function is a wrapper for skyfilter.c (IRDR)              
            
            INPUT
                list_file : a text file containing the suited structure in function of the observing_mode 
            
            OUTPUT
                The function generate a set of sky subtrated images (*.skysub.fits)             
            
            VERSION
                1.0, 20090909 by jmiguel@iaa.es
        
            TODO: extended objects !!!!
        """               
        
        # Skyfilter parameters
        halfnsky=4
        destripe='none'
        out_files=[]
            
        if obs_mode=='dither':
            skyfilter_cmd=self.m_papi_path+'/irdr/bin/skyfilter '+ list_file + '  ' + gain_file +' '+ str(halfnsky)+' '+ mask + '  ' + destripe 
        elif obs_mode=='dither_on_off':
            skyfilter_cmd=self.m_papi_path+'/irdr/bin/skyfilteronoff '+ list_file + '  ' + gain_file +' '+ str(halfnsky)+' '+ mask + '  ' + destripe
        elif obs_mode=='dither_off_on':
            skyfilter_cmd=self.m_papi_path+'/irdr/bin/skyfilteroffon '+ list_file + '  ' + gain_file +' '+ str(halfnsky)+' '+ mask + '  ' + destripe
        else:
            log.error("Observing mode not supported")
            raise
                  
        if misc.utils.runCmd( skyfilter_cmd )==1: # All was OK
            # Rename output sky-subtracted files
            for file in glob.glob(self.out_dir+'/*.fits.skysub'):
                shutil.move(file, file.replace('.fits.skysub', '.skysub.fits'))
                #out_files.append(file.replace('.fits.skysub', '.skysub.fits'))
            out_files=[line.split()[0].replace('.fits', '.skysub.fits') for line in fileinput.input(list_file)]
            if (obs_mode=='dither_on_off' or obs_mode=='dither_off_on'):
                out_files=glob.glob(self.out_dir+'/*.skysub.fits')
                out_files=self.sortOutData(out_files) 
            return out_files
        else: return []
                                  
    def getPointingOffsets (self, images_in=None,  p_min_area=5, p_mask_thresh=2, p_offsets_file='/tmp/offsets.pap'):
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
            
        log.info("Starting getPointingOffsets....")
        offsets_mat=None
           
        # STEP 1: Create SExtractor OBJECTS images
        suffix='_'+self.m_filter+'.skysub.fits'
        #p_min_area - SExtractor DETECT_MINAREA
        #p_mask_thresh - SExtractor DETECT_THRESH
        output_list_file=self.out_dir+"/gpo_objs.pap"
        
        log.debug("Creaing OBJECTS images (SExtractor)....")

        if images_in==None:
            makeObjMask( output_dir+'*'+suffix , p_min_area, p_mask_thresh, output_list_file)
        elif os.path.isfile(images_in):
            makeObjMask( images_in , p_min_area, p_mask_thresh, output_list_file)
        else:
            log.error("Option not recognized !!!")    
                           
        # STEP 2: Compute dither offsets (in pixles) using cross-correlation technique ==> offsets
        #>mosaic objfiles.nip $off_err > offsets1.nip
        search_box=10 # half_width of search box in arcsec (default 10)
        offsets_cmd=self.m_papi_path+'/irdr/bin/offsets '+ output_list_file + '  ' + str(search_box) + ' >' + p_offsets_file
        if misc.utils.runCmd( offsets_cmd )==0:
            log.error ("Some error while computing dither offsets")
            return 0
        else:
            try:
                offsets_mat = numpy.loadtxt(p_offsets_file, usecols = (1,2,3)) # columns => (xoffset, yoffset, match fraction) in PIXELS  
            except IOError:
                log.debug("Any offsets read. Check images ....")   
        
        log.debug("END of getPointingOffsets")                        
        return offsets_mat
                            
    def coaddStackImages(self, input='/tmp/stack.pap', gain='/tmp/gain.fits', output='/tmp/coadd.fits'):
        """  Coadd the stack of dithered FITS images listed in 'input_file' using dithercubemean from IRDR"""
                                                  
        log.info("Start coaddStackImages ...")                                          
        # STEP 1: Define parameters                                          
        prog = self.m_papi_path+"/irdr/bin/dithercubemean "
        input_file=input
        if input_file==None:
            log.error("Bad input file provided !")
            return
        
        if gain==None:
            gain_file=self.out_dir+"/gain_"+self.m_filter+".fits"
        else:
            gain_file=gain
        if output==None:
            output_file=self.out_dir+"/coadd_"+self.m_filter+".fits"     
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
                                    
    def createMasterObjMask( self, input_file, output_master_obj_mask ):
        """ Create a master object mask from a input file using SExtractor and then dilate the object file by a certain factor to remove also the  undetected object tails (dilate.c)"""
                                                             
                                                             
        log.info("Start createMasterObjMask....")
                                                             
        # STEP 1: create mask                                                     
        makeObjMask( input_file+"*", 5, 2.0)
        if os.path.exists(input_file+".objs"): 
            shutil.move(input_file+".objs", output_master_obj_mask)
            log.debug("New Object mask created : %s", output_master_obj_mask)
            
        # STEP 2: dilate mask
        log.info("Dilating image ....(NOT DONE by the moment)")
        prog = self.m_papi_path+"/irdr/bin/dilate "
        scale = 0.5 #mult. scale factor to expand object regions; default is 0.5 (ie, make 50%% larger)
        cmd  = prog + " " + input_file + " " + str(scale)
          
        return output_master_obj_mask
        
        """e=utils.runCmd( cmd )
        if e==0:
            log.debug("Some error while running command %s", cmd)
        else:
            log.debug("Succesful ending of createMasterObjMask")
        """        
                                        
    def makeAstrometry( self, input_file, catalog='2mass', re_grid=False):
        """Compute the initial astrometry of the given input_file"""
                                            
        if re_grid: regrid_str='regrid'
        else: regrid_str='noregrid'
                                            
        astrometry_cmd=self.m_papi_path+'/astrometry_scamp.pl '+ catalog + '  ' + regrid_str + ' ' + input_file
        if misc.utils.runCmd( astrometry_cmd )==0:
            log.error ("Some error while computing Astrometry")
            return 0
        else:
            return 1                                    
                                            
    def reduce(self):
        """ Main procedure for data reduction """
        
        log.info("###############################")
        log.info("#### Start data reduction #####")
        log.info("###############################")
        
        dark_flat = False
        
        # Clean old files 
        self.cleanUpFiles()
        
        # Copy/link source files (file or directory) to reduce to the working directory
        papi.linkSourceFiles(self.list_file, self.out_dir)
        files1=[line.replace( "\n", "") for line in fileinput.input(self.list_file)]
        self.m_LAST_FILES=[self.out_dir+"/"+os.path.basename(file_i) for file_i in files1]
        
        ######################################################
        # 0 - Some checks (filter, ....) 
        ######################################################
        try:
            self.checkData()
        except:
            raise 
        ######################################################
        # 00 - Sort out data by MJD (self.m_LAST_FILES)
        ######################################################
        try:
            self.m_LAST_FILES=self.sortOutData()
        except:
            raise
        
        ######################################################
        # 000 - Find out dither mode
        ######################################################
        try:
            self.obs_mode=self.getObsMode()  # overwrite initial given observing mode
        except:
            raise
        
        print "!!!!!!!!!!!!!!!!!!!!!!!!!!"
        print "OBS_MODE=", self.obs_mode
        print "!!!!!!!!!!!!!!!!!!!!!!!!!!"
        
        # keep a copy of the original file names
        self.m_rawFiles = self.m_LAST_FILES        
        
        ######################################
        # 1 - Apply dark, flat
        ######################################
        if self.master_dark!=None and self.master_flat!=None:
            log.debug("---> Applying dark and Flat")
            dark_flat=True
            res = reduce.ApplyDarkFlat(self.m_LAST_FILES, self.master_dark, self.master_flat, self.out_dir)
            self.m_LAST_FILES = res.apply()
            
        ######################################    
        # 2 - Compute Super Sky Flat-Field 
        ######################################
        # - Find out what kind of observing mode we have (dither, ext_dither, ...)
        log.debug('---> Now, compute Super-Sky Flat-Field ...')
        if self.obs_mode=="dither":
            misc.utils.listToFile(self.m_LAST_FILES, self.out_dir+"/files.list")
            superflat = reduce.SuperSkyFlat(self.out_dir+"/files.list", self.out_dir+"/superFlat.fits")
            superflat.create()
        elif self.obs_mode=="dither_on_off" or self.obs_mode=="dither_off_on":
            log.info("EXTENDED SOURCE !!!")
            sky_list=self.getSkyFrames()
            misc.utils.listToFile(sky_list, self.out_dir+"/files.list")
            superflat = reduce.SuperSkyFlat(self.out_dir+"/files.list", self.out_dir+"/superFlat.fits")
            superflat.create()                            
        else:
            log.error("Dither mode not supported")
            raise ("Error, dither mode not supported") 
              
              
        ######################################    
        # 3 - Compute Gain map and apply BPM
        ######################################
        log.debug("Computing gain-map .....")
        nxblock=16
        nyblock=16
        #The next values are to find out bad pixels 
        nsig=5
        mingain=0.7
        maxgain=1.3
        superflat_file = self.out_dir+"/superFlat.fits"
        gainfile = self.out_dir+'/gain_'+self.m_filter+'.fits'
        gain_cmd=self.m_papi_path+'/irdr/bin/gainmap '+ superflat_file + '  ' + gainfile +' '+ str(nsig) + '  ' + str(nxblock) + '  '+ str(nyblock) + '  ' + str(mingain) + '  ' + str(maxgain) 
        
        if misc.utils.runCmd( gain_cmd )==0:
            log.error("Some error while creating gainmap ....")
            return
            
        ########################################
        # Add external Bad Pixel Map to gainmap
        ########################################     
        if self.bpm !=None:
            if not os.path.exists( self.bpm ):
                print('No external Bad Pixel Mask found. Cannot find file : "%s"' %sef.bpm)
            else:
                iraf.imarith(operand1=gainfile,
                  operand2=self.bpm,
                  op='*',
                  result=gainfile.replace(".fits","_bpm.fits"),
                  verbose='yes'
                  )
                #replace old gain file
                os.rename(gainfile.replace(".fits","_bpm.fits"), gainfile)
                        
        #########################################
        # 4 - First Sky subtraction (IRDR) 
        #########################################
        misc.utils.listToFile(self.m_LAST_FILES, self.out_dir+"/skylist1.list")
        self.m_LAST_FILES=self.skyFilter(self.out_dir+"/skylist1.list", gainfile, 'nomask', self.obs_mode)      
                        
        #########################################
        # 4b -Divide by a normalised flat field image 
        #########################################
        # solo para probar a ver como sale (la condicion del IF no es importante)
        #if self.master_dark==None and self.master_flat!=None :
        #    res = reduce.ApplyDarkFlat(self.m_LAST_FILES, None, self.master_flat, self.out_dir)
        #    self.m_LAST_FILES = res.apply()                    
                              
        #########################################
        # 5 - Compute dither offsets from the first sky subtracted/filtered images using cross-correlation
        #########################################
        misc.utils.listToFile(self.m_LAST_FILES, self.out_dir+"/files_skysub.list")
        offset_mat=self.getPointingOffsets (self.out_dir+"/files_skysub.list", 5, 2, self.out_dir+'/offsets1.pap')                
                        
        
        #########################################
        # 6 - First pass coaddition using offsets
        #########################################
        fo=open(self.out_dir+'/offsets1.pap',"r")
        fs=open(self.out_dir+'/stack1.pap','w+')
        for line in fo:
          n_line = line.replace(".skysub.fits.objs", ".skysub.fits") 
          fs.write( n_line )
        fo.close()
        fs.close()    
        self.coaddStackImages(self.out_dir+'/stack1.pap', None, self.out_dir+'/coadd1.fits')
    
        ## END OF SINGLE REDUCTION  ##
        
        #########################################
        # 7 - Create master object mask
        #########################################
        obj_mask=self.createMasterObjMask(self.out_dir+'/coadd1.fits', self.out_dir+'/masterObjMask.fits')
    
        #########################################
        # 8 - Second Sky subtraction (IRDR) 
        #########################################
        #Compound masked sky list
        fs=open(self.out_dir+"/skylist2.pap","w+")
        i=0
        j=0
        for file in self.m_rawFiles:
            if dark_flat: 
                line = file.replace(".fits","_D_F.fits") + " " + obj_mask + " " + str(offset_mat[j][0]) + " " + str(offset_mat[j][1])
            else:
                line = file + " " + obj_mask + " " + str(offset_mat[j][0]) + " " + str(offset_mat[j][1])
            fs.write(line+"\n")
            if (self.obs_mode=='dither_on_off' or self.obs_mode=='dither_off_on') and i%2:
                j=j+1
            elif self.obs_mode=='dither':
                j=j+1
            i=i+1
                
        fs.close()
        self.m_LAST_FILES=self.skyFilter( self.out_dir+"/skylist2.pap", gainfile, 'mask', self.obs_mode)      
    
        #########################################
        # 9 - Create second coadded image of the dithered stack using new sky subtracted frames
        #########################################
        self.coaddStackImages(self.out_dir+'/stack1.pap', None, self.out_dir+'/coadd2.fits')
        reduce.imtrim.imgTrim(self.out_dir+'/coadd2.fits')
        
        #########################################
        # 10 - Make Astrometry
        #########################################
        self.makeAstrometry(self.out_dir+'/coadd2.fits', '2mass', 'noregrid') 
        
        
        
        #########################################
        # 11 - SWARP ???
        #########################################
        
        #Rename to named output file
        shutil.move(self.out_dir+'/coadd2.fits', self.out_file)
        log.info("Generated output file ==>%s", self.out_file)
        
        log.info("##################################")
        log.info("##### End of data reduction ######")
        log.info("##################################")
    
    def cleanUpFiles(self):
        """Clean up files from the working directory, probably from the last execution"""
        
        misc.fileUtils.removefiles(self.out_dir+"/*.fits", self.out_dir+"/c_*", self.out_dir+"/dc_*", self.out_dir+"/*.nip" )
        misc.fileUtils.removefiles(self.out_dir+"/coadd*", self.out_dir+"/*.objs", self.out_dir+"/uparm*", self.out_dir+"/*.skysub*" )
        misc.fileUtils.removefiles(self.out_dir+"/*.head", self.out_dir+"/*.list", self.out_dir+"/*.xml", self.out_dir+"/*.ldac", self.out_dir+"/*.png" )

                    
################################################################################
# main
################################################################################
if __name__ == "__main__":
    log.debug( 'Start PAPI....')
    # Get and check command-line options
    
    
    usage = "usage: %prog [options] arg1 arg2 ..."
    parser = OptionParser(usage)
    
                  
    parser.add_option("-s", "--source",
                  action="store", dest="source_file_list",
                  help="Source file list of data frames. It can be a file or directory name.")
    
    parser.add_option("-o", "--output",
                  action="store", dest="output_filename", help="final reduced output image")
    
    parser.add_option("-t", "--type",
                  action="store", dest="type", default="quick", help="type of reduction (quick|science)")
                  
    parser.add_option("-m", "--obs_mode",
                  action="store", dest="obs_mode", default="dither", help="observing mode (dither|ext_dither)")
    
    parser.add_option("-d", "--outdir",
                  action="store", dest="out_dir", default="/tmp",
                  help="output dir for intermidiate files")
    
    parser.add_option("-D", "--dark",
                  action="store", dest="dark",
                  help="master dark to subtract]")
    
    parser.add_option("-F", "--flat",
                  action="store", dest="flat",
                  help="master flat to divide by]")
    
    parser.add_option("-v", "--verbose",
                  action="store_true", dest="verbose", default=True,
                  help="verbose mode [default]")
    
    parser.add_option("-C", "--config_file",
                  action="store_true", dest="config_file", help="config file for the data reduction process")
                  
    parser.add_option("-S", "--show",
                  action="store_true", dest="show", default=False, help="show final reduced image")
                  
                  
    (options, args) = parser.parse_args()
    
    if not options.source_file_list or not options.output_filename or len(args)!=0: # args is the leftover positional arguments after all options have been processed
        parser.print_help()
        parser.error("incorrect number of arguments " )
    if options.verbose:
        print "reading %s ..." % options.source_file_list
    
    
    redSet = ReductionSet(options.source_file_list, options.out_dir, options.output_filename, options.obs_mode, options.dark, options.flat)
    redSet.reduce()
    
    if options.show==True:
        redSet.show()    
