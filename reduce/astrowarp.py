#! /usr/bin/env python
#encoding:UTF-8

# Copyright (c) 2010 Jose M. Ibanez All rights reserved.
# Institute of Astrophysics of Andalusia, IAA-CSIC
#
# This file is part of PAPI
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


import sys
import os
import shutil
import tempfile
import subprocess 
from optparse import OptionParser
import fileinput
import pyfits


# PAPI modules
import astromatic 
import datahandler
# Logging
from misc.paLog import log

# for initWCS (fk5prec)
#from PyWCSTools import wcscon

def initWCS( input_image ):
    """
    Call this routine to write rough WCS into FITS header and update RA,DEC
    coordinates to J2000.0 equinox; and this way allow SCAMP make astrometry
    with external Catalogs (in J2000??).
    
    Advert: The original file header (input_image) will be modified
    """
    
    try:
        f=datahandler.ClFits ( input_image )
    except Exception,e:
        raise e
    
    fits_file=pyfits.open(input_image, 'update')

    if f.isMEF(): # is a MEF
        raise Exception("Sorry, currently this function only works with simple FITS files with no extensions")
    else:  # is a simple FITS
        header=fits_file[0].header
        try:
            checkWCS(header)
        except Exception,e:
            log.debug("No WCS compliant header, trying to creating one ...")
            try:
                # Read some basic values
                naxis1=f.getNaxis1()
                naxis2=f.getNaxis2()
                ra=f.getRA()
                dec=f.getDec()
                equinox0=f.getEquinox()
                # 
                # Transform RA,Dec to J2000 -->fk5prec(epoch0, 2000.0, &ra, &dec);
                # EQUINOX precessing is DONE by SCAMP !!!
                WCS_J2000=1  #J2000(FK5) right ascension and declination
                WCS_B1950=2  #B1950(FK4) right ascension and declination
                #[new_ra, new_dec]=wcscon.wcscon(WCS_J2000, WCS_J2000, equinox0, 2000.0, ra, dec, 0)
                # Find out PIXSCALE
                if header.has_key("PIXSCALE"):
                    scale=header['PIXSCALE']
                    degscale=scale/3600.0
                else:
                    log.error("Cannot find out the PIXSCALE for the image")
                    fits_file.close()
                    raise Exception("Cannot find out PIXSCALE for the image")
                    
                #Create initial WCS
                #
                header.update("CRPIX1", naxis1/2.0, "Ref. pixel in <axis direction>")
                header.update("CRPIX2", naxis2/2.0, "Ref. pixel in <axis direction>")
                header.update("CRVAL1", ra, "Coordinate value of ref. pixel")
                header.update("CRVAL2", dec, "Coordinate value of ref. pixel")
                #header.update("RA", new_ra, "Coordinate value of ref. pixel")
                #header.update("DEC", new_dec, "Coordinate value of ref. pixel")
                header.update("CTYPE1", "RA---TAN", "Pixel coordinate system")
                header.update("CTYPE2", "DEC--TAN", "Pixel coordinate system")
                #header.update("RADECSYS","FK5","Coordinate reference frame")
                # CD matrix (the CDi_j elements) encode the sky position angle,
                # the pixel scale, and a possible flipping. 
                header.update("CD1_1",-degscale, "Translation matrix element")
                header.update("CD1_2",0.0, "Translation matrix element")
                header.update("CD2_1",0.0, "Translation matrix element")
                header.update("CD2_2",degscale, "Translation matrix element")
                header.update("SCALE",scale, "Image scale")
                #header.update("EQUINOX", 2000.0, "Standard FK5(years)")
                
                # clean imcompatible CDi_j and CDELT matrices
                if header.has_key("CDELT1"):
                    del header["CDELT1"]
                if header.has_key("CDELT2"):
                    del header["CDELT2"]
                
                log.debug("Successful WCS header created !")
                
            except Exception,e:
                log.error("Some error while creating initial WCS header...", str(e))
                fits_file.close(output_verify='ignore') # This ignore any FITS standar violation and allow write/update the FITS file
                raise e
            
        fits_file.close()
            
        log.debug("Right WCS info")
            
def checkWCS( header ):
    """
    Checks for a variety of WCS keywords and raise an Exception if the header lacks a proper
    combination of them.  This is needed because wcstools will not raise any
    sort of error if a WCS isn't present or is malformed, and SCAMP (E.Bertin)
    need a initial WCS information. This is probably 90% complete in terms of
    its checking for the types of FITS files that we are likely to be using.
    
    If you find any WCS keywords that cause wcstools to behave in an erratic
    manner without signaling errors, add them to this method.  Experience has
    shown that the astrophysical community has an uncanny ability to produce
    data sets that cause FITS readers and WCS projections to break.  It is
    important that we check for irregular cases and flag them before the code
    runs and produces confusing results.  Our only defense against this is
    experience with unusual data sets, so the more checks here the better.
    
    TODO(): Implement stricter checking under CDELT case, in
    particular for full PC matrices (both kinds), as well as LATPOLE and
    LONPOLE.  It would also be good to check for illegal values, but that's
    a lot of work.
    """
    
    keywords_to_check=['NAXIS1','NAXIS2','CTYPE1','CTYPE2','CRVAL1','CRVAL2',
                       'CRPIX1','CRPIX2']
    
    # Every header must have these keywords.
    for kw in keywords_to_check:
        if not header.has_key(kw):
            log.debug("Keyword %s not found",kw)
            raise Exception("Keyword %s not found",kw)
    
    # Check for the equinox, which can be specified in more than 1 way.
    if not header.has_key('EPOCH') and not header.has_key('EQUINOX'):
        log.debug("Missing keyword EPOCH or EQUINOX")
        raise Exception("Missing keyword EPOCH or EQUINOX")
        
    # Check some values
    if header['CTYPE1']=='PIXEL' or header['CTYPE2']=='PIXEL':
        log.debug("Wrong CTYPE value (PIXEL) for WCS header")
        raise Exception ("Wrong CTYPE value (PIXEL) for WCS header")
        
    # Check for CDi_j or CDELT matrix
    # CDELT matrix : Here we should probably be more rigorous and check
    # for a full PC matrix or CROTA value, but for now
    # this is pretty good.
    if not header.has_key('CD1_1') or not header.has_key('CD1_2') \
        or not header.has_key('CD2_1') or not header.has_key('CD2_2'):
            if not header.has_key('CDELT1') or not header.has_key("CDELT2"):
                log.debug("Couldn't find a complete set of CDi_j matrix or CDELT")
                raise Exception("Couldn't find a complete set of CDi_j matrix or CDELT")
                

    
def doAstrometry( input_image, output_image=None, catalog='2MASS', config_dict=None):
    """ Do the astrometric calibration to the input image (only one)
      
    The method does astrometry on the image entirely using Emmanuel Bertin's
    tools, namely SExtractor, SCAMP and SWarp. First, SExtractor is run on the
    image, saving the output catalog in the FITS_LDAC binary format that is then
    read by SCAMP, whose output is in turn a '.head' file. This file contains
    the updated astrometric information, which is finally merged with the 
    original image using the SWarp stacking tool. In case no output path is
    specified, the resulting, calibrated image is saved to a temporary
    directory, most likely /tmp.

    
    Note that the method does not receive any other parameters. How SExtractor,
    SCAMP and SWarp work is determined by their respective configuration files,
    located in the config file specified or with default values, and they are what 
    should be modified in case you need to adjust their operation so that they better suit
    your needs.

    Keyword arguments:
  
    input_image  - path to the FITS image whose astrometry is to be done.
    output_image - path to which the astrometically calibrated version of the
                  image. If no path is provided, it will be saved to a temporary
                  location.
    catalog      - external catalog to be used for the astrometric calibration. 
                   By default, 2MASS is used 

    """
    
    log.debug("*** Start Astrometric calibration ***")
    
    if not config_dict:
            raise Exception("Config dictionary not provided ...")
            
    # Save the resulting image to a temporary file if no path was specified
    if output_image is None:
        output_fd, output_path = tempfile.mkstemp(suffix='.fits')
        os.close(output_fd)

    ## STEP 0: Run IRDR::initwcs to initialize rough WCS header, thus modify the file headers
    #initwcs also converts to J2000.0 EQUINOX
    log.debug("***Doing WCS-header initialization ...")
    initWCS(input_image)
    """
    initwcs_path=config_dict['config_files']['irdr_bin']+"/initwcs" #os.environ['PAPI_HOME']+'/irdr/bin/initwcs'
    args = [initwcs_path, input_image]
    print "ARGS=", args
    ret_code = subprocess.call(args)
    if ret_code!=0:
        raise RuntimeError("There was an error while running 'initwcs'")
    """
    
    ## STEP 1: Create SExtractor catalog (.ldac)
    log.debug("*** Creating SExtractor catalog ....")
    sex = astromatic.SExtractor()
    #sex.config['CONFIG_FILE']="/disk-a/caha/panic/DEVELOP/PIPELINE/PANIC/trunk/config_files/sex.conf"
    #sex.ext_config['CHECKIMAGE_TYPE'] = "OBJECTS"
    sex.config['CATALOG_TYPE'] = "FITS_LDAC"
    sex.config['CATALOG_NAME'] = input_image + ".ldac"
    sex.config['DETECT_THRESH'] = config_dict['skysub']['mask_thresh']
    sex.config['DETECT_MINAREA'] = config_dict['skysub']['mask_minarea']
    try:
        sex.run(input_image, updateconfig=True, clean=False)
    except: 
        raise
                
    ## STEP 2: Make astrometric calibration 
    log.debug("Doing astrometric calibration....")
    scamp = astromatic.SCAMP()
    scamp.config['CONFIG_FILE']=config_dict['config_files']['scamp_conf']#"/disk-a/caha/panic/DEVELOP/PIPELINE/PANIC/trunk/config_files/scamp.conf"
    scamp.ext_config['ASTREF_CATALOG']=catalog
    scamp.ext_config['SOLVE_PHOTOM']="N"
    scamp.ext_config['CHECKPLOT_TYPE']="NONE"
    scamp.ext_config['WRITE_XML']="N"
    cat_file = input_image + ".ldac"   # xxxxx.fits.ldac
    #updateconfig=False means scamp will use the specified config file instead of the single config parameters
    #but, "ext_config" parameters will be used in any case
    try:
        scamp.run(cat_file, updateconfig=False, clean=False)
    except:
        raise
    
    ## STEP 3: Merge the astrometric parameters (.head keywords) with SWARP, and using .head files created by SCAMP
    log.debug("Merging astrometric calibration parameters and resampling ...")
    swarp = astromatic.SWARP()
    swarp.config['CONFIG_FILE']=config_dict['config_files']['swarp_conf']#"/disk-a/caha/panic/DEVELOP/PIPELINE/PANIC/trunk/config_files/swarp.conf"
    swarp.ext_config['IMAGEOUT_NAME']=output_image
    swarp.ext_config['COPY_KEYWORDS']='OBJECT,INSTRUME,TELESCOPE,IMAGETYP,FILTER,FILTER2,SCALE,MJD-OBS,RA,DEC'
    swarp.ext_config['WEIGHTOUT_NAME']=output_image.replace(".fits",".weight.fits")
    swarp.ext_config['HEADER_SUFFIX']=".head"
    #Rename the external header produced by SCAMP (.head) to a filename to be looked for by SWARP 
    #SWARP take into account for re-projection an externar header for file 'xxxxx.ext' if exists 
    #an 'xxxx.head' header ('ext' can be any string, i.e. fits, skysub, ....)
    shutil.move(input_image+".head", os.path.splitext(input_image)[0]+".head") # we use explitext() to remove the LAST suffix after, thus 'xxx.fits.skysub' ---> has as extension '.skysub'
    if (os.path.exists(input_image.replace(".fits",".weight.fits"))):
        swarp.ext_config['WEIGHT_TYPE']='MAP_WEIGHT'
        swarp.ext_config['WEIGHT_SUFFIX']='.weight.fits'
        swarp.ext_config['WEIGHT_IMAGE']=input_image.replace(".fits",".weight.fits")
        
    try:
        swarp.run(input_image, updateconfig=False, clean=False)
    except:
        raise
    
    return output_image 

class AstroWarp(object):
    """ Astrometric warping """

    def __init__(self, input_files, catalog='2MASS', coadded_file="/tmp/astrowarp.fits", config_dict=None):
        """ Instantiation method for AstroWarp class

        Keyword arguments:
        input_files  - the set of overlapping reduced frames, i.e., dark subtracted, flatted and sky subtracted 
        catalog      - the catalog to use for the astrometric calibration (by default, 2MASS)
        coadded_file - the output final coadded image
        
        """
        # TODO: I have to provide an alternate way to get a default config dictionary ...
        if not config_dict:
            raise Exception("Config dictionary not provided ...")
            
        self.input_files = input_files
        self.catalog = catalog
        self.coadded_file = coadded_file
        self.config_dict = config_dict # the config dictionary
        
    def run(self):
        """ Start the computing of the coadded image, following the next steps:
      
            0. Initialize rought WCS header
            1. Call up SExtractor to create pixel object list (ldac catalog)
            2. Call up SCAMP to make the astrometric calibration of the overlapped set of frames (.head calibration)
            3. Call up SWARP to make the coadd, using the .head files generated by SCAMP and containing the distortion model parameters
            4. Make the final astrometric calibration the the coadded frame (SExtractor+SCAMP)
            
        """
        
        log.debug("*** Start Astrowarp ***")

        ## STEP 0: Run IRDR::initwcs to initialize rough WCS header, thus modify the file headers
        # initwcs also converts to J2000.0 EQUINOX
        # TBD: re-implement in Python method the call to 'irdr:initwcs'
        log.debug("***Doing WCS-header initialization ...")
        #initwcs_path=self.config_dict['config_files']['irdr_bin']+"/initwcs"
        for file in self.input_files:
            initWCS(file)
            """
            args = [initwcs_path, file]
            #print "ARGS=", args
            ret_code = subprocess.call(args)
            if ret_code!=0:
                raise RuntimeError("There was an error while running 'initwcs'")
            """   
        ## STEP 1: Create SExtractor catalogs (.ldac)
        log.debug("*** Creating SExtractor catalog ....")
        for file in self.input_files:
            sex = astromatic.SExtractor()
            #sex.config['CONFIG_FILE']="/disk-a/caha/panic/DEVELOP/PIPELINE/PANIC/trunk/config_files/sex.conf"
            #sex.ext_config['CHECKIMAGE_TYPE'] = "OBJECTS"
            sex.config['CATALOG_TYPE'] = "FITS_LDAC"
            sex.config['CATALOG_NAME'] = file + ".ldac"
            sex.config['DETECT_THRESH'] = self.config_dict['skysub']['mask_thresh']
            sex.config['DETECT_MINAREA'] = self.config_dict['skysub']['mask_minarea']
            sex.run(file, updateconfig=True, clean=False)
                        
        ## STEP 2: Make the multi-astrometric calibration for each file (all overlapped-files together)
        log.debug("*** Doing multi-astrometric calibration ....")
        scamp = astromatic.SCAMP()
        scamp.config['CONFIG_FILE']=self.config_dict['config_files']['scamp_conf']#"/disk-a/caha/panic/DEVELOP/PIPELINE/PANIC/trunk/config_files/scamp.conf"
        scamp.ext_config['ASTREF_CATALOG']=self.catalog
        scamp.ext_config['SOLVE_PHOTOM']="N"
        cat_files = [(f + ".ldac") for f in self.input_files]
        #updateconfig=False means scamp will use the specified config file instead of the single config parameters
        scamp.run(cat_files, updateconfig=False, clean=False)
        
        
        ## STEP 3: Make the coadding with SWARP, and using .head files created by SCAMP
        # It requires the files are overlapped, i.e., have an common sky-area
        log.debug("*** Coadding overlapped files....")
        swarp = astromatic.SWARP()
        swarp.config['CONFIG_FILE']=self.config_dict['config_files']['swarp_conf']#"/disk-a/caha/panic/DEVELOP/PIPELINE/PANIC/trunk/config_files/swarp.conf"
        swarp.ext_config['COPY_KEYWORDS']='OBJECT,INSTRUME,TELESCOPE,IMAGETYP,FILTER,FILTER2,SCALE,MJD-OBS'
        swarp.ext_config['IMAGEOUT_NAME']=os.path.dirname(self.coadded_file)+ "/coadd_tmp.fits"
        swarp.ext_config['WEIGHTOUT_NAME']=os.path.dirname(self.coadded_file)+ "/coadd_tmp.weight.fits"
        if os.path.isfile(self.input_files[0].replace(".fits",".weight.fits")):
            swarp.ext_config['WEIGHT_TYPE']='MAP_WEIGHT'
            swarp.ext_config['WEIGHT_SUFFIX']='.weight.fits'
        swarp.run(self.input_files, updateconfig=False, clean=False)
        
        ## STEP 4: Make again the final astrometric calibration to the final coadd
        log.debug("*** Doing final astrometril calibration....")
        doAstrometry(os.path.dirname(self.coadded_file) + "/coadd_tmp.fits", self.coadded_file, self.catalog, self.config_dict)
        
        
        
################################################################################
# main
################################################################################
if __name__ == "__main__":
    log.debug( 'Start AstroWarp....')
    # Get and check command-line options
    
    
    usage = "usage: %prog [options] arg1 arg2 ..."
    parser = OptionParser(usage)
    
                  
    parser.add_option("-s", "--source",
                  action="store", dest="source_file_list",
                  help="Source file list of data frames. It can be a file or directory name.")
    
    parser.add_option("-o", "--output",
                  action="store", dest="output_filename", help="final coadded output image")
    
    parser.add_option("-v", "--verbose",
                  action="store_true", dest="verbose", default=True,
                  help="verbose mode [default]")
    
                                
    (options, args) = parser.parse_args()
    
    if not options.source_file_list or not options.output_filename or len(args)!=0: # args is the leftover positional arguments after all options have been processed
        parser.print_help()
        parser.error("incorrect number of arguments " )
    if options.verbose:
        print "reading %s ..." % options.source_file_list
    
    filelist=[line.replace( "\n", "") for line in fileinput.input(options.source_file_list)]
    astrowarp = AstroWarp(filelist, catalog="2MASS", coadded_file=options.output_filename)
    
    try:
        astrowarp.run()
    except:
        log.error("Some error while running Astrowarp....")
        raise
        sys.exit()
    
    