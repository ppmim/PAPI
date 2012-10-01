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
from optparse import OptionParser
import fileinput
import pyfits


# PAPI modules
import astromatic
import astromatic.ldac
import datahandler

# Logging
from misc.paLog import log
import misc.config

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
        f = datahandler.ClFits ( input_image )
    except Exception,e:
        raise e
    
    fits_file=pyfits.open(input_image, 'update', ignore_missing_end=True)

    if f.isMEF(): # is a MEF
        raise Exception("Sorry, currently this function only works with simple FITS files with no extensions")
    else:  # is a simple FITS
        header = fits_file[0].header
        try:
            checkWCS(header)
            log.debug("FITS looks having a right WCS header")
        except Exception,e:
            log.debug("No WCS compliant header, trying to creating one ...")
            try:
                # Read some basic values
                naxis1 = f.getNaxis1()
                naxis2 = f.getNaxis2()
                ra = f.ra
                dec = f.dec
                equinox0=f.getEquinox()
                # 
                # Transform RA,Dec to J2000 -->fk5prec(epoch0, 2000.0, &ra, &dec);
                # EQUINOX precessing is DONE by SCAMP !!!
                WCS_J2000 = 1  #J2000(FK5) right ascension and declination
                WCS_B1950 = 2  #B1950(FK4) right ascension and declination
                #[new_ra, new_dec]=wcscon.wcscon(WCS_J2000, WCS_J2000, equinox0, 2000.0, ra, dec, 0)
                # Find out PIXSCALE
                if "PIXSCALE" in header:
                    scale = header['PIXSCALE']
                    degscale = scale/3600.0
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
                if "CDELT1" in header:
                    del header["CDELT1"]
                if "CDELT2" in header:
                    del header["CDELT2"]
                
                log.debug("Successful WCS header created !")
                
            except Exception,e:
                log.error("Some error while creating initial WCS header...", str(e))
                fits_file.close(output_verify='ignore') # This ignore any FITS standar violation and allow write/update the FITS file
                raise e
            
        fits_file.close(output_verify='ignore')
            
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
        if kw not in header:
            log.debug("Keyword %s not found",kw)
            raise Exception("Keyword %s not found"%kw)
    
    # Check for the equinox, which can be specified in more than 1 way.
    if 'EPOCH' not in header and 'EQUINOX' not in header:
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
    if 'CD1_1' not in header or 'CD1_2' not in header \
        or 'CD2_1' not in header or 'CD2_2' not in header:
            if 'CDELT1' not in header or 'CDELT2' not in header:
                log.debug("Couldn't find a complete set of CDi_j matrix or CDELT")
                raise Exception("Couldn't find a complete set of CDi_j matrix or CDELT")
                

    
def doAstrometry( input_image, output_image=None, catalog='2MASS', 
                  config_dict=None, do_votable=False):
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

    config_dict  - dictionary with config arguments
    
    do_votable   - If True, build-out a VO-Table (SEx) with the fields required 
                   for the photometry calibration. 
    """
    
    log.debug("[doAstrometry] *** Start Astrometric calibration ***")
    
    if not config_dict:
            raise Exception("Config dictionary not provided ...")
            
    # Save the resulting image to a temporary file if no path was specified
    if output_image is None:
        output_fd, output_path = tempfile.mkstemp(suffix='.fits')
        os.close(output_fd)

    ## STEP 0: Run IRDR::initwcs to initialize rough WCS header, thus modify the file headers
    #initwcs also converts to J2000.0 EQUINOX
    log.debug("[doAstrometry] ***Doing WCS-header initialization ...")
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
    sex.config['DETECT_THRESH'] = config_dict['astrometry']['mask_thresh']
    sex.config['DETECT_MINAREA'] = config_dict['astrometry']['mask_minarea']
    try:
        sex.run(input_image, updateconfig=True, clean=False)
    except: 
        raise
    
    # A test to apply single point mask
    #filter_area(input_image + ".ldac",100)
    ## end of test
    
    ## STEP 2: Make astrometric calibration 
    log.debug("Doing astrometric calibration....")
    scamp = astromatic.SCAMP()
    scamp.config['CONFIG_FILE'] = config_dict['config_files']['scamp_conf']
    #"/disk-a/caha/panic/DEVELOP/PIPELINE/PANIC/trunk/config_files/scamp.conf"
    scamp.ext_config['ASTREF_CATALOG'] = catalog
    scamp.ext_config['SOLVE_PHOTOM'] = "N"
    #scamp.ext_config['CHECKPLOT_TYPE'] = "NONE"
    scamp.ext_config['WRITE_XML'] = "N"
    cat_file = input_image + ".ldac"   # xxxxx.fits.ldac
    #updateconfig=False means scamp will use the specified config file instead of the single config parameters
    #but, "ext_config" parameters will be used in any case
    try:
        scamp.run(cat_file, updateconfig=False, clean=False)
    except:
        raise
    
    ## STEP 3: Merge and Warp the astrometric parameters (.head keywords) with SWARP, and using .head files created by SCAMP
    log.debug("Merging astrometric calibration parameters and resampling ...")
    swarp = astromatic.SWARP()
    swarp.config['CONFIG_FILE'] = config_dict['config_files']['swarp_conf']
    #"/disk-a/caha/panic/DEVELOP/PIPELINE/PANIC/trunk/config_files/swarp.conf"
    swarp.ext_config['IMAGEOUT_NAME'] = output_image
    swarp.ext_config['COPY_KEYWORDS'] = 'OBJECT,INSTRUME,TELESCOPE,IMAGETYP,FILTER,FILTER1,FILTER2,SCALE,MJD-OBS,RA,DEC'
    basename_o, extension_o = os.path.splitext(output_image)
    swarp.ext_config['WEIGHTOUT_NAME'] = basename_o + ".weight" + extension_o
    basename, extension = os.path.splitext(input_image)
    swarp.ext_config['HEADER_SUFFIX'] = extension + ".head"
    if not os.path.isfile(input_image + ".head"):
        raise Exception ("Cannot find required .head file")
    
    #Rename the external header produced by SCAMP (.head) to a filename to be looked for by SWARP 
    #SWARP take into account for re-projection an external header for file 'xxxxx.ext' if exists 
    #an 'xxxx.head' header ('ext' can be any string, i.e. fits, skysub, ....)
    #We use splitext() to remove the LAST suffix after, thus 'xxx.fits.skysub' ---> has as extension '.skysub'
    # mmmmm, it depends on HEADER_SUFFIX config variable !!!! so, above comments are not completely true
    # But ...I'd rather to modify HEADER_SUFFIX than rename the file, efficiency !
    ##basename, extension = os.path.splitext(input_image)
    ##shutil.move(input_image+".head", basename +".head")  # very important !!
    
    if os.path.isfile(basename + ".weight" + extension):
        swarp.ext_config['WEIGHT_TYPE'] = 'MAP_WEIGHT'
        swarp.ext_config['WEIGHT_SUFFIX'] = '.weight' + extension
        swarp.ext_config['WEIGHT_IMAGE'] = basename + ".weight" + extension
        
    try:
        swarp.run(input_image, updateconfig=False, clean=False)
    except:
        raise

    ## STEP 4: Create SExtractor catalog (ascii.votable)
    if do_votable:
        log.debug("*** Creating SExtractor VOTable catalog ....")
        sex = astromatic.SExtractor()
        #sex.config['CONFIG_FILE']="/disk-a/caha/panic/DEVELOP/PIPELINE/PANIC/trunk/config_files/sex.conf"
        #sex.ext_config['CHECKIMAGE_TYPE'] = "OBJECTS"
        sex.config['CATALOG_TYPE'] = "ASCII_VOTABLE"
        sex.config['CATALOG_NAME'] = os.path.splitext(output_image)[0] + ".xml"
        sex.config['DETECT_THRESH'] = config_dict['skysub']['mask_thresh']
        sex.config['DETECT_MINAREA'] = config_dict['skysub']['mask_minarea']
        try:
            sex.run(output_image, updateconfig=True, clean=False)
        except: 
            raise
        
    return output_image 

def filter_area(cat_filename, max_size=200):
    """
    @summary: Filter the input catalog dropping out the objects with an area
    bigger than max_size

    @param cat_name: filename of the LDAC catalog
    @param max_size: maximun area size of objects; objects with a higher value
    are removed.
    
    @requires: the input catalog must be a LDAC catalog (SExtractor FITS_LDAC)
    LDAC = Leiden Data Analysis Center
    
    @note: the input filename will be overwritten
    
    @return: filename of the new catalog, otherwise None or exception
    
    @todo: test if works fine   ! not yet tested 
    
    """
    
    log.debug("Filtering catalog by ISOAREA_IMAGE")

    import numpy
    try:
        hdus = pyfits.open(cat_filename, "update")
        #mask = hdus[2].data.field('ISOAREA_IMAGE')<max_size
        mask = numpy.logical_and( hdus[2].data.field('ISOAREA_IMAGE')<max_size,
                                  hdus[2].data.field('FLAGS')==0)
        
        hdus[2].data = hdus[2].data[mask]
        hdus.writeto(cat_filename,output_verify='warn', clobber=True)
    except Exception,e:
        raise Exception("Error while filtering Catalog %s  : %s"%(cat_filename,str(e)))

    log.debug("Catalog filtered by Area (isoarea_image). File  : %s",cat_filename)
    
    return cat_filename

    """
    try:
        cat = astromatic.ldac.openObjectFile(cat_filename, table='LDAC_OBJECTS')
        mask = cat['ISOAREA_IMAGE']<max_size
        new_cat = cat.filter(mask)
        new_cat.saveas(cat_filename, clobber=True)
    except Exception,e:
        raise Exception("Error while filtering Catalog %s  : %s"%(cat_filename,str(e)))

    log.debug("Catalog filtered by Area (isoarea_image). File  : %s",cat_filename)
    
    return cat_filename
    """

    
class AstroWarp(object):
    """ Astrometric warping """

    def __init__(self, input_files, catalog=None, 
                 coadded_file="/tmp/astrowarp.fits", config_dict=None, do_votable=False):
        """ Instantiation method for AstroWarp class

        Keyword arguments:
        input_files  - the set of overlapping reduced frames, i.e., dark subtracted, flatted and sky subtracted 
        catalog      - the catalog to use for the astrometric calibration (by default, 2MASS)
        coadded_file - the output final coadded image
        
        """
        # TODO: I have to provide an alternate way to get a default config dictionary ...
        if not config_dict:
            raise Exception("Config dictionary not provided ...")
        else:
            self.config_dict = config_dict # the config dictionary
            
        self.input_files = input_files
        if catalog!=None:
            self.catalog = catalog
        else: 
            self.catalog = config_dict['astrometry']['catalog']
        
        self.coadded_file = coadded_file
        self.do_votable = do_votable
        
    def run(self):
        """ Start the computing of the coadded image, following the next steps:
      
            0. Initialize rought WCS header
            1. Call SExtractor to create pixel object list (ldac catalog)
            2. Call SCAMP to make the astrometric calibration of the overlapped set 
            of frames (.head calibration)
            3. Call SWARP to make the coadd, using the .head files generated by 
               SCAMP and containing the distortion model parameters
            4. Make the final astrometric calibration the the coadded frame (SExtractor+SCAMP)
            
            @todo: final astrometric calibration does not work fine ...(2011-09-21)
        """
        
        log.info("*** Start Astrowarp ***")

        ## STEP 0: Run IRDR::initwcs to initialize rough WCS header, thus modify the file headers
        # initwcs also converts to J2000.0 EQUINOX
        # TBD: re-implement in Python method the call to 'irdr:initwcs'
        log.debug("***Doing WCS-header initialization ...")
        #initwcs_path=self.config_dict['config_files']['irdr_bin']+"/initwcs"
        for file in self.input_files:
            log.debug("file: %s",file)
            initWCS(file)
            """
            args = [initwcs_path, file]
            #print "ARGS=", args
            ret_code = subprocess.call(args)
            if ret_code!=0:
                raise RuntimeError("There was an error while running 'initwcs'")
            """   
        ## STEP 1: Create SExtractor catalogs (.ldac)
        log.debug("*** Creating objects catalog (SExtractor)....")
        for file in self.input_files:
            sex = astromatic.SExtractor()
            #sex.config['CONFIG_FILE']="/disk-a/caha/panic/DEVELOP/PIPELINE/PANIC/trunk/config_files/sex.conf"
            #sex.ext_config['CHECKIMAGE_TYPE'] = "OBJECTS"
            sex.config['CATALOG_TYPE'] = "FITS_LDAC"
            sex.config['CATALOG_NAME'] = file + ".ldac"
            sex.config['DETECT_THRESH'] = self.config_dict['astrometry']['mask_thresh']
            sex.config['DETECT_MINAREA'] = self.config_dict['astrometry']['mask_minarea']
            try:
                sex.run(file, updateconfig=True, clean=False)
            except Exception,e:
                raise e

            
            #filter_area(file + ".ldac")    
                        
        ## STEP 2: Make the multi-astrometric calibration for each file (all overlapped-files together)
        log.debug("*** Doing multi-astrometric calibration (SCAMP)....")
        scamp = astromatic.SCAMP()
        scamp.config['CONFIG_FILE'] = self.config_dict['config_files']['scamp_conf']
        #"/disk-a/caha/panic/DEVELOP/PIPELINE/PANIC/trunk/config_files/scamp.conf"
        scamp.ext_config['ASTREF_CATALOG'] = self.catalog
        scamp.ext_config['SOLVE_PHOTOM'] = "N"
        cat_files = [(f + ".ldac") for f in self.input_files]
        #updateconfig=False means scamp will use the specified config file instead of the single config parameters
        
        print "CATFILES = ",cat_files
        try:
            scamp.run(cat_files, updateconfig=False, clean=False)
        except Exception,e:
            raise e
        
        
        ## STEP 3: Make the coadding with SWARP, and using .head files created by SCAMP
        # It requires the files are overlapped, i.e., have an common sky-area
        log.debug("*** Coadding overlapped files (SWARP)....")
        swarp = astromatic.SWARP()
        swarp.config['CONFIG_FILE'] = self.config_dict['config_files']['swarp_conf']
        #"/disk-a/caha/panic/DEVELOP/PIPELINE/PANIC/trunk/config_files/swarp.conf"
        basename, extension = os.path.splitext(self.input_files[0])
        swarp.ext_config['HEADER_SUFFIX'] = extension + ".head"  # very important !
        if not os.path.isfile(self.input_files[0]+".head"):
            raise Exception ("Cannot find required .head file")
            
        swarp.ext_config['COPY_KEYWORDS'] = 'OBJECT,INSTRUME,TELESCOPE,IMAGETYP,FILTER,FILTER1,FILTER2,SCALE,MJD-OBS,HISTORY'
        swarp.ext_config['IMAGEOUT_NAME'] = os.path.dirname(self.coadded_file) + "/coadd_tmp.fits"
        if os.path.isfile(basename + ".weight" + extension):
            swarp.ext_config['WEIGHT_TYPE'] = 'MAP_WEIGHT'
            swarp.ext_config['WEIGHT_SUFFIX'] = '.weight' + extension
            swarp.ext_config['WEIGHTOUT_NAME'] = os.path.dirname(self.coadded_file) + "/coadd_tmp.weight.fits"
        
        """
        #IMPORTANT: Rename the .head files in order to be found by SWARP
        #but it's much efficient modify the HEADER_SUFFIX value (see above)
        for aFile in self.input_files:
            basename, extension = os.path.splitext(aFile)
            shutil.move(aFile+".head", basename +".head")  # very important !!
        """    
        try:
            swarp.run(self.input_files, updateconfig=False, clean=False)
        except Exception,e:    
            raise e
        
        ## STEP 4: Make again the final astrometric calibration (only 
        ## if we coadded more that one file) to the final coadd)
        ## TODO: I am not sure if it is needed to do again ?????
        if (len(self.input_files)>1):
            log.debug("*** Doing final astrometric calibration....")
            doAstrometry(os.path.dirname(self.coadded_file) + "/coadd_tmp.fits", 
                         self.coadded_file, self.catalog, 
                         self.config_dict, self.do_votable)
        else:
            shutil.move(os.path.dirname(self.coadded_file) + "/coadd_tmp.fits", 
                        self.coadded_file)
        
        log.info("Lucky you ! file %s created", self.coadded_file)
        
################################################################################
# main
################################################################################
if __name__ == "__main__":
    
    log.debug( 'Start AstroWarp....')

    # Get and check command-line options
    usage = "usage: %prog [options] arg1 arg2 ..."
    parser = OptionParser(usage)
    
    parser.add_option("-c", "--config_file",
                  action="store", dest="config_file", help="config file")
                  
    parser.add_option("-s", "--source",
                  action="store", dest="source_file",
                  help="Source file list of data frames. It can be a file or directory name.")
    
    
    parser.add_option("-o", "--output",
                  action="store", dest="output_filename", help="final coadded output image")
    
    
    parser.add_option("-v", "--verbose",
                  action="store_true", dest="verbose", default=True,
                  help="verbose mode [default]")
    
                                
    (options, args) = parser.parse_args()
    
    # Read the default configuration file if none was specified by the user
    if not options.config_file:
        config_file = misc.config.default_config_file()       
    else:
        config_file = options.config_file
    
    # now, we "mix" the invokation parameter values with the values in the config file,
    # having priority the invokation values over config file values
    # note: only values of the 'general' section can be invoked
    cfg_options = misc.config.read_config_file(config_file)
    
    
    # args is the leftover positional arguments after all options have been processed
    if not options.source_file or not options.output_filename or len(args)!=0: 
        parser.print_help()
        parser.error("incorrect number of arguments " )
    if options.verbose:
        print "reading %s ..." % options.source_file
    
    # Check if source_file is a FITS file or a text file listing a set of files
    if os.path.exists(options.source_file):
        try:
            datahandler.fits_simple_verify(options.source_file)
            filelist=[options.source_file]
        except:
            filelist=[line.replace( "\n", "") for line in fileinput.input(options.source_file)]
            
        ##astrowarp = AstroWarp(filelist, catalog="2MASS", coadded_file=options.output_filename, config_dict=cfg_options)
        astrowarp = AstroWarp(filelist, catalog="2MASS", coadded_file=options.output_filename, config_dict=cfg_options)
    
        try:
            astrowarp.run()
        except Exception,e:
            log.error("Some error while running Astrowarp : %s",str(e))
    else:
        log.error("Source file %s does not exists",options.source_file)
        
    sys.exit()
        
    
