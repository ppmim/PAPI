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
import astropy.io.fits as fits


# PAPI modules
import astromatic
import astromatic.ldac
import datahandler

# Logging
from misc.paLog import log
from misc.version import __version__
import misc.config
import reduce.solveAstrometry

# for initWCS (fk5prec)
#from PyWCSTools import wcscon

def initWCS( input_image, pixel_scale):
    """
    Call this routine to write rough WCS into FITS header and update RA,DEC
    coordinates to J2000.0 equinox; and this way allow SCAMP make astrometry
    with external Catalogs (in J2000??).
    
    Warning: The original file header (input_image) will be modified
    """
    
    log.debug("Entering in initWCS!")
    
    try:
        f = datahandler.ClFits(input_image, check_integrity=False)
    except Exception,e:
        log.error("Error reading FITS %s : %s"%(f,str(e)))
        raise e
    
    try:
        fits_file = fits.open(input_image, 'update', ignore_missing_end=True)
    except Exception,e:
        log.error("Error reading FITS %s : %s"%(f,str(e)))
        raise e

    if f.isMEF(): # is a MEF
        #raise Exception("Sorry, currently this function only works with simple FITS files with no extensions")
        log.warning("MEF file detected !")
    
    for ext in range(0, len(fits_file)):
        if len(fits_file)>1 and ext==0:
            continue
    #else:  # is a simple FITS
        header = fits_file[ext].header
        try:
            checkWCS(header)
            log.debug("FITS looks having a right WCS header")
        except Exception,e:
            log.warning("No WCS compliant header, trying to create one ...")
            try:
                # Read some basic values
                naxis1 = header['NAXIS1']
                naxis2 = header['NAXIS2']
                ra = f.ra
                dec = f.dec 
                # 
                # Transform RA,Dec to J2000 -->fk5prec(epoch0, 2000.0, &ra, &dec);
                # EQUINOX precessing is DONE by SCAMP !!!
                #equinox0 = f.getEquinox()
                WCS_J2000 = 1  #J2000(FK5) right ascension and declination
                WCS_B1950 = 2  #B1950(FK4) right ascension and declination
                #[new_ra, new_dec]=wcscon.wcscon(WCS_J2000, WCS_J2000, equinox0, 2000.0, ra, dec, 0)
                # Find out PIXSCALE
                if "PIXSCALE" in header:
                    # scale must be in arcsec/pixel
                    scale = header['PIXSCALE']
                    # now convert to deg/pixel
                    degscale = scale/3600.0
                else:
                    scale = pixel_scale
                    degscale = scale/3600.0
                    log.warning("Cannot find out the PIXSCALE for the image.")
                    log.warning("Using Pixel scale = %s"%pixel_scale)
                    #fits_file.close()
                    #raise Exception("Cannot find out PIXSCALE for the image")

                create_wcs = True                
                if create_wcs:
                    #Create initial WCS
                    #
                    header.set("CRPIX1", naxis1/2.0, "RA and DEC reference pixel along axis 1")
                    header.set("CRPIX2", naxis2/2.0, "RA and DEC reference pixel along axis 2")
                    header.set("CRVAL1", ra, "[deg] RA Coordinate value of ref. pixel")
                    header.set("CRVAL2", dec, "[deg] DEC Coordinate value of ref. pixel")
                    #header.set("RA", new_ra, "Coordinate value of ref. pixel")
                    #header.set("DEC", new_dec, "Coordinate value of ref. pixel")
                    header.set("CUNIT1", "deg", "WCS units along axis 1")
                    header.set("CUNIT2", "deg", "WCS units along axis 2")
                    header.set("CTYPE1", "RA---TAN", "Pixel coordinate system")
                    header.set("CTYPE2", "DEC--TAN", "Pixel coordinate system")
                    #header.set("RADECSYS","FK5","Coordinate reference frame")
                    # CD matrix (the CDi_j elements) encode the sky position angle,
                    # the pixel scale, and a possible flipping.
                    # CD1_1 is <0 because East is supposed at Left = flipX
                    # CD2_2 is >0 because North is supposed at Up
                    # In addition, it must be noted that:
                    # CD1_1 = cos(r), CD1_2 = sin(r), CD2_1 = -sin(r), CD2_2 = cos(r)
                    # r = clockwise rotation_angle  
                    header.set("CD1_1", -degscale, "[deg/px] Translation WCS matrix element")
                    header.set("CD1_2", 0.0, "[deg/px] Translation WCS matrix element")
                    header.set("CD2_1", 0.0, "[deg/px] Translation WCS matrix element")
                    header.set("CD2_2", degscale, "[deg/px] Translation WCS matrix element")
                    header.set("SCALE", scale, "[arcsec/px] Image scale")
                    #header.set("EQUINOX", 2000.0, "Standard FK5(years)")
                else:
                    header.set("RA", ra, "Right Ascension (degree)")
                    header.set("DEC", dec, "Declination (degree)")
                    
                # clean incompatible CDi_j and CDELT matrices
                if "CDELT1" in header:
                    del header["CDELT1"]
                if "CDELT2" in header:
                    del header["CDELT2"]
                
                log.debug("Successful WCS header created !")
                
            except Exception,e:
                log.error("Some error while creating initial WCS header: %s"%str(e))
                fits_file.close(output_verify='ignore') # This ignore any FITS standar violation and allow write/update the FITS file
                raise e
        
        # Check whether CRPIXn need to be updated because of some kind of border
        # was added around the image during a previus coadd.
        try:
            if 'CRPIX1' in header and header['NAXIS1']>2048:
                log.info("Updating CRPIX1 taking into account border pixels due to coadd")
                log.debug("NAXIS1=%s  CRPIX1=%s"%(header['NAXIS1'],header['CRPIX1']))
                value = header['CRPIX1'] + (header['NAXIS1']-2048)/2
                header.set('CRPIX1', value )
                log.debug("VALUE1=%s"%value)
            if 'CRPIX2' in header and header['NAXIS2']>2048:
                log.info("Updating CRPIX2 taking into account border pixels due to coadd")
                value = header['CRPIX2'] + (header['NAXIS2']-2048)/2            
                header.set('CRPIX2', value )
                log.debug("VALUE2=%s"%value)
        except Exception,e:
            log.critial("[initWCS] Error updating header: %s"%str(e))
            raise e

    fits_file[0].header.set('PAPIVERS', __version__, 'PANIC Pipeline version') 
    
    try:        
        fits_file.close(output_verify='ignore')
    except Exception,e:
        log.critical("ERROR !")
        raise e

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

    ## pruebas con CCD Roper OSN
    ##raise Exception("Forzamos la creacion de un nuevo header")

    # Every header must have these keywords.
    for kw in keywords_to_check:
        if kw not in header:
            log.debug("Keyword %s not found",kw)
            raise Exception("Keyword %s not found"%kw)
    
    #
    # Check for the equinox, which can be specified in more than 1 way.
    # I am not sure if required by SCAMP ?
    #
    #if 'EPOCH' not in header and 'EQUINOX' not in header:
    #    log.debug("Missing keyword EPOCH or EQUINOX")
    #    raise Exception("Missing keyword EPOCH or EQUINOX")
        
    # Check some values
    if header['CTYPE1']=='PIXEL' or header['CTYPE2']=='PIXEL':
        log.debug("Wrong CTYPE value (PIXEL-Cartesian Coordinates) for WCS header")
        raise Exception ("Wrong CTYPE value (PIXEL-Cartesian Coordinates) for WCS header")
        
    # Check for CDi_j or CDELT matrix
    # CDELT matrix : Here we should probably be more rigorous and check
    # for a full PC matrix or CROTA value, but for now
    # this is pretty good.
    if 'CD1_1' not in header or 'CD1_2' not in header \
        or 'CD2_1' not in header or 'CD2_2' not in header:
            if 'CDELT1' not in header or 'CDELT2' not in header:
                log.debug("Couldn't find a complete set of CDi_j matrix or CDELT")
                raise Exception("Couldn't find a complete set of CDi_j matrix or CDELT")
             

    
def doAstrometry(input_image, output_image=None, catalog='2MASS', 
                  config_dict=None, do_votable=False,
                  resample=True, subtract_back=True):
    
    """ Do the astrometric calibration to an input image using SCAMP
      
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
    should be modified in case you need to adjust their operation so that they 
    better suit your needs.
    
    STEP 0) Run initwcs to initialize rough WCS header
    STEP 1) Create SExtractor catalog (.ldac)  (SExtractor)
    STEP 2) Make astrometric calibration       (SCAMP)
    STEP 3) Merge and Warp the astrometric parameters (.head keywords) with 
            SWARP, and using .head files created by SCAMP.
            This step imply the distortion correction of the final coadd.
    STEP 4) (optionaly) Create SExtractor catalog (ascii.votable)
    

    Note
    ----
    In priciple, this function supports MEF files because Astromatic software
    suite also supports it, but is has not been deeply tested.


    Parameters
    ----------
  
    input_image  - path to the FITS image whose astrometry is to be done.
    output_image - path to which the astrometically calibrated version of the
                  image. If no path is provided, it will be saved to a temporary
                  location.
    catalog      - external catalog to be used for the astrometric calibration. 
                   By default, 2MASS is used 

    config_dict  - dictionary with config arguments
    
    do_votable   - If True, build-out a VO-Table with the Sextractor catalog 
                including the fields required for the photometry calibration. 
    
    resample     - Resample image when SWARP is executed (RESAMPLE)
    
    subtract_back - Subtract sky background when scamp is executed (SUBTRACT_BACK)
    
    
    Returns
    -------
    
    Filename of astrometric calibrated image obtained.
     
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
    
    pixel_scale = config_dict['general']['pix_scale']
    
    try:
        initWCS(input_image, pixel_scale)
    except Exception,e:
        raise e
        
    
    ## STEP 1: Create SExtractor catalog (.ldac)
    log.debug("*** Creating SExtractor catalog ....")
    sex = astromatic.SExtractor()
    #sex.config['CONFIG_FILE']="/disk-a/caha/panic/DEVELOP/PIPELINE/PANIC/trunk/config_files/sex.conf"
    #sex.ext_config['CHECKIMAGE_TYPE'] = "OBJECTS"
    sex.config['CATALOG_TYPE'] = "FITS_LDAC"
    sex.config['CATALOG_NAME'] = input_image + ".ldac"
    sex.config['DETECT_THRESH'] = config_dict['astrometry']['mask_thresh']
    sex.config['DETECT_MINAREA'] = config_dict['astrometry']['mask_minarea']
    
    # SATUR_LEVEL and NCOADD
    try:
        dh = datahandler.ClFits(input_image, check_integrity=False)
        nc = dh.getNcoadds()
    except Exception,e:
        log.warning("Cannot read NCOADDS. Default value (=1) taken")
        nc = 1

    sex.config['SATUR_LEVEL'] = int(nc) * int(config_dict['astrometry']['satur_level'])
    
    try:
        sex.run(input_image, updateconfig=True, clean=False)
    except Exception,e:
        log.error("Error running SExtractor: %s"%str(e)) 
        raise e
    
    # A test to apply single point mask
    #filter_area(input_image + ".ldac",100)
    ## end of test
    
    # PAPI_HOME
    try:
        papi_home = os.environ['PAPI_HOME']
        if papi_home[-1]!='/':
            papi_home+='/'
    except Exception,e:
        log.error("Error, variable PAPI_HOME not defined.")
        raise e

    ## STEP 2: Make astrometric calibration 
    log.debug("Doing astrometric calibration....")
    scamp = astromatic.SCAMP()
    scamp.config['CONFIG_FILE'] = papi_home + config_dict['config_files']['scamp_conf']
    #"/disk-a/caha/panic/DEVELOP/PIPELINE/PANIC/trunk/config_files/scamp.conf"
    scamp.ext_config['ASTREF_CATALOG'] = catalog
    scamp.ext_config['SOLVE_PHOTOM'] = "N"
    # next parameters (POSANGLE_MAXERR, POSITION_MAXERR) are very important in 
    # order to be able to solve fields with large errors
    scamp.ext_config['MATCH_FLIPPED'] = 'Y'
    scamp.ext_config['POSANGLE_MAXERR'] = 5
    scamp.ext_config['POSITION_MAXERR'] = 5

    #scamp.ext_config['CHECKPLOT_TYPE'] = "NONE"
    scamp.ext_config['WRITE_XML'] = "N"
    cat_file = input_image + ".ldac"   # xxxxx.fits.ldac
    #updateconfig=False means scamp will use the specified config file instead 
    #of the single config parameters
    #but, "ext_config" parameters will be used in any case
    try:
        scamp.run(cat_file, updateconfig=False, clean=False)
    except Exception,e:
        log.error("Error running SCAMP: %s"%str(e))
        raise e
    
    ## STEP 3: Merge and Warp the astrometric parameters (.head keywords) with 
    ## SWARP, and using .head files created by SCAMP. Therefore, a field distortion
    ## correction is done if resample is YES.
    log.debug("Merging astrometric calibration parameters and re-sampling ...")
    swarp = astromatic.SWARP()
    swarp.config['CONFIG_FILE'] = papi_home + config_dict['config_files']['swarp_conf']
    swarp.ext_config['IMAGEOUT_NAME'] = output_image
    swarp.ext_config['COPY_KEYWORDS'] = 'OBJECT,INSTRUME,TELESCOPE,IMAGETYP,FILTER,FILTER1,FILTER2,SCALE,MJD-OBS,RA,DEC,HISTORY,NCOADDS,NDIT'
    basename_o, extension_o = os.path.splitext(output_image)
    #"Projected" weight-maps are created too, even if no weight-maps were given in input.
    swarp.ext_config['WEIGHTOUT_NAME'] = basename_o + ".weight" + extension_o
    basename, extension = os.path.splitext(input_image)
    swarp.ext_config['HEADER_SUFFIX'] = extension + ".head"
    if not os.path.isfile(input_image + ".head"):
        raise Exception ("Cannot find required .head file")
   
    if not resample:
        swarp.ext_config['RESAMPLE'] = 'N' # then, no field distortion removing is done

    if not subtract_back:
        swarp.ext_config['SUBTRACT_BACK'] = 'N'
   
     
    # Rename the external header produced by SCAMP (.head) to a filename to be 
    # looked for by SWARP.
    # SWARP take into account for re-projection an external header for 
    # file 'xxxxx.ext' if exists 
    # an 'xxxx.head' header ('ext' can be any string, i.e. fits, skysub, ....)
    # We use splitext() to remove the LAST suffix after, thus 
    # 'xxx.fits.skysub' ---> has as extension '.skysub'
    # mmmmm, it depends on HEADER_SUFFIX config variable !!!! so, above comments 
    # are not completely true.
    # But ...I'd rather to modify HEADER_SUFFIX than rename the file, efficiency !
    
    ##basename, extension = os.path.splitext(input_image)
    ##shutil.move(input_image+".head", basename +".head")  # very important !!
   
    # To avoid any problem concerning SWARP because by chance could exist
    # a old file about IMAGEOUT_NAME.head what would cause a bad resampling
    # and combination of files, we remove any IMAGEOUT_NAME.head
    if os.path.exists(swarp.ext_config['IMAGEOUT_NAME']+".head"):
        os.remove(swarp.ext_config['IMAGEOUT_NAME']+".head")
    
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
        # SATUR_LEVEL and NCOADD
        try:
            dh = datahandler.ClFits(input_image, check_integrity=False)
            nc = dh.getNcoadds()
        except:
            log.warning("Cannot read NCOADDS. Taken default (=1)")
            nc = 1
    
        sex.config['SATUR_LEVEL'] = int(nc) * int(config_dict['astrometry']['satur_level'])

        
        try:
            sex.run(output_image, updateconfig=True, clean=False)
        except: 
            raise
        
    return output_image 

def filter_area(cat_filename, max_size=200):
    """
    Filter the input catalog dropping out the objects with an area
    bigger than max_size.
    
    Parameters
    ----------

    cat_name: str
        filename of the LDAC catalog (SExtractor generated)
    
    max_size: int 
        maximun area size of objects; objects with a higher value are removed.
    
    LDAC = Leiden Data Analysis Center
    
    Notes
    -----
    The input filename will be overwritten.
    LDAC = Leiden Data Analysis Center
    
    Returns
    -------
    Filename of the new catalog, otherwise None or exception
    
    Todo
    ----
    Test if works fine, not yet tested ! 
    
    """
    
    log.debug("Filtering catalog by ISOAREA_IMAGE")

    import numpy
    try:
        hdus = fits.open(cat_filename, "update")
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
    """ 
    Astrometric warping, and in principle, taking into account field distortion
    during the coadding/warping.
    MEF files are not supported (at least with Astrometry.net)
    """

    def __init__(self, input_files, catalog=None, 
                 coadded_file="/tmp/astrowarp.fits", config_dict=None, 
                 do_votable=False, resample=True, subtract_back=True
                 ):
        """ 
        Instantiation method for AstroWarp class.

        Parameters
        ----------
        input_files  - the set of overlapping reduced frames, i.e., dark 
                       subtracted, flatted and sky subtracted. 
        catalog      - the catalog to use for the astrometric calibration (by default, 2MASS).
        coadded_file - the output final coadded image.
        resample     - Resample image when scamp is executed (RESAMPLE).
        subtract_back - Subtract sky background when scamp is executed (SUBTRACT_BACK)
        pixel_scale   - default pixel scale of the image used for initWCS().

        """
        

        # PAPI_HOME
        try:
            self.papi_home = os.environ['PAPI_HOME']
            if self.papi_home[-1]!='/':
                self.papi_home+='/'
        except Exception,e:
            log.error("Error, variable PAPI_HOME not defined.")
            raise Exception("Error, variable PAPI_HOME not defined")

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
        self.resample = resample
        self.subtract_back = subtract_back
        self.pix_scale = config_dict['general']['pix_scale']
        self.temp_dir = config_dict['general']['temp_dir']


    def run(self, engine='SCAMP'):
        """
        Run the astrometric calibration and coadding, taking into account the
        distortion.

        Parameters
        ----------
        engine    - default astrometric calibration engine to use. The possible
                    values are ('SCAMP', 'Astrometry.net').
        """

        if engine=='SCAMP':
            self.runWithSCAMP()
        else:
            self.runWithAstrometryNet()

    def runWithAstrometryNet(self):
        """ 
        Start the computing of the coadded image, following the next steps:
      
            0. Call solveAstrometry for first astrometric calibration
            1. Call SExtractor to create pixel object list (ldac catalog)
            2. Call SCAMP to make the distortion calibration of the 
            overlapped set of frames (.head calibration)
            3. Call SWARP to make the coadd and regriding, using the .head files 
            generated by SCAMP and containing the distortion model parameters
            4. Make the final astrometric calibration the the coadded frame (Astrometry.Net).
            
        
        Returns
        -------
        The name of the coadded output file generated.

        """
        
        log.info("*** Start Astrowarp (Astrometry.net engine) ***")

        ## STEP 0: Run solveAstrometry for first 
        log.debug("***Running solveAstrometry initialization ...")
        # Cannot call solveAstrometry.runMultiSolver because of next error:
        # "daemonic processes are not allowed to have children".
        # A pool cannot call a function that creates a new pool.
        # http://stackoverflow.com/questions/6974695/python-process-pool-non-daemonic
        # Then, we solve sequentially.
        solved_files = []
        
        # Check if input images have already been astrometically calibrated.
        for file in self.input_files:
            solved_msg = "--Start of Astrometry.net WCS solution--"
            if not solved_msg in fits.getheader(file)['COMMENT']:
                try:
                    solved = reduce.solveAstrometry.solveField( 
                                            file,
                                            self.temp_dir,
                                            self.config_dict['general']['pix_scale'])
                except Exception,e:
                    raise Exception("[runWithAstrometryNet] Cannot solve Astrometry for file: %s"%(file,str(e)))
                else:
                    solved_files.append(solved)
            else:
                log.warning("Image %s already astrometrically solved by Astrometry.net"%file)
                solved_files.append(file)

        ## STEP 1: Create SExtractor catalogs (.ldac)
        log.debug("*** Creating objects catalog (SExtractor)....")
        for file in solved_files:
            sex = astromatic.SExtractor()
            #sex.config['CONFIG_FILE']="/disk-a/caha/panic/DEVELOP/PIPELINE/PANIC/trunk/config_files/sex.conf"
            #sex.ext_config['CHECKIMAGE_TYPE'] = "OBJECTS"
            sex.config['CATALOG_TYPE'] = "FITS_LDAC"
            sex.config['CATALOG_NAME'] = file + ".ldac"
            sex.config['DETECT_THRESH'] = self.config_dict['astrometry']['mask_thresh']
            sex.config['DETECT_MINAREA'] = self.config_dict['astrometry']['mask_minarea']
            # SATUR_LEVEL and NCOADD
            try:
                dh = datahandler.ClFits(file, check_integrity=False)
                nc = dh.getNcoadds()
            except:
                log.warning("Cannot read NCOADDS. Taken default (=1)")
                nc = 1

            sex.config['SATUR_LEVEL'] = int(nc) * int(self.config_dict['astrometry']['satur_level'])
            
            try:
                sex.run(file, updateconfig=True, clean=False)
            except Exception,e:
                raise e

            
            #filter_area(file + ".ldac")    
                        
        ## STEP 2: Make the multi-astrometric calibration for each file (all overlapped-files together)
        log.debug("*** Doing multi-astrometric calibration (SCAMP)....")
        scamp = astromatic.SCAMP()
        scamp.config['CONFIG_FILE'] = self.papi_home + self.config_dict['config_files']['scamp_conf']
        scamp.ext_config['ASTREF_CATALOG'] = self.catalog
        scamp.ext_config['SOLVE_PHOTOM'] = "N"
        cat_files = [(f + ".ldac") for f in solved_files]
        #updateconfig=False means scamp will use the specified config file instead of the single config parameters
        
        try:
            scamp.run(cat_files, updateconfig=False, clean=False)
        except Exception,e:
            raise e
        
        
        ## STEP 3: Make the coadding with SWARP, and using .head files created by SCAMP
        # It requires the files are overlapped, i.e., have an common sky-area
        log.debug("*** Coadding overlapped files (SWARP)....")
        
        
        root, ext = os.path.splitext(os.path.basename(self.coadded_file))
        # Path to the temporary FITS file containing the WCS header
        kwargs = dict(prefix = '%s_coadd_' % root, suffix = ext, dir=self.temp_dir)
        with tempfile.NamedTemporaryFile(**kwargs) as fd:
            output_path = fd.name
        
        swarp = astromatic.SWARP()
        swarp.config['CONFIG_FILE'] = self.papi_home + self.config_dict['config_files']['swarp_conf']
        basename, extension = os.path.splitext(solved_files[0])
        swarp.ext_config['HEADER_SUFFIX'] = extension + ".head"  # very important !
        if not os.path.isfile(solved_files[0]+".head"):
            raise Exception ("Cannot find required .head file")
        
        
        swarp.ext_config['COPY_KEYWORDS'] = 'OBJECT,INSTRUME,TELESCOPE,IMAGETYP,FILTER,FILTER1,FILTER2,SCALE,MJD-OBS,RA,DEC,HISTORY,NCOADDS,NDIT'
        swarp.ext_config['IMAGEOUT_NAME'] = output_path
        swarp.ext_config['WEIGHTOUT_NAME'] = output_path.replace(".fits", ".weight.fits")
                                                                 
        # "Projected" weight-maps are created only if weight-maps were given in input.
        # That is not true !!!
        if os.path.isfile(basename + ".weight" + extension):
            swarp.ext_config['WEIGHT_TYPE'] = 'MAP_WEIGHT'
            swarp.ext_config['WEIGHT_SUFFIX'] = '.weight' + extension
            swarp.ext_config['WEIGHTOUT_NAME'] = output_path.replace(".fits", ".weight.fits")
        
        if not self.resample:
            swarp.ext_config['RESAMPLE'] = 'N' # then, no field distortion removing is done
  
        if not self.subtract_back:
            swarp.ext_config['SUBTRACT_BACK'] = 'N'
            
        # To avoid any problem concerning SWARP because by chance could exist
        # a old file about IMAGEOUT_NAME.head what would cause a bad resampling
        # and combination of files, we remove any IMAGEOUT_NAME.head
        if os.path.exists(swarp.ext_config['IMAGEOUT_NAME']+".head"):
            os.remove(swarp.ext_config['IMAGEOUT_NAME']+".head")
                
        """
        #IMPORTANT: Rename the .head files in order to be found by SWARP
        #but it's much efficient modify the HEADER_SUFFIX value (see above)
        for aFile in self.input_files:
            basename, extension = os.path.splitext(aFile)
            shutil.move(aFile+".head", basename +".head")  # very important !!
        """    
        try:
            swarp.run(solved_files, updateconfig=False, clean=False)
        except Exception,e:    
            raise e
        
        ## STEP 4: Make again the final astrometric calibration (only 
        ## if we coadded more that one file) to the final coadd.
        ## TODO: I am not sure if it is needed to do again ?????
        if (len(self.input_files)>1):
            log.debug("*** Doing final astrometric calibration....")
            try:
                solved = reduce.solveAstrometry.solveField(
                            output_path, 
                            self.temp_dir,
                            self.config_dict['general']['pix_scale'])
            except Exception,e:
                    raise Exception("[runWithAstrometryNet] Error doing Astrometric calibration: %s"%str(e))
            else:
                shutil.move(solved, self.coadded_file)
                shutil.move(output_path.replace(".fits", ".weight.fits"), 
                        self.coadded_file.replace(".fits", ".weight.fits"))
        else:
            shutil.move(output_path, self.coadded_file)
            shutil.move(output_path.replace(".fits", ".weight.fits"), 
                        self.coadded_file.replace(".fits", ".weight.fits"))
        
        log.info("Lucky you ! file %s created", self.coadded_file)

        return self.coadded_file


    def runWithSCAMP(self):
        """ Start the computing of the coadded image, following the next steps:
      
            0. Initialize rought WCS header
            1. Call SExtractor to create pixel object list (ldac catalog)
            2. Call SCAMP to make the astrometric calibration of the overlapped set 
            of frames (.head calibration)
            3. Call SWARP to make the coadd and regriding, using the .head files generated by 
               SCAMP and containing the distortion model parameters
            4. Make the final astrometric calibration the the coadded frame (SExtractor+SCAMP)
            
            TODO: final astrometric calibration does not work fine ...(2011-09-21)
        
            Returns
            -------
            The name of the coadded output file generated.
        
        """
        
        log.info("*** Start Astrowarp ***")

        ## STEP 0: Run IRDR::initwcs to initialize rough WCS header, thus modify the file headers
        # initwcs also converts to J2000.0 EQUINOX
        log.debug("***Doing WCS-header initialization ...")
        for file in self.input_files:
            log.debug("file: %s",file)
            initWCS(file, self.pix_scale)

        ## STEP 1: Create SExtractor catalogs (.ldac)
        log.debug("*** Creating objects catalog (SExtractor)....")
        for file in self.input_files:
            sex = astromatic.SExtractor()
            sex.config['CATALOG_TYPE'] = "FITS_LDAC"
            sex.config['CATALOG_NAME'] = file + ".ldac"
            sex.config['DETECT_THRESH'] = self.config_dict['astrometry']['mask_thresh']
            sex.config['DETECT_MINAREA'] = self.config_dict['astrometry']['mask_minarea']
            # SATUR_LEVEL and NCOADD
            try:
                dh = datahandler.ClFits(file, check_integrity=False)
                nc = dh.getNcoadds()
            except:
                log.warning("Cannot read NCOADDS. Taken default (=1)")
                nc = 1
        
            sex.config['SATUR_LEVEL'] = int(nc) * int(self.config_dict['astrometry']['satur_level'])

            
            try:
                sex.run(file, updateconfig=True, clean=False)
            except Exception,e:
                raise e

            
            #filter_area(file + ".ldac")    
                        
        ## STEP 2: Make the multi-astrometric calibration for each file (all overlapped-files together)
        log.debug("*** Doing multi-astrometric calibration (SCAMP)....")
        scamp = astromatic.SCAMP()
        scamp.config['CONFIG_FILE'] = self.papi_home + self.config_dict['config_files']['scamp_conf']
        scamp.ext_config['ASTREF_CATALOG'] = self.catalog
        scamp.ext_config['SOLVE_PHOTOM'] = "N"
        cat_files = [(f + ".ldac") for f in self.input_files]
        #updateconfig=False means scamp will use the specified config file instead of the single config parameters
        
        try:
            scamp.run(cat_files, updateconfig=False, clean=False)
        except Exception,e:
            raise e
        
        
        ## STEP 3: Make the coadding with SWARP, and using .head files created by SCAMP
        # It requires the files are overlapped, i.e., have an common sky-area
        log.debug("*** Coadding overlapped files (SWARP)....")
        swarp = astromatic.SWARP()
        swarp.config['CONFIG_FILE'] = self.papi_home + self.config_dict['config_files']['swarp_conf']
        basename, extension = os.path.splitext(self.input_files[0])
        swarp.ext_config['HEADER_SUFFIX'] = extension + ".head"  # very important !
        if not os.path.isfile(self.input_files[0]+".head"):
            raise Exception ("Cannot find required .head file")
            
        swarp.ext_config['COPY_KEYWORDS'] = 'OBJECT,INSTRUME,TELESCOPE,IMAGETYP,FILTER,FILTER1,FILTER2,SCALE,MJD-OBS,RA,DEC,HISTORY,NCOADDS,NDIT'
        swarp.ext_config['IMAGEOUT_NAME'] = os.path.dirname(self.coadded_file) + "/coadd_tmp.fits"
        #"Projected" weight-maps are created only if weight-maps were given in input.
        if os.path.isfile(basename + ".weight" + extension):
            swarp.ext_config['WEIGHT_TYPE'] = 'MAP_WEIGHT'
            swarp.ext_config['WEIGHT_SUFFIX'] = '.weight' + extension
            swarp.ext_config['WEIGHTOUT_NAME'] = os.path.dirname(self.coadded_file) + "/coadd_tmp.weight.fits"
        
        if not self.resample:
            swarp.ext_config['RESAMPLE'] = 'N' # then, no field distortion removing is done
  
        if not self.subtract_back:
            swarp.ext_config['SUBTRACT_BACK'] = 'N'

        #To avoid any problem concerning SWARP because by chance could exist
        #a old file about IMAGEOUT_NAME.head what would cause a bad resampling
        #and combination of files, we remove any IMAGEOUT_NAME.head
        if os.path.exists(swarp.ext_config['IMAGEOUT_NAME']+".head"):
            os.remove(swarp.ext_config['IMAGEOUT_NAME']+".head")
                
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
                         self.config_dict, self.do_votable,
                         self.resample, self.subtract_back)
        else:
            shutil.move(os.path.dirname(self.coadded_file) + "/coadd_tmp.fits", 
                        self.coadded_file)
        
        log.info("Lucky you ! file %s created", self.coadded_file)
        return self.coadded_file


################################################################################
# main
################################################################################
if __name__ == "__main__":
    
    log.debug( 'Start AstroWarp....')
    
    # Get and check command-line options
    usage = "usage: %prog [options]"
    desc = """Performs the alignment and warping of a set of images,
in principle previously reduced, but not mandatory.
"""
    parser = OptionParser(usage, description=desc)
    
    parser.add_option("-c", "--config_file",
                  action="store", dest="config_file", help="Mandatory PAPI configuration file.")
                  
    parser.add_option("-s", "--source",
                  action="store", dest="source_file",
                  help="Source input file. It can be a FITS file or"
                  "text file with a list of FITS files.")
    
    parser.add_option("-o", "--output",
                  action="store", dest="output_filename", 
                  help="final coadded output image")
    
    parser.add_option("-r", "--resample",
                  action="store_true", dest="resample", default=False,
                  help="Resample input image [default=%default]")
    
    parser.add_option("-b", "--subtract_back",
                  action="store_true", dest="subtract_back", default=False,
                  help="Subtract sky background [default=%default]")
    
    parser.add_option("-e", "--engine",
                  action="store", dest="engine", default="SCAMP",
                  help="Astrometric engine to use (SCAMP, AstrometryNet) [default=%default]")
    
    
                                
    (options, args) = parser.parse_args()
    
    if len(sys.argv[1:])<1:
       parser.print_help()
       sys.exit(0)

    # Read the default configuration file
    # If none was specified by the user, environment variable will be used
    if not options.config_file:
        try:
            config_file = os.environ['PAPI_CONFIG']
        except KeyError, error:
            print 'Environment variable PAPI_CONFIG not found!'
            sys.exit()
    else:
        config_file = options.config_file
        
    cfg_options = misc.config.read_config_file(config_file)
    
    # args is the leftover positional arguments after all options have been processed
    if not options.source_file or not options.output_filename or len(args)!=0: 
        parser.print_help()
        parser.error("incorrect number of arguments " )
    
    # Check if source_file is a FITS file or a text file listing a set of files
    if os.path.exists(options.source_file):
        try:
            datahandler.fits_simple_verify(options.source_file)
            filelist = [options.source_file]
        except:
            filelist = [line.replace( "\n", "") 
                      for line in fileinput.input(options.source_file)]
            
        if len(filelist)==1:
            # Astromatic
            if options.engine=="SCAMP":
                try:
                    log.debug("[Astrowarp] Solving with SCAMP engine")
                    doAstrometry(filelist[0], output_image=options.output_filename,
                              catalog='USNO-B1', config_dict=cfg_options, 
                              do_votable=False, resample=options.resample,
                              subtract_back=options.subtract_back)
                except Exception,e:
                    log.error("Some error while doing astrometric calibration")
            # AstrometryNet
            else: 
                try:
                    log.debug("[Astrowarp] Solving with Astrometry.Net engine")
                    solved = reduce.solveAstrometry.solveField( 
                                        filelist[0],
                                        cfg_options['general']['temp_dir'],
                                        cfg_options['general']['pix_scale'])
                except Exception,e:
                    raise Exception("Cannot solve Astrometry for file: %s"%(file,str(e)))
                else:
                    shutil.move(solved, options.output_filename)     
        else: 
            astrowarp = AstroWarp(filelist, catalog="2MASS", 
                                  coadded_file=options.output_filename, 
                                  config_dict=cfg_options,
                                  resample=options.resample,
                                  subtract_back=options.subtract_back)
        
            try:
                astrowarp.run(options.engine)
            except Exception,e:
                log.error("Some error while running Astrowarp : %s",str(e))
    else:
        log.error("Source file %s does not exists",options.source_file)
        
    sys.exit()
        
    
