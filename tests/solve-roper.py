#! /usr/bin/env python
#encoding:UTF-8

# Copyright (c) 2013 Jose M. Ibanez All rights reserved.
# Institute of Astrophysics of Andalusia, IAA-CSIC
#
# This file is part of ROPER CCD astrometry procedure
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
import logging
import subprocess
import glob

import numpy
import pywcs
import pyfits
import sys

def checkHeader(filename):
    """
    Check header from OSN FITS files and return scale and instrument.
    """
    
    hdulist = pyfits.open(filename)

    # Check is a single HDU
    if len(hdulist) >1:
        logging.error("MEF files not supported")
        raise Exception("MEF file not supported")
        
    # Check some basic keywords
    if 'INSTRUME' in hdulist[0].header:
        instrument = hdulist[0].header['INSTRUME']
    else:
        instrument = None
        
    if 'NAXIS1' in hdulist[0].header:
        naxis1 = hdulist[0].header['NAXIS1']
    else:
        naxis1 = -1
        
    if 'NAXIS2' in hdulist[0].header:
        naxis2 = hdulist[0].header['NAXIS2']
    else:
        naxis2 = -1
            
    if 'XBINNING' in hdulist[0].header:
        xbinning = hdulist[0].header['XBINNING']
    else:
        xbinning = -1
            
    if 'YBINNING' in hdulist[0].header:
        ybinning = hdulist[0].header['XBINNING']
    else:
        ybinning = -1

    
    if instrument=='Roper':
        # pixel scale for Roper@T150 (eaest focus) is 0.23 arcsec/pixel
        scale = 0.23 
    elif instrument=='RoperT90':
        # pixel scale for Roper@T90 is (west focus) 0.3857 arcsec/pixel
        scale = 0.3857
    else:
        #default scale (?)
        scale = 0.45
        
    pix_scale_x = scale * xbinning
    pix_scale_y = scale * ybinning
    
    
    #===========================================================================
    # try:
    #    checkWCS(hdulist[0].header)
    # except Exception,e:
    #    raise e
    #===========================================================================
    
    hdulist.close()
    
    return pix_scale_x
        
def loadWCS(filename):
    """
    Load WCS from input image
    """
    
    # Load the FITS hdulist using pyfits
    hdulist = pyfits.open(fileaname)

    # Parse the WCS keywords in the primary HDU
    wcs = pywcs.WCS(hdulist[0].header)

    # Print out the "name" of the WCS, as defined in the FITS header
    print wcs.wcs.name

    # Print out all of the settings that were parsed from the header
    wcs.wcs.print_contents()


def initWCS( input_image, pixel_scale):
    """
    Call this routine to write rough WCS into FITS header and update RA,DEC
    coordinates to J2000.0 equinox; and this way allow SCAMP make astrometry
    with external Catalogs (in J2000??).
    
    Warning: The original file header (input_image) will be modified
    """
    
    try:
        f = datahandler.ClFits ( input_image )
    except Exception,e:
        raise e
    
    fits_file = pyfits.open(input_image, 'update', ignore_missing_end=True)

    if f.isMEF(): # is a MEF
        raise Exception("Sorry, currently this function only works with simple FITS files with no extensions")
    else:  # is a simple FITS
        header = fits_file[0].header
        try:
            checkWCS(header)
            log.debug("FITS looks having a right WCS header")
        except Exception,e:
            log.debug("No WCS compliant header, trying to create one ...")
            try:
                # Read some basic values
                naxis1 = f.getNaxis1()
                naxis2 = f.getNaxis2()
                ra = f.ra
                dec = f.dec
                equinox0 = f.getEquinox()
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
                    # CD1_1 is <0 because East is supposed at Left = flipX
                    # CD2_2 is >0 because North is supposed at Up
                    # In addition, it must be noted that:
                    # CD1_1 = cos(r), CD1_2 = sin(r), CD2_1 = -sin(r), CD2_2 = cos(r)
                    # r = clockwise rotation_angle  
                    header.update("CD1_1", -degscale, "Translation matrix element")
                    header.update("CD1_2", 0.0, "Translation matrix element")
                    header.update("CD2_1", 0.0, "Translation matrix element")
                    header.update("CD2_2", degscale, "Translation matrix element")
                    header.update("SCALE", scale, "Image scale")
                    #header.update("EQUINOX", 2000.0, "Standard FK5(years)")
                else:
                    header.update("RA", ra, "Right Ascension (degree)")
                    header.update("DEC", dec, "Declination (degree)")
                    
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


    raise Exception("Forzamos la creacion de un nuevo header")

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
                

    
def solveField(filename, pix_scale, tmp_dir):
    """
    Do astrometric calibration to the given filename using Astrometry.net 
    function 'solve-field'
    """
    
    ## Create temporal directory
    if not os.path.isdir(tmp_dir):
        os.makedirs(tmp_dir)
    
    scale = checkHeader(filename)
    
    ## Create log file
    logging.debug("Starting to solve-Field to: %s  scale=%s  tmp_dir=%s"%(filename, scale, tmp_dir))
    
    ## Run calibration external command
    """solve-field --scale-units arcsecperpix --scale-low 0.4 --scale-high 0.5 /home/panic/as/150/PG1647+056-003R.fit -O -d 20 -p"""
    
    str_cmd = "solve-field -O -p -d 10,20,30,40 --scale-units arcsecperpix --scale-low %s --scale-high %s %s"%(pix_scale-0.1, pix_scale+0.1, filename)
    
    try:
        p = subprocess.Popen(str_cmd, bufsize=0, shell=True, stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             close_fds=True)
    except Exception,e:
        raise e

    #Warning
    #We use communicate() rather than .stdin.write, .stdout.read or .stderr.read 
    #to avoid deadlocks due to any of the other OS pipe buffers filling up and 
    #blocking the child process.(Python Ref.doc)

    (stdoutdata, stderrdata) = p.communicate()
    err = stdoutdata + " " + stderrdata

    print "STDOUTDATA=",stdoutdata
    print "STDERRDATA=",stderrdata
    
    if len(err)>1:
        print "[solveField]: STDOUT + STDERR = ", err
        logging.error(err)
    else:
        logging.debug("Field solved !")

################################################################################
# main
################################################################################
if __name__ == "__main__":
    
    
    # Get and check command-line options
    usage = "usage: %prog [options] arg1 arg2 ..."
    desc = """Performs the astrometric calibration of a set of images,
in principle previously reduced, but not mandatory.

"""
    parser = OptionParser(usage, description=desc)
    
    parser.add_option("-c", "--config_file",
                  action="store", dest="config_file", help="config file")
                  
    parser.add_option("-s", "--source",
                  action="store", dest="source_file",
                  help="Source file list of data frames. It can be a file or directory name.")
    
    parser.add_option("-o", "--output",
                  action="store", dest="output_filename", 
                  help="final coadded output image")
    
    parser.add_option("-p", "--pixel_scale",
                  action="store", dest="pixel_scale", type=float, default="0.45",
                  help="Pixel scale of the image")
    
    parser.add_option("-H", "--flip_horizontally",
                  action="store_true", dest="flipH", default=False,
                  help="Flip Horizontally the image")

    parser.add_option("-V", "--flip_vertically",
                  action="store_true", dest="flipH", default=False,
                  help="Flip Vertically the image")
                  
                                
    (options, args) = parser.parse_args()
    
    
    ## Logging setup
    FORMAT =  "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    logging_level = logging.DEBUG
    logging.basicConfig(format = FORMAT, level = logging_level)
    logging.debug("Logging setup done !")
    
    
    # args is the leftover positional arguments after all options have been processed
    if not options.source_file  or len(args)!=0: 
        parser.print_help()
        parser.error("incorrect number of arguments " )
    
    # Check if source_file is a FITS file or a text file listing a set of files
    if os.path.exists(options.source_file):
        if os.path.isfile(options.source_file):
            try:
                hdulist = pyfits.open(options.source_file)
                filelist = [options.source_file]
            except:
                filelist = [line.replace( "\n", "") 
                            for line in fileinput.input(options.source_file)]
        elif os.path.isdir(options.source_file):
            filelist = glob.glob(options.source_file+"/*.fit*")
                        
        for file in filelist:
            try:
                solveField(file, options.pixel_scale, "/tmp")
            except Exception,e:
                logging.error("Error solving file %s  [%s] "%(file,str(e)))
                    
    else:
        loggging.error("Source file %s does not exists",options.source_file)
        
    
    sys.exit()
        
    
