#!/usr/bin/env python

# Copyright (c) 2011-2012 IAA-CSIC  - All rights reserved. 
# Author: Jose M. Ibanez. 
# Instituto de Astrofisica de Andalucia, IAA-CSIC
#
# This file is part of PAPI (PANIC Pipeline)
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

################################################################################
#
# PANICtool
#
# applyDarkFlat.py
#
# Created    : 05/06/2009    jmiguel@iaa.es
# Last update: 05/06/2009    jmiguel@iaa.es
#              11/03/2010    jmiguel@iaa.es - Added out_dir for output files
#              18/03/2010    jmiguel@iaa.es - Modified to support only dark 
#                            subtraction, only flatfielding or both 
#              15/11/2010    jmiguel@iaa.es - Added normalization to flat-field
#              16/11/2010    jmiguel@iaa.es - Added support for MEF files
#              21/03/2014    jmiguel@iaa.es - Added support for BPM
#              18/03/2015    jmiguel@iaa.es - Added checking of NCOADD
################################################################################

################################################################################
# Import necessary modules

import getopt
import sys
import os
import fileinput
import time
from optparse import OptionParser
from scipy import interpolate, ndimage

import misc.fileUtils
import misc.statutils
import misc.utils as utils
import misc.robust as robust
import datahandler

# Logging
from misc.paLog import log
from misc.version import __version__

# Interact with FITS files
import astropy.io.fits as fits

import misc.cleanBadPix as cleanBadPix
import numpy


class ExError (Exception): 
    """ Next class if for a general execution Exception """
    pass

class ApplyDarkFlat(object):
    """
    Class used to subtract a master dark (or dark model), divide 
    by a master flat field and apply a BPM.
    
    Applies a master dark, master flat and BPM to a list of 
    non-calibrated science files.
    For each file processed a new file it is generated with the same 
    filename but with the suffix '_DF.fits' and/or '_BPM.fits'.
    
    Parameters
    ----------
    sci_raw_files: list
        A list of science raw files to calibrate 
    dark: str
        Master dark to subtract; it can be a master dark model to produce a 
        proper scaled master dark to subtract.
    mflat: str
        Master flat to divide by (not normalized !)
    bpm: str
        Input bad pixel mask or None

    out_dir: str
        Output directory where calibrated files will be created
    
    bpm_action: str
        Action to perform with BPM:
            - fix, to fix Bad pixels
            - grab, to set to 'NaN' Bad pixels.
            - none, nothing to do with BPM (default)
    
    force_apply: bool
        If true, no checking with data FITS header will be done (IMAGETYPE).  

    norm: bool
        If true, perform Flat-Field normalization (wrt median chip_1).

    Returns
    -------
    file_list
        If no error, return the list of files generated as result of the current 
        processing
  
    """
    def __init__(self, sci_raw_files, mdark = None, mflat = None,  bpm = None,
                 out_dir_ ="/tmp", bpm_action='none', force_apply=False,
                 norm = False):

        self.__sci_files = sci_raw_files  # list of files which apply dark and flat
        self.__mdark = mdark  # master dark (or master model) to apply
        self.__mflat = mflat  # master flat to apply (dome, twlight) - not normalized !!
        self.__bpm = bpm      # bad pixel mask file to apply
        self.__out_dir = out_dir_ # output directory of calibrated files
        self.__bpm_action = bpm_action
        self.__force_apply = force_apply 
        self.__norm = norm   
    
    def apply(self):
      
        """
        Applies masters DARK and/or FLAT and/or BPM to science file list. 
        Both master DARK and FLAT are optional,i.e., each one can be applied 
        even the other is not present.
        
        If dark_EXPTIME matches with the sci_EXPTIME, a straight subtraction is 
        done, otherwise, a Master_Dark_Model is required in order to compute a 
        scaled dark.

        Note: This routine works fine with MEF files and data cubes, of a with
        MEF-cubes. If means that if sci_data is a cube (3D array), the dark is 
        subtracted to each layer, and flat is applied also to each layer. Thus,
        the calibrations works correctly for data cubes of data, 
        no matter if they are MEF or single HDU fits.   
        """   

        log.debug("Start applyDarkFlat")
        log.debug(" + Master Dark: %s" % self.__mdark)
        log.debug(" + Master Flat: %s" % self.__mflat)
        log.debug(" + Master BPM: %s" % self.__bpm)
        log.debug(" + SCI Source images: %s" % str(self.__sci_files))
        
        
        start_time = time.time()
        t = utils.clock()
        t.tic()
    
        # some variables 
        cdark = None
        cflat = None
        out_suffix = ".fits"
        dmef = False # flag to indicate if dark is a MEF file or not
        fmef = False # flag to indicate if flat is a MEF file or not
        n_ext = 1 # number of extension of the MEF file (1=simple FITS file)
        median = 1 # median of the flat frame (if mef, mode of chip 0) 
        dark_time = None
        dark_ncoadd = None
        # default value used if normalization is not done
        
        # #######################
        # Master DARK reading:
        if self.__mdark != None:
            if not os.path.exists(self.__mdark): # check whether input file exists
                msg = "Master Dark '%s' does not exist"%self.__mdark
                log.error(msg)
                raise Exception(msg)
            else:
                dark = fits.open(self.__mdark)
                cdark = datahandler.ClFits (self.__mdark)
                dark_time = cdark.expTime()
                dark_ncoadd = cdark.getNcoadds()
                if (not self.__force_apply and cdark.getType() != 'MASTER_DARK' and 
                    cdark.getType() != 'MASTER_DARK_MODEL'):
                    log.error("File %s does not look a neither MASTER_DARK nor MASTER_DARK_MODEL",self.__mdark)
                    raise Exception("File %sdoes not look a neither MASTER_DARK nor MASTER_DARK_MODEL"%self.__mdark)
                
                n_ext = cdark.next  
                out_suffix = out_suffix.replace(".fits","_D.fits")
                log.debug("Master DARK to be subtracted : %s"%self.__mdark)
        # No master dark provided, then no dark to subtract
        else:
            log.warning("No master dark to be subtracted !")
            dark_data = 0
            dark_time = None
            dark_ncoadd = None
            
        # ######################
        # Master FLAT reading
        if self.__mflat != None:
            if not os.path.exists(self.__mflat): # check whether input file exists
                msg = "Master Flat '%s' does not exist"%self.__mflat
                log.error(msg)
                raise Exception(msg)
            else:
                flat = fits.open(self.__mflat)
                cflat = datahandler.ClFits (self.__mflat)
                
                if not cflat.isMasterFlat() and not self.__force_apply:
                    log.error("File %s does not look a neither MASTER_FLAT",self.__mflat)
                    raise Exception("File %s does not look MASTER_FLAT"%self.__mflat)
                
                flat_filter = cflat.getFilter()
                # check MEF compatibility
                if cdark != None and cdark.next != cflat.next:
                    raise Exception("Number of extensions does not match "
                        "in Dark and Flat files!")
                else: n_ext = cflat.next
                log.debug("Flat MEF file with %d extensions", n_ext)
                
                # Flat Field Normalization (in principle, input Flat must be already normalized):
                # compute mode for n_ext normalization (in case of MEF, we normalize 
                # all extension wrt chip 0)
                naxis1 = cflat.getNaxis1()
                naxis2 = cflat.getNaxis2()
                if cflat.next > 1: 
                    ext = 1
                else: 
                    ext = 0
                    
                # Take the center of the image
                off_naxis1 = int(naxis1 * 0.1)
                off_naxis2 = int(naxis2 * 0.1)
                dat = flat[ext].data[off_naxis2: naxis2 - off_naxis2, 
                                     off_naxis1: naxis1 - off_naxis1]
                # NaN values must not be replaced with 0.0 !!!
                # dat[numpy.isnan(dat)]= 0.0 #

                # Normalization is done with a robust estimator --> np.median()
                median = robust.r_nanmedian(dat)
                mean = robust.r_nanmean(dat)
                mode = 3 * median - 2 * mean
                
                log.info("Flat stats: MEDIAN= %f  MEAN=%f MODE(estimated)=%f ", \
                           median, mean, mode)
                
                if self.__norm: 
                    log.warning("Flat-field will be normalized by MEDIAN ( %f ) value", median)
                
                out_suffix = out_suffix.replace(".fits","_F.fits")
                log.debug("Found master FLAT to divide by: %s" % self.__mflat)
        else:
            log.warning("No master flat to be divided by !")
            flat_data = 1.0
            flat_filter = None
            modian = 1 
        
                
        # Add suffix whether BPM is to be applied
        if self.__bpm != None and self.__bpm_action != 'none':
            out_suffix = out_suffix.replace(".fits","_BPM.fits") 
            if self.__mdark == None and self.__mflat == None:
                with fits.open(self.__bpm) as bp:
                    n_ext = len(bp) - 1
            

        # Get the user-defined list of flat frames
        framelist = self.__sci_files
        
        # STEP 2: Check the TYPE and FILTER of each science file
        # If any frame on list missmatch the FILTER, then the procedure will be 
        # aborted. 
        # EXPTIME do not need be the same, so EXPTIME scaling will be done
        # However, if force_apply==True no checking will be done.
        n_removed = 0
        
        # List of files generated as result of this procedure and that will be returned
        result_file_list = [] 
        
        #
        # Start the applying of calibrations
        #
        for iframe in framelist:
            if not os.path.exists(iframe):
                log.error("File '%s' does not exist", iframe)
                continue  
            f = fits.open(iframe)
            cf = datahandler.ClFits(iframe)
            f_ncoadd = cf.getNcoadds()
            log.debug("Science frame %s, EXPTIME = %f, TYPE = %s, FILTER = %s, NCOADD = %s"\
                      %(iframe, cf.expTime(), cf.getType(), cf.getFilter(), f_ncoadd))
            
            # Check FILTER
            if (not self.__force_apply and flat_filter != None 
                and cf.getFilter() != flat_filter):
                log.error("Error: Task 'applyDarkFlat' found a frame with " 
                "different FILTER. Skipping frame...")
                f.close()
                n_removed = n_removed + 1
            else:
                # check Number of Extensions
                if (len(f) > 1 and (len(f) - 1) != n_ext):
                    raise Exception("File %s does not match the number of "
                    "extensions (%d)"%( iframe, n_ext))
                elif len(f) == 1 and n_ext != 1: 
                    raise Exception("File %s does not match the number of "
                    "extensions (%d)"%( iframe, n_ext))
                else:
                    log.debug("Good match of the number of extensions.")

                # Delete old files
                (path, name) = os.path.split(iframe)
                newpathname = (self.__out_dir + "/" + \
                             name.replace(".fits", out_suffix)).replace("//","/")
                misc.fileUtils.removefiles(newpathname)
                
                # Scale master DARK
                exp_time = float(cf.expTime()) # all extension have the same TEXP
                if self.__mdark != None and dark_time != None :
                    time_scale = float(exp_time / dark_time)
                else: time_scale = 1.0
                
                for chip in range(0, n_ext):
                    dark_data = None
                    # MEF
                    if n_ext > 1: # it means, MEF
                        # Get the right extension by name
                        chip_name = 'Q%d' % (chip + 1)
                        try:
                            f[chip_name].header
                        except KeyError:
                            chip_name = 'SG%i_1' % (chip + 1)
                        
                        log.info("Processing extension %s" % chip_name)
                        
                        # Get DARK
                        if self.__mdark != None: 
                            if (not self.__force_apply and 
                                (not numpy.isclose(time_scale, 1.0, atol=1e-02) 
                                 or f_ncoadd != dark_ncoadd)
                                ): # for dark_model time_scale==-1
                                log.debug("Dark EXPTIME or NCOADD mismatch ! looking for DarkModel ...")
                                if not cdark.isMasterDarkModel():
                                    log.error("Cannot find out a scaled dark to apply")
                                    raise Exception("Cannot find a scaled dark to apply")
                                else:
                                    log.debug("Scaling dark with DarkModel...")
                                    dark_data = dark[chip_name].data[1] * exp_time + dark[chip_name].data[0]
                            else:
                                dark_data = dark[chip_name].data
                        else: 
                            dark_data = 0
                    
                        # Get normalized FLAT
                        if self.__mflat != None: 
                            if self.__norm:
                                # normalization wrt chip 0
                                flat_data = flat[chip_name].data / median
                            else:
                                # suppose it's already normalized
                                flat_data = flat[chip_name].data 
                        else: 
                            flat_data = 1
                            
                        # Get SCI data 
                        sci_data = f[chip_name].data

                        # Get BPM
                        if self.__bpm != None: 
                            bpm_data = fits.getdata(self.__bpm, extname = chip_name,
                                                      header=False)

                    # Single
                    else:
                        # Get DARK
                        if self.__mdark != None:
                            if (not self.__force_apply and 
                                (not numpy.isclose(time_scale, 1.0, atol=1e-02) 
                                 or f_ncoadd != dark_ncoadd)
                                ): # for dark_model time_scale==-1
                                
                                if f_ncoadd != dark_ncoadd:
                                    log.warning("Dark NCOADD mismatch !. Checking if Dark is a dark model...")
                                else:
                                    log.warning("Dark EXPTIME mismatch (time_scale= %f)! Checking if Dark is a dark model ..." % time_scale)
                                
                                if not cdark.isMasterDarkModel():
                                    log.error("Dark is not a DarkModel, cannot find out a scaled dark to apply")
                                    raise Exception("Cannot find a scaled dark to apply")
                                else:
                                    log.debug("DarkModel found: Scaling dark with dark model...")
                                    dark_data = dark[0].data[1] * exp_time + dark[0].data[0]
                                    log.info("AVG(scaled_dark)=%s"% robust.r_nanmean(dark_data))
                            else:
                                dark_data = dark[0].data
                        else: dark_data = 0
                        
                        # Get normalized FLAT
                        if self.__mflat != None: 
                            if self.__norm:
                                log.debug("Normalizing FF...")
                                flat_data = flat[0].data / median  # normalization
                            else:
                                log.debug("No normalization will be done to FlatField.")
                                # we suppose it's already normalized
                                flat_data = flat[0].data

                        else: flat_data = 1
                        
                        ## Get RAW_SCI data    
                        sci_data = f[0].data
                        ## 
                        
                        # Get BPM
                        if self.__bpm != None:
                            # bpm_data: must be an array that is True or >0 
                            # where bad pixels
                            bpm_data = fits.getdata(self.__bpm, header=False)
                                                               
                    
                    # To avoid NaN values due to zero division by FLAT
                    __epsilon = 1.0e-20
                    flat_data = numpy.where(numpy.fabs(flat_data) < __epsilon, 
                                            1.0, flat_data)
                    # Other way to solve the zero division in FF
                    #sci_data = numpy.where(flat_data==0.0, 
                    #                       (sci_data - dark_data), 
                    #                       (sci_data - dark_data) / flat_data )
                    # Store back the new values

                    ###################################
                    # Finally, apply dark, Flat and BPM
                    # #################################
                    # Note: if sci_data is a cube (3D array), dark is subtracted
                    # to each layer, and flat is applied also to each layer. Thus,
                    # the calibrations works correctly for data cubes of data, 
                    # no matter if they are MEF or single HDU fits.  
                    sci_data = (sci_data - dark_data) / flat_data
                    
                    # TEST to convert NaNs to 0's !!
                    #sci_data = misc.statutils.nan2num(sci_data, 0)
                    # end-test
                    
                    # Now, apply BPM
                    if self.__bpm != None:
                        if sci_data.shape != bpm_data.shape:
                            print "SCI = ",sci_data.shape
                            print "BPM = ",bpm_data.shape
                            log.error("Source data and BPM do not match image shape")
                            raise Exception("Source data and BPM do not match image shape")
                        
                        if self.__bpm_action == 'fix':
                            log.debug("Fixing Bad Pixles...")
                            #sci_data = fixpix(sci_data, bpm_data, iraf=True)
                            sci_data = fixpix(sci_data, bpm_data)

                        elif self.__bpm_action == 'grab':
                            log.debug("Grabbing BPM")
                            sci_data[bpm_data == 1] = numpy.NaN

                        else:
                            log.debug("Nothing to do with BPM")

                    if n_ext > 1: # it means, MEF
                        f[chip_name].data = sci_data
                    else:
                        f[0].data = sci_data

                # Update header
                if self.__mdark != None: 
                    f[0].header.add_history('Dark subtracted %s' %self.__mdark)
                if self.__mflat != None: 
                    f[0].header.add_history('Flat-Field with %s' %self.__mflat)
                if self.__bpm != None and self.__bpm_action != 'none':
                    f[0].header.add_history('BPM with %s' %self.__bpm)

                f[0].header.set('PAPIVERS', __version__, 'PANIC Pipeline version')

                # Write output to outframe (data object actually still points 
                # to input data)
                try:
                    if n_ext == 1: f[0].scale('float32')
                    else:
                         for chip in range(0, n_ext):
                             f[chip + 1].scale('float32')
                    f.writeto(newpathname, output_verify='ignore')
                except IOError:
                    raise ExError('Cannot write output to %s' % newpathname)
                     
                f.close()
                result_file_list.append(newpathname)            
                log.debug('Saved new dark subtracted or/and flattened or/and BP Masked file  %s',
                          newpathname )
        
        log.debug(t.tac() )
        log.info("Successful end of applyDarkFlat !")
                
        return result_file_list
        
def fixpix( image_data, mask_data):
    """
    Clean masked (bad) pixels from an input image. Each masked pixel 
    is replaced by the median of unmasked pixels in a 2D window of ``size`` centered on
    it.  If all pixels in the window are masked, then the window is
    increased in size until unmasked pixels are found.
    
    Parameters
    ----------
    image_data: the image array to fix
    
    mask_data: an array that is True (or >0) where image contains bad pixels
    
    Returns
    -------
    The cleaned image array;otherwise an exception is raised.
    
        
    """
    
    try:
        #mask = numpy.where( mask_data == 0, 1, 0)
        #mask = numpy.logical_not(mask_data)
        return cleanBadPix._clean_masked_pixels(image_data, mask_data, size=5, exclude_mask=None)
    except Exception,e:
        log.error("Error cleanning bad pixels...")
        raise e
    
def fixpix_old(im, mask, iraf=False):
    """ 
    Applies a bad-pixel mask to the input image (im), creating an image with 
    masked values replaced with a bi-linear interpolation from nearby pixels.  
    Probably only good for isolated badpixels.

    Usage:
      fixed = fixpix(im, mask, [iraf=])

    Inputs:
      im = the image array
      mask = an array that is True (or >0) where im contains bad pixels
      iraf = True use IRAF.fixpix; False use numpy and a loop over 
             all pixels (extremelly low)

    Outputs:
      fixed = the corrected image

    v1.0.0 Michael S. Kelley, UCF, Jan 2008

    v1.1.0 Added the option to use IRAF's fixpix.  MSK, UMD, 25 Apr
           2011

    Notes
    -----
    - Non-IRAF algorithm is extremelly slow.

    """

    log.info("Fixpix - Bad pixel mask interpolation (iraf=%s)"%iraf)
    if iraf:
        # importing globally is causing trouble
        from pyraf import iraf as IRAF
        import tempfile

        badfits = tempfile.NamedTemporaryFile(suffix=".fits").name
        outfits = tempfile.NamedTemporaryFile(suffix=".fits").name
        fits.writeto(badfits, mask.astype(numpy.int16))
        fits.writeto(outfits, im)
        IRAF.fixpix(outfits, badfits)
        cleaned = fits.getdata(outfits)
        os.remove(badfits)
        os.remove(outfits)
        return cleaned

    
    #
    # Next approach is too slow !!!
    #

    # interp2d = interpolate.interp2d
    # x = numpy.array(im.shape)
    # y = numpy.array(im.shape)
    # z = im.copy()
    # good = (mask == False)
    # interp = interpolate.bisplrep(x[good], y[good],
    #                            z[good], kx=1, ky=1)
    # z[mask] = interp(x[mask], y[mask])

    # return z
  
    # create domains around masked pixels
    dilated = ndimage.binary_dilation(mask)
    domains, n = ndimage.label(dilated)

    # loop through each domain, replace bad pixels with the average
    # from nearest neigboors
    y, x = numpy.indices(im.shape, dtype=numpy.int)[-2:]
    #x = xarray(im.shape)
    #y = yarray(im.shape)
    cleaned = im.copy()
    for d in (numpy.arange(n) + 1):
        # find the current domain
        i = (domains == d)

        # extract the sub-image
        x0, x1 = x[i].min(), x[i].max() + 1
        y0, y1 = y[i].min(), y[i].max() + 1
        subim = im[y0:y1, x0:x1]
        submask = mask[y0:y1, x0:x1]
        subgood = (submask == False)

        cleaned[i * mask] = subim[subgood].mean()

    return cleaned
        
################################################################################
# main
if __name__ == "__main__":
    
    # Get and check command-line options
    usage = "usage: %prog [options]"
    desc = """
This module receives a series of FITS source images and applies a master Dark, Flat and
BPM (subtract dark, divide Flat, and fix bad pixel) using the given
calibration files (master dark, master flat-field, bad pixel mask). 
In principle, source raw files can be MEFs and data cubes, but cannot be
mixed MEF files and non-MEF files, and of course, images must have the
same image size.
"""
    parser = OptionParser(usage, description=desc)
    
                  
    parser.add_option("-s", "--source",
                  action="store", dest="source_file_list",
                  help="Source file listing the filenames of raw frames.")
    
    parser.add_option("-d", "--dark",
                  action="store", dest="dark_file", 
                  help="Master dark to be subtracted.")
    
    parser.add_option("-f", "--flat-field",
                  action="store", dest="flat_file",
                  help="Master flat-field to be divided by.")
    
    parser.add_option("-b", "--bpm",
                  action="store", dest="bpm_file",
                  help="Master Bad Pixel Mask to use (optional)")
    
    parser.add_option('-a', '--bpm_action',
                      type='choice',
                      action='store',
                      dest='bpm_action',
                      choices=['fix', 'grab', 'none',],
                      default='none',
                      help='Action (none,grab,fix) to perform with BPM [default: %default]')

    parser.add_option("-o", "--out_dir",
                  action="store", dest="out_dir", default="/tmp/",
                  help="Directory where output files will be saved [Default: %default]")

    parser.add_option("-F", "--force_apply",
                  action="store_true", dest="force_apply", default=False,
                  help="Forces operations withouth checking FITS data type [default: %default]")
   
    parser.add_option("-N", "--normalize_FF",
                  action="store_false", dest="normalize", default=False,
                  help="Performs Flat-Filed normalization [default: %default]")

    
    (options, args) = parser.parse_args()
    
    if len(sys.argv[1:]) < 1 :
       parser.print_help()
       sys.exit(0)

    if not options.source_file_list:
        parser.print_help()
        parser.error("Incorrect number of arguments " )
    
    
    if os.path.isdir(options.source_file_list):
        parser.print_help()
        parser.error("Source must be a file, not a directory")
    
    if options.dark_file is None and options.flat_file is None \
       and options.bpm_file is None:
        parser.print_help()
        parser.error("Incorrect number of arguments " )
    
    if datahandler.isaFITS(options.source_file_list):
        filelist = [options.source_file_list]
    else:
        filelist = [line.replace( "\n", "") 
                    for line in fileinput.input(options.source_file_list)]
        
    try:
        res = ApplyDarkFlat(filelist, options.dark_file, options.flat_file, 
                            options.bpm_file, options.out_dir, 
                            options.bpm_action, options.force_apply,
                            options.normalize)
        res.apply() 
    except Exception,e:
        log.error("Error running task: %s" % str(e))
        sys.exit(0)
    
    
