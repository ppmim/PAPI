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
# calBPM.py
#
# Created    : 25/09/2009    jmiguel@iaa.es
# Last update: 25/09/2009    jmiguel@iaa.es
#              19/04/2010    jmiguel@iaa.es - added master dark checking
#
# TODO:
#   - scale master dark before subtraction
################################################################################


# Import necessary modules
import os
import datetime
import sys
import fileinput
from optparse import OptionParser

from scipy import signal
#comentado el signal porque en el portatil da un error !!!!!

# Pyraf modules
import pyraf
from pyraf import iraf
from iraf import noao
from iraf import mscred

import pyfits
import numpy

# Logging
from misc.paLog import log

import datahandler
import misc.fileUtils
import misc.utils as utils

class ExError(Exception):
    pass


class BadPixelMask(object):
    """
    Build a bad pixel mask (hot and cold pixels) from a set of images 
    ( lamp_on and lamp_off frames and darks) using iraf.ccdmask task 
        
     
     Algorithm to create the BPM
     -------------------------- 
     1: classify/split the frames in 3 sets (DOME_FLAT_LAMP_ON, DOME_FLAT_LAMP_OFF, DARKS)
        and and check whether there are enough calib frames
     2: Check the master dark (Texp)
     3: Subtract the master dark to each dome flat
     4: Combine dome dark subtracted flats (on/off)
     5: Compute flat_low/flat_high
     6: Create BPM (iraf.ccdmask)

        
    Returns
    -------
    If no error, the bad pixel mask (fits file)
        
    """
    def __init__(self, i_file_list, master_dark, output_file=None, lsigma=20, 
                 hsigma=20, temp_dir="/tmp/"):
        
        self.i_file_list = i_file_list
        self.master_dark = master_dark
        # Default parameters values
        self.lsigma = float(lsigma)
        self.hsigma = float(hsigma)
        self.temp_dir = temp_dir
        
        if output_file==None:
            dt = datetime.datetime.now()
            self.output_file = self.temp_dir + 'BadPixMask'+dt.strftime("-%Y%m%d%H%M%S")
        else:
            self.output_file = output_file
                
        self.verbose = False
    
    def create(self):
        self.create_IRAF()
        
    def create_IRAF(self):
        
        """
        Generating a Bad Pixel Mask
        They are generated by taking flat field exposures at high and low 
        intensity levels, and identifying those pixels where the ratio of the 
        normalised flux is more than 15%. In detail (using IRAF commands):
        
            1-Take: 5 dark images, 5 Domeflats with low counts (~5000 counts), 
                    5 Domeflats with high counts (~40,000 counts)
            2-Zerocombine the darks images
            3-dark subtract the Domeflat images
            4-Flatcombine the domeflats with high counts - FLATHIGH
            5-Flatcombine the domeflats with low counts - FLATLOW
            6-iraf.imarith/iraf.mscred.mscarith FLATLOW.fits / FLATHIGH.fits FLATLHRATIO.fits
            7-ccdmask FLATLHRATIO.fits BADPIX_MASK_FILE with the following parameters set:
        
         
        
                            Image Reduction and Analysis Facility
        PACKAGE = ccdred
           TASK = ccdmask
        
        image   =                       Input image
        mask    =                       Output pixel mask
        (ncmed  =                    7) Column box size for median level calculation
        (nlmed  =                    7) Line box size for median level calculation
        (ncsig  =                   15) Column box size for sigma calculation
        (nlsig  =                   15) Line box size for sigma calculation
        (lsigma =                  20.) Low clipping sigma
        (hsigma =                  20.) High clipping sigma
        (ngood  =                    5) Minimum column length of good pixel seqments
        (linterp=                    2) Mask value for line interpolation 
        (cinterp=                    3) Mask value for column interpolation 
        (eqinter=                    2) Mask value for equal interpolation
        (mode   =                   ql)

        
        Returns
        -------

        The list of IRAF mask files given with a pixel list image (.pl).
        Pixel list images have the extension ".pl" signifying a special compact 
        file of integer values ideal for identifying sets of pixels. For a bad 
        pixel mask the good pixels have a value of zero and bad pixels have 
        positive integer values. 
        
        ==> Values are 0 for good pixels, non-zero for bad (masked) pixels.


        Notes
        -----
        See also:

        iraf.fixpix -- fix pixels identified by a bad pixel mask, image, or file

        """
        
        log.debug('createBadPixelMask started (iraf.ccdmask)')
        t = utils.clock()
        t.tic()
        
        flats_off_frames = []
        flats_on_frames = []
    
        
        # Check Master_Dark
        if os.path.exists(self.master_dark):
            f = datahandler.ClFits(self.master_dark)
            if not f.getType()=='MASTER_DARK':
                log.error("File %s does not look a MASTER_DARK! Check your data.", self.master_dark)
                raise ExError("File does not look a MASTER_DARK! Check your data.")              
        else:
            log.error("File %s does not exist", self.master_dark)
            raise ExError("File does not exist")    
        
        # Read the file list
        filelist = [line.replace( "\n", "").replace("//","/") 
                    for line in fileinput.input(self.i_file_list)]
        
        
        #STEP 1: classify/split the frames in 3 sets (DOME_FLAT_LAMP_ON, DOME_FLAT_LAMP_OFF)
        #and create string list for IRAF tasks
        for file in filelist:
            f = datahandler.ClFits(file)
            if f.isDomeFlatOFF():
                flats_off_frames.append(f.pathname)
            elif f.isDomeFlatON():
                flats_on_frames.append(f.pathname)
            else:
                # reject the frame
                log.error("DISCARTING: Frame %s is not dome_flat. Type = %s"%(f.pathname, f.getType()))
                
        # Check whether there are enough calib frames
        if (len(flats_off_frames)<1 or len(flats_on_frames)<1 or abs(len(flats_off_frames)-len(flats_off_frames))>10):
            log.error("There are not enough calib frames to create BPM !!")
            raise ExError("Not enough calib frames")
        

        #STEP 2: Subtrac the master dark to each dome flat
        for file in flats_off_frames:
            misc.fileUtils.removefiles(file.replace(".fits","_D.fits"))
            iraf.mscred.mscarith(operand1=file,
                        operand2=self.master_dark,
                        op='-',
                        result=file.replace(".fits","_D.fits"),
                        )
            flats_off_frames[flats_off_frames.index(file)] = file.replace(".fits","_D.fits")
        
        
        for file in flats_on_frames:
            misc.fileUtils.removefiles(file.replace(".fits","_D.fits"))
            iraf.mscred.mscarith(operand1=file,
                        operand2=self.master_dark,
                        op='-',
                        result=file.replace(".fits","_D.fits"),
                        )
            flats_on_frames[flats_on_frames.index(file)] = file.replace(".fits","_D.fits")                                                

        #STEP 4: Combine dome dark subtracted flats (OFF)
        flat_off_comb = self.temp_dir + "flats_comb_off.fits"
        misc.fileUtils.removefiles(flat_off_comb)
        # Due a bug in PyRAF that does not allow a long list of files separated 
        # with commas as 'input' argument, we need to build again a text file 
        # with the files.
        flats = self.temp_dir + "/flats_off.txt"
        misc.utils.listToFile(flats_off_frames, flats)

        iraf.mscred.flatcombine(input="@"+flats.replace('//','/'),
                        output=flat_off_comb,
                        combine='median',
                        ccdtype='',
                        process='no',
                        reject='sigclip',
                        subset='no',
                        scale='mode',
                        #verbose='yes'
                        #scale='exposure',
                        #expname='EXPTIME'
                        #ParList = _getparlistname ('flatcombine')
                        )

        flat_on_comb = self.temp_dir + "flats_comb_on.fits"
        misc.fileUtils.removefiles(flat_on_comb)
        flats = self.temp_dir + "flats_on.txt"
        misc.utils.listToFile(flats_on_frames, flats)
        
        #Combine dome dark subtracted flats (on)
        iraf.mscred.flatcombine(input="@"+flats.replace('//','/'),
                        output=flat_on_comb,
                        combine='median',
                        ccdtype='',
                        process='no',
                        reject='sigclip',
                        subset='no',
                        scale='mode',
                        #verbose='yes'
                        #scale='exposure',
                        #expname='EXPTIME'
                        #ParList = _getparlistname ('flatcombine')
                        )

        misc.fileUtils.removefiles(flats)
                                    
        #STEP 5: Compute flat_low/flat_high
        flat_ratio = self.temp_dir + 'flat_ratio.fits'
        misc.fileUtils.removefiles(flat_ratio)
        log.debug( 'Computing FlatRatio')
        iraf.mscred.mscarith(operand1=flat_off_comb,
                 operand2=flat_on_comb,
                 op='/',
                 result=flat_ratio,
                 )
        log.debug("Flat Ratio built")
                     
        #STEP 6: Create BPM
        log.debug( 'Now, computing bad pixel mask')
        misc.fileUtils.removefiles(self.output_file + ".pl")

        try:
            my_hdu = pyfits.open(flat_ratio)
            nExt = 1 if len(my_hdu)==1 else len(my_hdu)-1
        except Exception,e:
            log.error("Cannot read file: %"%flat_ratio)
            raise Exception("Cannot read file: %"%flat_ratio)

        # iraf.ccdmask does not support MEF files, so we have to split them and
        # run ccdmask for each extension. Values is 0 for good pixels, non-zero 
        # for bad (masked) pixels.
        list_outfitsnames = []
        if nExt==1:
            iraf.ccdmask(image=flat_ratio,
                     mask=self.output_file.partition(".fits")[0]+".pl",
                     ncmed = 7, # Column box size for median level calculation
                     nlmed = 7, # Line box size for median level calculation
                     ncsig = 15, # Column box size for sigma calculation
                     nlsig = 15, # Line box size for sigma calculation
                     lsigma = self.lsigma, # Low clipping sigma
                     hsigma = self.hsigma, # High clipping sigma
                     ngood = 5, # Minimum column length of good pixel seqments
                     linterp = 2, # Mask value for line interpolation
                     cinterp = 3, # Mask value for column interpolation
                     eqinter = 2 # Mask value for equal interpolation
                     )
            list_outfitsnames.append(self.output_file.partition(".fits")[0]+".pl")
        else:
            i = 0
            for i_nExt in range(0,nExt):
                mfnp = self.output_file.partition('.fits')
                outfitsname = mfnp[0] + ".Q%02d"%i + mfnp[1] + mfnp[2] 
                iraf.ccdmask(image=flat_ratio+"[%d]"%(i+1),
                         mask=outfitsname.partition(".fits")[0]+".pl",
                         ncmed = 7, # Column box size for median level calculation
                         nlmed = 7, # Line box size for median level calculation
                         ncsig = 15, # Column box size for sigma calculation
                         nlsig = 15, # Line box size for sigma calculation
                         lsigma = self.lsigma, # Low clipping sigma
                         hsigma = self.hsigma, # High clipping sigma
                         ngood = 5, # Minimum column length of good pixel seqments
                         linterp = 2, # Mask value for line interpolation
                         cinterp = 3, # Mask value for column interpolation
                         eqinter = 2 # Mask value for equal interpolation
                         )
                i = i +1 
                list_outfitsnames.append(outfitsname.partition(".fits")[0]+".pl")
            
        # Convert to IRAF (.pl) bpm to .fits and update header
        for mask_file in list_outfitsnames:
            new_file = mask_file.partition(".pl")[0]+".fits"
            iraf.imcopy(mask_file, new_file)
            pyfits.setval(new_file, keyword="PAPITYPE", 
                                             value="MASTER_BPM")            
            pyfits.setval(new_file, keyword="HISTORY",
                                        value="BPM created from %s"%filelist)
            log.info("Bad Pixel Mask file created : %s"%new_file)
        
        # Clean up tmp files
        #misc.fileUtils.removefiles(flat_off_comb, flat_on_comb, flat_ratio)
        # Change back to the original working directory
        iraf.chdir()
        
        log.debug("Time elapsed : [%s]" , t.tac() )

        return list_outfitsnames
        
    def create_simple(self):
        """
        Build a BPM following the algorithm described above.

        Returns
        -------

        - Image file (FITS) with 0 for good pixels, non-zero for bad (masked) 
        pixels.
        
        """
        
        log.debug('createBadPixelMask started (simple mode)')
        t = utils.clock()
        t.tic()
        
        flats_off_frames = []
        flats_on_frames = []

        # Check Master_Dark
        if os.path.exists(self.master_dark):
            f = datahandler.ClFits(self.master_dark)
            if not f.getType()=='MASTER_DARK':
                pass
                #log.error("File %s does not look a MASTER_DARK! Check your data.", self.master_dark)
                #raise ExError("File does not look a MASTER_DARK! Check your data.")              
        else:
            log.error("File %s does not exist", self.master_dark)
            raise ExError("File does not exist")    
        
        
        # Read the source file list
        filelist=[line.replace( "\n", "").replace("//","/") for line in fileinput.input(self.i_file_list)]
        
        #STEP 1: classify/split the frames in 3 sets (DOME_FLAT_LAMP_ON, DOME_FLAT_LAMP_OFF, DARKS)
        # and create string list for IRAF tasks
        for file in filelist:
            f = datahandler.ClFits(file)
            if f.getType()=='DOME_FLAT_LAMP_OFF':
                flats_off_frames.append(f.pathname)
            elif f.getType()=='DOME_FLAT_LAMP_ON':
                flats_on_frames.append(f.pathname)
            else:
                # reject the frame
                log.error("DISCARTING: Frame %s is neither dome_flat nor a dark frame", f.pathname)
                
        # Check whether there are enough calib frames
        if (len(flats_off_frames)<1 or len(flats_on_frames)<1 or abs(len(flats_off_frames)-len(flats_off_frames))>10 
            or self.master_dark==None):
            log.error("There are not enough calib frames for create BPM !!")
            raise ExError("Not enough calib frames")
        
        
        # STEP 2: Subtrac the master dark to each dome flat
        # NOTE: all dome flats (on and off) should have the same EXPTIME 
        for file in flats_off_frames:
            misc.fileUtils.removefiles(file.replace(".fits","_D.fits"))
            iraf.mscred.mscarith(operand1=file,
                        operand2=self.master_dark,
                        op='-',
                        result=file.replace(".fits","_D.fits"),
                        )
            flats_off_frames[flats_off_frames.index(file)]=file.replace(".fits","_D.fits")
        
        for file in flats_on_frames:
            misc.fileUtils.removefiles(file.replace(".fits","_D.fits"))
            iraf.mscred.mscarith(operand1=file,
                        operand2=self.master_dark,
                        op='-',
                        result=file.replace(".fits","_D.fits"),
                        )
            flats_on_frames[flats_on_frames.index(file)] = file.replace(".fits","_D.fits")                                                

        print "OFF=", flats_off_frames
        print "ON=", flats_on_frames
        
        # STEP 3: Combine dome dark subtracted flats (off & on)
        flat_off_comb = self.temp_dir + "flats_comb_off.fits"
        misc.fileUtils.removefiles(flat_off_comb)
        # Due to a bug in PyRAF that does not allow a long list of files separated with commas as 'input' argument
        # we need to build again a text file with the files
        flats = self.temp_dir + "/flats_off.txt"
        misc.utils.listToFile(flats_off_frames, flats)
        #Combine dome dark subtracted OFF-flats
        iraf.mscred.flatcombine(input="@"+flats.replace('//','/'),
                        output=flat_off_comb.replace('//','/'),
                        combine='median',
                        ccdtype='none',
                        process='no',
                        reject='sigclip',
                        subset='no',
                        scale='mode'
                        #verbose='yes'
                        #scale='exposure',
                        #expname='EXPTIME'
                        #ParList = _getparlistname ('flatcombine')
                        )            
        #misc.fileUtils.removefiles(flats)
        
        flat_on_comb = self.temp_dir + "flats_comb_on.fits"
        misc.fileUtils.removefiles(flat_on_comb)
        flats = self.temp_dir + "/flats_on.txt"
        misc.utils.listToFile(flats_on_frames, flats)

        # Combine dome dark subtracted ON-flats
        iraf.mscred.flatcombine(input="@"+flats.replace('//','/'),
                        output=flat_on_comb.replace('//','/'),
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
        misc.fileUtils.removefiles(flats)
                                    
        # STEP 4: Compute flat_comb_off/flat_comb_on
        flat_ratio = self.temp_dir + 'flat_ratio.fits'
        misc.fileUtils.removefiles(flat_ratio)
        log.debug( 'Computing FlatRatio')
        iraf.mscred.mscarith(operand1=flat_off_comb.replace('//','/'),
                 operand2=flat_on_comb.replace('//','/'),
                 op='/',
                 result=flat_ratio,
                 )
        log.debug("Flat Ratio built")
                     
        # STEP 5: Create BPM
        log.debug( 'Now, computing bad pixel mask')
        misc.fileUtils.removefiles(self.output_file)
        
        fr = pyfits.open(flat_ratio)
        
        # TBC: Bad pixel selection criterion 
        # Bad pixels >0, good_pixels=0
        bpm = numpy.where(fr[0].data<0.9, 1, 
                           numpy.where(fr[0].data>1.1, 1, 0))
        #bpm=numpy.where(bpm.isfinite(), 1, 0)

        # -- a test with a gauss kernel, it takes a long time ...(endless?) 
        #g = gauss_kern(5, sizey=2048)
        #log.debug( 'Gauss kernel done !')
        #improc = signal.convolve(bpm, g, mode='full')
        #log.debug( 'Convolution done !')
        #bpm2 = numpy.where (improc<0.55, 0, 1)
        fr.close()
        bpm2 = bpm
         
        # Save the BPM
        hdu = pyfits.PrimaryHDU()
        hdu.data = bpm2     
        hdu.scale('int16') # importat to set first data type
        hdulist = pyfits.HDUList([hdu])
        hdu.header.update('PAPITYPE','MASTER_BPM', 'TYPE of PANIC Pipeline generated file')
        hdulist.writeto(self.output_file)
        hdulist.close(output_verify='ignore')
        
        # Clean up tmp files
        #misc.fileUtils.removefiles(flat_off_comb, flat_on_comb, flat_ratio)
        # Change back to the original working directory
        iraf.chdir()
        
        log.info('Bad pixel mask created : %s', self.output_file)
        log.debug("Time elapsed : [%s]" , t.tac() )

        return self.output_file


def fixPix(image, mask):
    """
    Return an image with masked values replaced with a bi-linear
    interpolation from nearby pixels.  Probably only good for isolated
    badpixels.

    Usage:
      fixed = fixpix(im, mask, [iraf=])

    Inputs:
      image = the image 2D array
      mask = an array that is True where im contains bad pixels

    Outputs:
      fixed = the corrected image

    v1.0.0 Michael S. Kelley, UCF, Jan 2008

    v1.1.0 Added the option to use IRAF's fixpix.  MSK, UMD, 25 Apr
           2011
    """
    if iraf:
        # importing globally is causing trouble
        from pyraf import iraf as IRAF

        badfits = os.tmpnam() + '.fits'
        outfits = os.tmpnam() + '.fits'
        pyfits.writeto(badfits, mask.astype(np.int16))
        pyfits.writeto(outfits, im)
        IRAF.fixpix(outfits, badfits)
        cleaned = pyfits.getdata(outfits)
        os.remove(badfits)
        os.remove(outfits)
        return cleaned

    interp2d = interpolate.interp2d
    #     x = xarray(im.shape)
    #     y = yarray(im.shape)
    #     z = im.copy()
    #     good = (mask == False)
    #     interp = interpolate.bisplrep(x[good], y[good],
    #                                         z[good], kx=1, ky=1)
    #     z[mask] = interp(x[mask], y[mask])

    # create domains around masked pixels
    dilated = ndimage.binary_dilation(mask)
    domains, n = ndimage.label(dilated)

    # loop through each domain, replace bad pixels with the average
    # from nearest neigboors
    x = xarray(image.shape)
    y = yarray(image.shape)
    cleaned = image.copy()
    for d in (np.arange(n) + 1):
        # find the current domain
        i = (domains == d)

        # extract the sub-image
        x0, x1 = x[i].min(), x[i].max() + 1
        y0, y1 = y[i].min(), y[i].max() + 1
        subim = image[y0:y1, x0:x1]
        submask = mask[y0:y1, x0:x1]
        subgood = (submask == False)

        cleaned[i * mask] = subim[subgood].mean()

    return cleaned

#-----------------------------------------------------------------------

def gauss_kern(size, sizey=None):
    """ Returns a normalized 2D gauss kernel array for convolutions """
    size = int(size)
    if not sizey:
        sizey = size
    else:
        sizey = int(sizey)
    x, y = numpy.mgrid[-size:size+1, -sizey:sizey+1]
    g = numpy.exp(-(x**2/float(size) + y**2/float(sizey)))
    return g / g.sum()

def blur_image(im, n, ny=None) :
    """ blurs the image by convolving with a gaussian kernel of typical
        size n. The optional keyword argument ny allows for a different
        size in the y direction.
    """
    g = gauss_kern(n, sizey=ny)
    improc = signal.convolve(im, g, mode='valid')
    return(improc)



###############################################################################
usage = "usage: %prog [options] "
descr = u'''Build a bad pixel mask (hot and cold pixels) from a set of images '''+ \
        '''(lamp_ON and lamp_OFF frames and master dark)'''
exam = "Example: %prog -s /tmp/input_files.txt -o /tmp/bpm.fits -D /tmp/masterDark.fits"
parser = OptionParser(usage, description=descr, epilog = exam)
 
               
parser.add_option("-s", "--source",
              action="store", dest="source_file_list",
              help="list of input (optionally  corrected) dome ON and OFF flat images..")

parser.add_option("-o", "--output",
              action="store", dest="output_filename", 
              help="The output bad pixel mask.")

parser.add_option("-L", "--lthr",
              action="store", dest="lthr", type='float', default=20.0,
              help="The low rejection threshold in units of sigma [default 20]")

parser.add_option("-H", "--hthr",
              action="store", dest="hthr", type='float', default=20.0,
              help="The high rejection threshold in units of sigma [default 20]")

parser.add_option("-D", "--master_dark",
              action="store", dest="master_dark", type='str',
              help="[Optional] Master dark frame to subtract")    




################################################################################
# main
def main(arguments=None):
    
    if arguments is None:
        arguments = sys.argv[1:] # argv[0] is the script name
    (options, args) = parser.parse_args(args = arguments)

    if len(sys.argv[1:])<1 or len(args)!=0:
       parser.print_help()
       return 1

    if options.source_file_list is None or not os.path.exists(options.source_file_list):
        parser.print_help()
        parser.error("Wrong source file given")
        
    # Make sure we are not overwriting an existing file 
    if options.output_filename and os.path.exists(options.output_filename):
        print "Error. The output file '%s' already exists."  % \
              (options.output_filename)
        return 1
    if not options.master_dark or not os.path.exists(options.master_dark):
        print "Error, The master dark frame does not exists."
        return 1
    
    try:
        bpm = BadPixelMask(options.source_file_list, options.master_dark,
                           options.output_filename, 
                           options.lthr, options.hthr)
        bpm.create_IRAF()
    except Exception, e:
        log.error("Error creating BPM %s"%str(e))
        return 0
###############################################################################
if __name__ == "__main__":
    print 'Start BadPixelMap....'
    sys.exit(main())
