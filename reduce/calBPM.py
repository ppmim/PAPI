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
# Created    : 25/06/2009    jmiguel@iaa.es
# Last update: 25/06/2009    jmiguel@iaa.es
#              03/03/2010    jmiguel@iaa.es Added READMODE checking
#              21/02/2014    jmiguel@iaa.es Added support for MEF files 
# 
# TODO:
#  NOT USED !!!! 
# - include master dark subtraction !!!
# - NO FUNCIONA BIEN, muy conservador ??! no da mascaras buenas !! 
#   da pocos pixeles malos ????
################################################################################


# Import necessary modules

import datetime
import os
import sys
import fileinput
from optparse import OptionParser


# Pyraf modules
import pyraf
from pyraf import iraf
from iraf import noao
from iraf import imred
from iraf import ccdred
from iraf import mscred

import pyfits
import numpy

# Logging
from misc.paLog import log

import datahandler
import misc.fileUtils
import misc.utils as utils

class BadPixelMask(object):
    """
        Generate a bad pixel mask from a list of dark corrected dome flat images
        (extracted from VIRCAM pipeline, vircam_genbpm)
    
        A list of dark corrected dome flat images is given. A master flat 
        is created from all the input flats in the list. Each input flat is 
        then divided by the master. Bad pixels are marked on the new image as 
        those that are above or below the threshold (in sigma) in the new image. 
        Any pixel which has been marked as bad for more than a quarter of the 
        input images is defined as bad in the output mask. 
        
    """
    def __init__(self, input_file, outputfile=None, lthr=3.0, hthr=5.0, 
                 temp_dir="/tmp", raw_flag=False):
        
        self.input_file = input_file # file with the list of files to read and process 
        # Default parameters values
        self.lthr = float(lthr)
        self.hthr = float(hthr)
        self.temp_dir = temp_dir
        self.raw_flag = raw_flag

        if outputfile==None:
            dt = datetime.datetime.now()
            self.output = self.temp_dir + 'BadPixMask'+dt.strftime("-%Y%m%d%H%M%S")
        else:
            self.output = outputfile
        
    
    def create(self):
        """ 
         Algorith to create the BPM
         -------------------------- 
         1. Combine all of the dome flats into a master
         2. Divide the resulting image by its median -->normalized MASTER_FLAT
         3. Create and zero the rejection mask
         4. Loop for all input images and divide each by the master flat
            
            4.1 Divide the resulting image by its median
            4.2 Get the standard deviation of the image
            4.3 Define the bad pixels
            
         5. Go through the rejection mask and if a pixel has been marked bad 
            more than a set number of times, then it is defined as bad
        """
        
        t = utils.clock()
        t.tic()
        
        
        # Read the file list
        filelist = [line.replace( "\n", "") for line in fileinput.input(self.input_file)]
        
        # Here we could check if each frame is a good dome flat !!!
        good_flats = []
        f_readmode = -1
        if self.raw_flag:
            for iframe in filelist:
                fits = datahandler.ClFits(iframe)
                log.debug("Frame %s EXPTIME= %f TYPE= %s " %(iframe, fits.exptime, fits.type)) 
                # Check EXPTIME, TYPE (flat) and FILTER
                if not fits.isDomeFlat():
                    log.warning("Warning: Task 'create BPM' found a non-domeflat frame")
                else:
                    # Check READMODE
                    if ( f_readmode!=-1 and (f_readmode!= fits.getReadMode() )):
                        log.error("Error: Task 'calBPM' finished. Found a FLAT frame with different  READMODE")
                        raise Exception("Found a FLAT frame with different  READMODE") 
                    else: 
                        f_readmode  =fits.getReadMode()
                        good_flats.append(iframe)
        else: 
            good_flats = filelist                    
             
        if len(good_flats)<2:
            log.error('Not enough dome flats provided. At least 2 good flat frames are required')
            raise Exception("Not enough dome flats provided. At least 2 good flat frames are required")
        

        # Due a bug in PyRAF that does not allow a long list of files separated with commas as 'input' argument
        # we need to build again a text file with the good_files
        if len(good_flats)!=len(filelist):
            flats = self.temp_dir + "/flats.txt"
            ftemp = open(flats,"w")
            for flat in good_flats:
                ftemp.write(flat+"\n")
            ftemp.close() 
        else: flats = self.input_file
        
        # STEP 1: Make the combine of dome Flat frames
        # - Build the frame list for IRAF
        log.debug("Combining Flat frames...")
        flat_comb = self.temp_dir + '/flatcomb.fits'
        misc.fileUtils.removefiles(flat_comb)
        # Call IRAF task (it works with MEF or simple images)
        iraf.mscred.flatcombine(input=("'"+"@"+flats+"'").replace('//','/'), 
                        output=flat_comb, 
                        combine='median', 
                        ccdtype='', 
                        process='no', 
                        reject='sigclip', 
                        subsets='no', 
                        scale='mode'
                        )
        log.debug("Flatcombine created %s"%flat_comb)

        # STEP 2: Compute normalized flat 
        # - Divide the resulting combined flat by its median (robust estimator)
        log.debug("Divide the resulting combined flat by its median...")
        flat_comb_hdu = pyfits.open(flat_comb)
        nExt = 1 if len(flat_comb_hdu)==1 else len(flat_comb_hdu)-1
        if nExt ==1:
            median = numpy.median(flat_comb_hdu[0].data)
            flat_comb_hdu[0].data = flat_comb_hdu[0].data / median
        else:
            for i_nExt in range(0, nExt):
                median = numpy.median(flat_comb_hdu[i_nExt+1].data)
                flat_comb_hdu[i_nExt+1].data = flat_comb_hdu[i_nExt+1].data / median

        if nExt==1:            
            nx1 = flat_comb_hdu[0].header['NAXIS1']
            nx2 = flat_comb_hdu[0].header['NAXIS2']
        else:
            # we suppose all extension have the same shape
            nx1 = flat_comb_hdu[1].header['NAXIS1']
            nx2 = flat_comb_hdu[1].header['NAXIS2']


        # STEP 3: Create and zero the rejection mask
        bpm = numpy.zeros([nExt, nx1, nx2], dtype=numpy.uint8)


        # A TEST !!!!
        #low = self.lthr
        #high = self.hthr
        #bpm[ 0, (flat_comb_hdu[0].data < low) | (flat_comb_hdu[0].data > high)]+=1

        #import pdb
        #pdb.set_trace()
        

        # STEP 4: Loop for all input images and divide each by the master
        tmpf = numpy.zeros([nx1, nx2], dtype=numpy.float32)
        for flat in []:#good_flats:
            log.debug("Processing file %s"%flat)
            f_i = pyfits.open(flat)
            #ceros=(mflat[0].data==0).sum()
            #print "CEROS=", ceros
            for i_nExt in range(0, nExt):
                if nExt==1:
                    # to avoid zero division error
                    mydata = numpy.where(flat_comb_hdu[0].data==0, 0.0001, 
                                         flat_comb_hdu[0].data) 
                    tmpf = f_i[0].data/mydata
                else:
                    # to avoid zero division error
                    mydata = numpy.where(flat_comb_hdu[i_nExt+1].data==0, 0.0001, 
                                         flat_comb_hdu[i_nExt+1].data) 
                    tmpf = f_i[i_nExt+1].data/mydata

                mdat = numpy.ma.masked_array(tmpf, numpy.isnan(tmpf))
                std = numpy.std(mdat)
                median = numpy.median(tmpf)
                mad = numpy.median(numpy.abs(tmpf - median))
                mad*=1.4826
            
                print ">>MEDIAN=",median
                print ">>STD=",std
                print ">>MAD=",mad
            
                #log.debug("Divide each flatted flat by its median")
                tmpf = tmpf/median
                
                # Debug - save the normalized flat
                #misc.fileUtils.removefiles(flat.replace(".fits","_N.fits"))
                #hdulist = pyfits.HDUList()
                #hdr0 = pyfits.getheader(good_flats[0])
                #prihdu = pyfits.PrimaryHDU (data = tmpf, header = hdr0)
                #hdulist.append(prihdu)
      
                #hdulist.writeto(flat.replace(".fits","_N.fits"))
                #hdulist.close(output_verify='ignore')
                # End debug

                low = 1.0 - self.lthr*mad/median
                high = 1.0 + self.hthr*mad/median
                
                #low = self.lthr
                #high = self.hthr
                    
                print ">>LOW=", low
                print ">>HIGH=", high
            
                # STEP 4.3 Define the bad pixels
                tmpf.shape = nx1,nx2
                bpm[ i_nExt, (tmpf < low) | (tmpf > high)]+=1
            log.debug("BPM updated with current flat %s", flat)
            f_i.close()
            
        flat_comb_hdu.close()
        
        # STEP 5: Go through the rejection mask and if a pixel has been marked bad 
        # more than a set number of times (a quarter of number of images), 
        # then it is defined as bad.
        nbmax = numpy.maximum(2, len(good_flats)/4)
        #nbmax = 0 ### TEST !!!
        bpm = numpy.where(bpm>nbmax, 1, 0) # bad pixel set to 1

        # Show stats
        nbad = numpy.zeros(nExt)
        for i_nExt in range(0, nExt):
            nbad[i_nExt] = (bpm[i_nExt]==1).sum()
            badfrac = float(nbad[i_nExt])/float(nx1*nx2)
            log.info("# Bad pixels (extension %s): %f"%(i_nExt, nbad[i_nExt]))
            log.info("Fraction Bad pixel (extesion %s): %f"%(i_nExt,badfrac))
        
        # STEP 6: Save the BPM --- TODO MEF !!!!
        misc.fileUtils.removefiles(self.output)
        hdulist = pyfits.HDUList()     
        hdr0 = pyfits.getheader(good_flats[0])
        prihdu = pyfits.PrimaryHDU (data = None, header = None)
        try:
            prihdu.header.set('INSTRUME', hdr0['INSTRUME'])
            prihdu.header.set('TELESCOP', hdr0['TELESCOP'])
            prihdu.header.set('CAMERA', hdr0['CAMERA'])
            prihdu.header.set('MJD-OBS', hdr0['MJD-OBS'])
            prihdu.header.set('DATE-OBS', hdr0['DATE-OBS'])
            prihdu.header.set('DATE', hdr0['DATE'])
            prihdu.header.set('UT', hdr0['UT'])
            prihdu.header.set('LST', hdr0['LST'])
            prihdu.header.set('ORIGIN', hdr0['ORIGIN'])
            prihdu.header.set('OBSERVER', hdr0['OBSERVER'])
        except Exception,e:
            log.warning("%s"%str(e))

        prihdu.header.set('PAPITYPE','MASTER_BPM','TYPE of PANIC Pipeline generated file')
        prihdu.header.add_history('BPM created from %s' % good_flats)

        if nExt>1:
            prihdu.header.set('EXTEND', pyfits.TRUE, after = 'NAXIS')
            prihdu.header.set('NEXTEND', nExt)
            prihdu.header.set('FILENAME', self.output)
            hdulist.append(prihdu)
            for i_ext in range(0, nExt):
                hdu = pyfits.PrimaryHDU()
                hdu.scale('int16') # important to set first data type
                hdu.data = bpm[i_ext]
                hdulist.append(hdu)
                del hdu
        else:
            prihdu.scale('int16') # important to set first data type
            prihdu.data = bpm[0]
            hdulist.append(prihdu)
         
        
        # write FITS
        try:
            hdulist.writeto(self.output)
            hdulist.close(output_verify='ignore')
        except Exception,e:
            log.error("Error writing linearity model %s"%self.output)
            raise e

        # Remove temp files
        #misc.fileUtils.removefiles(flat_comb)
        
        log.debug('Saved Bad Pixel Mask  to %s' , self.output)
        log.debug("createBPM' finished %s", t.tac() )
        
    def create_JM(self):
            """ 
             Algorith to create the BPM
             -------------------------- 
             1. Combine all of the dome flats into a master
             2. Divide the resulting image by its median -->normalized MASTER_FLAT
             3. Create and zero the rejection mask
             4. Loop for all input images and divide each by the master flat
                
                4.1 Divide the resulting image by its median
                4.2 Get the standard deviation of the image
                4.3 Define the bad pixels
                
             5. Go through the rejection mask and if a pixel has been marked bad 
                more than a set number of times, then it is defined as bad
            """
            
            t = utils.clock()
            t.tic()
            
            
            # Read the file list
            filelist = [line.replace( "\n", "") for line in fileinput.input(self.input_file)]
            
            # Here we could check if each frame is a good dome flat !!!
            good_flats = []
            f_readmode = -1
            if self.raw_flag:
                for iframe in filelist:
                    fits = datahandler.ClFits(iframe)
                    log.debug("Frame %s EXPTIME= %f TYPE= %s " %(iframe, fits.exptime, fits.type)) 
                    # Check EXPTIME, TYPE (flat) and FILTER
                    if not fits.isDomeFlat():
                        log.warning("Warning: Task 'create BPM' found a non-domeflat frame")
                    else:
                        # Check READMODE
                        if ( f_readmode!=-1 and (f_readmode!= fits.getReadMode() )):
                            log.error("Error: Task 'calBPM' finished. Found a FLAT frame with different  READMODE")
                            raise Exception("Found a FLAT frame with different  READMODE") 
                        else: 
                            f_readmode  =fits.getReadMode()
                            good_flats.append(iframe)
            else: 
                good_flats = filelist                    
                 
            if len(good_flats)<2:
                log.error('Not enough dome flats provided. At least 2 good flat frames are required')
                raise Exception("Not enough dome flats provided. At least 2 good flat frames are required")
            

            # Due a bug in PyRAF that does not allow a long list of files separated with commas as 'input' argument
            # we need to build again a text file with the good_files
            if len(good_flats)!=len(filelist):
                flats = self.temp_dir + "/flats.txt"
                ftemp = open(flats,"w")
                for flat in good_flats:
                    ftemp.write(flat+"\n")
                ftemp.close() 
            else: flats = self.input_file
            
            # STEP 1: Make the combine of dome Flat frames
            # - Build the frame list for IRAF
            log.debug("Combining Flat frames...")
            flat_comb = self.temp_dir + '/flatcomb.fits'
            misc.fileUtils.removefiles(flat_comb)
            # Call IRAF task (it works with MEF or simple images)
            iraf.mscred.flatcombine(input=("'"+"@"+flats+"'").replace('//','/'), 
                            output=flat_comb, 
                            combine='median', 
                            ccdtype='', 
                            process='no', 
                            reject='sigclip', 
                            subsets='no', 
                            scale='mode'
                            )
            log.debug("Flatcombine created %s"%flat_comb)

            # STEP 2: Compute normalized flat 
            # - Divide the resulting combined flat by its median (robust estimator)
            log.debug("Divide the resulting combined flat by its median...")
            flat_comb_hdu = pyfits.open(flat_comb)
            nExt = 1 if len(flat_comb_hdu)==1 else len(flat_comb_hdu)-1
            if nExt ==1:
                median = numpy.median(flat_comb_hdu[0].data)
                flat_comb_hdu[0].data = flat_comb_hdu[0].data / median
            else:
                for i_nExt in range(0, nExt):
                    median = numpy.median(flat_comb_hdu[i_nExt+1].data)
                    flat_comb_hdu[i_nExt+1].data = flat_comb_hdu[i_nExt+1].data / median

            if nExt==1:            
                nx1 = flat_comb_hdu[0].header['NAXIS1']
                nx2 = flat_comb_hdu[0].header['NAXIS2']
            else:
                # we suppose all extension have the same shape
                nx1 = flat_comb_hdu[1].header['NAXIS1']
                nx2 = flat_comb_hdu[1].header['NAXIS2']

            # STEP 3: Create and zero the rejection mask
            bpm = numpy.zeros([nExt, nx1, nx2], dtype=numpy.uint8)

            # STEP 4: Loop for all input images and divide each by the master
            tmpf = numpy.zeros([nx1, nx2], dtype=numpy.float32)
            for flat in  good_flats:
                log.debug("Processing file %s"%flat)
                f_i = pyfits.open(flat)
                #ceros=(mflat[0].data==0).sum()
                #print "CEROS=", ceros
                for i_nExt in range(0, nExt):
                    if nExt==1:
                        # to avoid zero division error
                        mydata = numpy.where(flat_comb_hdu[0].data==0, 0.0001, 
                                             flat_comb_hdu[0].data) 
                        tmpf = f_i[0].data/mydata
                    else:
                        # to avoid zero division error
                        mydata = numpy.where(flat_comb_hdu[i_nExt+1].data==0, 0.0001, 
                                             flat_comb_hdu[i_nExt+1].data) 
                        tmpf = f_i[i_nExt+1].data/mydata
                    
                    std = numpy.std(tmpf)
                    median = numpy.median(tmpf)
                    mad = numpy.median(numpy.abs(tmpf - median))
                    mad*=1.4826
                
                    print ">>MEDIAN=",median
                    print ">>STD=",std
                    print ">>MAD=",mad
                
                    #log.debug("Divide each flatted flat by its median")
                    tmpf = tmpf/median
                    
                    low = self.lthr
                    high = self.hthr
                    #low = 1.0 - self.lthr*mad/median
                    #high = 1.0 + self.hthr*mad/median
                
                    print ">>LOW=", low
                    print ">>HIGH=", high
                
                    # STEP 4.3 Define the bad pixels
                    tmpf.shape = nx1,nx2
                    bpm[ i_nExt, (tmpf < low) | (tmpf > high)]+=1
                log.debug("BPM updated with current flat %s", flat)
                f_i.close()
                
            flat_comb_hdu.close()
            
            # STEP 5: Go through the rejection mask and if a pixel has been marked bad 
            # more than a set number of times (a quarter of number of images), 
            # then it is defined as bad.
            nbmax = numpy.maximum(2, len(good_flats)/4)
            bpm = numpy.where(bpm>nbmax, 1, 0) # bad pixel set to 1

            # Show stats
            nbad = numpy.zeros(nExt)
            for i_nExt in range(0, nExt):
                nbad[i_nExt] = (bpm[i_nExt]==1).sum()
                badfrac = float(nbad[i_nExt])/float(nx1*nx2)
                log.info("# Bad pixels (extension %s): %f"%(i_nExt, nbad[i_nExt]))
                log.info("Fraction Bad pixel (extesion %s): %f"%(i_nExt,badfrac))
            
            # STEP 6: Save the BPM --- TODO MEF !!!!
            misc.fileUtils.removefiles(self.output)
            hdulist = pyfits.HDUList()     
            hdr0 = pyfits.getheader(good_flats[0])
            prihdu = pyfits.PrimaryHDU (data = None, header = None)
            try:
                prihdu.header.set('INSTRUME', hdr0['INSTRUME'])
                prihdu.header.set('TELESCOP', hdr0['TELESCOP'])
                prihdu.header.set('CAMERA', hdr0['CAMERA'])
                prihdu.header.set('MJD-OBS', hdr0['MJD-OBS'])
                prihdu.header.set('DATE-OBS', hdr0['DATE-OBS'])
                prihdu.header.set('DATE', hdr0['DATE'])
                prihdu.header.set('UT', hdr0['UT'])
                prihdu.header.set('LST', hdr0['LST'])
                prihdu.header.set('ORIGIN', hdr0['ORIGIN'])
                prihdu.header.set('OBSERVER', hdr0['OBSERVER'])
            except Exception,e:
                log.warning("%s"%str(e))

            prihdu.header.set('PAPITYPE','MASTER_BPM','TYPE of PANIC Pipeline generated file')
            prihdu.header.add_history('BPM created from %s' % good_flats)

            if nExt>1:
                prihdu.header.set('EXTEND', pyfits.TRUE, after = 'NAXIS')
                prihdu.header.set('NEXTEND', nExt)
                prihdu.header.set('FILENAME', self.output)
                hdulist.append(prihdu)
                for i_ext in range(0, nExt):
                    hdu = pyfits.PrimaryHDU()
                    hdu.scale('int16') # important to set first data type
                    hdu.data = bpm[i_ext]
                    hdulist.append(hdu)
                    del hdu
            else:
                prihdu.scale('int16') # important to set first data type
                prihdu.data = bpm[0]
                hdulist.append(prihdu)
             
            
            # write FITS
            try:
                hdulist.writeto(self.output)
                hdulist.close(output_verify='ignore')
            except Exception,e:
                log.error("Error writing linearity model %s"%self.output)
                raise e

            # Remove temp files
            misc.fileUtils.removefiles(flat_comb)
            
            log.debug('Saved Bad Pixel Mask  to %s' , self.output)
            log.debug("createBPM' finished %s", t.tac() )

###############################################################################
usage = "usage: %prog [options] "
desc = """
Generate a bad pixel mask from a list of dark corrected dome flat images.
    
A list of dark corrected dome flat images is given. A master flat 
is created from all the input flats in the list. Each input flat is 
then divided by the master. Bad pixels are marked on the new image as 
those that are above or below the threshold (in sigma) in the new image. 
Any pixel which has been marked as bad for more than a quarter of the 
input images is defined as bad in the output mask.
This module receives a series of FITS images (darks) with increasing exposure 
time and creates the master dark model and computes several statistics.
"""
parser = OptionParser(usage, description=desc)
 
               
parser.add_option("-s", "--source",
              action="store", dest="source_file_list",
              help="list of input (optionally dark corrected) dome flat images..")

parser.add_option("-o", "--output",
              action="store", dest="output_filename", 
              help="The output bad pixel mask.")

parser.add_option("-L", "--lthr",
              action="store", dest="lthr", type='float', default=3.0,
              help="The low rejection threshold in units of sigma [default 3]")

parser.add_option("-H", "--hthr",
              action="store", dest="hthr", type='float', default=5.0,
              help="The high rejection threshold in units of sigma [default 5]")

parser.add_option("-r", "--raw",
              action="store_true", dest="raw_flag", default=False,
              help="Neither check FLAT type nor Readout-mode [default False]")    

#parser.add_option("-D", "--master_dark",
#              action="store", dest="master_dark", type='str',
#              help="[Optional] Master dark frame to subtract")    

#parser.add_option("-S", "--show_stats",
#              action="store_true", dest="show_stats", default=False,
#              help="Show statistics [default False]")    



################################################################################
# main
def main(arguments=None):
    
    if arguments is None:
        arguments = sys.argv[1:] # argv[0] is the script name
    (options, args) = parser.parse_args(args = arguments)

    if len(sys.argv[1:])<1:
       parser.print_help()
       return 2

    if len(args) !=0:
        parser.print_help()
        return 2 # used for command line syntax errors
    
    # Check mandatory arguments
    if not options.output_filename or not options.source_file_list:
        parser.print_help()
        parser.error("incorrect number of arguments " )
        return 2
        
    # Make sure we are not overwriting an existing file 
    if os.path.exists(options.output_filename):
        print "Error. The output file '%s' already exists."  % \
              (options.output_filename)
        return 1
    #if options.master_dark:
    #    print "Sorry, dark subtraction not yet implemented."
    #    return 1
    
    try:
        bpm = BadPixelMask(options.source_file_list, options.output_filename, 
                       options.lthr, options.hthr, raw_flag=options.raw_flag)
        bpm.create()
    except Exception, e:
        log.error("Error running BPM: %s"%str(e))
        return 0
        
###############################################################################
if __name__ == "__main__":
    print 'Start BadPixelMap....'
    sys.exit(main())
        
