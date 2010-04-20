#! /usr/bin/env python

# Copyright (c) 2009 Jose M. Ibanez. All rights reserved.
# Institute of Astrophysics of Andalusia, IAA-CSIC
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
# 
# TODO:
# - include master dark subtraction !!!
# - NO FUNCIONA BIEN, muy conservador ??!!!! no da mascaras buenas !!!! da pocos pixeles malos ????
################################################################################


################################################################################
# Import necessary modules

import datetime
import getopt
import os
import sys
import fileinput



# Pyraf modules
import pyraf
from pyraf import iraf
from iraf import noao
from iraf import imred
from iraf import ccdred

import pyfits
import numpy

# Logging
from misc.paLog import log

import datahandler
import misc.fileUtils
import misc.utils as utils

class BadPixelMask:
    """
    \brief 
        Generate a bad pixel mask from a list of dark corrected dome flat images
        (extracted from VIRCAM pipeline, vircam_genbpm)
    
    \par Class:
         BadPixelMask   
    \par Purpose:
        Work out the BPM
    \par Description:
            
    \par Language:
        Python
    \param data
        A input image
    \retval 0
        If no error, the number of bad pixels
    \author
        JMIbannez, IAA-CSIC
        
    """
    def __init__(self, input_file, output, lthr=1.0, hthr=3.0, verbose=False):
        
        self.input_file=input_file
        # Default parameters values
        self.lthr = float(lthr)
        self.hthr = float(hthr)
        self.output=output
        self.verbose= False
        
    
    def create(self):
        """ 
         Algorith to create the BPM
         -------------------------- 
         1. Combine all of the dome flats into a master
         2. Divide the resulting image by its median
         3. Create and zero the rejection mask
         4. Loop for all input images and divide each by the master flat
            
            4.1 Divide the resulting image by its median
            4.2 Get the standard deviation of the image
            4.3 Define the bad pixels
            
         5. Go through the rejection mask and if a pixel has been marked bad 
            more than a set number of times, then it is defined as bad
        """
        
        t=utils.clock()
        t.tic()
        
        
        # Read the file list
        filelist=[line.replace( "\n", "") for line in fileinput.input(self.input_file)]
        
        # Here we could check if each frame is a good dome flat !!!
        good_flats=[]
        f_readmode=-1
        for iframe in filelist:
            fits=datahandler.ClFits(iframe)
            log.debug("Frame %s EXPTIME= %f TYPE= %s " %(iframe, fits.exptime, fits.type)) 
            # Check EXPTIME, TYPE (flat) and FILTER
            if not fits.isDomeFlat():
                log.warning("Warning: Task 'create BPM' found a non dome flat frame")
            else:
                # Check READMODE
                if ( f_readmode!=-1 and (f_readmode!= fits.getReadMode() )):
                    log.error("Error: Task 'calBPM' finished. Found a FLAT frame with different  READMODE")
                    raise Exception("Found a FLAT frame with different  READMODE") 
                else: 
                    f_readmode  =fits.getReadMode()
                    good_flats.append(iframe)
                
            
        if len(good_flats)<2:
            log.error('Not enought dome flats provided. At least 2 good flat frames are required')
            sys.exit(1)
        
        # Due a bug in PyRAF that does not allow a long list of files separated with commas as 'input' argument
        # we need to build again a text file with the good_files
        if len(good_flats)!=len(filelist):
            ftemp=open("/tmp/flats.txt","w")
            for flat in good_flats:
                ftemp.write(flat+"\n")
            ftemp.close()       
            flats="/tmp/flats.txt"
        else: flats=self.input_file
             
        # STEP 1: Make the combine of dome Flat frames
        # - Build the frame list for IRAF
        log.debug("Combining Flat frames...")
        flat_comb='/tmp/flatcomb.fits'
        misc.fileUtils.removefiles(flat_comb)
        # Call IRAF task
        iraf.flatcombine(input='@'+flats, 
                        output=flat_comb, 
                        combine='median', 
                        ccdtype='none', 
                        process='no', 
                        reject='sigclip', 
                        subsets='no', 
                        scale='mode'
                        )
                        
        # STEP 2: Divide the resulting combined flat by its median
        # Compute the mean of the image
        log.debug("Divide the resulting combined flat by its median...")
        median = float(iraf.imstat (
            images=flat_comb,
            fields='midpt',Stdout=1)[1])
        
        # Compute normalized flat
        iraf.imarith(operand1=flat_comb,
                    operand2=median,
                    op='/',
                    result=flat_comb,
                    )
                    
        f=pyfits.open(flat_comb)
        nx1=f[0].header['NAXIS1']
        nx2=f[0].header['NAXIS1']
        f.close()
                    
        # STEP 3: Create and zero the rejection mask
        bpm=numpy.zeros([nx1, nx2])
        
        # STEP 4: Loop for all input images and divide each by the master
        mflat=pyfits.open(flat_comb)
        tmpf=numpy.zeros([nx1, nx2])
        for flat in  good_flats:
            f_i=pyfits.open(flat)
            #ceros=(mflat[0].data==0).sum()
            #print "CEROS=", ceros
            mydata=numpy.where(mflat[0].data==0, 0.0001, mflat[0].data) # to avoid zero division error
            tmpf=f_i[0].data/mydata
            std=numpy.std(tmpf)
            tmpf.shape=nx1*nx2
            median=numpy.median(tmpf)
            
            print "MEDIAN=",median
            print "STD=",std
            
            log.debug("Divide each flatted flat by its median")
            tmpf=tmpf/median
            
            low = 1.0 - self.lthr*std/median
            high = 1.0 + self.hthr*std/median
            
            print "LOW=", low
            print "HIGH=", high
            
            #STEP 4.3 Define the bad pixels
            tmpf.shape=nx1,nx2
            bpm[ (tmpf < low) | (tmpf > high)]+=1
            log.debug("BPM updated with current flat %s", flat)
            f_i.close()
            
        mflat.close()
        # STEP 5: Go through the rejection mask and if a pixel has been marked bad 
        # more than a set number of times, then it is defined as bad
        nbmax=numpy.max(2, len(good_flats)/4)
            
        bpm=numpy.where(bpm>nbmax,1,0) # bad pixel set to 1
        nbad=(bpm==1).sum()
        
        badfrac = float(nbad)/float(nx1*nx2)
        print "# bad pixels", nbad
        print " BAD_FRAC", badfrac
        
        # Save the BPM
        misc.fileUtils.removefiles( self.output )               
        hdu = pyfits.PrimaryHDU()
        hdu.scale('int16') # importat to set first data type
        hdu.data=bpm     
        hdulist = pyfits.HDUList([hdu])
        hdu.header.update('OBJECT','MASTER_PIXEL_MASK')
        hdu.header.add_history('BPM created from %s' % good_flats)
        hdulist.writeto(self.output)
        hdulist.close(output_verify='ignore')
        
        
        # Remove temp files
        #misc.fileUtils.removefiles(flat_comb)
        
        log.debug('Saved Bad Pixel Mask  to %s' , self.output)
        log.debug("createBPM' finished %s", t.tac() )
        

#-----------------------------------------------------------------------

def usage():
    print ''
    print 'NAME'
    print '       calBPM.py - Bad Pixel Mask creation\n'
    print 'SYNOPSIS'
    print '       calBPM.py [options] -f file.list\n'
    print 'DESCRIPTION'
    print '       Compute a BPM from a give dome-flat frames list'
    print ' '
    print 'OPTIONS'
    print '       -f --file          The list of input (dark corrected) dome flat images'
    print '       -d --dark          The master dark to subtract to dome flats'
    print '       -o --out bmp.fits  The output bad pixel mask'
    print '       -l --lthr 3        The low rejection threshold in units of sigma'
    print '       -h --hthr 5        The high rejection threshold in units of sigma'
    print '       -v : verbose debugging output\n'
    print 'VERSION'
    print '       25 June 2009'
    print ''
    raise SystemExit

#-----------------------------------------------------------------------




################################################################################
# main
if __name__ == "__main__":
    print 'Start BadPixelMask....'
    
    # Read command line parameters
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'f:d:o:l:h:v',['file=','dark=','out=','lthr=','hthr=','verbose'])
    except getopt.GetoptError:
        usage()
        sys.exit(1)
    
    nargs = len(sys.argv[1:])
    nopts = len(opts)
    if nargs<2:
        usage()
        sys.exit(2)  
    
    verbose= False
    inputfile=''
    dark=''
    outputfile='/tmp/bpm.fits'
    lthr=3
    hthr=5
    
            
    for option, par in opts:
        if option in ("-f", "--file"):
            inputfile=par
            print "inputfile=", inputfile
        if option in ("-d", "--dark"):
            dark=par
            print "dark=", dark
        if option in ("-o", "--out"):
            outputfile=par
            print "outputfile=",outputfile
        if option in ("-l", "--lthr"):
            lthr=par
            print "lthr=", lthr
        if option in ("-h", "--hthr"):
            hthr=par
            print "hthr=",hthr
        if option in ('-v','--verbose'):      # verbose debugging output
            verbose = True
            print "Verbose true"
                
    
    # Error checking:
    if not os.path.exists(inputfile):      # check whether input file exists
        print inputfile, 'does not exist'
        sys.exit(2)
    
    print '...reading', inputfile
    
    bpm = BadPixelMask(inputfile, outputfile, lthr, hthr, verbose)
    bpm.create()
    
    print 'ending BPM....'




    
