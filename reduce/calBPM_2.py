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
# Created    : 25/09/2009    jmiguel@iaa.es
# Last update: 25/09/2009    jmiguel@iaa.es
#              19/04/2010    jmiguel@iaa.es - added master dark checking
#
# TODO:
#  - escalar el master dark a restar !!!!
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
        Make a bad pixel mask (hot and cold pixels) from a set of images ( lamp_on and lamp_off frames and darks)
        using iraf.ccdmask task 
            
    \par Class:
         BadPixelMask   
    \par Purpose:
        Work out the BPM
    \par Description:
         
         Algorith to create the BPM
         -------------------------- 
         1: classify/split the frames in 3 sets (DOME_FLAT_LAMP_ON, DOME_FLAT_LAMP_OFF, DARKS)
            and and check whether there are enought calib frames
         2: Check the master dark (Texp)
         3: Subtrac the master dark to each dome flat
         4: Combine dome dark subtracted flats (on/off)
         5: Compute flat_low/flat_high
         6: Create BPM (iraf.ccdmask)

            
    \par Language:
        Python
    \param data
        A list of dome flat and a master dark
    \retval 0
        If no error, the bad pixel mask (fits file)
    \author
        JMIbannez, IAA-CSIC
        
    """
    def __init__(self, i_file_list, master_dark, output_file, lsigma=4, hsigma=4):
        
        self.i_file_list=i_file_list
        self.master_dark = master_dark
        # Default parameters values
        self.lsigma = float(lsigma)
        self.hsigma = float(hsigma)
        self.output_file=output_file
        self.verbose= False
        self.outputdir="/tmp/"
    
    def create(self):
        self.create_IRAF()
        
    def create_IRAF(self):
        
        
        log.debug('createBadPixelMask started (iraf.ccdmask)')
        t=utils.clock()
        t.tic()
        
        flats_off_frames=[]
        flats_on_frames=[]
    
        
        #Check Master_Dark
        if os.path.exists(self.master_dark):
            f = datahandler.ClFits(self.master_dark)
            if not f.getType()=='MASTER_DARK':
                log.error("File %s does not look a MASTER_DARK! Check your data.", self.master_dark)
                raise ExError, "File does not look a MASTER_DARK! Check your data."              
        else:
            log.error("File %s does not exist", self.master_dark)
            raise ExError, "File does not exist"    
        
        # Read the file list
        filelist=[line.replace( "\n", "") for line in fileinput.input(self.i_file_list)]
        
        #STEP 1: classify/split the frames in 3 sets (DOME_FLAT_LAMP_ON, DOME_FLAT_LAMP_OFF)
        #        and create string list for IRAF tasks
        #Also subtract master DARK
        for file in filelist:
            f = datahandler.ClFits(file)
            if f.getType()=='DOME_FLAT_LAMP_OFF':
                flats_off_frames.append(f.pathname)
            elif f.getType()=='DOME_FLAT_LAMP_ON':
                flats_on_frames.append(f.pathname)
            else:
                # reject the frame
                log.error("DISCARTING: Frame %s is not dome_flat", f.pathname)
                
        #Check whether there are enought calib frames
        if (len(flats_off_frames)<1 or len(flats_on_frames)<1 or abs(len(flats_off_frames)-len(flats_off_frames))>10):
            log.error("There are not enought calib frames to create BPM !!")
            raise ExError, "Not enought calib frames"
        
        
        #STEP 3: Subtrac the master dark to each dome flat
        for file in flats_off_frames:
            misc.fileUtils.removefiles(file.replace(".fits","_D.fits"))
            iraf.imarith(operand1=file,
                        operand2=self.master_dark,
                        op='-',
                        result=file.replace(".fits","_D.fits"),
                        )
            flats_off_frames[flats_off_frames.index(file)]=file.replace(".fits","_D.fits")
        
        for file in flats_on_frames:
            misc.fileUtils.removefiles(file.replace(".fits","_D.fits"))
            iraf.imarith(operand1=file,
                        operand2=self.master_dark,
                        op='-',
                        result=file.replace(".fits","_D.fits"),
                        )
            flats_on_frames[flats_on_frames.index(file)]=file.replace(".fits","_D.fits")                                                
        
        print "OFF=", flats_off_frames
        print "ON=", flats_on_frames
        
        
        #STEP 4: Combine dome dark subtracted flats (off)
        
        flat_off_comb=self.outputdir+"flats_comb_off.fits"
        misc.fileUtils.removefiles(flat_off_comb)
        # Due a bug in PyRAF that does not allow a long list of files separated with commas as 'input' argument
        # we need to build again a text file with the files
        flats="/tmp/flats_off.txt"
        ftemp=open(flats,"w")
        for flat in flats_off_frames:
            ftemp.write(flat+"\n")
        ftemp.close()       
        iraf.flatcombine(input="@"+flats,
                        output=flat_off_comb,
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
        #misc.fileUtils.removefiles(flats)
        
        
        flat_on_comb=self.outputdir + "flats_comb_on.fits"
        misc.fileUtils.removefiles(flat_on_comb)
        flats="/tmp/flats_on.txt"
        ftemp=open(flats,"w")
        for flat in flats_on_frames:
            ftemp.write(flat+"\n")
        ftemp.close()       
        #Combine dome dark subtracted flats (on)
        iraf.flatcombine(input="@"+flats,
                        output=flat_on_comb,
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
        #misc.fileUtils.removefiles(flats)
                                    
        #STEP 5: Compute flat_low/flat_high
        flat_ratio = self.outputdir + 'flat_ratio.fits'
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
        if self.output_file==None:
            outF = self.outputdir + 'BadPixMask'+dt.strftime("-%Y%m%d%H%M%S")
        else:
            outF = self.output_file
        
        misc.fileUtils.removefiles(outF+".pl")
        iraf.ccdmask(image=flat_ratio,
                 mask=outF,
                 lsigma=self.lsigma,
                 hsigma=self.hsigma,
                 )
        
        # Save the BPM
        #misc.fileUtils.removefiles( self.output )               
        #hdu = pyfits.PrimaryHDU()
        #hdu.scale('int16') # important to set first data type
        #hdu.data=bpm     
        #hdulist = pyfits.HDUList([hdu])
        #hdu.header.update('OBJECT','MASTER_PIXEL_MASK')
        #hdu.header.add_history('BPM created from %s' % good_flats)
        #hdulist.writeto(self.output)
        #hdulist.close(output_verify='ignore')
        
        # Clean up tmp files
        #misc.fileUtils.removefiles(flat_off_comb, flat_on_comb, flat_ratio)
        # Change back to the original working directory
        iraf.chdir()
        
        log.info('Bad pixel mask created : %s', outF+".pl")
        log.debug("Time elapsed : [%s]" , t.tac() )

        return outF+".pl"
        
    def create_simple(self):
        
        
        log.debug('createBadPixelMask started (simple mode)')
        t=utils.clock()
        t.tic()
        
        flats_off_frames=[]
        flats_on_frames=[]
    
    
        #Check Master_Dark
        if os.path.exists(self.master_dark):
            f = datahandler.ClFits(self.master_dark)
            if not f.getType()=='MASTER_DARK':
                log.error("File %s does not look a MASTER_DARK! Check your data.", self.master_dark)
                raise ExError, "File does not look a MASTER_DARK! Check your data."              
        else:
            log.error("File %s does not exist", self.master_dark)
            raise ExError, "File does not exist"    
        
        
        # Read the file list
        filelist=[line.replace( "\n", "") for line in fileinput.input(self.i_file_list)]
        
        #STEP 1: classify/split the frames in 3 sets (DOME_FLAT_LAMP_ON, DOME_FLAT_LAMP_OFF, DARKS)
        #        and create string list for IRAF tasks
        #Also subtract master DARK
        for file in filelist:
            f = datahandler.ClFits(file)
            if f.getType()=='DOME_FLAT_LAMP_OFF':
                flats_off_frames.append(f.pathname)
            elif f.getType()=='DOME_FLAT_LAMP_ON':
                flats_on_frames.append(f.pathname)
            else:
                # reject the frame
                log.error("DISCARTING: Frame %s is neither dome_flat nor a dark frame", f.pathname)
                
        #Check whether there are enought calib frames
        if (len(flats_off_frames)<1 or len(flats_on_frames)<1 or abs(len(flats_off_frames)-len(flats_off_frames))>10 
            or self.master_dark==None):
            log.error("There are not enought calib frames for create BPM !!")
            raise ExError, "Not enought calib frames"
        
        
        #STEP 3: Subtrac the master dark to each dome flat
        for file in flats_off_frames:
            iraf.imarith(operand1=file,
                        operand2=self.master_dark,
                        op='-',
                        result=file.replace(".fits","_D.fits"),
                        )
            flats_off_frames[flats_off_frames.index(file)]=file.replace(".fits","_D.fits")
        
        for file in flats_on_frames:
            iraf.imarith(operand1=file,
                        operand2=self.master_dark,
                        op='-',
                        result=file.replace(".fits","_D.fits"),
                        )
            flats_on_frames[flats_on_frames.index(file)]=file.replace(".fits","_D.fits")                                                
        print "OFF=", flats_off_frames
        print "ON=", flats_on_frames
        
        
        #STEP 4: Combine dome dark subtracted flats (off)
        
        flat_off_comb=self.outputdir+"flats_comb_off.fits"
        misc.fileUtils.removefiles(flat_off_comb)
        # Due a bug in PyRAF that does not allow a long list of files separated with commas as 'input' argument
        # we need to build again a text file with the files
        flats="/tmp/flats_off.txt"
        ftemp=open(flats,"w")
        for flat in flats_off_frames:
            ftemp.write(flat+"\n")
        ftemp.close()       
        iraf.flatcombine(input="@"+flats,
                        output=flat_off_comb,
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
        #misc.fileUtils.removefiles(flats)
        
        
        flat_on_comb=self.outputdir + "flats_comb_on.fits"
        misc.fileUtils.removefiles(flat_on_comb)
        flats="/tmp/flats_on.txt"
        ftemp=open(flats,"w")
        for flat in flats_on_frames:
            ftemp.write(flat+"\n")
        ftemp.close()       
        #Combine dome dark subtracted flats (on)
        iraf.flatcombine(input="@"+flats,
                        output=flat_on_comb,
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
        #misc.fileUtils.removefiles(flats)
                                    
        #STEP 5: Compute flat_low/flat_high
        flat_ratio = self.outputdir + 'flat_ratio.fits'
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
        if self.output_file==None:
            outF = self.outputdir + 'BadPixMask'+dt.strftime("-%Y%m%d%H%M%S")
        else:
            outF = self.output_file
        
        misc.fileUtils.removefiles(outF)
        
        fr=pyfits.open(flat_ratio)
        bpm=numpy.where (fr[0].data<0.9, 1, numpy.where(fr[0].data>1.0, 1, 0))
        fr.close()
         
        # Save the BPM
        hdu = pyfits.PrimaryHDU()
        hdu.data=bpm     
        hdu.scale('int16') # importat to set first data type
        hdulist = pyfits.HDUList([hdu])
        hdu.header.update('OBJECT','MASTER_PIXEL_MASK')
        hdulist.writeto(outF)
        hdulist.close(output_verify='ignore')
        
        # Clean up tmp files
        #misc.fileUtils.removefiles(flat_off_comb, flat_on_comb, flat_ratio)
        # Change back to the original working directory
        iraf.chdir()
        
        log.info('Bad pixel mask created : %s', outF)
        log.debug("Time elapsed : [%s]" , t.tac() )

        return outF
#-----------------------------------------------------------------------
 
def usage():
    print ''
    print 'NAME'
    print '       calBPM_2.py - Bad Pixel Mask creation\n'
    print 'SYNOPSIS'
    print '       calBPM_2.py [options] -f file.list\n'
    print 'DESCRIPTION'
    print '       Compute a BPM from a give dome-flat (on and off) frames list'
    print ' '
    print 'OPTIONS'
    print '       -f --file          The list of input dome flat images (on and off)'
    print '       -d --dark          The master dark to subtract to dome flats'
    print '       -o --out bmp.fits  The output bad pixel mask'
    print '       -l --lsig 20       Positive sigma factors to use for selecting pixels below and above the median level based on the local percentile sigma'
    print '       -h --hsig 20        '
    print '       -v : verbose debugging output\n'
    print 'VERSION'
    print '       25 Sep 2009'
    print ''
    raise SystemExit

#-----------------------------------------------------------------------




################################################################################
# main
if __name__ == "__main__":
    print 'Start creation of BadPixelMask....'
    
    # Read command line parameters
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'f:d:o:l:h:v',['file=','dark=','out=','lsig=','hsig=','verbose'])
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
    lsig=20
    hsig=20
    
            
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
        if option in ("-l", "--lsig"):
            lthr=par
            print "lsig=", lsig
        if option in ("-h", "--hsig"):
            hthr=par
            print "hsig=",hsig
        if option in ('-v','--verbose'):      # verbose debugging output
            verbose = True
            print "Verbose true"
                
    
    # Error checking:
    if not os.path.exists(inputfile):      # check whether input file exists
        print inputfile, 'does not exist'
        sys.exit(2)
    
    print '...reading', inputfile
    
    bpm = BadPixelMask(inputfile, dark, outputfile, lsig, hsig)
    bpm.create()
    #bpm.create_simple()
    
    print 'ending BPM....'