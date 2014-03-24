#!/usr/bin/env python

################################################################################
#
# symple test programm for PyRAF
#
# mkBadPix.py
#
# Last update 05/Feb/2008
#
################################################################################
################################################################################
#
# INPUT:
# Needs 3 files containing 
#                          - list of flats with high counts
#                          - list of flats with low counts
#
# DESCRIPTION:
# - Processes raw images to make 1 badpix files as follow:
#    >FLATLHRATIO.fits= FLATLOW.fits / FLATHIGH.fits
#    >ccdmask FLATLHRATIO.fits
################################################################################


################################################################################

# Import necessary modules


import os
import sys
import time
#import datetime

from pyraf import iraf
from iraf import noao
from iraf import imred
from iraf import ccdred

wk_dir='/disk-a/caha/panic/DATA/test1/'

def compute_badPixMask ( flat_files, threshold ):

    flat_frames=''
    nlow=0
    nhigh=0
    flats_low=''
    flats_high=''
    
    #Change to the source directory
    iraf.chdir(wk_dir)
    
    # Compound list of files for IRAF 
    for nframe in flat_files:
        #Compute the mean of the image
        mean=float(iraf.imstat (
            images=nframe,
            fields='mean',Stdout=1)[1])
        print 'Frame ' + nframe + ' MEAN value %f' % mean

        if ( mean < threshold ):
            flats_low+=nframe+ ' , '
            nlow+=1
        else:
            flats_high+=nframe+' , '
            nhigh+=1
            
    #Check whether there are enough flats files
    if ( nlow<1 or nhigh<1 or abs(nlow-nhigh)>40 ):
        #Error, not enough flat frames
        print 'Error in compute_badPixMask, not enough flats frames '
        return -1
    
    #Combine LOW flats
    print 'Combining LOW flats : ' + flats_low
    iraf.flatcombine(input=flats_low,
                     output=wk_dir+'FlatLow.fits',
                     combine='median',
                     ccdtype='none',
                     process='no',
                     reject='sigclip',
                     scale='exposure',
                     #expname='EXPTIME'
                     #ParList = _getparlistname ('flatcombine')
                     )
    
    #Combine HIGH flats
    print 'Combining HIGH flats : ' + flats_high
    iraf.flatcombine(input=flats_high,
                     output=wk_dir+'FlatHigh.fits',
                     combine='median',
                     ccdtype='none',
                     process='no',
                     reject='sigclip',
                     scale='exposure',
                     #expname='EXPTIME'
                     #ParList = _getparlistname ('flatcombine')
                     )

   
    #Make flat_low/flat_high
    print 'Computing FlatRatio'
    iraf.imarith(operand1='FlatLow.fits',
                 operand2='FlatHigh.fits',
                 op='/',
                 result='FlatRatio.fits',
                 )    

    #Determine mask
    print 'Determining bad pixel mask'
    t=datetime.datetime.now()
    iraf.ccdmask(image=wk_dir+'FlatRatio.fits',
                 mask='BadPixMask'+t.strftime("-%Y%m%d%H%M%S"),
                 lsigma=20,
                 hsigma=20,
                 )
    
    # Change back to the original working directory
    iraf.chdir()

    print 'Bad pixel mask created'

    
    






