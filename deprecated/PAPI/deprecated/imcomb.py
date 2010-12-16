#!/usr/bin/env python

################################################################################
#
#
# imcombine.py
#
# Last update 20/01/2009
#
################################################################################

import os
import glob
import shutil

# Pyraf modules
from pyraf import iraf
from iraf import noao
from iraf import imred
from iraf import ccdred

# Interact with FITS files
import pyfits

# Check the mask file exists
mask_file='mask.pl'
if not os.path.exists( mask_file ):
    print('No external Mask file found. Creating the mask file : "%s"' %mask_file)
    shutil.copy('gain.fits','mask.fits')
    #Replace all pixels values greater than or equal to 0.5 by 1
    iraf.imrepl(images='mask.fits',
                value='1',
                lower=0.5
               )
    shutil.copy('mask.fits','mask.pl')
    for filename in glob.glob('mask.fits*') : os.remove( filename ) 

#Edit headers to add the reference to the proper mask file
iraf.hedit(images='@im.list',
            fields='bmp',
            value='mask.pl',
            add='yes',
            update='yes',
            verify='no',
            show='yes'
            )

for filename in glob.glob('imc.fits*') : os.remove( filename ) 
            
iraf.imcombine(input='@im.list',
            output='imc.fits',
            combine='aver',
            offset='off.list',
            reject='sigclip',
            lsigma=2.5,
            hsigma=2.5,
            scale='none',
            zero='median',
            masktype='goodv',
            maskvalue=1
            #masktype='none'
            #scale='exposure',
            #expname='EXPTIME'
            #ParList = _getparlistname ('flatcombine')
            )
            
iraf.imarith(operand1='imc.fits',
            operand2=32768,
            op='-',
            result='imc.fits',
            verbose='yes'
            )
                
                  