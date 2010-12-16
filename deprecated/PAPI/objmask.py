#!/usr/bin/env python

################################################################################
#
#
# objmask.py
#
# Last update 21/01/2009
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


# Check the mask file exists
file='mask.pl.fits'

if not os.path.exists( file ):
    print('Warning, No external mask file %s found. Creating ...' %file)
    shutil.copy('gain.fits','mask.fits')
    iraf.imrepl(images='mask.fits',
                value='1',
                lower=0.5
               )
    shutil.copy('mask.fits', file)
    for filename in glob.glob('mask.fits*') : os.remove( filename )
    

iraf.imarith(operand1="@objfiles.nip",
            operand2=file,
            op='*',
            result="@objfiles.nip",
            verbose='yes'
            )
