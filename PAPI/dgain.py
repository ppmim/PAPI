#!/usr/bin/env python

################################################################################
#
#
# dgain.py
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


# Check the mask file exists
mask_file='dc_gain.fits'

if not os.path.exists( mask_file ):
    print('Error, No external Mask file found: "%s"' %mask_file)
else:
    iraf.imrepl(images=mask_file,
                value='0',
                lower='INDEF',
                upper=0.5
               )
