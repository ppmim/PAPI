#!/usr/bin/env python

################################################################################
#
#
# addmask.py
#
# Last update 20/01/2009
#
################################################################################

# Pyraf modules
import os
import sys
from pyraf import iraf
from iraf import noao
from iraf import imred
from iraf import ccdred


# Check the badmask file exists
bpm_file=sys.argv[1]
if not os.path.exists( bpm_file ):
    print('No external Bad Pixel Mask found. Cannot find file : "%s"' %bpm_file)
else:
    iraf.imarith(operand1='gain.fits',
                  operand2=bpm_file,
                  op='*',
                  result='gain_bpm.fits',
                  verbose='yes'
                  )
    