#!/usr/bin/env python

################################################################################
#
#
# applyFlat.py
#
# Last update 11/03/2009
#
################################################################################

import os
import pyraf
import iraf

# Pyraf modules
from pyraf import iraf
from iraf import noao
from iraf import imred
from iraf import ccdred

# Interact with FITS files
import pyfits




source= open( "files.nip", "r" )
destination= open( "flatened_files.nip", "w" )
for line in source:
    line=line.replace(".fits", "_flatened.fits")
    destination.write( line )
source.close()
destination.close()

iraf.imarith(operand1="@files.nip",
            operand2='flat.fits',
            op='/',
            result="@flatened_files.nip",
            verbose='yes'
            )
            

  
  
mean = float(iraf.imstat (
            images='flat.fits[100:900,100:900]',
            fields='mean',format='no',Stdout=1)[0])
