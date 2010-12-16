#!/usr/bin/env python

################################################################################
#
#
# imrep.py
#
# Last update 20/01/2009
#
################################################################################

import os
import sys

# Pyraf modules
from pyraf import iraf
from iraf import noao
from iraf import imred
from iraf import ccdred


# Check the mask file exists
file=sys.argv[1]

if not os.path.exists( file ):
    print('Error, No external file found: "%s"' %file)
else:
    values = iraf.imstat (
            images=file+"[400:1600,400:1600]",
            fields='mean,stddev',format='no',Stdout=1)[0].split()
            
    mean=float(values[0])
    sigma=float(values[1])
            
    print "MEAN=" ,mean
    print "SIGMA=" ,sigma
    iraf.imrepl(images=file,
                value=mean,
                lower='INDEF',
                upper=mean-4*sigma
               )

