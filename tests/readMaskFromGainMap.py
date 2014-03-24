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



def readMask( file ):
  
    f=pyfits.open(file)
    mask=numpy.where(f[0].data<0.0003, 0, 1)
    #mask=numpy.where( f[0].data<0.55, 0, numpy.where(f[0].data>1.44, 0, f[0].data) )
    #mask[(f[0].data < 0.55) | f[0].data> 1.45 ]=1
    
    os.remove("/tmp/mask.fits")
    hdu = pyfits.PrimaryHDU()
    hdu.data=mask     
    hdu.scale('int16') # importat to set first data type
    hdulist = pyfits.HDUList([hdu])
    hdu.header.update('OBJECT','MASTER_PIXEL_MASK')
    hdulist.writeto("/tmp/mask.fits")
    hdulist.close(output_verify='ignore')
    
    f.close()
    
################################################################################
#########   MAIN   #############################################################  
################################################################################

if __name__ == "__main__":
    readMask(sys.argv[1])
