#!/usr/bin/env python
################################################################################
#
# PAPI (PANIC PIpeline)
#
# mksuperflat.py
#
# Compute a super sky flat using the dither frames (IRAF implementation)
#
#
# 13/03/2009  : Created
# 15/04/2009  : created function and modified to accept command line arguments
#
################################################################################
#
# usage: mksuperflat.py dither_files.list superflat.fits
# 
# example:
# 
# input:
#
# output:
#
# method:
#
# notes:

import os
import pyraf
import iraf
import sys

# Pyraf modules
from pyraf import iraf
from iraf import noao
from iraf import imred
from iraf import ccdred

# Interact with FITS files
import pyfits

def makesuperflat(filelist, output_file):
  
    if os.path.exists(output_file): os.remove(output_file)
    
    iraf.imcombine(input=("'"+"@"+filelist+"'").replace('//','/'),
                output=output_file,
                combine='median',
                offset='none',
                reject='sigclip',
                lsigma=2.5,
                hsigma=2.5,
                scale='median',
                zero='none'
                #masktype='none'
                #scale='exposure',
                #expname='EXPTIME'
                #ParList = _getparlistname ('flatcombine')
               )
    
    median = float(iraf.imstat (
            images=("'"+output_file+"[100:1900,100:1900]'").replace('//','/'),
            fields='midpt',format='no',Stdout=1)[0])
                                              
    """iraf.imarith(operand1 = 'sflat.fits',
                    operand2 = median,
                    op = '/',
                    result ='sflatn.fits',
                    verbose = 'yes'
                    )                                              
    """
    flatframe = pyfits.open(output_file,'update')
    #Add a new keyword-->DATAMODE
    flatframe[0].header.update('DATAMODE',median,'Data mode of the frame')
    flatframe.close(output_verify='ignore') # This ignore any FITS standar violation and allow write/update the FITS file
                  
################################################################################            
#  Main
################################################################################
def usage():
    """Print help """
    print "\n\nUsage: \n\t>%s sci_files.list sflat.fits" %(sys.argv[0])
    
if __name__ == "__main__": 
    #log.debug("Starting the .....")
    
    if len(sys.argv)==3:
        try:
            makesuperflat(sys.argv[1], sys.argv[2])
            print 'Done !'
        except:
            print("Error while running ", sys.argv[0])
    else:
        usage()
        
        
        