#!/opt/local/bin/python2.6
###!/usr/bin/python2.6
# You may need to adapt the above path to call your preferred python executable
#
###########################################################################
#
# Copyright (c) 2009-2012 Max Planck Institute for Extraterrestrial Physics
# All Rights Reserved.
# JDK
#
# Do not distribute without contacting the author (kurk@mpe.mpg.de)
###########################################################################

import os, sys #, math, subprocess, shutil, time
#sys.path.insert(0,'/tera1/software/python-modules/lib/python2.6/site-packages')
sys.path.insert(0,'/opt2/local/lib/python2.6/site-packages')
import pyfits
import numpy

# first argument is fits file to be corrected
# read fits file to be corrected
# determine read mode
# read fits file containing polynomial for correction
# correct fits file
# save new fits file

# first argument is fits file to be corrected, e.g. ../raw_linearity/luci.20110407.0064.fits or ../raw_linearity/luci.20110407.0169.fits
if len(sys.argv) != 2:
    print sys.argv[0] + ' myfile.fits'
    exit(1)
else:
    myfitsname = sys.argv[1]

# read fits file to be corrected
print 'Reading '+myfitsname
myfits_hdu = pyfits.open(myfitsname)
myfits_hdr = myfits_hdu[0].header
# store scaling information for later use, as this is deleted when the array is updated
BZERO = myfits_hdr['BZERO']
BSCALE = myfits_hdr['BSCALE']
BITPIX = myfits_hdr['BITPIX']
bitpix_designation = pyfits.ImageHDU.NumCode[BITPIX]
# print BSCALE, BZERO, BITPIX, bitpix_designation
# read data array
myfits = myfits_hdu[0].data
myfits_hdu.close()

# determine read mode
if myfits_hdr['READMODE'] == 'multiple.endpoints':
    rmode = 'mer'
else:
    rmode = 'o2dcr'

# read fits file containing polynomial for correction
polyfile = 'nonlin_coeffs_'+rmode+'.fits'
print 'Reading '+polyfile
# for some weird reason I don't manage to read the cube using the standard functions
# but have to use the convenience function getdata
poly = pyfits.getdata(polyfile)
### read bad pixel file (currently not used and commented out)
# badfile = 'badpixel_'+rmode+'.fits'
# print 'Reading '+badfile
# bad = pyfits.getdata(badfile)

# correct fits file
# you may want to constrain this only to "positive" corrections
fac = 1 / (1 + poly[1,:,:] * myfits/1e4 + poly[0,:,:] * (myfits/1e4)**2)

# save new fits file
mfnp = myfitsname.partition('.fits')
# add _lincor before .fits extension, or at the end if no such extension present
outfitsname = mfnp[0] + '_lincor' + mfnp[1] + mfnp[2]
# keep same data type
myfits_hdu[0].data = (myfits * fac).astype(myfits.dtype)
if os.path.exists(outfitsname):
    os.unlink(outfitsname)
    print 'Overwriting '+outfitsname
else:
    print 'Writing '+outfitsname
# scale back data to original values
myfits_hdu[0].scale(bitpix_designation,'old')
myfits_hdu.writeto(outfitsname)



# in iraf
#imarith ../raw_linearity/luci.20110407.0064 / 1e4 myfile
#imfunc  myfile square square
#imarith nonlin_coeffs_o2dcr[*,*,1] * square fac2
#imarith nonlin_coeffs_o2dcr[*,*,2] * myfile fac1
#imarith fac1 + fac2 fac
#imarith fac + 1 fac
#imarith myfile / fac myfile_lincor
#imarith myfile_lincor * 1e4 myfile_lincor
