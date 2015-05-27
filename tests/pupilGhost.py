#! /usr/bin/env python
# WHIRC data reduction: Bad Pixels and Pupil Ghost
# https://kareemelbadry.wordpress.com/2014/06/03/whirc-data-reduction-bad-pixels-and-pupil-ghost/

import numpy as np
import pyfits
from astropy.convolution import convolve, Gaussian2DKernel

#READ IN FLATS
J_flat = pyfits.getdata('Calibs/J_master_flat.fits')

#READ IN BAD PIXEL MAP
bad_pix = pyfits.getdata('Calibs/bpix.whirc.fits')

#REMOVE BAD PIXELS IN FLATS TO BE SMOOTHED
for row in range(2048):
   for column in range(2048):
      if bad_pix[row][column] == 1:
           #bit sketchy, but close enough. 
           adjacent = [J_flat[row-1][column-1],J_flat[row-2][column-2],J_flat[row-3][column-3]]
           J_flat[row][column] = np.median(adjacent)

#MAKE KERNEL A GAUSSIAN WITH SIGMA = 10
gauss_kernel = Gaussian2DKernel(10)

#SMOOTH FLATS WITH 2D CONVOLUTION
print "smoothing J_flat with Gaussian2DKernel..."
smoothed_J = convolve(J_flat, gauss_kernel)

#NORMALIZE FLATS TO ONE
smoothed_J = smoothed_J/np.median(smoothed_J)
J_flat_bad = pyfits.getdata('Calibs/J_master_flat.fits')

#SUBTRACT OFF SMOOTHED SURFACE
#ADD ONE TO KEEP NORM OF FLATS AT ONE
J_flat = J_flat_bad - smoothed_J + 1

#WRITE TO FILE
J_file = pyfits.PrimaryHDU(J_flat)
J_file.writeto('Calibs/New_J_Flat.fits', clobber = True)
