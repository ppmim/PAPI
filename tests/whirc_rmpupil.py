# removes the pupil ghost from flats
# run in directory above Raw, Calibs, Final 
# run after having run whirc_sort() -- J_master_flat.fits 
# and K_master_flat.fits must be in same folder. 
# Must have 'bpix.whirc.fits' in Calibs/ 
# download from http://www.noao.edu/kpno/manuals/whirc/datared.html

import numpy as np
import pyfits
from astropy.convolution import convolve, Gaussian2DKernel

def whirc_rmpupil():
    #READ IN FLATS
    J_flat = pyfits.getdata('Calibs/J_master_flat.fits')
    K_flat = pyfits.getdata('Calibs/K_master_flat.fits')

    #READ IN BAD PIXEL MAP
    bad_pix = pyfits.getdata('Calibs/bpix.whirc.fits')
    
    #REMOVE BAD PIXELS IN FLATS TO BE SMOOTHED
    print "fixing bad pixels in J_flat..."
    
    for row in range(2048):
        for column in range(2048):
            if bad_pix[row][column] == 1:
                #bit sketchy, but close enough. 
                adjacent = [J_flat[row-1][column-1],J_flat[row-2][column-2],J_flat[row-3][column-3]]
                J_flat[row][column] = np.median(adjacent)
                
    #j_fixed = pyfits.PrimaryHDU(J_flat)
    #j_fixed.writeto('Calibs/j_good_pix.fits',clobber = True)
    
    print "fixing bad pixels in K_flat..."
    for row in range(2048):
        for column in range(2048):
            if bad_pix[row][column] == 1:
                adjacent = [K_flat[row-1][column-1],K_flat[row-2][column-2],K_flat[row-3][column-3]]
                K_flat[row][column] = np.median(adjacent)

    #MAKE KERNEL A GAUSSIAN WITH SIGMA = 10
    gauss_kernel = Gaussian2DKernel(10)

    #SMOOTH FLATS WITH 2D CONVOLUTION
    print "smoothing J_flat with Gaussian2DKernel..."
    smoothed_J = convolve(J_flat, gauss_kernel)
    print "smoothing K_flat..."
    smoothed_K = convolve(K_flat, gauss_kernel)
    
    smoothK = pyfits.PrimaryHDU(smoothed_K)
    smoothK.writeto('Calibs/smoothed_flat.fits',clobber = True)
    
    #NORMALIZE FLATS TO ONE
    print "normalizing and subtracting pupil ghost..."
        
    smoothed_J = smoothed_J/np.median(smoothed_J)
    smoothed_K = smoothed_K/np.median(smoothed_K)
    
    #READ IN AGAIN, NOT CORRECTED FOR BAD PIXELS
    J_flat_bad = pyfits.getdata('Calibs/J_master_flat.fits')
    K_flat_bad = pyfits.getdata('Calibs/K_master_flat.fits')
    
    J_flat_bad = J_flat_bad/np.median(J_flat_bad)
    K_flat_bad = K_flat_bad/np.median(K_flat_bad)
    
    #SUBTRACT OFF SMOOTHED SURFACE
    #ADD ONE TO KEEP NORM OF FLATS AT ONE
    J_flat = J_flat_bad - smoothed_J + 1
    K_flat = K_flat_bad - smoothed_K + 1
    
    #WRITE TO FILE
    print "writing to file..."
    
    J_file = pyfits.PrimaryHDU(J_flat)
    K_file = pyfits.PrimaryHDU(K_flat)
    
    J_file.writeto('Calibs/New_J_Flat.fits', clobber = True)
    K_file.writeto('Calibs/New_K_Flat.fits', clobber = True)