#!/opt/local/bin python
""" 
Filename : nirc2dewarp.py

Function : Apply distortion correction to NIRC2 Narrow camera
	   individual frames (e.g., n0123.fits) using distortion 
	   solution from Yelda et al. 2010 and IRAF drizzle function
           (see: http://www.astro.ucla.edu/~ghezgroup/distortion/)

Files : Distortion x and y fits map (set path in distXgeoim and distYgeoim)
	nirc2_X_distortion.fits	
	nirc2_Y_distortion.fits

Inputs : 
	in_frame = (string) The name of a single FITS frame to undistort
	out_frame = (string) The name of a single output FITS frame
	wgt_frame = (string) The name of a single output weight FITS frame
	log_file = (string) The name of the log file from drizzle (stdout)

Example :
	> import nirc2dewarp
	> nirc2dewarp.undistort('M101.fits', 'M101_dr.fits', 'M101_dr_wgt.fits', 'M101_drizzle.log')


Notes : (1) If output file (e.g., outfile.fits) already exists and overwrite=False, 
		IRAF drizzle will "drizzle" on top of this existing file! The
		default behavior of the code is to remove all output files and recreate.
	(2) This does not take account for Differential Atmospheric Refraction (DAR)
	(3) We recommend that the input file have bad pixels already corrected, 
		as a pixel in the raw image may be drizzled onto multiple pixels in 
		the drizzled frame.
	(4) If you need to run from command line or IDL, add sys.argv and add 
		variables to code and remove functions
	(5) If used, please cite Yelda et al. 2010,
		"Improving Galactic Center Astrometry by Reducing 
			the Effects of Geometric Distortion"
	 
"""

import os, sys
import pyfits
from pyraf import iraf as ir

distCoef = ''
distXgeoim = 'nirc2_X_distortion.fits'
distYgeoim = 'nirc2_Y_distortion.fits'

def undistort(in_frame, out_frame, wgt_frame, log_file, overwrite=True):

    # Check if the log file exists, if so, delete
    rmall([log_file])

    # Check if the output files already exist
    if overwrite == True:
        rmall([out_frame, wgt_frame])

    hdr = pyfits.getheader(in_frame)
    imgsizeX = float(hdr['NAXIS1'])
    imgsizeY = float(hdr['NAXIS2'])
    if (imgsizeX >= imgsizeY):
        imgsize = imgsizeX
    else:
        imgsize = imgsizeY

    ### Setup parameters for drizzle
    setup_drizzle(imgsize)

    ### Drizzle individual file ###
    print 'Drizzling input file %s' % in_frame
    ir.drizzle(in_frame, out_frame, outweig=wgt_frame, Stdout=log_file)


def setup_drizzle(imgsize):
    """Setup drizzle parameters for NIRC2 data.
    @param imgsize: The size (in pixels) of the final drizzle image.
    This assumes that the image will be square.
    @type imgsize: int
    @param type: str
    """
    print 'Setting up drizzle'
    # Setup the drizzle parameters we will use
    ir.module.load('stsdas', doprint=0, hush=1)
    ir.module.load('analysis', doprint=0, hush=1)
    ir.module.load('dither', doprint=0, hush=1)
    ir.unlearn('drizzle')
    ir.drizzle.outweig = ''
    ir.drizzle.in_mask = ''
    ir.drizzle.wt_scl = 1
    ir.drizzle.outnx = imgsize
    ir.drizzle.outny = imgsize
    ir.drizzle.pixfrac = 1
    ir.drizzle.kernel = 'lanczos3'
    ir.drizzle.scale = 1
    ir.drizzle.coeffs = distCoef
    ir.drizzle.xgeoim = distXgeoim
    ir.drizzle.ygeoim = distYgeoim
    ir.drizzle.shft_un = 'input'
    ir.drizzle.shft_fr = 'output'
    ir.drizzle.align = 'center'
    ir.drizzle.expkey = 'ITIME'
    ir.drizzle.in_un = 'counts'
    ir.drizzle.out_un = 'counts'


def rmall(files):
    """Remove list of files without confirmation."""
    for file in files:
        if os.access(file, os.F_OK): os.remove(file)
