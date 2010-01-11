#!/usr/bin/python
"""
Run SExtractor on an image.

Demonstrates how to remotely run a task and retrieve the results.

"""
import time, urllib
from astrogrid import acr
from astrogrid import Applications, MySpace
from astrogrid.utils import mkURI

# We need to login.
# Insert the credentials here if not already in ~/.python-acr
acr.login()

# Image
image='http://www2.astrogrid.org/science/documentation/workbench-advanced/advanced-usage/scripting-user-s-guide/image.fits'

# We define some configuration files here

config="""
#-------------------------------- Catalog ------------------------------------

CATALOG_NAME	test.cat	# name of the output catalog
CATALOG_TYPE	ASCII_HEAD	# "NONE","ASCII_HEAD","ASCII","FITS_1.0"
				# "FITS_LDAC" or "ASCII_SKYCAT"

PARAMETERS_NAME	int-wfs.param	
                                # name of the file containing catalog contents

#------------------------------- Extraction ----------------------------------

DETECT_TYPE	CCD		# "CCD" or "PHOTO" (*)
#FLAG_IMAGE	flag    	# filename for an input FLAG-image
DETECT_MINAREA	3		# minimum number of pixels above threshold
DETECT_THRESH	1.0		# <sigmas> or <threshold>,<ZP> in mag.arcsec-2
ANALYSIS_THRESH	1.0		# <sigmas> or <threshold>,<ZP> in mag.arcsec-2

FILTER		Y		# apply filter for detection ("Y" or "N")?
FILTER_NAME	gauss_4.0_7x7.conv	
                                # name of the file containing the filter

DEBLEND_NTHRESH	32		# Number of deblending sub-thresholds
DEBLEND_MINCONT	0.005		# Minimum contrast parameter for deblending

CLEAN		Y		# Clean spurious detections? (Y or N)?
CLEAN_PARAM	1.0		# Cleaning efficiency

MASK_TYPE	CORRECT		# type of detection MASKing: can be one of
				# "NONE", "BLANK" or "CORRECT"

#------------------------------ Photometry -----------------------------------

PHOT_APERTURES	3.5, 7.0, 9.9, 14.0, 19.8	# MAG_APER aperture diameter(s) in pixels
PHOT_AUTOPARAMS	2.5, 3.5	# MAG_AUTO parameters: <Kron_fact>,<min_radius>

SATUR_LEVEL	50000.0		# level (in ADUs) at which arises saturation

MAG_ZEROPOINT	0.0		# magnitude zero-point
MAG_GAMMA	4.0		# gamma of emulsion (for photographic scans)
GAIN		1.0		# detector gain in e-/ADU.
PIXEL_SCALE	0.333		# size of pixel in arcsec (0=use FITS WCS info).

#------------------------- Star/Galaxy Separation ----------------------------

SEEING_FWHM	1.2		# stellar FWHM in arcsec
STARNNW_NAME	default.nnw	
                                # Neural-Network_Weight table filename

#------------------------------ Background -----------------------------------

BACK_SIZE	64		# Background mesh: <size> or <width>,<height>
BACK_FILTERSIZE	3		# Background filter: <size> or <width>,<height>

BACKPHOTO_TYPE	GLOBAL		# can be "GLOBAL" or "LOCAL" (*)
BACKPHOTO_THICK	24		# thickness of the background LOCAL annulus (*)

#------------------------------ Check Image ----------------------------------

CHECKIMAGE_TYPE	NONE		# can be one of "NONE", "BACKGROUND",
				# "MINIBACKGROUND", "-BACKGROUND", "OBJECTS",
				# "-OBJECTS", "SEGMENTATION", "APERTURES",
				# or "FILTERED" (*)
#CHECKIMAGE_NAME	check	# Filename for the check-image (*)

#--------------------- Memory (change with caution!) -------------------------

MEMORY_OBJSTACK	2000		# number of objects in stack
MEMORY_PIXSTACK	100000		# number of pixels in stack
MEMORY_BUFSIZE	1024		# number of lines in buffer

#----------------------------- Miscellaneous ---------------------------------

VERBOSE_TYPE	NORMAL		# can be "QUIET", "NORMAL" or "FULL" (*)
"""

params="""NUMBER
X_WORLD
Y_WORLD
X_IMAGE
Y_IMAGE
MAG_ISO
MAGERR_ISO
MAG_AUTO
MAGERR_AUTO
CLASS_STAR
"""

filter="""CONV NORM
# 3x3 ``all-ground'' convolution mask with FWHM = 2 pixels.
1 2 1
2 4 2
1 2 1
"""

# We write the image to MySpace
m = MySpace()
m.savefile(image, '#sextractor/image.fits')

# Application ID. This is the name of the application we are going to run.
# It can be obtained from a registry search.
id = 'ivo://org.astrogrid/SExtractor'

# We fill in the application parameters
app = Applications(id)

app.inputs['ANALYSIS_THRESH']['value']=1.5
app.inputs['IMAGE_BAND']['value']='R'
app.inputs['MAG_ZEROPOINT']['value']=25.0
app.inputs['SEEING_FWHM']['value']=1.2

app.inputs['PARAMETERS_NAME']['value']=params
app.inputs['FILTER_NAME']['value']=filter
app.inputs['config_file']['value']=config

app.inputs['DetectionImage']['value'] = mkURI('#sextractor/image.fits')
app.inputs['DetectionImage']['indirect'] = True
app.inputs['PhotoImage']['value'] = mkURI('#sextractor/image.fits')
app.inputs['PhotoImage']['indirect'] = True

app.outputs['CATALOG_NAME']['value'] = mkURI('#sextractor/image_cat.fits')
app.outputs['CATALOG_NAME']['indirect'] = True

# We submit the application to the server at Cambridge.
# This is a temporal workaround because the other servers are not working.
task=app.submit(server='ivo://uk.ac.cam.ast/CEACommandline')

# Wait until completed
time.sleep(10)
while task.status() <> 'COMPLETED':
	time.sleep(10)
	
print 'Catalogue written to\n%s' % task.results()

# Now with Aladin running we can send the image and catalogue
from astrogrid.utils import broadcast
broadcast('#sextractor/image.fits')
broadcast('#sextractor/image_cat.fits')

