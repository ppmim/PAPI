#!/usr/bin/env python
#
# 2007 Jul 8  - Andrew W Stephens - alpha version
# 2007 Jul 9  - AWS - beta version
# 2007 Jul 10 - AWS - move most operations to cleanquad function
# 2007 Jul 11 - AWS - use stddev to decide whether to keep orig quad
# 2007 Jul 14 - AWS - generalized code to allow different pattern sizes
# 2007 Jul 18 - AWS - fix bug generating index arrays
# 2007 Jul 20 - AWS - add quadrant bias-level normalization
# 2007 Jul 23 - AWS - add force option
# 2007 Aug 06 - AWS - f/6 spectroscopy mode: use top & bottom for pattern
# 2007 Aug 22 - AWS - add -a flag to use all pixels
# 2008 Jan 11 - AWS - check for available image extensions
# 2008 Feb 05 - AWS - don't close input file until the end (req'd if next>2)

# To Do:

# automate quadrant bias value normalization
#
# specification of image section affected by pattern
#
# specification of image section to use to calculate pattern
#  (i.e. to ignore spectra when calculating pattern)
#
# look at stddev of each row to identify which have pattern noise

#-----------------------------------------------------------------------

import datetime
import getopt
import numpy
import os
import pyfits
import string
import sys

#-----------------------------------------------------------------------

def usage():
    print ''
    print 'NAME'
    print '       nirinoise.py - filter pattern noise out of NIRI frames\n'
    print 'SYNOPSIS'
    print '       nirinoise.py [options] infile\n'
    print 'DESCRIPTION'
    print '       This script assumes that the NIRI pattern noise in any quadrant'
    print '       can be represented by a fixed pattern which is repeated over the'
    print '       entire quadrant.  The default size for this pattern is 16 pixels'
    print '       wide and 4 pixels high, which may be changed via the -x and -y'
    print '       flags.  The pattern is determined for each quadrant by taking the'
    print '       median of the pixel value distribution at each position in the'
    print '       pattern.  Once the median pattern has been determined for a'
    print '       quadrant it is replicated to cover the entire quadrant and then'
    print '       subtracted from the quadrant, and the mean of the pattern is'
    print '       added back to preserve flux.  The standard deviation of all the'
    print '       pixels in the quadrant is compared to that before the pattern'
    print '       subtraction, and if no reduction was achieved the subtraction is'
    print '       undone.  The pattern subtraction may be forced via the -f flag.'
    print '       This process is repeated for all four quadrants of the frame, and'
    print '       the cleaned frame is written to <infile>.clean.fits.  The pattern'
    print '       derived for each quadrant is saved to <infile>.pattern.fits,'
    print '       where the output pattern has a size of 2*patternx by 2*patterny'
    print '       pixels; if no pattern was subtracted from a quadrant (because the'
    print '       standard deviation did not decrease) that quadrant of the pattern'
    print '       file will be filled with zeros.  Occasionally there will be an'
    print '       offset in the bias values between the four quadrants and one may'
    print '       want to use the -b flag to remove this offset.  In this case the'
    print '       median value of each quadrant is calculated in the central 10% of'
    print '       the image (changeable via the -c flag) and subtracted, and the'
    print '       median value of the original image is added back to each'
    print '       quadrant.\n'
    print '       Things are more difficult for spectroscopy since the vertical sky'
    print '       lines complicate the pattern mesurement.  By default f/6'
    print '       spectroscopy with the 2-pixel or blue slits (which does not fill'
    print '       the detector), uses the empty regions at the bottom (1-272) and'
    print '       top (720-1024) of the array for calculating the pattern.  This is'
    print '       not possible for other modes of spectroscopy as the spectrum'
    print '       fills the detector.  For these modes it is best to do sky'
    print '       subtraction first (or adjacent pair subtraction) and then run'
    print '       this pattern-removal script on the difference image.  If you want'
    print '       to use this method for spectra using the f/6 2-pixel or blue sits,'
    print '       the -a flag will force using all pixels for the pattern determination.\n'
    print 'OPTIONS'
    print '       -a : use all pixels for pattern determination'
    print '       -b : adjust bias values'
    print '       -c <frac> : use central <frac> of image for bias adjustment [0.1]'
    print '       -f : force cleaning of all quads even if stddev does not decrease'
    print '       -o <file> : write output to <file> [infile.clean.fits]'
    print '       -p <file> : write pattern to <file> [infile.pattern.fits]'
    print '       -x <size> : set pattern x size in pix [16]'
    print '       -y <size> : set pattern y size in pix [4]'
    print '       -v : verbose debugging output\n'
    print 'VERSION'
    print '       2008 February 5'
    print ''
    raise SystemExit

#-----------------------------------------------------------------------

try:
    opts, args = getopt.getopt(sys.argv[1:], 'abc:fo:p:x:y:v')
except getopt.GetoptError:
    usage()
    sys.exit(2)

nargs = len(sys.argv[1:])
nopts = len(opts)

allpixels = False
biasadjust = False
frac = 0.1
force = False
outputfile = 'default'
patternfile = 'default'
pxsize = 16
pysize = 4
verbose = False

for o, a in opts:
    if o in ('-a'):
        allpixels = True # force using all pixels for pattern determination
    if o in ('-b'):
        biasadjust = True
    if o in ('-c'):      # use central fraction for bias normalization
        frac = float(a)
        nopts += 1
    if o in ('-f'):      # force pattern subtraction in every quadrant
        force = True
    if o in ('-o'):      # specify cleaned output file
        outputfile = a
        nopts += 1
    if o in ('-p'):      # specify pattern file
        patternfile = a
        nopts += 1
    if o in ('-x'):      # specify pattern x-dimension
        pxsize = int(a)
        nopts += 1
    if o in ('-y'):      # specify pattern y-dimension
        pysize = int(a)
        nopts += 1
    if o in ('-v'):      # verbose debugging output
        verbose = True

if (nargs - nopts) != 1:
    usage()

inputfile  = sys.argv[nopts+1]

if (outputfile == 'default'):
    outputfile = string.rstrip(inputfile,'fits') + 'clean.fits'

if (patternfile == 'default'):
    patternfile = string.rstrip(inputfile,'fits') + 'pattern.fits'


#-----------------------------------------------------------------------
# Error checking

if not os.path.exists(inputfile):      # check whether input file exists
    print inputfile, 'does not exist'
    sys.exit(2)

if os.path.exists(outputfile):        # check whether output file exists
    print '...removing old', outputfile
    os.remove(outputfile)

if os.path.exists(patternfile):      # check whether pattern file exists
    print '...removing old', patternfile
    os.remove(patternfile)

if (frac < 0.1):
    print 'ERROR: central fraction must be > 0.1'
    sys.exit(2)

if (frac > 1):
    print 'ERROR: central fraction must be < 1'
    sys.exit(2)

#-----------------------------------------------------------------------

def cleanquad(quad,patternin):
    global qsize                                         # quadrant size
    global pxsize, pysize                                 # pattern size

    if (verbose):
        print '...mean of input quadrant =', quad.mean()
        print '...median of input quadrant =', numpy.median(quad.flatten())

    # create arrays of indices which correspond to the pattern tiled over
    # the region of the input quadrant to be used for pattern determination
    inpx = len(patternin[0])
    inpy = len(patternin)
    if (verbose):
        print '...size of pattern determination region =',inpx,'x',inpy
    indx = numpy.tile(numpy.arange(0,inpx,pxsize), inpy/pysize)
    indy = numpy.arange(0,inpy,pysize).repeat(inpx/pxsize)
    if (verbose):
        print '...indx:', indx
        print '...indy:', indy

    # create blank pattern array:
    pattern = numpy.zeros(pysize*pxsize).reshape(pysize,pxsize)
    origstddev = quad.std()
    print '...standard deviation of input quadrant =', origstddev

    for ix in range(0, pxsize):    # find median pattern across quadrant
        for iy in range(0, pysize):
            pattern[iy,ix] = numpy.median(patternin[indy+iy,indx+ix])

    # tile pattern over quadrant:
    quadpattern = numpy.tile(pattern, (qsize/pysize, qsize/pxsize))
    clean = quad - quadpattern                        # subtract pattern
    cleanstddev = clean.std()             # calculate standard deviation
    print '...standard deviation of clean quadrant =', cleanstddev

    if (force):
        print '...forcing pattern subtraction'
    else:
        # has subtracting the pattern reduced the standard deviation?
        if ( origstddev < cleanstddev ):
            print '...no improvement; using the original quadrant'
            clean = quad      # reset quadrant pixels to original values
            pattern = pattern * 0                 # set pattern to zeros
        else:
            print '...improvement!'

    patternmean = pattern.mean()
    print '...adding back pattern mean:', patternmean
    clean += patternmean

    print '----------'
    return clean, pattern

#-----------------------------------------------------------------------

print '...reading', inputfile
hdulist = pyfits.open(inputfile)
next = len(hdulist)
print '...number of extensions = ', next

if ( next == 1 ):
   sci = 0
else:
   sci = 1
print '...assuming the science data are in extension', sci

try:
    naxis1,naxis2 = hdulist[sci].header['naxis1'], hdulist[sci].header['naxis2']
except:
    print 'ERROR: cannot get the dimensions of extension ', sci
    pyfits.info(inputfile)
    sys.exit(2)

print '...image dimensions = ', naxis1, 'x', naxis2
if (naxis1 != naxis2):
    print 'ERROR: NAXIS1 != NAXIS2'
    sys.exit(2)

inputmedian = numpy.median( hdulist[sci].data.flatten() )
print '...median of input image =', inputmedian
print '...mean of input image =', hdulist[sci].data.mean()

print '...pattern size =', pxsize, 'x', pysize

qsize = naxis1 / 2                                       # quadrant size


#-----------------------------------------------------------------------
# set regions to be used for pattern determination:

try:
    fpmask = hdulist[0].header['FPMASK']
except:
    print 'ERROR: cannot find FPMASK header keyword'
    print '   Assuming that this is imaging...'
    fpmask = 'f6-cam_G5208'

if (verbose):
    print '...FPMASK = ', fpmask

if     allpixels or \
       fpmask == 'f6-cam_G5208'  or \
       fpmask == 'f14-cam_G5209' or \
       fpmask == 'f32-cam_G5210':
    print '...using whole image for pattern determination'
    patternin1 = hdulist[sci].data[    0:qsize,      0:qsize ]
    patternin2 = hdulist[sci].data[    0:qsize,  qsize:naxis1]
    patternin3 = hdulist[sci].data[qsize:naxis1,     0:qsize ]
    patternin4 = hdulist[sci].data[qsize:naxis1, qsize:naxis1]

elif     fpmask == 'f6-2pixBl_G5214' or \
         fpmask == 'f6-4pixBl_G5215' or \
         fpmask == 'f6-6pixBl_G5216' or \
         fpmask == 'f6-2pix_G5211':
    print '...using region above and below slit for pattern determination'
    patternin1 = hdulist[sci].data[   0:272,        0:qsize ]
    patternin2 = hdulist[sci].data[   0:272,    qsize:naxis1]
    patternin3 = hdulist[sci].data[ 720:naxis1,     0:qsize ]
    patternin4 = hdulist[sci].data[ 720:naxis1, qsize:naxis1]

elif     fpmask == 'f6-4pix_G5212'  or \
         fpmask == 'f6-6pix_G5213'  or \
         fpmask == 'f32-6pix_G5229' or \
         fpmask == 'f32-9pix_G5230':
    print '...using whole image for pattern determination'
    print 'WARNING: Sky lines may be altered by pattern removal!'
    patternin1 = hdulist[sci].data[    0:qsize,      0:qsize ]
    patternin2 = hdulist[sci].data[    0:qsize,  qsize:naxis1]
    patternin3 = hdulist[sci].data[qsize:naxis1,     0:qsize ]
    patternin4 = hdulist[sci].data[qsize:naxis1, qsize:naxis1]

#-----------------------------------------------------------------------
# check that supplied pattern size is acceptable:
if pxsize > qsize or pysize > qsize:
    print 'ERROR: pattern size is too large'
    sys.exit(2)
if qsize%pxsize != 0 or qsize%pysize != 0:
    print 'ERROR: image size is not a multiple of the pattern size'
    sys.exit(2)

#-----------------------------------------------------------------------
# set input regions to be pattern-subtracted:

quad1 = hdulist[sci].data[    0:qsize,      0:qsize ]
quad2 = hdulist[sci].data[    0:qsize,  qsize:naxis1]
quad3 = hdulist[sci].data[qsize:naxis1,     0:qsize ]
quad4 = hdulist[sci].data[qsize:naxis1, qsize:naxis1]

#-----------------------------------------------------------------------

print '----------'
print 'Lower left quad:'
clean1, pattern1 = cleanquad(quad1,patternin1)
print 'Lower right quad:'
clean2, pattern2 = cleanquad(quad2,patternin2)
print 'Upper left quad:'
clean3, pattern3 = cleanquad(quad3,patternin3)
print 'Upper right quad:'
clean4, pattern4 = cleanquad(quad4,patternin4)

#-----------------------------------------------------------------------
print '...updating header...'
timestamp = datetime.datetime.now().strftime("%Y.%m.%d %H:%M:%S")
if (verbose):
    print '...time stamp =', timestamp
hdulist[0].header.add_history('Cleaned with nirinoise.py ' + timestamp)

print '...reassembling image...'
hdulist[sci].data[     0:qsize,      0:qsize] = clean1
hdulist[sci].data[     0:qsize, qsize:naxis1] = clean2
hdulist[sci].data[qsize:naxis1,      0:qsize] = clean3
hdulist[sci].data[qsize:naxis1, qsize:naxis1] = clean4

#-----------------------------------------------------------------------
# Normalize each quadrant using the central quarter of the image:

if (biasadjust):
    print '----------'
    print '...normalizing by the central', frac,' of the image...'
    # set the median of each quadrant to zero:
    print '...lower left:  %8.4f' %                  numpy.median(hdulist[sci].data[(1-frac)*qsize:qsize, (1-frac)*qsize:qsize].flatten())
    hdulist[sci].data[     0:qsize,      0:qsize] -= numpy.median(hdulist[sci].data[(1-frac)*qsize:qsize, (1-frac)*qsize:qsize].flatten())
    print '...upper left:  %8.4f' %                  numpy.median(hdulist[sci].data[(1-frac)*qsize:qsize, qsize:(1+frac)*qsize].flatten())
    hdulist[sci].data[     0:qsize, qsize:naxis1] -= numpy.median(hdulist[sci].data[(1-frac)*qsize:qsize, qsize:(1+frac)*qsize].flatten())
    print '...lower right: %8.4f' %                  numpy.median(hdulist[sci].data[qsize:(1+frac)*qsize, (1-frac)*qsize:qsize].flatten())
    hdulist[sci].data[qsize:naxis1,      0:qsize] -= numpy.median(hdulist[sci].data[qsize:(1+frac)*qsize, (1-frac)*qsize:qsize].flatten())
    print '...upper right: %8.4f' %                  numpy.median(hdulist[sci].data[qsize:(1+frac)*qsize, qsize:(1+frac)*qsize].flatten())
    hdulist[sci].data[qsize:naxis1, qsize:naxis1] -= numpy.median(hdulist[sci].data[qsize:(1+frac)*qsize, qsize:(1+frac)*qsize].flatten())
    # add back the median value of the input image:
    hdulist[sci].data += inputmedian
    print '----------'

if (verbose):
    print '...median of output image =', numpy.median( hdulist[sci].data.flatten() )
    print '...mean of output image =', hdulist[sci].data.mean()

#-----------------------------------------------------------------------

print '...writing', outputfile
hdulist.writeto(outputfile)

print '...assembling pattern image...'
# create blank array to hold pattern from each quadrant:
fullpattern=numpy.zeros(4*pxsize*pysize).reshape(2*pysize,2*pxsize)
# assemble the quadrants into a full pattern image:
fullpattern[    0:pysize,        0:pxsize]     = pattern1
fullpattern[    0:pysize,    pxsize:2*pxsize]  = pattern2
fullpattern[ pysize:2*pysize,    0:pxsize]     = pattern3
fullpattern[ pysize:2*pysize, pxsize:2*pxsize] = pattern4

print '...writing', patternfile
hdu = pyfits.PrimaryHDU(fullpattern)
hdu.writeto(patternfile)

hdulist.close()

#-----------------------------------------------------------------------
