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
# 2008 Oct 02 - AWS - don't write pattern unless given the -p flag
# 2009 May 03 - AWS - use conformant default output file name
# 2009 May 13 - AWS - verify FITS header (old images have unquoted release date)
# 2009 May 22 - AWS - output full-frame pattern
# 2009 May 22 - AWS - improve quadrant bias normalization
# 2009 May 23 - AWS - add optional sky frame

# To Do:
# check if run before
# properly handle images < 1024 x 1024
# specification of image section to use to calculate pattern
#  (i.e. to ignore spectra when calculating pattern)
# look at stddev of each row to identify which have pattern noise
# split each (full-frame) quadrant into 4 bands

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
    print '       This script assumes that the NIRI pattern noise in a quadrant'
    print '       can be represented by a fixed pattern which is repeated over the'
    print '       entire quadrant.  The default size for this pattern is 16 pixels'
    print '       wide and 4 pixels high (which may be changed via the -x and -y'
    print '       flags).  The pattern is determined for each quadrant by taking the'
    print '       median of the pixel value distribution at each position in the'
    print '       pattern.  Once the median pattern has been determined for a'
    print '       quadrant it is replicated to cover the entire quadrant and'
    print '       subtracted, and the mean of the pattern is added back to preserve'
    print '       flux.  The standard deviation of all the pixels in the quadrant'
    print '       is compared to that before the pattern subtraction, and if no'
    print '       reduction was achieved the subtraction is undone.  The pattern'
    print '       subtraction may be forced via the -f flag.  This process is'
    print '       repeated for all four quadrants and the cleaned frame is written'
    print '       to c<infile> (or the file specified with the -o flag).  The'
    print '       pattern derived for each quadrant may be saved with the -p flag.'
    print ' '
    print '       Pattern noise is often accompanied by an offset in the bias'
    print '       values between the four quadrants.  One may want to use the'
    print '       -b flag to try to remove this offset.  This attempts to match'
    print '       the iteratively determined median value of each quadrant.'
    print '       This method works best with sky subtraction (i.e. with the -s'
    print '       flag), and does not work well if there are large extended objects'
    print '       in the frame.  By default the median is determined from the'
    print '       entire frame, although the -c flag will only use a central'
    print '       portion of the image.  Note that the derived quadrant offsets'
    print '       will be applied to the output pattern file.'
    print ' '
    print '       Removing the pattern from spectroscopy is more difficult because'
    print '       of many vertical sky lines.  By default f/6 spectroscopy with the'
    print '       2-pixel or blue slits (which do not fill the detector), uses the'
    print '       empty regions at the bottom (1-272) and top (720-1024) of the'
    print '       array for measuring the pattern.  This is not possible for other'
    print '       modes of spectroscopy where the spectrum fills the detector.'
    print '       For these modes it is best to do sky subtraction before pattern'
    print '       removal.  The quickest method is to pass a sky frame (or an offset'
    print '       frame) via the -s flag.  The manual method is to generate and'
    print '       subtract the sky, determine and save the pattern via the -p flag,'
    print '       then subtract the pattern from the original image.  One may use'
    print '       the -a flag to force using all of the pixels for the pattern'
    print '       determination.\n'
    print 'OPTIONS'
    print '       -a : use all pixels for pattern determination'
    print '       -b : adjust bias values for each quadrant'
    print '       -c <frac> : use central <fraction> of image for bias adjustment [1]'
    print '       -f : force cleaning of all quads even if stddev does not decrease'
    print '       -o <file> : write output to <file> (instead of c<infile>)'
    print '       -p <file> : write full-frame pattern to <file>'
    print '       -s <sky> : sky frame to help in pattern recognition'
    print '       -v : verbose debugging output'
    print '       -x <size> : set pattern x size in pix [16]'
    print '       -y <size> : set pattern y size in pix [4]\n'
    print 'VERSION'
    print '       2009 May 25'
    print ''
    raise SystemExit

#-----------------------------------------------------------------------

try:
    opts, args = getopt.getopt(sys.argv[1:], 'abc:fo:p:s:x:y:v')
except getopt.GetoptError:
    usage()
    sys.exit(2)

nargs = len(sys.argv[1:])
nopts = len(opts)

allpixels = False
biasadjust = False
cfrac = 1.0              # use whole image
force = False
outputfile = 'default'
savepattern = False
skysub = False
pxsize = 16
pysize = 4
verbose = False

for o, a in opts:
    if o in ('-a'):      # force using all pixels for pattern determination
        allpixels = True
    if o in ('-b'):      # try to adjust quadrant bias values
        biasadjust = True
    if o in ('-c'):      # use central fraction for bias normalization
        cfrac = float(a)
        nopts += 1
    if o in ('-f'):      # force pattern subtraction in every quadrant
        force = True
    if o in ('-o'):      # specify cleaned output file
        outputfile = a
        nopts += 1
    if o in ('-p'):      # write pattern file
        patternfile = a
        savepattern = True
        nopts += 1
    if o in ('-s'):      # sky frame
        skyfile = a
        skysub = True
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
    outputfile = 'c' + os.path.basename(inputfile)
else:
    if ( not outputfile.endswith('.fits') ):
        outputfile = outputfile + '.fits'

if (savepattern):
    if ( not patternfile.endswith('.fits') ):
        patternfile = patternfile + '.fits'

if (skysub):
     if ( not skyfile.endswith('.fits') ):
        skyfile = skyfile + '.fits'

#-----------------------------------------------------------------------
# Error checking

if not os.path.exists(inputfile):      # check whether input file exists
    print inputfile, 'does not exist'
    sys.exit(2)

if os.path.exists(outputfile):        # check whether output file exists
    print '...removing old', outputfile
    os.remove(outputfile)

if (savepattern):
    if os.path.exists(patternfile):  # check whether pattern file exists
        print '...removing old', patternfile
        os.remove(patternfile)

if (skysub):
    if not os.path.exists(skyfile):      # check whether sky file exists
        print skyfile, 'does not exist'
        sys.exit(2)

if (cfrac < 0.1):
    print 'ERROR: central fraction must be >= 0.1'
    sys.exit(2)

if (cfrac > 1):
    print 'ERROR: central fraction must be <= 1'
    sys.exit(2)

#-----------------------------------------------------------------------

print '...reading', inputfile
hdulist = pyfits.open(inputfile)

print '...verifying...'
if (verbose):
    hdulist.verify('fix')
else:
    hdulist.verify('silentfix')

next = len(hdulist)
if (verbose):
    print '...number of extensions = ', next

if ( next == 1 ):
   sci = 0
else:
   sci = 1
print '...assuming the science data are in extension', sci

image = numpy.array(hdulist[sci].data, dtype=numpy.double)

try:
    instrument = hdulist[0].header['INSTRUME']
    if (verbose):
        print '...instrument =', instrument
except:
    print 'WARNING: cannot determine instrument'

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

qsize = naxis1 / 2                                       # quadrant size

print '...pattern size =', pxsize, 'x', pysize

#-----------------------------------------------------------------------

if (skysub):
    print '...reading sky', skyfile
    sky = pyfits.open(skyfile)
    print '...verifying sky...'
    if (verbose):
        sky.verify('fix')
    else:
        sky.verify('silentfix')
    skyimage = numpy.array(sky[sci].data, dtype=numpy.double)

    # NEED ERROR CHECKING HERE! (extensions, image size, filter, instrument, etc.)

#-----------------------------------------------------------------------

def cleanquad(quad,patternin):
    # quad = region to be pattern-subtracted
    # patternin = region to use for pattern determination
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
    print '...standard deviation of input quadrant =%9.3f' % origstddev

    # find median pattern across quadrant:
    for ix in range(0, pxsize):    
        for iy in range(0, pysize):
            pattern[iy,ix] = numpy.median(patternin[indy+iy,indx+ix])

    # tile pattern over quadrant:
    quadpattern = numpy.tile(pattern, (qsize/pysize, qsize/pxsize))
    quadpattern -= pattern.mean()           # set the mean value to zero
    clean = quad - quadpattern                        # subtract pattern
    cleanstddev = clean.std()             # calculate standard deviation
    print '...standard deviation of clean quadrant =%9.3f' % cleanstddev

    if (force):
        print '...forcing pattern subtraction'
    else:
        # has subtracting the pattern reduced the standard deviation?
        if ( origstddev < cleanstddev ):
            print '...no improvement; using the original quadrant'
            clean = quad      # reset quadrant pixels to original values
            quadpattern = quadpattern * 0         # set pattern to zeros
        else:
            print '...improvement!'

    return clean, quadpattern

#-----------------------------------------------------------------------
# set regions to be used for pattern determination:

try:
    fpmask = hdulist[0].header['FPMASK']
except:
    print 'WARNING: cannot find FPMASK header keyword'
    print '   Assuming that this is imaging...'
    fpmask = 'f6-cam_G5208'

if (verbose):
    print '...FPMASK = ', fpmask

if     allpixels or \
       fpmask == 'f6-cam_G5208'  or \
       fpmask == 'f14-cam_G5209' or \
       fpmask == 'f32-cam_G5210':
    print '...using whole image for pattern determination'
    w1,x1,y1,z1 = 0,qsize,          0,qsize
    w2,x2,y2,z2 = 0,qsize,      qsize,naxis1
    w3,x3,y3,z3 = qsize,naxis1,     0,qsize
    w4,x4,y4,z4 = qsize,naxis1, qsize,naxis1
    lsigma = 3.0
    hsigma = 1.0      # set a very small upper threshold to reject stars

elif     fpmask == 'f6-2pixBl_G5214' or \
         fpmask == 'f6-4pixBl_G5215' or \
         fpmask == 'f6-6pixBl_G5216' or \
         fpmask == 'f6-2pix_G5211':
    print '...using region above and below slit for pattern determination'
    w1,x1,y1,z1 =   0,272,        0,qsize
    w2,x2,y2,z2 =   0,272,    qsize,naxis1
    w3,x3,y3,z3 = 720,naxis1,     0,qsize
    w4,x4,y4,z4 = 720,naxis1, qsize,naxis1
    lsigma = 3.0
    hsigma = 3.0

elif     fpmask == 'f6-4pix_G5212'  or \
         fpmask == 'f6-6pix_G5213'  or \
         fpmask == 'f32-6pix_G5229' or \
         fpmask == 'f32-9pix_G5230':
    print '...using whole image for pattern determination'
    print 'WARNING: Sky lines may be altered by pattern removal!'
    w1,x1,y1,z1 = 0,qsize,          0,qsize
    w2,x2,y2,z2 = 0,qsize,      qsize,naxis1
    w3,x3,y3,z3 = qsize,naxis1,     0,qsize
    w4,x4,y4,z4 = qsize,naxis1, qsize,naxis1
    lsigma = 3.0
    hsigma = 3.0

if (skysub):
    patternin1 = image[w1:x1, y1:z1] - skyimage[w1:x1, y1:z1]
    patternin2 = image[w2:x2, y2:z2] - skyimage[w2:x2, y2:z2]
    patternin3 = image[w3:x3, y3:z3] - skyimage[w3:x3, y3:z3]
    patternin4 = image[w4:x4, y4:z4] - skyimage[w4:x4, y4:z4]
else:
    patternin1 = image[w1:x1, y1:z1]
    patternin2 = image[w2:x2, y2:z2]
    patternin3 = image[w3:x3, y3:z3]
    patternin4 = image[w4:x4, y4:z4]

#-----------------------------------------------------------------------
# Calculate means and medians for reference:

inputmean = numpy.mean(image)                       # whole image
allpix = numpy.concatenate((patternin1,patternin2,patternin3,patternin4))
inputmedian = numpy.median(allpix)  # region used for pattern & offset determination
inputstddev = numpy.std(allpix)

if (verbose):
    print '...mean of input image =', inputmean
    print '...median of pixels used for pattern measurement =', inputmedian
    print '...stddev of pixels used for pattern measurement =', inputstddev

#-----------------------------------------------------------------------
# check that supplied pattern size is acceptable:
if pxsize > qsize or pysize > qsize:
    print 'ERROR: pattern size is too large'
    sys.exit(2)
if qsize%pxsize != 0 or qsize%pysize != 0:
    print 'ERROR: image size is not a multiple of the pattern size'
    sys.exit(2)

#-----------------------------------------------------------------------
# calculate and subtract pattern:

quad1 = image[    0:qsize,     0:qsize  ]
quad2 = image[    0:qsize,  qsize:naxis1]
quad3 = image[qsize:naxis1,     0:qsize ]
quad4 = image[qsize:naxis1, qsize:naxis1]

print '...lower left quadrant:'
clean1, pattern1 = cleanquad(quad1,patternin1)
print '...lower right quadrant:'
clean2, pattern2 = cleanquad(quad2,patternin2)
print '...upper left quadrant:'
clean3, pattern3 = cleanquad(quad3,patternin3)
print '...upper right quadrant:'
clean4, pattern4 = cleanquad(quad4,patternin4)

print '...reassembling image...'
image[     0:qsize,      0:qsize] = clean1
image[     0:qsize, qsize:naxis1] = clean2
image[qsize:naxis1,      0:qsize] = clean3
image[qsize:naxis1, qsize:naxis1] = clean4

print '...updating header...'
timestamp = datetime.datetime.now().strftime("%Y.%m.%d %H:%M:%S")
if (verbose):
    print '...time stamp =', timestamp
hdulist[0].header.add_history('Cleaned with nirinoise.py ' + timestamp)

#-----------------------------------------------------------------------
# Normalize each quadrant (with sigma clipping) to the input median:

def normquad(quad):
    global lsigma                             # low-rejection threshold
    global hsigma                             # high-rejection threshold
    if (verbose):
        print '...calculating median using low-sigma =',lsigma,' and high-sigma =',hsigma
    mincts = inputmedian - lsigma * inputstddev
    maxcts = inputmedian + hsigma * inputstddev
    if (verbose):
        print '...input median=',inputmedian,' min=',mincts,' max=',maxcts
   
    flatquad = quad.flatten()
    npix = numpy.size(flatquad)
    dn = 100
    while dn > 10:
        tmpquad = flatquad[(flatquad>mincts) & (flatquad<maxcts)]
        median = numpy.median(tmpquad)
        stddev = numpy.std(tmpquad)
        mincts = median - lsigma * stddev
        maxcts = median + hsigma * stddev
        dn = npix - numpy.size(tmpquad)
        npix = numpy.size(tmpquad)
        if (verbose):
            print '...median=',median,' stddev=',stddev,' min=',mincts,' max=',maxcts,' npix=',npix,' dn=',dn
    offset = inputmedian - median
    return offset

if (biasadjust):
    print '...normalizing the bias level of each quadrant...'
    # And apply the measured offset to the pattern output

    if (cfrac != 1):  # the user has specified a non-unity central fraction
        print '...using the central %4.0f pixels' % (2*cfrac*qsize)
        patternin1 = image[(1-cfrac)*qsize:qsize, (1-cfrac)*qsize:qsize]
        patternin2 = image[(1-cfrac)*qsize:qsize, qsize:(1+cfrac)*qsize]
        patternin3 = image[qsize:(1+cfrac)*qsize, (1-cfrac)*qsize:qsize]
        patternin4 = image[qsize:(1+cfrac)*qsize, qsize:(1+cfrac)*qsize]

    offset = normquad(patternin1)
    image[0:qsize, 0:qsize] += offset
    pattern1 -= offset
    print '...lower left quadrant offset =  %7.2f' % offset

    offset = normquad(patternin2)
    image[0:qsize, qsize:naxis1] += offset
    pattern2 -= offset
    print '...lower right quadrant offset = %7.2f' % offset

    offset = normquad(patternin3)
    image[qsize:naxis1, 0:qsize] += offset
    pattern3 -= offset
    print '...upper left quadrant offset =  %7.2f' % offset

    offset = normquad(patternin4)
    image[qsize:naxis1, qsize:naxis1] += offset
    pattern4 -= offset
    print '...upper right quadrant offset = %7.2f' % offset

    offset = inputmean - numpy.mean(image)
    print '...adjusting whole image by %7.2f to match input image...' % offset
    image += offset

if (verbose):
    print '...mean of input image =  %9.3f' % inputmean
    print '...mean of output image = %9.3f' % numpy.mean(image)
    print '...median of output image = %9.3f' % numpy.median(image)

#-----------------------------------------------------------------------
# Write cleaned output image

print '...writing', outputfile
hdulist[sci].data = image
hdulist.writeto(outputfile)

#-----------------------------------------------------------------------
# Write pattern image

if (savepattern):
    print '...assembling pattern image...'
    # create blank pattern array:
    fullpattern = numpy.zeros(naxis2*naxis1).reshape(naxis2,naxis1)
    # assemble the quadrants into a full pattern image:
    fullpattern[     0:qsize,      0:qsize] = pattern1
    fullpattern[     0:qsize, qsize:naxis1] = pattern2
    fullpattern[qsize:naxis1,      0:qsize] = pattern3
    fullpattern[qsize:naxis1, qsize:naxis1] = pattern4
    # normalize to zero:
    fullpattern -= fullpattern.mean()
    print '...writing', patternfile
    hdu = pyfits.PrimaryHDU(fullpattern)
    hdu.writeto(patternfile)

hdulist.close()

#-----------------------------------------------------------------------
