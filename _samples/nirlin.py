#! /usr/bin/env python

# Copyright (c) 2011 IAA-CSIC All rights reserved.
# Institute of Astrophysics of Andalusia, IAA-CSIC
#
# This file is part of PAPI (PANIC Pipeline)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

################################################################################
#
# PAPI (PANIC PIpeline)
#
# nirlin.py
#
# Created    : 03/02/2011    jmiguel@iaa.es
# Last update: 03/02/2011    jmiguel@iaa.es
# TODO : - add support to MEF files for PANIC FITS
#        - read coeffs from a external text file
#        - 

################################################################################

################################################################################


#-----------------------------------------------------------------------

import datetime
import getopt
import numpy
import os
import pyfits
import sys

#-----------------------------------------------------------------------

def usage():
    print ''
    print 'NAME'
    print '       nirlin.py - NIR linearization\n'
    print 'SYNOPSIS'
    print '       nirlin.py [options] infile\n'
    print 'DESCRIPTION'
    print '       Run on raw or nprepared PANIC NIRI data, this'
    print '       script calculates and applies a per-pixel linearity'
    print '       correction based on the counts in the pixel, the'
    print '       exposure time, the read mode, the bias level, the'
    print '       ROI, and the vertical position on the detector.'
    print '       Pixels over the maximum correctable value are set'
    print '       to BADVAL unless given the force flag.'
    print ' '
    print 'OPTIONS'
    print '       -b <badval> : value to assign to uncorrectable pixels [0]'
    print '       -f : force correction on all pixels'
    print '       -o <file> : write output to <file> [l<inputfile>]'
    print '       -v : verbose debugging output\n'
    print 'VERSION'
    print '       2011 Feb 03'
    print ''
    raise SystemExit

#-----------------------------------------------------------------------

try:
    opts, args = getopt.getopt(sys.argv[1:], 'b:fo:v')
except getopt.GetoptError:
    usage()
    sys.exit(2)

nargs = len(sys.argv[1:])
nopts = len(opts)

badval = 0
force = False
outputfile = 'default'
verbose = False

for o, a in opts:
    if o in ('-b'):      # value for bad pixels (over correction limit)
        badval = a
        nopts += 1
    if o in ('-f'):      # force correction on all pixels, even if over limit
        force = True
    if o in ('-o'):      # cleaned output file
        outputfile = a
        nopts += 1
    if o in ('-v'):      # verbose debugging output
        verbose = True

if (verbose):
    print "...nargs = ", nargs
    print "...nopts = ", nopts

if (nargs - nopts) != 1:
    usage()

inputfile  = sys.argv[nopts+1]

if ( not inputfile.endswith('.fits') ):
    inputfile = inputfile + '.fits'

if (outputfile == 'default'):
    outputfile = 'l' + os.path.basename(inputfile)
else:
    if ( not outputfile.endswith('.fits') ):
        outputfile = outputfile + '.fits'

# Error checking:
if not os.path.exists(inputfile):      # check whether input file exists
    print inputfile, 'does not exist'
    sys.exit(2)
if os.path.exists(outputfile):        # check whether output file exists
    print '...removing old', outputfile
    os.remove(outputfile)

print '...reading', inputfile
hdulist = pyfits.open(inputfile)

print '...verifying...'
if (verbose):
    hdulist.verify('fix')
else:
    hdulist.verify('silentfix')

# Check if this image has already been linearized:
try:
    history = hdulist[0].header['HISTORY']
except:
    history = ''
if verbose:
    print '...history = ', history
if (history.count("Linearized by nirlin.py") > 0):
    print "ERROR: ", history
    sys.exit(2)

# Get the number of extensions in the image:
next = len(hdulist)
if verbose:
    print '...number of extensions =', next
if ( next == 1 ):
   sci = 0
else:
   sci = 1
if verbose:
    print '...assuming science data are in extension', sci

# Get the image dimensions:
try:
    naxis1,naxis2 = hdulist[sci].header['NAXIS1'],hdulist[sci].header['NAXIS2']
    print '...image dimensions =', naxis1, 'x', naxis2
except:
    print 'ERROR: cannot get the dimensions of extension ', sci
    pyfits.info(inputfile)
    sys.exit(2)

if ( naxis2 < 1024 ):
    print 'ERROR: This script does not include coefficients for subarrays.'
    print 'Please contact CAHA Observatory for more information.'
    sys.exit(2)

exptime = hdulist[0].header['EXPTIME']
print '...exposure time =', exptime, 's'

# Check that exposure time is in range:
if ( exptime > 600 ):
    print 'WARNING: exposure time is outside the range used to derive correction.'

# Read science data:
counts = hdulist[sci].data
if verbose:
    print 'INPUT DATA:'
    print counts
print '...mean of input image =', numpy.mean(counts)

# Convert to counts / coadd:
coadds = hdulist[0].header['NCOADDS']
print '...number of coadds =', coadds
if ( coadds > 1 ):
    print '...converting to counts / coadd...'
    counts = counts / coadds

# Nprepare modifies the exposure time keyword value to be exptime * ncoadds
# so if nprepared, undo this operation to get the original exposure time:
try:
    hdulist[0].header['PREPARE']
    print '...image has been nprepared'
    if ( coadds > 1 ):
        print '...converting exptime to exptime / coadd...'
        exptime = exptime / coadds
        print '...exptime = ', exptime
except:
    if verbose:
        print '...image has not been nprepared (which is okay)'

# Read mode:
lnrs = hdulist[0].header['LNRS']
#lnrs=1
print '...number of low noise reads =', lnrs 
ndavgs = hdulist[0].header['NDAVGS']
#ndavgs=1
print '...number of digital averages =', ndavgs
if   ( lnrs == 1  and ndavgs == 1  ):
    readmode = 'high-noise'
elif ( lnrs == 1  and ndavgs == 16 ):
    readmode = 'medium-noise'
elif ( lnrs == 16 and ndavgs == 16 ):
    readmode = 'low-noise'
else:
    print 'ERROR: Unknown read mode'
    sys.exit(2)
print '...read mode =', readmode

# Bias level / Well-depth:
vdduc = hdulist[0].header['A_VDDUC']
vdet  = hdulist[0].header['A_VDET']
#vdduc=0.0
#vdet=0.6
biasvoltage = vdduc - vdet
print '...bias voltage =', biasvoltage, 'V'
shallowbias = -0.60   # shallow-well detector bias (V)
deepbias    = -0.87   # deep-well detector bias (V)
if ( abs(biasvoltage - shallowbias) < 0.05 ):
    welldepth = 'shallow'
elif ( abs(biasvoltage - deepbias) < 0.05 ):
    welldepth = 'deep'
else:
    print '...ERROR: can not determine well depth.'
    sys.exit(2)
print '...well depth =', welldepth

# Set non-linearity coefficients:
if   ( readmode == 'high-noise' and naxis2 == 1024 and welldepth == 'deep' ):
    maxcounts = 20000
    a1 =  0.991500E+00
    a2 =  0.362930E-02
    a3 =  0.256083E-04
    a4 = -0.638648E+01
    a5 =  0.726802E-05
    a6 =  0.131308E-04
    a7 = -0.283632E-05

elif ( readmode == 'high-noise' and naxis2 == 2048 and welldepth == 'shallow' ):
    maxcounts = 15000
    a1 =  0.991734E+00
    a2 =  0.338144E-02
    a3 =  0.182501E-04
    a4 = -0.547003E+01
    a5 =  0.109838E-04
    a6 =  0.115787E-04
    a7 = -0.335272E-05

elif ( readmode == 'medium-noise' and naxis2 == 1024 and welldepth == 'shallow' ):
    maxcounts = 12000
    a1 =  0.992382E+00
    a2 = -0.611100E-01
    a3 =  0.337251E-05
    a4 = -0.362465E+01
    a5 =  0.102578E-04
    a6 =  0.111977E-04
    a7 = -0.360413E-05

elif ( readmode == 'low-noise' and naxis2 == 1024 and welldepth == 'shallow' ):
    maxcounts = 11000
    a1 =  0.989357E+00
    a2 = -0.104807E+01
    a3 =  0.148414E-04
    a4 = -0.490751E+01
    a5 =  0.103698E-04
    a6 =  0.244230E-05
    a7 = -0.215586E-06

else:
    print 'ERROR: coefficients do not exist for this mode.'
    print 'Please contact CAHA Observatory for more information.'
    sys.exit(2)

if (force):
    maxcounts = 65000
    print '...forcing linearity correction on all pixels'
else:
    print '...upper limit for linearization =', maxcounts, 'ADU/coadd'

# Create array of pixel y-positions:
y1 = numpy.arange(naxis2/2,0,-1).repeat(naxis2).reshape(naxis2/2,naxis2)
y2 = numpy.arange(1,naxis2/2+1).repeat(naxis2).reshape(naxis2/2,naxis2)
ypos = numpy.concatenate((y1,y2))
if verbose:
    print 'YPOS:'
    print ypos

# Remember the true counts:
truecounts = counts.copy()

# Set the minimum number of counts for calculating a correction:
counts[counts<10] = 10.

# Calculate correction factor:
correction = a1 + \
             a2 / exptime + a3 * exptime + \
             a4 / counts  + a5 * counts + \
             a6 * ypos * numpy.log10(exptime) + \
             a7 * ypos * numpy.log10(counts)

if verbose:
    print 'CORRECTION:'
    print correction

# Calculate statistics ignoring out-of-range values:
print '...', numpy.min(correction[counts<maxcounts]), \
      '< correction factors <', numpy.max(correction[counts<maxcounts])
print '...mean correction factor =', numpy.mean(correction[counts<maxcounts])
print '...median correction factor =', numpy.median(correction[counts<maxcounts])

# Apply correction:
newcounts = truecounts * correction * coadds

# Set out-of-range pixels to BADVAL:
if (not force):
    print '...setting out-of-range pixels to', badval
    newcounts[counts>maxcounts] = badval

# Write FITS data:
hdulist[sci].data = newcounts

if verbose:
    print 'OUTPUT DATA:'
    print hdulist[sci].data

print '...mean of output image =', hdulist[sci].data.mean()

print '...updating header...'
timestamp = datetime.datetime.now().strftime("%Y.%m.%d %H:%M:%S")
if (verbose):
    print '...time stamp =', timestamp
hdulist[0].header.add_history('Linearized by nirlin.py ' + timestamp)

print '...writing', outputfile
hdulist.writeto(outputfile)
hdulist.close()

#-----------------------------------------------------------------------
