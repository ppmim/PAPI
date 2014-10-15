#!/usr/bin/env python
#########################################################################
# PANIC data processing
#########################################################################
# Example of nonlinearity correction for full frame MEF files
# Includes many sanity checks, the actual calculation is very simple
# Incorrectable and saturated pixels are set to NaN. Therefore the ouptut
# data is float.
#
# 1.0 14/07/2014 BD Creation
#
_name = 'correct_nonlinearity.py'
_version = '1.0'
#########################################################################
import numpy as np
from astropy.io import fits
import os
import dateutil
import dateutil.parser

# path to data and files
inputpath = '/data1/PANIC/testNLC'
nonlinpath = '/data1/Calibs'
outputpath = '.'
inputfile = 'Test_MEF_0001.fits'

# inputpath = '/Users/dorner/work/PANIC/Exposures'
# nonlinpath = 'output'

if not os.access(outputpath,os.F_OK):
	os.mkdir(outputpath)

def polyval_map(poly, map):
	'''Evaluate individual polynomials on an array. Looping over each pixel
	is stupid, therefore we loop over the order and calculate the
	polynomial directly.
	Note: The output is a float array!
	
	Input
	-----
	poly : array_like
		   Polynomial coefficients without constant offset. The order
		   must be along the last axis.
	map : array_like
		  Data array.
		  
	Returns
	-------
	polymap : array_like
			  Result of evaluation, same shape as map, dtype float
	'''
	order = poly.shape[-1]
	polymap = map * 0.
	for io in range(order):
		polymap += poly[Ellipsis, -io-1] * map**(io+1)
	return polymap

# load image data and header
hdulist = fits.open(os.path.join(inputpath, inputfile))
dataheader = hdulist[0].header
# load the appropriate correction file
if dataheader['READMODE'] == 'line.interlaced.read':
# 	nonlinfile = '29_nonlinearity_data_lir/mNONLIN_LIR_00.01.fits'
	nonlinfile = 'mNONLIN_LIR_00.01.fits'
elif dataheader['READMODE'] == 'fast-reset-read.read':
# 	nonlinfile = '29_nonlinearity_data_rrr-mpia/mNONLIN_RRR-MPIA_00.01.fits'
	nonlinfile = 'mNONLIN_RRR-MPIA_00.01.fits'
else:
	raise ValueError('No nonlinearity correction available for read mode %s' %dataheader['READMODE'])
nlhdulist = fits.open(os.path.join(nonlinpath, nonlinfile))
nlheader = nlhdulist[0].header

# check header data
if nlheader['PAPITYPE'] != 'MASTER_LINEARITY':
	raise ValueError('Wrong type of nonlinearity correction file')
datadate = dateutil.parser.parse(dataheader['DATE-OBS'])
nldate = dateutil.parser.parse(nlheader['USE_AFT'])
if datadate < nldate:
	raise ValueError('Nonlinearity calibration too new for input data')
keys = ['INSTRUME', 'PREAD', 'PSKIP', 'LSKIP', 'READMODE', 'IDLEMODE', 'IDLETYPE']
for key in keys:
	if dataheader[key] != nlheader[key]:
		raise ValueError('Mismatch in header data for keyword \'%s\'' %key)
# some may not be present in old data
keys = ['DETROT90', 'DETXYFLI']
for key in keys:
	if not key in dataheader:
		print 'Warning: key \'%s\' not in data header' %key
	elif dataheader[key] != nlheader[key]:
		raise ValueError('Mismatch in header data for keyword \'%s\'' %key)
keys = ['B_EXT', 'B_DSUB', 'B_VREST', 'B_VBIAG']
for key in keys:
	for i in range (1, 5):
		if dataheader[key + '%i'%i] != nlheader[key + '%i'%i]:
			raise ValueError('Mismatch in header data for keyword \'%s%i\'' %(key, i))

# output fits HDU
linhdu = fits.PrimaryHDU()
linhdu.header = dataheader.copy()
hdus = []
# loop over detectors
for iSG in range(1, 5):
	extname = 'SG%i_1' %iSG
	# check detector sections
	# another way would be to loop until the correct one is found
	datadetsec = hdulist[extname].header['DETSEC']
	nldetsec = nlhdulist['LINMAX%i' %iSG].header['DETSEC']
	if datadetsec != nldetsec:
		raise ValueError('Mismatch of detector sections for SG%i' %iSG)
	# or check SG IDs (as long as they are reliable)
	datadetid = hdulist[extname].header['DET_ID']
	nldetid = nlhdulist['LINMAX%i' %iSG].header['DET_ID']
	if datadetid != nldetid:
		raise ValueError('Mismatch of detector IDs for extension' %extname)

	# load file data
	data = hdulist[extname].data
	nlmaxs = nlhdulist['LINMAX%i' %iSG].data
	nlpolys = np.rollaxis(nlhdulist['LINPOLY%i' %iSG].data, 0, 3)
	# calculate linear corrected data
	lindata = polyval_map(nlpolys, data)
	# mask saturated inputs - to use nan it has to be a float array
	lindata[data > nlmaxs] = np.nan
	# mask where max range is nan
	lindata[np.isnan(nlmaxs)] = np.nan
	exthdu = fits.ImageHDU(lindata.astype('float32'), header=hdulist[extname].header.copy())
	# this may rearrange the MEF extensions, otherwise loop over extensions
	hdus.append(exthdu)

# add some info in the header
linhdu.header['HISTORY'] = 'Nonlinearity correction applied'
linhdu.header['HISTORY'] = 'Nonlinearity data: %s' %nlheader['ID']
linhdu.header['HISTORY'] = '<-- The German team made this on 2014/07/13'
linhdulist = fits.HDUList([linhdu] + hdus)
linhdulist.writeto(os.path.join(outputpath, inputfile), clobber=True)
	
