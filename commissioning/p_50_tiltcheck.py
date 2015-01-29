#!/usr/bin/env python
#########################################################################
# PANIC commissioning
# Instrument tilt adjustment: coordinate calculation, plane fit, and shim
# changes
#
# It requires the input files named 'FILTERNAME_SGX_X.txt' for the
# individual filters and the four detectors+sub-quadrants (e.g.
# Ks_SG1_3.txt). They have to be placed in the input path (optional
# argument). In case one file is not found, the area is skipped.
# These files contain the output from the starfocus log (Best focus
# estimate ) and must have two header lines of the form
# "# Object: XXX"
# "# Average best focus of 22.9332 with GFWHM of 2.43" (line from starfocus log)
# The output results are created in the optional output path, one logfile
# for each filter selection. New executions add log entries in existing
# files.
#
# usage: p_50_tiltcheck.py [-h] [-i INPUT] [-o OUTPUT] [--nolog] {1,2,3,4,5,6,7}
# 
# positional arguments:
#   {1,2,3,4,5,6,7}       Select filter: (1)Ks (2)H (3)J (4)Y (5)Z (6)H2 (7)All
# 
# optional arguments:
#   -h, --help            show this help message and exit
#   -i INPUT, --input INPUT
#                         Input path, default: '.'
#   -o OUTPUT, --output OUTPUT
#                         Output path, default: '.'
#   --nolog               Do not write output log file
# 
# 0.1 23/01/2015 BD Creation from p_49_2.2_tilt
#	1) Missing all filters and shim calculation
# 0.2 23/01/2015 BD Update
#	1) Added tilt angle and shim change calculation
#	2) Added safety check for no input files of a single filter
#	3) Added plane tilt in plot
# 1.0 23/01/2015 BD Update
#	1) Added fit of all filters
# 1.1 26/01/2015 BD Update
#	1) Added shim calculation for all filters
#
_name = 'p_50_tiltcheck.py'
_version = '1.1'
#########################################################################
import numpy as np
import os
import argparse
import matplotlib
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as axes3d
import datetime

# parse command line arguments
parser = argparse.ArgumentParser()
arghelp = 'Select filter: (1)Ks (2)H (3)J (4)Y (5)Z (6)H2 (7)All'
parser.add_argument('filter', type=int, choices=range(1,8), help=arghelp)
parser.add_argument('-i', '--input', help="Input path, default: '.'", type=str, default='.')
parser.add_argument('-o', '--output', help="Output path, default: '.'", type=str, default='.')
parser.add_argument('--nolog', help='Do not write output log file', action='store_true')
args = parser.parse_args()
ifilter = args.filter
inputpath = args.input
outputpath = args.output
writelog = not args.nolog

# FPA gap size / px
FPAgap = 167
# Radial magnification FS-FPA
FPAscale = 0.49647
# Scale M2 to FS
M2scale = 5.35

filters = ['Ks', 'H', 'J', 'Y', 'Z', 'H2', 'All']
filter = filters[ifilter - 1]
logfilename = 'Log_tiltcheck_%s.txt' %filter

# create output if necessary
if not os.access(outputpath,os.F_OK):
	os.mkdir(outputpath)
	
def nowtime():
	return datetime.datetime.now().isoformat()

def pr(msg):
	'''Print messages and write to log file'''
	print msg
	if writelog:
		logfile = open(os.path.join(outputpath, logfilename), 'a')
		logfile.write(msg + '\n')
		logfile.close()

# data arrays indices: SG1_1, 1_2, 1_3, ..., 4_4
dataindex = lambda iSG, iq: (iSG - 1) * 4 + iq -1

def loadfiledata(filter):
	'''Load data of input file of a single filter
	Input
	-----
	filter : str
			 Name of filter
	Returns
	-------
	obj : str
		  Name of the object
	subquads : list
			   List of the subquadrant names, missing ones as '--'
	pxi, pxj, focus : masked arrays
			   		  Average pixel coordinates (px) and focus data (mm) per subquadrant
	'''
	subquads = ['--'] * 16
	pxi = np.ma.masked_all(16)
	pxj = np.ma.masked_all(16)
	focus = np.ma.masked_all(16)
	obj = '--'
	for iSG in range(1, 5):
		for iq in range(1, 5):
			datafilename = '%s_SG%1i_%1i.txt' %(filter, iSG, iq)
			subquads[dataindex(iSG, iq)] = 'SG%1i_%1i' %(iSG, iq)
			try:
				datafile = open(os.path.join(inputpath, datafilename), 'r')
			except IOError:
				pr('WARNING: input file %s not found!' %datafilename)
			else:
				lines = datafile.readlines()
				# read header data: object, best focus
				tokens = lines[0].split(': ')
				obj = tokens[1][:-1]
				tokens = lines[1].split()
				focus[dataindex(iSG, iq)] = float(tokens[5])
				# read table data
				convfunc = lambda s: s.strip('(),:m=')
				data = np.loadtxt(os.path.join(inputpath, datafilename), usecols=[4, 5, 7], converters={4: convfunc, 5: convfunc, 7: convfunc})		
				# calculate relative flux and weighted pixel positions	
				flux = 10**(data[:, 2] / -2.5)
				pxi[dataindex(iSG, iq)] = (data[:, 0] * flux).sum() / flux.sum()
				pxj[dataindex(iSG, iq)] = (data[:, 1] * flux).sum() / flux.sum()
	return obj, subquads, pxi, pxj, focus

def convertdata(pxi, pxj, focus):
	'''Convert IRAF data to fieldstop coordinates
	Input
	-----
	pxi, pxj, focus : masked arrays
			   		  Average pixel and focus data per subquadrant
	Returns
	-------
	FPAx, FPAy, FSx, FSy, FSz : masked arrays
					Field coordinates and equivalent focus (mm) in the FPA and field stop
	'''
	# FPA xy in mmm
	FPAx = pxi.copy()
	FPAx[pxi <= 2048.5] = (-2048.5 + pxi[pxi <= 2048.5] - 0.5 * FPAgap) * 0.018
	FPAx[pxi > 2048.5] = (-2048.5 + pxi[pxi > 2048.5] + 0.5 * FPAgap) * 0.018
	FPAy = pxj.copy()
	FPAy[pxj <= 2048.5] = (-2048.5 + pxj[pxj <= 2048.5] - 0.5 * FPAgap) * 0.018
	FPAy[pxj > 2048.5] = (-2048.5 + pxj[pxj > 2048.5] + 0.5 * FPAgap) * 0.018
	FSx = FPAx / FPAscale
	FSy = FPAy / FPAscale
	FSz = focus * M2scale
	return FPAx, FPAy, FSx, FSy, FSz	


# Tilt of a single filter
if ifilter in range(1, 7):
	pr('#########################')
	pr('# New script executiuon #')
	pr('#########################')
	pr('# %s, %s' %(_name, _version))
	pr('# ' + nowtime())
	pr('# Analyze filter: %s' %filter)
	
	# load data from text files
	pr('# Loading input file data')
	obj, subquads, pxi, pxj, focus = loadfiledata(filter)
	if pxi.count() == 0:
		raise IOError('No input files found, aborting')
	pr('# Object: %s' %obj)
	
	# Convert to field stop coordinates via FPA xy
	pr('# Converting to field stop coordinates')
	FPAx, FPAy, FSx, FSy, FSz = convertdata(pxi, pxj, focus)
	pr('# Parameters:')
	pr('  * FPA gap size (px): %i' %FPAgap)
	pr('  * Radial magnification FS-FPA: %f' %FPAscale)
	pr('  * Focus scale M2 to FS: %f' %M2scale)
	pr('# Data table:')
	labels = ['Field quad', 'Pixel i', 'Pixel j', 'M2 focus/mm', 'FPA x/mm', 'FPA y/mm', 'FS x/mm', 'FS y/mm', 'FS focus/mm']
	header = '# '
	for label in labels: header += label + ', '
	pr(header[:-2])
	for i in range(16):
		pr('%s, %f, %f, %7.4f, %f, %f, %f, %f, %f' %(subquads[i], pxi[i], pxj[i], focus[i], FPAx[i], FPAy[i], FSx[i], FSy[i], FSz[i]))
	
	# fit plane to field stop data
	pr('# Fitting plane in field stop')
	# fit plane to points, only non-masked data
	# Where z = ax + by + d
	# Rearranging: ax + by - z + d = 0
	# Gives normal (a, b, -1)
	npts = FSx.count()
	x = FSx.compressed()
	y = FSy.compressed()
	z = FSz.compressed()
	G = np.ones((npts, 3))
	G[:, 0] = x
	G[:, 1] = y
	(a,b,d), resid, rank, s = np.linalg.lstsq(G, z)
	c = -1.
	# calculate residuals at each point
	zplane = a * x + b * y + d
	zdiff = np.abs(z - zplane)
	pr('# Fitted plane parameter A, B, C, D for Ax + By + Cz + D = 0:')
	pr(str(a))
	pr(str(b))
	pr(str(c))
	pr(str(d))
	# tilt angle in arcmin, position angle
	tilt = np.rad2deg(np.arctan(np.sqrt(a**2 + b**2) / -c)) * 60
	pr('# Plane tilt (arcmin): %5.2f' %tilt)
	angle = np.rad2deg(np.arctan2(b, a))
	pr('# Position angle (deg): %5.1f' %angle)

	# calculate the change of the shim thicknesses
	pr('# Calculating shim changes')
	# screw positions in image / mm
	TAdia = 680
	shims = np.arange(1, 13)
	shimsx = TAdia / 2 * np.sin(np.radians(15 + (shims - 1) * 30))
	shimsy = -TAdia / 2 * np.cos(np.radians(15 + (shims - 1) * 30))
	# plane positions and normalized to maximum (material to remove)
	shimsz = a * shimsx + b * shimsy
	dz = shimsz - shimsz.max()
	pr('# Shim nr, Thickness change/mm')
	for i in range(12):
		pr('%2i, %6.4f' %(i+1, dz[i]))

	# 3D plot
	pr('# Plotting points and plane fit')
	fig = plt.figure()
	# dummy contour for size colorbar
	levels = np.linspace(zdiff.min(), zdiff.max(), 256)
	ct = plt.contourf([[0,0], [0,0]], levels)
	plt.clf()
	ax = fig.add_subplot(111, projection='3d')
	FSmax = 77.261
	xlim = [-FSmax, FSmax]
	ylim = [-FSmax, FSmax]
	# scale data to color values
	cmap = plt.get_cmap()
	normalize = matplotlib.colors.Normalize()
	colorvals = normalize(zdiff)
	# plot points
	for ipt in range(npts):
		ax.plot([x[ipt]], [y[ipt]], [z[ipt]], 'o', color=cmap(colorvals[ipt]), zorder=10)
	# fitted plane
	xx, yy = np.meshgrid(xlim, ylim)
	zz = a * xx + b * yy + d
	ax.plot_surface(xx, yy, zz, alpha=0.1, color=[0,1,0])
	zmin, zmax = ax.get_zlim()
	# z drop lines
	for ipt in range(npts):
		ax.plot([x[ipt], x[ipt]], [y[ipt], y[ipt]], [z[ipt], zplane[ipt]], ':', color='0.5', zorder=1)
		ax.plot([x[ipt]], [y[ipt]], [zplane[ipt]], '.', color='0.5', zorder=1)
	# side projections
	for ipt in range(npts):
		# to back x
		ax.plot([-FSmax], [y[ipt]], [z[ipt]], '.', color='0.7', zorder=1)
		# to back y
		ax.plot([x[ipt]], [FSmax], [z[ipt]], '.', color='0.7', zorder=1)
		# at bottom
		ax.plot([x[ipt]], [y[ipt]], [zmin], '.', color='0.7', zorder=1)
	# SG cross
	ax.plot([0, 0], ylim, [zmin, zmin], 'g', zorder=0)
	ax.plot(xlim, [0, 0], [zmin, zmin], 'g', zorder=0)
	
	ax.set_xlim(xlim)
	ax.set_ylim(ylim)
	plt.suptitle(r'FS focus positions %s %s' %(obj, filter))
	ax.set_xlabel('FS x / mm')
	ax.set_ylabel('FS y / mm')
	ax.set_zlabel('FS focus / mm')
	# colorbar for sizes
	cbar = plt.colorbar(ct)
	cbar.set_label('z residuum / mm')
	plt.tight_layout()
	plt.figtext(0.02, 0.03, "Plane tilt angle: %5.2f'" %tilt, ha='left', style='italic', size=8)
	plt.figtext(0.02, 0.01, "Plane position angle: %5.1f'" %angle, ha='left', style='italic', size=8)
	filename = '%s_%s_tiltcheck' %(obj, filter)
	plt.savefig(os.path.join(outputpath, filename + '.png'))
	plt.close()

# Plane of all filters
elif ifilter in [7]:
	pr('#########################')
	pr('# New script executiuon #')
	pr('#########################')
	pr('# %s, %s' %(_name, _version))
	pr('# ' + nowtime())
	pr('# Analyze filter: %s' %filter)
	pr('# Parameters:')
	pr('  * FPA gap size (px): %i' %FPAgap)
	pr('  * Radial magnification FS-FPA: %f' %FPAscale)
	pr('  * Focus scale M2 to FS: %f' %M2scale)
	# loop over filters, create single data set and plot
	# all points
	fulldata = np.empty((0, 9))
	# of single filter
	filterdata = dict()
	# list of available filters
	fullfilters = []
	# tilt and position angles
	fulltilts = np.empty((0, 2))
	# plane parameters
	fullplanes = np.empty((0, 4))
	# 3D plot of single filters
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	FSmax = 77.261
	xlim = [-FSmax, FSmax]
	ylim = [-FSmax, FSmax]
	
	# load data from text files
	pr('# Loading input file data')
	obj = '--'
	for sfilter in filters[0:6]:	
		pr('# Filter: %s' %sfilter)
		objfile, subquads, pxi, pxj, focus = loadfiledata(sfilter)
		if pxi.count() == 0:
			pr('WARNING: No input files found for filter %s' %sfilter)
			continue
		obj = objfile
		pr('  * Object: %s' %obj)
		# Convert to field stop coordinates via FPA xy
		pr('  * Converting to field stop coordinates')
		FPAx, FPAy, FSx, FSy, FSz = convertdata(pxi, pxj, focus)

		# fit plane to field stop data
		pr('  * Fitting plane in field stop')
		# fit plane to points, only non-masked data
		# Where z = ax + by + d
		# Rearranging: ax + by - z + d = 0
		# Gives normal (a, b, -1)
		npts = FSx.count()
		x = FSx.compressed()
		y = FSy.compressed()
		z = FSz.compressed()
		G = np.ones((npts, 3))
		G[:, 0] = x
		G[:, 1] = y
		(a,b,d), resid, rank, s = np.linalg.lstsq(G, z)
		c = -1.
		# tilt angle in arcmin, position angle
		tilt = np.rad2deg(np.arctan(np.sqrt(a**2 + b**2) / -c)) * 60
		angle = np.rad2deg(np.arctan2(b, a))
		# Filter data with subtracted offset
		data = np.ma.array([pxi, pxj, focus, FPAx, FPAy, FSx, FSy, FSz, FSz - d]).transpose()
		# plot points
		ax.plot(x, y, z-d, 'o', zorder=10, label='%s' %(sfilter))
		# result collections
		fullfilters.append(sfilter)
		fulldata = np.vstack((fulldata, np.ma.compress_rows(data)))
		filterdata[sfilter] = data
		fulltilts = np.vstack((fulltilts, [tilt, angle]))
		fullplanes = np.vstack((fullplanes, [a, b, c, d]))
	
	pr('# Data tables of all filters')
	for sfilter in fullfilters:
		pr('  * Filter: %s' %sfilter)
		labels = ['Field quad', 'Pixel i', 'Pixel j', 'M2 focus/mm', 'FPA x/mm', 'FPA y/mm', 'FS x/mm', 'FS y/mm', 'FS focus/mm', 'FS focus - D/mm']
		header = '# '
		for label in labels: header += label + ', '
		pr(header[:-2])
		for i in range(16):
			pr('%s, %f, %f, %7.4f, %f, %f, %f, %f, %f, %f' %((subquads[i],) + tuple(nmbr for nmbr in filterdata[sfilter][i, :])))

	pr('# Fitting plane to all points')
	[rows,cols] = fulldata.shape
	x = fulldata[:, -4]
	y = fulldata[:, -3]
	z = fulldata[:, -1]
	G = np.ones((rows, 3))
	G[:, 0] = x
	G[:, 1] = y
	(a,b,d), resid, rank, s = np.linalg.lstsq(G, z)
	c = -1.
	# tilt angle in arcmin, position angle
	tilt = np.rad2deg(np.arctan(np.sqrt(a**2 + b**2) / -c)) * 60
	angle = np.rad2deg(np.arctan2(b, a))
	# add to result collection
	fullfilters.append(filter)
	fulltilts = np.vstack((fulltilts, [tilt, angle]))
	fullplanes = np.vstack((fullplanes, [a, b, c, d]))
	# relative offsets at M2, relative to first (Ks?)
	DM2 = fullplanes[:, 3] / M2scale
	DM2 -= DM2[0]
	
	pr('# Plane results of all filters')
	labels = ['Filter', 'Plane A', 'Plane B', 'Plane C', 'Plane D', "Tilt/'", 'Pos angle/deg', 'rel D at M2/um']
	header = '# '
	for label in labels: header += label + ', '
	pr(header[:-2])
	for i, sfilter in enumerate(fullfilters):
		pr('%s, %s, %s, %s, %s, %5.2f, %5.1f, %5.1f' %((sfilter,) + tuple(nmbr for nmbr in fullplanes[i, :]) + (fulltilts[i, 0], fulltilts[i, 1]) + (DM2[i]*1000,)))
				
	# fitted plane
	xx, yy = np.meshgrid(xlim, ylim)
	zz = a * xx + b * yy + d
	ax.plot_surface(xx, yy, zz, alpha=0.1, color=[0,1,0])
	zmin = ax.get_zlim()[0]
	zplane = a * x + b * y + d
	# z drop lines
	zmin = ax.get_zlim()[0]
	for ipt in range(rows):
		ax.plot([x[ipt], x[ipt]], [y[ipt], y[ipt]], [z[ipt], zplane[ipt]], ':', color='0.5', zorder=1)
		ax.plot([x[ipt]], [y[ipt]], [zplane[ipt]], 'x', color='0.5', zorder=11)
	# side projections
	for ipt in range(rows):
		# at bottom
		ax.plot([x[ipt]], [y[ipt]], [zmin], '.', color='0.7', zorder=1)
	# SG cross
	ax.plot([0, 0], ylim, [zmin, zmin], 'g', zorder=0)
	ax.plot(xlim, [0, 0], [zmin, zmin], 'g', zorder=0)	

	ax.set_xlim(xlim)
	ax.set_ylim(ylim)
	# 	ax.set_zlim(zmin=zmin)
	plt.suptitle(r'FS focus positions %s all filters' %(obj))
	ax.set_xlabel('FS x / mm')
	ax.set_ylabel('FS y / mm')
	ax.set_zlabel('Rel. FS focus / mm')
	plt.legend(loc='upper left', numpoints=1)
	plt.tight_layout()
	plt.figtext(0.02, 0.03, "Plane tilt angle: %5.2f'" %tilt, ha='left', style='italic', size=8)
	plt.figtext(0.02, 0.01, "Plane position angle: %5.1f'" %angle, ha='left', style='italic', size=8)
	filename = '%s_%s_tiltcheck' %(obj, filter)
	plt.savefig(os.path.join(outputpath, filename + '.png'))
	plt.close()

	# calculate the change of the shim thicknesses
	pr('# Calculating shim changes')
	# screw positions in image / mm
	TAdia = 680
	shims = np.arange(1, 13)
	shimsx = TAdia / 2 * np.sin(np.radians(15 + (shims - 1) * 30))
	shimsy = -TAdia / 2 * np.cos(np.radians(15 + (shims - 1) * 30))
	# plane positions and normalized to maximum (material to remove)
	shimsz = a * shimsx + b * shimsy
	dz = shimsz - shimsz.max()
	pr('# Shim nr, Thickness change/mm')
	for i in range(12):
		pr('%2i, %6.4f' %(i+1, dz[i]))
