#!/usr/bin/env python

################################################################################
#
#
# test.py
#
# Last update 11/03/2009
#
################################################################################

import os
import pyraf
import iraf

# Pyraf modules
from pyraf import iraf
from iraf import noao
from iraf import imred
from iraf import ccdred

# Interact with FITS files
import pyfits

lista='/disk-a/caha/panic/DATA/aa0001.fits, /disk-a/caha/panic/DATA/aa0002.fits,'
lista='/disk-a/caha/panic/DATA/data_mat/dark0001.fits , /disk-a/caha/panic/DATA/data_mat/dark0002.fits , /disk-a/caha/panic/DATA/data_mat/dark0003.fits , /disk-a/caha/panic/DATA/data_mat/dark0004.fits , /disk-a/caha/panic/DATA/data_mat/dark0005.fits'

sky_bkg='sflat.fits'
master_flat='/tmp/master_dflat.fits'
sky_bkg_n='sflatn.fits'
#########################################################
##### Option 1 # substract sky divide by twilight flat
#########################################################
## Substract sky
source= open( "files.nip", "r" )
destination= open( "files_sub_sky.nip", "w" )
for line in source:
    line=line.replace(".fits", "sub_sky.fits")
    destination.write( line )
source.close()
destination.close()

iraf.imarith(operand1 = '@files.nip',
                    operand2 = sky_bkg,
                    op = '-',
                    result ='@files_sub_sky.nip',
                    verbose = 'yes'
                    )
                    
## Divide by flat
source= open( "files_sub_sky.nip", "r" )
destination= open( "files_sub_sky_flatted.nip", "w" )
for line in source:
    line=line.replace(".fits", "_flatted.fits")
    destination.write( line )
source.close()
destination.close()

iraf.imarith(operand1 = '@files_sub_sky.nip',
                    operand2 = master_flat,
                    op = '/',
                    result ='@files_sub_sky_flatted.nip',
                    verbose = 'yes'
                    )
                               
#########################################################
##### Option 2 # divide by sky-flat (superflat)
#########################################################
## Substract sky
source= open( "files.nip", "r" )
destination= open( "files_flat_sky.nip", "w" )
for line in source:
    line=line.replace(".fits", "flat_sky.fits")
    destination.write( line )
source.close()
destination.close()

iraf.imarith(operand1 = '@files.nip',
                    operand2 = sky_bkg_n,
                    op = '/',
                    result ='@files_flat_sky.nip',
                    verbose = 'yes'
                    )
                    
"""
## Divide by flat
source= open( "files_sub_sky.nip", "r" )
destination= open( "files_sub_sky_flatted.nip", "w" )
for line in source:
    line=line.replace(".fits", "_flatted.fits")
    destination.write( line )
source.close()
destination.close()

iraf.imarith(operand1 = '@files_sub_sky.nip',
                    operand2 = master_flat,
                    op = '/',
                    result ='@files_sub_sky_flatted.nip',
                    verbose = 'yes'
                    )
"""  
                       
