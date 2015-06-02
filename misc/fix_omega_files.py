#!/usr/bin/env python


import datetime
import os
import sys
import fileinput
import astropy.io.fits as fits

input_file = "/tmp/omega_raw.txt"
filelist = [line.replace( "\n", "") 
                    for line in fileinput.input(input_file)]

for i in xrange(0, len(filelist), 4):
    # Modify 1st file
    with fits.open(filelist[i], 'update') as f:
        print "Updating file :",filelist[i]
        f[0].header.set("IMAGETYP", "SKY")
    if (i + 1) < len(filelist):
        # 2nd file
        print "Updating file :",filelist[i + 1]
        with fits.open(filelist[i + 1], 'update') as f:
            f[0].header.set("IMAGETYP", "SKY")
