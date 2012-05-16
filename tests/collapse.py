#!/usr/bin/env python

# Copyright (c) 2012 IAA-CSIC  - All rights reserved. 
# Author: Jose M. Ibanez. 
# Instituto de Astrofisica de Andalucia, IAA-CSIC
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



# Import necessary modules
from optparse import OptionParser
import sys
import os
import pyfits

# Logging
from misc.paLog import log
import misc.utils

def collapse(   frame_list, out_filename_suffix="_sum"):
    """
    Collapse (sum) the data cube into a single 2D image

    Return a list with the new collapsed frames
    """

    log.debug("Starting collapse() method ....")

    new_frame_list = [] 
    n = 0

    if frame_list==None or len(frame_list)==0 or frame_list[0]==None:
        return []

    for frame_i in frame_list:
        f = pyfits.open(frame_i)
        # First, we need to check if we have MEF files
        if len(f)>1:
            log.error("MEF files cannot be collapsed. First need to be splitted !")
            raise Exception("MEF files cannot be collapsed. First need to be splitted !")
            new_frame_list.append(frame_i)
        elif len(f[0].data.shape)!=3:
            log.info("No collapse required. It is not a cube image")
            new_frame_list.append(frame_i)
        else:            
            #Suppose we have single CUBE file ...
            #sum=numpy.zeros([2048,2048],dtype='float32')
            out_hdulist = pyfits.HDUList()               
            prihdu = pyfits.PrimaryHDU (data = f[0].data.sum(0), header = f[0].header)
            prihdu.scale('float32') 
            # Start by updating PRIMARY header keywords...
            #prihdu.header.update ('EXTEND', pyfits.FALSE, after = 'NAXIS')
            
            out_hdulist.append(prihdu)    
            #out_hdulist.verify ('ignore')
            # Now, write the new MEF file
            new_filename = frame_i.replace(".fits","_sum.fits")
            out_hdulist.writeto (new_filename, output_verify = 'ignore', clobber=True)
            
            out_hdulist.close(output_verify = 'ignore')
            del out_hdulist
            new_frame_list.append(new_filename)
            log.info("FITS file %s created" % (new_frame_list[n]))
            n+=1
            

if __name__ == "__main__":

    frames = ['/home/panic/DATA/IAA_test_0008.fits']
    
    log.info("Start-of-collapse")

    collapse(frames)
    
    log.info("End-of-collapse")