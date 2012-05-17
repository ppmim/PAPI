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

def collapse(frame_list, out_filename_suffix="_c"):
    """
    Collapse (sum) a (list) of data cube into a single 2D image.

    Return a list with the new collapsed frames.
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
            log.debug("No collapse required. It is not a FITS-cube image")
            new_frame_list.append(frame_i)
        else:            
            #Suppose we have single CUBE file ...
            #sum=numpy.zeros([2048,2048],dtype='float32')
            out_hdulist = pyfits.HDUList()               
            prihdu = pyfits.PrimaryHDU (data = f[0].data.sum(0), header = f[0].header)
            prihdu.scale('float32') 
            # Updating PRIMARY header keywords...
            prihdu.header.update('NCOADDS', f[0].data.shape[0])
            prihdu.header.update('EXPTIME', f[0].header['EXPTIME']*f[0].data.shape[0])
            
            
            out_hdulist.append(prihdu)    
            #out_hdulist.verify ('ignore')
            # Now, write the new MEF file
            new_filename = frame_i.replace(".fits", out_filename_suffix+".fits")
            out_hdulist.writeto (new_filename, output_verify = 'ignore', clobber=True)
            
            out_hdulist.close(output_verify = 'ignore')
            del out_hdulist
            new_frame_list.append(new_filename)
            log.info("FITS file %s created" % (new_frame_list[n]))
            n+=1
     
    return new_frame_list

def collapse_distinguish(frame_list, out_filename_suffix="_c"):
    """
    Collapse (sum) a set of distinguish files (not cubes) into a single 2D image.

    Return a list with the new collapsed frames.
    """

    log.debug("Starting collapse_distinguish() method ....")
    
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
            log.debug("No collapse required. It is not a FITS-cube image")
            new_frame_list.append(frame_i)
        else:            
            #Suppose we have single CUBE file ...
            #sum=numpy.zeros([2048,2048],dtype='float32')
            out_hdulist = pyfits.HDUList()               
            prihdu = pyfits.PrimaryHDU (data = f[0].data.sum(0), header = f[0].header)
            prihdu.scale('float32') 
            # Updating PRIMARY header keywords...
            prihdu.header.update('NCOADDS', f[0].data.shape[0])
            prihdu.header.update('EXPTIME', f[0].header['EXPTIME']*f[0].data.shape[0])
            
            
            out_hdulist.append(prihdu)    
            #out_hdulist.verify ('ignore')
            # Now, write the new MEF file
            new_filename = frame_i.replace(".fits", out_filename_suffix+".fits")
            out_hdulist.writeto (new_filename, output_verify = 'ignore', clobber=True)
            
            out_hdulist.close(output_verify = 'ignore')
            del out_hdulist
            new_frame_list.append(new_filename)
            log.info("FITS file %s created" % (new_frame_list[n]))
            n+=1
     
    return new_frame_list
    
################################################################################
# main
################################################################################
if __name__ == "__main__":

    log.info("Start-of-collapse")

    # Get and check command-line options
        
    usage = "usage: %prog [options] arg1 arg2 ..."
    parser = OptionParser(usage)
    
    parser.add_option("-i", "--input_image",
                  action="store", dest="input_image", 
                  help="input cube image to collapse to a 2D image")
                  
    parser.add_option("-s", "--suffix",
                  action="store", dest="suffix", 
                  help="output filename suffix to add (default = %default)",
                  default="_col")
                                
    (options, args) = parser.parse_args()
    
    
    if not options.input_image or len(args)!=0: 
    # args is the leftover positional arguments after all options have been processed
        parser.print_help()
        parser.error("wrong number of arguments " )
    

    if not os.path.exists(options.input_image):
        log.error ("Input image %s does not exist", options.input_image)
        sys.exit(0)
        
    try:
        frames = [options.input_image]
        print collapse(frames, options.suffix)
    except Exception,e:
        log.info("Some error while collapsing image to 2D: %s"%str(e))
        sys.exit(0)
        
    log.info("End-of-collapse")
    