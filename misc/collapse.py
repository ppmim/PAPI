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
import astropy.io.fits as fits
import fileinput
import numpy

# Logging
from misc.paLog import log
from misc.version import __version__


def collapse(frame_list, out_dir="/tmp"):
    """
    Collapse (add them up arithmetically) a (list) of data cubes into a single 
    2D image.

    Return a list with the new collapsed  frames.
    """

    log.debug("Starting collapse() method ....")
    
    new_frame_list = [] 
    n = 0

    if frame_list == None or len(frame_list) == 0 or frame_list[0] == None:
        return []

    for frame_i in frame_list:
        f = fits.open(frame_i)
        # First, we need to check if we have MEF files
        if len(f)>1 and len(f[1].data.shape)==3:
            try:
                log.info("Collapsing a MEF cube %s"%frame_i)
                outfile = out_dir + "/" + os.path.basename(frame_i).replace(".fits", "_coadd.fits")
                out = collapse_mef_cube(frame_i, outfile)
                new_frame_list.append(out)
            except Exception,e:
                log.error("Some error collapsing MEF cube: %s"%str(e))
                raise e
        elif len(f)>1 and len(f[1].data.shape)==2:
            log.debug("MEF file has no cubes, no collapse required.")
            new_frame_list.append(frame_i)
        elif len(f[0].data.shape)!=3: # 2D !
            log.debug("It is not a FITS-cube image, no collapse required.")
            new_frame_list.append(frame_i)
        else:            
            # Suppose we have single CUBE file ...
            out_hdulist = fits.HDUList()               
            prihdu = fits.PrimaryHDU (data = f[0].data.sum(0), header = f[0].header)
            prihdu.scale('float32') 
            # Updating PRIMARY header keywords...
            prihdu.header.set('NCOADDS', f[0].data.shape[0])
            prihdu.header.set('EXPTIME', f[0].header['EXPTIME']*f[0].data.shape[0])
            prihdu.header.set('PAPIVERS', __version__, "PANIC Pipeline version")

            out_hdulist.append(prihdu)    
            #out_hdulist.verify ('ignore')
            # Now, write the new collapsed file
            t_filename = out_dir + "/" + os.path.basename(frame_i).replace(".fits", "_coadd.fits")
            out_hdulist.writeto (t_filename, output_verify = 'ignore', 
                                 clobber=True)
            
            out_hdulist.close(output_verify = 'ignore')
            del out_hdulist
            new_frame_list.append(t_filename)
            log.info("FITS file %s created" % (new_frame_list[n]))
            n+=1
     
    return new_frame_list

def collapse_mef_cube(inputfile, out_filename=None):
    """
    Collapse each of the extensions of a MEF file
    """

    f = fits.open(inputfile)

    out_hdulist = fits.HDUList()
    prihdu = fits.PrimaryHDU (data = None, header = f[0].header)
    prihdu.header.set('NCOADDS', f[1].data.shape[0])
    prihdu.header.set('EXPTIME', f[0].header['EXPTIME']*f[1].data.shape[0])
    prihdu.header.set('PAPIVERS', __version__, "PANIC Pipeline version")
    out_hdulist.append(prihdu)    
 
    # Sum each extension
    for ext in range(1,len(f)):
        hdu = fits.ImageHDU (data = f[ext].data.sum(0), header = f[ext].header)
        hdu.scale('float32') 
        out_hdulist.append(hdu)    
    
    # Now, write the new collapsed file
    if out_filename==None:
        outfile = inputfile.replace(".fits", "_coadd_%s.fits"%str(f[1].data.shape[0]).zfill(3))
    else:
        outfile = out_filename 
    out_hdulist.writeto (outfile, output_verify = 'ignore', 
                             clobber=True)
        
    out_hdulist.close(output_verify = 'ignore')
    del out_hdulist
    log.info("FITS file %s created" % (outfile))

    return outfile

def collapse_distinguish(frame_list, out_filename="/tmp/collapsed.fits"):
    """
    Collapse (sum) a set of distinguish files (not cubes) into a single 2D image.

    Return the name of the output file created.
    
    Curretly not used from PAPI, **only** from command-line
    """

    log.debug("Starting collapse_distinguish() method ....")
    
    new_frame_list = [] 
    if frame_list == None or len(frame_list) == 0 or frame_list[0] == None:
        return []

    for frame_i in frame_list:
        f = fits.open(frame_i)
        # First, we need to check if we have MEF files
        if len(f)>1 and len(f[1].data.shape)==3:
            log.error("MEF-cubes files cannot be collapsed. First need to be split !")
            raise Exception("MEF-cubes files cannot be collapsed. First need to be split !")
        elif len(f)>1 and len(f[1].data.shape)==2:
            log.error("Not implemented yet.")
            raise Exception("MEF-2D files cannot be collapsed. Not implemented yet.")
            ## TO BE COMPLETED !!! ##
            log.debug("Found a MEF with a 2D-image: %s:"%frame_i)
            new_frame_list.append(frame_i)
            if len(new_frame_list)==1:
                sum = numpy.zeros(len(f)-1, [f[1].header['NAXIS1'], 
                                   f[1].header['NAXIS2']], dtype='float32')
                header1 = f[0].header
            for i in range(len(f)):
                sum[i] += f[i+1].data
        elif len(f)==1 and len(f[0].data.shape)==2:
            log.debug("Found a 2D-image: %s:"%frame_i)
            new_frame_list.append(frame_i)
            if len(new_frame_list)==1:
                sum = numpy.zeros([f[0].header['NAXIS1'], 
                                   f[0].header['NAXIS2']], dtype='float32')
                header1 = f[0].header
            sum += f[0].data
            f.close()
        
    # Now, save the collapsed set of files in a new single file        
    out_hdulist = fits.HDUList()
                   
    prihdu = fits.PrimaryHDU (data = sum, header = header1)
    prihdu.scale('float32') 
        
    # Updating PRIMARY header keywords...
    prihdu.header.set('NCOADDS', len(new_frame_list))
    prihdu.header.set('EXPTIME', header1['EXPTIME']*len(new_frame_list))
    prihdu.header.set('PAPIVERS', __version__, "PANIC Pipeline version")
    
    out_hdulist.append(prihdu)    
    #out_hdulist.verify ('ignore')
    # Now, write the new single FITS file
    out_hdulist.writeto(out_filename, output_verify = 'ignore', clobber=True)
    
    out_hdulist.close(output_verify = 'ignore')
    del out_hdulist
    log.info("FITS file %s created" % (out_filename))
     
    return out_filename
    
################################################################################
# main
################################################################################
if __name__ == "__main__":

    log.info("Start-of-collapse")

    # Get and check command-line options
        
    USAGE = "usage: %prog [options] arg1 arg2 ..."
    desc = "Collapse (add them up arithmetically) each cube of a list files into a single 2D image"
    
    parser = OptionParser(USAGE, description=desc)
    
    parser.add_option("-i", "--input_image",
                  action="store", dest="input_image", 
                  help="input cube image to collapse into a 2D image")

    parser.add_option("-l", "--input_image_list",
                  action="store", dest="input_image_list", 
                  help="input image list to collapse into a single 2D image")

    parser.add_option("-o", "--output_file",
                  action="store", dest="output_file", 
                  help="output filename (default = %default)",
                  default="/tmp/out.fits")

    parser.add_option("-d", "--output_dir",
                  action="store", dest="output_dir", 
                  help="output directory (default = %default)",
                  default="/tmp")
                                
    (options, args) = parser.parse_args()
    
    if len(sys.argv[1:])<1:
       parser.print_help()
       sys.exit(0)

    if (not options.input_image and not options.input_image_list) or len(args)!=0: 
    # args is the leftover positional arguments after all options have been processed
        parser.print_help()
        parser.error("Wrong number of arguments " )
    
    if options.input_image and options.input_image_list:
    # only one option can be executed 
        parser.print_help()
        parser.error("Only one option can be used")
        
    if options.input_image:
        if not os.path.exists(options.input_image):
            log.error("Input image %s does not exist", options.input_image)
            sys.exit(0)
        if not options.output_dir or not os.path.exists(options.output_dir):
            parser.print_help()
            parser.error("Wrong number of arguments " )

        try:
            frames = [options.input_image]
            print collapse(frames, options.output_dir)
        except Exception, e:
            log.info("Some error while collapsing image to 2D: %s"%str(e))
            sys.exit(0)
    elif options.input_image_list:
        try:
            frames = [line.replace("\n", "") for line in 
                      fileinput.input(options.input_image_list)]
            print collapse_distinguish(frames, options.output_file)
        except Exception, e:
            log.info("Some error while collapsing image to 2D: %s"%str(e))
            sys.exit(0)
        
    log.info("End-of-collapse")
    