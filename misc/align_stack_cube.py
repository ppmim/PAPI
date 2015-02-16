# #! /usr/bin/env python
#encoding:UTF-8

# Copyright (c) 2015 Jose M. Ibanez All rights reserved.
# Institute of Astrophysics of Andalusia, IAA-CSIC
#
# This file is part of PAPI
#
# PAPI is free software: you can redistribute it and/or modify
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


import sys
import os
import errno
import shutil
import tempfile
from optparse import OptionParser
import fileinput
import astropy.io.fits as fits
import numpy

# PAPI modules
import astromatic
import datahandler
from reduce.makeobjmask import *

# Logging
from misc.paLog import log
from misc.version import __version__
import misc.config
import misc.utils



def align_stack_cube(input_file, output_file=None):
    """
    Align and stack (arithmetic sum) the planes of a cube of frames.
    
    Parameters
    ----------
    input_file: str
        Filename of FITS cube file to be aligned and stacked.
    
    output_file: str
        Filename of the aligned and stacked image (FITS).
        
    Returns
    -------
    output_file: str
        Filename of the aligned and stacked image (FITS).
        
    """
    
    log.debug("Start align_stack_cube...")
    # Split planes into temporal fits files
    sp_files = split_fits_cube(input_file, os.path.split(output_file)[0])
    
    log.debug("Lets get offsets....")
    # Get offsets
    offsets = getPointingOffsets(sp_files)
    
    print "OFFSETS: \n %s", offsets
    
    # Aling and coadd the images
    # TODO
    
def split_fits_cube(filename, out_path, overwrite=True):
    """
    Split up the slices of a FITS cube into individual files
    
    Parameters
    ----------
    filename: str
        Name of the FITS cube file to split.
    
    overwrite: bool
        If True, overwrite output filename if it exits.
    
    Returns
    -------
    List of individual filenames generated.
    
    """
    
    print "OUT_PATH",out_path
    imfiledir = os.path.join( out_path, os.path.splitext(os.path.split(filename)[1])[0] )
    print "IMFILEDIR",imfiledir
    mkdir_p( imfiledir )
    
    out_filenames = []
    with fits.open(filename) as fcube:
        if len(fcube)>1:
            msg = "MEF files are currently not supported !"
            log.error(msg)
            raise Exception(msg)
        
        numslices = fcube[0].data.shape[0]
        print "FILE", filename
        for i in range(numslices):
            imfileout = os.path.join( imfiledir, os.path.splitext(os.path.split(filename)[1])[0] + \
                '_' + ( "%05d" % i ) + '.fits' )
            if overwrite and (os.path.exists(imfileout) or os.path.islink(imfileout)):
                os.remove(imfileout)
            hdu = fits.PrimaryHDU( header=fcube[0].header, data=fcube[0].data[i] )
            hdu.writeto(imfileout)
            out_filenames.append(imfileout)
    
    return out_filenames

def getPointingOffsets (image_list, 
                            p_offsets_file='/tmp/offsets.pap',
                            out_path='/tmp'):
        """
        Derive pointing offsets between each image using SExtractor OBJECTS 
        (makeObjMask) and offsets (IRDR)
        
        Parameters
        ----------
                
        image_list: list
            list of FITS files.
            
        p_offsets_file: str
            filename of file where offset will be saved.
        
        Returns
        -------    
        offsets: narray          
                two dimensional array with offsets
                
        Notes:
            It assumed that North is up and East is left.
        """
            
        log.info("Starting getPointingOffsets....")
        
        offsets_mat = None
        out_dir =  out_path
        m_irdr_path = '/home/panic/DEVELOP/papi/irdr/bin'
        
        # STEP 1: Create SExtractor OBJECTS images
        output_fd, output_list_file = tempfile.mkstemp(suffix='.pap', 
                                                       dir=out_dir)
        os.close(output_fd)
        os.unlink(output_list_file) # we only need the name
        
        log.debug("Creating OBJECTS images (SExtractor)....")
        
        mask_minarea = 5
        mask_maxarea = 200
        mask_thresh = 1.5
        satur_level = 300000
        single_p = False
        MIN_CORR_FRAC = 0.1
            
        # In order to set a real value for satur_level, we have to check the
        # number of coadds of the images (NCOADDS or NDIT keywords).
        try:
            pf = datahandler.ClFits(image_list[0])
            satur_level = int(satur_level) * int(pf.getNcoadds())
            log.critical("SAT_LEVEL=%s"%satur_level)
        except:
            log.warning("Error read NCOADDS value. Taken default value (=1)")
            satur_level = satur_level
            log.critical("SAT_LEVEL(def)=%s"%satur_level)

        # Make mask
        try:
            misc.utils.listToFile(image_list, out_dir + "/makeObjMask.list")
            makeObjMask( out_dir + "/makeObjMask.list", mask_minarea, mask_maxarea, mask_thresh, 
                         satur_level, output_list_file, single_point=single_p)
        except Exception,e:
            log.error("Error making object mask")
            raise e
        
        # STEP 2: Compute dither offsets (in pixles) using cross-correlation technique ==> offsets
        #>mosaic objfiles.nip $off_err > offsets1.nip
        search_box = 10 # half_width of search box in arcsec (default 10)
        offsets_cmd = m_irdr_path + '/offsets '+ output_list_file + '  ' + str(search_box) + ' >' + p_offsets_file
        if misc.utils.runCmd( offsets_cmd ) == 0:
            log.critical("Some error while computing dither offsets")
            raise Exception("Some error while computing dither offsets")
        else:
            try:
                offsets_mat = numpy.loadtxt(p_offsets_file, usecols = (1,2,3)) # columns => (xoffset, yoffset, match fraction) in PIXELS
                # check if correlation overlap fraction is good enough for all offsets computed
                if (offsets_mat[:,2] < MIN_CORR_FRAC).sum() > 0:
                    log.critical("Some error while computing dither offsets. Overlap correlation fraction is < %f", MIN_CORR_FRAC)
                    raise Exception("Wrong overlap correlation fraction for translation offsets")
                    
            except IOError:
                log.critical("No offsets read. There may be some problem while computing translation offsets ....")
                raise Exception("No offsets read")
        
        log.debug("END of getPointingOffsets")                        
        return offsets_mat

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise
    
################################################################################
# main
################################################################################
if __name__ == "__main__":
    
    log.debug( 'Start AlignStackCube....')
    
    # Get and check command-line options
    usage = "usage: %prog [options]"
    desc = """Performs the alignment and Stack (sum) the planes of a given cube,
in principle, with small drift of the planes. No resampling is done.
"""
    parser = OptionParser(usage, description=desc)
    
    parser.add_option("-c", "--config_file",
                  action="store", dest="config_file", 
                  help="Mandatory PAPI configuration file.")
                  
    parser.add_option("-s", "--source",
                  action="store", dest="source_file",
                  help="Source input file. It can be a FITS file or "
                  "text file with a list of FITS files.")
    
    parser.add_option("-o", "--output",
                  action="store", dest="output_filename", 
                  help="final coadded output image",
                  default="/tmp/coadd.fits")
    
                                
    (options, args) = parser.parse_args()
    
    if len(sys.argv[1:])<1:
       parser.print_help()
       sys.exit(0)

    # Read the default configuration file
    # If none was specified by the user, environment variable will be used
    if not options.config_file:
        try:
            config_file = os.environ['PAPI_CONFIG']
        except KeyError, error:
            print 'Environment variable PAPI_CONFIG not found!'
            sys.exit()
    else:
        config_file = options.config_file
        
    cfg_options = misc.config.read_config_file(config_file)
    
    # args is the leftover positional arguments after all options have been processed
    if not options.source_file or not options.output_filename or len(args)!=0: 
        parser.print_help()
        parser.error("incorrect number of arguments " )
    
    # Check if source_file is a FITS file or a text file listing a set of files
    if os.path.exists(options.source_file):
        try:
            datahandler.fits_simple_verify(options.source_file)
            filelist = [options.source_file]
        except:
            filelist = [line.replace( "\n", "") 
                      for line in fileinput.input(options.source_file)]
        
        for ifile in filelist:
            try:
                align_stack_cube(ifile, options.output_filename)
                sys.exit(0)
            except Exception,e:
                msg = "Error processing file %s - %s"%(ifile,str(e))
                log.error(msg)
                raise e
            
            