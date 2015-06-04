#! /usr/bin/env python
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


import math
import sys
import os
import shutil
import fileinput
from optparse import OptionParser


import numpy
import astropy.io.fits as fits
from astropy import wcs

from misc.paLog import log
from misc.version import __version__

import datahandler


def createDataSeq(input_files, seq_type, overwrite=False, output_dir=None):
    """
    Given a list of FITS files, update header to create  Data Sequence 
    compliant with PAPI, adding/updating the next keywords:
    
    OBS_TOOL,
    PROG_ID,
    OB_ID,
    OB_NAME,
    OB_PAT,
    PAT_NAME,
    PAT_EXPN,
    PAT_NEXP,
    IMAGETYP
    
    Parameters
    ----------
    input_files: list
        list of FITS files
    seq_type: str
        Type of sequece to create (IMAGETYP):
            DARK, SKY_FLAT, DOME_FLAT, FOCUS, SCIENCE
    overwrite: bool
        True if input files are overwritten, otherwise new files are 
        created in 'output_dir'.
    output_dir: str
        Dir name where new files are created when overwrite=False
        
    """
    
    # Check of Seq. type
    if seq_type not in ['DARK', 'SKY_FLAT', 'DOME_FLAT', 'FOCUS', 'SCIENCE']:
        msg = "Error, wrong Sequece type specified"
        log.error(msg)
        raise Exception(msg)
    
    # First of all, check if we need to create  new FITS files or overwrite
    new_files = []
    if not overwrite:
        for ifile in input_files:
            new_file = output_dir + "/" + os.path.basename(ifile).replace(".fits", "_seq.fits")
            shutil.copyfile(ifile, new_file)
            new_files.append(new_file)
    else:
        new_files = input_files
        
    pat_nexp = len(new_files)
    if pat_nexp < 2:
        msg = "Cannot create a sequence, not enought number of files (<2)"
        log.error(msg)
        raise Exception(msg)
    
    prog_id = '       '
    ob_id = '1'
    ob_name = 'OB_dummy'
    ob_pat = 'unknown'
    pat_name = 'unknown'
    pat_expn = 0
    
    # Sort out files by MJD !!
    dataset = []
    sorted_files = []
    for ifile in new_files:
        mjd = fits.getval(ifile, "MJD-OBS", ext=0)
        dataset.append((ifile, mjd))
    
    # Sort out the files
    dataset = sorted(dataset, key=lambda data_file: data_file[1])
    
    # Created sorted list
    for l_tuple in dataset:
        sorted_files.append(l_tuple[0])

    print "SORTED", sorted_files
      
    for ifile in sorted_files:
        pat_expn += 1
        with fits.open(ifile, 'update') as hdu:
            log.debug("Updating file %s"%ifile)
            try:
                if not 'OBS_TOOL' in hdu[0].header:
                    hdu[0].header.set('OBS_TOOL', 'createDS',
                                'PANIC Observing Tool Software version')
                if not 'PROG_ID' in hdu[0].header:
                    hdu[0].header.set('PROG_ID', prog_id,
                                'PANIC Observing Program ID')
                if not 'OB_ID' in hdu[0].header:
                    hdu[0].header.set('OB_ID', ob_id,
                                'PANIC Observing Block ID')
                if not 'OB_NAME' in hdu[0].header:
                    hdu[0].header.set('OB_NAME', ob_name,
                                'PANIC Observing Block Name')
                if not 'OB_PAT' in hdu[0].header:
                    hdu[0].header.set('OB_PAT', ob_pat,
                                'PANIC Observing Block Pattern Type')
                if not 'PAT_NAME' in hdu[0].header:
                    hdu[0].header.set('PAT_NAME', pat_name,
                                'PANIC Observing Sequence Pattern Name')
                # next ones need to be updated in anycase
                hdu[0].header.set('PAT_EXPN', pat_expn,
                                'PANIC Observing Exposition Number')
                hdu[0].header.set('PAT_NEXP', pat_nexp,
                                'PANIC Observing Number of Expositions')
                hdu[0].header.set('IMAGETYP', seq_type,
                                'PANIC Image type')
            except Exception,e:
                msg = "Error updating FITS header: %s" %str(e)
                log.error(msg)
                raise Exception(msg)
    
    log.info("New data sequence of type %s created."%seq_type)
    log.info("Files: %s"%sorted_files)
    
    return sorted_files

    
#################
# MAIN #
#################
if __name__ == "__main__":
    
    usage = "usage: %prog [options]"
    desc = """Modify headers of a set of FITS files to create a Data Sequece compliant with PAPI""" 
    
    parser = OptionParser(usage, description=desc)
    
    parser.add_option("-s", "--source", type="str",
                  action="store", dest="source_file",
                  help="Input text file with the list of FITS file to be analized.")
                  
    parser.add_option("-t", "--seq_type", type="str",
                  action="store", dest="seq_type", default='SCIENCE',
                  help="Tyep of the Data Sequence: DARK, SKY_FLAT, DOME_FLAT, FOCUS, SCIENCE")
    
    parser.add_option("-o", "--output_dir", type="str",
                  action="store", dest="output_dir", default="/tmp",
                  help="Output file to write the new FITS files of the Data Sequence.")
    
    parser.add_option("-O", "--overwrite",
                  action="store_true", dest="overwrite", default=False,
                  help="Overwrite input files [default=%default]")
    
    (options, args) = parser.parse_args()
    
    if len(sys.argv[1:]) < 1:
       parser.print_help()
       sys.exit(0)
        
    # args is the leftover positional arguments after all options have been processed 
    if not options.source_file or len(args)!=0: 
        parser.print_help()
        parser.error("incorrect number of arguments " )
        
    
    if os.path.isfile(options.source_file):
        files = [line.replace("\n", "").replace('//','/')
            for line in fileinput.input(options.source_file)]
    else:
        print "Error, cannot read file : ", options.source_file
        sys.exit(0)
    try:    
        res_files = createDataSeq(files, options.seq_type, options.overwrite, options.output_dir)
    except Exception,e:
        log.error("Error, cannot create Date Sequence: %s",str(e))
    
    
    sys.exit(0)