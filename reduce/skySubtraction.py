#! /usr/bin/env python
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


################################################################################
#
#
# PAPI (PANIC PIpeline)
#
# skySubtration.py
#
# Created    : 20/11/2012    jmiguel@iaa.es -
# Last update: 
# TODO
################################################################################
"""
TO BE DONE !!!! - NOT IMPLEMENTED YET !!!
"""

# Import necessary modules
from optparse import OptionParser
import os
import fileinput
import sys

# Logging
from misc.paLog import log
import reduce.reductionset as RS
import misc.config


def do_skySubtration(in_filelist, out_dir=None, n_nearest=1, 
             overwrite=False, config_file=None):
    """
    Given a set of images, the sky subtraction is done to each one.
    The sky background to be subtracted is computed with the (2*N)-nearest 
    frames around each frame (N>1). 
    
    Parameters
    ----------
    in_filelist : list
        input list of images to be stacked
    
    out_dir : str
        Output directory for sky subtracted images
    
    n_nearest: int
        Number of (2*n)-nearest frames to be used for the sky background 
        computation. (n>1)
        
    overwrite: boolean
        If true, the 'out_image' filename will be overwritten if it exists,
        otherwise, the 'out_image' will be created 
    
    Returns
    -------
    If all was successful, the name of the output file is returned
        
    """
    
    # read the input file
    if os.path.exists(in_filelist):
        filelist = [line.replace( "\n", "") for line in fileinput.input(in_filelist)]
    else:
        log.error("Can not read input file list")
        raise Exception("Can not read input file list")
    
    if out_dir==None:
        out_dir = os.path.abspath(os.path.join(in_filelist, os.pardir))

    if os.path.exists(config_file):
        conf = misc.config.read_config_file(config_file)
        
    try:    
        task = RS.ReductionSet( filelist, 
                                out_dir,
                                out_file=os.path.basename(in_filelist) + ".fits",
                                obs_mode="dither", 
                                dark=None, 
                                flat=None,
                                bpm=None, 
                                red_mode="quick",
                                group_by='filter', 
                                check_data=True, 
                                config_dict=conf)
                
        task.subtractNearSky()
    except Exception, e:
        log.error("Error in ReductionSet")
        raise e
    
# main
if __name__ == "__main__":
    
    usage = "usage: %prog [options]"
    desc = """Subtract sky-background from the input image

TO BE COMPLETED !!!!
"""
    parser = OptionParser(usage, description=desc)
    
    parser.add_option("-i", "--input_file",
                  action="store", dest="input_file", 
                  help="Input file list")
                  
    parser.add_option("-o", "--output_dir",
                  action="store", dest="output_dir", 
                  help="Output filename (default = CWD)")
    
    parser.add_option("-c", "--config_file",
                  action="store", dest="config_file", 
                  help="Configuration file required")
    
    parser.add_option("-n", "--n_nearest",
                  action="store_true", dest="n_nearest", default=1,
                  help="Number of (2*n)-nearest frames for sky computation")
    
    parser.add_option("-O", "--overwrite",
                  action="store_true", dest="overwrite", default=False,
                  help="Overwrite the output image if it exits")
                                
    (options, args) = parser.parse_args()
    
    if len(sys.argv[1:])<1:
       parser.print_help()
       sys.exit(0)

    if not options.input_file or len(args)!=0: 
    # args is the leftover positional arguments after all options have been processed
        parser.print_help()
        parser.error("wrong number of arguments " )
    if not options.output_dir:
        options.output_image = os.path.abspath(os.path.join(options.input_file, os.pardir))

    try:
        do_skySubtration(options.input_file, options.output_dir, 
                         options.n_nearest, options.overwrite, options.config_file)
    except Exception, e:
        log.error("Fail of sky subtraction !")
        raise e
    print "\nWell done !"
    
