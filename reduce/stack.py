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
# stack.py
#
# Created    : 20/11/2012    jmiguel@iaa.es -
# Last update: 
# TODO
################################################################################
"""
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


""""  TO BE COMPLETED !!!! """

def do_stack(in_filelist, out_image=None, s_type='median', run_offsets=True, 
             overwrite=False, config_file=None):
    """
    Do a stacking of images computing the offsets, doing the alignment and
    finally the combination using the mean or median of the values.
     
    Notes
    -----
    It is supposed that image background is already subtracted from images
    
    Parameters
    ----------
    in_filelist : list
        input list of images to be stacked
    
    out_image : str
        Output filename of resulted image
    
    s_type: str
        Type of image combination to be done (mean, median)
        
    run_offsets: boolean
        Whether true, the offsets between the images is computed and then
        image alignment is done previous combination
        
    overwrite: boolean
        If true, the 'out_image' filename will be overwritten if it exists,
        otherwise, the 'out_image' will be created 
    
    Returns
    -------
    If all was successful, the name of the output file is returned
        
    """
    
    # read the input file
    if in_filelist and os.path.exists(in_filelist):
        filelist = [line.replace( "\n", "") for line in fileinput.input(in_filelist)]
    else:
        log.error("Can not read input file list")
        raise Exception("Can not read input file list")
    
    if config_file and os.path.exists(config_file):
        conf = misc.config.read_config_file(config_file)
    
    out_dir = os.path.abspath(os.path.join(out_image, os.pardir))  
    
    task = RS.ReductionSet( filelist, 
                            out_dir,
                            out_file=os.path.basename(out_image),
                            obs_mode="dither", 
                            dark=None, 
                            flat=None,
                            bpm=None, 
                            red_mode="quick",
                            group_by='filter', 
                            check_data=True, 
                            config_dict=conf)
    
            
    
    # Coaddition using offsets
    if run_offsets:
        log.info("**** Coaddition of sky subtracted frames ****")
        try:
            offset_mat = task.getPointingOffsets(out_dir+"/files_skysub.list", 
                                                     out_dir+'/offsets1.pap')                
        except Exception,e:
            log.error("Erron while getting pointing offsets. Cannot continue with data reduction...")
            raise e
        
    # Coaddition using offsets
    log.info("**** Coaddition of sky subtracted frames ****")
    fo = open(out_dir + '/offsets1.pap',"r")
    fs = open(out_dir + '/stack1.pap','w+')
    for line in fo:
        n_line = line.replace(".fits.objs", ".fits") 
        fs.write( n_line )
    fo.close()
    fs.close()
    gainmap = '???'  # TODO    
    task.coaddStackImages(out_dir+'/stack1.pap', gainmap, 
                          out_dir+'/coadd1.fits','average')
    
    
# main
if __name__ == "__main__":
    
    usage = "usage: %prog [options]"
    desc = """Do a stacking of images computing the offsets, doing the alignment 
and finally the combination using the mean or median of the values.

TO BE COMPLETED !!!!

"""
    parser = OptionParser(usage, description=desc)
    
    parser.add_option("-i", "--input_file",
                  action="store", dest="input_file", 
                  help="Input file list")
                  
    parser.add_option("-o", "--output",
                  action="store", dest="output_image", 
                  help="output filename (default = %default)",
                  default="stack.fits")
    
    parser.add_option("-c", "--config_file",
                  action="store", dest="config_file", 
                  help="Configuration file required")
    
    parser.add_option("-f", "--run_offsets",
                  action="store_true", dest="run_offsets", default=False,
                  help="Compute image offsets and alignment")
    
    parser.add_option("-O", "--overwrite",
                  action="store_true", dest="overwrite", default=False,
                  help="overwrite the output image if exits")
                                
    (options, args) = parser.parse_args()
    
    if len(sys.argv[1:])<1:
       parser.print_help()
       sys.exit(0)

    if not options.input_file or len(args)!=0: 
    # args is the leftover positional arguments after all options have been processed
        parser.print_help()
        parser.error("wrong number of arguments " )
    if not options.output_image:
        options.output_image = None

    try:    
        do_stack(options.input_file, options.output_image, 
                         options.run_offsets, options.overwrite, options.config_file)
    except Exception, e:
        log.error("Fail of stacking !")
        raise e
    print "\nWell done !"
    
