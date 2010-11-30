#! /usr/bin/env python

# Copyright (c) 2009 Jose M. Ibanez. All rights reserved.
# Institute of Astrophysics of Andalusia, IAA-CSIC
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
# PAPI (PAnic PIpeline)
#
# papi.py
#
# Last update 30/Nov/2010
#
################################################################################


    ####################################################################
    #                                                                  #
    #                   This is where it all starts.                   #
    #         The body of this program is at the end of this file      #
    #                                                                  #
    ####################################################################

    
#From system
import sys
import os
import os.path
from optparse import OptionParser
from optparse import IndentedHelpFormatter
import fileinput
import glob
import shutil
import tempfile
import dircache


#Log
import misc.paLog
from misc.paLog import log    

#PAPI packages 
import reduce.reductionset as RS
import misc.config

################################################################################
# main
################################################################################
def main(arguments = None):
    
    if arguments is None:
        arguments = sys.argv[1:] # ignore argv[0], the script name
    
    # Get and check command-line options
    
    description = \
    "This module in the main application for the PANIC data reduction system "

    wider_format = IndentedHelpFormatter(max_help_position = 50, width = 79)
    parser = OptionParser(description = description, 
                                   formatter = wider_format, 
                                   usage = "%prog [OPTION]... DIRECTORY...", 
                                   version = "%prog 1.0")

    parser.add_option("--config", action="store", type="str",
                      dest="config_file",
                      help="if not specified, '%s' is used" \
                            % misc.config.default_config_file())
                  
    parser.add_option("-s", "--source",
                  action="store", dest="source",
                  help="Source file list of data frames. It can be a file or directory name.")
    
    parser.add_option("-o", "--output_file",
                  action="store", dest="output_filename", 
                  help="final reduced output image")
    
    parser.add_option("-d", "--outdir",
                  action="store", dest="out_dir", default="/tmp",
                  help="output dir for intermidiate files")
    
    parser.add_option("-t", "--type",
                  action="store", dest="type", default="quick", 
                  help="type of reduction (quick|science)")
                  
    parser.add_option("-m", "--obs_mode",
                  action="store", dest="obs_mode", 
                  default="dither", help="observing mode (dither|ext_dither)")
    
    parser.add_option("-D", "--dark",
                  action="store", dest="dark",
                  help="master dark to subtract")
    
    parser.add_option("-F", "--flat",
                  action="store", dest="flat",
                  help="master flat to divide by")

    parser.add_option("-b", "--bpm",
                  action="store", dest="bpm", 
                  help="bad pixel mask")
    
    parser.add_option("-C", "--config_file",
                  action="store", dest="config_file", 
                  help="config file for the data reduction process")
                  
    parser.add_option("-1", "--single",
                  action="store_true", dest="single", default=False,
                  help="make a single reduction")              

    parser.add_option("-n", "--no_class",
                  action="store_true", dest="no_class", default=False,
                  help="not try to do dataset classification")

    parser.add_option("-k", "--no_check",
                  action="store_true", dest="no_check", default=False, 
                  help="not check data properties match (type, expt, filter, ncoadd, mjd)")
    
    parser.add_option("-v", "--verbose",
                  action="store_true", dest="verbose", default=True,
                  help="verbose mode [default]")

        
    (init_options, args) = parser.parse_args(args=arguments)
    
    # Read the default configuration file if none was specified by the user
    if not init_options.config_file:
        config_file = misc.config.default_config_file()       
    else:
        config_file = init_options.config_file
    
    options = misc.config.read_options(init_options, 'general', config_file)


    if not args:
        parser.print_help()
        sys.exit(2) 
    
    # Read the options from the configuration file
    options = config.read_config_file()
    general_opts  = options['general'] 
  
        
    if not options.source or not options.output_filename or len(args)!=0: # args is the leftover positional arguments after all options have been processed
        parser.print_help()
        parser.error("incorrect number of arguments " )
    if options.verbose:
        print "reading %s ..." % options.source
    
    
    sci_files=[]
    
    # Read file or list-directory 
    if os.path.isfile(options.source):
        sci_files=[line.replace("\n", "").replace('//','/') for line in fileinput.input(options.source)]
    elif os.path.isdir(options.source):
        for file in dircache.listdir(options.source):
            if file.endswith(".fits") or file.endswith(".fit"):
                sci_files.append((options.source+"/"+file).replace('//','/'))
                
    if options.single: red_m="single"
    else: red_m="full"
    
    log.debug("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    log.debug(">Start PAPI....")
    log.debug("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    
    # Create the RS (it works both simple FITS as MEF files)            
    rs = RS.ReductionSet( sci_files, options.out_dir, out_file=options.output_filename, obs_mode="dither", \
                                dark=options.dark, flat=options.flat, bpm=options.bpm, red_mode=red_m, \
                                classf=not(options.no_class), check_data=not(options.no_check))
    
    try:
        rs.reduceSet(red_mode=red_m)
    except Exception,e:
        print "Cannot reduce the Data Set, check error log...."
        raise e
    else:
        print "Well done (i hope) !!!"
        return 0
    
######################################################################
if __name__ == "__main__":
    sys.exit(main())