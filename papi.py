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
import dircache


#Log
#import misc.paLog
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

    #wider_format = IndentedHelpFormatter(max_help_position = 50, width = 79)
    parser = OptionParser(description = description, 
                                   #formatter = wider_format, 
                                   usage = "%prog [OPTION]... DIRECTORY...", 
                                   version = "%prog 1.0")
    # general options

    parser.add_option("-c", "--config", type="str",
                      action="store", dest="config_file",
                      help="config file for the PANIC Pipeline application. If not specified, '%s' is used" \
                      % misc.config.default_config_file())
                  
    parser.add_option("-s", "--source", type="str",
                  action="store", dest="source",
                  help="Source file list of data frames. It can be a file or directory name.")
    
    parser.add_option("-o", "--output_file", type="str",
                  action="store", dest="output_filename", 
                  help="final reduced output image")
    
    parser.add_option("-t", "--temp_dir", type="str",
                  action="store", dest="temp_dir",
                  help="directory for temporal files")              
    
    parser.add_option("-d", "--out_dir", type="str",
                  action="store", dest="output_dir", 
                  help="output dir for product files")
    
    parser.add_option("-r", "--red_mode", type="str",
                  action="store", dest="reduction_mode", 
                  help="Mode of data reduction to do (quick|science)")
                  
    parser.add_option("-m", "--obs_mode", type="str",
                  action="store", dest="obs_mode", 
                  default="dither", help="observing mode (dither|ext_dither|other)")
    
    parser.add_option("-v", "--verbose",
                  action="store_true", dest="verbose", default=True,
                  help="verbose mode [default]")
                  
    # file calibration options
    
    parser.add_option("-D", "--master_dark",
                  action="store", dest="master_dark",
                  help="master dark to subtract")
    
    parser.add_option("-F", "--master_flat", type="str",
                  action="store", dest="master_flat",
                  help="master flat to divide by")

    parser.add_option("-b", "--bpm_file", type="str",
                  action="store", dest="bpm_file", 
                  help="bad pixel mask file")

    parser.add_option("-g", "--group_by", type="str",
                  action="store", dest="group_by",
                  help="kind of data grouping (based on) to do with the dataset files (ot |filter)")

    parser.add_option("-k", "--check_data", 
                  action="store_true", dest="check_data", 
                  help="if true, check data properties matching (type, expt, filter, ncoadd, mjd)")

        
    (init_options, i_args) = parser.parse_args(args=arguments)
    
    # Read the default configuration file if none was specified by the user
    if not init_options.config_file:
        config_file = misc.config.default_config_file()       
    else:
        config_file = init_options.config_file
    
    # now, we "mix" the invokation parameter values with the values in the config file,
    # having priority the invokation values over config file values
    # note: only values of the 'general' section can be invoked
    options = misc.config.read_options(init_options, 'general', config_file)

    print "options=",options

    if  i_args:    # i_args is the leftover positional arguments after all options have been processed
        parser.print_help()
        sys.exit(2) 
    
    general_opts  = options['general'] 
    print "GEN_OPTS",general_opts
  
        
    if not general_opts['source'] or not general_opts['output_file'] \
        or not general_opts['output_dir'] or not general_opts['temp_dir']   \
        or len(i_args)!=0: # args is the leftover positional arguments after all options have been processed
        parser.print_help()
        parser.error("incorrect number of arguments " )
    if general_opts['verbose']:
        print "reading %s ..." % general_opts['source']
    
    
    sci_files=[]
    
    # Read file or list-directory 
    if os.path.isfile(general_opts['source']):
        sci_files=[line.replace("\n", "").replace('//','/') for line in fileinput.input(general_opts['source'])]
    elif os.path.isdir(general_opts['source']):
        for file in dircache.listdir(general_opts['source']):
            if file.endswith(".fits") or file.endswith(".fit"):
                sci_files.append((general_opts['source']+"/"+file).replace('//','/'))
    
    log.debug("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    log.debug(">> Starting PAPI....")
    log.debug("   + source  : %s",general_opts['source'])
    log.debug("   + out_dir : %s",general_opts['output_dir'])
    log.debug("   + reduction_mode: %s",general_opts['reduction_mode'])
    log.debug("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    
    # Create the RS (it works both simple FITS as MEF files)            
    rs = RS.ReductionSet( sci_files, general_opts['output_dir'], out_file=general_opts['output_file'], \
                          obs_mode=general_opts['obs_mode'], \
                          dark=general_opts['master_dark'], flat=general_opts['master_flat'], \
                          bpm=general_opts['master_bpm'], red_mode=general_opts['reduction_mode'], \
                          group_by=general_opts['group_by'], check_data=general_opts['check_data'], \
                          config_dict=options \
                        )
    
    try:
        rs.reduceSet(red_mode=general_opts['reduction_mode'])
    except Exception,e:
        print "Cannot reduce the Data Set, check error log...."
        raise e
    else:
        print "Well done (i hope) !!!"
        return 0
    
######################################################################
if __name__ == "__main__":
    sys.exit(main())
    
    