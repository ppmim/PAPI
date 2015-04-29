#! /usr/bin/env python
#
# Copyright (c) 2008-2015 IAA-CSIC  - All rights reserved. 
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
#

################################################################################
#
# PAPI (PAnic PIpeline)
#
# papi.py
#
# Last update 17/May/2011
#
################################################################################

'''User command line interface of PAPI.'''

    ####################################################################
    #                                                                  #
    #                   This is where it all starts.                   #
    #         The body of this program is at the end of this file      #
    #                                                                  #
    ####################################################################

    
#From system
import sys
#import os
import os.path
from optparse import OptionParser
#from optparse import IndentedHelpFormatter
import fileinput
import dircache


#Log
#import misc.paLog
from misc.paLog import log    

#PAPI packages 
import reduce.reductionset as RS
import misc.config
import misc.utils as utils
import misc.genLogsheet as gls
import misc.check_papi_modules
from misc.version import __version__


################################################################################
# main
################################################################################
def main(arguments = None):
    
    # Clock to measure the CPU time used
    t = utils.clock()
    t.tic()
        
    if arguments is None:
        arguments = sys.argv[1:] # ignore argv[0], the script name
    
    # Get and check command-line options
    
    description = \
    "This is PAPI, the PANIC PIpeline data reduction system "\
    "- IAA-CSIC - Version %s"%__version__

    parser = OptionParser(description = description, 
                          usage = "%prog [OPTION]... DIRECTORY...", 
                          version = "%%prog version %s"%__version__)
    # general options

    parser.add_option("-c", "--config", type = "str",
                      action = "store", dest = "config_file",
                      help = "Config file for the PANIC Pipeline application."
                       "If not specified, '%s' is used."
                       % misc.config.default_config_file())

    parser.add_option("-s", "--source", type = "str",
                  action = "store", dest = "source",
                  help = "Source file list of data frames. It can be a file" 
                  "or directory name.")
    
    parser.add_option("-d", "--out_dir", type = "str",
                  action = "store", dest = "output_dir", 
                  help = "Output dir for product files")
    
    parser.add_option("-o", "--output_file", type = "str",
                  action = "store", dest = "output_file", 
                  help = "Final reduced output image")
    
    parser.add_option("-t", "--temp_dir", type = "str",
                  action = "store", dest = "temp_dir",
                  help = "Directory for temporal files")              
   
    parser.add_option("-r", "--rows", nargs=2,
                  action="store", dest="rows", type=int,
                  help="Use _only_ files of the source file-list in the range" 
                  "of rows specified (0 to N, both included)")
    
    parser.add_option("-R", "--recursive",
                  action="store_true", dest = "recursive", default = False,
                  help="Does recursive search for files in source directory")
    
    parser.add_option("-l", "--list",
                  action="store_true", dest="list", default = False,
                  help="Generate a list with all the source files read from" 
                  "the source and sorted by MJD")

    parser.add_option("-M", "--red_mode", type = "str",
                  action = "store", dest = "reduction_mode", 
                  help = "Mode of data reduction to do "
                  "(quick|science|lab|lemon|quick-lemon).")
                  
    parser.add_option("-m", "--obs_mode", type = "str",
                  action = "store", dest = "obs_mode", 
                  default = "dither", 
                  help = "Observing mode (dither|ext_dither|other)")
    
    
    parser.add_option("-S", "--seq_to_reduce", type=int, nargs=2,
                  action = "store", dest = "seq_to_reduce",
                  default = -1, help = "Sequence number to reduce. By default, " 
                  "all sequences found will be reduced.")

    parser.add_option('-W', '--window_detector',
                      type='choice',
                      action='store',
                      dest='detector',
                      choices=['Q1', 'Q2', 'Q3', 'Q4', 'Q123', 'all'],
                      default='all',
                      help="Specify which detector to process:"
                      "Q1(SG1), Q2(SG2), Q3(SG3), Q4(SG4), Q123(all except SG4), all [default: %default]")
    
    parser.add_option("-p", "--print",
                  action = "store_true", dest = "print_seq", default = False,
                  help = "Print all detected sequences in the Data Set")
    
    parser.add_option('-T', '--sequences_type',
                      type='choice',
                      action='store',
                      dest='seq_type',
                      choices=['DARK', 'FLAT', 'DOME_FLAT', 'SKY_FLAT', 'FOCUS', 'SCIENCE', 'CAL', 'all'],
                      default='all',
                      help="Specify the type of sequences to show: "
                      "DARK, FLAT(all), DOME_FLAT, SKY_FLAT, FOCUS, SCIENCE, CAL, all [default: %default]")
    
    parser.add_option("-b", "--build_calibrations",
                  action = "store_true", dest = "build_calibrations", default = False,
                  help = "Build all the master calibrations files")
    
    # file calibration options
    
    parser.add_option("-C", "--ext_calibration_db", type = "str",
                  action = "store", dest = "ext_calibration_db", 
                  default = None, 
                  help = "External calibration directory "
                  "(library of Dark & Flat calibrations)")
    
    parser.add_option("-D", "--master_dark",
                  action = "store", dest = "master_dark",
                  help = "Master dark to subtract")
    
    parser.add_option("-F", "--master_flat", type="str",
                  action = "store", dest = "master_flat",
                  help = "Master flat to divide by")

    parser.add_option("-B", "--bpm_file", type = "str",
                  action = "store", dest = "bpm_file", 
                  help = "Bad pixel mask file")

    parser.add_option("-g", "--group_by", type = "str",
                  action = "store", dest = "group_by",
                  help = "kind of data grouping (based on) to do with the" 
                  "dataset files (ot |filter)")

    parser.add_option("-k", "--check_data", 
                  action = "store_true", dest = "check_data", 
                  help = "if true, check data properties matching (type, expt, "
                  "filter, ncoadd, mjd)")
    
    parser.add_option("-e", "--Check",
                      action = "store_true", dest = "check_modules",
                      help = "Check if versions of PAPI modules are right.",
                       default = False)
    
        
    (init_options, i_args) = parser.parse_args(args = arguments)
    
    # If no arguments, print help
    if len(arguments) < 1:
       parser.print_help()
       sys.exit(0)

    # Check required modules and versions
    if init_options.check_modules:
        misc.check_papi_modules.check_modules()
        return
    
    # Read the default configuration file
    # If none was specified by the user, environment variable will be used
    if not init_options.config_file:
        try:
            config_file = os.environ['PAPI_CONFIG']
        except KeyError, error:
            print 'Environment variable PAPI_CONFIG not found!'
            sys.exit()
    else:
        config_file = init_options.config_file
    
    # Now, we "mix" the invokation parameter values with the values in the 
    # config file, having priority the invokation values over config file values
    # Note: only values of the 'general' section can be invoked
    options = misc.config.read_options(init_options, 'general', config_file)

    # Manually mix bpm_file, having priority the invokation values over config
    # file value.
    if init_options.bpm_file != None:
      options['bpm']['bpm_file'] = init_options.bpm_file
    
    # Add the configuration filename as an extra value to the dictionary
    options['general']['config_filename'] = config_file
    
    #print "options = ",options
    
    # 'i_args' is the leftover positional arguments after all options have been 
    # processed.
    if i_args: 
        parser.print_help()
        sys.exit(2) 
    
    general_opts  = options['general'] 
    #print "GEN_OPTS",general_opts
  
        
    if not general_opts['source'] or not general_opts['output_file'] \
        or not general_opts['output_dir'] or not general_opts['temp_dir']   \
        or len(i_args) != 0: # args is the leftover positional arguments after all options have been processed
        parser.error("Incorrect number of arguments " )
        parser.print_help()
    
    rs_files = []
    
    # Read file or list-directory 
    if os.path.isfile(general_opts['source']):
        rs_files = [line.replace("\n", "").replace('//','/') 
                     for line in fileinput.input(general_opts['source'])]
    elif os.path.isdir(general_opts['source']):
        if init_options.recursive:
            for root, dirs, files in os.walk(general_opts['source'], topdown=True):
                for name in files:
                    if name.endswith(".fits") or name.endswith(".fit"):
                        rs_files.append(os.path.join(root, name))
        else:
            for file in dircache.listdir(general_opts['source']):
                if file.endswith(".fits") or file.endswith(".fit"):
                    rs_files.append((general_opts['source']+"/"+file).replace('//','/'))
    
    # Check for list files sorted by MJD
    if init_options.list == True:
        log.debug("Creating logsheet file ....")
        logSheet = gls.LogSheet(rs_files, "/tmp/logsheet.txt", [0,len(rs_files)], True)
        logSheet.create()
        return
   
    # Take only the rows(files) required
    if (os.path.isfile(general_opts['source']) and init_options.rows!=None):
        if (init_options.rows[0]<0) or (init_options.rows[1]>len(rs_files)-1):
            parser.error("wrong rows index values (0,%s)"%(len(rs_files)-1))
            parser.print_help()
        i = 0
        tmp_sci_files = []
        for file in rs_files:
            if i>=init_options.rows[0] and i<=init_options.rows[1]:
                tmp_sci_files.append(file)
            i = i+1
        rs_files = tmp_sci_files 
    
    
    log.debug("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    log.debug(">> Starting PAPI....")
    log.debug("   + source  : %s",general_opts['source'])
    log.debug("   + out_dir : %s",general_opts['output_dir'])
    log.debug("   + Master Dark : %s",general_opts['master_dark'])
    log.debug("   + Master Flat : %s",general_opts['master_flat'])
    log.debug("   + Master BPM  : %s",options['bpm']['bpm_file'])
    log.debug("   + Calibration library: %s",general_opts['ext_calibration_db'])
    log.debug("   + Reduction_mode: %s",general_opts['reduction_mode'])
    log.debug("   + Astrometric catalog: %s", options['astrometry']['catalog'])
    log.debug("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    
    # Create the RS (it works both simple FITS as MEF files)            
    try:
        rs = RS.ReductionSet( rs_files, general_opts['output_dir'], 
                out_file=general_opts['output_file'],
                obs_mode=general_opts['obs_mode'],
                dark=general_opts['master_dark'], 
                flat=general_opts['master_flat'],
                bpm=general_opts['master_bpm'], 
                red_mode=general_opts['reduction_mode'],
                group_by=general_opts['group_by'], 
                check_data=general_opts['check_data'],
                config_dict = options )
    
                    
        if init_options.print_seq:
            rs.getSequences(show=True, stype=init_options.seq_type)
        elif init_options.build_calibrations:
            rs.buildCalibrations()
        else:
            if init_options.seq_type == 'CAL':
                stype = ['DARK', 'DOME_FLAT', 'SKY_FLAT']
            elif init_options.seq_type == 'FLAT':
                stype = ['DOME_FLAT', 'SKY_FLAT']
            else:
                stype = [init_options.seq_type]
            if init_options.seq_to_reduce == -1: #all
                rs.reduceSet(red_mode=general_opts['reduction_mode'],
                             types_to_reduce=stype)
            else:
                m_seqs_to_reduce = []
                m_seqs_to_reduce = range(init_options.seq_to_reduce[0], 
                                         init_options.seq_to_reduce[1]+1)
                #for elem in init_options.seq_to_reduce: 
                #    m_seqs_to_reduce.append(int(elem))
                #return
                rs.reduceSet(red_mode=general_opts['reduction_mode'], 
                             seqs_to_reduce=m_seqs_to_reduce,
                             types_to_reduce=stype)
                
    except RS.ReductionSetException, e:
        log.error("Error during data reduction: %s " % str(e))
    
    # The following handles ctrl-c. We need to do it this way due to a
    # bug in multiprocessing, see:
    # http://stackoverflow.com/questions/1408356/keyboard-interrupts-with-pythons-multiprocessing-pool
    # and http://bugs.python.org/issue8296.
    except (KeyboardInterrupt, SystemExit):
        log.error("Ctrl-c, KeyboardInterrupt !")
        sys.exit()
    except Exception, e:
        print "Cannot reduce the Data Set, check error log...."
        print str(e)
    else:
        print "\n\nWell done (I hope) -  %s!!!"%t.tac()
        return 0
    
######################################################################
if __name__ == "__main__":
    sys.exit(main())
    
    
