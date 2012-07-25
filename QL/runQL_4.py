#! /usr/bin/env python

# Copyright (c) 2010 Jose M. Ibanez. All rights reserved.
# Institute of Astrophysics of Andalusia, IAA-CSIC
#
# This file is part of PAPI QL (PANIC Quick Look)
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
# runQL (run the QL for PANIC pipeline)
#
# runQL.py
#
# Last update 09/Dic/2010
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
#import misc.paLog
from misc.paLog import log    

#PAPI packages
from PyQt4 import QtCore, QtGui

import mainGUI_4
import misc.config

################################################################################
# main
################################################################################
def main(arguments = None):
    
    if arguments is None:
        arguments = sys.argv[1:] # ignore argv[0], the script name
    
    # Get and check command-line options
    
    description = \
    "This module in the main application for the PANIC Quick Loook (PQL) data reduction system "

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
    
    parser.add_option("-v", "--verbose",
                  action="store_true", dest="verbose", default=True,
                  help="verbose mode [default]")
    
    parser.add_option("-s", "--source", type="str",
                  action="store", dest="source",
                  help="Source directory of data frames. It has to be a fullpath file name")
    
    parser.add_option("-o", "--output_dir", type="str",
                  action="store", dest="output_dir", help="output directory to write products")
    
    parser.add_option("-t", "--temp_dir", type="str",
                  action="store", dest="temp_dir", help="temporary directory to write")
    
    # parser.add_option("-c", "--config",
    #                action="store", dest="config_file", help="Quick Look config file")
              
        
    
    (init_options, i_args) = parser.parse_args(args=arguments)
    
    # Read the default configuration file if none was specified by the user
    if not init_options.config_file:
        config_file = misc.config.default_config_file()
    else:
        config_file = init_options.config_file
    
    # now, we "mix" the invokation parameter values with the values in the config file,
    # having priority the invokation values over config file values
    # note: only values of the 'quicklook' section can be invoked
    options = misc.config.read_options(init_options, 'quicklook', config_file)

    #print "options=",options

    if  i_args: # i_args is the leftover positional arguments after all options have been processed
        parser.print_help()
        sys.exit(2) 
    
    ql_opts  = options['quicklook'] 
    
    if not ql_opts['source'] or not ql_opts['output_dir'] or not ql_opts['temp_dir']: 
        parser.print_help()
        parser.error("Incorrect number of arguments " )
    
    log.debug("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    log.debug(">> Starting Quick-Look....")
    log.debug("   + source  : %s",ql_opts['source'])
    log.debug("   + out_dir : %s",ql_opts['output_dir'])
    log.debug("   + temp_dir : %s",ql_opts['temp_dir'])
    log.debug("   + run_mode: %s",ql_opts['run_mode'])
    log.debug("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    
    app = QtGui.QApplication(sys.argv)
    try:
        f = mainGUI_4.MainGUI(ql_opts['source'], ql_opts['output_dir'], 
                            ql_opts['temp_dir'], config_opts=options)
        f.show()
        app.setMainWidget(f)
        app.exec_loop()
    except:
        log.debug("Some error while running mainGUI")
        raise

######################################################################
if __name__ == "__main__":
    sys.exit(main())
