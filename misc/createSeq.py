#! /usr/bin/env python

# Copyright (c) 2011 JM Ibanez All rights reserved.
# Institute of Astrophysics of Andalusia, IAA-CSIC
#
# This file is part of PAPI.
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


import sys
import fileinput
import pyfits
from optparse import OptionParser
import os.path

# PAPI modules
import style

def create_obs_sequence (filelist, ob_id, ob_pat, suffix='OS', overwrite=False):
    """ Function to modify the FITS headers in order to simulate an 
    observing sequence of PANIC. The keywords to write are:
     
        - OB_ID (given by the user)
        - OB_NAME (not used)
        - OB_PAT (given by the user)
        - PAT_NAME (not used)
        - PAT_EXPN (1...N)
        - PAT_NEXP (N)
        - END_SEQ (T/F)
        
    If overwrite is True, the files are overwriten, otherwise, a copy of files
    is created with the new keywords values    
    """
    
    number_files = len(filelist)
    pat_expn = 1
    pat_nexp = number_files
    
    for file in filelist:
        hdus = pyfits.open(file,"update")
        hdus[0].header.update("OB_ID", ob_id, "Observing Block ID (simulated)")
        hdus[0].header.update("OB_NAME", "OB_DEMO", "Observing Block Name (simulated)")
        hdus[0].header.update("OB_PAT", ob_pat, "OB Pattern used (simulated)")
        hdus[0].header.update("PAT_NAME", "PAT_DEMO", "Pattern name (simulated)")
        hdus[0].header.update("PAT_EXPN", pat_expn, "Exposition number of the current (simulated)")
        hdus[0].header.update("PAT_NEXP", pat_nexp, "Total  number of the expositions of current pattern (simulated)")
        if pat_expn==pat_nexp:
            hdus[0].header.update("END_SEQ", "True", "end of dither sequence (simulated)")
        else:
            hdus[0].header.update("END_SEQ", "False", "end of dither sequence (simulated)")
        
        hdus.close()
        pat_expn = pat_expn + 1
        
        
        
 

################################################################################
# main
if __name__ == "__main__":
    
    # Get and check command-line options
    usage = "usage: %prog [options] arg1 arg2"
    parser = OptionParser (usage)
    
    
    parser.add_option ("-i", "--input",
                  action = "store", dest = "input_file_list", type="str",
                  help = "Source file list of data frames to modify \
                  It has to be a fullpath file name")
    
    parser.add_option ("-s", "--suffix", type="str",
                  action = "store", dest = "out_suffix", \
                  help = "suffix added to new out files (default .%02d.fits)")
    
    parser.add_option ("-b", "--ob_id", type="str",
                  action = "store", dest = "ob_id", \
                  help = "Observing Block ID of sequence" , 
                  default = False)
                                 
    parser.add_option ("-p", "--ob_pat", type="str",
                  action = "store", dest = "ob_pat", \
                  help = "OB pattern of sequence")
    
    parser.add_option("-o", "--overwrite input files",
                  action = "store_true", dest = "overwrite", default = False,
                  help = "overwrite input files (false)")
    
    parser.add_option("-v", "--verbose",
                  action = "store_true", dest = "verbose", default = True,
                  help = "verbose mode [default]")
                  
    
    (options, args) = parser.parse_args()
    
    if options.input_file_list:
        filelist = [line.replace( "\n", "") for line in fileinput.input(options.input_file_list)]
        print filelist
    else:
        parser.print_help()
        parser.error("incorrect number of arguments " )
        
        
    create_obs_sequence(filelist, options.ob_id, options.ob_pat, 
                        options.out_suffix, options.overwrite)
