#! /usr/bin/env python
#encoding:UTF-8

# Copyright (c) 2010 Jose M. Ibanez All rights reserved.
# Institute of Astrophysics of Andalusia, IAA-CSIC
#
# This file is part of PAPI
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

# ======================================================================
#
# SWARP module.
#
# Author: Jose Miguel Ibanez Mengual <jmiguel@iaa.es>
#
# ======================================================================

"""
A wrapper for SWARP

A wrapper for SWARP from Astromatic (E.Bertin).
by Jose M. Ibanez Mengual
version: 0.1 - last modified: 2010-06-09

This wrapper allows you to configure SWARP, run it and get
back its outputs without the need of editing SWARP
configuration files. by default, configuration files are created
on-the-fly, and SWARP is run silently via python.

Tested on SWARP versions 2.17.x


Example of use:

-----------------------------------------------------------------

TODO   

-----------------------------------------------------------------


"""

# ======================================================================

import __builtin__

import sys
import os
import popen2
import exceptions
import re
import copy
import fileinput
import types


#PAPI 
import misc.utils

# ======================================================================

__version__ = "0.1 (2010-06-09)"

# ======================================================================

class SWARPException(Exception):
    pass

# ======================================================================

class SWARP:
    """
    A wrapper class to transparently use SWARP.

    """

    _SW_config = { 

#----------------------------- Output--------- ---------------------------------
    
        "IMAGEOUT_NAME":
        {"comment": "Output filename",
         "value": "coadd.fits"},
        
        "WEIGHTOUT_NAME":
        {"comment": "Output weight-map filename",
         "value": "coadd.weight.fits"},
        
        "HEADER_ONLY":
        {"comment": "Only a header as an output file (Y/N)?",
         "value": "N"},
        
        "HEADER_SUFFIX":
        {"comment": 'Filename extension for additional headers',
         "value": ".head"},
         
#------------------------------- Input Weights --------------------------------        
        
        "WEIGHT_TYPE":
        {"comment": "BACKGROUND,MAP_RMS,MAP_VARIANCE or MAP_WEIGHT",
         "value": "MAP_WEIGHT"},
        
        "WEIGHT_SUFFIX":
        {"comment": "Suffix to use for weight-maps",
         "value": ".weight.fits"},
        
        "WEIGHT_IMAGE":
        {"comment": "Weightmap filename if suffix not used (all or for each weight-map)",
         "value": ""},
        
        "WEIGHT_THRESH":
        {"comment": "Bad pixel weight-threshold",
         "value": ""},

#------------------------------- Co-addition ----------------------------------
        
        "COMBINE":
        {"comment": 'Combine resampled images (Y/N)?',
         "value": 'Y'},
        
        "COMBINE_TYPE":
        {"comment": "MEDIAN,AVERAGE,MIN,MAX,WEIGHTED,CHI2 or SUM",
         "value": "MEDIAN"},
        
        "BLANK_BADPIXELS":
        {"comment": "Set to 0 pixels having a weight of 0",
         "value": "N"},

#-------------------------------- Astrometry ----------------------------------

        "CELESTIAL_TYPE":
        {"comment": "NATIVE, PIXEL, EQUATORIAL,GALACTIC,ECLIPTIC, or SUPERGALACTIC",
         "value": "NATIVE"},

        "PROJECTION_TYPE":
        {"comment": "Any WCS projection code or NONE",
         "value": "TAN"},

        "PROJECTION_ERR":
        {"comment": "Maximum projection error (in output pixels), or 0 for no approximation",
         "value": 0.001 },

        "CENTER_TYPE":
        {"comment": "MANUAL, ALL or MOST",
         "value": "ALL"},

        "CENTER":
        {"comment": "Coordinates of the image center",
         "value": "00:00:00.0, +00:00:00.0"},

        "PIXELSCALE_TYPE":
        {"comment": "MANUAL,FIT,MIN,MAX or MEDIAN",
         "value": "MEDIAN"},
        
        "PIXEL_SCALE":
        {"comment": "Pixel scale",
         "value": 0},

        "IMAGE_SIZE":
        {"comment": "Image size (0 = AUTOMATIC)",
         "value": "N"},

#-------------------------------- Resampling ----------------------------------
       
        "RESAMPLE":
        {"comment": "Resample input images (Y/N)?",
         "value": "Y"},

        "RESAMPLE_DIR":
        {"comment": "Directory path for resampled images",
         "value": "."},

        "RESAMPLE_SUFFIX":
        {"comment": "filename extension for resampled images",
         "value": ".resamp.fits"},

        "RESAMPLING_TYPE":
        {"comment": "NEAREST,BILINEAR,LANCZOS2,LANCZOS3or LANCZOS4 (1 per axis)",
         "value": "LANCZOS3"},

        "OVERSAMPLING":
        {"comment": "Oversampling in each dimension 0 = automatic)",
         "value": 0},

        "INTERPOLATE":
        {"comment": "Interpolate bad input pixels (Y/N)? (all or for each image)",
         "value": "N"},

        "FSCALASTRO_TYPE":
        {"comment": "NONE,FIXED, or VARIABLE",
         "value": "FIXED"},

        "FSCALE_KEYWORD":
        {"comment": "FITS keyword for the multiplicative factor applied to each input image",
         "value": "FLXSCALE"},

        "FSCALE_DEFAULT":
        {"comment": "Default FSCALE value if not in header",
         "value": 1.0},

        "GAIN_KEYWORD":
        {"comment": "FITS keyword for effect. gain (e-/ADU)",
         "value": "GAIN"},

        "GAIN_DEFAULT":
        {"comment": "Default gain if no FITS keyword found 0 = infinity (all or for each image)",
         "value": 0.0},

        "SATLEV_KEYWORD":
        {"comment": "FITS keyword for saturation level (ADU)",
         "value": "SATURATE"},

        "SATLEV_DEFAULT":
        {"comment": "Default saturation if no FITS keyword",
         "value": 50000.0},

#--------------------------- Background subtraction ---------------------------

        "SUBTRACT_BACK":
        {"comment": "Subtraction sky background (Y/N)? (all or for each image)",
         "value": "Y"},

        "BACK_TYPE":
        {"comment": "AUTO or MANUAL (all or for each image)",
         "value": "AUTO"},

        "BACK_DEFAULT":
        {"comment": "Default background value in MANUAL (all or for each image)",
         "value": 0.0},

        "BACK_SIZE":
        {"comment": "Background mesh size (pixels) all or for each image)",
         "value": 128},
         
        "BACK_FILTERSIZE":
        {"comment": "Background map filter range (meshes) (all or for each image)",
         "value": 3},

        "BACK_FILTTHRESH":
        {"comment": "Threshold above which the background-map filter operates",
         "value": 0.0},
         
#------------------------------ Memory management -----------------------------

        "VMEM_DIR":
        {"comment": "Directory path for swap files",
         "value": "."},

        "VMEM_MAX":
        {"comment": "Maximum amount of virtual memory (MB)",
         "value": 2047},
         
        "MEM_MAX":
        {"comment": "Maximum amount of usable RAM (MB)",
         "value": 128},

        "COMBINE_BUFSIZE":
        {"comment": "RAM dedicated to co-addition(MB)",
         "value": 64},
        
#------------------------------ Miscellaneous ---------------------------------
        
        "DELETE_TMPFILES":
        {"comment": "Delete temporary resampled FITS files (Y/N)?",
         "value": "Y"},
        
        "COPY_KEYWORDS":
        {"comment": 'List of FITS keywords to propagate from the input to the output headers',
         "value": "OBJECT,RA,DEC,OBJECT,PIXSCALE,ROT-RTA,INSTRUME,TELESCOPE"},
        
        "WRITE_FILEINFO":
        {"comment": "Write information about each input file in the output image header?",
         "value": "N"},
        
        "WRITE_XML":
        {"comment": 'Write XML file (Y/N)?',
         "value": "N"},
        
        "XML_NAME":
        {"comment": "Filename for XML output",
         "value": "swarp.xml"},
        
        "XSL_URL":
        {"comment": "Filename for XSL style-sheet",
         "value": "file:///usr/local/Terapix//share/swarp/swarp.xsl"},
        
        "VERBOSE_TYPE":
        {"comment": "QUIET,NORMAL or FULL",
         "value": "NORMAL"},
        
        "NNODES":
        {"comment": "Number of nodes (for clusters)",
         "value": 1},
        
        "NODE_INDEX":
        {"comment": "Node index (for clusters)",
         "value": 0},

        "NTHREADS":
        {"comment": "1 single thread",
         "value": 1},
        

        # -- Extra-keys (will not be saved in the main configuration file

        "CONFIG_FILE":
        {"comment": '[Extra key] name of the main configuration file',
         "value": "SWARP.conf"}
    }
        
    
    # -- Special config. keys that should not go into the config. file.

    _SW_config_special_keys = ["CONFIG_FILE"]



    def __init__(self):
        """
        SWARP class constructor.
        If a specific config_file is provided, it is used 
        """

        self.config = (
            dict([(k, copy.deepcopy(SWARP._SW_config[k]["value"]))\
                  for k in SWARP._SW_config.keys()]))

        # Extra config parameters that will be added/updated to the current values of the config file
        self.ext_config = {}       

        self.program = None
        self.version = None


    def setup(self, path=None):
        """
        Look for SWARP program ('swarp').
        If a full path is provided, only this path is checked.
        Raise a SWARPException if it failed.
        Return program and version if it succeed.
        """

        # -- Finding SWARP program and its version
        # first look for 'swarp'

        candidates = ['swarp']

        if (path):
            candidates = [path]
        
        selected=None
        for candidate in candidates:
            try:
                (_out_err, _in) = popen2.popen4(candidate)
                versionline = _out_err.read()
                if (versionline.find("SWarp") != -1):
                    selected=candidate
                    break
            except IOError:
                continue
                
        if not(selected):
            raise SWARPException, \
                  """
                  Cannot find SWARP program. Check your PATH,
                  or provide the SWARP program path in the constructor.
                  """

        _program = selected

        #print versionline
        _version_match = re.search("[Vv]ersion ([0-9\.])+", versionline)
        if not _version_match:
            raise SWARPException, \
                  "Cannot determine SWARP version."

        _version = _version_match.group()[8:]
        if not _version:
            raise SWARPException, \
                  "Cannot determine SWARP version."

        # print "Use " + self.program + " [" + self.version + "]"

        return _program, _version



    def update_config(self):
        """
        Update the configuration files according to the current
        in-memory SWARP configuration.
        """
        

        # -- Write main configuration file

        main_f = __builtin__.open(self.config['CONFIG_FILE'], 'w')

        for key in self.config.keys():
            if (key in SWARP._SW_config_special_keys):
                continue

            if (key == "CENTER"): # tuple instead of a single value
                value = " ".join(map(str, self.config[key]))
            else:
                value = str(self.config[key])
            
            print >>main_f, ("%-16s       %-16s # %s" %
                             (key, value, SWARP._SW_config[key]['comment']))

        main_f.close()


    def run(self, file_list, updateconfig=True, clean=False, path=None):
        """
        Run SWARP for a given list of files (fits files), and it can be one single file

        If updateconfig is True (default), the configuration
        files will be updated before running SWARP.

        If clean is True (default: False), configuration files 
        (if any) will be deleted after SWARP terminates.

        """

        if updateconfig:
            self.update_config()

        # Try to find SWARP program
        # This will raise an exception if it failed

        self.program, self.version = self.setup(path)
        
        # check how many files in the input
        if type(file_list) == types.ListType:
            my_files=""
            for file in file_list:
               my_files = my_files + " " + file
        else:
            my_files=file_list # a single file
                
        # Compound extra config command line args
        ext_args=""
        for key in self.ext_config.keys():
            ext_args=ext_args + " -" + key+ " " + str(self.ext_config[key])
               
        commandline = (
            self.program + " -c " + self.config['CONFIG_FILE'] + " " + ext_args + " " + my_files)
        
        #print commandline

        #rcode = os.system(commandline)
        rcode = misc.utils.runCmd(commandline)
        
        if (rcode==0):
            raise SWARPException, \
                  "SWARP command [%s] failed." % commandline
            
        if clean:
            self.clean()



    def clean(self, config=True):
        """
        Remove the generated SWARP files (if any).
        If config is True, remove generated configuration files.
        """

        try:
            if (config):
                os.unlink(self.config['CONFIG_FILE'])
                
        except OSError:
            pass



# ======================================================================
if __name__ == "__main__":


    swarp = SWARP()
    # Using a specific config file (updateconfig=False)
    swarp.config['CONFIG_FILE']="/disk-a/caha/panic/DEVELOP/PIPELINE/PANIC/trunk/config_files/swarp.conf"
    files = [line.replace( "\n", "") for line in fileinput.input(sys.argv[1])]
    swarp.run(files, updateconfig=False, clean=False)
    
    # Using and creating internal default config file (updateconfig=True)
    #swarp2 = SWARP()
    #cat_files = [line.replace( "\n", "") for line in fileinput.input(sys.argv[1])]
    #swarp2.run(cat_files, updateconfig=True, clean=False)
    
    sys.exit()
