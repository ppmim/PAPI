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
        
        "ASTREF_CATALOG":
        {"comment": "NONE, FILE, USNO-A1, USNO-A2, USNO-B1, GSC-1.3, GSC-2.2, UCAC-1, UCAC-2,NOMAD-1, 2MASS, DENIS-3, SDSS-R3, SDSS-R5 or SDSS-R6 ",
         "value": "2MASS"},
        
        "ASTREF_BAND":
        {"comment": "Photom. band for astr.ref.magnitudes or DEFAULT, BLUEST, or REDDEST",
         "value": "DEFAULT"},
        
        "ASTREFCAT_NAME":
        {"comment": "Local astrometric reference catalogs",
         "value": "astrefcat.cat"},
        
        "ASTREFCENT_KEYS":
        {"comment": "Local ref.cat.centroid parameters",
         "value": "X_WORLD,Y_WORLD"},
        
        "ASTREFERR_KEYS":
        {"comment": 'Local ref.cat.error ellipse parameters',
         "value": 'ERRA_WORLD, ERRB_WORLD, ERRTHETA_WORLD'},
        
        "ASTREFMAG_KEY":
        {"comment": "Local ref.cat.magnitude parameter",
         "value": "MAG"},
        
        "SAVE_REFCATALOG":
        {"comment": "Save ref catalogs in FITS-LDAC format?",
         "value": "N"},

#--------------------------- Merged output catalogs ---------------------------
       
        "MERGEDOUTCAT_NAME":
        {"comment": "Merged output catalog filename",
         "value": "SWARP.cat"},
        
        "MERGEDOUTCAT_TYPE":
        {"comment": "NONE, ASCII_HEAD, ASCII, FITS_LDAC",
         "value": 'NONE'},
        
#----------------------------- Pattern matching -------------------------------
        
        "MATCH":
        {"comment": "Do pattern-matching (Y/N) ?",
         "value": "Y"},
        
        "MATCH_NMAX":
        {"comment": 'Max.number of detections for MATCHing(0=auto)',
         "value": 0},
        
        "PIXSCALE_MAXERR":
        {"comment": "Max scale-factor uncertainty",
         "value": 1.2},
        
        "POSANGLE_MAXERR":
        {"comment": 'Max position-angle uncertainty (deg)',
         "value": 5.0},
        
        "POSITION_MAXERR":
        {"comment": "Max positional uncertainty (arcmin)",
         "value": 1.0},
        
        "MATCH_RESOL":
        {"comment": "Matching resolution (arcsec); 0=auto",
         "value": 0},
        
        "MATCH_FLIPPED":
        {"comment": "Allow matching with flipped axes?",
         "value": "N"},
        
        "FIXFOCALPLANE_NMIN":
        {"comment": "Min number of dets for FIX_FOCALPLANE",
         "value": 1},

#---------------------------- Cross-identification ----------------------------
        
        "CROSSID_RADIUS":
        {"comment": "Cross-id initial radius (arcsec)",
         "value": 2.0},

#---------------------------- Astrometric solution ----------------------------
        
        "SOLVE_ASTROM":
        {"comment": "Compute astrometric solution (Y/N) ?",
         "value": "Y"},
        
        "ASTRINSTRU_KEY":
        {"comment": "FITS keyword(s) defining the astrom",
         "value": "INSTRID, INSFLNAM"},
        
        "STABILITY_TYPE":
        {"comment": "EXPOSURE, GROUP, INSTRUMENT or FILE",
         "value": "INSTRUMENT"},
        
        "CENTROID_KEYS":
        {"comment": "Cat. parameters for centroiding",
         "value": "XWIN_IMAGE,YWIN_IMAGE"},
        
        "CENTROIDERR_KEYS":
        {"comment": 'Cat. params for centroid err ellipse',
         "value": "ERRAWIN_IMAGE,ERRBWIN_IMAGE,ERRTHETAWIN_IMAGE"},
        
        "DISTORT_KEYS":
        {"comment": 'Cat. parameters or FITS keywords',
         "value": "XWIN_IMAGE,YWIN_IMAGE"},
        
        "DISTORT_GROUPS":
        {"comment": "Filename for the check-image",
         "value": "1,1"},
        
        "DISTORT_DEGREES":
        {"comment": "Polynom degree for each group",
         "value": 3},
        
        "ASTREF_WEIGHT":
        {"comment": "Relative weight of ref.astrom.cat.",
         "value": 1.0},
        
        "ASTRCLIP_NSIGMA":
        {"comment": "Astrom. clipping threshold in sigmas",
         "value": 3.0},
        
        "CORRECT_COLOURSHIFTS":
        {"comment": 'Correct for colour shifts (Y/N)?',
         "value": "N"},



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

        # print self.config

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
        
        my_files=""
        for file in file_list:
           my_files = my_files + " " + file
            
        commandline = (
            self.program + " -c " + self.config['CONFIG_FILE'] + " " + my_files)
        print commandline

        rcode = os.system(commandline)
        
        if (rcode):
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
