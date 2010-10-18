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
# scamp module.
#
# Author: Jose Miguel Ibanez Mengual <jmiguel@iaa.es>
#
# ======================================================================

"""
A wrapper for SCAMP

A wrapper for SCAMP, the Source Extractor.
by Jose M. Ibanez Mengual
version: 0.1 - last modified: 2010-06-08

This wrapper allows you to configure SCAMP, run it and get
back its outputs without the need of editing SCAMP
configuration files. by default, configuration files are created
on-the-fly, and SCAMP is run silently via python.

Tested on SCAMP versions 1.4.6 and 1.7.0


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

# ======================================================================

__version__ = "0.1 (2010-06-08)"

# ======================================================================

class SCAMPException(Exception):
    pass

# ======================================================================

class SCAMP:
    """
    A wrapper class to transparently use SCAMP.

    """

    _SC_config = { 

#----------------------------- Field grouping ---------------------------------
    
        "FGROUP_RADIUS":
        {"comment": "Max dist (deg) between field groups",
         "value": 1.0},

#---------------------------- Reference catalogs ------------------------------
        
        "REF_SERVER":
        {"comment": "Internet addresses of catalog servers",
         "value": "cocat1.u-strasbg.fr"},
        
        "REF_PORT":
        {"comment": "Ports to connect to catalog servers",
         "value": 1660},
        
        "CDSCLIENT_EXEC":
        {"comment": 'CDSclient executable',
         "value": "aclient"},
        
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
         "value": "scamp.cat"},
        
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
         "value": [1,1]},
        
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

#---------------------------- Photometric solution ----------------------------
        "SOLVE_PHOTOM":
        {"comment": 'Compute photometric solution (Y/N) ?',
         "value": "N"},
         
        "MAGZERO_OUT":
        {"comment": 'Magnitude zero-point(s) in output',
         "value": 0.0},
 
        "MAGZERO_INTERR":
        {"comment": 'Internal mag.zero-point accuracy',
         "value": 0.01},
         
        "MAGZERO_REFERR":
        {"comment": 'Photom.field mag.zero-point accuracy',
         "value": 0.03},

        "PHOTINSTRU_KEY":
        {"comment": 'FITS keyword(s) defining the photom.',
         "value": "INSFLNAM"},

        "MAGZERO_KEY":
        {"comment": 'FITS keyword for the mag zero-point',
         "value": "PHOT_C"},

        "EXPOTIME_KEY":
        {"comment": 'FITS keyword for the exposure time (s)',
         "value": "EXPTIME"},

        "AIRMASS_KEY":
        {"comment": 'FITS keyword for the airmass',
         "value": "AIRMASS"},

        "EXTINCT_KEY":
        {"comment": 'FITS keyword for the extinction coeff',
         "value": "PHOT_K"},

        "PHOTOMFLAG_KEY":
        {"comment": 'FITS keyword for the photometry flag',
         "value": "PHOTFLAG"},
         
        "PHOTFLUX_KEY":
        {"comment": 'Catalog param. for the flux measurement',
         "value": "FLUX_AUTO"},         

        "PHOTFLUXERR_KEY":
        {"comment": 'Catalog parameter for the flux error',
         "value": "FLUXERR_AUTO"},
         
         "PHOTCLIP_NSIGMA":
        {"comment": 'Photom.clipping threshold in sigmas',
         "value": 3.0},
         
#------------------------------- Check-plots ----------------------------------
         
         "CHECKPLOT_CKEY":
        {"comment": 'FITS keyword for PLPLOT field colour',
         "value": "NULL"},
         
         "CHECKPLOT_DEV":
        {"comment": 'NULL, XWIN, TK, PS, PSC, XFIG, PNG, JPEG, AQT, PDF or SVG',
         "value": "PNG"},
         
         "CHECKPLOT_RES":
        {"comment": 'Check-plot resolution (0 = default)',
         "value": 0},
         
         "CHECKPLOT_ANTIALIAS":
        {"comment": 'Anti-aliasing using convert (Y/N) ?',
         "value": "Y"},
         
         "CHECKPLOT_TYPE":
        {"comment": 'Check-plots to perform',
         "value": "FGROUPS,DISTORTION,ASTR_INTERROR2D,ASTR_INTERROR1D,ASTR_REFERROR2D,ASTR_REFERROR1D,ASTR_CHI2,PHOT_ERROR"},
         
         "CHECKPLOT_NAME":
        {"comment": 'Check-plot filename(s)',
         "value": "fgroups,distort,astr_interror2d,astr_interror1d,astr_referror2d,astr_referror1d,astr_chi2,psphot_error"},
         
#------------------------------- Check-images ---------------------------------
         "CHECKIMAGE_TYPE":
        {"comment": 'NONE, AS_PAIR, AS_REFPAIR, or AS_XCORR',
         "value": "NONE"},
         
         "CHECKIMAGE_NAME":
        {"comment": 'Check-image filename(s)',
         "value": "check.fits"},
         

#------------------------------ Miscellaneous ---------------------------------

         "SN_THRESHOLDS":
        {"comment": 'S/N thresholds (in sigmas) for all and high-SN sample',
         "value": [10.0, 100.0] },
         
         "FWHM_THRESHOLDS":
        {"comment": 'FWHM thresholds (in pixels) for sources',
         "value": [0.0, 100.0]},
         
         "FLAGS_MASK":
        {"comment": 'Rejection mask on SEx FLAGS',
         "value": "0x00f0"},
         
         "WEIGHTFLAGS_MASK":
        {"comment": 'Rejection mask on SEx FLAGS_WEIGHT',
         "value": "0x00ff"},
         
         "IMAFLAGS_MASK":
        {"comment": 'Rejection mask on SEx IMAFLAGS_ISO',
         "value": "0x0"},
         
         "AHEADER_GLOBAL":
        {"comment": 'Filename of the global INPUT header',
         "value": "scamp.ahead"},
         
         "AHEADER_SUFFIX":
        {"comment": 'Filename extension for additional INPUT headers',
         "value": ".ahead"},
         
         "HEADER_SUFFIX":
        {"comment": 'Filename extension for OUTPUT headers',
         "value": ".head"},
         
         "HEADER_TYPE":
        {"comment": 'NORMAL or FOCAL_PLANE',
         "value": "NORMAL"},
         
         "VERBOSE_TYPE":
        {"comment": 'QUIET, NORMAL, LOG or FULL',
         "value": "NORMAL"},
         
         "WRITE_XML":
        {"comment": 'Write XML file (Y/N)?',
         "value": "N"},
         
         "XML_NAME":
        {"comment": 'Filename for XML output',
         "value": "scamp.xml"},
         
         "XSL_URL":
        {"comment": 'Filename for XSL style-sheet',
         "value": "file:///usr/share/scamp/scamp.xsl"},
         
         "NTHREADS":
        {"comment": 'umber of simultaneous threads for the SMP version of SCAMP (0=automatic)',
         "value": 0},
         


        # -- Extra-keys (will not be saved in the main configuration file

        "CONFIG_FILE":
        {"comment": '[Extra key] name of the main configuration file',
         "value": "scamp.conf"}
    }
        
    
    # -- Special config. keys that should not go into the config. file.

    _SC_config_special_keys = ["CONFIG_FILE"]



    def __init__(self):
        """
        SCAMP class constructor.
        If a specific config_file is provided, it is used 
        """

        self.config = (
            dict([(k, copy.deepcopy(SCAMP._SC_config[k]["value"]))\
                  for k in SCAMP._SC_config.keys()]))

        # Extra config parameters that will be added/updated to the current values of the config file
        self.ext_config = {}
              
        self.program = None
        self.version = None


    def setup(self, path=None):
        """
        Look for SCAMP program ('scamp').
        If a full path is provided, only this path is checked.
        Raise a SCAMPException if it failed.
        Return program and version if it succeed.
        """

        # -- Finding scamp program and its version
        # first look for 'scamp'

        candidates = ['scamp']

        if (path):
            candidates = [path]
        
        selected=None
        for candidate in candidates:
            try:
                (_out_err, _in) = popen2.popen4(candidate)
                versionline = _out_err.read()
                if (versionline.find("SCAMP") != -1):
                    selected=candidate
                    break
            except IOError:
                continue
                
        if not(selected):
            raise SCAMPException, \
                  """
                  Cannot find SCAMP program. Check your PATH,
                  or provide the SCAMP program path in the constructor.
                  """

        _program = selected

        #print versionline
        _version_match = re.search("[Vv]ersion ([0-9\.])+", versionline)
        if not _version_match:
            raise SCAMPException, \
                  "Cannot determine SCAMP version."

        _version = _version_match.group()[8:]
        if not _version:
            raise SCAMPException, \
                  "Cannot determine SCAMP version."

        # print "Use " + self.program + " [" + self.version + "]"

        return _program, _version



    def update_config(self):
        """
        Update the configuration files according to the current
        in-memory SCAMP configuration.
        """
        

        # -- Write main configuration file

        main_f = __builtin__.open(self.config['CONFIG_FILE'], 'w')

        for key in self.config.keys():
            if (key in SCAMP._SC_config_special_keys):
                continue

            if (key == "DISTORT_GROUPS" or key == "SN_THRESHOLDS" or key == "FWHM_THRESHOLDS"): # tuple instead of a single value
                value = " ".join(map(str, self.config[key]))
            else:
                value = str(self.config[key])
            
            print >>main_f, ("%-16s       %-16s # %s" %
                             (key, value, SCAMP._SC_config[key]['comment']))

        main_f.close()


    def run(self, catalog_list, updateconfig=True, clean=False, path=None):
        """
        Run SCAMP for a given list of catalog (.ldac files), and it can be one single catalog list

        If updateconfig is True (default), the configuration
        files will be updated before running SCAMP.

        If clean is True (default: False), configuration files 
        (if any) will be deleted after SCAMP terminates.

        """

        if updateconfig:
            self.update_config()

        # Try to find SCAMP program
        # This will raise an exception if it failed

        self.program, self.version = self.setup(path)
        
        # check how many files in the input
        my_catalogs=""
        if type(catalog_list) == types.ListType:
            for file in catalog_list:
               my_catalogs = my_catalogs + " " + file
        else:
            my_catalogs=catalog_list # a single file
        
            
        # Compound extra config command line args
        ext_args=""
        for key in self.ext_config.keys():
            ext_args=ext_args + " -" + key+ " " + str(self.ext_config[key])     
            
        commandline = (
            self.program + " -c " + self.config['CONFIG_FILE'] + " " + ext_args + " " + my_catalogs)
        print commandline

        rcode = os.system(commandline)
        
        if (rcode):
            raise SCAMPException, \
                  "SCAMP command [%s] failed." % commandline
            
        if clean:
            self.clean()



    def clean(self, config=True, catalog=False, check=False):
        """
        Remove the generated SCAMP files (if any).
        If config is True, remove generated configuration files.
        If catalog is True, remove the output catalog.
        If check is True, remove output check image.
        """

        try:
            if (config):
                os.unlink(self.config['CONFIG_FILE'])
            if (check):
                os.unlink(self.config['CHECKIMAGE_NAME'])
                
        except OSError:
            pass



# ======================================================================
if __name__ == "__main__":


    scamp = SCAMP()
    # Using a specific config file (updateconfig=False)
    scamp.config['CONFIG_FILE']="/disk-a/caha/panic/DEVELOP/PIPELINE/PANIC/trunk/config_files/scamp.conf"
    scamp.config['ASTREF_CATALOG']="2MASS"
    cat_files = [line.replace( "\n", "") for line in fileinput.input(sys.argv[1])]
    scamp.run(cat_files, updateconfig=False, clean=False)
    
    # Using and creating internal default config file (updateconfig=True)
    #scamp2 = SCAMP()
    #scamp2.config['ASTREF_CATALOG']="2MASS"
    #cat_files = [line.replace( "\n", "") for line in fileinput.input(sys.argv[1])]
    #scamp2.run(cat_files, updateconfig=True, clean=False)
    
    sys.exit()
