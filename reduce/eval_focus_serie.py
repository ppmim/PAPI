#!/usr/bin/env python

# Copyright (c) 2011 IAA-CSIC  - All rights reserved. 
# Author: Jose M. Ibanez. 
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


# ==============================================================================
# script to calculate 'best' focus value taking
# into account the FWHM of a focus exposures.
# ==============================================================================

import string
import math
import optparse 
import sys
import os
import dircache
import fileinput

import numpy as np
import pyfits
import matplotlib.pyplot as plt


import checkQuality

class FocusSerie(object):
    """
    @summary:  
        Class used to estimate the best focus value of a focus exposures 
    
    @author: 
        JMIbannez, IAA-CSIC - 2011
        
    """
    
    def __init__(self, input_files, output, pix_size, sat_level, show=False, *a, **k):
        """
        Init method
        """
        
        super (FocusSerie, self).__init__ (*a,**k)
        #import pdb
        #pdb.set_trace()

        if type(input_files)!=type(list()) and os.path.isfile(input_files):
            files = [line.replace("\n", "").replace('//','/')
                     for line in fileinput.input(input_files)]

            self.input_files = files
        else:
            self.input_files = input_files

        self.output = output
        self.pix_size = pix_size
        self.sat_level = sat_level
        self.show = show # whether to show or not the pdf plot file

    def eval_serie(self):
        """
        Do the focus evaluation
        """
        
        fwhm_values = []
        focus_values = []
        
        # Find out the FWHM of each image
        for file in self.input_files:
            try:
                print "Evaluating %s\n"%file    
                cq = checkQuality.CheckQuality(file, 
                                               pixsize=self.pix_size, 
                                               sat_level=self.sat_level)
                fwhm = cq.estimateFWHM()[0]
                fwhm_values.append(fwhm)
                focus = self.get_t_focus(file)
                focus_values.append(focus)
                print " >> FWHM =%f, T-FOCUS =%f <<\n"%(fwhm,focus)
            except Exception,e:
                sys.stderr.write("Some error while processing file %s\n >>Error: %s\n"%(file,str(e)))
                #log.debug("Some error happened") 
    
        # Fit the the values to a 2-degree polynomial
        sys.stdout.write("Focus values : %s \n"%str(focus_values))
        sys.stdout.write("FWHM values : %s \n"%str(fwhm_values))
        if len(focus_values)>0 and len(focus_values)==len(fwhm_values):
            print "Lets do the fit...."
            z = np.polyfit(focus_values, fwhm_values, 2)
            print "Fit = %s  \n"%str(z)
            pol = np.poly1d(z)
            xp = np.linspace(np.min(focus_values), np.max(focus_values), 2000)
            best_focus = xp[pol(xp).argmin()]

            # Plotting
            plt.plot(focus_values, fwhm_values, '.', xp, pol(xp), '-')
            plt.title("Focus serie - Fit: %f X^2 + %f X + %f\n Best Focus=%f" 
                      %(pol[0],pol[1],pol[2],best_focus))
            plt.xlabel("T-FOCUS (mm)")
            plt.ylabel("FWHM (pixels)")
            plt.xlim(np.min(focus_values),np.max(focus_values))
            plt.ylim(pol(xp).min(), np.max(fwhm_values))
            plt.savefig(self.output)
            if self.show:
                plt.show(block=True)

        else:
            print "Not enough data for fitting"
            best_focus = np.NaN

        log.info("Plot generated: %s"%self.output)
        
        return best_focus, self.output
           
    def get_t_focus(self, file):
        """
        @summary: Look for the focus value into the FITS header keyword "T-FOCUS"
        
        @return: the "T-FOCUS" keyword value
        """ 
                
        focus = -1           
        try:
            fits = pyfits.open(file)
            if "T-FOCUS" in fits[0].header:
                focus = fits[0].header["T-FOCUS"]
            elif "T_FOCUS" in fits[0].header:
                focus = fits[0].header["T_FOCUS"]
            else:
                raise Exception("Canno find the FOCUS value")
        except Exception,e:
            sys.stderr.write("Cannot find the T-FOCUS value in file %s\n"%file)
            raise e
        
        return focus
        
        
def check_python_env():
    """    
    check for Python 2.X with X >= 5; the 'optparse' module needs
    Python 2.5 for the used 'epilog' feature (see below):
    """
    version = string.split(string.split(sys.version) [0], ".")
    # well, Python version 3 just gives us a syntax error at the
    # first print statement :-)
    if map(int, version) >= [3, 0, 0] or map(int, version) < [2, 5, 0]:
        sys.stderr.write("This script needs Python 2.Y.X (with Y >= 5)\n\n")
        sys.stderr.write("You have Python V%s.%s.%s\n" \
                          % (version[0], version[1], version[2]))
        sys.exit(1)
        
    # import non-standard modules:
    try:
        import pyfits, numpy
    except ImportError:
        sys.stderr.write("This script needs the modules 'pyfits' and 'numpy'!\n")
        sys.stderr.write("see http://www.stsci.edu/resources/software_hardware/pyfits\n")
        sys.stderr.write("and/or http://numpy.scipy.org/\n")
        sys.exit(1)



################################################################################
# main
if __name__ == "__main__":

    check_python_env()

    # The following class allows us to define a usage message epilogue which does
    # not skip newline characters:
    class MyParser(optparse.OptionParser):
        def format_epilog(self, formatter):
            return self.epilog
    
    usage = "%prog [Options] "
    parser = MyParser(usage=usage, epilog=
    """
    Description:
       Given a focus exposures serie as a set of FITS files images, the FWHM is 
       computed for each of the images, and the best focus is determined by 
       interpolation.
       The FWHM is computation is based on the SExtractor catalog generated.
    
    Example:
       - eval_focus_serie.py -s /data-dir/focus -o fwhm_values.pdf
    
       - eval_focus_serie.py -s /data-dir/focus.list -o fwhm_values.pdf
    
    Known Bugs/Shortcomings:
    
    Author:
       J.M. Ibanez       (jmiguel@iaa.es)
    
    """)
    
    parser.add_option("-i", "--input", dest="input",
                      help="name of input directory or list file (default: %default)",
                      default="")
    parser.add_option("-o", "--output", dest="output",
                      help="name of output [pdf] file with results (default: %default)",
                      default="output.pdf")
    
    parser.add_option("-p", "--pix_scale", dest="pix_scale",
                      help="Pixel scale (default: %default)",
                      default=0.45)
    
    parser.add_option("-s", "--satur_level", dest="satur_level",
                      help="Saturation level in ADUs. NCOADD is not taken into account.(default: %default)",
                      default=50000)
    
    (options, args) = parser.parse_args()
    

    if len(sys.argv[1:])<1 or len(args)!=0:
       parser.print_help()
       sys.exit(0)

    # Read input files
    files = []    
    if os.path.isfile(options.input):
        files = [line.replace("\n", "").replace('//','/')
                     for line in fileinput.input(options.input)]
    elif os.path.isdir(options.input):
        for file in dircache.listdir(options.input):
            files.append(options.input+"/"+file)
    else:
        parser.print_help()
        sys.exit(1)
    
    # Eval the focus exposures
    try:
        focus_serie = FocusSerie( files , options.output, 
            options.pix_scale, options.satur_level )
        best_focus = focus_serie.eval_serie()
        print "BEST_FOCUS =",best_focus
    except Exception,e:
        sys.stderr.write("Could not process focus serie. --> '%s'\n" %str(e))
        sys.exit(1)
    
        
    # and bye:
    sys.exit(0)

