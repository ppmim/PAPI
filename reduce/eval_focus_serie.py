#!/usr/bin/env python

# script to calculate 'best' focus value taking
# into account the FWHM of a focus exposures.

# module 're' is for regular expressions:
import re
import string
import math
import optparse 
import sys
import os

import checkQuality

class EvalFocusSerie:
    """
    @summary:  
        Class used to estimate the best focus value of a focus esposures 
    
    @author: 
        JMIbannez, IAA-CSIC
        
    """
    
    def __init__(self, input_files):
        self.input_files = input_files

    def eval_serie(self):
        """
        Do the focus evaluation
        """
        
        fwhm_values = []
        focus_values = []
        # Find out the FWHM of each image
        for file in self.input_files:
            try:
                cq = CheckQuality(file)
                fwhm_values.append(fwhm)
                fwhm = cq.extimateFWHM()
                focus = ClFits(file).getFocus()
            except Exception,e:
                sys.stderr.writhe("Some error in CheckQuality")
                #log.debug("Some error happened") 
                
                

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

 eval_focus_serie():


################################################################################
# main
if __name__ == "__main__":

    check_python_env()

    # The following class allows us to define a usage message epilogue which does
    # not skip newline characters:
    class MyParser(optparse.OptionParser):
        def format_epilog(self, formatter):
            return self.epilog
    
    usage = "%prog [Options] - compute the best focus value of a focus serie"
    parser = MyParser(usage=usage, epilog=
    """
    Description:
    Given a focus exposures serie as a set of FITS files images, the FWHM is 
    computed for each of the images, and the best focus is determined by 
    interpolation.
    The FWHM is computation is based on the SExtractor catalog generated.
    
    Example:
    - eval_focus_serie.py -s /data-dir/focus -o fwhm_values.txt
    
    - eval_focus_serie.py -s /data-dir/focus.list -o fwhm_values.txt
    
    Known Bugs/Shortcomings:
    
    Author:
      J.M. Ibanez       (jmiguel@iaa.es)
    
    """)
    
    parser.add_option("-i", "--input", dest="input",
                      help="name of input directory or list file (default: %default)",
                      default=os.curdir)
    parser.add_option("-o", "--output", dest="output",
                      help="name of output file with results (default: %default)",
                      default="output.txt")
    
    (options, args) = parser.parse_args()
    
    try:
        hdu = pyfits.open(options.input)
    except:
        sys.stderr.write("Could not read input catalogue '%s'\n" % (options.input))
        parser.print_help()
        sys.exit(1)
    
        
    # and bye:
    sys.exit(0)

