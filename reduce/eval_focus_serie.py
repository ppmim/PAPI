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
import os.path
import dircache
import fileinput

import numpy as np
import pyfits
import matplotlib.pyplot as plt


import checkQuality

class TFOCUSNotFound(StandardError):
    pass

class FocusSerie(object):
    """
    Class used to estimate the best focus value of a focus exposures 
    
    Author: 
        JMIbannez, IAA-CSIC - 2011
        
    """
    
    def __init__(self, input_files, output, pix_size, sat_level, show=False, 
                    window='all', *a, **k):
        """
        Init method.

        Parameters
        ----------
        input_files: str
            File list of files to be processed.
        
        output: str
            Filename of the output pdf plot.
        
        pix_size: float
            Pixel size of input files.

        sat_level: int
            Saturation level of pixels (in ADUs).
        
        show: bool
            whether or not to show the pdf plot file generated.
        
        window: str
            Window/detector to be processed.

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
        self.show = show # whether or not to show the pdf plot file genetated
        self.window = window

    def eval_serie(self):
        """
        Performs the focus evaluation and write the results in a pdf plot and 
        the best focus in the ~/tmp/ql_focus.

        Returns
        -------
        Two options:

        (1)    best_focus, output_file: TFOCUS is present 
                best focus value (mm) and plot file with the fitted curve.

        (2)    best_fwhm, filename:    TFOCUS is not present or incorrect
                    best FWHM (min) found and the FITS file for this min FWHM
        """
        
        fwhm_values = []
        focus_values = []
        good_files = []

        # Find out the FWHM of each image
        for file in self.input_files:
            try:
                print "Evaluating %s\n"%file    
                cq = checkQuality.CheckQuality(file, 
                                               pixsize=self.pix_size, 
                                               sat_level=self.sat_level,
                                               window=self.window)
                try:
                    fwhm = cq.estimateFWHM()[0]
                except Exception,e:
                    sys.stderr.write("Error while computing FWHM for file %s"%file)
                    sys.stderr.write(str(e))
                    continue
                else:
                    fwhm_values.append(fwhm)

                # Try to read Telescope Focus (T-FOCUS)
                try:
                    focus = self.get_t_focus(file)
                except TFOCUSNotFound, e:
                    # Because we could be interested in knowing the best FWHM
                    # of the files even if the do not have the TFOCUS, we use
                    # a special value (-1) for this purpose. Obviously, no
                    # Poly fit will be done, but only show the FWHM values
                    # obtained for each file.
                    focus = -1
                focus_values.append(focus)
                good_files.append(file)
                print " >> FWHM =%f, T-FOCUS =%f <<\n"%(fwhm, focus)
            except Exception,e:
                sys.stderr.write("Some error while processing file %s\n"
                    " >>Error: %s\n"%(file,str(e)))
                raise Exception("Some error while processing file %s"%file)
                #log.debug("Some error happened") 
    
        # First, check if we have good values (!=-1) for T-FOCUS
        good_focus_values = [ v for v in focus_values if v!=-1 ]
        
        # Secondly, we remove duplicated values of T-FOCUS
        # (not sure if it should be done...)
        #seen = set()
        #seen_add = seen.add
        #good_focus_values = [ x for x in good_focus_values 
        #                            if not (x in seen or seen_add(x))]
        
        if len(good_focus_values)>1 and len(good_focus_values)==len(fwhm_values):
            # Fit the the values to a 2-degree polynomial
            sys.stdout.write("Focus values : %s \n"%str(good_focus_values))
            sys.stdout.write("FWHM values : %s \n"%str(fwhm_values))
            print "Lets do the fit...."
            z = np.polyfit(good_focus_values, fwhm_values, 2)
            print "Fit = %s  \n"%str(z)
            pol = np.poly1d(z)
            xp = np.linspace(np.min(good_focus_values), 
                    np.max(good_focus_values), 20000) # nro. puntos a interpolar
            best_focus = xp[pol(xp).argmin()]

            # Plotting
            plt.plot(good_focus_values, fwhm_values, '.', xp, pol(xp), '-')
            plt.title("Focus serie - Fit: %f X^2 + %f X + %f\n Best Focus=%f Detector=%s" 
                      %(pol[0],pol[1],pol[2],best_focus, self.window))
            plt.xlabel("T-FOCUS (mm)")
            plt.ylabel("FWHM (pixels)")
            plt.xlim(np.min(good_focus_values),np.max(good_focus_values))
            plt.ylim(pol(xp).min(), np.max(fwhm_values))
            # Annotate filenames on points
            for indx, i_file in enumerate(good_files):
                plt.annotate(os.path.basename(i_file), 
                    xy=(good_focus_values[indx], fwhm_values[indx]),  
                    xycoords='data',
                    xytext=(-50, 30), textcoords='offset points',
                    arrowprops=dict(arrowstyle="->")
                    )
            plt.savefig(self.output)
            if self.show:
                plt.show(block=True)

            sys.stdout.write("\nPlot generated: %s\n"%self.output)
            
            # Print out Values
            for idx, i_file in enumerate(good_files):
                sys.stdout.write("\nFile: %s   -->  FWHM: %s"%(i_file, fwhm_values[idx]))

            #
            # Write focus value into text file for OT
            #
            from os.path import expanduser
            home = expanduser("~")
            ql_focus_text_file = home + "/tmp/ql_focus"

            if not os.path.isdir(home + "/tmp"):
                msg = "tmp directory %s not found. Using %s directory"
                sys.stderr.write(msg % (home + "/tmp", home ))
                ql_focus_text_file = home + "/ql_focus"

            with open(ql_focus_text_file, "w") as text_file:
                # best_focus [mm] are converted to [microns] and written to
                # text file ready to be read and used by OT.
                text_file.write("%d"%int(round(best_focus*1000)))
            #
            # End-of-focus-file-writing
            #
            return best_focus, self.output

        else:
            if len(fwhm_values)==0:
                raise Exception("Empty list of FWHM values.")

            # Because cannot read TFOCUS values, only FWHM per file are shown. 
            sys.stdout.write("** Because cannot read properly the TFOCUS values,"
                             " only the FWHM per file is shown **")
            for idx, i_file in enumerate(good_files):
                sys.stdout.write("\nFile: %s   -->  FWHM: %s"%(i_file, fwhm_values[idx]))

            min_fwhm = np.min(fwhm_values)
            min_filename = good_files[np.argmin(fwhm_values)]
            return min_fwhm, min_filename
           
    def get_t_focus(self, file):
        """
        Look for the focus value into the FITS header keyword "T-FOCUS"
        
        Parameters
        ----------
        file: str
            Name of FITS file to look for 'T-FOCUS' keyword

        Returns
        ------- 
        The "T-FOCUS" keyword value read.
        
        """ 
                
        focus = -1           
        try:
            fits = pyfits.open(file)
            if "T-FOCUS" in fits[0].header:
                focus = fits[0].header["T-FOCUS"]
            elif "T_FOCUS" in fits[0].header:
                focus = fits[0].header["T_FOCUS"]
            else:
                sys.stderr.write("Cannot find the T-FOCUS value")
                raise TFOCUSNotFound("Cannot find the T-FOCUS value")
        finally:
            fits.close()
        
        return focus
        
        
def check_python_env():
    """    
    Check for Python 2.X with X >= 5; the 'optparse' module needs
    Python 2.5 for the used 'epilog' feature (see below).
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
    
    """)
    
    parser.add_option("-i", "--input", dest="input",
                      help="Name of input directory or list file "
                      "[default: %default].",
                      default="")
    parser.add_option("-o", "--output", dest="output",
                      help="Name of output [pdf] file with results " 
                      "[default: %default]",
                      default="output.pdf")
    
    parser.add_option("-p", "--pix_scale", dest="pix_scale",
                      help="Pixel scale [default: %default].",
                      default=0.45)
    
    parser.add_option("-s", "--satur_level", dest="satur_level",
                      help="Saturation level in ADUs. NCOADD is not taken "
                      "into account [default: %default].",
                      default=50000)

    parser.add_option('-W', '--window',
                      type='choice',
                      action='store',
                      dest='window',
                      choices=['Q1', 'Q2', 'Q3', 'Q4', 'all'],
                      default='all',
                      help="When input is a MEF, it means the "
                      "window/dectector/extension to process: "
                      "Q1, Q2, Q3, Q4, full [default: %default]")

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
                        options.pix_scale, options.satur_level,
                        False, options.window )

        best_focus = focus_serie.eval_serie()
        
        print "\nBEST_FOCUS =",best_focus
    except Exception,e:
        sys.stderr.write("Could not process focus serie. --> '%s'\n" %str(e))
        sys.exit(1)
    
        
    # and bye:
    sys.exit(0)

