#!/usr/bin/env python

# Copyright (c) 2009-2012 IAA-CSIC  - All rights reserved. 
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


################################################################################
#
# PAPI (PANIC PIpeline)
#
# spatial_noise.py
#
# Created    : 08/11/2012    jmiguel@iaa.es -
# Last update: 
# TODO
#       
################################################################################

# System modules
from optparse import OptionParser
import sys
import itertools
import tempfile
import os

import numpy
import math
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib.cm as cm
import pylab
            

# PAPI modules
from misc.paLog import log
from misc.utils import *
from misc.print_table import print_table
    
def run_spatial_noise ( input_catalog, area, window, gain, 
                        out_filename="/tmp/hc_out.pdf"):
    """ 
    Takes a input catalog (ascii file) listing all the files to be used in the
    spatial noise analysis and performs the computation required for it. 
    
    Parameters
    ----------
    input_catalot: str
        Text file listing the files to use for noise computation 
    area: str
        Area in frame (Q1, Q2, Q3, Q4, full-detector, central) to use
    window: list of int coordinates
        Window coordinates(x1,y1,x2,y2) to use for computation
    gain: float
        Detector gain
    out_filename: str
        filename where results will be saved

    Returns
    -------
    out_filename: str
        Filename where results where saved
    
    TODO
    ----
      - print pretty stats in output, not only plots
      - Do computations per channel and/or detector
      - Allow a custom window size (coordinates)
      
    """
    
    
    # Read the file list from input_file
    if os.path.exists(input_catalog):
        filelist = [line.replace( "\n", "") for line in fileinput.input(input_catalog)]
        n_files = len(filelist)
    else:
        raise Exception("Input file %s does not exist"%input_catalog)
    
    #Read image shape
    with pyfits.open(filelist[0]) as pf:
        shape = pf[0].data.shape
    
    #Define window
    if window:
        x1 = window[0]
        y1 = window[1]
        x2 = window[2]
        y2 = window[3]
        s_area = [x1, y1, x2, y2] 
    elif area == 'Q1':
        x1 = 10
        y1 = 10
        x2 = 2030
        y2 = 2030
        s_area = [10, 10, 2030, 2030]
    elif area == 'Q2':
        x1 = 2100
        y1 = 10
        x2 = 4080
        y2 = 2100
        s_area = [2100, 10, 4080, 2100]    
    elif area == 'Q3':
        x1 = 2100
        y1 = 2100
        x2 = 4080
        y2 = 4080
        s_area = [2100, 2100, 4080, 4080]    
    elif area == 'Q4':
        x1 = 10
        y1 = 2100
        y2 = 2100
        x2 = 4080
        s_area = [10, 2100, 2100, 4080]
    elif area == 'central':
        x1 = 512
        y1 = 512
        x2 = 3584
        y2 = 3584
        s_area = [512, 512, 3584, 3584]
    else: # or 'full-frame'
        x1 = 10
        x2 = shape[0]-10
        y1 = 10
        y2 = shape[1]-10
        s_area = [x1, y1, x2, y2] 
    
    print "Selected window = ",s_area
    print "Gain = ", gain
    print "Files (%d) = %s" %(n_files, filelist)
    
    # check window-shape
    pf = pyfits.open(filelist[0])
    if not (x1 < shape[0] and x2 < shape[0] and 
        y1 < shape[1] and y2 < shape[1]):
        raise Exception("Wrong window definition; check image and window size")
    pf.close() 

    std_sub = numpy.zeros([n_files/2], dtype=numpy.float32)
    itime = numpy.zeros([n_files/2], dtype=numpy.float32)
    n = 0
    key_time = 'EXPTIME' # 'ITIME'
    for packet in grouper(2, filelist):
        if len(packet)==2 and not (None in packet):
            try:
                pf1 = pyfits.open(packet[0])
                pf2 = pyfits.open(packet[1])
                if (key_time in pf1[0].header and 
                    key_time in pf2[0].header and
                    pf1[0].header[key_time] == pf2[0].header[key_time]): 
                    itime[n] = pf1[0].header[key_time]
                else:
                    log.error("ITIME mismatch, frames skipped %s"%str(packet))
                    continue
                if pf1[0].data.shape == pf2[0].data.shape:
                    # Do subtraction
                    std_sub[n] = numpy.std(pf1[0].data[x1:x2, y1:y2] - 
                                           pf2[0].data[x1:x2, y1:y2])
                else:
                    log.error("Shape mismatch; frames %s skipped"%(str(packet)))
                    continue
            except Exception, e:
                raise e
            pf1.close()
            pf2.close()
            n = n + 1


    noise = gain * std_sub
    mean_noise = numpy.mean(noise)
    
    t = [
        ["Packet", "ITime  ", "Std     ", "Noise"]
        ]
    i = 0        
    for row in zip(itime, std_sub, noise):
        t.append([str(i), str(row[0]), str(row[1]), str(row[2]) ])
        i += 1
    print_table(t)
    print "------------------------"
    #print "Std", std_sub
    #print "Itimes", itime
    print "Gain = ", gain
    print "Mean noise", mean_noise
    print "------------------------"
    
    # Compute the linear fit Itime VS stddev ( std = a*itime+b )
    # ==============================================================
    res = numpy.polyfit(itime, std_sub, 1, None, True)
    a = res[0][1] # intercept
    b = res[0][0] # slope
    r = res[3][1] # regression coeff ??? not exactly
    
    #print "Coeffs =", res
    
    # Plot the Signal(x) VS Variance(y)
    pol = numpy.poly1d(res[0])
    plt.plot(itime, std_sub, '.', itime, pol(itime), '-')
    plt.title("Poly fit: %f X + %f  r=%f Mean_Noise=%f" %(b, a, r, mean_noise))
    plt.xlabel("Itime / s")
    plt.ylabel("Noise  / e-")
    plt.savefig(out_filename)
    plt.show()
    

    return out_filename

def grouper(group_size, iterable, fillvalue=None):
    "grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * group_size
    return itertools.izip_longest(fillvalue=fillvalue, *args)

################################################################################
# main
################################################################################
if __name__ == "__main__":

    # Get and check command-line options
        
    usage = "usage: %prog [options] arg1 arg2 ..."
    desc = "Compute the Spatial Noise from a set of dark images grouped in \
    pairs with the same Integration Time"
    
    parser = OptionParser(usage, description=desc)
    
    parser.add_option("-i", "--input_images",
                  action="store", dest="input_images", 
                  help="input image list (darks and flats")
    
    parser.add_option("-a", "--area",
                  action="store", dest="area", type=str,
                  help="Area to use for computations")
    
    parser.add_option("-w", "--window",
                  action="store", dest="window", type=int, 
                  metavar='x1 y1 x2 y2', nargs=4,
                  help="Window to use for computations (exclusive with area)")
    
    parser.add_option("-g", "--gain",
                  action="store", dest="gain", type=float, default=1,
                  help="Detector gain (default=%default)")
                  
    parser.add_option("-o", "--output",
                  action="store", dest="output_file", 
                  help="output plot filename (default = %default)",
                  default="noise.pdf")
    
    parser.add_option("-v", "--verbose",
                  action="store_true", dest="verbose", default=True,
                  help="verbose mode [default=%default]")
    
                                
    (options, args) = parser.parse_args()
    
    if options.area and options.window:
        parser.print_help()
        parser.error("<Area> and <Window> arguments are exclusive")
        
    if not options.input_images or len(args)!=0: 
    # args is the leftover positional arguments after all options have been processed
        parser.print_help()
        parser.error("Wrong number of arguments")

    if not os.path.exists(options.input_images):
        log.error ("Input image %s does not exist", options.input_images)
        sys.exit(0)
    
    log.debug( 'Spatial noise computation starts')
    try:
        run_spatial_noise(options.input_images, options.area, options.window,
                         options.gain, options.output_file)
    except Exception, e:
        log.error("Some error while running routine: %s"%str(e))
        sys.exit(0)
        
