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
# health.py
#
# Tools used:
#
#    IRAF
#
# Created    : 29/10/2012    jmiguel@iaa.es -
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
import datahandler

def stack_frames_2 (frames, f_from, f_to, outframe, sigmaframe, 
                  data_range_low, data_range_high, kappa):
    """
    Computes average, median or mode for stack of frames
    
    Parameters
    ----------
    frames: list
        list of frames to be stacked
    f_from: int
        frame from start the stacking
    f_to: int 
        frame to end the stacking
    outframe: str
         filename where stack frame is saved to
    sigmaframe: str 
        filename where sigma frame is saved to
    data_range_low: long 
        min. value in range of accepted values in frame
    data_range_high: long 
        max. value in range of accepted values in frame
    kappa: float 
        number of sigmas to use in clipping during the stack
    
    Returns
    -------
    outframe, sigmaframe 
    """
    
    # First, some checks
    if len(frames)<2:
        raise Exception("Not enough number of frames")
    
    #cube = numpy.zeros([n_stripes, height_st, width_st], dtype=numpy.float)
    #data_out = numpy.zeros([n_stripes*height_st*2, width_st*2], dtype=numpy.float32)
    #median = numpy.zeros([4], dtype=numpy.float32)
    
    for i_frame in frames:
        try:
            pf = pyfits.open(i_frame)
        except Exception, e:
            raise e
        
        if len(pf)>1:
            raise Exception("MEF files not yet supported.")
        
        #cube[i] = f[0].data
        
            
    return outframe, sigmaframe

def stack_frames(file_of_frames, type_comb, f_from, f_to, outframe, sigmaframe, 
                  data_range_low, data_range_high, kappa):
    """
    Computes average, median or mode for stack of frames
    
    Parameters
    ----------
    file_of_frames: str
        text file listing the frames to stack
    type_comb: str
        type of combination to make (mean, median, sigma)
    f_from: int
        frame from start the stacking [0,...,N-1], where N = number of files
    f_to: int [1,...,N], where N = number of files
        frame to end the stacking
    (Both indexes follow python 0-index)
    outframe: str
         filename where stack frame is saved to
    sigmaframe: str 
        filename where sigma frame is saved to
    data_range_low: long 
        min. value in range of accepted values in frame
    data_range_high: long 
        max. value in range of accepted values in frame
    kappa: float 
        number of sigmas to use in clipping during the stack
    
    Returns
    -------
    outframe, sigmaframe 
    
    TODO
    ----
     - Implement the use of data_range and kappa for clipping -->cubemean
    """
    
    # First, some checks
    if not os.path.exists(file_of_frames):
        raise Exception("File does not exists")

    irdr_path = "/home/panic/DEVELOP/PIPELINE/PANIC/trunk/irdr/bin"
    outweight = "/tmp/weight.fits"
    
    # Compute the outframe
    mode = 'mean'
    if type_comb == 'mean': 
        mode = 'mode'
    elif type_comb == 'median': 
        mode = 'median'
    else:
        log.error("Type of stacking not supported")
        raise Exception("Type of stacking not supported")

    # Limit the number of frames
    tmp_list = fileToList(file_of_frames)
    if f_from != 0 or f_to != len(tmp_list):
        new_file_of_frames = file_of_frames + ".ranged"
        if f_from >= 0 and f_to <=len(tmp_list):
            listToFile(tmp_list[f_from : f_to], new_file_of_frames)
        else:
            print "Len=", len(tmp_list)
            print "f_to", f_to
            raise Exception("list index out of range")
    else:
        new_file_of_frames = file_of_frames
    
    # (use IRDR::cubemean)
    prog = irdr_path+"/cubemean "
    cmd  = prog + " " + new_file_of_frames + " " + outframe + " " + outweight + \
            " " + "offset " + mode + " noweight float"
    e = runCmd( cmd )
    if e==0:
        log.debug("Some error while running command %s", cmd)
        raise Exception("Some error while running command %s"%cmd)
        
    # Now, compute the sigmaframe
    cmd  = prog + " " + new_file_of_frames + " " + sigmaframe + " " + \
            outweight + " " + "offset sigma noweight float"
            
    e = runCmd( cmd )
    if e==0:
        log.debug("Some error while running command %s", cmd)
        raise Exception("Some error while running command %s"%cmd)
    
    log.debug("Successful ending of stack_frames")
    
    return (outframe, sigmaframe)
            
    
def run_health_check ( input_file, f_from, t_to, packet_size, window='full-frame',
                        out_filename="/tmp/hc_out.pdf"):
    """ 
    Takes a input catalog (ascii file) listing all the files to be used in the
    check-health analysis and performs the computation required for it. 
    
    Parameters
    ----------
    input_file: str
        Text file listing the files to use for health computation 
    f_from: int [0,...,N-1], where N = number of files
        File number where start the packet 
    t_to: int [1,...,N], where N = number of files
        File number where end the packet 
    packet_size: int
        Size of the packet
    window: str
        Window in frame (Q1, Q2, Q3, Q4)
    out_filename: str
        filename where results will be saved

    Returns
    -------
    out_filename: str
        Filename where results where saved
    
    Notes
    -----
    It is based on JWF MIDAS routine.
    
    TODO
    ----
      - Implement usage of f_from, f_to 

    """
    
    
    # Read the file list from input_file
    if os.path.exists(input_file):
        filelist = [line.replace( "\n", "") for line in fileinput.input(input_file)]
    else:
        raise Exception("Input file %s does not exist"%input_file)
    
    if window == 'Q1':
        x1 = 10
        y1 = 10
        x2 = 2030
        y2 = 2030
        area = [10, 10, 2030, 2030]
    elif window == 'Q2':
        x1 = 2100
        y1 = 10
        x2 = 4080
        y2 = 2100
        area = [2100, 10, 4080, 2100]    
    elif window == 'Q3':
        x1 = 2100
        y1 = 2100
        x2 = 4080
        y2 = 4080
        area = [2100, 2100, 4080, 4080]    
    elif window == 'Q4':
        x1 = 10
        y1 = 2100
        y2 = 2100
        x2 = 4080
        area = [10, 2100, 2100, 4080]
    elif window == 'central':
        x1 = 512
        y1 = 512
        x2 = 3584
        y2 = 3584
        area = [512, 512, 3584, 3584]
    else: # or 'full-frame'
        x1 = 10
        x2 = 10
        y1 = 4080
        y2 = 4080
        area = [10, 10, 4080, 4080] 
    
    print "Selected area = ",area
    print "Packet size", packet_size
    print "Files", filelist    
    #fileToList()
    #tmp_file, _ = tempfile.mkstemp()
    tmp_file = "/tmp/packet.txt"
    n = 0
    for packet in grouper(packet_size, filelist):
        if len(packet)==packet_size and not (None in packet):
            listToFile(packet, tmp_file)
            try:
                stack_frames(tmp_file, 'median', 1, 3, "/tmp/stack%02d.fits"%n,
                             "/tmp/stack_sigma%02d.fits"%n, 10, 100000, 3)    
            except Exception, e:
                raise e
            n = n + 1

            
    #Get stats from images
    stat_values = {}
    itime = numpy.zeros([n], dtype=numpy.float32)
    signal = numpy.zeros([n], dtype=numpy.float32)
    std = numpy.zeros([n], dtype=numpy.float32)
    
    for i in range(n):
        try:
            print "Reading file %s"%("/tmp/stack%02d.fits"%i) 
            pf = pyfits.open("/tmp/stack%02d.fits"%i)
            pf2 = pyfits.open("/tmp/stack_sigma%02d.fits"%i)
        except Exception, e:
            log.error("Cannot open file %s"%("/tmp/stack%02d.fits"%i))
            continue
        signal[i] = numpy.mean(pf[0].data[x1:x2, y1:y2])
        std[i] = numpy.std(pf2[0].data[x1:x2, y1:y2])
        if 'ITIME' in pf[0].header: 
            itime[i] = pf[0].header['ITIME']
            
        else:
            itime[i] = NaN
        stat_values[i] = [itime[i], signal[i], std[i]]
        pf.close()
        pf2.close()
        
    
    print "Signal",signal
    print "Std",std
    print "Itimes",itime
    
    # Compute the linear fit signal VS variance ( var = a*signal+b )
    res = numpy.polyfit(signal, std**2, 1, None, True)
    a = res[0][1] # intercept
    b = res[0][0] # slope
    r = res[3][1] # regression coeff ??? not exactly
    
    print "Coeffs =", res
    
    
    # Plot the Signal(x) VS Variance(y)
    pol = numpy.poly1d(res[0])
    plt.plot(signal, std**2, '.', signal, pol(signal), '-')
    #Note: gain = 1/slope
    #Note: intercept = (RON/gain)**2 => RON = gain*sqrt(intercept)
    gain = 1.0 / b
    ron = gain/math.sqrt(math.fabs(a))
    print "Gain = %s [e-/ADU]"%gain
    print "RON = %s e-"%ron
    plt.title("Poly fit: %f X + %f  r=%f Gain=%s RON=%s" %(b, a, r, 
                                                           gain, ron))
    plt.xlabel("Signal / ADU")
    plt.ylabel("Variance  / ADU")
    plt.savefig(out_filename)
    plt.show()
    
    ##
    # Fit and Plot the ITime(x) VS Signal(y)
    ##
    res = numpy.polyfit(itime, signal, 1, None, True)
    a = res[0][1] # intercept
    b = res[0][0] # slope
    r = res[3][1] # regression coeff ??? not exactly
    
    pol = numpy.poly1d(res[0])
    plt.clf()
    plt.plot(itime, signal, '.', itime, pol(itime), '-')
    plt.title("Poly fit: %f X + %f  r=%f " %(b, a, r))
    plt.xlabel("ITime / s")
    plt.ylabel("Signal  / ADU")
    plt.savefig(out_filename+"_2.pdf")
    plt.show()
    
    
    # remove tmp files
    os.unlink(tmp_file)

    return 0

def grouper(group_size, iterable, fillvalue=None):
    "grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * group_size
    return itertools.izip_longest(fillvalue=fillvalue, *args)

################################################################################
# main
################################################################################
if __name__ == "__main__":

    log.debug( 'Health-Check routines for PANIC')
    
    # Get and check command-line options
        
    usage = "usage: %prog [options] arg1 arg2 ..."
    parser = OptionParser(usage)
    
    parser.add_option("-i", "--input_images",
                  action="store", dest="input_images", 
                  help="input image list (darks and flats")
    
    parser.add_option("-w", "--window",
                  action="store", dest="window", type=str,
                  help="Window to use for computations (default = %default)",
                  default='full-detector')
    
    parser.add_option("-s", "--start_packet",
                  action="store", dest="start_packet", type=int, default=1,
                  help="File where start the packet (default=%default)")
    
    parser.add_option("-e", "--end_packet",
                  action="store", dest="end_packet", type=int, default=10,
                  help="File where end the packet (default=%default)")
    
    parser.add_option("-p", "--packet-size",
                  action="store", dest="packet_size", type=int, default=10,
                  help="Number of files per packet (default=%default)")
                  
    parser.add_option("-o", "--output",
                  action="store", dest="output_file", 
                  help="output plot filename (default = %default)",
                  default="health_check.pdf")
    
    parser.add_option("-v", "--verbose",
                  action="store_true", dest="verbose", default=True,
                  help="verbose mode [default=%default]")
    
                                
    (options, args) = parser.parse_args()
    
    
    if not options.input_images or len(args)!=0: 
    # args is the leftover positional arguments after all options have been processed
        parser.print_help()
        parser.error("Wrong number of arguments")

    if not os.path.exists(options.input_images):
        log.error ("Input image %s does not exist", options.input_images)
        sys.exit(0)
        
    try:
        #stack_frames(options.input_images, 'median', 1, 5, "/tmp/stack.fits",
        #             "/tmp/stack_sigma.fits", 10, 100000, 3)
        
        run_health_check(options.input_images, options.start_packet, 
                         options.end_packet, options.packet_size, 
                         options.window, options.output_file)
    except Exception, e:
        log.error("Some error while running Health-Check routine: %s"%str(e))
        sys.exit(0)
        
