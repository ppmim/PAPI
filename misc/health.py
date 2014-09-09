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
#    PyRAF
#
# Created    : 29/10/2012    jmiguel@iaa.es -
# Last update: 03/06/2014    jmiguel@iaa.es - Some improvements
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
import astropy.io.fits as fits
import matplotlib
import matplotlib.pyplot as plt
            

# PAPI modules
from misc.paLog import log
from misc.utils import *
import misc.robust as robust
from misc.print_table import print_table



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

    papi_home = os.getenv("PAPI_HOME")
    if papi_home==None:
        raise Exception("Cannot find PAPI_HOME environment")
    else:
        irdr_path = papi_home + "/irdr/bin"
    
    # next file really is not used nor created
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
            
    
def run_health_check ( input_file, packet_size, f_from, f_to,  window='full-frame',
                        out_filename="/tmp/hc_out.pdf", temp_dir="/tmp"):
    """ 
    Takes a input catalog (ascii file) listing all the files (flat_fields) to 
    be used in the gain and noise computation. 
    
    Parameters
    ----------
    input_file: str
        Text file listing the files to use for health computation 
    packet_size: int
        Size of the packet; define how files are grouped
    f_from: int [0,...,N-1], where N = number of files
        File number inside the packet where start the computation 
    f_to: int [1,...,N], where N = number of files
        File number inside the packet where end the computation
    window: str
        Window in frame (Q1, Q2, Q3, Q4, full-detector, central)
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
      - print pretty stats in output, not only plots
      - Do dark subtraction to input files (flat-fields)
      - Do computations per channel and/or detector
      - Allow a custom coordinates for window definition
      
    """
    
    
    # Read the file list from input_file
    if os.path.exists(input_file):
        filelist = [line.replace( "\n", "") for line in fileinput.input(input_file)]
    else:
        raise Exception("Input file %s does not exist"%input_file)
    
    if window == 'Q1':
        x1, y1, x2, y2 = 10, 10, 2030, 2030
    elif window == 'Q2':
        x1, y1, x2, y2 = 2100, 10, 4080, 2100
    elif window == 'Q3':
        x1, y1, x2, y2 = 2100, 2100, 4080, 4080    
    elif window == 'Q4':
        x1, y1, x2, y2 = 10, 2100, 2100, 4080
    elif window == 'central':
        x1, y1, x2, y2 = 512, 512, 3584, 3584
    elif window == 'full': 
        x1, y1, x2, y2 = 10, 10, 4080, 4080 
    else:
        msg = "Wrong window definition."
        log.error(msg)
        raise Exception(msg)

    print "Selected area = [%d:%d, %d:%d]"%(x1,x2,y1,y2)
    print "Packet size = ", packet_size
    
    # check packet-range
    if not (f_from >= 0 and f_from < packet_size and 
            f_to >0 and f_to <= packet_size):
        raise Exception("Wrong values of packet file range")
      
    # check window-shape
    pf = fits.open(filelist[0])
    if len(pf)>1:
        msg = "MEF files are not supported. Expected a 4kx4k single image."
        log.error(msg)
        raise Exception(msg)
    
    if pf[0].shape!=(4096,4096):
        msg = "Expected a 4kx4k single FITS image."
        log.error(msg)
        raise Exception(msg)

    if not (x1 < pf[0].data.shape[0] and x2 < pf[0].data.shape[0] and 
        y1 < pf[0].data.shape[1] and y2 < pf[0].data.shape[1]):
        raise Exception("Wrong window definition; check image and window size")
     
    #tmp_file, _ = tempfile.mkstemp()
    tmp_file = temp_dir + "/packet.txt"
    n = 0
    for packet in grouper(packet_size, filelist):
        if len(packet)==packet_size and not (None in packet):
            listToFile(packet, tmp_file)
            try:
                stack_frames(tmp_file, 'median', f_from, f_to, 
                    temp_dir+"/stack%02d.fits"%n,
                    temp_dir+"/stack_sigma%02d.fits"%n, 
                    10, 100000, 3)    
            except Exception, e:
                raise e
            n = n + 1

            
    # Get stats from images
    stat_values = {}
    itime = numpy.zeros([n], dtype=numpy.float32)
    signal = numpy.zeros([n], dtype=numpy.float32)
    std = numpy.zeros([n], dtype=numpy.float32)
    kw_time = 'ITIME'
    for i in range(n):
        try:
            log.debug("Reading file %s"%(temp_dir+"/stack%02d.fits"%i)) 
            pf = fits.open(temp_dir + "/stack%02d.fits"%i)
            pf2 = fits.open(temp_dir + "/stack_sigma%02d.fits"%i)
        except Exception, e:
            log.error("Cannot open file %s"%(temp_dir+"/stack%02d.fits"%i))
            continue
        signal[i] = robust.mean(pf[0].data[x1:x2, y1:y2])
        std[i] = robust.mean(pf2[0].data[x1:x2, y1:y2])
        if kw_time in pf[0].header: 
            itime[i] = pf[0].header[kw_time]
        else:
            itime[i] = NaN

        stat_values[i] = [itime[i], signal[i], std[i]]
        pf.close()
        pf2.close()
        
    
    #print "Signal", signal
    #print "Std", std
    #print "Itimes", itime
    t = [
        ["Packet", "Signal  ", "Std     ", "ITime   "]
        ]
    i = 0        
    for row in zip(signal, std, itime):
        t.append([str(i), str(row[0]), str(row[1]), str(row[2])])
        i += 1
    print_table(t)
        
    # Gain
    # Compute the linear fit signal VS variance ( var = a*signal+b )
    # ==============================================================
    res = numpy.polyfit(signal, std**2, 1, None, True)
    a = res[0][1] # intercept
    b = res[0][0] # slope
    r = res[3][1] # regression coeff ??? not exactly
    
    #print "Coeffs =", res
    
    # Plot the Signal(x) VS Variance(y)
    pol = numpy.poly1d(res[0])
    plt.plot(signal, std**2, '.', signal, pol(signal), '-')
    #Note: gain = 1/slope
    #Note: intercept = (RON/gain)**2 => RON = gain*sqrt(intercept)
    gain = 1.0 / b
    ron = gain/math.sqrt(math.fabs(a))
    log.info("Gain = %f [e-/ADU]"%gain)
    log.info("RON = %f e-"%ron)

    plt.title("Poly fit: %f X + %f  r=%f Gain=%s RON=%s" %(b, a, r, 
                                                           gain, ron))
    plt.xlabel("Signal / ADU")
    plt.ylabel("Variance  / ADU")
    plt.savefig(out_filename)
    plt.show()
    
    # Linearity and full-well
    # Fit and Plot the ITime(x) VS Signal(y)
    # ======================================
    res = numpy.polyfit(itime, signal, 1, None, True)
    a = res[0][1] # intercept
    b = res[0][0] # slope
    r = res[3][1] # regression coeff ??? not exactly
    
    pol = numpy.poly1d(res[0])
    plt.clf()
    plt.plot(itime, signal, '.', itime, pol(itime), '-')
    
    full_well = numpy.max(signal[numpy.where( ((signal - pol(itime))/signal) < 0.05)]) 
    log.info("Full-well = %f"%full_well)
    
    plt.title("Poly fit: %f X + %f  r=%f Full-well=%s" %(b, a, r, full_well))
    plt.xlabel("ITime / s")
    plt.ylabel("Signal  / ADU")
    plt.savefig(out_filename+"_2.pdf")
    plt.show()
    

    # remove tmp files
    os.unlink(tmp_file)

    return out_filename, out_filename + "_2.pdf"

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
    desc = """Compute the Gain and Noise from a set of Flat images grouped in
packets and with increased level of Integration Time (ITIME). Flat files should
be dark corrected and non MEF 4kx4k files."""
    
    parser = OptionParser(usage, description=desc)
    
    parser.add_option("-i", "--input_images",
                  action="store", dest="input_images", 
                  help="Input image list (single non-integrated flats).")
    
    parser.add_option('-w', '--window',
                      type='choice',
                      action='store',
                      dest='window',
                      choices=['Q1', 'Q2', 'Q3', 'Q4', 'central', 'full'],
                      default='full',
                      help="Window/dectector to process: "
                      "Q1, Q2, Q3, Q4, central, full [default: %default]")

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
                  help="Output plot filename (default = %default)",
                  default="health_check.pdf")

    parser.add_option("-t", "--temp_dir",
                  action="store", dest="temp_dir", default="/tmp",
                  help="Temporal directory (default = %default)")
    
                                
    (options, args) = parser.parse_args()
    
    if not options.input_images or len(args)!=0: 
    # args is the leftover positional arguments after all options have been processed
        parser.print_help()
        parser.error("Wrong number of arguments")

    if not os.path.exists(options.input_images):
        log.error ("Input image %s does not exist", options.input_images)
        sys.exit(0)
    
    log.debug( 'Health-Check routines for PANIC')
    try:
        run_health_check(options.input_images, options.packet_size,
                         options.start_packet, options.end_packet, 
                         options.window, options.output_file,
                         options.temp_dir)
    except Exception, e:
        log.error("Some error while running Health-Check routine: %s"%str(e))
        sys.exit(0)
        
