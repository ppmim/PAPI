#! /usr/bin/env python

""" Run IRAF.starfocus, get the best Focus for a focus sequence. """

# Author: Jose M. Ibanez (c) 2014
# Email: jmiguel@iaa.es
# License: GNU GPLv3

import os
import re
from optparse import OptionParser
import sys
import fileinput

# To avoid conflict with QL-Focus_Evaluation
try:
    if os.path.exists(os.path.expanduser("~") + "/iraf/focus_seq.txt"):
        os.unlink(os.path.expanduser("~") + "/iraf/focus_seq.txt")
        print("Deleted file ~/iraf/focus_seq.txt")
except Exception,e:
    print("Error, cannot delete ~/iraf/focus_seq.txt")

    

import astropy.io.fits as fits
import numpy as np
import matplotlib.pyplot as plt
from pyraf import iraf



class IrafError(StandardError):
    """ Raised if some IRAF error happens """
    pass
  
def getBestFocusfromStarfocus(images, coord_file, log_file):
    """
    Calculate the average Full Width Half Max for the objects in image
    at the coords specified in coord_file calling iraf.obsutil.psfmeasure.

    Paramaters
    ----------
    images: str
        Filename of the file listing the files to be analyzed. It must start
        with '@' for iraf formatting.
        
    coord_file: str
        Filename of the file with the coordinages (x, y) of the stars to be
        analyzed by iraf.starfocus.
    
    log_file: str
        Filename of the log file where iraf.starfocus will write the results.
       
      
    The coordinates in coord_file should be in the same world coordiantes
    as the WCS applied to the image.
    
    Exam. coord_file:
    1024  1024
    
    
    Returns
    -------
    The best focus obtained by iraf.starfocus.
    
    
    """
    
    # Read NCOADDS of the images to set the SATURATION limit
    if os.path.isfile(images[1:]):
        files = [line.replace("\n", "").replace('//','/')
                     for line in fileinput.input(images[1:])]
        with fits.open(files[0]) as hdu:
            if 'NCOADDS' in hdu[0].header:
                satur_level = hdu[0].header['NCOADDS'] * 50000
            else:
                satur_level = 50000
                
    print "SATUR_LEVEL =", satur_level
    
    # Config and launch the iraf.starfocus task
    try:
        iraf.noao(_doprint=0)
        iraf.obsutil(_doprint=0)
        iraf.unlearn("starfocus")
      
        starfocus = iraf.obsutil.starfocus
        starfocus.focus = "T_FOCUS"
        starfocus.fstep = "" 
        starfocus.nexposures = 1 
        starfocus.coords = "mark1"
        starfocus.wcs = "world" 
        starfocus.size = "GFWHM"
        starfocus.scale = 1
        starfocus.radius = 10
        starfocus.sbuffer = 10 
        starfocus.swidth= 10
        starfocus.saturation = satur_level
        #starfocus.saturation = "INDEF"
        starfocus.ignore_sat = "no"
        starfocus.imagecur = coord_file
        starfocus.display = "no"
        starfocus.frame = 1
        #starfocus.graphcur = "/dev/null" 
        starfocus.logfile = log_file
        res = starfocus(images, Stdout=1)[-1] # get last linet of output
        numMatch = re.compile(r'(\d+(\.\d+)?)')
        match = numMatch.search(res)
        
        return float (match.group(1))
    
    except Exception,e:
        print "Error running IRAF.starfocus: %s"%str(e)
        
def writeDataFile(log_file, data_file, target):
    """
    Read iraf.starfocus log file and write a data file to be used
    later for the Tilt analysis (p_50_tiltcheck.py).
    """
    
    if data_file:
        # write log output to data file
        # parse backwards: look for best focus line and block with single results
        with open(log_file, "r") as f:
            f.seek (0, 2)           # Seek @ EOF
            fsize = f.tell()        # Get Size
            f.seek (max (fsize-2**15, 0), 0) # Set pos @ last chars
            lines = f.readlines()       # Read to end
        lines.reverse()
        fo = open(data_file, 'w')
        if target:
            obj = target
        else:
            print 'WARNING: Object name not provided'
            obj = 'Unknowm'
        fo.write('# Object: %s\n' %obj)
        while not lines[0].strip().startswith('Average'):
            lines.pop(0)
        line = lines.pop(0)
        fo.write('#%s' %line)
        while not lines[0].strip().startswith('Best'):
            lines.pop(0)
        k = lines.index('\n')
        for i in range(k):
            if lines[k -i -1].strip().startswith('Best'):
                fo.write(lines[k - i - 1])
        fo.close()
        print 'Data file written: %s' %data_file
    else:
        print 'Error, no data file given'
        
def readStarfocus(log_file):
    """
    Read the results from the iraf.starfocus log file and compute
    the best focus for that execution. Only non-saturated data
    are read.
    """
    
    # Read the last lines 
    with open(log_file, "r") as f:
        f.seek (0, 2)           # Seek @ EOF
        fsize = f.tell()        # Get Size
        f.seek (max (fsize-2**15, 0), 0) # Set pos @ last chars
        lines = f.readlines()       # Read to end
    
    lines.reverse()
    my_lines = []
    
    # Look for the last execution (It must start with
    # 'NOAO/IRAF') 
    while not lines[0].strip().startswith('NOAO/IRAF'):
        my_lines.append(lines[0])
        lines.pop(0)

    my_lines.reverse()
    data = []
    # start to read the right columns
    for line in my_lines:
        if line.strip().startswith('Image'):
            # Heading line
            continue
        elif line.strip().startswith('Best'):
            # End reading
            break
        elif len(line.split()) == 0: 
            # Blank line
            continue
        else:
            # Read columns, only non-saturated data
            if line.split()[-1] != '*':
                if line.split()[0].startswith("/"):
                    data.append(line.split()[1:9])
                else:
                    data.append(line.split()[0:8])
     
    return data

def getBestFocus(data, output_file):
    """
    Fit the input data read from starfocus log file,
    to a parabola and find out the minimin (best focus estimation).
    
    Parameters
    ----------
    data: list
        A list with N rows x M (8) columns with the next correspondence:
        Column    Line     Mag   Focus   MFWHM Beta   Ellip      PA  
        
    output_file: str
        Filename of out plot with the fitting computed.
        
    Returns
    -------
    If success, returns the best focus computed.
    
    """
    
    d = np.array(data, dtype=np.float32)
    good_focus_values = d[:, 3] # focus
    fwhm_values = d[:, 4] # PSF-value (MFWHM, GFWHM, FWHM, ...)
    
    
    print "\n---------"
    print "N_POINTS: ", len(fwhm_values)
    print "----------\n"
    
    
    m_foc = good_focus_values.mean()
    good_focus_values = good_focus_values - m_foc
    z = np.polyfit(good_focus_values, fwhm_values, 2)
    print "Fit = %s  \n"%str(z)
    # Note that poly1d returns polynomials coefficients, in increasing powers !
    pol = np.poly1d(z)
    
    xp = np.linspace(np.min(good_focus_values ) - 0.5, 
                     np.max(good_focus_values ) + 0.5, 500) # number or points to interpolate

    # best focus is derivative of parabola = 0
    # but check if it is correctly curved
    if pol[2] < 0:
        print "ERROR: Parabola fit unusable!"
    best_focus = - pol[1] / (2. * pol[2])
    if best_focus + m_foc < 10 or best_focus + m_foc > 30:
        print "ERROR: Best focus out of range!"
    #print "BEST_FOCUS=", best_focus + m_foc
    
    # Plotting
    plt.plot(good_focus_values + m_foc, fwhm_values, '.')
    plt.plot(xp + m_foc, pol(xp), '-')
    plt.axvline(best_focus + m_foc, ls='--', c='r')
    plt.title("Focus serie - Fit: %f X^2 + %f X + %f\n Best Focus=%6.3f mm" 
        %(pol[2], pol[1], pol[0], best_focus + m_foc))
    
    plt.xlabel("T-FOCUS (mm)")
    plt.ylabel("FWHM (pixels)")
    plt.xlim(np.min(good_focus_values + m_foc) - 0.1, np.max(good_focus_values + m_foc) + 0.1)
    plt.ylim(np.min(fwhm_values) - 1, np.max(fwhm_values) + 1 )
    if pol[2] < 0:
        plt.figtext(0.5, 0.5, 'ERROR: Parabola fit unusable!', size='x-large', color='r', weight='bold', ha='center', va='bottom')
    if best_focus + m_foc < 10 or best_focus + m_foc > 30:
        plt.figtext(0.5, 0.5, 'ERROR: Best focus out of range!', size='x-large', color='r', weight='bold', ha='center', va='top')
    plt.grid()
    
    
    plt.savefig(output_file)
    
    print 'Image saved: ', output_file
    show = True
    if show:
        plt.show(block=True)
    
    # Print out Values
    #for idx, foc_value in enumerate(good_focus_values):
    #    sys.stdout.write("\nFoc. value: %s   -->  FWHM: %s"%(foc_value, fwhm_values[idx]))
    
    return (best_focus + m_foc)

if __name__ == "__main__":
    
    usage = "usage: %prog [options] "
    desc = """Run IRAF.starfocus for a focus sequecen and return the best focus"""

    parser = OptionParser(usage, description=desc)
    
                  
    parser.add_option("-s", "--source",
                  action="store", dest="source_file",
                  help="Source file listing the filenames of input images.")
    
    parser.add_option("-c", "--coordiantes",
                  action="store", dest="coord_file",
                  help="Coordinates file listing the x,y coordiantes of stars in input images.")
    
    parser.add_option("-o", "--output_log",
                  action="store", dest="log_file",default="starfocus.log", 
                  help="Output log file generated [default: %default]")
    
    parser.add_option("-d", "--data_file",
                  action="store", dest="data_file",
                  help="Output data file for analysis")
    
    parser.add_option("-t", "--target",
                  action="store", dest="target",
                  help="Object name for output data")

    (options, args) = parser.parse_args()
    
    if len(sys.argv[1:])<1:
       parser.print_help()
       sys.exit(0)

    if not options.source_file or not options.coord_file:
        parser.print_help()
        parser.error("incorrent number of arguments")
        
    best_focus = getBestFocusfromStarfocus("@" + options.source_file, 
                                           options.coord_file, 
                                           options.log_file)
    print "\n\nBest Focus (IRAF)= ",best_focus
    
    # Read log file and write values into data file for the Tilt analysis.
    writeDataFile(options.log_file, options.data_file, options.target)
    
    # Compute our own BEST_FOCUS value and plot the fittting
    print "Now, our own fitting...\n"
    data = readStarfocus(options.log_file)
    
    # I do not know why if IRF window is open, then matplotlib crash with a SF !!!
    print "Please, close IRAF window and press any key to continue..."
    raw_input('')
    my_best_focus = getBestFocus(data, "starfocus.pdf")
    print "\n\nBest Focus (OWN):", my_best_focus
    
    print 'The end !'
    sys.exit(0)
