#! /usr/bin/env python

""" Run IRAF.psfmeasure, get field FWHM of given stars """

# Author: Jose M. Ibanez (c) 2014
# Email: jmiguel@iaa.es
# License: GNU GPLv3

from pyraf import iraf
#import pyraf.iraf as iraf
#from iraf import obsutil
import re

from optparse import OptionParser
import sys


class IrafError(Exception):
    """ Raised if some IRAF error happens """
    pass
  
def getAverageFWHMfromPsfmeasure(images, coord_file, log_file):
    """
    Calculate the average Full Width Half Max for the objects in image
    at the coords specified in coord_file calling iraf.obsutil.psfmeasure.

    The coordinates in coord_file should be in the same world coordiantes
    as the WCS applied to the image.
    
    Exam. coord_file:
    1024  1024
    
    """
    try:
        iraf.noao(_doprint=0)
        iraf.obsutil(_doprint=0)
        iraf.unlearn("psfmeasure")
      
        psfmeasure = iraf.obsutil.psfmeasure
        # setup all paramaters
        psfmeasure.coords = "mark1"
        psfmeasure.wcs = "world"
        psfmeasure.display = "no"
        psfmeasure.size = "GFWHM"
        psfmeasure.radius = 15
        psfmeasure.sbuffer = 5 
        psfmeasure.swidth = 5
        psfmeasure.iterations = 3
        psfmeasure.saturation = "INDEF"
        psfmeasure.ignore_sat = "no"
        psfmeasure.imagecur = coord_file
        #psfmeasure.graphcur = 'q.txt' #'myfile.mc' #'/dev/null' #file that is empty by definition
        psfmeasure.logfile = log_file
        res = psfmeasure(images, Stdout=1)[-1] # get last linet of output
        numMatch = re.compile(r'(\d+(\.\d+)?)')
        match = numMatch.search(res)
        
        return float (match.group(1))
    
    except Exception as e:
        print("Error running IRAF.psfmeasure: %s"%str(e))
        


if __name__ == "__main__":

    
    usage = "usage: %prog [options] "
    desc = """Run IRAF.psfmeasure and get field FWHM of given stars"""

    parser = OptionParser(usage, description=desc)
    
                  
    parser.add_option("-s", "--source",
                  action="store", dest="source_file",
                  help="Source file listing the filenames of input images.")
    
    parser.add_option("-c", "--coordiantes",
                  action="store", dest="coord_file",
                  help="Coordinates file listing the x,y coordiantes of stars in input images.")
    
    parser.add_option("-o", "--output_log",
                  action="store", dest="log_file",default="psfmeasure.log", 
                  help="Output log file generated [default: %default]")
    
    (options, args) = parser.parse_args()
    
    if len(sys.argv[1:])<1:
       parser.print_help()
       sys.exit(0)

    if not options.source_file or not options.coord_file:
        parser.print_help()
        parser.error("incorrent number of arguments")
        
    fwhm = getAverageFWHMfromPsfmeasure("@"+options.source_file, options.coord_file, options.log_file)
    print("Mean FWHM = ", fwhm)
    sys.exit()
