#! /usr/bin/env python

""" Run IRAF.starfocus, get the best Focus for a focus sequence. """

# Author: Jose M. Ibanez (c) 2014
# Email: jmiguel@iaa.es
# License: GNU GPLv3

from pyraf import iraf
#import pyraf.iraf as iraf
#from iraf import obsutil
import re

from optparse import OptionParser
import sys


class IrafError(StandardError):
    """ Raised if some IRAF error happens """
    pass
  
def getBestFocusfromStarfocus(images, coord_file, log_file):
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
        starfocus.saturation = "INDEF"
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
    
    ############### Insert start
    parser.add_option("-d", "--data_file",
                  action="store", dest="data_file",
                  help="Output data file for analysis")
    
    parser.add_option("-t", "--target",
                  action="store", dest="target",
                  help="Object name for output data")
    ############### Insert end

    (options, args) = parser.parse_args()
    
    if len(sys.argv[1:])<1:
       parser.print_help()
       sys.exit(0)

    if not options.source_file or not options.coord_file:
        parser.print_help()
        parser.error("incorrent number of arguments")
        
    best_focus = getBestFocusfromStarfocus("@"+options.source_file, options.coord_file, options.log_file)
    print "Best Focus = ",best_focus
    
    ############### Insert start
    if options.data_file:
        # write log output to data file
        # parse backwards: look for best focus line and block with single results
        with open(options.log_file, "r") as f:
            f.seek (0, 2)           # Seek @ EOF
            fsize = f.tell()        # Get Size
            f.seek (max (fsize-2**15, 0), 0) # Set pos @ last chars
            lines = f.readlines()       # Read to end
        lines.reverse()
        fo = open(options.data_file, 'w')
        if options.target:
            obj = options.target
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
        print 'Data file written: %s' %options.data_file
    ############### Insert end
    
    sys.exit()
