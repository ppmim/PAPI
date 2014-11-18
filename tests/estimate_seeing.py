#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import numpy as np
import pyraf.iraf as iraf
import astropy.io.fits as fits  
import argparse
import sys
import os
import astropy.wcs.wcs as wcs 

os.environ['PYRAF_NO_DISPLAY'] = '0'

def calculate_seeing():
    """ Program to estimate the seeing from an image and a list of estimates
    for the positions of stars. After calculating the seeing, some of the 
    stars might get recalculated centers. The list will be updated with the 
    """
    
    output = "output.txt"
    ignore = "ignore.txt"

    iraf.noao()
    iraf.obsutil()
    iraf.module.psfmeasure("@/tmp/focus.txt", coords = "mark1", size = "GFWHM", 
                           sbuffer = 10, swidth=10, radius=10,
                           satura = 55000, ignore_sat="no", scale=0.45,
                           imagecur = "", display = "yes", 
                           graphcur = "", 
                           wcs ="logical")
        
        
    print "[psfmeasure] Y ahora que .....???"


def calculate_best_focus():
    """ 
    Program to estimate the best focus from a sequence of images taken 
    at diferent focus values. 
    The routine uses iraf.obsutil.starfocus() task.
    """
    
    output = "output.txt"
    ignore = "ignore.txt"

    iraf.noao()
    iraf.obsutil()
    iraf.module.starfocus("@/tmp/focus.txt", 
                            focus = "T_FOCUS", fstep = "", nexposures = "1",
                            coords = "mark1", wcs = "logical", 
                            size = "MFWHM", scale = 1,  
                            sbuffer = 10, swidth=10, radius=5,
                            satura = 55000, ignore_sat="no",
                            imagecur = "", display = "yes", frame = 1,
                            graphcur = "", logfile = "starfocus.log"
                            )
        
        
    print "[starfocus] Y ahora que .....???"
############################################################################


# Create parser
parser = argparse.ArgumentParser(description='Program to find stars in an image. ')

# Add necessary arguments to parser
parser.add_argument("input", metavar='input', action='store', help='list of ' +\
                    'input images for which to estimate the FWHM.',
                    nargs="+", type=str)
parser.add_argument("--cat", metavar='cat', action='store', dest="cat",
                    help='list of ' +\
                    'catalogs of the position of stars for the input images.', \
                    nargs=1, type=str)
parser.add_argument("--wcs", metavar="wcs_in", action="store", dest="wcs", 
                    default="logical",
                    help = "System in which the input coordinates are. Can "
                    "be 'logical', 'tv', 'physical' and 'world' accordind to "
                    "IRAF. If your coordinates are, for example, in RA and DEC "
                    "you should provide 'world'. Default: logical.")


def main(arguments = None):
  # Pass arguments to variable args
  if arguments == None:
      arguments = sys.argv[1:]

  args = parser.parse_args(arguments)
  
  if len(args.input) != len(args.cat):
    sys.exit("\n\n number of star catalogues and input images do not coincide \n ")      
      
  calculate_seeing(args)  
  return None    
     
if __name__ == "__main__":
    calculate_best_focus()
