#! /usr/bin/env python

# Copyright (c) 2010 Jose M. Ibanez. All rights reserved.
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
################################################################################
#
# PANICtool
#
# checkAstrometry.py
#
# Created    : 07/10/2010    jmiguel@iaa.es
#
################################################################################

# Import necessary modules
import os
import numpy as np
import tempfile
import subprocess
from optparse import OptionParser

def checkAstrometry(image, catalog="2mass"):
    """ 
        Evaluate the astrometric accuracy of an image comparing with an reference catalog (using the results from wcstool::immatch)
        
        'image',   must have a good enought WCS header
        'catalog', must be 2mass, usnob, usnoa, gsc and conresponding environtment variable must to be exported accordingly
    """
  
    if catalog=="2mass":
        my_cat = "tmc"
    elif catalog=="usnob":
        my_cat = "ub"
    elif catalog=="usnoa":
        my_cat = "ua"
    elif catalog=="gsc":
        my_cat = catalog
    else:
        print "\nCatalog %s not found, using 2MASS"%catalog
        my_cat = "tmc"
    
    # The output of 'immatch' will be saved to a temporary file, which will
    # then be parsed in order to extract the stars that were returned and the position errors.
    temp_fd, temp_path = tempfile.mkstemp()
    max_stars = 9999999 # allow for lots of stars to be extracted
    args = "/disk-a/caha/panic/SOFTWARE/wcstools-3.8.1/bin/immatch -c %s -h %d %s" \
            %(my_cat, max_stars, image)
    retcode = subprocess.call(args, shell = True, stdout = temp_fd)
    os.close(temp_fd)

    # We could use subprocess.check_call(), but this error message provides
    # more useful information
    if retcode != 0:    
            raise RuntimeError("The client 'immatch' encountered an error " \
                               "while querying the catalog server.")     

    # Now we load from the temp file the columns with needed data (7=x_err, 8=y_err, 9=radii_err)
    m=np.loadtxt(temp_path, comments="#", usecols=(7,8,9))
    if len(m)==0:
        print "Error checking astormetry, any star matched with catalog"
        return -1,-1,-1
    
    # Compute the Root-mean-squeare-deviation (http://en.wikipedia.org/wiki/Root_mean_square_error)
    x_rmsd = np.sqrt(np.sum(np.power(m[:,0],2))/len(m[:,0]))
    y_rmsd = np.sqrt(np.sum(np.power(m[:,1],2))/len(m[:,1]))
    r_rmsd = np.sqrt(np.sum(np.power(m[:,2],2))/len(m[:,2]))
    print "x_std=",np.std(m[:,0])
    print "y_std=",np.std(m[:,1])
    print "r_std=",np.std(m[:,2])
    
  
    return x_rmsd, y_rmsd, r_rmsd

################################################################################
# main
if __name__ == "__main__":
    # Get and check command-line options
    
    usage = "usage: %prog [options] arg1 arg2"
    parser = OptionParser(usage)
    
    parser.add_option("-v", "--verbose",
                  action="store_true", dest="verbose", default=True,
                  help="verbose mode [default]")
                  
    parser.add_option("-s", "--source",
                  action="store", dest="source",
                  help="Image with WCS to evaluate astrometrically. It has to be a fullpath file name (required)")
                  
    parser.add_option("-c", "--catalog",
                  action="store", dest="catalog", help="name of the catalog to use for evaluation (2MASS, USNO-B1,GSC2)")
    
   

    (options, args) = parser.parse_args()
    if not options.source or not options.catalog or len(args)!=0: # args is the leftover positional arguments after all options have been processed
        parser.print_help()
        parser.error("incorrect number of arguments " )
    
    result = checkAstrometry(options.source, options.catalog)
    print "\nRMS x= %f arcsec ,  y= %f arcsec, rad= %f arcsec\n"%(result[0],result[1], result[2])