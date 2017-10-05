#! /usr/bin/env python
#encoding:UTF-8

# Copyright (c) 2015 Jose M. Ibanez All rights reserved.
# Institute of Astrophysics of Andalusia, IAA-CSIC
#
# This file is part of PAPI
#
# PAPI is free software: you can redistribute it and/or modify
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


import math
import sys
import os
import shutil
import fileinput
from optparse import OptionParser


import numpy
import astropy.io.fits as fits

from misc.paLog import log
from misc.version import __version__
import misc.mef 

import datahandler


def createBPM(input_nlc, output_dir=None, joined=False):
    """
    Given a NLC MEF file, a BPM is produced reading NaNs from LINMAXi extensions. 
    The Bad Pixels will be saved as 1's on the BPM file created.
    
    Bernhard Dorner's comment:
    I suggest to take them from the nonlinearity file. The LINMAXi extensions 
    has uncorrectable bad pixels as nan. I've not cross-checked with dark and 
    flatfield origins, but from how they are determined, they should be similar 
    (the fractions are actually). Besides, these are the ones which cannot be 
    corrected, and therefore must be left out in the processing anyway.
    
    The structure of the NLC files is as follow:
    
    No.    Name         Type      Cards   Dimensions   Format
    0    PRIMARY     PrimaryHDU      46   ()              
    1    LINMAX1     ImageHDU        12   (2048, 2048)   float32   
    2    LINMAX2     ImageHDU        12   (2048, 2048)   float32   
    3    LINMAX3     ImageHDU        12   (2048, 2048)   float32   
    4    LINMAX4     ImageHDU        12   (2048, 2048)   float32   
    5    LINPOLY1    ImageHDU        12   (2048, 2048, 4)   float32   
    6    LINPOLY2    ImageHDU        12   (2048, 2048, 4)   float32   
    7    LINPOLY3    ImageHDU        12   (2048, 2048, 4)   float32   
    8    LINPOLY4    ImageHDU        12   (2048, 2048, 4)   float32

    Parameters
    ----------
    input_nlc: list
        MEF FITS files
        
    output_dir: str
        Dir name where BPM is created
        
    joined: bool
        Create a Single Extension FITS file in addition to the MEF.
        
    """
    
    new_bpm = output_dir + "/" + os.path.basename(input_nlc).replace(".fits","_BPM.fits")
    # Read NLC
    f = fits.open(input_nlc)
    if len(f) != 9:
        msg = "Error, NLC file has wrong number of extensions"
        log.error(msg)
        raise Exception(mgs)
    
    # Init the BPM FITS file
    hdulist = fits.HDUList()
    primaryHeader = fits.getheader(input_nlc)
    prihdu = fits.PrimaryHDU (data = None, header = primaryHeader)
    prihdu.header.set('IMAGETYP', 'MASTER_BPM')
    prihdu.header.set('PAPITYPE', 'MASTER_BPM')
    prihdu.header.rename_keyword('ID', 'NLC_FILE')
    prihdu.header.add_history("BPM file created from NLC file:%s"%input_nlc)
    
    hdulist.append(prihdu)
    
    # Create BPM arrays (NaNs will be 1's in BPM)
    for iExt in range(1, 5):
        # Check for compatitibility of extension naming
        try:
            extName = 'LINMAX%d'%iExt
            f[extName].header
        except KeyError, e:
            extName = 'SG%i_1' %iExt
        # create temporal array
        tmp = numpy.zeros([2048, 2048], dtype=numpy.int16)
        p = numpy.isnan(f[extName].data)
        # Set NaNs (bad pixels) to 1
        tmp[p] = 1 
        hdui = fits.ImageHDU(data = tmp, header = None)
        hdui.header.set('EXTNAME', 'Q%d' % iExt)
        hdui.header.set('DET_ID',  'SG%d' % iExt)
        
        
        hdulist.append(hdui)
    
    hdulist.writeto(new_bpm)
    hdulist.close(output_verify='ignore')
    
    # In addition to the MEF file, if selected, we create a SEF file
    # for the BPM.
    if joined:
        mef = misc.mef.MEF([new_bpm])
        res = mef.doJoin(".join.fits", output_dir=output_dir)[1]
        log.info("Joined file created %s"%res)
    
    return new_bpm

    
#################
# MAIN #
#################
if __name__ == "__main__":
    
    usage = "usage: %prog [options]"
    desc = """Create the BPM file from the NonLinearity correction MEF file.
The bad pixels will be saved as 1's""" 
    
    parser = OptionParser(usage, description=desc)
    
    parser.add_option("-s", "--source", type="str",
                  action="store", dest="source_file",
                  help="NonLinearity correction MEF file.")
    
    parser.add_option("-o", "--output_dir", type="str",
                  action="store", dest="output_dir", default="/tmp",
                  help="Output file to write the new BPM FITS file obtained.")
    
    parser.add_option("-j", "--joined_bpm",
                  action="store_true", dest="joined", default=False,
                  help="Create a SEF (single extension FITS) BPM in addition to MEF BPM [default=%default]")
    
    
    (options, args) = parser.parse_args()
    
    if len(sys.argv[1:]) < 1:
       parser.print_help()
       sys.exit(0)
        
    # args is the leftover positional arguments after all options have been processed 
    if not options.source_file or len(args)!=0: 
        parser.print_help()
        parser.error("incorrect number of arguments " )
        
    try:
        r = createBPM(options.source_file, 
                      output_dir=options.output_dir,
                      joined=options.joined)
        log.info("BPM file created: %s"%r)
    except Exception,e:
        log.error("Error, cannot create BPM: %s",str(e))
    
    
    sys.exit(0)