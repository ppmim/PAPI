#! /usr/bin/env python
#encoding:UTF-8

# Copyright (c) 2013 Jose M. Ibanez All rights reserved.
# Institute of Astrophysics of Andalusia, IAA-CSIC
#
# This file is part of PAPI
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

###############################################################################
#
# PAPI (PANIC PIpeline)
#
# imtrim.py
#
# Created    : 21/01/2009    jmiguel@iaa.es
# Last update: 28/07/2014    jmiguel@iaa.es: remove import noao,ccdred, imred
#
###############################################################################

import sys
import os
import shutil
from optparse import OptionParser
import fileinput
import astropy.io.fits as fits


# Pyraf modules
from pyraf import iraf


# Interact with FITS files
from misc.paLog import log
from misc.version import __version__


                
def imgTrim(inputfile, outputfile=None, p_step=128):
    """
    Crop/cut the input image edges
       
    Parameters
    ----------
    inputfile: str
        Filename of file to be trimmed. It can be overwriten.
        
    outputfile: str
        Filename of the output file. If outputfile name exists and is not the 
        same than inputfile name, the funtion gives and error.
        It means, that inputfile name can be overwritten.
        For more details, see iraf.imcopy() function.

    p_step: int
        Highest value or row/column to be looked for border pixels
        
    
    Returns
    -------
        Whether success, return the crop values
          
    """ 

    if outputfile == None:
        # inputfile name is overwritten        
        outputfile = inputfile
    elif os.path.exists(outputfile) and outputfile != inputfile:
        raise Exception("Error, output file name already exists.")
        
    # clean double slash (//) due to problems in IRAF
    file = inputfile.replace("//", "/")
    outputfile = outputfile.replace("//", "/")
    
    log.debug("Start imgTrim ....")
    
    try:
        indata = fits.open(file, ignore_missing_end=True)
        if len(indata) > 1:
            raise Exception("MEF files currently not supported")
        indata[0].verify()
    except Exception, e:
        log.error('Could not open frame - something wrong with input data')
        raise e
    
    try:
        nx = indata[0].header['NAXIS1']
        ny = indata[0].header['NAXIS2']
        log.debug("NX=%d  NY=%d" % (nx, ny))
    except Exception,e:
        raise e
    
    indata.close()

    xmin = 0
    lasti = 0
    start = 1
    step = int(p_step)
    std = 0.0
    
    
    ####### 1st loop (xmin) #############################
    i = start
    while i<=nx:
        if lasti==i:
            std = 1.0
        else:
            #print "FILE =", file+"["+str(i)+",*]"
            std = float(iraf.imstat (
                images=file+"["+str(i)+",*]",
                fields='stddev', format='no', Stdout=1)[0])
        if (std!=0.0):
            if (i==1):
                xmin = 1
                break
            else:
                if (step==1):
                    xmin = i
                    break
                else:
                    lasti = i
                    step = step/2
                    i = i -2*step
        else:
            pass
        i+=step
        #print "DEBUG_xmin (i,std,lasti):" , i, std, lasti
         
    if (xmin==0):
        log.error("No data in file %s"%file)
        raise Exception("No data in file %s"%file)
    
    
    #### 2nd loop (xmax) ###############################  
    
    lasti = 0
    start = nx
    step = int(p_step)
    i = start
    while i>=1:
        if (lasti==i):
            std = 1.0
        else:
            std = float(iraf.imstat (
                images=file+"["+str(i)+",*]",
                fields='stddev',format='no',Stdout=1)[0])
        
        if (std!=0.0):
            if (i==nx):
                xmax = nx
                break
            else:
                if (step==1):
                    xmax = i
                    break
                else:
                    lasti = i
                    step = step/2
                    i = i+2*step
            
        else:
            pass
        i-=step    
        #print "DEBUG_xmax (i,std,lasti):" , i, std, lasti
        
    ##### 3rd loop (ymin) ####################  
    
    lasti = 0
    start = 1
    step = int(p_step)
    i = start
    while i<=ny:
        if (lasti==i):
            std = 1.0
        else:
            std = float(iraf.imstat (
                images=file+"[*,"+str(i)+"]",
                fields='stddev',format='no',Stdout=1)[0])
        
        if (std!=0.0):
            if (i==1):
                ymin = 1
                break
            else:
                if (step==1):
                    ymin = i
                    break
                else:
                    lasti = i
                    step = step/2
                    i=i-2*step
        else:
            pass
        i+=step
        #print "DEBUG_ymin (i,std,lasti):" , i, std, lasti
    

    ##### 4th loop (ymax) ####################  
    
    lasti = 0
    start = ny
    step = int(p_step)
    i = start
    while i >= 1:
        if (lasti == i):
            std = 1.0
        else:
            std = float(iraf.imstat (
                images=file+"[*,"+str(i)+"]",
                fields='stddev',format='no',Stdout=1)[0])
        
        if (std != 0.0):
            if (i == ny):
                ymax = ny
                break
            else:
                if (step == 1):
                    ymax = i
                    break
                else:
                    lasti = i
                    step = step / 2
                    i= i + 2 * step
        else:
            pass
        i-=step
        #print "DEBUG_ymax (i,std,lasti):" , i, std, lasti
    
    
    # Add a certain number or row/colums to avoid the noise
    xmin += 10
    ymin += 10
    xmax -= 10
    ymax -= 10
    
    #
    # Finaly, update (overwriting) the original images (.fits, 
    # .weight.fits, .objs.fits)
    # Note that iraf.imcopy makes the WCS update accordingly to the image crop,
    # and modify WCS keywords on header (but not always as we'd wish, e.g.,
    # remove CDi_j values equal to 0)
    log.debug("Trimming image %s [ %d : %d, %d : %d ]" % (file.replace("//", "/"), xmin, xmax, ymin, ymax))
    iraf.imcopy(input=file.replace("//", "/") + "[" + str(xmin) + ":" + str(xmax) + "," +
            str(ymin) + ":" + str(ymax) + "]", output=outputfile)
    
    fits.setval(outputfile, keyword='HISTORY', value='Image trimmed', ext=0)
    fits.setval(outputfile, keyword='PAPIVERS', value=__version__, 
                    comment='PANIC Pipeline version', ext=0)

    # 
    # Look for weight and objs images
    #
    try:
        ima_sec = file.replace(".fits", ".weight.fits").replace("//", "/")
        if os.path.exists(ima_sec):
            log.debug("Trimming image %s [ %d : %d, %d : %d ]" % (ima_sec, xmin, xmax, ymin, ymax))
            iraf.imcopy(input=ima_sec + "[" + str(xmin) + ":" + str(xmax) + "," +
                str(ymin) + ":" + str(ymax) + "]", output=ima_sec)
        ima_objs = file.replace(".fits", "objs.fits").replace("//", "/")
        if os.path.exists( ima_objs ):
            log.debug("Trimming image %s [ %d : %d, %d : %d ]" % (ima_objs, xmin, xmax, ymin, ymax))
            iraf.imcopy(input=ima_objs + "[" + str(xmin) + ":" + str(xmax) + "," +
                str(ymin) + ":" + str(ymax) + "]", output=ima_objs)
                
    except Exception,e:
        log.debug("Some error trimming image %s . Probaby input image is wrong." % ima_sec)

    log.debug("....End of imgTrim --> XMIN= %d YMIN=%d XMAX= %d YMAX=%d"
            %(xmin, ymin, xmax, ymax))
    
    return (xmin, ymin, xmax, ymax)
        
###############################################################################
# main
if __name__ == "__main__":
        
    
    # Get and check command-line options
    usage = "usage: %prog [options]  arg2 ..."
    desc = """Performs image trimming looking for a constant frame around the \
image.
"""
    parser = OptionParser(usage, description=desc)
    
    parser.add_option("-i", "--input_image",
                  action="store", dest="input_image", 
                  help="Input image to trim")
    
    parser.add_option("-o", "--output",
                  action="store", dest="output_image", 
                  help="Output file for trimmed image")
    
    parser.add_option("-s", "--step",
                  action="store", dest="step", default=128, type=int,
                  help="Max. step to look for border pixels [default=%default]")
    
                                
    (options, args) = parser.parse_args()
    
    
    if options.input_image and os.path.exists(options.input_image) and \
    options.output_image:
        try:    
            imgTrim(options.input_image, 
                    options.output_image,
                    options.step)
        except Exception,e:
            log.error("Error while trimming image %s"%options.input_image)
            log.error("%s"%e)
            raise e
    else:
        parser.error("Input or output file does not exist.")
