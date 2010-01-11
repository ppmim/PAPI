#!/usr/bin/env python

################################################################################
#
# simple_reduce
#
# simple_reduce.py
#
# Last update 30/04/2008
#
################################################################################


    ####################################################################
    #                                                                  #
    #                   This is where it all starts.		       #
    #         The body of this program is at the end of this file      #
    #								       #
    ####################################################################


"""
   Main routines for the simple automatic data reduction
"""

_version     = "1.0.0"
_date        = "30-04-2008"
_author      = "Jose M. Ibanez (jmiguel@iaa.es)"


_minversion_numpy   = "1.0.1"
_minversion_pyfits  = "1.1"
_minversion_pyraf   = "1.4"
_minversion_biggles = "1.6.4"


################################################################################
# Import necessary modules

import getopt
import sys
import os
import fileUtils

import utils
import fits

# Pyraf modules
from pyraf import iraf
from iraf import noao
from iraf import imred
from iraf import ccdred

# Interact with FITS files
import pyfits

# Import Pyro core
import Pyro.core
import Pyro.naming



###############################################################################
def run ( args ):

    #Initialize variables
    source_frame=""
    master_dark=""
    master_flat=""
    frame_out=""
    
    # Get and check command-line options
    try:
        opts, args = getopt.getopt(args, "-s -d -f -o", ['source=','dark=','flat=','out='])
    except getopt.GetoptError:
        # print help information and exit:
        usage()
        sys.exit(2)

    #if (opts==[]):
    #    usage()
    #    sys.exit(2)
    
    print opts
    
    for option, parameter in opts:
        if option in ("-s", "--source"):
            source_frame = parameter
            print "SOURCE frame =", source_frame
        if option in ("-d", "--dark"):
            master_dark = parameter
            print "master_dark =", master_dark
        if option in ("-f", "--flat"):
            master_flat = parameter
        if option in ("-o", "--out"):
            frame_out = parameter
            
            
    if  source_frame=="" or master_dark=="" or master_flat=="" or frame_out=="":
        usage()
        sys.exit(2)

    if simple_reduce ( source_frame, master_dark, master_flat, frame_out, True )==0:
        print "Simple reduce was sucessfull !! "
         
     

######################################################################################
def usage ():
     print "Unknown command line parameter. Required parameters are : "
     print "-s / --source=		Source data frame in"
     print "-d / --dark=                Master DARK frame"
     print "-f / --flat=                Master FLAT frame"
     print "-o / --out=                 Output frame"

     
######################################################################################

def simple_reduce(frame_in, master_dark, master_flat, frame_out, bdPixel):
    
  """Subtract sky from an image using only its own background for sky computing
  """
  
  
  print('Simple reduce started ...')
  t=utils.clock()
  t.tic()
  local="/usr/local/bin/"
  
  # Get the user-defined list of science frames
  #framelist = config.options['sciencelist'].value
  #frame=framelist[0]  # For now, only the first science list image
  #frameoutD=frame.replace(".fits","_D.fits")
  #frameoutD=frameoutD.replace(config.options['inpath'].value, config.options['outpath'].value)
  
  frameoutD=frame_in.replace(".fits","_D.fits")
  #frameoutD=frameoutD.replace(config.options['inpath'].value, config.options['outpath'].value)
  
  fc=fits.FitsFile(frame_in)
  ff=fits.FitsFile(master_flat)
  if ( fc.filter!=ff.filter ):
      messageLog.put_error("Error reducing frame %s. Filter missmatch with calibration frames" %frame_in)
      return -1
  
  # STEP 1: Subtract masterdark
  iraf.imarith(operand1=frame_in,
               operand2=master_dark,
               op='-',
               result=frameoutD,
               verbose='yes'
               )
  
  # STEP 2: Divide by normalized master FlatField
  frameoutF=frameoutD.replace("_D.fits","_D_F.fits")
  iraf.imarith(operand1=frameoutD,
               operand2=master_flat,
               op='/',
               result=frameoutF,
               verbose='yes'
               )
  
  # STEP 3: Compute sky background with only one frame (itself)
  frame_out=frameoutF.replace("_D_F.fits","_D_F_S.fits")
  sex(frameoutF, frame_out )
  
  
  # STEP 4: Apply pixel mask
  if bdPixel:
      applyPixelMask( frame_out, "/disk-a/caha/panic/DATA/data_alh/output/BadPixMask-20080409192131.pl")
      
      
  # STEP 5: Cleanup
  fileUtils.removefiles(frameoutD, frameoutF)
  
  print("Simple reduce finished. %s" % t.tac() )
  
  return 0


################################################################################
def applyPixelMask( frame_in, mask ):
  """
    Apply a previously created pixel mask to a frame
  """

  print ('Start applyPixelMask; applying mask %s to frame %s' %(mask,frame_in))
  t=utils.clock()
  t.tic()
  
  base=os.path.dirname(frame_in)
  iraf.chdir(base)

  iraf.fixpix(images=frame_in,
              masks=mask,
              verbose='yes'
              )

  print ( "applyPixelMask finished. %s" % t.tac() )

################################################################################

def sex (image_in, image_out, params=None, sconfig=None):
  
  """
    Set SExtractor parameters and configuration and run.
    
    image_in: Detection (and measurement) image
    image_out: (optional) Measurement image
    params: (optional) extra output parameters (Parameters object)
    config: (optional) configuration options (Config object)
    
    NOTE: If either params or config is not specified, the defaults will be
    used
  """

  sextractor_dir="/disk-a/caha/panic/DATA/test1/bin/"


  #cmd_line=config.options['sextractor_dir'].value+"/sex -c sex.conf" + image_in
  #os.system("d")
  params="-PARAMETERS_NAME " + sextractor_dir + "/default.param" + \
          " -FILTER N " + \
          " -FITS_UNSIGNED Y " + \
          " -CATALOG_NAME catalog.cat " + \
          " -DETECT_MINAREA  5" + \
          " -DETECT_THRESH 5.0 " + \
          " -CHECKIMAGE_TYPE -BACKGROUND " + \
          " -CHECKIMAGE_NAME " +  image_out

  command = sextractor_dir + "/sex -c sex.conf1 " + image_in +" "+ params
  print "SEX COMMAND = " , command 

  (child_stding, child_stdout, child_stderr ) = os.popen3( command )
  #child_stding.write("-c sex.conf " + image_in)
  output = child_stdout.readline()
  while output!="":
    print output
    output = child_stdout.readline()


################################################################################      
if __name__=="__main__":

    run(sys.argv[1:])

        
################################################################################
