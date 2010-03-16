################################################################################
#
# PAPI (PANIC PIpeline)
#
# makeobjmask.py
#
# Created    : 04/08/2008    jmiguel@iaa.es
# Last update: 04/08/2009    jmiguel@iaa.es
#
################################################################################
#
# Create object masks (SExtractor OBJECTS images) for a list of FITS images.
# Expects the command "sex" (SExtractor Version 2+) in path.  If weight maps
# exist they will be used (assume weight map filename given by replacing .fits
# with .weight.fits).
#

################################################################################
# Import necessary modules

import getopt
import sys
import os
import glob
import fileinput

import misc.utils as utils
from misc.paLog import log


#-----------------------------------------------------------------------

def usage():
    print ''
    print 'NAME'
    print '       makeobjmask.py - Object mask creation with SExtractor\n'
    print 'SYNOPSIS'
    print '       makeobjmask.py [options] -f file.list\n'
    print 'DESCRIPTION'
    print '       Create object mask (SExtractor OBJECTS images) for a list of FITS images'
    print '       If weight maps exist they will be used.'
    print ' '
    print 'OPTIONS'
    print '       -v : verbose debugging output\n'
    print '       -m --minarea    5     SExtractor DETECT_MINAREA (int)'
    print '       -t --threshold  2.0   SExtractor DETECT_THRESH  (float)'
    print 'VERSION'
    print '       04 August 2009'
    print ''
    raise SystemExit

#-----------------------------------------------------------------------

                
def makeObjMask (inputfile, minarea=5,  threshold=2.0, outputfile="/tmp/out.txt"):
    """DESCRIPTION
                Create an object mask of the inputfile/s based on SExtractor
           
           INPUTS
                inputfile       it can be a regular expresion or a filename of a filelist
                
           OPTIONAL INPUTS
                
                outputfile      filename of file list with the object file/s created by SExtractor
                
                p_min_area      SExtractor DETECT_MINAREA (minimun object area)
                           
                p_mask_thresh   SExtractor DETECT_THRESH
                
                remove_objs     Flat to indicate if remove or not tha .obj files created after masks creation
                
           OUTPUTS
                masks           Objects mask files ending with '.objs' suffix
              
      """
         
    # Some other settings
    irdr_basedir=''
    try:
        irdr_basedir=os.environ['IRDR_BASEDIR']
    except KeyError:
        log.error("Please, setenv IRDR_BASEDIR")
        sys.exit(0)
        
    sex_config=irdr_basedir+"/src/config/default.sex"        
    if not os.path.exists(sex_config):      # check whether input file exists
        log.error( 'File %s does not exist', sex_config)
        sys.exit(1) 
           
    files = []       
    # Check if inputfile is a filename of a file list OR a regular expresion        
    if os.path.isfile(inputfile):
        files=[line.replace( "\n", "") for line in fileinput.input(inputfile)]
    else:
        files = glob.glob(inputfile)
        files.sort()
        
    f_out = open(outputfile,"w")
    
    for fn in files:
        print 'File : ', fn
        wfn = fn
        wfn=wfn.replace(".fits",".weight.fits")
        
        if not os.path.exists(fn):      # check whether input file exists
            log.error( 'File %s does not exist', fn)
            sys.exit(2) 
        
        if os.path.exists(wfn):
            wstr=" -WEIGHT_IMAGE " + wfn +" -WEIGHT_TYPE MAP_WEIGHT "
        else: wstr=""
        
        cmd="sex %s -c %s -FITS_UNSIGNED Y %s -DETECT_MINAREA %s  -DETECT_THRESH %s -PIXEL_SCALE 0.45 -GAIN 4.15 -CHECKIMAGE_TYPE OBJECTS \
        -CHECKIMAGE_NAME %s" % (fn, sex_config, wstr, str(minarea), str(threshold), fn+".objs")  
        
        #print cmd
        e=utils.runCmd( cmd )
        if e==0:
            log.debug("Some error while running command %s", cmd)
        else:
            log.debug("Succesful ending of makeObjMask")
            f_out.write(fn+".objs"+"\n")
            
    f_out.close()
    print 'ending makeObjMask....'
    
################################################################################
# main
if __name__ == "__main__":
    print 'Start makeObjMask....'
    
    # Read command line parameters
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'f:m:t:o:v',['file=','minarea=','threshold=','out=','verbose'])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    
    nargs = len(sys.argv[1:])
    nopts = len(opts)
      
    
    minarea = 5
    threshold = 2.0
    verbose= False
    inputfile=''
    outputfile='/tmp/objs.txt'
            
    for option, par in opts:
        if option in ('-v','--verbose'):      # verbose debugging output
            verbose = True
            print "Verbose true"
        if option in ("-f", "--file"):
            inputfile=par
            print "inputfile=", inputfile
        if option in ("-m", "--minarea"):
            isomin=par
            print "minarea=", minarea
        if option in ("-t", "--threshold"):
            ellipmax=par
            print "threshold=",threshold
        if option in ("-o", "--out"):
            outputfile=par
            print "out_file=",outputfile
            

    makeObjMask( inputfile, minarea, threshold, outputfile)
    
