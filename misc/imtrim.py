#!/usr/bin/env python
################################################################################
#
# PAPI (PANIC PIpeline)
#
# imtrim.py
#
# Created    : 21/01/2009    jmiguel@iaa.es
# Last update: 10/09/2009    jmiguel@iaa.es
#
# TODO : 
#    - update AR,DEC after image trimming !!!!
################################################################################

# Import necessary modules
import os
import sys
import getopt

# Pyraf modules
from pyraf import iraf
from iraf import noao
from iraf import imred
from iraf import ccdred

# Interact with FITS files
import pyfits
from misc.paLog import log



#-----------------------------------------------------------------------

def usage():
    print ''
    print 'NAME'
    print '       imtrim.py - Cut/crop edges of the input image\n'
    print 'SYNOPSIS'
    print '       imtrim.py -i file.fits\n'
    print 'DESCRIPTION'
    print '       Crop/cut the edges of the input image'
    print ' '
    print 'OPTIONS'
    print '       -v : verbose debugging output\n'
    print 'VERSION'
    print '       10 September 2009'
    print ''
    raise SystemExit

#-----------------------------------------------------------------------

                
def imgTrim (inputfile):
  
        """DESCRIPTION
                Crop/cut the input image edges
           
           INPUTS
                -i inputfile       input FITS file to trim
                
           OPTIONAL INPUTS
                
           OUTPUTS
                                The trimmed input image
              
        """ 
        ### Start the script #####
        file = inputfile.replace("//","/") # clean double slash (//) due to problems in IRAF
        
        log.debug("Start imgTrim ....")
        
        
        
        try:
            indata = pyfits.open(file, ignore_missing_end=True)
            indata[0].verify()
        except:
            print('Could not open frame - something wrong with input data')
            raise
        try:
            nx = indata[0].header['NAXIS1']
            ny = indata[0].header['NAXIS2']
        except:
            raise
        
                    
        xmin = 0
        lasti = 0
        start = 1
        step = 128
        std = 0.0
        
        ####### 1st loop #############################
        i = start
        while i<=nx:
            if lasti==i:
                std = 1.0
            else:
                #print "FILE =", file+"["+str(i)+",*]"
                std = float(iraf.imstat (
                    images=file+"["+str(i)+",*]",
                    fields='stddev', format='no', Stdout=1)[0])
            if (std!=0.0) :
                if (i==1):
                    xmin = 1
                    break
                else:
                    if (step==1):
                        xmin = 1
                        break
                    else:
                        lasti = i
                        step = step/2
                        i = i -2*step
            else:
                pass
            i+=step
            
        if (xmin==0 ):
            print "No data in file " , file
            sys.exit()
        
        #print "XMIN= ", xmin
        
        #### 2nd loop ###############################  
        
        lasti = 0
        start = nx
        step = 128
        i = start
        while i>=1:
            if (lasti==i):
                std = 1.0
            else:
                print "DEBUG (i,std,lasti):" , i, std, lasti
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
            
        #print "XMAX= ", xmax
        ##### 3rd loop  ####################  
        
        lasti = 0
        start = 1
        step = 128
        i = start
        while i<=ny:
            if (lasti==i):
                std = 1.0
            else:
                std = float(iraf.imstat (
                    images=file+"[*,"+str(i)+"]",
                    fields='stddev',format='no',Stdout=1)[0])
                print "debud --> I, STD, LASTI", i, std, lasti
            
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
        
        #print "YMIN= ", ymin
        ##### 4th loop  ####################  
        
        lasti = 0
        start = ny
        step = 128
        i = start
        while i>=1:
            if (lasti==i):
                std=1.0
            else:
                print "debud --> I, STD, LASTI", i, std, lasti
                std = float(iraf.imstat (
                    images=file+"[*,"+str(i)+"]",
                    fields='stddev',format='no',Stdout=1)[0])
            
            if (std!=0.0):
                if (i==ny):
                    ymax = ny
                    break
                else:
                    if (step==1):
                        ymax = i
                        break
                    else:
                        lasti = i
                        step = step/2
                        i=i+2*step
            else:
                pass
            i-=step
                
        #print "YMAX= ", ymax          
        
        
        
        xmin+=10
        ymin+=15
        xmax-=10
        ymax-=15
        
        ####### Finally ##########
        iraf.imcopy(input=file+"["+str(xmin)+":"+str(xmax)+","+str(ymin)+":"+str(ymax)+"]", 
                    output=file)
        
        
        ima_sec = file.replace(".fits", ".weight.fits")
        if os.path.exists(ima_sec):
            iraf.imcopy(input=ima_sec+"["+str(xmin)+":"+str(xmax)+","+str(ymin)+":"+str(ymax)+"]", output=ima_sec)
        
        ima_objs = file.replace(".fits", "objs.fits")
        if os.path.exists( ima_objs ):
            iraf.imcopy(input=ima_objs+"["+str(xmin)+":"+str(xmax)+","+str(ymin)+":"+str(ymax)+"]", 
                    output=ima_objs)
                
        
        log.debug("....End of  imgTrim ....")
        #### End of the function #######
        
################################################################################
# main
if __name__ == "__main__":
    print 'Start imtrim....'
    
    # Read command line parameters
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'i:v',['input=','verbose'])
    except getopt.GetoptError:
        usage()
        sys.exit(1)
    
    nargs = len(sys.argv[1:])
    nopts = len(opts)
      
    
    verbose = False
    inputfile=''
            
    for option, par in opts:
        if option in ('-i','--input'):
            inputfile = par
        elif option in ('-v','--verbose'):      # verbose debugging output
            verbose = True
            
    if inputfile=='':
        usage()
        sys.exit(2)
    else:
        imgTrim( inputfile )
