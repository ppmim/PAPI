################################################################################
#
# PANICtool
#
# display.py
#
# Last update 17/04/2008
#
################################################################################
"""
  Interface and basic operations with the Display tool (DS9)
"""



# Import necessary modules

# System modules
import os

# Other modules with which communication is needed (the logging facility
# and the helper to autoload configuration files)
#import messageLog

# Utilities
import datetime
import sys
import time

#Log
from misc.paLog import log
import datahandler

ds9_path="/usr/local/bin/"

################################################################################
# Launch the DS9 display
def startDisplay():
  """ funtion to launch the DS9 display, checking previusly if it is already started"""

  #First all, check if already is running
  stdout_handle=os.popen("/sbin/pidofproc ds9","r")
  if stdout_handle.read() =='':
    print "DS9 not running "
    # DS9 is not running, so we start it  
    os.system(("%s/ds9 &" % ds9_path))
    time.sleep(2)
    stdout_handle=os.popen("/sbin/pidofproc ds9","r")
    if stdout_handle.read() =='':
        time.sleep(1)
    time.sleep(1)
  else:
    # DS9 is already running...
    print "Warning, display is already running"
    #os.system(("%s/xpaset -p ds9 frame delete all" % ds9_path))
      
def startDisplay2():
   
  """ NOT USED, because Esta funcion da muchos problemas cuando se llama desde el QL, por eso no la usamos"""

  #messageLog.put('Starting DS9 display...')
  # First all, check if already is running
  pid = os.fork()
  if pid ==0:
    ## The child 
    #First all, check if already is running
    stdout_handle=os.popen("/sbin/pidofproc ds9","r")
    if stdout_handle.read() =='':
      print "DS9 not running "
      #messageLog.put_warning("Display is not running, so we start it")
      # DS9 is not running, so we start it  
      #os.execl(("%sds9" % ds9_path), "ds9","-cmap","Heat")  ## ojo, reemplaza al proceso actual, hace un fork
      os.system(("%s/ds9 &" % ds9_path))
      time.sleep(1)
      #os.execl(("/usr/local/bin/ds9"), "ds9", "-cmap", "Heat")
      # Now, we create 4 new frames and enable #1  
      os.system(("%s/xpaset -p ds9 frame new" % ds9_path)) #2
      os.system(("%s/xpaset -p ds9 frame new " % ds9_path)) #3
      os.system(("%s/xpaset -p ds9 frame new" % ds9_path)) #4
      os.system(("%s/xpaset -p ds9 frame frameno 1" % (ds9_path)))
    else:
      # DS9 is already running...
      #messageLog.put_warning("Display ALREADY is running")
      print "Warning, display is already running"
      #os.system(("%s/xpaset -p ds9 frame delete all" % ds9_path))
    sys.exit(0)
  else:
    # The parent 
    os.wait()[0]
    
################################################################################
# Show the current frame into DS9 display

def showFrame(frame):

  f=datahandler.ClFits(frame)

  #Check display
  startDisplay()
  
  if (f.mef==True):
        # Multi-Extension FITS files
        if (f.isDark()):
            # (Hawki) Dark files don't have WCS information required by ds9 mosaicimage
            os.system(("%s/xpaset -p ds9 frame delete all" % (ds9_path)))
            os.system(("%s/xpaset -p ds9 frame new" % ds9_path))
            os.system(("%s/xpaset -p ds9 cmap Heat" % ds9_path))
            os.system(("%s/xpaset -p ds9 scale zscale" % ds9_path ))
            os.system(("%s/xpaset -p ds9 file multiframe %s" % (ds9_path, frame)))
        else:
            # Beware, 'mosaicimage' ds9 facility require WCS information 
            os.system(("%s/xpaset -p ds9 file mosaicimage %s" % (ds9_path, frame)))
            os.system(("%s/xpaset -p ds9 cmap Heat" % ds9_path))
            os.system(("%s/xpaset -p ds9 scale zscale" % ds9_path ))
            os.system(("%s/xpaset -p ds9 zoom to fit" % ds9_path))
  else:
        # Single FITS files
        os.system(("%s/xpaset -p ds9 frame delete all" % (ds9_path)))
        os.system(("%s/xpaset -p ds9 frame new" % ds9_path))
        os.system(("%s/xpaset -p ds9 single" % ds9_path))
        os.system(("%s/xpaset -p ds9 cmap Heat" % ds9_path))
        os.system(("%s/xpaset -p ds9 file %s" %(ds9_path, frame)))
        os.system(("%s/xpaset -p ds9 scale zscale" % ds9_path ))
        os.system(("%s/xpaset -p ds9 zoom to fit" % ds9_path))


################################################################################
# Show the current frame into DS9 display

def showSingleFrames(framelist):
  """ Display a single frame, not supposing a MEF file """
  
  global next_frameno

  #Check display
  startDisplay()
  
  nframes=len(framelist)
  for i in range(nframes):
      #os.system(("%s/xpaset -p ds9 frame new" % ds9_path))
      os.system(("%s/xpaset -p ds9 frame frameno %d" % (ds9_path, next_frameno)))
      os.system(("%s/xpaset -p ds9 frame reset" % ds9_path))
      os.system(("%s/xpaset -p ds9 tile" % ds9_path))
      os.system(("%s/xpaset -p ds9 cmap Heat" % ds9_path))
      os.system(("%s/xpaset -p ds9 file %s" %(ds9_path, framelist[i])))
      os.system(("%s/xpaset -p ds9 scale zscale" % ds9_path ))
      os.system(("%s/xpaset -p ds9 zoom to fit" % ds9_path))
      
      if next_frameno == 4:
          next_frameno=1
      else:
          next_frameno=next_frameno+1
          
    
################################################################################
#
def stopDisplay():
    """Kill the actual ds9 display, if already exists"""
    
    os.system("/usr/local/bin/xpaset -p ds9 exit")

################################################################################
#
if __name__ == "__main__":
   startDisplay()
   
   
