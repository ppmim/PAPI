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
  Interface and basic operations with the Display tool (DS9
"""



# Import necessary modules

# System modules
import os
import config



# Other modules with which communication is needed (the logging facility
# and the helper to autoload configuration files)
import messageLog
import autoLoader

# Utilities
import frameUtils
import utils
import datetime

next_frameno=1

################################################################################
# Launch the DS9 display

def startDisplay(local):

  messageLog.put('Starting DS9 display...')
  # First all, check if already is running
  pid = os.fork()
  if pid ==0:
    #First all, check if already is running
    stdout_handle=os.popen("/sbin/pidofproc ds9","r")
    if stdout_handle.read() =='':
      messageLog.put_warning("Display is not running, so we start it")
      # DS9 is not running, so we start it  
      os.execl(("%s/ds9" % local), "ds9","-cmap","Heat")
      # Now, we create 4 new frames and enable #1  
      os.system(("%s/xpaset -p ds9 frame new" % local)) #2
      os.system(("%s/xpaset -p ds9 frame new " % local)) #3
      os.system(("%s/xpaset -p ds9 frame new" % local)) #4
      os.system(("%s/xpaset -p ds9 frame frameno 1" % (local)))
    else:
      # DS9 is already running...
      messageLog.put_warning("Display ALREADY is running")
      #os.system(("%s/xpaset -p ds9 frame delete all" % local))
    os.wait()[0]
    
################################################################################
# Show the current frame into DS9 display

def showFrame(local, frame):
  
  global next_frameno
  messageLog.put('Showing frame %s nro.: %d'  %( frame, next_frameno))
  
  #os.system(("%s/xpaset -p ds9 frame new" % local))
  os.system(("%s/xpaset -p ds9 frame frameno %d" % (local, next_frameno)))  
  os.system(("%s/xpaset -p ds9 frame reset" % local))
  os.system(("%s/xpaset -p ds9 tile" % local))
  os.system(("%s/xpaset -p ds9 cmap Heat" % local))
  os.system(("%s/xpaset -p ds9 file %s" %(local, frame)))
  os.system(("%s/xpaset -p ds9 scale zscale" % local ))
  os.system(("%s/xpaset -p ds9 zoom to fit" % local))

  if next_frameno == 4:
    next_frameno=1
  else:
    next_frameno=next_frameno+1

################################################################################
# Show the current frame into DS9 display

def showFrames(local, framelist):
  
  global next_frameno

  print "FRAMELIST TO SHOW ", framelist
  nframes=len(framelist)
  for i in range(nframes):
      #os.system(("%s/xpaset -p ds9 frame new" % local))
      os.system(("%s/xpaset -p ds9 frame frameno %d" % (local, next_frameno)))  
      os.system(("%s/xpaset -p ds9 frame reset" % local))
      os.system(("%s/xpaset -p ds9 tile" % local))
      os.system(("%s/xpaset -p ds9 cmap Heat" % local))
      os.system(("%s/xpaset -p ds9 file %s" %(local, framelist[i])))
      os.system(("%s/xpaset -p ds9 scale zscale" % local ))
      os.system(("%s/xpaset -p ds9 zoom to fit" % local))
      messageLog.put('Showing frame %s nro.: %d'  %( framelist[i], next_frameno))
      
      if next_frameno == 4:
          next_frameno=1
      else:
          next_frameno=next_frameno+1
          
    
################################################################################
#
def stopDisplay():
    
    os.system("/usr/local/bin/xpaset -p ds9 exit")
