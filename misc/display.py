################################################################################
#
# PANICtool
#
# display.py
#
# Last update 09/10/2019
#
################################################################################
"""
  Interface and basic operations with the Display tool (DS9)
"""

# System modules
import os
from distutils import spawn
import time

# papi
from misc.paLog import log
import datahandler


try:
    ds9_path = os.path.dirname(spawn.find_executable("ds9"))
except Exception:
    msg = "Cannot find DS9 binary. Check your ~/.papirc config file"
    log.error(msg)
    raise Exception(msg)

try:
    xpa_path = os.path.dirname(spawn.find_executable("xpaaccess"))
except Exception:
    msg = "Cannot find XPA bin directory. Check your ~/.papirc config file"
    log.error(msg)
    raise Exception(msg)

frame_no = 0
created_frames = 0
MAX_FRAMES_NO = 10

################################################################################
# Launch the DS9 display
def startDisplay():
    """
    funtion to launch the DS9 display, checking previusly if it is already 
    started
    """

    # First all, check if already is running
    # stdout_handle = os.popen("/sbin/pidofproc ds9","r")
    stdout_handle = os.popen("%s/xpaaccess ds9"%xpa_path, "r")
    if stdout_handle.read() == 'no\n':
        # print "DS9 not running "
        # DS9 is not running, so we start it  
        os.system(("%s/ds9 &" % ds9_path))
        time.sleep(2)
        stdout_handle = os.popen("%s/xpaaccess ds9"%xpa_path, "r")
        if stdout_handle.read() == 'no\n':
            time.sleep(3)
        time.sleep(1)
        #for i in range(0, MAX_FRAMES_NO-1): # when ds9 start, it has already one frame
        #    os.system(("%s/xpaset -p ds9 frame new" % xpa_path))
            
    else:
        pass
        # DS9 is already running...
        #print "Warning, display is already running"
        #os.system(("%s/xpaset -p ds9 frame delete all" % xpa_path))
      
    
################################################################################
# Show the current frame into DS9 display

def showFrame(frame, del_all=False):
  """
  Show in DS9 display the file/s given in the input
  """
  
  # Check display
  startDisplay()
  
  global frame_no
  global created_frames
  
  # frame could be a single file or a file list 
  if type(frame)==type(list()): 
        fileList = frame
        os.system(("%s/xpaset -p ds9 frame delete all" % (xpa_path)))
        delete_all = False
  elif os.path.isfile(frame):
        delete_all = del_all
        fileList = [frame]
  
  for file in fileList:
        f = datahandler.ClFits(file)
        if f.mef:
                # Multi-Extension FITS files
                if f.getInstrument() == 'hawki':
                    # HAWK-I Dark files don't have WCS information required by ds9 mosaicimage
                    # PANIC  Dark files do have WCS information.  
                    if delete_all: os.system(("%s/xpaset -p ds9 frame delete all" % (xpa_path)))
                    os.system(("%s/xpaset -p ds9 frame new" % xpa_path))
                    os.system(("%s/xpaset -p ds9 cmap Heat" % xpa_path))
                    os.system(("%s/xpaset -p ds9 scale zscale" % xpa_path ))
                    os.system(("%s/xpaset -p ds9 file multiframe %s" % (xpa_path, file)))
                    #os.system(("%s/xpaset -p ds9 medatacube multiframe %s" % (xpa_path, file)))
                else:
                    # Beware, 'mosaicimage' ds9 facility require WCS information
                    if delete_all: os.system(("%s/xpaset -p ds9 frame delete all" % (xpa_path)))
                    if frame_no<MAX_FRAMES_NO:
                        if created_frames<(frame_no+1):
                            os.system(("%s/xpaset -p ds9 frame new" % xpa_path))
                            created_frames+=1
                        os.system(("%s/xpaset -p ds9 frame frameno %d" % (xpa_path, frame_no+1)))
                        frame_no+=1
                    else:
                        frame_no = 1
                        os.system(("%s/xpaset -p ds9 frame frameno %d" % (xpa_path, frame_no)))
                    os.system(("%s/xpaset -p ds9 single" % xpa_path)) 
                    os.system(("%s/xpaset -p ds9 file mosaicimage %s" % (xpa_path, file)))
                    os.system(("%s/xpaset -p ds9 cmap Heat" % xpa_path))
                    os.system(("%s/xpaset -p ds9 scale zscale" % xpa_path ))
                    os.system(("%s/xpaset -p ds9 zoom to fit" % xpa_path))
        else:
                # Single FITS files
                if delete_all: os.system(("%s/xpaset -p ds9 frame delete all" % (xpa_path)))
                if frame_no<MAX_FRAMES_NO:
                    #os.system(("%s/xpaset -p ds9 frame new" % xpa_path))
                    os.system(("%s/xpaset -p ds9 frame frameno %d" % (xpa_path, frame_no+1)))
                    frame_no+=1
                else:
                    frame_no=1
                    os.system(("%s/xpaset -p ds9 frame frameno %d" % (xpa_path, frame_no))) 
                os.system(("%s/xpaset -p ds9 single" % xpa_path))
                os.system(("%s/xpaset -p ds9 cmap Heat" % xpa_path))
                os.system(("%s/xpaset -p ds9 file %s" %(xpa_path, file)))
                os.system(("%s/xpaset -p ds9 scale zscale" % xpa_path ))
                os.system(("%s/xpaset -p ds9 zoom to fit" % xpa_path))


################################################################################
# Show the current frame into DS9 display

def showSingleFrames(framelist):
  """ Display a single frame, not supposing a MEF file """
  
  global next_frameno

  #Check display
  startDisplay()
  
  nframes=len(framelist)
  for i in range(nframes):
      #os.system(("%s/xpaset -p ds9 frame new" % xpa_path))
      os.system(("%s/xpaset -p ds9 frame frameno %d" % (xpa_path, next_frameno)))
      os.system(("%s/xpaset -p ds9 frame reset" % xpa_path))
      os.system(("%s/xpaset -p ds9 tile" % xpa_path))
      os.system(("%s/xpaset -p ds9 cmap Heat" % xpa_path))
      os.system(("%s/xpaset -p ds9 file %s" %(xpa_path, framelist[i])))
      os.system(("%s/xpaset -p ds9 scale zscale" % xpa_path ))
      os.system(("%s/xpaset -p ds9 zoom to fit" % xpa_path))
      
      if next_frameno == 4:
          next_frameno=1
      else:
          next_frameno=next_frameno+1
          
    
################################################################################
#
def stopDisplay():
    """Kill the actual ds9 display, if already exists"""
    
    os.system(("%s/xpaset -p ds9 exit" % xpa_path))

################################################################################


def showPDF(filename):
    """
    Display a PDF file with the default PDF viewer (mupdf)
    """
    
    if filename.endswith(".pdf"):
        # try to find out a pdf client
        try:
            pdf_path = os.path.dirname(spawn.find_executable("mupdf"))
            os.system(pdf_path + "/mupdf %s" % filename)
        except Exception as e:
            msg = "Cannot display PDF file, cannot find PDF client mupdf"
            log.error(msg)
            raise e

    else:
        log.info("Give file %s does not look a PDF file" % filename)

################################################################################
if __name__ == "__main__":
    # test
    startDisplay()
    f = "../data_test/panic.20191008_0001.fits"
    showFrame(f)

