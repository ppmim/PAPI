################################################################################
#
# PANICtool
#
# fitsClass.py
#
# Last update 07/04/2008
#
################################################################################

"""
   FITS class for data structuring
"""

# import external modules

#import config
#import messageLog
import os.path
import pyfits


from time import strftime

detected_files=[]
################################################################################

# Define local constants

class FitsFile:

    """
      A class for PANIC FITS file classification and recognition. Really it is a
      wrapper for pyfits.

      TYPE of files to be considered :
         UNKNOW,
         DARK,
         FLAT_SKY,
         FLAT_SKY_DUSK,
         FLAT_SKY_DAWN,
         FLAT_DOME_ON,
         FLAT_DOME_OFF,
         FOCUS,
         SCIENCE_RAW,
         MASTER_DARK,
         MASTER_FLAT_SKY,
         MASTER_FLAT_DOME,
         SCIENCE_RED
         
    """
    

    # Class initialization
    def __init__(self, pathname ):

        self.pathname  = pathname
        self.filename  = os.path.basename(pathname)
        self.type      = "UNKNOW"
        self.size      = -1
        self.runID     = "-1"
        self.dataSetID = -1
        self.processed = False
        self.exptime   = -1
        self.filter    = ""
        self.mef       = False
        self.next      = 0
        self.date      = ""
        self.detectorID= ""
        self.date_obs  = ""
        self.time_obs  = ""
        
        self.recognize()

    def recognize(self):

        # Check the file exists
        if not os.path.exists( self.pathname ):
            print 'Cannot find frame : "%s"' % self.filename
            return
        if not self.pathname.endswith(".fits"):
            print 'File %s does not seem FITS file'  %self.pathname
            raise NameError, 'NotAFitsFile'
        
        try:
            indata = pyfits.open(self.pathname)
        except:
            print 'Could not open frame - something wrong with input data'
            raise
        
        # First, find out the type of frame ( DARK, DOME_FLAT_LAMP_ON/OFF, SKY_FLAT, SCIENCE , UNKNOW)  
        if indata[0].header['OBJECT'].count('dark'):
            self.type="DARK"
        elif indata[0].header['OBJECT'].count('lamp off'):
            self.type="DOME_FLAT_LAMP_OFF"
        elif indata[0].header['OBJECT'].count('lamp on'):
            self.type="DOME_FLAT_LAMP_ON"
        elif indata[0].header['OBJECT'].count('dusk'):
            self.type="SKY_FLAT_DUSK"
        elif indata[0].header['OBJECT'].count('dawn'):
            self.type="SKY_FLAT_DAWN"
        else:
            self.type="SCIENCE"

        #Is pre-reduced the image ? by default, no
        self.processed=False

        #EXTEND
        try:    
            if indata[0].header['EXTEND']==True:
                self.mef=True
                next=indata[0].header['NEXTEND']
            else:
                self.mef=False
        except KeyError:
            #print 'Error, EXTEND keyword not found'
            self.mef=False
            pass

        #print "File :"+ self.pathname
        #Filter
        try:
            self.filter  = indata[0].header['FILTER']
        except KeyError:
            print 'Error, FILTER keyword not found'
            self.filter  = 'UNKNOWN'
            pass
        
        #Exposition Time
        try:
            self.exptime=indata[0].header['EXPTIME']
        except KeyError:
            print 'Error, EXPTIME keyword not found'
            self.exptime  = -1
            pass   

        #UT-date of observation
        try:
            self.datetime_obs = indata[0].header['DATE-OBS']
            self.date_obs, self.time_obs = self.datetime_obs.split('T')
        except KeyError:
            print 'Error, DATE-OBS keyword not found'
            self.date_obs  = ''
            self.time_obs  = ''
            self.datetime_obs = ''
            pass
        
        #DetectorID
        self.detectorID='O2k'
        #RunID
        self.runID='0'
               
        indata.close()
            
            

        
