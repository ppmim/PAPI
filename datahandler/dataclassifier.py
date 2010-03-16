################################################################################
#
# PANICtool
#
# DataClassifier.py
#
# Created    : 07/04/2008    jmiguel@iaa.es
# Update     : 25/05/2009    jmiguel@iaa.es   Added object field
#              14/12/2009    jmiguel@iaa.es   Renamed SKY_FLAT by TW_FLAT
#                                             Reanamed isSkyFlat()  by isTwFlat()
#              02/03/2010    jmiguel@iaa.es   Added READMODE checking
# 
################################################################################

"""
   Class for FITS data classifying
"""

# import external modules


import os.path
import pyfits
import sys

from time import strftime

#Log
from misc.paLog import log


################################################################################

################################################################################

class ClFits:

    """
      Classified FITS
      A class for PANIC FITS file classification and recognition. Really it is a
      wrapper for pyfits.

      TYPE of files to be considered :
         UNKNOW,
         DARK,
         TW_FLAT,
         TW_FLAT_DUSK,
         TW_FLAT_DAWN,
         DOME_FLAT_LAMP_OFF,
         DOME_FLAT_LAMP_ON,
         FOCUS,
         SCIENCE, (RAW)
         MASTER_DARK,
         MASTER_SKY_FLAT,
         MASTER_DOME_FLAT,
         MASTER_TW_FLAT,
         SCIENCE_REDUCED
         
    """

    # Class initialization
    def __init__(self, full_pathname ):
      
        """
        \param pathname The full filename of the fits file to classify
        """
      
        self.pathname  = full_pathname # path_name + root_name
        self.filename  = os.path.basename(self.pathname)  # root_name
        self.type      = "UNKNOW"
        self.naxis1    = -1
        self.naxis2    = -1
        self.runID     = "-1"
        self.dataSetID = -1
        self.processed = False
        self.exptime   = -1
        self.filter    = ""
        self.mef       = False
        self.next      = 0
        self.datetime_obs = ""
        self.detectorID= ""
        self.date_obs  = ""
        self.time_obs  = ""
        self.ra        = ""
        self.dec       = ""
        self.mdj       = "" # modified julian date of observation
        self.object    = ""
        self.chipcode  = 1
        self.ncoadds   = 0
        self.itime     = 0.0
        self.readmode  = ""
        self.data      = None
        self.my_header = None
        
        
        self.recognize()
		
    def getType(self):
        return self.type
    
    def getFilter(self):
        return self.filter
    
    def getNaxis1(self):
        return self.naxis1
    
    def getNaxis2(self):
        return self.naxis2
    
    def pathname(self):
        return self.pathname
    
    def isDark(self):
		return (self.type=="DARK")

    def isDomeFlat(self):
        return (self.type=="DOME_FLAT_LAMP_ON" or self.type=="DOME_FLAT_LAMP_OFF" or self.type=="DOME_FLAT")
    
    def isDomeFlatON(self):
        return (self.type=="DOME_FLAT_LAMP_ON")
    
    def isDomeFlatOFF(self):
        return (self.type=="DOME_FLAT_LAMP_OFF")
    
    def isTwFlat(self):
        return (self.type=="TW_FLAT_DUSK" or self.type=="TW_FLAT_DAWN" or self.type=="TW_FLAT")
	
    def isScience(self):
        return (self.type.count("SCIENCE"))
    
	def isRawScience(self):
		return (self.type=="SCIENCE_RAW")
    
    def isReducedScience(self):
        return (self.type=="SCIENCE_REDUCED")
    
    def isMasterDark(self):
        return (self.type=="MASTER_DARK")
    
    def expTime(self):
        return (self.exptime)
    
    def getRA(self):
        return self.ra
    
    def getDec(self):
        return self.dec
    
    def getNcoadds(self):
        return self.ncoadds
    
    def getItime(self):
        return self.itime
    
    def getReadMode(self):
        return self.readmode
    
    def getDateTimeObs(self):
        return self.datetime_obs
    
    def getData(self):
        """No tested, to check !!!"""
        try:
            indata = pyfits.open(self.pathname)
            temp=indata[0].data
            indata.close()
            return temp
        except:
            log.error('Could not open frame - something wrong with input data')
            raise
    

    def recognize(self):
     
        log.debug("Recognizing file %s" %self.filename)

        # Check the file exists
        if not os.path.exists( self.pathname ):
            log.error('Cannot find frame : "%s"' % self.pathname)
            raise "File not found"
        
        # Verify fits file
        #if not self.pathname.endswith(".fits") or not self.pathname.endswith("*.fit"):
        #    log.('File %s does not seem FITS file'  %self.pathname)
        #    raise NameError, 'NotAFitsFile'
        
        try:
            indata = pyfits.open(self.pathname)                                 
            #indata = pyfits.open(self.pathname, 'update')
            #indata[0].verify()
        except:
            log.error('Could not open frame - something wrong with input data')
            raise
        
        self.naxis1=indata[0].header['NAXIS1']
        self.naxis2=indata[0].header['NAXIS2']
        self.my_header = indata[0].header
        
        # First, find out the type of frame ( DARK, DOME_FLAT_LAMP_ON/OFF, SKY_FLAT, SCIENCE , UNKNOW)  
        try:
            if indata[0].header['OBJECT'].count('dark'):
                self.type="DARK"
            elif indata[0].header['OBJECT'].count('lamp off'):
                self.type="DOME_FLAT_LAMP_OFF"
            elif indata[0].header['OBJECT'].count('lamp on'):
                self.type="DOME_FLAT_LAMP_ON"
            elif indata[0].header['OBJECT'].count('dusk'):
                self.type="TW_FLAT_DUSK"
            elif indata[0].header['OBJECT'].count('dawn'):
                self.type="TW_FLAT_DAWN"
            elif indata[0].header['OBJECT'].count('skyflat'):
                self.type="TW_FLAT"
            elif indata[0].header['OBJECT'].count('MASTER'):
                self.type=indata[0].header['PAPI.TYPE']
            elif indata[0].header['OBJECT'].count('focus'):
                self.type="SCIENCE"  #por una razon que desconozco, CAHA le asigna el id 'focus' !!!
            else:
                self.type="SCIENCE"
        except KeyError:
            log.error('Error, keyword not found')
            self.type='UNKNOW'
        
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

        #print "File :"+ self.pathname
        #Filter
        try:
            self.filter  = indata[0].header['FILTER']
        except KeyError:
            log.error('Error, FILTER keyword not found')
            self.filter  = 'UNKNOWN'
        
        #Exposition Time
        try:
            self.exptime=indata[0].header['EXPTIME']
        except KeyError:
            log.error('Error, EXPTIME keyword not found')
            self.exptime  = -1

        #Integration Time
        try:
            self.itime=indata[0].header['ITIME']
        except KeyError:
            log.error('Error,  ITIME keyword not found')
            self.itime  = -1
            
        #Number of coadds
        try:
            self.ncoadds=indata[0].header['NCOADDS']
        except KeyError:
            log.error('Error, NCOADDS keyword not found')
            self.ncoadda  = -1
                 
        #Read-Mode
        try:
            self.readmode=indata[0].header['READMODE']
        except KeyError:
            log.error('Error, READMODE keyword not found')
            self.readmode  = ""
                     
        #UT-date of observation
        try:
            self.datetime_obs = indata[0].header['DATE-OBS']
            if self.datetime_obs.count('T'):
                self.date_obs, self.time_obs = self.datetime_obs.split('T')
            else:
                self.date_obs = self.datetime_obs
                self.time_obs = ''
        except KeyError:
            log.error('Error, DATE-OBS keyword not found')
            self.date_obs  = ''
            self.time_obs  = ''
            self.datetime_obs = ''
        
        #RA-coordinate (in degrees)
        try:
            self.ra = indata[0].header['RA']
        except KeyError:
            log.error('Error, RA keyword not found')
            self.ra  = ''
        
        #Dec-coordinate (in degrees)
        try:
            self.dec = indata[0].header['DEC']
        except KeyError:
            log.error('Error, DEC keyword not found')
            self.dec  = ''

        #MJD-Modified julian date 'days' of observation
        try:
            self.mjd = indata[0].header['MJD-OBS']
        except KeyError:
            log.error('Error, MJD-OBS keyword not found')
            self.mjd  = ''
           
        #OBJECT
        try:
            self.object = indata[0].header['OBJECT']
        except KeyError:
            log.error('Error, OBJECT keyword not found')
            self.object  = ''   
        
        #CHIPCODE
        try:
            self.chipcode = indata[0].header['CHIPCODE']
        except KeyError:
            log.warning('Warning, CHIPCODE keyword not found, setting default (1)')
            self.chipcode  = 1  # default
        
        #DetectorID
        self.detectorID='O2k'
        #RunID
        self.runID='0'
        
        
        # PRESS1 and PRESS2 wrong keyword values
        try:
            if indata[0].header['INSTRUME']=='Omega2000':
                indata[0].header.update('PRESS1', 0.0)
                indata[0].header.update('PRESS2', 0.0)
        except KeyError:
            log.error('Warning,no INSTRUME keyword not found')
        
        indata.close()
			
    def addHistory(self, string_history):
        """ Add a history keyword to the header
            NOTE: The update is only done in memory-header(my_header). To flush to disk, we should
            to write/update to a new file. 
            
            STILL NOT USED !!!!
        """
        log.warning("Header not updated nicely !!")
        try:
            #t=pyfits.open(self.pathname,'update')
            self.my_header.add_history(string_history)
            #t.close(output_verify='ignore')
            log.warning('History added')
        except:
            log.error('Error while adding history to header')
            
        
        
    def printClass(self):
        print "---------------------------------"
        print "Fichero   : ", self.pathname
        print "Data class: ", self.type
        print "Filter    : ", self.filter
        print "Processed : ", self.processed 
        print "TEXP      : ", self.exptime
        print "MEF       : ", self.mef
        print "Date-Obs  : ", self.date_obs
        print "Time-Obs  : ", self.time_obs
        print "RA        : ", self.ra
        print "Dec       : ", self.dec
        print "MJD       : ", self.mjd
        print "---------------------------------"
		
        
      
################################################################################            
#  Useful function to check data integrity
################################################################################
def checkDataProperties( file_list, c_type=True, c_filter=True, c_texp=True, c_ncoadds=True, c_readmode=True):
    """This function will check all the files in the file_list have the same properties required as True in the parameters"""
    
    m_type=''
    m_filter=''
    m_texp=''
    m_ncoadds=''
    m_readmode=''
    
    # First file as reference
    if file_list[0]:
            f=ClFits ( file_list[0] )
            m_type=f.getType()
            m_filter=f.getFilter()
            m_texp=f.expTime()
            m_ncoadds=f.getNcoadds()
            m_readmode=f.getReadMode()
    
    # Check all files        
    for file in file_list:
            f=ClFits ( file )
            if (  (c_type and m_type!=f.getType()) or (c_filter and m_filter!=f.getFilter()) or (c_texp and m_texp!=f.expTime()) or (c_ncoadds and m_ncoadds!=f.getNcoadds()) or (c_readmode and m_readmode!=f.getReadMode())):
                log.debug("Missmath on some property")
                return False
        
    log.debug("Successful properties checking")
    return True      
      
		
################################################################################            
#  Main for Testing
################################################################################
def usage():
    """Print help """
    print "Usage: %s filename", sys.argv[0]
    
if __name__ == "__main__": 
    log.debug("Starting the tests.....")
    
    if len(sys.argv)>1:
        try:
            dr = ClFits(sys.argv[1]) #"/disk-a/caha/panic/DATA/data_mat/skyflat0020.fits")
            dr.printClass()
            print dr.isTwFlat()
        except:
            log.error("Error reading fits File")
    else:
        usage()

        
