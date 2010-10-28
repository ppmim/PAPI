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
      Represents a Classified FITS of PANIC
      A class for PANIC FITS file classification and recognition. Actually it is a
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
         SKY_FOR,              # Sky for extended objects during dither sequence
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
        self.ncoadds   = -1
        self.itime     = 0.0
        self.readmode  = ""
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
    
    def isSky(self):
        return (self.type=="SKY_FOR")
    
    def isObject(self):
        return (self.isScience())
    
    def isMEF(self):
        return (self.mef)
    
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
    
    def getMJD(self):
        return self.mjd
    
    def getData(self,ext=0):
        """No tested, to check !!!"""
        if ext!=0 and not f.mef:
            log.error("Wrong extension number especified. Not a MEF file.")
        else:
            try:
                myfits = pyfits.open(self.pathname)
                temp=myfits[0].data
                myfits.close()
                return temp
            except:
                log.error('Could not open frame - something wrong with input data')
                raise
    

    def recognize(self):
     
        log.info("Recognizing file %s" %self.filename)

        # Check the file exists
        if not os.path.exists( self.pathname ):
            log.error('Cannot find frame : "%s"' % self.pathname)
            raise Exception("File %s not found", self.pathname)
                    
        # Open the file            
        try:
            myfits = pyfits.open(self.pathname)                                 
            #myfits = pyfits.open(self.pathname, 'update')
            #myfits[0].verify()
        except:
            log.error('Could not open frame - something wrong with input data')
            raise
        
        #Check if is a MEF file 
        if len(myfits)>1:
            self.mef=True
            self.next=len(myfits)-1
            log.debug("Found a MEF file with %d extensions", self.next)
        else:
            self.mef=False
            self.next=1
            log.debug("Found a simple FITS file")
        
        # If file is a MEF, some values will be read from the header extensions      
        if self.mef:
            # we suppose all extension have the same dimensions
            self.naxis1=myfits[1].header['NAXIS1']
            self.naxis2=myfits[1].header['NAXIS2']
        else:
            self.naxis1=myfits[0].header['NAXIS1']
            self.naxis2=myfits[0].header['NAXIS2']
        # pointer to the primary-main header
        self.my_header = myfits[0].header
          
        # First, find out the type of frame ( DARK, DOME_FLAT_LAMP_ON/OFF, SKY_FLAT, SCIENCE , UNKNOW)  
        try:
            if myfits[0].header['OBJECT'].lower().count('master'):
                self.type=myfits[0].header['PAPI.TYPE']
            elif myfits[0].header['OBJECT'].lower().count('dark') :
                self.type="DARK"
            elif myfits[0].header['OBJECT'].lower().count('lamp off'):
                self.type="DOME_FLAT_LAMP_OFF"
            elif myfits[0].header['OBJECT'].lower().count('lamp on'):
                self.type="DOME_FLAT_LAMP_ON"
            elif myfits[0].header['OBJECT'].lower().count('dusk'):
                self.type="TW_FLAT_DUSK"
            elif myfits[0].header['OBJECT'].lower().count('dawn'):
                self.type="TW_FLAT_DAWN"
            elif myfits[0].header['OBJECT'].lower().count('skyflat') or myfits[0].header['OBJECT'].lower().count('flat'):
                self.type="TW_FLAT"
            elif myfits[0].header['OBJECT'].lower().count('sky for'):
                self.type="SKY_FOR"
            elif myfits[0].header['OBJECT'].lower().count('focus'):
                self.type="SCIENCE"  #por una razon que desconozco, CAHA le asigna el id 'focus' !!!
            else:
                self.type="SCIENCE"
            log.debug("Image type: %s", self.type)
        except KeyError:
            log.warning('OBJECT keyword not found')
            self.type='UNKNOW'
        
        #Is pre-reduced the image ? by default, no
        self.processed=False
        

        #print "File :"+ self.pathname
        #Filter
        try:
            self.filter  = myfits[0].header['FILTER']
        except KeyError:
            log.warning('FILTER keyword not found')
            self.filter  = 'UNKNOWN'
        
        #Exposition Time
        try:
            self.exptime=myfits[0].header['EXPTIME']
        except KeyError:
            log.warning('EXPTIME keyword not found')
            self.exptime  = -1

        #Integration Time
        try:
            self.itime=myfits[0].header['ITIME']
        except KeyError:
            log.warning('ITIME keyword not found')
            self.itime  = -1
            
        #Number of coadds
        try:
            self.ncoadds=myfits[0].header['NCOADDS']
        except KeyError:
            log.warning('NCOADDS keyword not found')
            self.ncoadda  = -1
                 
        #Read-Mode
        try:
            self.readmode=myfits[0].header['READMODE']
        except KeyError:
            log.warning('READMODE keyword not found')
            self.readmode  = ""
                     
        #UT-date of observation
        try:
            self.datetime_obs = myfits[0].header['DATE-OBS']
            if self.datetime_obs.count('T'):
                self.date_obs, self.time_obs = self.datetime_obs.split('T')
            else:
                self.date_obs = self.datetime_obs
                self.time_obs = ''
        except KeyError:
            log.warning('DATE-OBS keyword not found')
            self.date_obs  = ''
            self.time_obs  = ''
            self.datetime_obs = ''
        
        #RA-coordinate (in degrees)
        try:
            self.ra = myfits[0].header['RA']
        except KeyError:
            log.warning('RA keyword not found')
            self.ra  = ''
        
        #Dec-coordinate (in degrees)
        try:
            self.dec = myfits[0].header['DEC']
        except KeyError:
            log.warning('DEC keyword not found')
            self.dec  = ''

        #MJD-Modified julian date 'days' of observation
        try:
            self.mjd = myfits[0].header['MJD-OBS']
        except KeyError:
            log.warning('MJD-OBS keyword not found')
            self.mjd  = ''
           
        #OBJECT
        try:
            self.object = myfits[0].header['OBJECT']
        except KeyError:
            log.warning('OBJECT keyword not found')
            self.object  = ''   
        
        #CHIPCODE
        try:
            self.chipcode = myfits[0].header['CHIPCODE']
        except KeyError:
            log.warning('CHIPCODE keyword not found, setting default (1)')
            self.chipcode  = 1  # default
        
        #DetectorID
        self.detectorID='O2k'
        #RunID
        self.runID='0'
        
        
        # PRESS1 and PRESS2 wrong keyword values
        try:
            if myfits[0].header['INSTRUME']=='Omega2000':
                myfits[0].header.update('PRESS1', 0.0)
                myfits[0].header.update('PRESS2', 0.0)
        except KeyError:
            log.warning('INSTRUME keyword not found')
        
        myfits.close()
			
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
    """
        This function will check all the files in the file_list have the same properties required as True in the parameters
        Note that that properties should be available in the main header(0) in case of MEF files
    """
    
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
            debug=0
            if debug:
                print 'FILE=',file
                print '------------------------'
                print 'TYPE=',f.getType()
                print 'FILTER=',f.getFilter()
                print 'TEXP=',f.expTime()
                print 'NCOADDS=',f.getNcoadds()
                print 'READMODE=',f.getReadMode()
            
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

        
