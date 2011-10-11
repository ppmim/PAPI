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

import pywcs 
import numpy

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
         SKY,              # Sky for extended objects during dither sequence
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
        self.runID     = -1 # exposition ID for the current night (unique for that night)
        self.obID      = -1 # Observation Block ID (unique) provided by the OT
        self.obPat     = -1 # Observation (dither) Pattern 
        self.pat_expno = -1 # Exposure number within (dither) pattern
        self.pat_noexp = -1 # Number of exposures within pattern
        self.instrument = ""
        self.processed = False
        self.exptime   = -1
        self.filter    = ""
        self.mef       = False
        self.next      = 0
        self.datetime_obs = ""
        self.detectorID= ""
        self.date_obs  = ""
        self.time_obs  = ""
        self._ra       = -1
        self._dec      = -1
        self.equinox   = -1
        self.mdj       = -1 # modified julian date of observation
        self.object    = ""
        self.chipcode  = 1
        self.ncoadds   = -1
        self.itime     = 0.0
        self.readmode  = ""
        self.my_header = None
        self.obs_tool  = False
        self._softwareVer = ''
        
        
        self.recognize()

    def getType(self):
        return self.type 
    
    def getFilter(self):
        return self.filter
    
    def getOBId(self):
        """return the Observation Block ID"""
        return self.obID
    
    def getNoExp(self):
        return self.pat_noexp
    
    def getExpNo(self):
        return self.pat_expno
    
    def getOBPat(self):
        return self.obPat
    
    def getNaxis1(self):
        return self.naxis1
    
    def getNaxis2(self):
        return self.naxis2
    
    def pathname(self):
        return self.pathname
    
    def isDark(self):
        return (self.type=="DARK")

    def isDomeFlat(self):
        return (self.type=="DOME_FLAT_LAMP_ON" or 
                self.type=="DOME_FLAT_LAMP_OFF" or self.type=="DOME_FLAT")
    
    def isDomeFlatON(self):
        return (self.type=="DOME_FLAT_LAMP_ON")
    
    def isDomeFlatOFF(self):
        return (self.type=="DOME_FLAT_LAMP_OFF")
    
    def isTwFlat(self):
        return (self.type=="TW_FLAT_DUSK" or self.type=="TW_FLAT_DAWN" 
                or self.type=="TW_FLAT" or self.type=="SKY_FLAT")
    
    def isFocusSerie(self):
        return False # not yet implemented
    
    def isScience(self):
        return (self.type.count("SCIENCE") or self.count("STD"))
    
    def isRawScience(self):
        return (self.type=="SCIENCE_RAW")

    def isReducedScience(self):
        return (self.type=="SCIENCE_REDUCED")
    
    def isMasterDark(self):
        return (self.type=="MASTER_DARK")
    
    def isSky(self):
        return (self.type=="SKY")
    
    def isObject(self):
        return (self.isScience())
    
    def isMEF(self):
        return (self.mef)
    
    def isFromOT(self):
        return (self.obs_tool)
    
    def expTime(self):
        return (self.exptime)
    
    @property
    def ra(self):
        return self._ra
    
    @property
    def dec(self):
        return self._dec
    @property
    def softwareVer(self):
        return self._softwareVer
    
    #def getRA(self):
    #    return self.ra
    
    #def getDec(self):
    #    return self.dec
    
    def getEquinox(self):
        return self.equinox
    
    def getNcoadds(self):
        return self.ncoadds
    
    def getItime(self):
        return self.itime
    
    def getInstrument(self):
        return self.instrument
    
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
                myfits = pyfits.open(self.pathname, 
                                     ignore_missing_end=True) # since some problems with O2k files   
                temp=myfits[0].data
                myfits.close(output_verify='ignore')
                return temp
            except:
                log.error('Could not open frame - something wrong with input data')
                raise
    

    def recognize(self):
     
        ###log.info("Recognizing file %s" %self.filename)

        # Check the file exists
        if not os.path.exists( self.pathname ):
            log.error('Cannot find frame : "%s"' % self.pathname)
            raise Exception("File %s not found"%self.pathname)
                    
        # Open the file            
        try:
            myfits = pyfits.open(self.pathname, 
                                 ignore_missing_end=True) # since some problems with O2k files                                
            #myfits = pyfits.open(self.pathname, 'update')
            #myfits[0].verify()
        except Exception,e:
            log.error('Could not open frame %s - something wrong with input data : %s',self.pathname, str(e))
            raise e
        
        #Check if is a MEF file 
        if len(myfits)>1:
            self.mef=True
            self.next=len(myfits)-1
            log.debug("Found a MEF file with %d extensions", self.next)
        else:
            self.mef=False
            self.next=1
            ###log.debug("Found a simple FITS file")
        
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
          
        
         
        # INSTRUMENT
        try:
            if myfits[0].header.has_key('INSTRUME'):
                self.instrument = myfits[0].header['INSTRUME']
            else:
                self.instrument = "Unknown"
        except Exception,e:
            log.warning("INSTRUME keyword not found")
            self.instrument = "Unknown"
        
        # Find out the how data file were"observed"
        if myfits[0].header.has_key('OBS_TOOL') or self.instrument=='HAWKI':
            self.obs_tool=True
        else:
            self.obs_tool=False
        
        #Software Version (GEIRS Version)
        if myfits[0].header.has_key('SOFTWARE'):
            self._softwareVer = myfits[0].header['SOFTWARE']
        
        
        # IMAGE TYPE
            
        # Some temporal to allow work with diff instrument data files
        inst = myfits[0].header['INSTRUME'].lower()
        try:
            if inst=='omega2000' and myfits[0].header.has_key('OBJECT'):
                keyword_with_frame_type = 'OBJECT'
            elif inst=='hawki' and myfits[0].header.has_key('IMAGETYP'):
                keyword_with_frame_type = 'IMAGETYP'
            elif inst=='hawki' and myfits[0].header.has_key('OBJECT'):
                keyword_with_frame_type = 'OBJECT'
            elif inst=='panic': # current ID in GEIRS for PANIC
                if self.obs_tool:
                    keyword_with_frame_type = 'IMAGETYP'
                    #keyword_with_frame_type = 'OBJECT'
                else:
                    keyword_with_frame_type = 'OBJECT'
            else: keyword_with_frame_type = 'OBJECT' # default, even for 'Panic'    
        except KeyError:
            log.warning('INSTRUME keyword not found')
            keyword_with_frame_type = 'OBJECT' #default
            
        # First, find out the type of frame ( DARK, DOME_FLAT_LAMP_ON/OFF, SKY_FLAT, SCIENCE , MASTER_calibration, UNKNOW)     
        panic_types = {'dark':'DARK', 
                       'lamp_on_flat':'DOME_FLAT_LAMP_ON', 
                       'lamp_off_flat':'DOME_FLAT_LAMP_OFF',
                       'tw_flat_dawn':'TW_FLAT_DAWN', 
                       'tw_flat_dusk':'TW_FLAT_DUSK', 
                       'sky_flat':'SKY_FLAT',
                       'sky': 'SKY',
                       'science': 'SCIENCE',
                       'std':'STD',
                       'focus':'FOCUS'
                       }
        #IMAGETYP = [DARK, LAMP_ON_FLAT, LAMP_OFF_FLAT, TW_FLAT_DUSK, TW_FLAT_DAWN, SKY_FLAT, SCIENCE, SKY, STD, FOCUS ]
        #FIELDTYP = [POINTLIKE, SPARSE_FIELD, CROWDED_FIELD, EXT_OBJECT ]
        if inst =='panic':
            try:
                if myfits[0].header.has_key('PAPITYPE'):
                    self.type = myfits[0].header['PAPITYPE']
                else:
                    if myfits[0].header.has_key('IMAGETYP'):
                        ltype = myfits[0].header['IMAGETYP'].lower()
                    else:
                        ltype = myfits[0].header[keyword_with_frame_type].lower()
                    try:
                        self.type = panic_types[ltype]
                    except KeyError:
                        log.error("Frame type '%s' does not match any kind !"%ltype)
            except KeyError:
                log.error('PAPITYPE/OBJECT/IMAGETYP keyword not found')
                self.type = 'UNKNOW'
                raise Exception("Cannot classify (dark,flat, science) FITS image")
        else: #o2000, hawk-i?
            try:
                if myfits[0].header.has_key('PAPITYPE'):
                    self.type = myfits[0].header['PAPITYPE']
                elif myfits[0].header.has_key('IMAGETYP'):
                    self.type = myfits[0].header['IMAGETYP']
                elif myfits[0].header[keyword_with_frame_type].lower().count('dark') :
                    self.type = "DARK"
                elif myfits[0].header[keyword_with_frame_type].lower().count('lamp off'):
                    self.type = "DOME_FLAT_LAMP_OFF"
                elif myfits[0].header[keyword_with_frame_type].lower().count('lamp on'):
                    self.type = "DOME_FLAT_LAMP_ON"
                elif myfits[0].header[keyword_with_frame_type].lower().count('dusk'):
                    self.type = "TW_FLAT_DUSK"
                elif myfits[0].header[keyword_with_frame_type].lower().count('dawn'):
                    self.type = "TW_FLAT_DAWN"
                elif myfits[0].header[keyword_with_frame_type].lower().count('sky_flat') or \
                     myfits[0].header[keyword_with_frame_type].lower().count('flat'): 
                    self.type = "SKY_FLAT"
                elif myfits[0].header[keyword_with_frame_type].lower().count('sky'):
                    self.type = "SKY"
                elif myfits[0].header[keyword_with_frame_type].lower().count('focus'):
                    self.type = "SCIENCE"  
                    #por una razon que desconozco, CAHA le asigna el id 'focus' en algunas images, 
                    #pero tiene pinta que fue  un despiste del operador !!!
                elif myfits[0].header[keyword_with_frame_type].lower().count('science'):
                    self.type = "SCIENCE"
                else:
                    self.type = "SCIENCE"
                    #log.debug("DEFAULT Image type: %s"%self.type)
            except KeyError:
                log.error('PAPITYPE/OBJECT/IMAGETYP keyword not found')
                self.type = 'UNKNOW'
                raise Exception("Cannot classify (dark,flat, science) FITS image")
        
        #Is pre-reduced the image ? by default, no
        self.processed = False
        
        #print "File :"+ self.pathname
        #Filter
        try:
            if myfits[0].header['INSTRUME']=='HAWKI':
                if myfits[0].header.has_key('ESO INS FILT1 NAME'): 
                    self.filter = myfits[0].header['ESO INS FILT1 NAME']
                elif myfits[0].header.has_key('FILTER1'):
                    self.filter = myfits[0].header['FILTER1']
                elif myfits[0].header.has_key('FILTER2'):
                    self.filter = myfits[0].header['FILTER2']
            else: # PANIC, O2000, ...
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
            # WCS-coordinates are prefered than RA,DEC
            if (myfits[0].header.has_key('CTYPE1') 
                            and myfits[0].header['CTYPE1']=='RA---TAN'):
                
                wcs = pywcs.WCS(myfits[0].header)
                center_pix = numpy.array([[self.naxis1/2,self.naxis2/2]], numpy.float_)
                self._ra = wcs.wcs_pix2sky(center_pix, 1)[0][0] # ups, we are supposing naxis1 is RA axis
                log.debug("Read RA-WCS coordinate =%s", self._ra)
            elif myfits[0].header.has_key('RA'):
                self._ra = myfits[0].header['RA']
            else:
                raise Exception("No valid RA coordinate found")
        except Exception,e:
            log.warning('Error reading RA keyword :%s',str(e))
            self._ra  = -1
            
        #Dec-coordinate (in degrees)
        try:
            # WCS-coordinates are prefered than RA,DEC
            if (myfits[0].header.has_key('CTYPE2') 
                            and myfits[0].header['CTYPE2']=='DEC--TAN'):
                wcs = pywcs.WCS(myfits[0].header)
                center_pix = numpy.array([[self.naxis1/2,self.naxis2/2]], numpy.float_)
                self._dec = wcs.wcs_pix2sky(center_pix, 1)[0][1] # ups, we are supposing naxis2 is Declination axis
                log.debug("Read Dec-WCS coordinate =%s", self._dec)
            elif myfits[0].header.has_key('DEC'):
                self._dec = myfits[0].header['DEC']
            else:
                raise Exception("No valid DEC coordinates found")
        except Exception,e:
            log.warning('Error reading DEC keyword : %s', str(e))
            self._dec  = -1
    
        try:
            self.equinox = myfits[0].header['EQUINOX']
        except KeyError:
            log.warning("EQUINOX keyword not found")
            self.equinox=-1
            
        #MJD-Modified julian date 'days' of observation
        try:
            self.mjd = myfits[0].header['MJD-OBS']
        except KeyError:
            log.warning('MJD-OBS keyword not found')
            self.mjd  = -1
           
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
            ###log.warning('CHIPCODE keyword not found, setting default (1)')
            self.chipcode  = 1  # default
        
        #DetectorID
        try:
            if myfits[0].header['INSTRUME']=='HAWKI': 
                self.detectorID=myfits[0].header['HIERARCH ESO DET CHIP NAME']
            else: 
                self.detectorID='O2k'
        except:
              log.warning("Cannot find HIERARCH ESO DET CHIP NAME")
              self.detectorID='unknown'
               
        #RunID
        self.runID=-1
        
        #OB_ID : Observation Block Id
        try:
            if myfits[0].header['INSTRUME']=='HAWKI':
                self.obID = myfits[0].header['HIERARCH ESO OBS ID']
            elif myfits[0].header['INSTRUME']=='Omega2000':
                self.obID = myfits[0].header['POINT_NO'] # for O2000
            elif myfits[0].header['INSTRUME']=='Panic':
                # check how was observed
                if self.obs_tool:
                    self.obID = myfits[0].header['OB_ID'] # for PANIC using OT
                else:
                    self.obID = myfits[0].header['POINT_NO'] # for PANIC using MIDAS or whatever
            else:
                self.obID = -1
        except Exception,e:
            log.warning("Cannot find OB_ID keyword : %s:",str(e))
            self.obID = -1
               
        #OB_PAT : Observation Block Pattern
        try:
            if myfits[0].header['INSTRUME']=='HAWKI':
                self.obPat = myfits[0].header['HIERARCH ESO TPL ID']
            elif myfits[0].header['INSTRUME']=='Omega2000':
                self.obPat = myfits[0].header['POINT_NO'] # for O2000
            elif myfits[0].header['INSTRUME']=='Panic':
                if self.obs_tool:
                    self.obPat = myfits[0].header['OB_PAT'] # for PANIC using OT
                else:
                    self.obPat = myfits[0].header['POINT_NO'] # for PANIC using MIDAS or whatever
            else:
                self.obPat = -1
        except Exception,e:
            log.warning("Cannot find keyword : %s:",str(e))
            self.obPat = -1
                   
        #PAT_EXPN : Pattern Exposition Number (expono of noexp)
        try:
            if myfits[0].header['INSTRUME']=='HAWKI':
                self.pat_expno = myfits[0].header['HIERARCH ESO TPL EXPNO']
            elif myfits[0].header['INSTRUME']=='Omega2000':
                self.pat_expno = myfits[0].header['DITH_NO']
            elif myfits[0].header['INSTRUME']=='Panic':
                if self.obs_tool:
                    self.pat_expno = myfits[0].header['PAT_EXPN'] # for PANIC using OT
                else:
                    self.pat_expno = myfits[0].header['DITH_NO'] # for PANIC using MIDAS or whatever
            else:
                self.pat_expno = -1
        except Exception,e:
            log.warning("Cannot find keyword : %s:",str(e))
            self.pat_expno = -1
            
        #PAT_NEXP : Number of Expositions of Pattern (expono of noexp)
        try:
            if myfits[0].header['INSTRUME']=='HAWKI':
                self.pat_noexp = myfits[0].header['HIERARCH ESO TPL NEXP']
            elif myfits[0].header['INSTRUME']=='Omega2000':
                self.pat_noexp = -1 # not available
            elif myfits[0].header['INSTRUME']=='Panic':
                if self.obs_tool:
                    self.pat_noexp = myfits[0].header['PAT_NEXP'] # for PANIC using OT
                else:
                    # we could try to parse OBJECT key
                    self.pat_noexp = -1 # for PANIC using MIDAS or whatever
            else:
                self.pat_noexp = -1
        except Exception,e:
            log.warning("Cannot find keyword : %s:",str(e))
            self.pat_noexp = -1
            
        # To Fix PRESS1 and PRESS2 wrong keyword values of Omega2000 headers
        try:
            if myfits[0].header['INSTRUME']=='Omega2000':
                myfits[0].header.update('PRESS1', 0.0)
                myfits[0].header.update('PRESS2', 0.0)
        except Exception,e:
            log.warning("Keyword not found : %s ->", str(e))
        
        try:
            myfits.close(output_verify='ignore')
        except Exception, e:
            log.error("Error while closing FITS file %s   : %s",self.pathname, str(e))
        
            
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
            
        
        
    def print_info(self):
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
            dr.print_info()
            print dr.isTwFlat()
        except Exception, e:
            log.error("Error reading fits File %s", str(e))
    else:
        usage()

        
