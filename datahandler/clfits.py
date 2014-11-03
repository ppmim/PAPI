#
# PANICtool
#
# DataClassifier.py
#
# Created    : 07/04/2008    jmiguel@iaa.es
# Update     : 25/05/2009    jmiguel@iaa.es   Added object field
#              14/12/2009    jmiguel@iaa.es   Renamed SKY_FLAT by TW_FLAT
#                                             Reanamed isSkyFlat() by isTwFlat()
#              02/03/2010    jmiguel@iaa.es   Added READMODE checking
#              20/08/2013    jmiguel@iaa.es   Addapted to support OSN CCDs. 
#
###############################################################################

"""
Module containing all basic operation with FITS files.
"""

# import external modules


import os.path
import sys
import time
from time import strftime

from astropy import wcs
import astropy.io.fits as fits

# Logging (PAPI or Python built-in)
try:
    # PAPI logging module
    from misc.paLog import log
except ImportError:
    import logging as log


import warnings


###############################################################################
class FitsTypeError(ValueError):
    """Raised when trying to classify a FITS file which is
    not supported or a special file (e.g. a non PANIC generated file)"""


###############################################################################
class ClFits (object):

    """
      Represents a Classified FITS of PANIC
      A class for PANIC FITS file classification and recognition. Actually it is 
      a wrapper for astropy.io.fits.

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
         MASTER_DARK_MODEL,
         MASTER_SKY_FLAT,
         MASTER_DOME_FLAT,
         MASTER_TW_FLAT,
         SCIENCE_REDUCED
         
    """

    # Class initialization
    def __init__(self, full_pathname, check_integrity=True, *a,**k):
      
        """
        Init the object
          
        Parameters
        ----------
        full_pathname: str
            pathname The full filename of the fits file to classify.
        check_integrity: bool
            When True, the FITS integrity is done to check the file is complete.
            Mainly used on QL the know whether file writting finished. 
        """
        
        super (ClFits, self).__init__ (*a,**k)
      
        self.check_integrity = check_integrity
        self.pathname = full_pathname # path_name + root_name
        self.filename = os.path.basename(self.pathname)  # root_name
        self.type = "UNKNOWN"
        self.naxis1 = -1
        self.naxis2 = -1
        
        # shape of first HDU of the FITS =(naxis3, naxis2, naxis1)
        self._shape = ()
        
        # exposition ID for the current night (unique for that night)
        self.runID = -1
        
        # Observation Block ID (unique) provided by the OT
        self.obID = -1
        
        # Observation (dither) Pattern 
        self.obPat = -1 
        
        # Exposure number within (dither) pattern
        self.pat_expno = -1 
        
        # Number of exposures within pattern
        self.pat_noexp = -1 

        self.instrument = ""
        self.processed = False
        self.exptime = -1
        self.filter = ""
        self.mef = False
        self.next = 0
        self.datetime_obs = ""
        self.detectorID = ""
        self.date_obs = ""
        self.time_obs = ""
        self._ra = -1
        self._dec = -1
        self.equinox = -1
        
        # modified julian date of observation
        self.mdj = -1 

        self.object = ""
        self.chipcode = 1
        self.ncoadds = -1
        self.itime = 0.0
        self.readmode = ""
        self.my_header = None
        self.obs_tool = False
        self.pix_scale = 1.0
        self.binning = 1
        self.telescope = "unknown"
        self._softwareVer = ''
        
        
        self.recognize()

    def getType(self, distinguish_domeflat=True):
        
        # Because in some special case (on-line PQL detecting sequences)
        # we might need to not distinguish between DOME_FLAT_ON and DOME_FLAT_OFF
        if distinguish_domeflat and self.isDomeFlat():
            return "DOME_FLAT"
        else:
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
        return (self.type=="FOCUS") 
    
    def isScience(self):
        return (self.type.count("SCIENCE") or self.type.count("STD"))
    
    def isRawScience(self):
        return (self.type=="SCIENCE_RAW")

    def isReducedScience(self):
        return (self.type=="SCIENCE_REDUCED")
    
    def isMasterDark(self):
        return (self.type=="MASTER_DARK")
    
    def isMasterDarkModel(self):
        return (self.type=="MASTER_DARK_MODEL")
    
    def isMasterFlat(self):
        return (self.type=="MASTER_DOME_FLAT" or self.type=="MASTER_TW_FLAT"
                or self.type=="MASTER_SKY_FLAT")
    
    def isSky(self):
        return (self.type=="SKY")
    
    def isObject(self):
        return (self.isScience())
    
    def isMEF(self):
        return (self.mef)

    def getNExt(self):
        return (self.next)
    
    def isFromOT(self):
        return (self.obs_tool)
    
    def isFromGEIRS(self):
        if self.softwareVer.count("GEIRS"):
            return True
        else:
            return False
    
    def isFromPANIC(self):
        return self.instrument=='panic'

    def isPANICFullFrame(self):
        return (self.instrument=='panic' and self.naxis1>=4096 and self.naxis2>=4096)
                
    def expTime(self):
        return (self.exptime)
    
    @property
    def shape(self):
        return self._shape
    
    @property
    def ra(self):
        return self._ra
    
    @property
    def dec(self):
        return self._dec
    
    @property
    def pixScale(self):
        return pix_scale

    @property
    def Telescope(self):
        return self.telescope
        
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
        if ext!=0 and not self.mef:
            log.error("Wrong extension number especified. Not a MEF file.")
        else:
            try:
                myfits = fits.open(self.pathname, 
                                     ignore_missing_end=True) # since some problems with O2k files   
                temp = myfits[0].data
                myfits.close(output_verify='ignore')
                return temp
            except:
                log.error('Could not open frame - something wrong with input data')
                raise
    

    def recognize(self, retries=5):
     
        log.info("Recognizing file %s" %self.filename)

        # Check the file exists
        if not os.path.exists( self.pathname ):
            log.error('Cannot find frame : "%s"' % self.pathname)
            raise Exception("File %s not found"%self.pathname)
        
        # Read FITS and check FITS-file integrity
        nTry = 0
        found_size = 0

        if self.check_integrity: 

            # We change that behavior of fits warnings with a filter.
            # See http://bit.ly/1etvfJC
            # It is very useful because some fits warnings are quite usual
            # when reading non completely written files. Non captured warnings 
            # slow down the reading process.
            # Ex. 
            #    "Header block contains null bytes instead of spaces for padding, and is not 
            #     FITS-compliant. Nulls may be replaced with spaces upon writing."
            # 
            #log.debug("Cheking FITS integrity")
            # Turn matching warnings into exceptions
            warnings.simplefilter('error', UserWarning)
            while True:
                log.debug("FITS integrity check. FILE=%s ITER=%d"%(self.pathname,nTry))
                try:
                    # First level of checking
                    # found_size = fits_simple_verify(self.pathname)
                    # Now, try to read the whole FITS file
                    myfits = fits.open(self.pathname, mode='readonly', memmap=True,
                                     ignore_missing_end=False) # since some problems with O2k files 
                except Exception, e:
                    log.warning("Error reading FITS : %s"%self.pathname)
                    if nTry<retries:
                        nTry +=1
                        time.sleep(nTry*1.0)
                        log.warning("Error reading FITS. Trying to read again file : %s\n %s"%(self.pathname, str(e)))
                    else:
                        log.error("Finally, FITS-file could not be read with data integrity:  %s\n %s"%(self.pathname, str(e)))
                        log.error("File discarded : %s"%self.pathname)
                        raise e
                else:
                    break

            # Undo the warning filter to avoid Exceptions forever even when 
            # non-hard warnings!
            warnings.resetwarnings()
        else:
            myfits = fits.open(self.pathname, 
                                     ignore_missing_end=True) # since some problems with O2k files     

        # Check if is a MEF file 
        if len(myfits)>1:
            self.mef = True
            self.next = len(myfits)-1
            log.debug("Found a MEF file with %d extensions", self.next)
        else:
            self.mef = False
            self.next = 1
            ###log.debug("Found a simple FITS file")
        
        # If file is a MEF, some values will be read from the header extensions      
        if self.mef:
            # we suppose all extension have the same dimensions (subwindows ???)
            self.naxis1 = myfits[1].header['NAXIS1']
            self.naxis2 = myfits[1].header['NAXIS2']
            #self._shape = myfits[1].data.shape # (naxis3, naxis2, naxis1)
            # Because the fits.data.shape makes extremealy slow read
            # of the FITS-file (mainly because shape need to read the whole image),
            # we build a shape tuple from the header values (NAXIS1, NAXIS2, and NAXIS3 )
            if 'NAXIS3' in myfits[1].header:
                self._shape = (myfits[1].header['NAXIS3'], self.naxis2, self.naxis1)
            else:
                self._shape = (self.naxis2, self.naxis1)
        else:
            self.naxis1 = myfits[0].header['NAXIS1']
            self.naxis2 = myfits[0].header['NAXIS2']
            #self._shape = myfits[0].data.shape # (naxis2, naxis1)
            #see comments about fits.data.shape above 
            if 'NAXIS3' in myfits[0].header:
                self._shape = (myfits[0].header['NAXIS3'], self.naxis2, self.naxis1)
            else:
                self._shape = (self.naxis2, self.naxis1)

        # pointer to the primary-main header
        self.my_header = myfits[0].header
          
        
         
        # INSTRUMENT
        try:
            if 'INSTRUME' in myfits[0].header:
                self.instrument = myfits[0].header['INSTRUME'].lower()
            else:
                self.instrument = "Unknown"
        except Exception,e:
            log.warning("INSTRUME keyword not found")
            self.instrument = "Unknown"
        
        # Find out the how data file were"observed"
        if 'OBS_TOOL' in myfits[0].header or self.instrument=='hawki':
            self.obs_tool = True
        else:
            self.obs_tool = False
        
        # Software Version (GEIRS Version):
        # Old versions of GEIRS used 'SOFTWARE'
        if 'SOFTWARE' in myfits[0].header:
            self._softwareVer = myfits[0].header['SOFTWARE']
        # New versions of GEIRS moved to CREATOR keyword for software version
        if 'CREATOR' in myfits[0].header:
            self._softwareVer = myfits[0].header['CREATOR']
        
        
        # IMAGE TYPE
        try:
            if self.instrument=='omega2000' and 'PAPITYPE' in myfits[0].header:
                keyword_with_frame_type = 'PAPITYPE'
                # It happens if the image is product of PAPI 
            elif self.instrument=='omega2000' and 'IMAGETYP' in myfits[0].header:
                keyword_with_frame_type = 'IMAGETYP'
            elif self.instrument=='omega2000' and 'OBJECT' in myfits[0].header:
                keyword_with_frame_type = 'OBJECT'
            elif self.instrument=='hawki' and 'IMAGETYP' in myfits[0].header:
                keyword_with_frame_type = 'IMAGETYP'
            elif self.instrument=='hawki' and 'OBJECT' in myfits[0].header:
                keyword_with_frame_type = 'OBJECT'
            elif self.instrument=='panic': # current ID in GEIRS for PANIC
                if self.obs_tool:
                    keyword_with_frame_type = 'IMAGETYP'
                    #keyword_with_frame_type = 'OBJECT'
                else:
                    keyword_with_frame_type = 'OBJECT'
            else: keyword_with_frame_type = 'IMAGETYP' #default, even for Roper    
        except KeyError:
            log.warning('IMAGETYP or OBJECT keywords not found')
            keyword_with_frame_type = 'IMAGETYP' #default
            
        # Find out the type of frame ( BIAS, DARK, DOME_FLAT_LAMP_ON/OFF, 
        # SKY_FLAT, SCIENCE , MASTER_calibration, UNKNOW)     
        # IMAGETYP = [BIAS, DARK, LAMP_ON_FLAT, LAMP_OFF_FLAT, TW_FLAT_DUSK, 
        #               TW_FLAT_DAWN, SKY_FLAT, SCIENCE, SKY, STD, FOCUS ]
        # FIELDTYP = [POINTLIKE, SPARSE_FIELD, CROWDED_FIELD, EXT_OBJECT ]
        if self.instrument =='panic':
            try:
                # Self-typed file, created by PAPI
                if 'PAPITYPE' in myfits[0].header:
                    self.type = myfits[0].header['PAPITYPE']
                else:
                    #
                    if 'IMAGETYP' in myfits[0].header:
                        ltype = myfits[0].header['IMAGETYP'].lower()
                    else:
                        ltype = myfits[0].header[keyword_with_frame_type].lower()
                        
                    if ltype.count('dark') :
                        self.type = "DARK"
                    elif ltype.count('lamp_off'):
                        self.type = "DOME_FLAT_LAMP_OFF"
                    elif ltype.count('lamp_on'):
                        self.type = "DOME_FLAT_LAMP_ON"
                    elif ltype.count('dusk'):
                        self.type = "TW_FLAT_DUSK"
                    elif ltype.count('dawn'):
                        self.type = "TW_FLAT_DAWN"
                    elif ltype.count('sky_flat'): 
                        self.type = "SKY_FLAT"
                    elif ltype.count('flat'):
                        self.type = "DOME_FLAT"
                    elif ltype.count('sky'):
                        self.type = "SKY"
                    elif ltype.count('focus'):
                        self.type = "FOCUS"
                    elif ltype.count('object'):
                        self.type = "SCIENCE"
                    else:
                        # By default, the image is classified as SCIENCE object
                        self.type = "SCIENCE"
            except KeyError:
                log.warning('PAPITYPE/OBJECT/IMAGETYP keyword not found')
                self.type = 'UNKNOW'
        elif self.instrument=='hawki':
            try:
                if myfits[0].header[keyword_with_frame_type].lower().count('dark') :
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
                elif myfits[0].header[keyword_with_frame_type].lower().count('object'):
                    self.type = "SCIENCE"
                else:
                    #By default, the image is classified as SCIENCE object
                    self.type = "SCIENCE"
            except KeyError:
                log.warning('PAPITYPE/OBJECT/IMAGETYP keyword not found')
                self.type = 'UNKNOW'
                
        else: #o2000, Roper or  ??
            try:
                if myfits[0].header[keyword_with_frame_type].lower().count('master'):
                    self.type = myfits[0].header[keyword_with_frame_type]
                elif myfits[0].header[keyword_with_frame_type].lower().count('bias'):
                    self.type = "BIAS"
                elif myfits[0].header[keyword_with_frame_type].lower().count('dark'):
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
                    self.type = "FOCUS"  
                elif myfits[0].header[keyword_with_frame_type].lower().count('science'):
                    self.type = "SCIENCE"
                # CCD Roper (OSN)
                elif myfits[0].header[keyword_with_frame_type].lower().count('light'):
                    self.type = "SCIENCE"
                else:
                    #By default, the image is classified as SCIENCE object
                    self.type = "SCIENCE"
                    #log.debug("DEFAULT Image type: %s"%self.type)
            except KeyError:
                log.warning('PAPITYPE/OBJECT/IMAGETYP keyword not found')
                self.type = 'UNKNOW'
        
        #Is pre-reduced the image ? by default, isn't
        self.processed = False
        
        #print "File :"+ self.pathname
        #FILTER
        try:
            if self.instrument=='hawki':
                if 'ESO INS FILT1 NAME' in myfits[0].header: 
                    self.filter = myfits[0].header['ESO INS FILT1 NAME']
                elif 'FILTER1' in myfits[0].header:
                    self.filter = myfits[0].header['FILTER1']
                elif 'FILTER2' in myfits[0].header:
                    self.filter = myfits[0].header['FILTER2']
            else: # PANIC, O2000, Roper, ...
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
            self.itime = myfits[0].header['ITIME']
        except KeyError:
            log.warning('ITIME keyword not found')
            self.itime  = -1
            
        #Number of coadds
        try:
            if 'NCOADDS' in myfits[0].header:
                self.ncoadds = myfits[0].header['NCOADDS']
            elif 'NDIT' in myfits[0].header:
                self.ncoadds = myfits[0].header['NDIT']
            elif 'HIERARCH ESO DET NDIT' in myfits[0].header:
                self.ncoadds = myfits[0].header['HIERARCH ESO DET NDIT']
            else:
                self.ncoadds = 1
        except KeyError:
            log.warning('NCOADDS keyword not found. Taken default value (=1)')
            self.ncoadds  = 1
                 
        #Read-Mode
        try:
            self.readmode = myfits[0].header['READMODE']
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
        
        # ############################
        # RA-coordinate (in degrees)
        # ############################
        try:
            # WCS-coordinates are preferred than RA,DEC (both in degrees)
            if ('CTYPE1' in myfits[0].header and
                     'TAN' in myfits[0].header['CTYPE1'] ): #'RA---TAN' or 'RA---TAN--SIP'
                m_wcs = wcs.WCS(myfits[0].header)
                #self._ra, self._dec = wcs.image2sky( self.naxis1/2, self.naxis2/2, True)
                # No SIP or Paper IV table lookup distortion correction is applied.
                self._ra = m_wcs.wcs_pix2world([[self.naxis1/2, self.naxis2/2]], 1)[0][0]
                #log.debug("Read RA-WCS coordinate =%s", self._ra)
            elif 'RA' in myfits[0].header:
                self._ra = myfits[0].header['RA'] # degrees supposed
            elif 'OBJCTRA' in myfits[0].header:
                # Mainly for Roper CCDs
                a = myfits[0].header['OBJCTRA']
                self._ra = float(a.split()[0]) + float(a.split()[1])/60.0 + float(a.split()[2])/3600.0
                # convert to degrees                
                self._ra = self._ra * 360.0 / 24.0
            else:
                raise Exception("No valid RA coordinate found")
        except Exception,e:
            log.warning('Error reading RA keyword :%s',str(e))
            self._ra  = -1
        finally:
            #log.debug("RA = %s"%str(self._ra))
            pass

        # ############################
        # Dec-coordinate (in degrees)
        # ############################
        try:
            # WCS-coordinates are preferred than RA,DEC
            if ('CTYPE2' in myfits[0].header and
                     'TAN' in myfits[0].header['CTYPE2']): #=='DEC--TAN' or 'DEC--TAN--SIP'
                m_wcs = wcs.WCS(myfits[0].header)
                #self._ra, self._dec = wcs.image2sky( self.naxis1/2, self.naxis2/2, True)
                # No SIP or Paper IV table lookup distortion correction is applied.
                self._dec = m_wcs.wcs_pix2world([[self.naxis1/2, self.naxis2/2]], 1)[0][1]
                # log.debug("Read Dec-WCS coordinate =%s", self._dec)
            elif 'DEC' in myfits[0].header:
                self._dec = myfits[0].header['DEC']
            elif 'OBJCTDEC' in myfits[0].header:
                a = myfits[0].header['OBJCTDEC']
                if a.split()[0][0]=='-':
                    self._dec = (-1.0) * float(a.split()[0]) + float(a.split()[1])/60.0 + float(a.split()[2])/3600.0
                    self._dec*=-1.0
                else:
                    self._dec =  float(a.split()[0]) + float(a.split()[1])/60.0 + float(a.split()[2])/3600.0
            else:
                raise Exception("No valid DEC coordinates found")
        except Exception,e:
            log.warning('Error reading DEC keyword : %s', str(e))
            self._dec  = -1
        finally:
            #log.debug("DEC = %s"%str(self._dec))
            pass

        try:
            self.equinox = myfits[0].header['EQUINOX']
        except KeyError:
            log.warning("EQUINOX keyword not found")
            self.equinox = -1
            
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
            if self.instrument=='hawki': 
                self.detectorID = myfits[0].header['HIERARCH ESO DET CHIP NAME']
            elif self.instrument=='panic':
                self.detectorID = 'HAWAII-2RG'
            else:
                self.detectorID = 'UNKNOWN'
        except:
            log.warning("Cannot find HIERARCH ESO DET CHIP NAME")
            self.detectorID='unknown'
               
        #RunID
        self.runID=-1
        
        #OB_ID : Observation Block Id
        try:
            if self.instrument=='hawki':
                self.obID = myfits[0].header['HIERARCH ESO OBS ID']
            elif self.instrument=='omega2000':
                self.obID = myfits[0].header['POINT_NO'] # for O2000
            elif self.instrument=='panic':
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
            if self.instrument=='hawki':
                self.obPat = myfits[0].header['HIERARCH ESO TPL ID']
            elif self.instrument=='omega2000':
                self.obPat = myfits[0].header['POINT_NO'] # for O2000
            elif self.instrument=='panic':
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
            if self.instrument=='hawki':
                self.pat_expno = myfits[0].header['HIERARCH ESO TPL EXPNO']
            elif self.instrument=='omega2000':
                self.pat_expno = myfits[0].header['DITH_NO']
            elif self.instrument=='panic':
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
            if self.instrument=='hawki':
                self.pat_noexp = myfits[0].header['HIERARCH ESO TPL NEXP']
            elif self.instrument=='omega2000':
                self.pat_noexp = -1 # not available
            elif self.instrument=='panic':
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
            
        # Try Fix PRESS1 and PRESS2 wrong keyword values of Omega2000 headers
        # Has no sense, because file is not opened with 'update' flag
        try:
            if self.instrument=='omega2000':
                myfits[0].header.update('PRESS1', 0.0)
                myfits[0].header.update('PRESS2', 0.0)
        except Exception,e:
            log.warning("Keyword not found : %s ->", str(e))


        # Telescope
        if 'TELESCOPE' in myfits[0].header:
            self.telescope = myfits[0].header['TELESCOPE']
        else:
            self.telescope = "unknown"
        
        # Binning 
        if 'BINNING' in myfits[0].header:
            self.binning = int(myfits[0].header['BINNING'])
        elif 'XBINNING' in myfits[0].header:
            self.binning = int(myfits[0].header['XBINNING'])
        else:
            self.binning = 1

        # Pixel scale (updated with binning factor)
        if 'PIXSCALE' in myfits[0].header:
            self.pix_scale = float(myfits[0].header['PIXSCALE']) * self.binning
        else:
            if self.instrument=='omega2000':
                self.pix_scale = 0.45 * self.binning
            elif self.instrument=='panic':
                self.pix_scale = 0.45 * self.binning
            elif self.instrument=='roper':  # roperT150
                self.pix_scale = 0.23 * self.binning
            elif self.instrument=='ropert90':
                self.pix_scale = 0.386 * self.binning
            elif self.instrument=='hawki':
                self.pix_scale = 0.106 * self.binning
            else:
                # default scale ?
                self.pix_scale = 0.45 * self.binning

        # 
        # Close the file. Some updates can be done (PRESS1, ...) but it mustn't
        #            
        try:
            myfits.close(output_verify='ignore')
        except Exception, e:
            log.error("Error while closing FITS file %s   : %s",
                      self.pathname, str(e))
        
        log.debug("End of FITS recognition: %s"%self.pathname)
        

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
        print "OB_ID     : ", self.getOBId()
        print "NoExp     : ", self.getNoExp()
        print "ExpNo     : ",self.getExpNo()
        print "---------------------------------"
      
################################################################################            
#  Useful function to check data integrity
################################################################################
def checkDataProperties( file_list, c_type=True, c_filter=True, c_texp=True, 
                         c_ncoadds=True, c_readmode=True):
    """
    This function will check all the files in the file_list have the same 
    properties required as True in the parameters.
    Note that that properties should be available in the main header(0) in case 
    of MEF files.
    """
    
    m_type = ''
    m_filter = ''
    m_texp = ''
    m_ncoadds = ''
    m_readmode = ''
    
    # First file as reference
    if len(file_list)>0 and file_list[0]:
            f = ClFits( file_list[0] )
            m_type = f.getType()
            m_filter = f.getFilter()
            m_texp = f.expTime()
            m_ncoadds = f.getNcoadds()
            m_readmode = f.getReadMode()
    
    # Check all files        
    for file in file_list:
            f = ClFits ( file )
            debug = 0
            if debug:
                print 'FILE=',file
                print '------------------------'
                print 'TYPE=', f.getType()
                print 'FILTER=', f.getFilter()
                print 'TEXP=', f.expTime()
                print 'NCOADDS=', f.getNcoadds()
                print 'READMODE=', f.getReadMode()
            
            if (  (c_type and m_type!=f.getType()) or 
                  (c_filter and m_filter!=f.getFilter()) or 
                  (c_texp and m_texp!=f.expTime()) or 
                  (c_ncoadds and m_ncoadds!=f.getNcoadds()) or 
                  (c_readmode and m_readmode!=f.getReadMode())):
                
                log.debug("Missmath on some property (FILTER, EXPTIME, NCOADDS or  READMODE)")
                return False
        
    log.debug("Successful properties checking")
    
    return True      

################################################################################
####### fits tools #############################################################
################################################################################
def isaFITS(filepath):
    """
    Check if a given filepath is a FITS file
    """
    
    if os.path.exists(filepath):
        try:
            fd = fits.open(filepath, ignore_missing_end=True)
            if fd[0].header['SIMPLE']==True:
                return True
            else:
                return False
        except Exception:
            return False
    else:
        return False

def fits_simple_verify(fitsfile):
    """
    Performs 2 simple checks on the input fitsfile, which is a string
    containing a path to a FITS file.  First, it checks that the first card is
    SIMPLE, and second it checks that the file 2880 byte aligned.
    
    This function is useful for performing quick verification of FITS files.
    
    Raises:
      ValueError:  if either of the 2 checks fails
      IOError:     if fitsfile doesn't exist

    Returns
    -------
    file_size: the current size of the fitsfile just read.  
    """
    
    if not os.path.exists(fitsfile):
        raise IOError("file '%s' doesn't exist" % fitsfile)


    f = open(fitsfile, "readonly")
        
    FITS_BLOCK_SIZE = 2880

    try:
        # check first card name
        card = f.read(len("SIMPLE"))
        if card != "SIMPLE":
            raise ValueError("input file is not a FITS file")

        # check file size
        file_size = os.stat(fitsfile).st_size
        #time.sleep(0.1)
        #while file_size != os.stat(fitsfile).st_size:
        #    time.sleep(0.1)
        #    file_size = os.stat(fitsfile).st_size

        # check that file_size>fits_block_size*5 to be sure all the header/s content can be read     
        if file_size<FITS_BLOCK_SIZE*4 or file_size % FITS_BLOCK_SIZE != 0:
            log.warning("FITS file is not 2880 byte aligned (corrupted?) or file_size too small")
            raise ValueError("FITS file is not 2880 byte aligned (corrupted?) or file_size too small")
    # Exceptions are re-raised after the finally clause has been executed
    finally:
        f.close()
            
    return file_size
		
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
            dr = ClFits(sys.argv[1])
            dr.print_info()
            print dr.isTwFlat()
        except Exception, e:
            log.error("Error reading fits File %s", str(e))
    else:
        usage()

        
