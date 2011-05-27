#!/usr/bin/env python
"""Module to do query to on-line catalogs"""

################################################################################
#
#
# PAPI (PANIC PIpeline)
#
# mef.py
#
# Multi Extension FITS file basic operations
#
# Created    : 23/05/2011    jmiguel@iaa.es -
# Last update: 
# TODO
#       
################################################################################

################################################################################
# Import necessary modules

import urllib, os
from optparse import OptionParser


# Logging
from misc.paLog import log



class ICatalog (object):
    """ Class to query on-line catalogs through Gator, the IRSA's 
        Catalog Search engine 
    """
    
    cat_names = {'2MASS':'fp_psc', 
                     'USNOB1': 'usno_b1',
                     'IRAS': 'iraspsc'
                     }
    outfmt = {'votable': 3,
              'ascii': 1
              }
    
    url = "http://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query"
    
    def __init__(self, *a, **k):
        """ The constructor """
        super (ICatalog, self).__init__ (*a, **k)
        
    
    def queryCatalog(self, ar, dec, sr=0.1, cat_name=None, 
                     out_filename=None, out_format='votable'):    
        """
        @summary: Query the catalog and return the output file format selected
        @param ra,dec: Right Ascension & Declination of center of conesearch
        @param sr: search radious (arcsecs)
        @param cat_name: catalog name where search will be done
        @param out_filename: filename where results will be saved;if absent, 
                the location will be a tempfile with a generated name
        @param out_format: format of the output generated; current options available are:
            - VO Table (XML) (votable) (default)
            - SVC (Software handshaking structure) message (svc)
            - ASCII table (ascii)
        @return: filename where results where saved (VOTABLE, ASCII_TABLE, ...)
        """
        
        params={}
        params['outfmt'] = ICatalog.outfmt[out_format]
        params['objstr'] = str(ar) + "+" 
        if float(dec)>0:
            params['objstr'] = str(ar) + "+" + str(dec)
        else:
            params['objstr'] = str(ar) + str(dec) 
        params['spatial'] = 'Cone'
        params['radius'] = sr
        params['catalog'] = cat_name
        
        query = urllib.urlencode(params)
        get_url = ICatalog.url + "?" + query

        #query = self.url + "?" + out_fmt + "&" +  radius + "&" + objstr + "&" \
        #        + spatial + "&" + catalog 
        
        
        log.debug("Query: %s", get_url)        

        # check if file exist
        if os.path.exists(out_filename):
            log.debug("Warning, overwriting file %s" %out_filename)
        
        try:
            urllib.urlcleanup()    
            r = urllib.urlretrieve(get_url, out_filename)
        except urllib.ContentTooShortError, e:
            log.error("Amount of data available was less than the expected \
            amount; might download was interrupted : %s", e)
            raise e
        else:
            return r
        

################################################################################
# main
################################################################################
if __name__ == "__main__":
    log.debug( 'Testing ICatalog')
    
    icat = ICatalog ()
    res_file = None
    # RA      =           243.298750 / (deg) R.A.:  16:13:11.7
    # DEC     =            54.600278 / Dec.:  54:36:01.0

    try:
        res_file = icat.queryCatalog(243.298750, +54.600278, 500, ICatalog.cat_names['2MASS'], 
                      "/tmp/prueba.xml", 'votable')[0]
        log.debug("Output file generated : %s", res_file) 
    except Exception,e:
        log.error("Sorry, cann't solve the query")
        raise e
            
