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

# if you comment out this line, it will download to the directory from which you run the script.
os.chdir('/directory/to/save/the/file/to')

url = 'http://www.mydomain.com/myfile.txt'
urllib.urlretrieve(url)



class ICatalog (object):
    """ Class to query on-line catalogs through Gator, the IRSA's 
        Catalog Search engine 
    """
    def __init__(self, *a, **k):
        """ The constructor """
        super (ICatalog, self).__init__ (*a, **k)
        
        self.cat_names = {'2MASS':'fp_psc', 
                     'USNOB1': 'usno_b1',
                     'IRAS': 'iraspsc'
                     }
        self.url = "http://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query"
    
    def queryCatalog(self, ar, dec, sr=0.1, cat_name=None, 
                     out_filename=None, out_format='votable',):    
        """Query the catalog and return the output file format selected"""
        
        
        os.chdir(os.path.dirname(out_filename))
        outmft = "outfmt=" + out_format
        objstr = str(ar) + str(dec)
        spatial = "spatial=Cone"
        radius = "radius=1" #arcsecs
        catalog = "catalog=" + cat_name
        query = self.url + "?" + outfmt + "&" + objstr + "&" + spatial + "&" + \
                catalog 

        try:
            urllib.urlretrieve(url, out_filename)
        except ContentTooShortError, e:
            log.error("Amount of data available was less than the expected \
            amount; might download was interrupted : %s", e)
            raise e
        else:
            return out_filename
        
        
