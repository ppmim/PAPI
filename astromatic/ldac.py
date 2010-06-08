###############
# @file ldac.py
# @author Douglas Applegate
# @date 9/2/2008
#
# @brief Utilities to make accessing LDAC cats easier
###############

from __future__ import with_statement
import pyfits

class LDACCat(object):

    def __init__(self, hdu):
        self.hdu = hdu
        
    def __len__(self):
        return len(self.hdu.data)

    def __getitem__(self, key):

        if type(key) == type(5) or \
                type(key) == type(slice(5)):
            return self.hdu.data[key]

        if type(key) == type("a"):
            try:
                return self.hdu.data.field(key)
            except AttributeError:
                raise KeyError(key)

        raise TypeError

    def __setitem__(self, key, val):
        raise NotImplementedError

    def __delitem__(self, key):
        raise NotImplementedError

    def keys(self):
        return self.hdu.columns.names

    def __iter__(self):
        return self.hdu.data.__iter__()

    def __contains__(self, item):
        return item in self.keys()

    def has_key(self, key):
        return self.__contains__(key)

    def filter(self, mask):
        return LDACCat(pyfits.BinTableHDU(data=self.hdu.data[mask],
                                          header=self.hdu.header))

    def saveas(self, file):
        self.hdu.writeto(file)


def openObjects(hdulist, table='OBJECTS'):

    for hdu in hdulist:
        try:
            if table == hdu.header['EXTNAME']:
                return LDACCat(hdu)
        except KeyError:
            pass

    return None
    

def openObjectFile(filename, table='OBJECTS'):
    hdulist = pyfits.open(filename)
    if hdulist is None:
        return None
    return openObjects(hdulist, table)

        
        
