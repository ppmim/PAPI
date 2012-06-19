###############
# @file ldac.py
# @author Douglas Applegate & Thomas Erben
# @date 9/2/2008
#
# @brief Utilities to make accessing LDAC cats easier
###############

# HISTORY INFORMATION:
# ====================
#
# 01.09.2010:
# I included treatment of slices through vectorkeys in
# LDACCat.__getitem__
#
# 09.09.2010:
# I made the module more robust against the non-existence
# of necessary libraries

# standard-library includes: 
import sys

# import non-standard modules:
try:
    import pyfits, numpy
except ImportError:
    sys.stderr.write("This script needs the modules 'pyfits' and 'numpy'!\n")
    sys.stderr.write("see http://www.stsci.edu/resources/software_hardware/pyfits\n")
    sys.stderr.write("and/or http://numpy.scipy.org/\n")
    sys.exit(1)

#from __future__ import with_statement

class LDACCat(object):

    def __init__(self, hdu):
        self.hdu = hdu
        self.current=0
        
    def __len__(self):
        return self.hdu.header['NAXIS2'] # this way we prevent the error when there aren't any data row (naxis2=0)
        #return len(self.hdu.data) # throw an error if data

    def __getitem__(self, key):

        if type(key) == type(5) or \
                type(key) == type(slice(5)):
            return self.hdu.data[key]

        if type(key) == type("a"):
            # we need to deal with slices through vector keys
            # such as 'MAG_APER(2)'
            startind = key.find("(")
            endind = key.find(")")

            if startind > 0 and endind > 0:
                keyname = key[:startind]
                keyindex = int(key[startind + 1:endind]) - 1

                try:
                    return self.hdu.data.field(keyname)[:,keyindex]
                except AttributeError:
                    raise KeyError(key) 
            else:
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
        return self
        #return self.hdu.data.__iter__()

    def __contains__(self, item):
        return item in self.keys()

    def has_key(self, key):
        return self.__contains__(key)

    def filter(self, mask):
        return LDACCat(pyfits.BinTableHDU(data=self.hdu.data[mask],
                                          header=self.hdu.header))
    def next(self):
        if self.current >= self.__len__(): # len(self.hdu.data):
            raise StopIteration
        else:
            self.current += 1
            return self.hdu.data[self.current -1]

    def saveas(self, file, clobber=False):
        self.hdu.writeto(file, output_verify='warn', clobber=clobber)


def openObjects(hdulist, table='LDAC_OBJECTS'):

    for hdu in hdulist:
        try:
            if table == hdu.header['EXTNAME']:
                return LDACCat(hdu)
        except KeyError:
            pass

    return None
    

def openObjectFile(filename, table='LDAC_OBJECTS'):
    hdulist = pyfits.open(filename)
    if hdulist is None:
        return None
    return openObjects(hdulist, table)

        
        
