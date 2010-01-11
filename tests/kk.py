# Pyraf modules
import pyraf
from pyraf import iraf
from iraf import noao
from iraf import imred
from iraf import ccdred


import pyfits
import numpy


f=pyfits.open('/tmp/c1.fits',memmap=0)
bpm=numpy.zeros([2048,2048],dtype='int16')
#cdata=f[0].data.copy()
cdata=f[0].data
f.close()
print 'start...'
bpm[ (cdata<1) | (cdata>1.e6) ] = 1

#bpm = numpy.where(cdata<1,1,numpy.where(cdata>100000,1,cdata))

hdu = pyfits.PrimaryHDU()
#hdu.scale('int16') # importat to set first data type
hdu.data=bpm     
hdulist = pyfits.HDUList([hdu])
hdu.header.update('OBJECT','MASTER_PIXEL_MASK')
hdulist.writeto("/tmp/a.fits")
hdulist.close(output_verify='ignore')

"""for i in range(0,2048):
    for j in range(0,2048):
        if cdata[i,j]<1 or cdata[i,j]>100000:
            bpm[i,j]=1    

"""
"""
lista='/disk-a/caha/panic/DATA/SIMU_PANIC_2/cali_0041.fits , /disk-a/caha/panic/DATA/SIMU_PANIC_2/cali_0042.fits'
iraf.flatcombine(input=lista, 
                        output='/tmp/c1.fits', 
                        combine='median', 
                        ccdtype='none', 
                        process='no', 
                        reject='sigclip', 
                        subsets='no', 
                        scale='mode'
                        )
"""                 