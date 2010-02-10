#!/usr/bin/env python
################################################################################
#
# PANICtool
#
# quadSplit.py
#
# Created    : 20/01/2010    jmiguel@iaa.es
# Last update: 
#
################################################################################

# Import necessary modules
import pyfits
import numpy as np
import scipy.ndimage as nd

import os
import fileinput
from math import cos

def makeSplit(files, dst_path, gap):
  
    os.mkdir(dst_path+"/Q1")
    os.mkdir(dst_path+"/Q2")
    os.mkdir(dst_path+"/Q3")
    os.mkdir(dst_path+"/Q4")
    
    for file in files:
        src_img=pyfits.open(file)
        src=src_img[0].data
    

        w, h = src.shape
        #dest = np.empty(src_w, src_h)
      
        try:
            print "\nSplitting image with a gap of %d pixels"%(2*gap)
            print "Original RA=", src_img[0].header['RA']
            
            #gap=32 # in fact, 32*2
            scale = 0.45 # arcsec/pixel
            hdu = pyfits.PrimaryHDU()
            #hdu.scale('float32') # importat to set first data type
            hdu.data=src[0:h/2-gap,0:w/2-gap].copy()     
            hdulist = pyfits.HDUList([hdu])
            hdu.header=src_img[0].header.copy()
            hdu.verify('fix')
            print "RA1=",src_img[0].header['RA']+((scale*(w+2*gap)/4)/3600.0)/cos(src_img[0].header['DEC']*2*3.1415/360.0)
            print "COS=", cos(src_img[0].header['DEC']*2*3.1415/360.0)
            hdu.header.update("RA", src_img[0].header['RA']+((scale*(w+2*gap)/4)/3600.0)/cos(src_img[0].header['DEC']*2*3.1415/360.0))
            hdu.header.update("DEC", src_img[0].header['DEC']-((scale*(h+2*gap)/4)/3600.0))
            hdu.header.add_history('Warp image created from %s' % file)
            hdulist.writeto(dst_path+"/Q1/"+os.path.basename(file).replace(".fits","_q1.fits"))
            hdulist.close(output_verify='ignore')
            
            hdu = pyfits.PrimaryHDU()
            hdu.scale('float32') # importat to set first data type
            hdu.data=src[0:h/2-gap,w/2+gap:w]     
            hdulist = pyfits.HDUList([hdu])
            hdu.header=src_img[0].header.copy()
            hdu.verify('fix')
            print "RA2=",src_img[0].header['RA']-((scale*(w+2*gap)/4)/3600.0)/cos(src_img[0].header['DEC']*2*3.1415/360.0)
            hdu.header.update("RA", src_img[0].header['RA']-((scale*(w+2*gap)/4)/3600.0)/cos(src_img[0].header['DEC']*2*3.1415/360.0))
            hdu.header.update("DEC", src_img[0].header['DEC']-((scale*(h+2*gap)/4)/3600.0))
            hdu.header.add_history('Warp image created from %s' % file)
            hdulist.writeto(dst_path+"/Q2/"+os.path.basename(file).replace(".fits","_q2.fits"))
            hdulist.close(output_verify='ignore')
            
            
            hdu = pyfits.PrimaryHDU()
            hdu.scale('float32') # importat to set first data type
            hdu.data=src[h/2+gap:h,w/2+gap:w]     
            hdulist = pyfits.HDUList([hdu])
            hdu.header=src_img[0].header.copy()
            hdu.verify('fix')
            print "RA3=",src_img[0].header['RA']-((scale*(w+2*gap)/4)/3600.0)/cos(src_img[0].header['DEC']*2*3.1415/360.0)
            hdu.header.update("RA", src_img[0].header['RA']-((scale*(w+2*gap)/4)/3600.0)/cos(src_img[0].header['DEC']*2*3.1415/360.0))
            hdu.header.update("DEC", src_img[0].header['DEC']+((scale*(h+2*gap)/4)/3600.0))
            hdu.header.add_history('Warp image created from %s' % file)
            hdulist.writeto(dst_path+"/Q3/"+os.path.basename(file).replace(".fits","_q3.fits"))
            hdulist.close(output_verify='ignore')
        
        
            hdu = pyfits.PrimaryHDU()
            hdu.scale('float32') # importat to set first data type
            hdu.data=src[h/2+gap:h,0:w/2-gap]     
            hdulist = pyfits.HDUList([hdu])
            hdu.header=src_img[0].header.copy()
            hdu.verify('fix')
            print "RA4=",src_img[0].header['RA']+((scale*(w+2*gap)/4)/3600.0)/cos(src_img[0].header['DEC']*2*3.1415/360.0)
            hdu.header.update("RA", src_img[0].header['RA']+((scale*(w+2*gap)/4)/3600.0)/cos(src_img[0].header['DEC']*2*3.1415/360.0))
            hdu.header.update("DEC", src_img[0].header['DEC']+((scale*(h+2*gap)/4)/3600.0))
            hdu.header.add_history('Warp image created from %s' % file)
            hdulist.writeto(dst_path+"/Q4/"+os.path.basename(file).replace(".fits","_q4.fits"))
            hdulist.close(output_verify='ignore')
        except:
            raise
        

if __name__ == "__main__":
    import sys
    try:
        source_file_list = sys.argv[1]
        dst_path = sys.argv[2]
        gap = int(sys.argv[3])
    except IndexError:
        print "<source image path>, <destination image path> <gap>"
    else:
        filelist=[line.replace( "\n", "") for line in fileinput.input(source_file_list)]        
        makeSplit(filelist, dst_path, gap)
        
        
