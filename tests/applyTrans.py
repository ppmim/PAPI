#!/usr/bin/env python
################################################################################
#
# PANICtool
#
# applyTrans.py
#
# Created    : 19/01/2010    jmiguel@iaa.es
# Last update: 26/01/2010
#
################################################################################

# Import necessary modules
import pyfits
import numpy as np
import scipy.ndimage as nd
import math

def applyTransformation(src_path, dst_path, rot=0, xc=0, yc=0, scale=1, shear=0 ):
    """ Apply an Affine transformation to an image (a rotation, a translation, a scale, or a shear)"""
    
    src_img=pyfits.open(src_path)
    src=src_img[0].data

    src_h, src_w = src.shape
    #dest = np.empty(src_w, src_h)
        
          
    #def shift_func(output_coordinates):
    #    return (output_coordinates[0]-10.5, output_coordinates[1]-10.5)        
    
    #dst=nd.geometric_transform(src, shift_func)
    cos_r=math.cos(2*3.141592*rot/360) 
    sin_r=math.sin(2*3.141592*rot/360)
    # anticlockwise rotation arround point xc,yc
    H=np.array([[cos_r, -sin_r, xc-cos_r*xc+sin_r*yc],[sin_r, cos_r, yc-sin_r*xc-cos_r*yc],[0,0,1]])
    #H=np.array([[1, 0, -512],[0, 1, -512],[0,0,1]])
    dst=nd.affine_transform(src, H[:2,:2],(H[0,2],H[1,2]))
        
    
    hdu = pyfits.PrimaryHDU()
    hdu.header=src_img[0].header.copy()
    hdu.data=dst     
    hdulist = pyfits.HDUList([hdu])
    #hdu.header.update('OBJECT','MASTER_PIXEL_MASK')
    hdu.header.add_history('Warp image created from %s' % src_path)
    hdulist.writeto(dst_path)
    hdulist.close(output_verify='ignore')


if __name__ == "__main__":
    import sys
    try:
        src_path = sys.argv[1]
        dst_path = sys.argv[2]
        rot = sys.argv[3]
    except IndexError:
        print "<source image path>, <destination image path> <rot>"
    else:
        applyTransformation(src_path, dst_path, int(rot), 1024, 1024)
