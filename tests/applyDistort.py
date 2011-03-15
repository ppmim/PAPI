#! /usr/bin/env python

# Copyright (c) 2010 Jose M. Ibanez. All rights reserved.
# Institute of Astrophysics of Andalusia, IAA-CSIC
#
# This file is part of PAPI (PANIC Pipeline)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
################################################################################
#
# PANICtool
#
# applyDistort.py
#
# Created    : 04/10/2010    jmiguel@iaa.es
#
# References: http://en.wikipedia.org/wiki/Bilinear_interpolation
################################################################################

# Import necessary modules
import pyfits
import numpy as N
import scipy.ndimage
import math
from scipy import interpolate
import os

def applyDistort(src_image, dst_path, dist_matrix ):
    """ Apply (remove/add) to an image a distorion using a matrix distortion for a grid of points, so   
        firstly, we need to compute the full distorion matrix to all the pixels of the image,
        and then map the new coordinates to remove/add the distorion. Note how matrix interpolation
        is done (prediction-->real), because it caused to me some nightmares, but finally I caught it!
    """
    
    src_img=pyfits.open(src_path)
    src=src_img[0].data

    src_h, src_w = src.shape
         
    #read the distorion matrix from file     
    dm = readDistMat(dist_matrix)
         
    SX=2048#4263
    SY=2048#4263
    #generate full-image coordinates in mm, having into account a pixel is 18um
    #Xp,Yp = N.meshgrid(N.arange(0,SX)*0.018, N.arange(0,SY)*0.018) 
    
    # ojo, me costo entender como hacer la interpolacion (prediction-->real),
    # pero es asi por como luego mapeamos con map_coordinates para corregir la distorsion
    x = dm[:,0] + 1*38.267
    y = dm[:,1] + 1*38.267
    zx = dm[:,2] + 1*38.767 # interpol values for x-axis
    zy = dm[:,3] + 1*38.767 # interpol values for y-axis
    ipx = interpolate.interp2d(x,y,zx, kind='linear')
    ipy = interpolate.interp2d(x,y,zy, kind='linear')
    # another option would be use 'interpolate.griddata' instead of interp2d,
    # but I can't try it because need  Scipy_version>=0.98.3
    
    #vx=N.arange(-38.267,38.367,0.036)
    #vy=N.arange(-38.267,38.367,0.036)
    vx=N.arange(0, 37.782*2,0.036)
    vy=N.arange(0, 37.296*2,0.036)
    
    Ix=ipx(vx,vy)/(0.036*1) # interpolated value with call
    Iy=ipy(vx,vy)/(0.036*1)
    
    new_image = scipy.ndimage.map_coordinates(src.reshape(src_h, src_w), \
                                              N.array([[Iy],[Ix]]), order=3) # ojo con la pos de los ejes !! 
    
    #remove old file
    if os.path.exists(dst_path): os.remove(dst_path)    
    # Save the new image in a FITS file
    hdu = pyfits.PrimaryHDU()
    hdu.header=src_img[0].header.copy()
    hdu.data=new_image[0,:,:]     
    hdulist = pyfits.HDUList([hdu])
    hdu.header.add_history('Warp image created from %s and distortion matrix %s' %(src_path, dist_matrix))
    hdulist.writeto(dst_path)
    hdulist.close(output_verify='ignore')

def readDistMat(dist_mat_path):
    """ Read distortion matrix from file an return a MxN matrix"""
    a=N.loadtxt(dist_mat_path, skiprows=1, usecols=(5,6,7,8))
    
    return a
    
    
if __name__ == "__main__":
    import sys
    try:
        src_path = sys.argv[1]
        dst_path = sys.argv[2]
        mat_path = sys.argv[3]
    except IndexError:
        print "<source image path>, <destination image path> <dist_mat>"
    else:
        applyDistort(src_path, dst_path, mat_path)
