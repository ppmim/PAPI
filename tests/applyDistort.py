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
from optparse import OptionParser

def applyDistort(src_image, dst_path, dist_matrix , scale_factor=1.0, inverse=False):
    """ Apply (remove/add) to an image a distorion using a matrix distortion for a grid of points, so   
        firstly, we need to compute the full distorion matrix to all the pixels of the image,
        and then map the new coordinates to remove/add the distorion. Note how matrix interpolation
        is done (prediction-->real), because it caused to me some nightmares, but finally I caught it!;
        the point it that with that interpolation (predic-->real), we'll get a corrected image, because what
        we get is for each point in a generated-corrected image, what point of the distorted image should be,
        what is the input for the map_coordinates routine !
        So, if we wist to apply a artificial distortion, we should do the inverse interpolation, thus real--->predic
        
        
        Predicted : where points should be when no distortion
        Real      : where points are due to the optical distortion, then real means distorted
        
        # NOTA: me costo entender como hacer la interpolacion (prediction-->real),
        # pero es asi por como luego mapeamos con map_coordinates para corregir la distorsion
    """
    
    src_img=pyfits.open(src_image)
    src=src_img[0].data

    src_h, src_w = src.shape
         
    #read the distorion matrix from file     
    dm = readDistMat(dist_matrix)
         
    #SX=src_h
    #SY=src_w
    #generate full-image coordinates in mm, having into account a pixel is 18um
    #Xp,Yp = N.meshgrid(N.arange(0,SX)*0.018, N.arange(0,SY)*0.018) 
    
    
    t=1 # offset factor
    pixel_size=0.0168  # pixel size in mm, it means the sampling rate of the matrix interpolation
    offset=abs(min(dm[:,0]))
    
    x = dm[:,0] + t*offset # predicted-X (undistorted)
    y = dm[:,1] + t*offset  # predicted-Y (undistorted)
    zx = dm[:,2]*scale_factor + t*offset # real-X ; interpol values for x-axis
    zy = dm[:,3]*scale_factor + t*offset # real-Y ; interpol values for y-axis
    
    if not inverse: # correct distortion
        ipx = interpolate.interp2d(x,y,zx, kind='linear')
        ipy = interpolate.interp2d(x,y,zy, kind='linear')
    else:   # distort the image
        ipx = interpolate.interp2d(zx,zy,x, kind='linear')
        ipy = interpolate.interp2d(zx,zy,y, kind='linear')
        
    # another option would be use 'interpolate.griddata' instead of interp2d,
    # but I can't try it because it's needed  Scipy_version>=0.98.3
    
    #vx=N.arange(-38.267,38.367,0.036)
    #vy=N.arange(-38.267,38.367,0.036)
    vx=N.arange(0, offset*2, pixel_size)
    vy=N.arange(0, offset*2, pixel_size)
    
    Ix=ipx(vx,vy)/(pixel_size) # we divide by a sampling rate to have same points as pixels in the final image; 0.018 interpolated value with call
    Iy=ipy(vx,vy)/(pixel_size)
    
    # Apply distortion mapping the real(distorted) grid of points  previously interpolated
    new_image = scipy.ndimage.map_coordinates(src.reshape(src_h, src_w), \
                                              N.array([[Iy],[Ix]]), order=3) # ojo con la pos de los ejes !! 
    
    ## save results
    #remove old file
    if os.path.exists(dst_path): os.remove(dst_path)    
    # Save the new image in a FITS file
    hdu = pyfits.PrimaryHDU()
    hdu.header=src_img[0].header.copy()
    hdu.data=new_image[0,:,:]     
    hdulist = pyfits.HDUList([hdu])
    hdu.header.add_history('Warp image created from %s and distortion matrix %s' %(src_image, dist_matrix))
    hdulist.writeto(dst_path)
    hdulist.close(output_verify='ignore')

def readDistMat(dist_mat_path):
    """ Read distortion matrix from a file (given by ZEMAX) and return a MxN matrix"""
    a=N.loadtxt(dist_mat_path, skiprows=1, usecols=(5,6,7,8))
    
    return a
    
    
if __name__ == "__main__":
    import sys
    
    usage = "usage: %prog [options] arg1 arg2 ..."
    parser = OptionParser(usage)
    
    parser.add_option("-s", "--source_file",
                  action="store", dest="source_file", help="souce file (FITS) to apply matrix distortion")
                  
    parser.add_option("-o", "--output_file",
                  action="store", dest="output_file",
                  help="Output file generated after distortion matrix is applied")
    
    parser.add_option("-m", "--matrix_file",
                  action="store", dest="matrix_file", help="matrix file describing the distortion to apply to the source file")
    
    parser.add_option("-v", "--verbose",
                  action="store_true", dest="verbose", default=True,
                  help="verbose mode [default]")
    
    # optional values
    parser.add_option("-f", "--scale_factor", type="float", default=1.0,
                  action="store", dest="scale_factor", 
                  help="scale factor to multiply by the distortion matrix")
    
    parser.add_option("-i", "--inverse",
                  action="store_true", dest="inverse", default=False,
                  help="apply inverse distortion (correction?)")
    
                                
    (options, args) = parser.parse_args()
    
    
    if not options.source_file or not options.output_file or not options.matrix_file:
        parser.print_help()
        parser.error("incorrect number of arguments " )
    else:
        applyDistort(options.source_file, options.output_file, options.matrix_file, \
                     options.scale_factor, options.inverse)
        
