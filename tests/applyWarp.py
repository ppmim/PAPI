#! /usr/bin/env python

# Copyright (c) 2009 Jose M. Ibanez. All rights reserved.
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
# Created    : 14/01/2010    jmiguel@iaa.es
# Last update: 25/01/2010    jmiguel@iaa.es
#
################################################################################

# Import necessary modules
import pyfits
import numpy as np
import scipy.ndimage


class CRadialDistortion:
  
    def __init__(self, src_path, dst_path, xc, yc, k):
        """ Constructor """
        self.src_path=src_path
        self.dst_path=dst_path
        self.cx = xc
        self.cy = yc
        self.k  = k
        
        
        self.xscale = -1
        self.yscale = -1
        self.xshift = -1
        self.yshift = -1
        
        
    def getRadialX( self, x, y):
    
        x = (x*self.xscale+self.xshift)
        y = (y*self.yscale+self.yshift)
        res = x+((x-self.cx)*self.k*((x-self.cx)*(x-self.cx)+(y-self.cy)*(y-self.cy)))
        return res
             
    def getRadialY( self, x, y):
    
        x = (x*self.xscale+self.xshift)
        y = (y*self.yscale+self.yshift)
        res = y+((y-self.cy)*self.k*((x-self.cx)*(x-self.cx)+(y-self.cy)*(y-self.cy)))
        return res         
    
    def calc_shift(self,x1, x2, cx):
      
        x3 = x1+(x2-x1)*0.5
        res1 = x1+((x1-cx)*self.k*((x1-cx)*(x1-cx)))
        res3 = x3+((x3-cx)*self.k*((x3-cx)*(x3-cx)))

        thresh = 1
        if (res1>-thresh and res1 < thresh):
            return x1
        if (res3<0):
            return self.calc_shift(x3,x2,cx)
        else:
            return self.calc_shift(x1,x3,cx)
        
    def sampleImage( self, img, i,  j):
        """ Make image point interpolation (bi-linear, first order interpolation"""
        
        h, w = img.shape
        # avoid the borders of the image
        if(i<=0 or j<=0 or i>=(h-1) or j>=(w-1)):
            res=0
        else:
            p=i-int(i)
            q=j-int(j)
            res = (img[i+1,j]-img[i,j])*p + (img[i,j+1]-img[i,j])*q + (img[i+1,j+1]+img[i,j]-img[i+1,j]-img[i,j+1])*p*q + img[i,j]
        return res
    
  
    def makeDistort(self):
  
        src_img=pyfits.open(self.src_path)
        src=src_img[0].data
    
            
        centerX = self.cx
        centerY = self.cy
        height,width = src.shape
        img_res=np.empty((width, height),float)
    
        """self.xshift = self.calc_shift(0,centerX-1,centerX)
        newcenterX = width-centerX
        xshift_2 = self.calc_shift(0,newcenterX-1,newcenterX)
    
        self.yshift = self.calc_shift(0,centerY-1,centerY)
        newcenterY = height-centerY
        yshift_2 = self.calc_shift(0,newcenterY-1,newcenterY)
        #scale = (centerX-xshift)/centerX
        self.xscale = (width-self.xshift-xshift_2)/width
        self.yscale = (height-self.yshift-yshift_2)/height
        """
        self.xscale = 1
        self.yscale = 1
        self.xshift = 0
        self.yshift = 0
    
        #print "xshift: %f   yshift: %f   xscale: %f   yscale: %f  "%(self.xshift,self.yshift,self.xscale,self.yscale)
    
        for j in range(0,height):
            for i in range(0,width):
                x = self.getRadialX(float(i),float(j))
                y = self.getRadialY(float(i),float(j))
                res=self.sampleImage(src,y,x)
                #print "(%d,%d)=%f,%f ---> V= %f"%(j,i,y,x,res)
                img_res[j,i]=res
    

        hdu = pyfits.PrimaryHDU()
        hdu.scale('float32') # importat to set first data type
        hdu.data=img_res     
        hdulist = pyfits.HDUList([hdu])
        #hdu.header.update('OBJECT','MASTER_PIXEL_MASK')
        hdu.header.add_history('Warp image created from %s' % self.src_path)
        hdulist.writeto(self.dst_path)
        hdulist.close(output_verify='ignore')


    
    

if __name__ == "__main__":
    import sys
    try:
        src_path = sys.argv[1]
        dst_path = sys.argv[2]
    except IndexError:
        print "<source image path>, <destination image path>"
    else:
        dist = CRadialDistortion (src_path, dst_path, 1024,1024, -0.00000001)
        dist.makeDistort()
        
        