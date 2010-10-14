#! /usr/bin/env python
#encoding:UTF-8

#euclideand warp

import pyfits
import numpy as N
import scipy.interpolate as interpol
from scipy import ndimage


# load the image 
myfits=pyfits.open("/tmp/image_orig.fits")

(h,w)=myfits[0].data.shape

# show input image
#TBD

#Transformation matrix (T)
s = 1 # scale
a = 5*N.pi/180.0 # rotation
tx = 0 # x-translation
ty = 0 # y-translation

T = N.matrix([[ s*N.cos(a), s*N.sin(a), tx], [ -s*N.sin(a), s*N.cos(a), ty], [0,0,1]])

# warp incoming corners to determine the size of the output image

cp = N.dot(T, N.array([[1,1,w,w], [1,h,1,h], [1,1,1,1]]))
Xpr = N.arange(N.min(cp[0,:]), N.max(cp[0,:]))
Ypr = N.arange(N.min(cp[1,:]), N.max(cp[1,:])) 

Xp,Yp = N.meshgrid(Xpr, Ypr)
wp,hp = Xp.shape
print "WP=%d,HP=%d"%(wp,hp)
n = wp*hp
#X = T.transpose() \ [ Xp  Yp N.ones((n,1))].transpose()
M=N.array([ Xp.flatten(1),  Yp.flatten(1), N.ones((n,1))])    
#X = T.inv*M' 
X = N.linalg.solve(T, M )
print "T.shape=",T.shape
print "M.shape=",M.shape

xI = X[0,:].reshape((wp,hp))
yI = X[1,:].reshape((wp,hp))
print "xI size=",xI.size
print "yI size=",yI.size
    
#Ip= interpol.interp2d( xI.flatten(1), yI.flatten(1), myfits[0].data.flatten(1),'linear')
#re-sample pixel values with 'linear' interpolation
new_image = ndimage.map_coordinates(myfits[0].data.reshape((h,w)), N.array([[xI],[yI]]))

#write new image in FITS file
hdu = pyfits.PrimaryHDU()
hdu.scale('float32')
hdu.data=new_image[0]
hdulist = pyfits.HDUList([hdu])
hdulist.writeto("/tmp/rota.fits")
hdulist.close(output_verify='ignore')    






