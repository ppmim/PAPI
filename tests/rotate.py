import numpy as N
import pylab as P
import scipy.ndimage as ndi

# Create a test pattern
I = N.zeros((200,200,3),dtype=N.uint8)
I[0:50,:,0] = 255
I[:,0:50,1] = 255
I[150:200,:,2] = 255
I[:,150:200,1] = 128

# Apply rotation
I_rot_1 = ndi.rotate(I,30,axes=(0,1),order=1)
I_rot_2 = ndi.rotate(I,30,axes=(0,1),order=2)
I_rot_3 = ndi.rotate(I,30,axes=(0,1),order=2,prefilter=False)

# Display
P.subplot(141)
P.title('Original')
P.imshow(I)
P.subplot(142)
P.title('Rotation\n1st order spline')
P.imshow(I_rot_1)
P.subplot(143)
P.title('Rotation\n2nd order spline')
P.imshow(I_rot_2)
P.subplot(144)
P.title('Rotation\n2nd order spline\nNo prefilter')
P.imshow(I_rot_3)
P.show()
P.close()
