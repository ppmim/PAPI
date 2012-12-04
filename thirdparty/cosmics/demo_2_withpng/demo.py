# The following two lines are only needed as cosmic.py is not in this directory nor in the python path.
# They would not be required if you copy cosmics.py in this directory.
import sys
sys.path.append("../.") # The directory that contains cosmic.py


import cosmics, f2n

# Read the FITS :
(array, header) = cosmics.fromfits("euler.fits", verbose = True)

# We can of course crop the numpy array :
#array = array[100:500,210:525]

z1=170
z2=5000
upsample = 1

# Build the object :
c = cosmics.cosmicsimage(array, pssl = 0.0, gain=2.2, readnoise=9.5,
	sigclip = 5.0, sigfrac = 0.2, objlim = 3.0,
	satlevel = 50000.0, verbose = True)

# Run :
c.run(maxiter = 4, verbose = False)

# And now we use f2n.py to make several PNG images :


# The raw input array :
im = f2n.f2nimage(c.getrawarray(), verbose=False)
im.setzscale(z1, z2)
im.makepilimage("log")
im.upsample(upsample)
im.writetitle("Input image", colour = (0,255,0))
im.tonet("0_raw.png")


# The mask of saturated stars, upon the raw image :
im = f2n.f2nimage(c.getrawarray(), verbose=False)
im.setzscale(z1, z2)
im.makepilimage("log")
im.drawmask(c.getsatstars(), colour=(255, 0, 255))
im.upsample(upsample)
im.writetitle("Saturated stars", colour = (0,255,0))
im.tonet("1_satstars.png")

# We output a list of the positions of detected cosmic ray hits.
# This is made on purpose to be fed into f2n's drawstarslist :
labeldict = c.labelmask()
im = f2n.f2nimage(c.getrawarray(), verbose=False)
im.setzscale(z1, z2)
im.makepilimage("log")
im.drawmask(c.getsatstars(), colour=(255, 0, 255))
im.upsample(upsample)
im.drawstarslist(labeldict, colour=(255,0,0))
im.writetitle("Cosmic ray hits", colour = (0,255,0))
im.tonet("2_labels.png")

# One png with the precise mask in green and a wider version in blue :
im = f2n.f2nimage(c.getrawarray(), verbose=False)
im.setzscale(z1, z2)
im.makepilimage("log")
im.drawmask(c.getdilatedmask(size=5), colour=(0, 0, 255))
im.drawmask(c.getmask(), colour=(0, 255, 0))
im.upsample(upsample)
im.writetitle("Mask", colour = (0,255,0))
im.tonet("3_mask.png")

# And of course one png with the clean array :
im = f2n.f2nimage(c.getcleanarray(), verbose=False)
im.setzscale(z1, z2)
im.makepilimage("log")
im.upsample(upsample)
im.writetitle("Cleaned image", colour = (0,255,0))
im.tonet("4_clean.png")


