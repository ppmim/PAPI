#!/usr/bin/env python

# A collection of tasks I find useful for reducing echelle spectra.  These
# can all be called within pyraf after an 'from myiraf import *'.  

from pyraf import iraf
import string
import re
import math
import os
import numarray as num
from ReadSex import readsex
import stats
import pygplot

# a file patter, omitting the trailing [#]...
file_pat = re.compile(r'(^.*?)\[.*')

def display_triplet(root):
   # Display the original, template and difference images in 3 frames
   SN = root+"SN.fits"
   diff = root+"diff.fits"
   temp = root+"temp.fits"
   iraf.display(SN, 1)
   iraf.display(temp,2)
   iraf.display(diff,3)

def get_image_size(images):
   result = iraf.imheader(images=images, longheader=0, Stdout=1)
   if result[0] == '':
      return None
   array = {}
   for item in result:
      res = re.search(r'(.*)\[([0-9]+),([0-9]+)\]', item)
      array[res.group(1)] = [int(res.group(2)),int(res.group(3))]

   if len(array.keys()) == 1:
      return array[array.keys()[0]]
   else:
      return array

def get_mask(image):
   nmasks = get_header(image, 'NMASKS', float)
   im_masks = []
   if nmasks is not None:
      for i in range(nmasks):
         im_masks.append(map(int,string.split(\
               get_header(image, "MASK%d" % i,str))))
   else:
      nmasks = 0
   return (im_masks)


def sex(image, output, thresh=2.0, fwhm=4, gain=None, zmag=None,
        sexdir='/Users/burns/CSP/template_subtraction/sex', scale=0.125,
        sigma_map=None):
   '''Construct a sextractor command and run it.'''
   if sexdir[-1] != '/':  sexdir += '/'
   com = ["sex", image, "-c "+sexdir+"default.sex",
          "-PARAMETERS_NAME "+sexdir+"default.param",
          "-DETECT_MINAREA 5", "-DETECT_THRESH " + str(thresh),
          "-ANALYSIS_THRESH " + str(thresh), "-PIXEL_SCALE %f" % (scale),
          "-SEEING_FWHM "+str(fwhm), "-CATALOG_NAME "+output,
          "-FILTER_NAME "+sexdir+"gauss_3.0_5x5.conv",
          "-STARNNW_NAME "+sexdir+"default.nnw",
          "-SATUR_LEVEL "+str(30000), "-VERBOSE_TYPE QUIET"]
   if gain is not None:
      com += ["-GAIN",str(gain)]
   if zmag is not None:
      com += ["-MAG_ZEROPOINT", str(zmag)]
   if sigma_map is not None:
      com += ["-WEIGHT_IMAGE",sigma_map,"-WEIGHT_TYPE MAP_WEIGHT"]
   com = string.join(com)
   print com
   res = os.system(com)
   return res


def rectify_image(image, reference, output, refout=None, fitgeometry='general', 
                   function='polynomial',
                   xxorder=3, yyorder=3, xyorder=3, yxorder=3, clobber=1,
                   objects=None, interactive=1, maxiter=5, reject=3, xrange=None,
                   yrange=None, thresh=2.0, scale=0.125, fwhm=4, image_sigma=None,
                   reference_sigma=None, debug=0, aoffset=0):
   '''This task will take an image and a reference image and apply geomap/
   geotran to aligned them.  image is registered and output to output.
   You can specify the fit geometry (shift, xyscale, rotate, rscale,
   rxyscale, general).  You can also choose the geenral function
   (legendre, chebyshev, polynomial).  You can specfy the order of all
   the terms (xxorder, yyorder, xyorder, yxorder).  If you choose
   clobber=1 (the default), new images overwrite older ones.  Lastly, if
   you want to use a pre-existing list of objects, specify them as the
   objects keyword.'''
   
   # First, if we already have objects to use, save time by not finding them
   #  again!

   if objects is None:
      sex(image, output='image.cat', thresh=thresh, scale=scale, fwhm=fwhm,
            sigma_map=image_sigma)
      sex(reference, output='reference.cat', thresh=thresh, scale=scale,
            fwhm=fwhm, sigma_map=reference_sigma)

      im_masks1 = get_mask(image)
      im_masks2 = get_mask(reference)
      if aoffset == 0:
         # Check the ROTANG header field to see if we need to give a rotation
         # offset to the object matcher.
         rotang1 = get_header(reference, "ROTANG", float)
         rotang2 = get_header(image, "ROTANG", float)
         if rotang1 is not None and rotang2 is not None:
            if abs(rotang1-rotang2) > 0.1:
               aoffset = rotang1-rotang2
 
      # Run the matchum routine to get corresponding objects:
      (x2,y2,x1,y1,o1,o2) = matchum('image.cat', 'reference.cat', 
            im_masks1=im_masks1, im_masks2=im_masks2, xrange=xrange,
            yrange=yrange, debug=debug, aoffset=aoffset)

 
      f = open('geomap.objects', 'w')
      for i in range(len(x1)):
         f.write('%.2f %.2f %.2f %.2f\n' % (x1[i],y1[i],x2[i],y2[i]))
      f.close()
      objects = "geomap.objects"
   else:
      data = pygplot.columns(objects)
      x1 = data[0]
      y1 = data[1]
      x2 = data[2]
      y2 = data[3]

   xshift = int(stats.bwt(x2-x1)[0])
   yshift = int(stats.bwt(y2-y1)[0])

   # We now have some selected objects to use with geomap.  Let's figure out
   #  the area of validity for the data using ths shifts computed above.
   (r_xmax,r_ymax) = get_image_size(reference)
   if xshift < 0:
      xmin = -xshift
      xmax = r_xmax
   else:
      xmin = 1
      xmax = r_xmax - xshift
   if yshift < 0:
      ymin = -yshift
      ymax = r_ymax
   else:
      ymin = 1
      ymax = r_ymax - yshift
   #xmin = min(x1);  xmax = max(x1)
   #ymin = min(y1);  ymax = max(y1)
   print xshift,yshift,xmin,xmax,ymin,ymax
   
   iraf.images()
   if os.path.isfile('geomap.db'):  os.unlink('geomap.db')

   # We may not have enough points to do a general fit, so we do them in order of
   #   complexity.  Hopefully, we have at least ONE point (for the shift)
   fitgeometries = [fitgeometry, 'rxyscale', 'rotate', 'shift']

   success = 0
   for fitgeo in fitgeometries:
      print "using fitgeometry = ", fitgeo
      res = iraf.geomap(input=objects,database='geomap.db', xmin=1,
               ymin=1, xmax=r_xmax, ymax=r_ymax, fitgeometry=fitgeo,
               function=function, xxorder=xxorder, yyorder=yyorder, 
               xyorder=xyorder, yxorder=yxorder, results='', 
               transforms=objects, interactive=interactive, 
               maxiter=maxiter, reject=reject, Stdout=1)
      # assume all went well
      success = 1
      for line in res:
         print line
         if re.search('Too few points',line):
            success = 0
      if success:
         break

   if clobber and os.path.isfile(output):  iraf.imdel(output)
   iraf.geotran(input=image, output='temp.rectify', database='geomap.db', 
         transforms=objects, geometry='geometric', interpolant="linear")
   section = "[%d:%d,%d:%d]" % (xmin, xmax, ymin, ymax)
   print 'Valid section of data:', section
   iraf.imcopy(input='temp.rectify',output=output)
   iraf.imdel('temp.rectify')
   return(xmin,xmax,ymin,ymax)

def mygeotran(input, output, database, transforms, geometry='geometric', 
      interpolant="linear"):

   # First, we get the shifts from the database
   f = open(database, 'r')
   for line in f.readlines():
      if re.search('xshift', line):
         xshift = float(string.split(line)[1])
      elif re.search('yshift', line):
         yshift = float(string.split(line)[1])
   f.close()
   f = open('offsets.dat', 'w')
   print >>f, 0,0
   print >>f, xshift+1,yshift+1
   f.close()

   # Now we make a blank holder for the mapped image
   nuke('placeholder.fits')
   iraf.imcombine(input+','+input, 'placeholder.fits', offsets='offsets.dat')
   iraf.imarith('placeholder.fits', '*', 0.0, 'placeholder.fits')

   xsize,ysize = get_image_size(input)
   # next, copy the input into this place holder
   iraf.imcopy(input+'[*,*]', 'placeholder.fits[1:%d,1:%d]' % (xsize,ysize))

   # now run geogran on this padded image
   iraf.geotran(input='placeholder.fits', output=output, database=database,
         transforms=transforms, geometry=geometry, interpolant=interpolant)
#   nuke('placeholder.fits')

def combine_images(images, output, reference=None, method='average', zscale=1,
      fitgeometry='general', function='polynomial', xxorder=3, yyorder=3, 
      xyorder=3, yxorder=3, objects=None, trim=0, interactive=1, maxiter=5,
      reject=3, sigma=None, xrange=None, yrange=None, rectify=1, thresh=2.0, scale=0.125,
      fwhm=4.0, debug=0, aoffset=0, creject=None):
   '''This Routine will take a list of images and use rectify_image to rectify
   each with respect to referece.  It will then use imcombine to combine the
   images into one and output to "output".  Default method is 'average', but
   you can use any available to imcombine.  The rest of the arguments are sent
   to rectify_images.  If zscale=1, all the images are scaled to a common
   zmag (30).  Also, if a file with _sigma.fits is found for all input files,
   we combine them too.  If you specify a sigma, that file will hold the
   standard deviations from the files.'''
   names = get_files(images)
   if reference is None:
      reference = names[0]
      names = names[1:]

   pclip = -0.5
   lsigma = 3.
   hsigma = 3.
   nkeep = -2
   if creject:
      creject = "pclip"
   else:
      creject = "none"

   # Check to see if we have _sigma files for each input file
   if os.path.isfile(reference.replace('.fits','_sigma.fits')):
      do_sigmas = 1
   else:
      do_sigmas = 0

   sigmas = []
   for name in names:
      if os.path.isfile(name.replace('.fits','_sigma.fits')):
         sigmas.append(name.replace('.fits','_sigma.fits'))
      else:
         do_sigmas = 0
         break

   # Put the reference image first, to ensure that all positional header keywords are
   #  valid.
   if do_sigmas:
      ref_sig = reference.replace('.fits','_sigma.fits')
   update = 0
   if zscale:
      zmag = get_header(reference, "ZMAG", float)
      if zmag is not None and zmag != 26:
         update = 1
         zmagfac = num.power(10, (zmag-30)/2.5)
         if zmagfac != 1:
             print 'zmagfac = ',zmagfac
             iraf.imarith(reference, '/', zmagfac, reference)
             iraf.hedit(reference, "ZMAG", 30.0, add=1, verify=0)
             if do_sigmas:
                iraf.imarith(ref_sig, '/', zmagfac, ref_sig)
                iraf.hedit(ref_sig, "ZMAG", 30.0, add=1, verify=0)

   print "Using %s as reference" % (reference)

   if rectify:
      nuke(reference.replace('.fits','_rec.fits'))
      iraf.imcopy(reference, reference.replace('.fits','_rec.fits'))
      if do_sigmas:
         nuke(reference.replace('.fits','_rec_sigma.fits'))
         iraf.imcopy(reference.replace('.fits','_sigma.fits'),
                     reference.replace('.fits','_rec_sigma.fits'))
      temps = [reference.replace('.fits','_rec.fits')]
      if do_sigmas:  sig_temps = [reference.replace('.fits','_rec_sigma.fits')]
   else:
      temps = [reference]
      if do_sigmas:  sig_temps = [ref_sig]
   
   i = 0
   for name in names:
      print name
      if rectify:
         temp = name.replace('.fits','_rec.fits')
         nuke(temp)
         if do_sigmas:
            (x0,x1,y0,y1) = rectify_image(name, reference, output=temp, 
                 fitgeometry=fitgeometry, 
                 function=function, xxorder=xxorder, yyorder=yyorder,
                 xyorder=xyorder, yxorder=yxorder, objects=objects,
                 interactive=interactive, maxiter=maxiter, reject=reject, 
                 xrange=xrange,
                 yrange=yrange, thresh=thresh, scale=scale, fwhm=fwhm,
                 image_sigma=sigmas[i], reference_sigma=ref_sig, debug=debug,
                 aoffset=aoffset)
         else:
            (x0,x1,y0,y1) = rectify_image(name, reference, output=temp, 
                 fitgeometry=fitgeometry, 
                 function=function, xxorder=xxorder, yyorder=yyorder,
                 xyorder=xyorder, yxorder=yxorder, objects=objects,
                 interactive=interactive, maxiter=maxiter, reject=reject, 
                 xrange=xrange,
                 yrange=yrange, thresh=thresh, scale=scale, fwhm=fwhm, debug=debug,
                 aoffset=aoffset)
      else:
         temp = name
      bpm = get_header(name, 'BPM', str)
      if bpm and rectify:
         bpm_fits = bpm.replace('.pl','.fits')
         temp_bpm = bpm.replace('.pl','_rec.fits')
         temp_bpm_pl = bpm.replace('.pl','_rec.pl')
         nuke(bpm_fits);  nuke(temp_bpm);  nuke(temp_bpm_pl)
         iraf.imcopy(bpm, bpm_fits)
         iraf.imarith(bpm_fits, "*", 10.0, bpm_fits)
         if objects is None:  sigmaobjects = 'geomap.objects'
         iraf.geotran(input=bpm_fits, output=temp_bpm, database='geomap.db', 
            transforms=sigmaobjects, geometry='geometric', 
            interpolant="linear")
         iraf.imreplace(temp_bpm, lower=0.02, upper="INDEF", value=1.)
         iraf.imcopy(temp_bpm, temp_bpm_pl)
         iraf.hedit(temp, 'BPM', temp_bpm_pl, update=1, verify=0)
      if do_sigmas:
         if rectify:
            temp_sig = sigmas[i].replace('_sigma.fits','_rec_sigma.fits')
            nuke(temp_sig)
            if objects is None:  
               sigmaobjects = 'geomap.objects'
            else:
               sigmaobjects = objects

            iraf.geotran(input=sigmas[i], output=temp_sig, database='geomap.db', 
               transforms=sigmaobjects, geometry='geometric', 
               interpolant="linear")
         else:
            temp_sig = sigmas[i]
      if rectify:
         if i == 0:
            xmin = x0;  ymin = y0;   xmax = x1;  ymax = y1
         else:
            xmin = max(xmin,x0);  ymin = max(ymin,y0)
            xmax = min(xmax,x1);  ymax = min(ymax,y1)
      if zscale:
         zmag = get_header(name, "ZMAG", float)
         if zmag is not None and zmag != 26:
            zmagfac = num.power(10, (zmag-30)/2.5)
            if zmagfac != 1:
               print 'zmagfac = ',zmagfac
               iraf.imarith(temp, '/', zmagfac, temp)
               iraf.hedit(temp, 'ZMAG', 30.0, add=1, verify=0)
               if do_sigmas:
                  iraf.imarith(temp_sig, '/', zmagfac, temp_sig)
                  iraf.hedit(temp_sig, 'ZMAG', 30.0, add=1, verify=0)

      temps.append(temp)
      if do_sigmas:  sig_temps.append(temp_sig)
      i = i + 1
   Nimages = len(temps)
   temps = string.join(temps, ',')
   print "Combining ", temps
   nuke(output)
   if sigma is not None:
      iraf.imcombine(temps, output, combine=method, sigma=sigma, 
            masktyp='badvalue', maskval=1., reject=creject,
            lsigma=lsigma, hsigma=hsigma, nkeep=nkeep)
   else:
      iraf.imcombine(temps, output, combine=method, masktyp='badvalue', 
            maskval=1., reject=creject, lsigma=lsigma, hsigma=hsigma, nkeep=nkeep)

   if do_sigmas:
      Ns = []
      for this in sig_temps:
         N = this.replace('.fits','_N.fits')
         nuke(N)
         iraf.imcopy(this, N)
         iraf.imreplace(N, lower=0.01, upper="INDEF", value=1)
         Ns.append(N)
      sig_temps = string.join(sig_temps, ',')
      Ns = string.join(Ns,',')
      iraf.imfunction(sig_temps, sig_temps, 'square')
      file = output.replace('.fits','_sigma.fits')
      nuke(file)
      iraf.imcombine(sig_temps, file, combine='sum')
      if method == 'average':
         nuke("N.fits")
         iraf.imcombine(Ns, 'N.fits', combine='sum')
         iraf.imfunction('N.fits','N.fits','square')
         iraf.imarith(file, '/', 'N.fits', file)
      iraf.imfunction(file, file, 'sqrt')
      iraf.imfunction(sig_temps, sig_temps, 'sqrt')
   if update:
      iraf.hedit(output, "ZMAG", 30.0, verify=0, add=1)
   #nuke('temp?.fits')

   # Now, let's clean up the gain and rdnoise headers, as they don't seem to
   #  be handled in the PANIC pipeline.  NOTE:  we really need a noise map, but
   #  for now, we'll deal with just global values.
   gain = get_header(output, "gain", float)
   rdnoise = get_header(output, "rdnoise", float)
   if gain is None:
      # have to compute it from scratch
      egain = get_header(output, "egain", float)
      enoise = get_header(output, "enoise", float)
      nloop = get_header(output, "nloops", int)
      i = 1
      key = ""
      # This seems to be the only safe way to compute this.  Can't trust
      #   NDITHERS
      while key is not None:
         key = get_header(output, "imcmb%03d" % (i), str)
         i = i + 1
      ndither = i - 1
      print "N, ndither, nloop, egain, enoise = ", Nimages, ndither, nloop, egain, enoise
      gain = egain*Nimages*nloop*ndither
      rdnoise = enoise*num.sqrt(Nimages*nloop*ndither)
   else:
      gain = gain*Nimages
      rdnoise = rdnoise*num.sqrt(Nimages)
   iraf.hedit(output, "gain", gain, verify=0, add=1)
   iraf.hedit(output, "rdnoise", rdnoise, verify=0, add=1)

   # If requested, trim to the common regions of the combined images
   if rectify and trim and Nimages > 1:
      iraf.imcopy(output+"[%d:%d,%d:%d]" % (xmin,xmax,ymin,ymax), 
            output.replace('.fits','_trim.fits'))

   
def matchum(file1, file2, tol=10, perr=4, aerr=1.0, nmax = 40, 
      im_masks1=[], im_masks2=[], debug=0, domags=0, xrange=None,
      yrange=None, sigma=4, aoffset=0):
   '''Take the output of two sextractor runs and match up the objects with
   each other (find out which objects in the first file match up with
   objects in the second file.  The routine considers a 'match' to be any 
   two objects that are closer than tol pixels (after applying the shift).  
   Returns a 6-tuple:  (x1,y1,x2,y2,o1,o2).  o1 and o2 are the ojbects
   numbers such that o1[i] in file 1 corresponds to o2[i] in file 2.'''
   NA = num.NewAxis

   sexdata1 = readsex(file1)
   sexdata2 = readsex(file2)

   # Use the readsex data to get arrays of the (x,y) positions
   x1 = num.asarray(sexdata1[0]['X_IMAGE'])
   y1 = num.asarray(sexdata1[0]['Y_IMAGE'])
   x2 = num.asarray(sexdata2[0]['X_IMAGE'])
   y2 = num.asarray(sexdata2[0]['Y_IMAGE'])
   m1 = num.asarray(sexdata1[0]['MAG_BEST'])
   m2 = num.asarray(sexdata2[0]['MAG_BEST'])
   o1 = num.asarray(sexdata1[0]['NUMBER'])
   o2 = num.asarray(sexdata2[0]['NUMBER'])
   f1 = num.asarray(sexdata1[0]['FLAGS'])
   f2 = num.asarray(sexdata2[0]['FLAGS'])

   # First, make a cut on the flags:
   gids = num.where(f1 < 4)
   x1 = x1[gids];  y1 = y1[gids];  m1 = m1[gids]; o1 = o1[gids]
   gids = num.where(f2 < 4)
   x2 = x2[gids];  y2 = y2[gids];  m2 = m2[gids]; o2 = o2[gids]

   # next, if there is a range to use:
   if xrange is not None and yrange is not None:
      cond = num.greater(x1, xrange[0])*num.less(x1,xrange[1])*\
            num.greater(y1, yrange[0])*num.less(y1,yrange[1])
      gids = num.where(cond)
      x1 = x1[gids];  y1 = y1[gids];  m1 = m1[gids]; o1 = o1[gids]
      cond = num.greater(x2, xrange[0])*num.less(x2,xrange[1])*\
            num.greater(y2, yrange[0])*num.less(y2,yrange[1])
      gids = num.where(cond)
      x2 = x2[gids];  y2 = y2[gids];  m2 = m2[gids]; o2 = o2[gids]

   # Use the user masks
   for m in im_masks1:
      print "applying mask (%d,%d,%d,%d)" % tuple(m)
      condx = num.less(x1, m[0]) + num.greater(x1, m[1])
      condy = num.less(y1, m[2]) + num.greater(y1, m[3])
      gids = num.where(condx + condy)
      x1 = x1[gids];  y1 = y1[gids];  m1 = m1[gids]; o1 = o1[gids]

   for m in im_masks2:
      print "applying mask (%d,%d,%d,%d)" % tuple(m)
      condx = num.less(x2, m[0]) + num.greater(x2, m[1])
      condy = num.less(y2, m[2]) + num.greater(y2, m[3])
      gids = num.where(condx + condy)
      x2 = x2[gids];  y2 = y2[gids];  m2 = m2[gids];  o2 = o2[gids]

   if nmax:
      if len(x1) > nmax:
         ids = num.argsort(m1)[0:nmax]
         x1 = x1[ids]
         y1 = y1[ids]
         m1 = m1[ids]
         o1 = o1[ids]
      if len(x2) > nmax:
         ids = num.argsort(m2)[0:nmax]
         x2 = x2[ids]
         y2 = y2[ids]
         m2 = m2[ids]
         o2 = o2[ids]
   if debug:
      print "objects in frame 1:"
      print o1
      print "objects in frame 2:"
      print o2
      mp = pygplot.MPlot(2,1,device='/XWIN')
      p = pygplot.Plot()
      p.point(x1,y1)
      [p.label(x1[i],y1[i],"%d" % o1[i]) for i in range(len(x1))]
      mp.add(p)
      p = pygplot.Plot()
      p.point(x2,y2)
      [p.label(x2[i],y2[i],"%d" % o2[i]) for i in range(len(x2))]
      mp.add(p)
      mp.plot()
      mp.close()

   # Now, we make 2-D arrays of all the differences in x and y between each pair
   #  of objects.  e.g., dx1[n,m] is the delta-x between object n and m in file 1 and
   #  dy2[n,m] is the y-distance between object n and m in file 2.
   dx1 = x1[NA,:] - x1[:,NA];   dx2 = x2[NA,:] - x2[:,NA]
   dy1 = y1[NA,:] - y1[:,NA];   dy2 = y2[NA,:] - y2[:,NA]
   # Same, but with angles
   da1 = num.arctan2(dy1,dx1)*180/num.pi
   da2 = num.arctan2(dy2,dx2)*180/num.pi
   # Same, but with absolute distances
   ds1 = num.sqrt(num.power(dx1,2) + num.power(dy1,2))
   ds2 = num.sqrt(num.power(dx2,2) + num.power(dy2,2))

   # Here's the real magic:  this is a matrix of matrices (4-D).  Consider 4 objects:
   #  objects i and j in file 1 and objects m and n in file 2.  dx[i,j,m,n] is the
   #  difference between delta-xs for objects i,j in file 1 and m,n in file 2.  If object
   #  i corresponds to object m and object j corresponds to object n, this should be a small
   #  number, irregardless of an overall shift in coordinate systems between file 1 and 2.
   dx = dx1[::,::,NA,NA] - dx2[NA,NA,::,::]
   dy = dy1[::,::,NA,NA] - dy2[NA,NA,::,::]
   da = da1[::,::,NA,NA] - da2[NA,NA,::,::] + aoffset
   ds = ds1[::,::,NA,NA] - ds2[NA,NA,::,::]
   # pick out close pairs.
   #use = num.less(dy,perr)*num.less(dx,perr)*num.less(num.abs(da),aerr)
   use = num.less(ds,perr)*num.less(num.abs(da), aerr)
   use = use.astype(num.Int32)

   #use = num.less(num.abs(da),perr)
   suse = num.add.reduce(num.add.reduce(use,3),1)
   print suse[0]

   guse = num.greater(suse,suse.flat.max()/2)
   i = [j for j in range(x1.shape[0]) if num.sum(guse[j])]
   m = [num.argmax(guse[j]) for j in range(x1.shape[0]) if num.sum(guse[j])]
   xx0,yy0,oo0,mm0 = num.take([x1,y1,o1,m1],i,1)
   xx1,yy1,oo1,mm1 = num.take([x2,y2,o2,m2],m,1)
   if debug:
      mp = pygplot.MPlot(2,1,device='/XWIN')
      p = pygplot.Plot()
      p.point(xx0,yy0)
      [p.label(xx0[i],yy0[i],"%d" % oo0[i]) for i in range(len(xx0))]
      mp.add(p)
      p = pygplot.Plot()
      p.point(xx1,yy1)
      [p.label(xx1[i],yy1[i],"%d" % oo1[i]) for i in range(len(xx1))]
      mp.add(p)
      mp.plot()
      mp.close()
   xshift,xscat = stats.bwt(xx0-xx1)
   xscat = max([1.0,xscat])
   yshift,yscat = stats.bwt(yy0-yy1)
   yscat = max([1.0,yscat])
   mshift,mscat = stats.bwt(mm0 - mm1)
   print "xscat = ", xscat
   print "yscat = ", yscat
   print "xshift = ", xshift
   print "yshift = ", yshift
   print "mshift = ", mshift
   print "mscat = ", mscat
   keep = num.less(num.abs(xx0-xx1-xshift),sigma*xscat)*\
         num.less(num.abs(yy0-yy1-yshift),sigma*yscat)
   # This is a list of x,y,object# in each file.
   xx0,yy0,oo0,xx1,yy1,oo1 = num.compress( keep, [xx0,yy0,oo0,xx1,yy1,oo1], 1)

   if debug:
      print file1, oo0
      print file2, oo1
      mp = pygplot.MPlot(2,1, device='temp.ps/CPS')
      p1 = pygplot.Plot()
      p1.point(xx0,yy0, symbol=25, color='red')
      for i in range(len(xx0)):
         p1.label(xx0[i],yy0[i], " %d" % oo0[i], color='red')
      mp.add(p1)
      p2 = pygplot.Plot()
      p2.point(xx1,yy1, symbol=25, color='green')
      for i in range(len(xx1)):
         p2.label(xx1[i],yy1[i], " %d" % oo1[i], color='green')
      mp.add(p2)
      mp.plot()
      mp.close()


   if domags:
      return(xx0,yy0,mm0,xx1,yy1,mm1,mshift,mscat,oo0,oo1)
   else:
      return(xx0,yy0,xx1,yy1,oo0,oo1)
   
def rename_files(files, key, keypattern, oldpattern, newpattern, exact=0):
   '''Quick little function that lets you search a list of files for
      a particular key (e.g., OBJECT) and if it matches the keypattern
      (e.g., 'flat'), it applies imrename using the provided oldpattern
      and newpattern to do the renaming (e.g. 'a' and 'comp').  If you 
      want an exact match for the key pattern, specify exact=1'''
   list = get_header(files, key, str)
   if exact:  keypattern = r"\A" + keypattern + r"\Z"
   for name in list:
      if re.search(keypattern, list[name]):
         iraf.imrename(oldnames=name, newnames=re.sub(oldpattern,newpattern,
                       name))
         
def nuke(fileparm):
   '''Delete a list of images/files, but only try it if they exist. USE WITH CAUTION.
   With FITS files, this is a bit over-cautious, but with imh files, you reallly need
   to use imdel to get rid of the .pix files properly.'''
   files = get_files(fileparm)
   if files is None: return
   for file in files:
      if os.path.isfile(file):
         if re.search(r'(fits|imh)$', file):
            iraf.imdel(file)
         else:
            os.unlink(file)

def get_files(fileparm):
   '''In iraf cl, sometimes a pattern or a file list is given instead of a 
   proper name.  Python doesn't handle that very well, so here's a cludge to
   get it to work using imhead.'''

   result = iraf.imhead(images=fileparm, longheader=0, Stdout=1)
   if len(result) == 0:
      return None
   if result[0] == '' or result[0] == 'no images found':
      return None

   out = []
   for res in result:
      out.append(file_pat.search(res).group(1))

   return out

def label_objects(images, xfield, yfield):
   '''Given a list of images, allow the user to pick an object in the
   display and have the x and y coordinates saved as header keywords
   xfield and yfield.'''
   files = get_files(images)
   print "Please use the 'a' key to select objects"
   for file in files:
      iraf.display(file,1)
      res = iraf.imexamine(file, Stdout=1)
      #res = iraf.rimcursor(file, Stdout=1)
      x,y = map(float,(string.split(res[2])[0:2]))

      iraf.hedit(images=file, fields=xfield, value=x, add=1, verify=0)
      iraf.hedit(images=file, fields=yfield, value=y, add=1, verify=0)


def get_header(images,field,type):
   '''Get a header keyword (field) from the images (images) and apply
      the given type conversion (type).  Returns an associated array
      with the image name as key, unless there's only one image, in which
      case, just the value of the key is returned.'''
   result = iraf.hedit(images=images, fields=field, value='.', Stdout=1)
   if len(result) == 0:  return None
   if result[0] == '':  return None
   array = {}
   for res in result:
      res = re.sub('"', '', res)
      value = type(string.split(res, '=')[1])
      name = string.split(res, ',')[0]
      array[name] = value
   if len(array.keys()) == 1:
      return array[array.keys()[0]]
   else:
      return array

def closest_time(images, reference, keyword='JD', printit=0):
   '''This function takes a list of images (in the iraf sense) and determines
   which one is closest in time to the reference (good to figure out 
   which comparison spectra to use.'''

   todo = get_files(reference)
   if todo is None:  
      raise IOError, "Could not access any images:  %s" % images

   to_return = {}
   for this_ref in todo:
      result = get_header(images, keyword, float)
      list = result.keys()
      ref = get_header(this_ref, keyword, float)
      ref = float(ref[this_ref])
 
      diff = abs(ref - result[list[0]])
      best = list[0]
      for obj in list[1:]:
         if abs(ref - result[obj]) < diff:
            diff = abs(ref - result[obj])
            best = obj
      to_return[this_ref] = best

   if printit:
      list = to_return.keys()
      list.sort()
      for item in list:
         print item,': ',to_return[item]
      return
         

   if len(to_return.keys()) == 1:
      return to_return[to_return.keys()[0]]
   else:
      return to_return

def sort_dict(dict, cmp=lambda x, y: x > y):
   '''A function that takes a dictionary and sorts the values using the 
      proviced cmp function.  It then returns the associated sorted keys.
      '''
   keys = dict.keys()
   indices = range(len(keys))
   indices.reverse()
   for j in indices:
      for i in range(j):
         if cmp(dict[keys[i]], dict[keys[i+1]]): 
            (keys[i], keys[i+1]) = (keys[i+1], keys[i])
   return keys
               
def closest_time_pos(images, reference, time_key='JD', ra_key='RA', 
      dec_key='DEC', printit=0):
   '''This function takes a list of images (in the iraf sense) and determines
   which one is closest in time and position to each reference (good to 
   figure out which comparison spectra to use.  It returns a dictionary
   where the key is the reference image and the value is the image closest
   in time/position to it.  This can be fed directly into extractComps.'''


   # This is a bias multiplied to the coordinate difference in order to
   #  make it "worse" to be off target rather than off in time
   pos_bias = 100.0
   
   #  Get a list of all references we need to do.
   todo = get_files(reference)
   if todo is None:  
      raise IOError, "Could not access any images:  %s" % images

   # The the times, RA's and DEC's of all the comparison spectra:
   obs_times = get_header(images, time_key, iraf.real)
   obs_ras = get_header(images, ra_key, iraf.real)
   obs_decs = get_header(images, dec_key, iraf.real)

   # Just in case the reference files are included in the images list,
   #  get rid of them (otherwise, this is just the identity)
   for this_ref in todo:
      if this_ref in obs_times.keys():
         del obs_times[this_ref]
         del obs_ras[this_ref]
         del obs_decs[this_ref]

   to_return={}

   # Now go through each file and find the closest in time and position.
   for this_ref in todo:
      this_time = get_header(this_ref, time_key, iraf.real)[this_ref]
      this_ra = get_header(this_ref, ra_key, iraf.real)[this_ref]
      this_dec = get_header(this_ref, dec_key, iraf.real)[this_ref]

      # Get a list of differences in observing times:
      deltas = {}
      for obj in obs_times:
         deltas[obj] = abs(this_time - obs_times[obj])

      # Get a list of differences in observing positions
      dists = {}
      for obj in obs_ras:
         deltas[obj] = deltas[obj] + pos_bias*math.sqrt((this_ra - \
               obs_ras[obj])**2 + ((this_dec - obs_decs[obj])/15.0)**2)

      # deltas is now a dictionary to which each objects is associated
      # the difference in position plus the difference in position (all
      # in units of hours).  We now sort on this value and the smallest
      # will be the comparison we want.
      sorted_objs = sort_dict(deltas)
      
      # Now find the closest position
      to_return[this_ref] = sorted_objs[0]
      if printit:
         obj = to_return[this_ref]
         print '%s: %s  delta_t=%.4f sec, delta_position=%.2f sec' % \
            (this_ref,obj, 3600.0*(this_time - obs_times[obj]),
             3600.0*math.sqrt((this_ra - obs_ras[obj])**2 + \
                  ((this_dec - obs_decs[obj])/15.0)**2))



   if len(to_return.keys()) == 1:
      return to_return[to_return.keys()[0]]
   else:
      return to_return

def doDarks(dark, readnoise, images):
   '''This function will take the provided name of the dark frame (dark),
      determine the dark signal (median of the frame) and compute
      the dark current expected for each file name in images.  This is then
      compared to the readnoise.  A list of all images in which the dark 
      current is expected to be greater than the readnoise is printed,
      separated by commas, so you can cut-n-paste into ccdproc's epar.'''

   if not iraf.imaccess(dark):
      print "Error:  can't open dark frame %s" % (dark)
      return

   # First off, get the dark current and exposure time
   dc = float(iraf.imstat(images=dark, fields='mean', Stdout=1)[1])
   dt = get_header(dark, 'darktime', float)[dark]
   print 'dark:  count=%f   exptime=%f   dark current rate=%f' % (dc,dt,dc/dt)

   # Next, get the exposure times for each image
   res = get_header(images=images, field='exptime', type=float)

   # a list to be printed out at the end (good for cut-n-paste)
   names = []
   print "expected dark currents:"
   for name in res.keys():
      print name,"%.2f" % (res[name]*dc/dt)
      if res[name]*dc/dt > readnoise:
         names.append(name)
   print "expected dark current greater than read noise:"
   print string.join(names, ',')


def extractComps(assignments, dorefspec=1):
   '''This task takes a dictionary 'assignments' which has keywords equal
      to the objects spectrum and value equal to the comparison spectrum
      which should be used for the object spectrum.  You could use 
      closest_time_pos() to figure this out, for example.  The apertures
      of each comp are then extracted using apall and the object's spectrum
      as a reference.  Finally, the object's headers are updated using
      the REFSPEC1 keyword (unless dorefspec=0).'''
   iraf.noao();  iraf.imred();  iraf.echelle()
   # First things first, we save the old paramter list.
   iraf.apall.saveParList(filename='temp.apall.par')

   # iraf.apall.unlearn()
   iraf.apall.format = 'echelle'

   # now, get the associated list of references and comps:
   stuff = assignments

   for ref in stuff.keys():
      print 'extracting %s using %s as reference' % (stuff[ref], ref)
      comp = stuff[ref]
      if iraf.imaccess(re.sub('\.fits','.ec.fits',comp)):
         print 'Extracted comparison exists, skipping...'
      else:
         iraf.apall(input=comp,
                 references=re.sub('\.ec','',ref),
                 interactive=0,
                 find=0,
                 recenter=0,
                 resize=0,
                 edit=0,
                 trace=0,
                 fittrace=0,
                 extract=1,
                 extras=0,
                 review=0)

      # Now, since we've used them to extract the spectra, we can assign
      # them using refspec:
      # print "now assigning reference spectrum"
      if dorefspec:  iraf.hedit(images=ref, fields='REFSPEC1', 
                 value=re.sub('\.fits','.ec',comp), add=1, addonly=1,
                 verify=0, show=1)


   iraf.apall.setParList(ParList='temp.apall.par')


def precess_coords(file, epoch=2000.0):
   '''This function takes a file, accesses its header keywords (RA,DEC and
      EPOCH), precesses the coordinates to the proviced epoch and returns
      the RA and DEC as a tuple of string representations.'''

   iraf.noao();  iraf.astutil()
   if not iraf.imaccess(file):
      print "Error:  can't open image frame %s" % (file)
      return None

   ra = get_header(file, 'RA', str)[file]
   dec = get_header(file, 'DEC', str)[file]
   in_epoch = get_header(file, 'EPOCH', str)[file]

   f = open('iraf.tempfile', 'w')
   f.write('%s %s %s %.1f\n' % (ra,dec,in_epoch,epoch))
   f.close()
   result = iraf.precess(input='iraf.tempfile', startyear=1950, endyear=2000, 
            Stdout=1)
   #os.unlink('iraf.tempfile')
   return string.split(result[0])[0:2]

def closer_than(ra1,dec1,ra2,dec2, tol=5):
   '''Returns 1 if ra2,dec2 is closer than tol arc-sec from ra1,dec1 and
   0 otherwise.  ra?  is assumed to be in hours, dec? in degrees.'''
   ra1 = iraf.real(ra1)
   dec1 = iraf.real(dec1)
   ra2 = iraf.real(ra2)
   dec2 = iraf.real(dec2)

   dist = math.sqrt(((ra1 - ra2)*15.0)**2 + (dec1 - dec2)**2)*3600.0

   if dist < tol:
      return 1
   else:
      return 0

