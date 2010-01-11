################################################################################
#
# FIEStool
#
# frameUtils.py
#
# Last update 08/08/2006
#
################################################################################

"""
   Frame-related utility routine(s) that make life easier.
"""

# Import external modules

import os
import numpy
import pyfits
import config
import messageLog

import irafDbNames

################################################################################

# Helper functions :

def _flipaxis(arr, axis=0):

  # Might be done more efficiently with numpy.fliplr() and flipud()
  return arr.take(numpy.arange(arr.shape[axis], 0, -1) - 1, axis=axis)

################################################################################


def WCSvector(frame, extension=0):

  """
     Return vectors of the WCS (World Coordinate System) x and y pixel
     coordinates in this frame. Particularly useful for binned or windowed
     frames.
  """

  # 'frame' is assumed to be an object created by pyfits.open(). Failure
  # to read any of the header keywords generates a KeyError exception, which
  # is supposed to be handled by the calling routine. However, for properly
  # defined frames, exceptions should not need to occur

  # Get the number of pixels for the two axes from the FITS headers
  # No try/except needed, because these headers MUST be defined.
  try:
    xsize = frame[extension].header['NAXIS1']
    ysize = frame[extension].header['NAXIS2']
  except KeyError:
    messageLog.put_error('NAXIS1 or NAXIS2 keyword not defined. Are you working with the correct extension?')
    return

  # Get the remaining vector values, in a truly failsafe way.

  try:
    x0 = frame[extension].header['CRVAL1']
  except KeyError:
    x0 = 1

  try:
    y0 = frame[extension].header['CRVAL2']
  except KeyError:
    y0 = 1
  
  try:
    xref = frame[extension].header['CRPIX1']
  except KeyError:
    xref = 1

  try:
    yref = frame[extension].header['CRPIX2']
  except KeyError:
    yref = 1

  try:
    xstep, dummy = frame[extension].header['CCDSUM'].split()
    xstep = int(xstep)
  except:
    raise
    try:
      xstep = frame[extension].header['CDELT1']
    except KeyError:
      xstep = 1

  try:
    dummy, ystep = frame[extension].header['CCDSUM'].split()
    ystep = int(ystep)
  except:
    raise
    try:
      ystep = frame[extension].header['CDELT2']
    except KeyError:
      ystep = 1
  

  # Construct the actual vectors with x and y pixel coordinates.
  # numpy.arange(n) gives a vector with n numbered elements
  xvec = (x0 + ( (numpy.arange(xsize) - xref + 1) * xstep))
  yvec = (y0 + ( (numpy.arange(ysize) - yref + 1) * ystep))

  # Return a tuple with the two vectors
  return (xvec, yvec)


################################################################################


def createframe(object="", xbin=1, ybin=1, nodata=0, extension=0):

  """
     Create a new frame from scratch using the default frame size and a
     given binning. Returns a PyFITS object.
  """

  # Determine the default frame size
  xmin = config.options['x_start'].value
  xmax = config.options['x_end'  ].value
  ymin = config.options['y_start'].value
  ymax = config.options['y_end'  ].value

  xsize = (xmax - xmin)/xbin + 1
  ysize = (ymax - ymin)/ybin + 1

  messageLog.put('Creating new frame: %s (%i, %i) pixels' % (object, xsize, ysize), 7)


  # Create a list of HDUs (Header Data Units).
  newframe = pyfits.HDUList()

  # ...and create the primary HDU
  primhdu  = pyfits.PrimaryHDU()

  # Set OBJECT header
  primhdu.header.update('OBJECT', object)

  # Attach the primary HDU to the frame
  newframe.append(primhdu)

  # Return if frame should not contain data block
  if nodata: return newframe

  # Also create an image HDU
  if (extension != 0):

    for extno in range(1, extension+1) :

      imhdu = pyfits.ImageHDU()

      # Attach the image HDU to the frame
      newframe.append(imhdu)

      # Create EXTNAME header (to satisfy IRAF)
      newframe[extno].header.update('EXTNAME', 'im%i' % extno)


  # Fill the data unit with an empty frame of requested size
  newframe[extension].data = numpy.zeros((ysize, xsize), dtype=numpy.float64)

  # Update header values accordingly
  newframe[extension].header.update('CRVAL1', xmin)
  newframe[extension].header.update('CRPIX1', 1)
  newframe[extension].header.update('CDELT1', xbin)
  newframe[extension].header.update('CTYPE1', 'PIXEL')

  newframe[extension].header.update('CRVAL2', ymin)
  newframe[extension].header.update('CRPIX2', 1)
  newframe[extension].header.update('CDELT2', ybin)
  newframe[extension].header.update('CTYPE2', 'PIXEL')

  newframe[extension].update_header()


  # Return the frame object
  return newframe


################################################################################


def flipframe(frame, direction=0, extension=0):

  # Retrieved the x and y coordinate vectors
  (xvec, yvec) = WCSvector(frame, extension)


  # MMM... MAYBE NOT TRUE... 
  # Comment : flipping the xvec and yvec below is not really
  # necessary, because the xnewmin and ynewmin will still be the
  # lowest pixel number in the (flipped) xvec and yvec. However,
  # exchanging the xvec and yvec is _very_ relevant.

  # NB: xaxis is axis nr 1, and yaxis is axis nr 0

  flippeddata = frame[extension].data

  # Rotate 'n' times 90 degrees
  # x ->  x, y ->  y
  if   direction == 0:
    flippeddata = flippeddata

  # x -> -y, y ->  x
  elif direction == 1:
    flippeddata = _flipaxis(flippeddata, axis=1)
    flippeddata.transpose()
    yvec = _flipaxis(yvec, axis=0)

  # x -> -x, y -> -y
  elif direction == 2:
    flippeddata = _flipaxis(flippeddata, axis=1)
    flippeddata = _flipaxis(flippeddata, axis=0)
    xvec = _flipaxis(xvec, axis=0)
    yvec = _flipaxis(yvec, axis=0)

  # x ->  y, y -> -x
  elif direction == 3:
    flippeddata = _flipaxis(flippeddata, axis=0)
    flippeddata.transpose()
    xvecnew = yvec
    yvec = _flipaxis(xvec, axis=0)
    xvec = xvecnew
    del(xvecnew)

  # Rotate 'n' times 90 degrees, and transpose
  # x ->  y, y ->  x
  elif direction == 4:
    flippeddata.transpose()
    xvecnew = yvec
    yvec = xvec
    xvec = xvecnew
    del(xvecnew)

  # x -> -x, y ->  y
  elif direction == 5:
    flippeddata = _flipaxis(flippeddata, axis=1)
    xvec = _flipaxis(xvec, axis=0)

  # x -> -y, y -> -x
  elif direction == 6:
    flippeddata = _flipaxis(flippeddata, axis=1)
    flippeddata = _flipaxis(flippeddata, axis=0)
    flippeddata.transpose()

  # x ->  x, y -> -y
  elif direction == 7:
    flippeddata = _flipaxis(flippeddata, axis=0)
    yvec = _flipaxis(yvec, axis=0)

  else:
    messageLog.put_warning('Frame orientation parameter out of range - ignored')


  # Determine 'new' vector minima and stepsize
  xnewbin = abs(xvec[2] - xvec[1])
  ynewbin = abs(yvec[2] - yvec[1])

  xnewmin = xvec.min()
  ynewmin = yvec.min()


  # Update header values accordingly
  frame[extension].header.update('CRVAL1', xnewmin)
  frame[extension].header.update('CRPIX1', 1)
  frame[extension].header.update('CDELT1', xnewbin)
  frame[extension].header.update('CTYPE1', 'PIXEL')

  frame[extension].header.update('CRVAL2', ynewmin)
  frame[extension].header.update('CRPIX2', 1)
  frame[extension].header.update('CDELT2', ynewbin)
  frame[extension].header.update('CTPYE2', 'PIXEL')

  # Store the rotated and flipped data in fits object
  frame[extension].data = flippeddata

  # And adjust the image headers to reflect the data
  frame[extension].update_header()


################################################################################


def clipframe(frame, extension=0):

  """
     Remove under- and overscan regions from a frame using the WCS
     pixel coordinates. If a frame is already within the limits,
     leave it unchanged (no padding).
  """

  # 'frame' is assumed to be an object created by pyfits.open().

  # Get the maximum allowed frame size from the current configuration
  xabsmin = config.options['x_start'].value
  xabsmax = config.options['x_end'  ].value
  yabsmin = config.options['y_start'].value
  yabsmax = config.options['y_end'  ].value
 
 
  # Retrieved the x and y coordinate vectors
  (xvec, yvec) = WCSvector(frame, extension)


  # Create new x and y vectors of the pixels coordinates that should be
  # in the output frame. Use boolean logic to select the pixels within bounds
  # Here, 'newxvec' and 'newyvec' contain only 0 or 1.
  newxvec = (xvec >= xabsmin) * (xvec <= xabsmax)
  newyvec = (yvec >= yabsmin) * (yvec <= yabsmax)
 

  # Get a subset of the data for the non-zero elements of the new x and y
  # vectors. Then, remove all non-zero elements of the new x and y vectors.
  #
  # The [0] after nonzero() is there because nonzero() returns a tuple
  # that contains an array (don't know why - bug?)
  #
  # Do not change the order of the statements below! Slicing is done twice.

  # Get the subset of the data along the x axis
  subset  = frame[extension].data.take(newxvec.nonzero()[0], axis=1)
  # Make 'newxvec' contain the actual pixel positions, rather than 0 or 1
  newxvec = xvec.take(newxvec.nonzero()[0])
  # Force the new vector to be one-dimensional (take returns 2D)
  newxvec.shape = (-1,)

  # Idem for y
  subset  = subset.take(newyvec.nonzero()[0], axis=0)
  newyvec = yvec.take(  newyvec.nonzero()[0])
  newyvec.shape = (-1,)

  # 'subset' only contains pointers (!) to (part of) the data in 'data'.
  # Therefore, make a real copy of the data, that can be stored in an output
  # file
  newdata = subset.copy()

  # Adjust the 2D shape of newdata
  newdata.shape = (len(newyvec), len(newxvec))

  # Determine if there is binning from the pixel positions
  newxbin = newxvec[2] - newxvec[1]
  newybin = newyvec[2] - newyvec[1]

  # Store the sliced data in the original frame
  frame[extension].data = newdata

  # And update the headers accordingly
  frame[extension].header.update('CRVAL1', newxvec[0])
  frame[extension].header.update('CRPIX1', 1)
  frame[extension].header.update('CDELT1', newxbin)
  frame[extension].header.update('CTYPE1', 'PIXEL')

  frame[extension].header.update('CRVAL2', newyvec[0])
  frame[extension].header.update('CRPIX2', 1)
  frame[extension].header.update('CDELT2', newybin)
  frame[extension].header.update('CTYPE2', 'PIXEL')

  frame[extension].update_header()


################################################################################


def padframe(frame, extension=0):

  """
     Extend an existing frame to the default image size. Returns 0 if the frame
     is already larger or equal to the default size, and 1 if the frame was
     padded.
  """

  # Get pixel numbering of this frame (from WCS)
  (xvec, yvec) = WCSvector(frame, extension)

  # Get current binning factor
  xbin = xvec[2] - xvec[1]
  ybin = yvec[2] - yvec[1]

  # Determine current frame size
  xmin = min(xvec)
  xmax = max(xvec)
  ymin = min(yvec)
  ymax = max(yvec)
  xsize = (xmax - xmin)/xbin + 1
  ysize = (ymax - ymin)/ybin + 1

  # Determine the maximum allowed frame size
  xnewmin = config.options['x_start'].value
  xnewmax = config.options['x_end'  ].value
  ynewmin = config.options['y_start'].value
  ynewmax = config.options['y_end'  ].value

  # And calculate the expected frame size from these numbers
  # (Use same algorithm as in createframe and createcube)
  xnewsize = (xnewmax - xnewmin)/xbin + 1
  ynewsize = (ynewmax - ynewmin)/ybin + 1


  # Check if padding is really necessary. If not, return with exit code 0
  if ( xsize == xnewsize and ysize == ynewsize ):
    return 0


  # Create new x and y index vectors
  # Note: arange starts with 0
  newxvec = numpy.arange(xnewsize)*xbin + xnewmin
  newyvec = numpy.arange(ynewsize)*ybin + ynewmin


  # Find which part of the new frame should contain the data
  # from the new frame (because old frame is smaller than new frame)
  #
  # The [0] after nonzero is there because nonzero return a tuple
  # that contains an array (don't know why - bug?)

  selxvec = (newxvec >= xmin) * (newxvec <= xmax)
  selxvec = selxvec.nonzero()[0]

  selyvec = (newyvec >= ymin) * (newyvec <= ymax)
  selyvec = selyvec.nonzero()[0]


  # Create a new zero-filled full-sized frame
  newdata = numpy.zeros((ynewsize, xnewsize), dtype=numpy.float64)

  # Transfer the data into the new frame. Note that the upper boundary of a
  # slice in a Python array is exclusive! Therefore +1...
  newdata[min(selyvec):max(selyvec)+1,
          min(selxvec):max(selxvec)+1] = frame[extension].data


  # Store the sliced data in the original frame
  frame[extension].data = newdata

  # Update header values accordingly
  frame[extension].header.update('CRVAL1', xnewmin)
  frame[extension].header.update('CRPIX1', 1)
  frame[extension].header.update('CDELT1', xbin)
  frame[extension].header.update('CTYPE1', 'PIXEL')

  frame[extension].header.update('CRVAL2', ynewmin)
  frame[extension].header.update('CRPIX2', 1)
  frame[extension].header.update('CDELT2', ybin)
  frame[extension].header.update('CTPYE2', 'PIXEL')

  frame[extension].update_header()

  # Return that padding was succesful
  return 1


################################################################################


def createdatacube(framelist, extension=0):

  """
     Read data for a list of images, and store this in a 3-dimensional array
     (numpy type). A 3-D array is practical for performing efficient
     calculations on a large number of images (for example, averaging).
  """

  # Determine the number of frames in the list (will be z-size of cube)
  nframes = len(framelist)

  # Counter to keep track of position in the datacube
  n = 0

  # Read first frame in list to determine frame size (x- and y-size of cube)
  messageLog.put('Reading frame %i out of %i' % (n+1, nframes))
  try:
    image = pyfits.open(framelist[0])
  except IOError:
    messageLog.put_error('Cannot open %s' % framelist[0])
    return None

  image = extractMEF(image, extension=extension)

  # THIS IS NOT THE BEST PLACE TO PUT THESE, BUT IT SHOULD GO SOMEWHERE
  #flipframe(image, config.options['frameorientation'].value)
  #clipframe(image)
  data = image[0].data
  framesize = data.shape

  # Create am empty 3D dataset and insert first frame
  # (Use .copy() to keep the data in memory after the image is closed)
  datacube = numpy.zeros((nframes, framesize[0], framesize[1]), data.dtype)
  datacube[n, :, :] = data.copy()

  # Close the file
  image.close()

  # Next frame...
  n = n + 1

  # Loop over remaining frames (start from index 1, not 0)
  for frame in framelist[1:]:

    messageLog.put('Reading frame %i out of %i' % (n+1, nframes))

    # Read a frame
    try:
      image = pyfits.open(frame)
    except IOError:
      messageLog.put_error('Cannot open %s' % frame)
      return None


    # THIS IS NOT THE BEST PLACE TO PUT THESE, BUT IT SHOULD GO SOMEWHERE
    image = extractMEF(image, extension=extension)
    #flipframe(image, config.options['frameorientation'].value)
    #clipframe(image)

    data = image[0].data

    # Check that the frame size is compatible with first frame
    if data.shape != framesize:
      messageLog.put_error('Not all frames have identical sizes')
      return None

    # Insert the frame in the cube
    # (Use .copy() to keep the data in memory after the image is closed)
    datacube[n, :, :] = data.copy()

    # Close the file
    image.close()

    # Next frame...
    n = n + 1


  # Return the cube
  return datacube


################################################################################


def getorderdef(frame):

  """
     Not so nice routine to read the x and y positions of spectral orders from
     an IRAF definition file, and store these in a dictionary object.
  """

  # Separate the directory and the file name
#  base, frame = os.path.split(frame)

  # Split off the filename extension
#  frame, ext = os.path.splitext(frame)

  # Construct the filename of the IRAF order definition data file that
  # belongs to 'frame'
  databasefilename = irafDbNames.apname(frame)

  # Initialize parameters
  norders = 0
  orderarray = []

  # Open the definition file for reading
  fh = open(databasefilename)


  try:

    # Loop until broken
    while 1 :

      # Read a single line
      line = fh.readline()
      
      # If this was unsuccessful (EOF) then break the loop
      if not line: break

      # Split the line into words (using whitespace)
      words = line.split()

      # Empty line?
      if len(words) == 0 : continue

      # Does this line define the start of a new order definition?
      if words[0] == 'begin'  : 
				# Increase total number of orders by 1
				norders = norders + 1
				# Append empty order definition dictionary object
				orderarray.append( {} )

      # Does this line define the center position of the order?
      if words[0] == 'center' : 
    				# Store position of the center pixel and it's
				# corresponding (x-)value in keywords in the
				# dictionary 'orderarray'. Python array numbering
				# starts from 0, so therefore 'norders-1'.
				orderarray[norders-1]['centerpix'] = float(words[2])
				orderarray[norders-1]['centerval'] = float(words[1])

      # Is this line the first of several lines that define the coefficients
      # of the order location curve?
      if words[0] == 'curve' :  
				# Get the number of coefficients of this order
				ncoefs = int(words[1])
				orderarray[norders-1]['ncoefs'] = ncoefs
				
				# Create an array of coefficients
    				coefs = []

    				# Read the coefficients from the file
				for j in range(ncoefs):
			          value = fh.readline()
				  # Store in the coefficient array
			          coefs.append(float(value.strip()))

				# Store the coefficient array in 'orderarray'
				orderarray[norders-1]['coefs'] = coefs


  # Cope with i/o-errors (not implemented)
  except IOError: raise

  # Close the input file
  fh.close
  
  # Return the set of order definition parameters
  return orderarray

################################################################################

def extractMEF(frame, extension=1):

  """
     Extract a FITS Header Data Unit (HDU) from a Multi-Extension FITS (MEF)
     file and return non-MEF data
  """

  if (extension == 0) : return frame

  try: frame[extension]
  except IndexError: return None

  # (NB: frame is an existing PyFITS object)

  # Create a list of HDUs (Header Data Units).
  newframe = pyfits.HDUList()

  # Create and attach a primary HDU to the frame
  newframe.append(pyfits.PrimaryHDU())

  # Copy headers from primary data unit
  for (headername, headervalue) in frame[0].header.items() :
    if (headername == "COMMENT") :
      newframe[0].header.add_comment(headervalue)
      continue
    if (headername == "HISTORY") :
      newframe[0].header.add_history(headervalue)
      continue
    newframe[0].header.update(headername, headervalue)

  # Copy headers from selected image data unit
  for (headername, headervalue) in frame[extension].header.items() :
    if (headername == "COMMENT") :
      newframe[0].header.add_comment(headervalue)
      continue
    if (headername == "HISTORY") :
      newframe[0].header.add_history(headervalue)
      continue
    newframe[0].header.update(headername, headervalue)

  # Remove old headers that may contain incorrect values or that are
  # not allowed in the primary HDU
  del(newframe[0].header['NAXIS'])
  del(newframe[0].header['NAXIS1'])
  del(newframe[0].header['NAXIS2'])
  del(newframe[0].header['BZERO'])
  del(newframe[0].header['BSCALE'])
  del(newframe[0].header['BITPIX'])
  del(newframe[0].header['EXTNAME'])
  del(newframe[0].header['XTENSION'])
  del(newframe[0].header['PCOUNT'])
  del(newframe[0].header['GCOUNT'])

  # Attach the new dataset to the primary HDU
  newframe[0].data = frame[extension].data.copy()
  
  # And fix the FITS output. This means reintroducing some of the keywords that
  # were just deleted, but now with values that agree with the new dataset
  newframe.verify('silentfix')
  newframe[0].update_header()

  return newframe

################################################################################

def splitMEF(fnameMEF, out_filenames):

 """
   Split a MEF file into individual files (one per extension).
   The name of the output files will be given in 'output_filenames'
   The new output filename will be 'filenameorig_N.fits' where N goes 1 to NEXTENSIONS

   Input 
      fnameMEF:  filename of MEF
   Output
      out_filenames : list of output filenames 
      Return: the number of extension extracted to outputfiles
 """

 next = 0
 
 try:
   hdulist = pyfits.open(fnameMEF)
 except IOError:
   print 'Error, can not open file %s' %(fnameMEF)
   return 0

 try:
   if hdulist[0].header['EXTEND']!=True:
     print 'Error, file %s is not a MEF file' %(fnameMEF)
     return 0
 except KeyError:
   print 'Error, file %s is not a MEF file' %(fnameMEF)
   return 0

 try:
   next=hdulist[0].header['NEXTEND']
 except KeyError:
   print 'Error, card NEXTEND not found'

 for i in range(1,next+1):
   sufix="_%d.fits" %i
   out_filenames.append(fnameMEF.replace('.fits', sufix))
   out_hdulist = pyfits.HDUList([pyfits.PrimaryHDU(header=hdulist[i].header, data=hdulist[i].data)])
   out_hdulist.verify('silentfix')
   out_hdulist.writeto(out_filenames[i-1], output_verify='ignore', clobber=True)
   out_hdulist.close(output_verify='ignore')
   del out_hdulist
   print "File %s created " %(out_filenames[i-1])

 
 return next 
 print "End of createMEF"
   
 
 
   

################################################################################

def fixheaders(frame):

  # Make sure header cards GAIN and RDNOISE are defined
  # Gain = 1.0 and RDNOISE = 0.0 are the only feasible assumptions one can
  # make if nothing else is known.

  data = pyfits.open(frame, 'update')
  try:
    gain = data[0].header['GAIN']
  except KeyError:
    data[0].header.update('GAIN', 1.0)
    messageLog.put_warning('Could not find GAIN header: assumed GAIN = 1.0')
    data[0].header.add_history('Inserted missing GAIN header card, assumed GAIN = 1.0')
  try:
    rdnoise = data[0].header['RDNOISE']
  except KeyError:
    data[0].header.update('RDNOISE', 0.0)
    messageLog.put_warning('Could not find RDNOISE header: assumed RDNOISE = 0.0')
    data[0].header.add_history('Inserted missing RDNOISE header card, assumed RDNOISE = 0.0')
  data.close()


################################################################################
