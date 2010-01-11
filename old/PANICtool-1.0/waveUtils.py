################################################################################
#
# FIEStool
#
# waveUtils.py
#
# Last update 29/09/2005
#
################################################################################

"""
   Wavelength-definition-related utility routines that make life easier.
"""

# Import external modules

import numpy
import messageLog

import irafWrappers

################################################################################


def getwavedef1D(frame):


  """
     Extract and return the wavelength definition vector [ordernumber, zero
     wavelength offset, wavelength step] using the FITS header keywords of a
     1-dimensional spectrum.
  """

  # Read and store the appropriate headers
  try:
    npixels  = frame[0].header['NAXIS1']
    zerowave = frame[0].header['CRVAL1']
    refpix   = frame[0].header['CRPIX1']
    wstep    = frame[0].header['CDELT1']
  except KeyError:
    # Error when reading headers!
    messageLog.put_warning('Could not find FITS headers defining the wavelength solution')
    return None

  # Determine starting (reference) wavelength
  woffset = zerowave + (refpix - 1.0) * wstep

  # And create the wavelength vector
  wavearray   = numpy.zeros((1, 4), dtype=numpy.float64)
  wavearray[0, :]   = [0, 0, woffset, wstep]

  # Return the vector to the caller
  return wavearray


def getwavedef2D(frame):

  """
     Extract and return the wavelength definition vector [ordernumber, zero
     wavelength offset, wavelength step] using the FITS header keywords of a
     2-dimensional frame, following the IRAF definition. It is a rather dirty
     routine.
  """

  # Make sure I can find the appropriate headers
  if not (frame[0].header.get('WAT0_001', None)):
    messageLog.put_warning('Could not find FITS headers defining the wavelength solution')
    return None

  # Get the number of pixels per order and the number of orders from the headers
  npixels = frame[0].header['NAXIS1']
  norders = frame[0].header['NAXIS2']

  # Create an (empty) output array
  wavearray   = numpy.zeros((norders, 4), dtype=numpy.float64)

  # Get all the FITS header keywords in an array
  fullhlist =  frame[0].header.items()
  wavedef = ""

  # And loop over the headers to extract any header that contains (part of)
  # the wavelength definition, and store these in one long string object
  for name, value in fullhlist:
    # Only 68 characters can be part of FITS header value
    if name[0:4] == 'WAT2': wavedef = wavedef + value.ljust(68)

  # Now it gets very IRAF-specific! Bin the first 16 characters (garbage)
  wavedef = wavedef[16:]
  # and split the string object into smaller pieces, one for each order
  wavedef = wavedef.split('spec')

  # Bin the first entry (also garbage)  
  del(wavedef[0])

  # Check that we got the correct number of order definitions from the headers
  if len(wavedef) != norders:
    messageLog.put_warning('No wavelength solution found', 5)
    return None

  # Initialize order number counter
  linecounter = 0

  # Get rid of " and = characters in the string and isolate the parameter values
  for line in wavedef:
    line = line.replace('"', '')
    line = line.replace('=', '')
    line = line.split()

    # According to IRAF's definition the following parameters can be extracted
    # from the line

    number  = int(line[1])
    orderno = int(line[2])
    woffset = float(line[4])
    wstep   = float(line[5])

    # Store the order number, zero pixel wavelength and wavelength step in the
    # output array
    wavearray[linecounter, :]   = [number, orderno, woffset, wstep]
    linecounter = linecounter + 1

  # If an incompatible format was used to define the order locations, complain!
  if int(line[3]) != 0:
    messageLog.put_error('Error in interpreting spectrum definition')
    messageLog.put_warning('Dispersion solution format seems to be non-linear. Trying to make ')
    messageLog.put_warning('the best of it. Wavelength scale may be terribly wrong.')

  # Finally, return the result to the caller
  return wavearray


################################################################################


def helcorr(frame):

  # NB: This routine assumes that the file containes the 'DATE-AVG' keyword,
  #     which is currently only added by the MEF post-processing scripts
  #
  #     To be checked : which epoch are RA and DEC given in, and which role
  #     does EQUINOX play in this?

  """
     Determine the heliocentric velocity correction from the frame headers
  """

  # Read the appropriate headers
  try:
    dateavg  = frame[0].header['DATE-AVG']
  except KeyError:
    # Error when reading headers!
    messageLog.put_warning('Could not find DATE-AVG header: skipping heliocentric velocity calculation')
    return

  try:
    ra       = float(frame[0].header['RA'])
    dec      = float(frame[0].header['DEC'])
    epoch    = float(frame[0].header['EQUINOX'])
  except KeyError:
    # Error when reading headers!
    messageLog.put_warning('Could not read object coordinates: skipping heliocentric velocity calculation')
    return

  # Example of header card :
  # DATE-AVG= '2003-12-13T16:52:12.0' / Midpoint of observation

  # Extract useful data from header
  year  = int(dateavg[0:4])
  month = int(dateavg[5:7])
  day   = int(dateavg[8:10])
  ut    = dateavg[11:]

  # Obtain the heliocentric correction from IRAF
  vhelio = irafWrappers.rvcorrect(year, month, day, ut, ra / 15.0, dec, epoch)

  # Write result to frame header
  frame[0].header.update('VHEL', vhelio)
  frame[0].header.add_history('Added header card for heliocentric velocity')
  messageLog.put('Added header card VHEL with heliocentric velocity', 5)

