################################################################################
#
# FIEStool
#
# plotUtils.py
#
# Last update 03/10/2005
#
################################################################################

"""
   A collection of utility routines that make plotting life easier....
"""

# Import necessary modules

import numpy
import biggles

################################################################################

# Helper function :

def _chebev(c, x):
  """
     Calculate Chebyshev polynomial expansion with coefficients c for points
     (or vector of points) x.
     
     Adopted from Numerical Recipes.
  """

  a = min(x)
  b = max(x)
  m = len(c)
    
  d  = 0.0
  dd = 0.0
  
  y = (2.0 * x - a - b) / (b - a)
  y2 = 2.0 * y
  
  rng = range(1, m)
  rng.reverse()
  
  for j in rng :
  
    sv = d
    d = y2 * d - dd + c[j]
    dd = sv
  
  result = y * d - dd + 0.5 * c[0]
  
  return result


################################################################################

def plot(yvecarr, xvecarr=None, color="black", title="", 
         xsize=600, ysize=300, xrange=None, yrange=None, psfile=None):

  """
     Use biggles to plot the spectrum in the 'yvecarr' array. 'yvecarr' may be
     a 1-dimensional or a 2-dimensional array containing a single spectrum
     or a set of spectra, respectively
  """

  # Set the plotting screen to the demanded size
  biggles.configure('screen', 'width',  xsize)
  biggles.configure('screen', 'height', ysize)

  # Prepare the frame object and give it a title
  frame = biggles.FramedPlot()
  frame.aspect_ratio = 0.5
  frame.title = title


  # Make sure yvecarr is defined and make sense of it
  if yvecarr is None:
    return
  else:
    try:
      nyvec, nypixels = numpy.shape(yvecarr)
    except ValueError:
      # Apparently, yvecarr is a one-dimensional array
      nyvec = 1
      nypixels = len(yvecarr)
      yvecarr = numpy.array([yvecarr])

  # If no xvecarr is defined, fill the x-axis array with increasing
  # integers
  if xvecarr is None:
    xvecarr = numpy.zeros((nyvec, nypixels), dtype=numpy.float64)
    nxpixels = nypixels
    for i in range(nyvec):
      xvecarr[i, :] = numpy.arange(nypixels)
  else :
    # If not, it is either a 1- or 2-dimensional array (similar to yvecarr)
    try:
      nxvec, nxpixels = numpy.shape(xvecarr)
    except ValueError:
      nxvec = 1
      nxpixels = len(xvecarr)
      xvecarr = numpy.array([xvecarr])

  # Make sure that x and y arrays have identical length.
  if (nxpixels != nypixels):
    print 'Unequal length of x- and y-data'
    return


  # Set the plotting x-range if this was specified
  if xrange:
    frame.xrange = xrange
  else :
    frame.xrange = (xvecarr.min(), xvecarr.max())

  # And similar for the y-range
  if yrange:
    frame.yrange = yrange
#  else :
#    frame.yrange = (yvecarr.min(), yvecarr.max())



  # Now, start the loop over the spectra (or spectrum)
  for i in range(nyvec) :

    # Select the x- and y-vectors

    try:
      xvec = xvecarr[i]
    except IndexError:
      xvec = range(0, len(yvec))

    yvec = yvecarr[i]

    # And plot these in the 'frame' object
    frame.add(biggles.Curve(xvec, yvec, color=color))

  # Also add a line indicating the zero-level
  frame.add(biggles.LineY(0, color='blue'))
  
  # And display this plot on screen
  frame.show()

  # Save to encapsulated postscript
  # CURRENTLY NOT IMPLEMENTED. FOR LARGE PLOTS THIS MAY RESULT IN HEAVY DISK I/O
  if psfile is not None:
    frame.write_img(xsize, ysize, psfile)

  # Return the object to the caller, in case it would like to do more with it...
  return frame


################################################################################

def plotorders(orders, interlacedorders=None, xsize=2048, ysize=2048,
               npoints=50, title=''):

  "Plot an overview of the order locations"

  # Set the plotting screen size (square!)
  biggles.configure('screen', 'width',  600)
  biggles.configure('screen', 'height', 600)

  # Find the number of orders and determine if there are any interlaced
  # orders to be plotted
  norders = len(orders)
  if interlacedorders:
    ninterlacedorders = len(interlacedorders)


  # Create the 'frame' object that will hold the plot and give it a title
  frame = biggles.FramedPlot()
  frame.title = title

  # Set the x- and y-ranges
  frame.xrange = (0, xsize)
  frame.yrange = (0, ysize)

  # Loop over orders
  for order in orders:

    # Determine the x and y positions of the orders on the detector
    (xvec, yvec) = _getorderxy(order, npoints=npoints)
    # And plot these into 'frame'
    frame.add(biggles.Curve(xvec, yvec))


  if interlacedorders:

    # Do similar procedure for the interlaced orders
    for order in interlacedorders:

      (xvec, yvec) = _getorderxy(order, npoints=npoints)
      frame.add(biggles.Curve(xvec, yvec, color='blue'))


  # Display the plot on screen
  frame.show()

################################################################################

def _getorderxy(order, npoints=50):

  """
     Simple hack to convert the IRAF order definition to x and y positions
     on the detector. 'npoints' determines how coarse the sampling along the
     orders is allowed to be.
  """

  # This routine can only handle order definitions using Chebychev polynoms!
  if order['coefs'][0] != 1.0 : 
    raise TypeError, 'Wrong type of polynom used in order fitting (not Chebychev)'

  # Get min and max pixel values
  ymin = order['coefs'][2]
  ymax = order['coefs'][3]

  # Determine the stepsize needed to put in 'npoints' points.
  stride = (ymax - ymin) / npoints

  # Get the central pixel offset
  centerpix = int((order['centerpix'] - ymin)/stride)
  centerval = order['centerval']

  # Take the coefficients...
  chebcoef = order['coefs'][4:]

  # ...and use these coefficients to determine the x-positions from the
  # chebychev polynom.
  yvec = numpy.arange(ymin, ymax, stride)
  xvec = _chebev(chebcoef, yvec)

  # Move the x-vector to the correct offset of the central pixel
  try:
    xvec = xvec - xvec[centerpix] + centerval
  except IndexError:
    print "Error in order definition?"
    pass

  # return two sets of arrays that can be plotted
  return (xvec, yvec)
