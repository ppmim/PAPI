Modules
=======

The PAPI pipeline includes N processing steps from basic calibration 
to generating final co-added registered mosaics (See table below). Among these
steps are steps for successfully handling three of the most common PANIC 
anomalies: electronic ghosting, the pedestal effect and cosmic-ray persistence. 
Electronic ghosts occur when a bright object is observed in one of the four 
quadrants on a PANIC detector. This results in an echo of the bright object 
in the other three quadrants. The module undopuft attempts to remove these 
artifacts. The pedestal effect is the result of variable biases in each of the four 
PANIC quadrants - leaving a significant pedestal signature in the processed 
data. 

.. index:: modules

.. tabularcolumns:: |r|l|

====================     ===========
Module                   Description
====================     ===========
``papi``                 Main pipeline module to start the entire data reduction process 
``calDark``              Creates a master dark by combination of a dark sequence
``calDarkModel``         Creates a master dark model from a dark serie
``calBPM``               Creates a master Bad Pixel Mask from a set of darks and flats calibration files
``calDomeFlat``          Creates a master Dome Flat 
``calTwFlat``            Creates a master Twilight Flat
``calSuperFlat``         Creates a master Super Flat from a set of object or sky frames
``calGainMap``           Creates a Gain Map from any master flat
``calNonLinearity``      Corrects the images pixel values for non-linearity
``checkQuality``         Computes some quality values from the image (FWHM, STD, RMS)
``applyDarkFlat``        Finds out the best Focus value from a focus serie
``eval_focus_serie``     Finds out the best Focus value from a focus serie
``skyfilter``            Subtracts sky background to a dither sequence of frames
``astrowarp``            Creates final aligned and coadded frame using SEx,SCAMP and SWARP 
``photometry``           Performs a photometric calibration comparison with 2MASS
``genLogsheet``          Creates a log sheet from a set of FITS files
``makeobjmask``          Creates a objects mask (SExtractor OBJECTS images) for a list of FITS images.
====================     ===========

.. index:: setup, sqlite

``papi``
********

The ``papi`` module is the main PAPI routine that start the data reduction.  
It starts by creating a subdirectory in the ``PAPI_DIR`` using the run name 
give on the command line.  Within the run directory the following
subdirectories are created:


.. automodule:: papi
   :members:


.. tabularcolumns:: |r|l|

=========   ===========  
Directory   Description
=========   ===========
``CALIB``   A copy of all the uncalibrated input data and the output processed data products for modules ``undupuft`` through ``nonlincor``
``ALIGN``   Output processed data products for modules ``weightmap`` through ``mdrizzle``
``FINAL``   The final data products (final image mosiacs, weightmaps and context images)
=========   ===========

PAPI creates a SQLite_ database to store the uncalibrated input data fits headers and pipeline metadata:

.. index:: log, logging, status, FITS, headers

.. tabularcolumns:: |r|l|

==============   ===========  
Table            Description
==============   ===========
``headers``      Select FITS header keywords for all input images 
``run_log``      Runtime log messages
``run_pars``     Value of each runtime option
``run_status``   Runtime status information for each module
``raw``          A *VIEW* of the header table listing only the raw FITS image headers
==============   ===========

.. _SQLite: http://www.sqlite.org

.. index:: papi

``calDark``
***********

The ``calDark`` module creates a master dark image from a set of dark frames.
In addition computes the Read-Out Noise of the detectors along with several
statistics.

Input

.. index:: calped, calnica, pedsky, cridcalc, multiaccum, calibration

``calDarkModel``
****************

The ``calDarkModel`` module performs a dark model. To do that, a input dark serie
exposures with a range of exposure times is given. Then a linear fit is done at 
each pixel position of data number versus exposure time. A each pixel position 
in the output map represents the slope of the fit done at that position and is 
thus the dark current expressed in units of data numbers per second.

.. index:: , SAA, pyraf

``calBPM``
**********

The ``calBPM`` module 

``calDomeFlat``
***************

The ``calDomeFlat`` module creates a master flat field from dome flat observations,
a bad pixel map an various statistics.

``calDomeFlat`` uses the following criteria for determining which super-median image to use:

1. Same camera.
2. Same sample sequence.
3. Same filter.
4. Same PANIC proposal ID (PROP_ID fits header keyword). Or..
5. The super-median reference image with and observation date nearest the observation date of the input image.



``calTwFlat``
*************

The ``calTwFlat`` module creates a master flat field from twilight observations,
a bad pixel map an various statistics.

.. index:: flat-field, twilight 

``calSuperFlat``
****************

The ``calSuperFlat`` module creates a master super flat field from science observations,
a bad pixel map an various statistics.

.. index:: flat-field, super-flat 


``applyDarkFalt``
*****************

The ``applyDarkFalt`` module subtract a master dark and divide by the list of images
given.

.. index:: flat, dark, calibration 


``calGainMap``
****************

The ``calGainMap`` module creates a master gain map from a flat field 

.. index:: flat-field, super-flat 


``calNonLinearity``
*******************

The ``calNonLinearity`` module corrects PANIC images for their count-rate dependent 
non-linearity. It used the header keywords READMODE and FILTER to determine the 
non-linearity parameter. It corrects the first image, and in the case of a 
multi-extension image, the second image as well, with the appropriate power law. 
For details see `Correcting the PANIC count-rate 
dependent non-linearity <http://www.iaa.es/PANIC/papi/documents/nonlinearity.pdf>`_

.. index:: mask, masking, applymask, ds9

``checkQuality``
****************
The ``checkQuality`` module computes some initial image quality estimations using 
SExtractor.

.. index:: fwhm, seeing, sextractor

    
``skyfilter``
*************

The ``skyfilter`` module uses the external package ``irdr_skyfilter`` to perform the
sky background subtraction from a dither sequence of science frames. It works
with almost all kind of dither sequences, even with sequences used for extended
objects (T-S-T-S- ...., T-T-S-T-T-S-T-T-....)

For more details on ``skyfilter`` see the Appendix section :ref:`skyfilter`. 

.. index:: sky-background, irdr, sky

``eval_focus_serie``
********************

The ``eval_focus_serie`` module computes the best focus estimation for a focus
exposure serie. It is done according to the FWHM value estimated for each
frame, fitting a curve the the values pair values (FWHM,focus) and finding out the 
minimun.

- Requirements

    - T-FOCUS (telescope focus) keyword value present in the header 
    - (Raw) Images with enought number of stars
    - A serie of images taken with covering a range of telescope focus values including the best focus value.
 

.. index:: focus, fwhm, seeing

``astrowarp``
*************

The ``astrowarp`` module performs alingment and warping of a set of images previosuly reduced. 
That module uses the Astromatic_ packages sextractor_ , scamp_ and swarp_
to perform the...

.. _astromatic: http://www.astromatic.net/
.. _sextractor: http://www.astromatic.net/software/sextractor
.. _scamp: http://www.astromatic.net/software/scamp
.. _swarp: http://www.astromatic.net/software/swarp


