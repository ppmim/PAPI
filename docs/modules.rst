Modules
=======

The PAPI pipeline consists of a set of processing modules which implement from 
basic calibration to generating final co-added registered mosaics (See table below).
These modules can be run as a stand alone routines depending of your needs, e.g. 
to create a master dark or flat-field for calibration, or you can use them as a
pipeline calling the main script ``papi``, which will use each of the modules 
as they are needed in order to accomplish a complete data reduction of a set of raw images.   
 

.. index:: modules

.. tabularcolumns:: |r|l|

====================     ===========
Module                   Description
====================     ===========
``papi``                 Main pipeline script to start the entire data reduction process 
``applyDarkFlat``        Finds out the best Focus value from a focus series
``astrowarp``            Creates final aligned and coadded frame using SEx, SCAMP and SWARP 
``calBPM``               Creates a master Bad Pixel Mask from a set of darks and flats calibration files
``calCombineFF``         Combine a dome Flat-field and a sky Flat-field into a new Flat-field
``calDark``              Creates a master dark by combination of a dark sequence
``calDarkModel``         Creates a master dark model from a dark series
``calDomeFlat``          Creates a master Dome Flat 
``calSuperFlat``         Creates a master Super Flat from a set of object or sky frames
``calTwFlat``            Creates a master Twilight Flat
``calGainMap``           Creates a Gain Map from any master flat
``calNonLinearity``      Corrects the images pixel values for non-linearity (TBC)
``dxtalk``               Removes cross-talk spots from input images
``makeobjmask``          Creates a objects mask (SExtractor OBJECTS images) for a list of FITS images.
``photometry``           Performs a photometric calibration comparison with 2MASS
====================     ===========

.. tabularcolumns:: |r|l|

=====================    ===========
Utilities                Description
=====================    ===========
``checkQuality``         Computes some quality values from the image (FWHM, STD, RMS)
``check_papi_modules``   Check whether all python modules required by PAPI are installed
``collapse``             Collapse (sum) each cube of a list of files into a single 2D image
``eval_focus_serie``     Finds out the best Focus value from a focus series
``genLogsheet``          Creates a log sheet from a set of FITS files
``health``               Compute the Gain and Noise from a set of flat images grouped in packets and with increased level of Integration Time
``imtrim``               Cut/crop edges of the input image
``modFits``              Modifies a keyword inside a FITS header
``skyfilter``            Subtracts sky background to a dither sequence of frames
``spatial_noise``        Compute the Spatial Noise from a set of dark images grouped in pairs with the same Integration Time
=====================    ===========

.. index:: setup, sqlite


``papi``
********

The ``papi`` module is the main PAPI script that run the data reduction.  
It starts by creating a subdirectory in the ``PAPI_PROD`` directory using the 
run name give on the command line.  Within the run directory the following
sub-directories are created:


.. automodule:: papi
   
   :members:


.. tabularcolumns:: |r|l|

=========   ===========  
Directory   Description
=========   ===========
``CALIB``   A copy of all the calibration created (master dark, master flat, ...)
``TMP``     Temporal data products needed for data reduction, normally purged
``FINAL``   The final data products (final image mosiacs, weightmaps and context images)
=========   ===========

PAPI creates a in-memory SQLite_ database to store the uncalibrated input data fits 
headers and pipeline metadata. 


.. _SQLite: http://www.sqlite.org

.. index:: papi

``calDark``
***********

The ``calDark`` module creates a master dark image from a set of dark frames.
In addition computes several statistics of the detectors.

Options::


   Usage: calDark.py [options] arg1 arg2 ...
   
   Options:
     -h, --help            show this help message and exit
     -s SOURCE_FILE_LIST, --source=SOURCE_FILE_LIST
                           Source file list of data frames. It can be a file or
                           directory name.
     -o OUTPUT_FILENAME, --output=OUTPUT_FILENAME
                           final coadded output image
     -n, --normalize       normalize master dark to 1 sec [default False]
     -e, --scale           scale raw frames by TEXP [default False]
     -v, --verbose         verbose mode [default]

Example::


   $ calDark.py -s /data/PANIC_V0/dark_seq_1/ -o /data/out


.. index:: dark, calibration

``calDarkModel``
****************

The ``calDarkModel`` module performs a dark model. To do that, a input dark series
exposures with a range of exposure times is given. Then a linear fit is done at 
each pixel position of data number versus exposure time. A each pixel position 
in the output map represents the slope of the fit done at that position and is 
thus the dark current expressed in units of data numbers per second.



.. index:: , dark, calibration

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
**************

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

``eval_focus_series``
*********************

The ``eval_focus_series`` module computes the best focus estimation for a focus
exposure series. It is done according to the FWHM value estimated for each
frame, fitting a curve the the values pair values (FWHM,focus) and finding out the 
minimun.

- Requirements

    - T-FOCUS (telescope focus) keyword value present in the header 
    - (Raw) Images with enought number of stars
    - A series of images taken with covering a range of telescope focus values including the best focus value.
 

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


