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
``papi``                 Ingest uncalibrated data and builds a SQLite database containing fits header data
``calDark``              Removes electronic ghosts (a.k.a. the "Mr. Staypuft" effect)
``calDarkModel``         Basic calibration, sky subtraction/pedestal effect removal, cosmic-ray rejection
``calBPM``               Removes cosmic-ray persistence signal if target was observed shortly after the SAA
``calDomeFlat``          Removes residual instrument signatures by subtracting a "super-median" reference image
``calTwFlat``            Uses bicubic spline to flatten the background
``calSuperFlat``         Corrects for count-rate non-linearity
``calGainMap``           Apply user-defined masks (optional)
``calNonLinearity``      Uses object matching algorithms to improve image alignment and registration
``checkQuality``         Creates accurate *RMS* maps for use with MultiDrizzle
``eval_focus_serie``     Creates final CR-cleaned, distortion-free drizzled image mosaics using multidrizzle
``astrowarp``            Create final aligned and coadded frame using SEx,SCAMP and SWARP 
====================     ===========

.. index:: setup, sqlite

``papi``
*********

The ``papi`` module in the main PAPI routine that start the data reduction.  
It starts by creating a subdirectory in the ``PAPI_DIR`` using the run name 
give on the command line.  Within the run directory the following
subdirectories are created:


.. autmodule:: papi
   :members:


.. tabularcolumns:: |r|l|

=========   ===========  
Directory   Description
=========   ===========
``CALIB``   A copy of all the uncalibrated input data and the output processed data products for modules ``undupuft`` through ``nonlincor``
``ALIGN``   Output processed data products for modules ``weightmap`` through ``mdrizzle``
``FINAL``   The final data products (final image mosiacs, weightmaps and context images)
=========   ===========

NICRED creates a SQLite_ database to store the uncalibrated input data fits headers and pipeline metadata:

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

.. index:: undopuft, staypuft

``undopuft``
************

The ``undopuft`` module attempts to remove the electronic ghosts that can appear when observing
a bright source. For details see `Electronic Ghosts: Mr. Staypuft, Ringing, and Streaking 
<http://www.stsci.edu/hst/nicmos/performance/anomalies/staypuft.html>`_. 

.. index:: calped, calnica, pedsky, cridcalc, multiaccum, calibration

``calped``
**********

The ``calped`` module performs basic instrumental calibration (dark current subtraction, flat fielding, 
conversion to count rates, and cosmic ray identification and rejection) and attempts to remove 
the NICMOS pedestal effect. These task are performed by the STSCI IRAF package tasks calnica_ and pedsky_. 

The NICMOS pedestal effect is the result of variable biases in each of the four NICMOS detector quadrants these 
varying bias levels can leave a significant pedestal signature in the processed data. For details see the 
NICMOS anomaly page `Residual Bias (Pedestal) <http://www.stsci.edu/hst/nicmos/performance/anomalies/pedestal.html>`_

NICRED runs all of the calibration steps provided by calnica_ in the default sequence with the exception of one 
additional step. Before the calnica_ cosmic ray identification and removal step *CRIDCALC* is run NICRED runs an
additional step to improve the cosmic ray rejection. For NICMOS MultiAccum mode observations, *CRIDCALC* assumes 
that accumulating background counts over the entire observation is a linear function. This assumption may not 
be the true for all observations. Depending on circumstances of the observation the background 
count rate may vary over the duration of the observation. In order to determine if the background count rate 
is sufficiently non-linear, NICRED computes the median of the first and last three readouts of the MultiAccum 
observation.  If the NIRCED finds the count rate has varied it applies the additional step of running pedsky_ 
on each of the individual readouts in the MultiAccum observation. This additional step assures the background 
count rate is linear before running the *CRIDCALC* step. 

.. _calnica: http://stsdas.stsci.edu/cgi-bin/gethelp.cgi?calnica

.. _pedsky: http://stsdas.stsci.edu/cgi-bin/gethelp.cgi?pedsky.hlp

.. index:: saaclean, SAA, pyraf

``saaclean``
************

The ``saaclean`` module removes cosmic ray persistence due to observations following an HST transit of 
the South Atlantic Anomaly (SAA). See `Removing Post-SAA Persistence in NICMOS Data 
<http://www.stsci.edu/hst/nicmos/documents/isrs/isr_2007_001.pdf>`_. 
NICRED uses the PyRAF_ task ``saaclean`` to perform this processing.

.. _PyRAF: http://www.stsci.edu/resources/software_hardware/pyraf


``medsub``
**********

The ``medsub`` module attempts to remove any residual instrument signature left over after basic calibration 
by subtracting a "super-median" reference image. These super-median images are created by median stacking a large 
number of images that have been processed by the NICRED modules ``undopuft``, ``calped`` and ``saaclean``. 
Many super-median reference images (based on various camera, sample-sequence, observation window or HST proposal ID, 
and filter combinations) have been pre-generated and are provided in `nicred_reffiles.tgz 
<http://www.firstgalaxies.org/downloads/nicred/nicred_reffiles.tgz>`_. 

``medsub`` uses the following criteria for determining which super-median image to use:

1. Same camera.
2. Same sample sequence.
3. Same filter.
4. Same HST proposal ID (PROP_ID fits header keyword). Or..
5. The super-median reference image with and observation date nearest the observation date of the input image.


``flatten``
***********

The ``flatten`` module attempts to remove any discontinuities between the four quadrants of a NICMOS camera 2 or 3 image. 
Discontinuities between quadrants can occur when an exposure contains a large bright object in one of the quadrants. 


.. index:: nonlincor, non-linearity, count-rate

``nonlincor``
*************

The ``nonlincor`` module corrects NICMOS images for their count-rate dependent 
non-linearity. It used the header keywords CAMERA and FILTER to determine the 
non-linearity parameter. It corrects the first image, and in the case of a 
multi-extension image, the second image as well, with the appropriate power law. 
For details see `Correcting the NICMOS count-rate 
dependent non-linearity <http://www.stsci.edu/hst/nicmos/documents/isrs/isr_2006_003.pdf>`_

.. index:: mask, masking, applymask, ds9

``applymask``
*************

NICRED has the ability to mask any residual artifacts that may occur in one's 
data (e.g., as may occur when satellites pass through the HST focal plane). 
Masks are easily generated using `SAO’s DS9 <http://http://hea-www.harvard.edu/RD/ds9>`_ 
image display tool using the following procedure:

    1. Display all ``_cal4.fits`` images in DS9.
    2. Marked artifacts on each image with the DS9 polygon region tool.
    3. A script is run that saves a DS9 region file for each image which has a marked artifact.
    4. A second script is run that applies the marked regions in these region files to the associated image’s data-quality array.


.. image:: _static/applymask.png
    :align: center
    
``align``
*********

The ``align`` module uses the external package ``superalign`` to determine the internal
shifts and rotations for an arbitrary number of (overlapping)
contiguous images from a set of (distortion free) catalogs.  It
requires good initial guesses for the shifts and rotations (within 2.5
arcsec and 0.5 degrees of the true solution, respectively), and thus
is ideal for use with NICMOS HST data where these quantities are
approximately known.  It offers several useful advantages relative to
other alignment programs:

    1. It does not require that all images be contiguous with a single reference image. 
       This allows one to construct arbitrarily large mosaics out of individual images.

    2. Input catalogs can include substantial (>80%) contamination from cosmic rays.

For more details on ``superalign`` see the Appendix section :ref:`superalign`. 

.. index:: multidrizzle, weightmap, variance

``weightmap``
*************

The ``weightmap`` module generates an inverse variance weigh map image of each input image as input to MultiDrizzle.

.. math::spee

        Var\; =\; \frac{\left( \left( {D}\; +\; {A} \right)\times {G}\; +\; \left( {B}\times {f} \right)\times {G}\; +\; {\sigma_{read}}^{2} \right)}{\left( {f}^{2}\; \times \; {t}^{2} \right)}

        W\; =\; \frac{1}{Var}

Where *D* is the dark image; *A* is the amplifier glow image; *G* is the gain; *B* is the average background as computed by ``calnica``; *sigma* is the readnoise; *f* is the inverse flatfield image; and *t* is the exposure time. 


.. index:: multidrizzle, dithering, weightmap, variance

``mdrizzle``
************

The ``mdrizzle`` module performs cosmic ray rejection and combination of dithered observations using the STSCi software package MultiDrizzle. 
For a complete discussion of MultiDrizzle and the Drizzle alorgithm for combining dithered imaging data see the `MultiDrizzle Handbook Wiki  
<http://incubator.stsci.edu/mediawiki/index.php/Telescopedia:Multidrizzle>`_.