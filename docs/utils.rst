Utilities
=========
Besides the modules for data reduction, PAPI has a set of utilities than
can be used as tools for preparing the data reduction execution; they are:


.. tabularcolumns:: |r|J|

======================   ===========
Utilities                Description
======================   ===========
``check_papi_modules``   Check whether all python modules required by PAPI are installed
``checkQuality``         Computes some quality values from the image (FWHM, STD, RMS)
``collapse``             Collapse (sum) each cube of a list of files into a single 2D image
``eval_focus_serie``     Finds out the best Focus value from a focus series
``genLogsheet``          Creates a log sheet from a set of FITS files
``health``               Compute the Gain and Noise from a set of flat images grouped in packets and with increased level of Integration Time
``imtrim``               Cut/crop edges of the input image
``modFits``              Modifies a keyword inside a FITS header
``skyfilter``            Subtracts sky background to a dither sequence of frames
``spatial_noise``        Compute the Spatial Noise from a set of dark images grouped in pairs with the same Integration Time
======================   ===========


``check_papi_modules``
**********************
Check whether all Python modules required by PAPI are installed. The modules
currently required are:

.. tabularcolumns:: |l|l|

======================   ===========
Module                   Version
======================   ===========
Numpy                    1.6
PyRaf                    1.1
PyFITS                   3.0
Matplotlib               0.98.1
Scipy                    0.10
PyQt4.QtCore             4.8
PyWCS                    1.11
vo                       0.7
Atpy                     0.95
======================   ===========

Example::

    $ check_papi_modules.py 
    
    PAPI Python checking tool
    =========================
    
    Checking Python Version:
    PAPI needs Python Version 2.Y with Y>=2.7
    Your Python version 2.7.3 is fine!
    
    Testing Python module installation for module 'atpy':
    PAPI needs at least version 0.9.5
    Your version 0.9.6 of 'atpy' is fine!
    
    Testing Python module installation for module 'scipy':
    PAPI needs at least version 0.10
    Your version 0.10.1 of 'scipy' is fine!
    
    Testing Python module installation for module 'pywcs':
    PAPI needs at least version 1.11
    Your version 1.11-4.10 of 'pywcs' is fine!
    
    Testing Python module installation for module 'PyQt4.QtCore':
    PAPI needs at least version 4.8
    Your version 4.9.1 of 'PyQt4.QtCore' is fine!
    
    Testing Python module installation for module 'vo':
    PAPI needs at least version 0.7
    Your version 0.8 of 'vo' is fine!
    
    Testing Python module installation for module 'numpy':
    PAPI needs at least version 1.6
    Your version 1.6.2 of 'numpy' is fine!
    
    Testing Python module installation for module 'pyraf':
    PAPI needs at least version 1.1
    Your version 2.0 of 'pyraf' is fine!
    
    Testing Python module installation for module 'matplotlib':
    PAPI needs at least version 0.98.1
    Your version 1.1.0 of 'matplotlib' is fine!
    
    Testing Python module installation for module 'pyfits':
    PAPI needs at least version 3.0
    Your version 3.1.0 of 'pyfits' is fine!


``collapse``
************
Sum the planes of each cube of a list files into a single plane 2D-image.

::

    Usage: collapse.py [options] arg1 arg2 ...

    Options:
    -h, --help            show this help message and exit
    -i INPUT_IMAGE, --input_image=INPUT_IMAGE
                            input cube image to collapse into a 2D image
    -l INPUT_IMAGE_LIST, --input_image_list=INPUT_IMAGE_LIST
                            input image list to collapse into a single 2D image
    -o OUTPUT_FILE, --output_file=OUTPUT_FILE
                            output filename (default = /tmp/out.fits)

Example::

    $ collapse -i /data/mycube.fits -o /data/anymore_a_cube.fits
    
    $ collapse -l /data/list.txt -o /data/anymore_a_cube.fits
    

``checkQuality``
****************
The ``checkQuality`` module computes some initial image quality estimations using 
SExtractor.

.. index:: fwhm, seeing, sextractor


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


``genLogsheet``
***************

``health``
**********

``imtrim``
**********

``modFits``
***********
    
``skyfilter``
*************

The ``skyfilter`` module uses the external package ``irdr_skyfilter`` to perform the
sky background subtraction from a dither sequence of science frames. It works
with almost all kind of dither sequences, even with sequences used for extended
objects (T-S-T-S- ...., T-T-S-T-T-S-T-T-....)

For more details on ``skyfilter`` see the Appendix section :ref:`skyfilter`. 

.. index:: sky-background, irdr, sky


``spatial_noise`` 
*****************


.. _astromatic: http://www.astromatic.net/
.. _Sextractor: http://www.astromatic.net/software/sextractor
.. _scamp: http://www.astromatic.net/software/scamp
.. _swarp: http://www.astromatic.net/software/swarp


