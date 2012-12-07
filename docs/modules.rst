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

.. index:: setup, sqlite


``papi``
********
The papi module (see :ref:`papi`) is the main PAPI script that run the data reduction.  
It starts by creating a subdirectory in the ``PAPI_PROD`` directory using the 
run name give on the command line.  Within the run directory the following
sub-directories are created:


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

.. index:: papi

``applyDarkFlat``
*****************
This module receives a series of FITS images and applies basic calibration: 
subtract and divide by the given calibration files (master dark and master flat-field).

Options::

      -h, --help            show this help message and exit
      -s SOURCE_FILE_LIST, --source=SOURCE_FILE_LIST
                            Source file listing the filenames of raw frames
      -d DARK_FILE, --dark=DARK_FILE
                            Master dark to be subtracted
      -f FLAT_FILE, --flat-field=FLAT_FILE
                            Master flat-field to be divided by
      -o OUT_DIR, --out_dir=OUT_DIR
                            Directory where output files will be saved

``astrowarp``
*************

The ``astrowarp`` module performs the alignment and warping of a set of images,
in principle previously reduced, but not mandatory. 
The module uses the Astromatic_ packages sextractor_ , scamp_ and swarp_
to accomplish this task.

Usage::

    Options:
      -h, --help            show this help message and exit
      -c CONFIG_FILE, --config_file=CONFIG_FILE
                            config file
      -s SOURCE_FILE, --source=SOURCE_FILE
                            Source file list of data frames. It can be a file or directory name.
      -o OUTPUT_FILENAME, --output=OUTPUT_FILENAME
                            final coadded output image
      -v, --verbose         verbose mode [default]


Example::

    $ astrowarp.py -c papi.cfg -s /tmp/test_files.txt -o /tmp/astrowarp.fits

``calBPM``
**********

This module creates a master Bad Pixel Map (.pl iraf file) from a set of dome (on and off) flats.

The algorithm followed to create the BPM is the next:

     1. Classify/split the frames in 3 sets (DOME_FLAT_LAMP_ON, DOME_FLAT_LAMP_OFF, DARKS) and and check whether there are enough calib frames
     2. Check the master dark (Texp)
     3. Subtract the master dark to each dome flat
     4. Combine dome dark subtracted flats (on/off)
     5. Compute flat_low/flat_high
     6. Create BPM (iraf.ccdmask)

Usage::

    Options:
      -h, --help            show this help message and exit
      -s SOURCE_FILE_LIST, --source=SOURCE_FILE_LIST
                            list of input (optionally  corrected) dome ON and OFF flat images..
      -o OUTPUT_FILENAME, --output=OUTPUT_FILENAME
                            The output bad pixel mask.
      -L LTHR, --lthr=LTHR  The low rejection threshold in units of sigma [default 20]
      -H HTHR, --hthr=HTHR  The high rejection threshold in units of sigma [default 20]
      -D MASTER_DARK, --master_dark=MASTER_DARK
                            [Optional] Master dark frame to subtract
      -S, --show_stats      Show statistics [default False]
      -v, --verbose         verbose mode [default]
    

Example::
    
    $ calBPM.py -s /tmp/domesF.txt -D /tmp/masterDark.fits -o /tmp/masterBPM.pl
    
    

``calCombineFF``         
****************
Combine a master dome Flat-field and a master sky Flat-field into a combined
master Flat-field. The procedure followed is :

The procedure for taking advantage of the facts that the large-scale flat-field
variation of the dark-sky flat match that of the program frames and the dome 
flats have very high S/N in each pixel goes as follows:
 
(a) Median smooth the combined, dark-sky flat -this improves the S/N and
preserves the large-scale features of the flat.

(b) Median smooth the combined dome flats using the same filter size as was
used for the dark-sky flat.

(c) Divide the combined dome flat by it's median smoothed-version. The result is
a frame that is flat on large scales but contains all the high spatial frequency
flat-field information.

(d) Now multiply the smoothed dark-sky frame and the result of the division in
the previous step.


As result a flat-field with the low spatial frequency properties of the dark-sky 
flat combined with the high S/N, high spatial frequency properties of the dome 
flat is obtained.

Usage::

    $ calCombineFF.py [options] arg1 arg2 ...
    
    Module to combine a dome Flat-field and a sky Flat-field.
    
    Options:
      --version             show program's version number and exit
      -h, --help            show this help message and exit
      -d DOMEFF, --domeFF=DOMEFF
                            input dome Flat-Field
      -s SKYFF, --skyFF=SKYFF
                            input sky Flat-Field
      -o OUTPUT_IMAGE, --output=OUTPUT_IMAGE
                            output filename of combined Flat-Field (default = combinedFF.fits)

Example::

    $ calCombineFF.py -d /data/masterDF.fits -s /data/masterSF.fits -o /data/masterFF.fits
                   
``calDark``
***********

The ``calDark`` module receives a series of FITS images (master darks) and
create the master dark and computer several statistics.

Usage::

    Usage: calDark.py [options] arg1 arg2 ...
    
    Options:
      -h, --help            show this help message and exit
      -s SOURCE_FILE_LIST, --source=SOURCE_FILE_LIST
                            Source file listing the filenames of dark frames.
      -o OUTPUT_FILENAME, --output=OUTPUT_FILENAME
                            final coadded output image
      -n, --normalize       normalize master dark to 1 sec [default False]
      -e, --scale           scale raw frames by TEXP [default False]
      -S, --show_stats      Show frame stats [default False]
      -t, --no_type_checking
                            Do not make frame type checking [default False]
      -v, --verbose         verbose mode [default]
    
       Usage: calDark.py [options] arg1 arg2 ...
   

Example::

   $ calDark.py -s /data/PANIC_V0/dark_seq.txt -o /data/out/masterDark.fits


.. index:: dark, calibration

``calDarkModel``
****************

The ``calDarkModel`` module performs a dark model. To do that, a input dark series
exposures with a range of exposure times is given. Then a linear fit is done at 
each pixel position of data number versus exposure time. A each pixel position 
in the output map represents the slope of the fit done at that position and is 
thus the dark current expressed in units of data numbers per second.
The dark model obtained will be a FITS files with two planes (extensions): 
    
    * plane 0 = dark current in DN/sec
    * plane 1 = bias
        
    DARKCURRENT value: The median dark current in data numbers per second found 
    from the median value of the output dark current map.



Usage::

    Usage: calDarkModel.py [options] arg1 arg2 ...

    Options:
      -h, --help            show this help message and exit
      -s SOURCE_FILE_LIST, --source=SOURCE_FILE_LIST
                            Source file listing the filenames of dark frames.
      -o OUTPUT_FILENAME, --output=OUTPUT_FILENAME
                            final coadded output image
      -S, --show_stats      Show frame stats [default False]

Example::

    $ calDarkModel.py -s /tmp/darkModel.txt -o /tmp/darkModel.fits

.. index:: dark, calibration


``calDomeFlat``
***************

The ``calDomeFlat`` module creates a master flat field from dome flat observations,
a bad pixel map an various statistics.


Usage::

    Options:
      -h, --help            show this help message and exit
      -s SOURCE_FILE_LIST, --source=SOURCE_FILE_LIST
                            Source file list of data frames. It can be a file or directory name.
      -o OUTPUT_FILENAME, --output=OUTPUT_FILENAME
                            final coadded output image
      -n, --normalize       normalize master flat by median. If image is multi-detector,                  then normalization wrt chip 1 is done) [default False]
      -m, --median_smooth   Median smooth the combined flat-field [default False]
      -v, --verbose         verbose mode [default]


Example::

    $ calDomeFlat -s /tmp/domeFlats.txt -o /tmp/masterDF.fts -n
    

``calSuperFlat``
****************

The ``calSuperFlat`` module creates a master super flat field from science observations,
a bad pixel map an various statistics.


Usage::

    Options:
      -h, --help            show this help message and exit
      -s SOURCE_FILE_LIST, --source=SOURCE_FILE_LIST
                            Source file list of data frames. It has to be a fullpath file name
      -o OUTPUT_FILENAME, --output=OUTPUT_FILENAME
                            output file to write SuperFlat
      -b BPM, --bpm=BPM     bad pixel map file (default=none)
      -N, --norm            normalize output SuperFlat. If image is multi-chip, normalization wrt chip 1 is done (default=True)
      -m, --median_smooth   Median smooth the combined flat-field (default=False)
    
  
Example::

    $ calSuperFlat.py -s /tmp/test_files.txt  -o /tmp/superFlat.fits -N

.. index:: flat-field, super-flat 


``calTwFlat``
*************

This module receives a series of FITS images (twilight flats) and a master dark 
model and creates the master twilight flat-field.


Usage::


    Options:
      -h, --help            show this help message and exit
      -s SOURCE_FILE_LIST, --source=SOURCE_FILE_LIST
                            Source file list of data frames. It can be a file or directory name.
      -d MASTER_DARK, --master_dark_model=MASTER_DARK
                            Master dark model to subtract each raw flat (it will be scaled by TEXP)
      -o OUTPUT_FILENAME, --output=OUTPUT_FILENAME
                            final coadded output image
      -b MASTER_BPM, --master_bpm=MASTER_BPM
                            Bad pixel mask to be used (optional)
      -n, --normalize       normalize master flat by median. If image is multi-detector,then normalization wrt chip 1 is done)[default False]
      -m, --median_smooth   Median smooth the combined flat-field [default False]
      -L MINLEVEL, --low=MINLEVEL
                            flats with median level bellow (default=1000) are rejected
      -H MAXLEVEL, --high=MAXLEVEL
                            flats with median level above (default=100000) are rejected
      -v, --verbose         verbose mode [default]


Example::

    $ calTwFlat.py -s /tmp/twflats.txt -d /tmp/darkModel.fits  -o /tmp/masterTF.fits -n


.. index:: flat-field, twilight 


``calGainMap``
**************

The ``calGainMap`` module creates a master gain map from a master flat field (dome, twilight or superflat)
previously created.

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


``dxtalk``               
**********
Removes cross-talk spots from input images

``makeobjmask``          
***************
Creates a objects mask (SExtractor OBJECTS images) for a list of FITS images.


``photometry``
**************
This module receives a reduced image of any known NIR filter and match to 2MASS 
catalog performing a fit in order to get a estimation of the  Zero Point.
It is based on the method followed here ::

    http://www.ast.cam.ac.uk/ioa/research/vdfs/docs/reports/2masscal.pdf

Usage::

    Options:
      -h, --help            show this help message and exit
      -i INPUT_IMAGE, --input_image=INPUT_IMAGE
                            Input image to calibrate to do photometric comparison with
      -c BASE_CATALOG, --base_catalog (2MASS, USNO-B)=BASE_CATALOG
                            Name of base catalog to compare with (2MASS, USNO-B) -- not used !!! (default = 2MASS)
      -S SNR, --snr=SNR     Min SNR of stars used for linear fit (default = 10.0)
      -z ZERO_POINT, --zero_point=ZERO_POINT
                            Initial Magnitude Zero Point estimation [25.0]; used for SExtractor
      -o OUTPUT_FILENAME, --output=OUTPUT_FILENAME
                            Output plot filename (default = photometry.pdf)

Example::

    $ photometry.py -i /data/reduced.fits -o /tmp/calibration.pdf
    

.. _astromatic: http://www.astromatic.net/
.. _sextractor: http://www.astromatic.net/software/sextractor
.. _scamp: http://www.astromatic.net/software/scamp
.. _swarp: http://www.astromatic.net/software/swarp
.. _SQLite: http://www.sqlite.org


