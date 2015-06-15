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

=======================     ===========
Main Modules                Description
=======================     ===========
``papi``                    Main pipeline script to start the entire data reduction process 
``applyDarkFlat``           Finds out the best Focus value from a focus series
``astrowarp``               Creates final aligned and coadded frame using SEx, SCAMP and SWARP 
``calBPM``                  Creates a master Bad Pixel Mask from a set of darks and flats calibration files
``calCombineFF``            Combine a dome Flat-field and a sky Flat-field into a new Flat-field
``calDark``                 Creates a master dark by combination of a dark sequence
``calDarkModel``            Creates a master dark model from a dark series
``calDomeFlat``             Creates a master Dome Flat 
``calSuperFlat``            Creates a master Super Flat from a set of object or sky frames
``calTwFlat``               Creates a master Twilight Flat
``calGainMap``              Creates a Gain Map from any master flat
``correctNonLinearity``     Corrects the images pixel values for non-linearity
``dxtalk``                  Removes cross-talk spots from input images
``makeobjmask``             Creates a objects mask (SExtractor OBJECTS images) for a list of FITS images.
``photometry``              Performs a photometric calibration comparison with 2MASS
``solveAstrometry``         Performs a astrometric calibration using Astrometry.net and 42xx index files
``remove_cosmics``          Detects and clean cosmic ray hits on images based on Pieter van Dokkum's L.A.Cosmic algorithm.
``eval_focus_serie``        Estimates the best focus value of a focus exposures
``cleanBadPix``             Cleans masked (bad) pixels from an input image. 
=======================     ===========


.. tabularcolumns:: |r|l|

=======================     ===========
Utilities                   Description
=======================     ===========
``createDataSeq``           Modifies headers of a set of FITS files to create a Data Sequece compliant with PAPI
``getBPM``                  Creates the BPM file from the NonLinearity correction MEF file. The bad pixels will be saved as 1's
``mef``                     Tool to convert from SEF to MEF and viceversa; also allows to give splits of the extensions or join SEFs.
``collapse``                Collapses (add them up arithmetically) each cube of a list files into a single 2D image.
``genLogsheet``             Generates a text file as a log sheet from a set of images.
``imtrim``                  Crops/cuts the input image edges
``modFITS``                 Allows to perfom the modification of any FITS keyword
``runStarfocus``            Run IRAF.starfocus for a focus sequece and return the best focus value and a plot of the fit.
``runPsfmeasure``           Run IRAF.psfmeasure for a focus sequece and get field FWHM of given stars
``getDarks``                Gives the unique values of [read_mode, itime, ncoadd, save_mode] of a set of files of a given directory. 
                            Used to know the DARKS required from them.
``getImageOffsets``         Gives the image offsets (arcsecs) based on the WCS of the image headers
=======================     ===========



``papi``
********

**Description:**

The papi module (see :ref:`papi`) is the main PAPI module to run the data reduction.
It starts by creating a subdirectory in the ``output_dir`` directory using the 
name give on the command line or in the $PAPI_CONFIG file.  Within the run directory 
a [Q1-Q4] subdirectories, one for each detector, will be created. The temporal files
will be saved (and deleted at the end) in the ``temp_dir`` directory.


**Syntax:**

::

    Usage: papi.py [OPTION]... DIRECTORY...
    
    This is PAPI, the PANIC PIpeline data reduction system - IAA-CSIC - Version 1.2.20150508064845

    Options:
    --version             show program's version number and exit
    -h, --help            show this help message and exit
    -c CONFIG_FILE, --config=CONFIG_FILE
                            Config file for the PANIC Pipeline application.If not
                            specified, './config_files/papi.cfg' is used.
    -s SOURCE, --source=SOURCE
                            Source file list of data frames. It can be a fileor
                            directory name.
    -d OUTPUT_DIR, --out_dir=OUTPUT_DIR
                            Output dir for product files
    -o OUTPUT_FILE, --output_file=OUTPUT_FILE
                            Final reduced output image
    -t TEMP_DIR, --temp_dir=TEMP_DIR
                            Directory for temporal files
    -r ROWS, --rows=ROWS  Use _only_ files of the source file-list in the
                            rangeof rows specified (0 to N, both included)
    -R, --recursive       Does recursive search for files in source directory
    -l, --list            Generate a list with all the source files read fromthe
                            source and sorted by MJD
    -M REDUCTION_MODE, --red_mode=REDUCTION_MODE
                            Mode of data reduction to do (quick|science|lab|lemon
                            |quick-lemon).
    -m OBS_MODE, --obs_mode=OBS_MODE
                            Observing mode (dither|ext_dither|other)
    -S SEQ_TO_REDUCE, --seq_to_reduce=SEQ_TO_REDUCE
                            Sequence number to reduce. By default, all sequences
                            found will be reduced.
    -W DETECTOR, --window_detector=DETECTOR
                            Specify which detector to process:Q1(SG1), Q2(SG2),
                            Q3(SG3), Q4(SG4), Q123(all except SG4), all [default:
                            all]
    -p, --print           Print all detected sequences in the Data Set
    -T SEQ_TYPE, --sequences_type=SEQ_TYPE
                            Specify the type of sequences to show: DARK,
                            FLAT(all), DOME_FLAT, SKY_FLAT, FOCUS, SCIENCE, CAL,
                            all [default: all]
    -b, --build_calibrations
                            Build all the master calibrations files
    -C EXT_CALIBRATION_DB, --ext_calibration_db=EXT_CALIBRATION_DB
                            External calibration directory (library of Dark & Flat
                            calibrations)
    -D MASTER_DARK, --master_dark=MASTER_DARK
                            Master dark to subtract
    -F MASTER_FLAT, --master_flat=MASTER_FLAT
                            Master flat to divide by
    -B BPM_FILE, --bpm_file=BPM_FILE
                            Bad pixel mask file
    -g GROUP_BY, --group_by=GROUP_BY
                            kind of data grouping (based on) to do with thedataset
                            files (ot |filter)
    -k, --check_data      if true, check data properties matching (type, expt,
                            filter, ncoadd, mjd)
    -e, --Check           Check if versions of PAPI modules are right.


PAPI creates a in-memory SQLite_ database to store the uncalibrated input data fits 
headers and pipeline metadata. 

**Results:**

FITS file/s with coadd as result of the reduction and calibration of the specified sequences; otherwise,
the error will be shown in the console and log file.


**Examples:**

The following example reduce, in quick mode, all the sequences of the given directory:

::
   
   $papi.py -s /my/raw_data/directory -d /my/output/directory -M quick

   
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
NOT normalized and previously created. 
The flatfield will be normalized to make a gainmap and set bad pixels to 0.

Usage::

Options:
  -h, --help            show this help message and exit
  -s SOURCE_FILE, --source=SOURCE_FILE
                        Flat Field image NOT normalized. It has to be a fullpath file name (required)
  -o OUTPUT_FILENAME, --output=OUTPUT_FILENAME
                        output file to write the Gain Map
  -L MINGAIN, --low=MINGAIN
                        pixel below this gain value  are considered bad (default=0.5)
  -H MAXGAIN, --high=MAXGAIN
                        pixel above this gain value are considered bad (default=1.5)
  -x NXBLOCK, --nx=NXBLOCK
                        X dimen. (pixels) to compute local bkg (even) (default=16)
  -y NYBLOCK, --ny=NYBLOCK
                        Y dimen. (pixels) to compute local bkg (even) (default=16)
  -n NSIGMA, --nsigma=NSIGMA
                        number of (+|-)stddev from local bkg to be bad pixel (default=5)
  -N, --normal          if true, the input flat-field will be normalized before build the gainmap (default=True)


Example::

    $ calGainMap.py -s /tmp/masterTF.fits -o /tmp/masterGain.fits
    $ calGainMap.py -s /tmp/masterTF.fits -o /tmp/masterGain.fits -L 0.7 -H 1.2
    
    
.. index:: flat-field, super-flat 


``calNonLinearity``
*******************
In the moment of writing this manual is not known if PANIC detectors (HAWAII-2RG_) 
will need a non-linearity correction, thus that procedure is not completed yet.
During the commissioning phase of the instrument will be decided if the correction
is needed and what is the correction to apply.  

The ``calNonLinearity`` module corrects PANIC images for their count-rate dependent 
non-linearity. It used the header keywords READMODE and FILTER to determine the 
non-linearity parameter. It corrects the first image, and in the case of a 
multi-extension image, the second image as well, with the appropriate power law. 
For details see `Correcting the PANIC count-rate 
dependent non-linearity <http://www.iaa.es/PANIC/papi/documents/nonlinearity.pdf>`_

.. index:: mask, masking, applymask, ds9


``dxtalk``               
**********
In the moment of writing this manual is not known if PANIC detectors (HAWAII-2RG_) 
will have a crosstalk effect between the data transfer lines of the different
channels. However, because PAPI can be also used to reduce Omega2000 images which show
this effect and in order to be ready to remove this crosstalk effect if it appears
in PANIC, a de-crosstalk routine has been included in the pipeline. It can be activated 
or deactivated in the :ref:`config` (remove_crosstalk=True|False).

During the commissioning phase of the instrument will be checked if there is any
crosstalk effect and this routine could be debugged and tuned.



"Characterization, Testing and Operation of Omega2000 Wide Field Infrared
Camera", Zoltan Kovacs et.al.

Although bright stars can saturate the detector, resetting of the full array
prevents this excess in the pixel values from causing any residual image 
effects in the following image of the dithering. Nevertheless, the satured
pixels generate a crosstalk between the data transfer lines of the different
channels of the quadrant in which they are situated. The data lines of the 
channels are organized in parallel and there might be an interference between 
the data lines transferring the high video signal and the neighbour ones. As a 
result of this crosstalk, a series of spots with the distances of 128 pixels 
from each other appeares in the whole quadrant, corresponding to each channel. 
The average values of the spots were lower than the background signal and their
difference was few percent, which is large enough to degrade the photometric
correctness at the place they are situated. These spots could not be measured
in the raw images but they were well discernible in the reduced frames (Fig. 9). 
This effect was a general feature of the operation of all the  HAWAII-2 detectors 
we tested and should be considered for the choice of pointing positions in any 
field of next observations.  

Usage::

    Options:
      -h, --help            show this help message and exit
      -i INPUT_IMAGE, --input_image=INPUT_IMAGE
                            input image to remove crosstalk
      -o OUTPUT_IMAGE, --output=OUTPUT_IMAGE
                            output filename (default = dxtalk.fits)
      -O, --overwrite       overwrite the original image with the corrected one

Example::
    
    $ ./dxtalk.py -i /tmp/pruebaDC.fits -O
    $ ./dxtalk.py -i /tmp/pruebaDC.fits -o /tmp/pruebaDC_dx.fits
    
``makeobjmask``          
***************
Creates object masks (SExtractor_ OBJECTS images) for a list of FITS images or a 
single FITS image.
Expects the command "sex" (SExtractor Version 2+) in path.  If weight maps
exist they will be used (assume weight map filename given by replacing .fits
with .weight.fits).

The module can produce single poing masks,i.e, a single pixel set to 1 per each
detected object if `single_poing` option is true.

Usage::

    Options:
      -h, --help            show this help message and exit
      -s INPUTFILE, --file=INPUTFILE
                            It can be a source file listing data frames or a single FITS file to process.
      -o OUTPUTFILE, --output=OUTPUTFILE
                            Output text file including the list of objects mask files created by SExtractor ending with '.objs' suffix
      -m MINAREA, --minarea=MINAREA
                            SExtractor DETECT_MINAREA (default=5)
      -t THRESHOLD, --threshold=THRESHOLD
                            SExtractor DETECT_THRESH (default=2.0)
      -l SATURLEVEL, --saturlevel=SATURLEVEL
                            SExtractor SATUR_LEVEL (default=300000)
      -1, --single_point    Create a single point object mask (default=False)



  
Example::
    $ ./makeobjmask.py -s /tmp/reduced_SEQ.fits -o /tmp/obj_mask.txt
    $ ./makeobjmask.py -s /tmp/reduced_SEQ.fits -o /tmp/obj_mask.txt -1 -l 100000 -m 10


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

``correctNonLinearity``
***********************
Performs the non-linearity correction of the PANIC raw data files using the 
proper NL-Model (FITS file).

Usage::
  
  Options:
    -h, --help            show this help message and exit
    -m MODEL, --model=MODEL
                          FITS MEF-cube file of polinomial coeffs (c4, c3, c2, c1).
    -s SOURCE_FILE_LIST, --source=SOURCE_FILE_LIST
                          Source file list of FITS files to be corrected.
    -o OUT_DIR, --out_dir=OUT_DIR
                          filename of out data file (default=/tmp)
    -S SUFFIX, --suffix=SUFFIX
                          Suffix to use for new corrected files.
    -f, --force           Force Non-linearity correction with no check of headervalues (NCOADD, DATE-OBS, DETROT90, ...



``solveAstrometry``
*******************
Performs the astrometric calibration of a set of images, in principle previously 
reduced, but not mandatory; Astromety.net tool is used.

Options:

  -h, --help            show this help message and exit
  -s SOURCE_FILE, --source=SOURCE_FILE
                        Source file list of data frames. It can be a file or directory name.
  -o OUTPUT_DIR, --output_dir=OUTPUT_DIR
                        Place all output files in the specified directory [default=/tmp]
  -p PIXEL_SCALE, --pixel_scale=PIXEL_SCALE
                        Pixel scale of the images
  -r, --recursive       Recursive subdirectories (only first level)


``remove_cosmics``
******************
Remove the cosmic ray hits in the input image.

Options:

  -h, --help            show this help message and exit
  -i INPUT_IMAGE, --input_image=INPUT_IMAGE
                        input image to remove cosmics
  -o OUTPUT_IMAGE, --output=OUTPUT_IMAGE
                        output filename (default = without_cosmics.fits)
  -O, --overwrite       overwrite the original image with the corrected one
  -m, --mask            If true, the mask with cosmics detected and removed is written into a FITS file.

.. _astromatic: http://www.astromatic.net/
.. _SExtractor: http://www.astromatic.net/software/sextractor
.. _scamp: http://www.astromatic.net/software/scamp
.. _swarp: http://www.astromatic.net/software/swarp
.. _SQLite: http://www.sqlite.org
.. _HAWAII-2RG: http://panic.iaa.es/detectors

