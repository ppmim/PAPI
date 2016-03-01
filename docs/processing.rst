.. _processing:

**********************
Processing description
**********************
This section gives a description of each step of the pipeline in a greater 
detail and list the parameters that can be changed if needed.

Next figure shows the main steps that are involved in the PANIC pipeline:


.. image:: _static/PAPI_flowchart.jpg
   :align: center
   :scale: 90%

Outline
-------

    * the non-linearity is corrected;
    * a master flat-field in computed by combining all the appropriate 
      images without offsets; if the mosaic is of an Sky-Target type, 
      only the sky frames are used;
    * a bad pixel mask, is computed by detecting, in the master flat 
      field, the pixels with deviant values;
    * if provided, an external bad pixel mask is also used, adding the 
      bad pixel in the previus one;
    * for each image, a sky frame is computed by combining a certain 
      number of  the closest images; 
    * this sky frame is subtracted by the image and the result is 
      divided by the master flat;
    * bright objects are detected by SExtractor_ in these cleaned images 
      to measure the offsets among the images; the object mask are 
      multiplied by the bad pixel mask to remove false detections;
    * a cross-correlation algorithm among the object masks is used to 
      measure the relative offsets. It also works if no object is 
      common to all the images; 
    * the cleaned images are combined using the offsets, creating the 
      "quick" image;
    * to remove the effect of faint obejcts on the estimate of the sky 
      frames, SExtractor_ is used on the combined image to create a master 
      object mask;
    * the object mask is dilatated by a certain factor to remove also 
      the undetected object tails;
    * for each image a new sky is computed by taking into account 
      this object mask;
    * if field distortion can be neglected, these images are combined 
      by using the old offsets, creating the "science" image;
    * field distortion is removed from the cleaned images by using 
      SCAMP computed distortion model
    * the pixels containing deviant pixels are identified and flagged;
    * the old offsets could be effected by field distortion, therefore 
      new offsets are computed for the undistorted images;
    * finally, the cleaned corrected images are combined.

Main configuration file
***********************
See :ref:`Main config file <config>`


Data-set classification
***********************

One of the main features of PAPI is that the software is able to do an automatic
data reduction. While most of the pipelines are run interactively, PAPI is able
to run without human interaction. It is done because of the classificaton algorithm
that is implemented in PAPI and that allow an automatic identification of the 
data sets grouping the files according to the observation definition with the OT.

1 - The data grouping algorithm
2 - Sky finding algorithm for extended objects


In case of not using the OT during the observation, also a data grouping is possible,
althouth with some limitations. Let's see how it works:

[...]

Data Preparation
****************
Firstly, each FITS file is linearity corrected if it was enabled in the configuration 
file (nonlinearity:apply). If integrations where done with repetitions >1 and saved as
a cube with N-planes, then the FITS cube is collapsed doing a simple arithmetich sum of
N-planes.

Then the image is divided into the number of chips in the FPA (which constitutes 4 chips 
in a mosaic). From this step on, the pipeline works on individual chips rather than whole 
images, thereby enhancing the speed and enabling us to do multi-chip processing on multi CPUs.


Calibrations
************
In next sections we describe the main calibration to be done by PAPI.

Computing the master dark
-------------------------
TBC

Computing the master flat-field
-------------------------------
TBC

Computing the Bad Pixel Mask
----------------------------

The map of all bad pixels (hot, low QE) are derived from the non-linearity tests. However, also
the nonlinearity analysis provides a list of non-correctable pixels, which always will be
considered invalid. 

So, currently there is no procedure in PAPI to compute the right bad pixel mask (BPM).



First pass sky subtraction
**************************

Sky model
---------
TBC

Object detection
****************

Offset computation
******************

First pass coaddition
*********************

Master object mask
******************
SExtractor_ is again used to find objects in this first-pass coadded image in 
order to mask then during next sky estimation. This time the parameters controlling
the detection threshold should be set to have deeper detections and mask faint
objects. The parameters involved nad ther default values are:

mask_minarear = 10
mask_thresh = 1.5

The resulting object mask is extended by a certain fraction to reject also 
the undetected object tails. 


Non-Linearity
*************

HAWAII-2RG near-IR detectors exhibit an inherent non-linear response. 
It is caused by the change of the applied reverse bias voltage due to the 
accumulation of generated charge.
The effect increases with signal levels, so that the measured signal deviates stronger 
from the incident photon number at higher levels, and eventually levels out when 
the  pixel well reaches saturation.

The  common  approach  is  to  extrapolate  the  true  signal Si(t) from measurements
with low values, and fit it as a function of the measured data S(t) with a polynomial of 
order n:


For the correction, PAPI uses a master Non-Linearity FITS file that store the fit to be
applied to the raw images. There is file for each readout mode. The filename is composed
as::

    mNONLIN_<readmode>_<version>.fits

The FITS file has a primary header with no data, and two data extensions for each detector.
They are labeled LINMAX<i> and LINPOLY<i> with i=1...4 being the quadrant index, numbered
similar to the scheme for MEF data files from GEIRS. Note that the indices do not
necessarily correspond to SG hardware IDs, which are written in the header instead.

The extension LINMAX<i> is a 32bit float 2048x2048 data array containing the maximum
correctable signal for each detector. Uncorrectable pixels have a NaN instead of a 
numerical value.
The extension and LINPOLY<i> is a 32bit float 2048x2048x4 data cube containing the
polynomial coefficients c[1...4] in reverse order. The first slice in the cube 
is c[4], the second c[3], etc.

The module used to correct the non-linearity is ``correctNonLinearity.py``; in adition
the non-linearity correction can be enable in the configuration file $PAPI_CONFIG setting
in the *nonlinearity* section the keyword *apply = True*.




Crosstalk
*********

HAWAII2 sensors with multiple parallel readout sections can show crosstalk 
in form of compact positive and negative ghost images whose amplitude varies between 
readout sections. PAPI has a optional de-crosstalk module that assumes that the 
amplitude is the same, therefore the correction will only partially remove the 
effect (if at all). If you know in advance that this will be a problem for your 
science case, then consider choosing different camera rotator angles for your 
observations.


The first effort at characterizing and removing the cross-talks made use of 
the "Medamp" technique. By this we mean isolating then subtracting what is 
common to all 32 amplifiers. This effectively seems to remove the edge and 
negative cross-talks which both affect all 32 amplifiers. But it does not 
remove the positive crosstalk. Note that the assumption is that the amplitude 
of the edge and negative cross-talks is the same ona ll 32 channels. We tried 
inconclusively to prove/disprove that assumption. If amplifier-dependant, the 
amplitude variations must be less than 10%.

We experimented doing the medamp at various stages of the processing and found 
the best results when removing the crosstalk as the very last step, after sky 
subtraction. Rigorously, it should actually be the very first step since 
crosstalk effects are produced in the very last stages of image generation.

The module used to correct the crosstalk is ``dxtalk.py.py``; in adition
the crosstalk correction can be enable in the configuration file $PAPI_CONFIG setting
in the *general* section the keyword *remove_crosstalk = True*.




Extended Objects
****************
If your targets are really extended and/or very faint, then you should seriously 
consider observing blank SKY fields. They will be recognized and automatically 
used in the correct manner once identified by PAPI. No additional settings 
have to be made. You should check though that the images have correct header keys.


.. _astromatic: http://www.astromatic.net/
.. _SExtractor: http://www.astromatic.net/software/sextractor
.. _scamp: http://www.astromatic.net/software/scamp
.. _swarp: http://www.astromatic.net/software/swarp
.. _HAWAII-2RG: http://w3.iaa.es/PANIC/index.php/gb/workpackages/detectors

