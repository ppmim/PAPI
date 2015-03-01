++++++++
PAPI FAQ
++++++++

------------
Installation
------------

Where can I get PAPI ?
----------------------
The stable version can be dowloaded from https://github.com/ppmim/PAPI/.

Does PAPI work with python 3?
-----------------------------
It has not been tested yet, but probably it does not work.

What GUI toolkit does PAPI-QL use?
----------------------------------
PAPI Quick-Look uses Qt toolkit.

Can PAPI-QL work with PyQt5?
--------------------------
It has not been tested yet, but probably it does not work.

---------------
Data reduction
---------------

Are there any known problems with PANIC data?
---------------------------------------------
Please, go to http://www.iaa.es/PANIC to get the list of last known problems.


Where can I find tutorial on PANIC data reduction?
--------------------------------------------------
Please, visit http://www.iaa.es/~jmiguel/PANIC/PAPI/html/index.html for a detaild
description.


I have a problem with the PAPI pipeline. What should I do?
----------------------------------------------------------
The PANIC pipeline is continuosly under development; please, follow the changes
and publish issues on https://github.com/ppmim/PAPI/issues.


What is the best way to reduce PAPI data?
-----------------------------------------
The recommend to use the OT and execute the OBs. That way the headers will include
meta-data about the observation, and thus the pipeline can group the data and
find the required calibrations for a successful reduction. Then, you only have to
type:

> papi -s /data1/PANIC/my_program/ 


How can we reach hundredth of magnitude accuracy in photometry ?
----------------------------------------------------------------
The best way to accurately photometrically calibrate PANIC images is to use 2MASS 
stars in the field itself to derive the photometric solution. The accuracy 
strongly depends on the number of bright 2MASS stars within the filed of view, 
but ranges from a few 1/100th of a magnitude to 0.1 magnitudes if only faint 
stars are contained in the field. Additionally, observing supplementary standard
star fields can be asked for when preparing the observations. To perform the 2MASS 
photometric calibration on an image you should use the 'photometry' command as 
follow:

>  photometry -i /directory/prereducedField.fits -o test.pdf


How good is PAPI astrometry and how are PSF variations corrected ?
------------------------------------------------------------------
At present the pipeline applies a correction for PSF distortions based on a 
distortion map derived during the astrometric calibration done with SCAMP (a 
software developped by Emmanuel Bertin) and 2MASS.
