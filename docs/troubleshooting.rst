***************
Troubleshooting
***************

This section gives a description of some problems that can be obtained and how 
they can be solved.

PAPI does not recognize a group of files as a well formed sequence.
===================================================================

It maybe because the sequence is missing some files of the sequence, or some 
of them have no proper headers (ie. they were not observed with OT).


Solutions
---------
* You can try to reduce the sequence using the groping mode 'filter' instead 
  of 'ot'. It means PAPI will not use any OT header information, but RA, DEC, 
  DATE-OBS and FILTER.

* Use the 'build_sequence' tool which will create a new sequence adding the
  required header keywords.  

PAPI can not run the astrometric calibration
============================================

Verify that:
------------
* You have cdsclient installed
* You have an internet connection and no firewall is blocking CDSClient
* You can config CDSClient to run with a proxy.

PAPI can not run the photometric calibration
============================================

Verify that:
------------

* You have an Internet connection. It is used to query the on-line 2MASS or 
  USNO-B catalog.
* Check that the input reduced image is astrometrically calibrated and it the 
  filter name FILTER keyword match some of the 2MASS or USNO-B catalog.


PAPI says some files in the input file list does not look a FITS file.
======================================================================

Verify that:
------------
* The input file is not ending with a blank/empy line.
* The file has Unix text format. Text files created on DOS/Windows machines have 
  different line endings than files created on Unix/Linux. DOS uses carriage 
  return and line feed ("\r\n") as a line ending, which Unix uses just line feed ("\n"). 


Problem 4 ...
=============

.. _astromatic: http://www.astromatic.net/
.. _sextractor: http://www.astromatic.net/software/sextractor
.. _scamp: http://www.astromatic.net/software/scamp
.. _swarp: http://www.astromatic.net/software/swarp
.. _HAWAII-2RG: http://w3.iaa.es/PANIC/index.php/gb/workpackages/detectors

