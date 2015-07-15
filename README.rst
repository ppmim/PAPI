|logo| PAPI
===========

PAPI (PANIC Pipeline) is the automatic image processing pipeline for data obtained 
with the PAnoramic Near Infrared Camera (PANIC_) for the 2.2m and 3.5m Telescopes at 
Calar Alto Observatory (CAHA_). The pipeline is written in Python and developed 
at the `Institute of Astrophysics of Andalusia (CSIC) <http://www.iaa.es/>`_. 
The automated processing steps include basic calibration (removeing instrumental 
signature), cosmic-ray removal, treatment for electronic ghosts (cross-talk), 
sky subtraction, non-linear count-rate correction, robust alignment and 
registration.


PANIC_ is a general purpose Panoramic Near Infrared camera for Calar Alto. 
It is optimized for use at the 2.2m telescope, but can also be installed 
at the 3.5m telescope. It works in the nIR bands Z, J, H and K. 



Installation
============

PAPI has the following strict requirements:
 
 - Python 2.7
 - Numpy 1.6 or later 

and also depends on next packages:

 - `SciPy <http://www.scipy.org>`_ (> v0.12.2)
 - `Astropy <http://www.astropy.org/>`_ (> v0.3.1)
 - `Matplotlib <http://matplotlib.org/>`_ (> v1.3.0)
 - `PyQt4 <http://www.riverbankcomputing.co.uk/software/pyqt/download>`_
 - `IRAF <http://iraf.noao.edu/>`_ with STSDAS and MSCRED (< v2.16 or higher)
 - `stsci_python <http://www.stsci.edu/resources/software_hardware/pyraf/stsci_python>`_ (> v2.14)
 - `CDSClient <http://cdsarc.u-strasbg.fr/doc/cdsclient.html>`_
 - `SExtractor <http://astromatic.iap.fr/software/sextractor/>`_ (> v2.8.6)
 - `SCAMP <http://www.astromatic.net/software/scamp>`_ (> v1.7.0)
 - `SWarp <http://www.astromatic.net/software/swarp>`_ (> v2.19.1)
 - `Astrometry.net <http://astrometry.net/>`_ with 42xx index files (optional to SCAMP)
 - `SAO DS9 and XPA <http://hea-www.harvard.edu/RD/ds9>`_ (> v7.3b5)
 - `Montage <http://montage.ipac.caltech.edu/download/Montage_v3.3.tar.gz>`_ (v3.3)
 


Note that, for PyRAF_ you have to install IRAF_(v2.16 or later), what can be a 
tricky task. However, is has been simplified in recent versions.


To install PAPI, follow the next steps:

1. Clone the PAPI files in a suitable location. Note that, it is a development 
version:

	* ``git clone https://github.com/ppmim/PAPI ~/papi``

#. Build the sources:

    * cd papi/irdr/src
    * make all

#. Edit papi_setup.sh script:

    * Modify the PAPI_HOME and PAPI_PROD in the papi_setup.[sh] file in the papi 
    directory
    * Run the papi_setup.sh 
    * Re-load your new profile (.bashrc or .cshrc ) 

        - Bash: . ~/.bashrc
        - CSH: source ~/.cshrc

#. Go to config_files/ directory to setup the config file to use.


Supported Platforms
===================
Currently PAPI has only be tested under openSuSE12.x and openSuSE13.1, but it
should work on any 64-bit Linux box with the software packages required above.


Documentation
=============
You can browse the latest release documentation_ online.



Webpage: http://www.iaa.es/PANIC
Maintainer: jmiguel@iaa.es


.. links:
.. |logo| image:: ./QL4/images/logo_PANIC_100.jpg
          :width: 127 px
          :alt: PANIC icon

.. _PANIC: http://www.iaa.es/PANIC
.. _CAHA: http://www.caha.es
.. _iaa_web: http://www.iaa.es
.. _mpia_web: http://www.mpia.de
.. _source code: http://github.com/ppmim/PAPI
.. _documentation: http://www.iaa.es/~jmiguel/PANIC/PAPI/html/index.html
.. _SciPy: http://www.scipy.org
.. _PyFITS: http://www.stsci.edu/resources/software_hardware/pyfits
.. _PyRAF: http://www.stsci.edu/institute/software_hardware/pyraf
.. _PyQt4: http://www.riverbankcomputing.co.uk/software/pyqt/download
.. _Astropy: http://www.astropy.org/
.. _Astrometry.net: http://astrometry.net/
.. _Astromatic: http://www.astromatic.net/
.. _Sphinx: http://sphinx-doc.org/
.. _IRAF: http://www.iraf.net
