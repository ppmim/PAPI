
.. _installation:

Installation & Configuration  
****************************

.. index:: prerequisites, requirements, pipeline


Requirements and Supported Platforms
------------------------------------

Because PAPI is written mostly in Python_ and ANSI C, it can run on any platform
that has the required Python modules and GCC compilier. However, it has been developed
and deeply tested under `openSuSE`_ 12.x/13.x x86_64 Linux OS.  
`Python 2.7.x <http://www.python.org>`_ or higher and the following packages 
are required:

    * `NumPy <http://numpy.scipy.org/>`_ (> v1.6)
    * `SciPy <http://www.scipy.org>`_ (> v0.12.2)
    * `Astropy <http://www.astropy.org/>`_ (> v0.3.1)
    * `Matplotlib <http://matplotlib.org/>`_ (> v1.3.0)
    * `PyQt4 <http://www.riverbankcomputing.co.uk/software/pyqt/download>`_
    * `IRAF <http://iraf.noao.edu/>`_ with STSDAS and MSCRED (v2.16)
    * `x11iraf <http://iraf.noao.edu/iraf/ftp/iraf/x11iraf/x11iraf-v2.0BETA-bin.linux.tar.gz>`_ for xgterm
    * `stsci_python <http://www.stsci.edu/resources/software_hardware/pyraf/stsci_python>`_ (> v2.14)
    * `CDSClient <http://cdsarc.u-strasbg.fr/doc/cdsclient.html>`_
    * `SExtractor <http://astromatic.iap.fr/software/sextractor/>`_ (> v2.8.6)
    * `SCAMP <http://www.astromatic.net/software/scamp>`_ (> v1.7.0)
    * `SWarp <http://www.astromatic.net/software/swarp>`_ (> v2.19.1)
    * `Astrometry.net <http://astrometry.net/>`_ with `42xx index files <http://broiler.astrometry.net/~dstn/4200/>`_
    * `SAO DS9 and XPA <http://hea-www.harvard.edu/RD/ds9>`_ (> v7.3b5)
    * `Montage <http://montage.ipac.caltech.edu/download/Montage_v3.3.tar.gz>`_ (v3.3)
    * `montage_wrapper <https://pypi.python.org/pypi/montage-wrapper>`_ (0.9.8)
 
Additional packages are optionally required:
    * `sphinx`_  to build the documentation

.. note::
    
    If you are using a SCAMP version <= 2.0.4 (lastest stable version), then you need to install the CDSClient. Otherwise, if you are using SCAMP version > 2.0.4, then you need **libcurl**. 

    Anycase, if you are behind a proxy, you need to set the proxy server in your system::
    
    http_proxy=http//your_proxy:your_port; export http_proxy

    
.. index:: installing, building, source, downloading

Download
--------

The latest stable version of PAPI can be downloaded from `GitHub repository <https://github.com/ppmim/PAPI>`_ .

Building and Installation
-------------------------
PAPI installation is thought to be done as a 'personal user' (non-root), however it should work
under any system directory (ie., /usr/local/). 

1. To install PAPI as a "personal user" (non-root), follow the next steps:

Once you have installed the required packages described above, you are ready to install
PAPI; for this, follow the next steps::

    $ git clone https://github.com/ppmim/PAPI.git ~/papi
    $ cd papi
    $ ./papi_setup.sh


2. To install PAPI as root on your system, follow the next steps::

    $ cd /usr/local
    $ git clone https://github.com/ppmim/PAPI.git papi
    $ cd papi
    
    Edit the papi_setup.sh and set the right values to PAPI_HOME and PAPI_BIN variables, and then run the script as an user:
    
    $ ./papi_setup.sh


.. warning::
    
    The script papi_setup.sh is currently implemented **only** for the Bash shell, and will modify your .bashrc file adding a new line at the end.

    

Building the documentation
--------------------------

The PAPI documentation is base on `sphinx`_. With the package installed, the 
html documentation can be built from the `doc` directory::

  $ cd papi/doc
  $ make html
  
The documentation will be copied to a directory under `build/sphinx`.
  
The documentation can be built in different formats. The complete list will appear
if you type `make`.

Bug reports
-----------

Please submit issues with the `issue tracker <https://github.com/ppmim/PAPI/issues>`_ on github.


Release Notes
-------------

* 1.2.x
    - Support for new MEF structure (Qi); old format (SGi_1) also supported
    - Bug Fixes
* 1.0.x
    - First version
    
    
.. _PANIC: http://www.iaa.es/PANIC
.. _CAHA: http://www.caha.es
.. _Omega2000: http://www.caha.es/CAHA/Instruments/O2000/index.html
.. _HAWK-I: http://www.eso.org/sci/facilities/paranal/instruments/hawki/
.. _sphinx: http://sphinx.pocoo.org
.. _pdf: http://www.iaa.es/~jmiguel/PANIC/PAPI/PAPI.pdf
.. _openSuSE: http://www.opensuse.org/
.. _issue tracker
.. _Python: http://www.python.org
