
.. _installation:

Installation & Configuration  
============================

.. index:: prerequisites, requirements, pipeline


Prerequisites
*************

Because PAPI is written mostly in Python, `Python 2.7.x <http://www.python.org>`_  
or  higher is required. In addition, PAPI requires the following packages 
installed in order to work properly:

    * `NumPy <http://numpy.scipy.org/>`_ (> v1.6)
    * `SciPy <http://www.scipy.org>`_ (> v0.12.2)
    * `Astropy <http://www.astropy.org/>`_ (> v0.3.1)
    * `Matplotlib <http://matplotlib.org/>`_ (> v1.3.0)
    * `PyQt4 <http://www.riverbankcomputing.co.uk/software/pyqt/download>`_
    * `IRAF <http://iraf.noao.edu/>`_ with STSDAS and MSCRED (< v2.16 or higher)
    * `stsci_python <http://www.stsci.edu/resources/software_hardware/pyraf/stsci_python>`_ (> v2.14)
    * `CDSClient <http://cdsarc.u-strasbg.fr/doc/cdsclient.html>`_
    * `SExtractor <http://astromatic.iap.fr/software/sextractor/>`_ (> v2.8.6)
    * `SCAMP <http://www.astromatic.net/software/scamp>`_ (> v1.7.0)
    * `SWarp <http://www.astromatic.net/software/swarp>`_ (> v2.19.1)
    * `Astrometry.net <http://astrometry.net/>`_ with 42xx index files (optional to SCAMP)
    * `SAO DS9 and XPA <http://hea-www.harvard.edu/RD/ds9>`_ (> v7.3b5)

Additional packages are optionally required:
    * `sphinx`_  to build the documentation

.. index:: installing, building, source, downloading

Download
********
The latest stable version of PAPI can be downloaded from `here <https://github.com/ppmim/PAPI>`_

Installing
**********
Once you have installed the packages described above, you are ready to install
PAPI; for this, follow the next steps::

    $ git clone https://github.com/ppmim/PAPI.git ~/papi
    $ cd papi
    $ ./papi_setup.sh


 * Modify the PAPI_HOME in the papi_setup.[sh][.csh] file in the PAPI bin directory
 * Add the papi_setup.[sh][.csh] to your .bash_profile or .cshrc (.tcshrc)

    	* Bash: ``. $PAPI_HOME/bin/papi_setup.sh``
    	* CSH: ``source $PAPI_HOME/bin/papi_setup.csh``


Example papi_setup.sh::
	
    #---------------------------------------------------------------------------
    # User Configurable Settings
    #---------------------------------------------------------------------------

    # path to PAPI directory
    export PAPI_HOME=${HOME}/pipelines/papi

    # path to PAPI output data products
    export PAPI_PROD=${HOME}/DataProd

    #---------------------------------------------------------------------------
    # Fixed Settings
    #---------------------------------------------------------------------------
    # path to PAPI reference files
    export PAPI_CONFIG=${PAPI_HOME}/config_files
    export PATH=${PATH}:${PAPI_HOME}/bin
    export PYTHONPATH=${PYTHONPATH}:${PAPI_HOME}/lib


Building the documentation
**************************
The PAPI documentation is base on `sphinx`_. With the package installed, the 
html documentation can be built from the `doc` directory::

  $ cd doc
  $ make html
  
The documentation will be copied to a directory under `build/sphinx`.
  
The documentation can be built in different formats. The complete list will appear
if you type `make`.
 

.. _PANIC: http://www.iaa.es/PANIC
.. _CAHA: http://www.caha.es
.. _Omega2000: http://www.caha.es/CAHA/Instruments/O2000/index.html
.. _HAWK-I: http://www.eso.org/sci/facilities/paranal/instruments/hawki/
.. _sphinx: http://sphinx.pocoo.org
.. _pdf: http://www.iaa.es/~jmiguel/PANIC/PAPI/PAPI.pdf
  