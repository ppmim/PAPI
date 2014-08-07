
.. _installation:

Installation & Configuration  
============================

PAPI is the automatic image processing pipeline for data taken with the 
`PAnoramic Near Infrared Camera (PANIC) <http://www.iaa.es/PANIC>`_ for the 2.2m 
and 3.5m Telescopes at `Calar Alto Observatory (CAHA) <http://www.caha.es>`_. 
The pipeline is written in Python and developed at the `Institute of Astrophysics 
of Andalusia (CSIC) <http://www.iaa.es/>`_. The automated processing steps 
include basic calibration (removeing instrumental signature), cosmic-ray removal, 
treatment for electronic ghosts (cross-talk), sky subtraction, non-linear 
count-rate correction, robust alignment and registration for 
large mosaics. 


This manual is a complete description of the data reduction recipes implemented 
by the PAPI pipeline, showing the status of the current pipeline version
and describing data reduction process of the PANIC data using PAPI.

Although PAPI was developed for the PANIC camera, it also works with data from 
Omega2000_ NIR camera at the 3.5m CAHA telescope and HAWK-I_ camera at the VLT.

In addition to this html version of the manual, there is also a pdf_ version to download.


**Development Team:** José M. Ibáñez (IAA-CSIC)

Caveats
*******

Currently PAPI it is able to reduce data taken with the Observing Tool (OT) 
defining the required observing blocks (OB), or manually through GEIRS scripts.
PAPI was primarily developed and optimized for reducing broad-band imaging data 
of extragalactic sources (such as imaging data taken for field galaxy surveys and 
galaxy cluster surveys). Other types imaging data have been reduced with PAPI 
but results can not be as good as desired. (See :ref:`troubleshooting` for tips).
PAPI is **not** designed to reduce any kind of field taken with PANIC.  

.. index:: prerequisites, requirements, pipeline


Prerequisites
*************

`Python 2.7 <http://www.python.org>`_ is required. In addition, PAPI requires the following packages 
installed in order to work properly:

    * `NumPy <http://numpy.scipy.org/>`_ 
    * `SciPy <http://www.scipy.org>`_
    * `IRAF <http://iraf.noao.edu/>`_
    * `STSDAS/TABLES <http://www.stsci.edu/institute/software_hardware/stsdas/download-stsdas/>`_
    * `stsci_python <http://www.stsci.edu/resources/software_hardware/pyraf/stsci_python>`_ (> v2.2)
    * `CDSClient <http://cdsarc.u-strasbg.fr/doc/cdsclient.html>`_
    * `SExtractor <http://astromatic.iap.fr/software/sextractor/>`_ (> v2.3.2)
    * `SCAMP <http://www.astromatic.net/software/scamp>`_
    * `SWarp <http://www.astromatic.net/software/swarp>`_
    * `SAO DS9 and XPA <http://hea-www.harvard.edu/RD/ds9>`_

Additional packages are optionally required:
    * `sphinx`_  to build the documentation

.. index:: installing, building, source, downloading

Download
********
The latest stable version of PAPI can be downloaded from `here <http://www.iaa.es/~jmiguel/software/papi.tgz>`_

Installing
**********

To install PAPI, use the standard installation procedure:::

    $ tar zxvf papi-X.Y.Z.tar.gz
    $ cd papi-X.Y.Z
    $ python setup.py install


Edit ``papi_setup`` script

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

Development version
*******************

The development version can be checked out with:::

    $ git clone https://github.com/ppmim/PAPI.git

And then installed following the next procedure:::

    $ cd papi
    $ cd irdr 
    $ make all
    

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
  