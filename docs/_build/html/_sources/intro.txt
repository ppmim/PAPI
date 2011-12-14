
Introduction
============

PAPI is an automatic image processing pipeline for data taken with the 
PAnoramic Near Infrared Camera (PANIC) on the 2.2m Telescope at Calar Alto Observatory. The pipeline
currently supports imaging data from camera and is written in Python 
and C making it portable across many platforms. The automated processing 
steps include basic calibration (removeing instrumental signature), cosmic-ray 
removal, treatment for post-SAA cosmic ray persistence and electronic ghosts, 
sky subtraction, non-linear count-rate correction, 
artifact masking, robust alignment and registration for large mosaics, 
weight map generation, and drizzling onto a final image mosaic. 

**Development Team:** Jose M. Ibanez (IAA-CSIC)

Caveats
*******

Currently it is only possible for PANIC to reduce data taken with the
Observing Tool (OT), but not manually with GEIRS.
PAPI was primarily developed and optimized for reducing broad-band imaging data of
extragalactic sources (such as imaging data taken for field galaxy surveys and galaxy cluster surveys). 
Other types imaging data have been reduced with PAPI but YMMV (See :ref:`troubleshooting` for tips).
PAPI is **not** designed to reduce any kind of field taken with PANIC.  

.. index:: prerequisites, requirements

Prerequisites
*************

These software must be install for PAPI to run:

	* `python <http://www.python.org>`_ (2.4 or 2.5)
	* `sqlite <http://www.sqlite.org>`_ (v3.0 > if using Python 2.4)
	* `pysqlite <http://initd.org/tracker/pysqlite>`_ (v2.2.0 > if using Python 2.4)
	* `stsci_python <http://www.stsci.edu/resources/software_hardware/pyraf/stsci_python>`_ (v2.2 >)
	* `SExtractor <http://astromatic.iap.fr/software/sextractor/>`_ (v2.3.2 >)
	* `SAO DS9 and XPA <http://hea-www.harvard.edu/RD/ds9>`_ (if applying user defined masks)

.. index:: installing, building, source, downloading

Download
********

Download the papi_latest.tgz file and the papi_reffiles.tgz file.


    * `papi_latest.tgz <http://code.google.com/p/panicdrs/files/papi_latests.tgz>`_
    * `papi_reffiles.tgz <http://code.google.com/p/panicdrs/files/papi_reffiles.tgz>`_


Installing
**********

1. Unzip and untar the PAPI files in a suitable location.


    * tar zxf papi_X.X.X.tgz
    * ln -s papi_X.X.X papi
    * cd papi
    * tar zxf papi_reffiles.tgz


2. Build the source


	* cd src
	* make
	* make install


3. Edit ``papi_setup`` script


    * Modify the PAPI_DIR and PAPI_PIPE in the papi_setup.[sh][.csh] file in the PAPI bin directory
    * Add the papi_setup.[sh][.csh] to your .bash_profile or .cshrc (.tcshrc)

    	* Bash: ``. $PAPI_DIR/bin/papi_setup.sh``
    	* CSH: ``source $PAPI_DIR/bin/papi_setup.csh``


Example papi_setup.sh::
	
	#------------------------------------------------------------------------------
	# User Configurable Settings
	#------------------------------------------------------------------------------

	# path to PAPI directory
	export PAPI_DIR=${HOME}/pipelines/papi

	# path to PAPI output data products
	export PAPI_PIPE=${HOME}/Data

	#------------------------------------------------------------------------------
	# Fixed Settings
	#------------------------------------------------------------------------------
	# path to PAPI reference files
	export PAPI_REF=${PAPI_DIR}/PAPI_REF
	export PATH=${PATH}:${PAPI_DIR}/bin
	export PYTHONPATH=${PYTHONPATH}:${PAPI_DIR}/lib
