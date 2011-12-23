PANIC Quick-Look Tool (PQL)
===========================

.. index:: quick-look, running


Running on-line
***************

Running PQL can be as simple as executing the following command in a terminal::
	
	runPQL.py config_file

Where ``config_file`` is the name of the configuration file to use.

.. index:: config, quicklook, papi, on-line

Getting PAPI Data
*****************

The PAPI pipeline requires the full set of uncalibrated data products 
and best reference files for each observation in the input image set. These files 
can be readily obtained through the CAHA archive. When
requesting data from CAHA you need to specify:
	
	* Science Files Requested: **Uncalibrated** 
	* Reference Files: **Best Reference Files**

.. image:: _static/caha_archive.jpg
   :align: center
   :height: 300 px
   :width: 565 px

.. CAHA: http://www.caha.es/APPS/ARCHIVE/

.. index:: options

Optional Commands
*****************

For most image sets PAPI can be run in the default configuration with no 
additional interaction required. If the default settings are insufficient for 
processing a particular data set, there are a number of run-time options which 
may be applied to help improve the reductions:

	* Modules can be run manually step-by-step allowing for the inspection of the output at each step.
	* Modules can be skipped.
	* A pipeline run can be stopped, restarted or rerun at any stage of the reduction after the initial setup.
	* An external reference image can be used to improve the internal alignment of the reduced NICMOS frames.


Here's a listing of the PAPI command line options::

   Usage: runQL.py [OPTION]... DIRECTORY...

   This module in the main application for the PANIC Quick Loook (PQL) data
   reduction system
   
   Options:
     --version             show program's version number and exit
     -h, --help            show this help message and exit
     -c CONFIG_FILE, --config=CONFIG_FILE
                           config file for the PANIC Pipeline application. If not
                           specified, './config_files/papi_portatil.cfg' is used
     -v, --verbose         verbose mode [default]
     -s SOURCE, --source=SOURCE
                           Source directory of data frames. It has to be a
                           fullpath file name
     -o OUTPUT_DIR, --output_dir=OUTPUT_DIR
                           output directory to write products
     -t TEMP_DIR, --temp_dir=TEMP_DIR
                           temporary directory to write
      
   
   
	
.. index:: quicklook, off-line

Running off-line
****************

Run PQL in off-line mode means that data were already taken and are in a specific
directory that we wish to inspect in quick way.



.. _troubleshooting:

Troubleshooting
***************

As we stated previously, PAPI was developed primarily for reducing imaging data of extragalactic sources. 
Here are some tips for reducing other types of data:

*Add tips here*