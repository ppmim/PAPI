Running PAPI
============

.. index:: quickstart, running

Quickstart
**********

Running PAPI can be as simple as executing the following command in a terminal::
	
	papi.py run_name input_data_directory

Where ``run_name`` is the name of the dataset (e.g. A1689_NIC3) and 
``input_data_directory`` is the path to the uncalibrated data directory.

.. index:: uncalibrated, data, mast

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

   Usage: papi.py [OPTION]... DIRECTORY...

   This module in the main application for the PANIC data reduction system
   
   Options:
     --version             show program's version number and exit
     -h, --help            show this help message and exit
     -c CONFIG_FILE, --config=CONFIG_FILE
                           config file for the PANIC Pipeline application.
                           If not specified, './config_files/papi_portatil.cfg'
                           is used
     -s SOURCE, --source=SOURCE
                           Source file list of data frames.                   It
                           can be a file or directory name.
     -o OUTPUT_FILENAME, --output_file=OUTPUT_FILENAME
                           final reduced output image
     -t TEMP_DIR, --temp_dir=TEMP_DIR
                           directory for temporal files
     -d OUTPUT_DIR, --out_dir=OUTPUT_DIR
                           output dir for product files
     -r REDUCTION_MODE, --red_mode=REDUCTION_MODE
                           Mode of data reduction to do (quick|science)
     -m OBS_MODE, --obs_mode=OBS_MODE
                           observing mode (dither|ext_dither|other)
     -p, --print           print detected sequences in the Data Set
     -S SEQ_TO_REDUCE, --seq_to_reduce=SEQ_TO_REDUCE
                           Sequence number to reduce.                   Be
                           default, all sequences found will be reduced.
     -v, --verbose         verbose mode [default]
     -D MASTER_DARK, --master_dark=MASTER_DARK
                           master dark to subtract
     -F MASTER_FLAT, --master_flat=MASTER_FLAT
                           master flat to divide by
     -b BPM_FILE, --bpm_file=BPM_FILE
                           bad pixel mask file
     -g GROUP_BY, --group_by=GROUP_BY
                           kind of data grouping (based on) to do with the
                           dataset files (ot |filter)
     -k, --check_data      if true, check data properties matching
                           (type, expt, filter, ncoadd, mjd)
   
	
.. index:: mask, masking, applymask

Applying User Defined Masks
***************************

Run PAPI with the ``--applymask`` option. PAPI will stop processing after the ``nonlincor`` module and give
instructions on how to run the masking tools on your data.

.. _troubleshooting:

Troubleshooting
***************

As we stated previously, PAPI was developed primarily for reducing imaging data of extragalactic sources. 
Here are some tips for reducing other types of data:

*Add tips here*