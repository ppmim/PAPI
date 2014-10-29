.. PAPI documentation master file, created by
   sphinx-quickstart on Wed Dec 14 16:53:57 2011.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

PANIC Data Pipeline Documentation
=================================
.. warning::

   This "Documentation" is still a work in progress; some of the material
   is not organized, and several aspects of PAPI are not yet covered
   with sufficient detail.

Welcome! This is the Documentation for PAPI, the Data Reduction Pipeline of PANIC
instrument (PANIC Pipeline, version |version|, |today|). 


PAPI is the automatic image processing pipeline for data taken with the 
`PAnoramic Near Infrared Camera (PANIC) <http://www.iaa.es/PANIC>`_ for the 2.2m 
and 3.5m Telescopes at `Calar Alto Observatory (CAHA) <http://www.caha.es>`_. 
The pipeline is written in Python and developed at the `Institute of Astrophysics 
of Andalusia (CSIC) <http://www.iaa.es/>`_. The automated processing steps 
include basic calibration (removeing instrumental signature), cosmic-ray removal, 
treatment for electronic ghosts (cross-talk), sky subtraction, non-linear 
count-rate correction, robust alignment and registration. 


This manual is a complete description of the data reduction recipes implemented 
by the PAPI pipeline, showing the status of the current pipeline version
and describing data reduction process of the PANIC data using PAPI.

Although PAPI was developed for the PANIC camera, it also works with data from 
Omega2000_ NIR camera at the 3.5m CAHA telescope and HAWK-I_ camera at the VLT.

In addition to this html version of the manual, there is also a pdf_ version to download.


**Development Team:** José M. Ibáñez (IAA-CSIC)

Caveat
******

Currently PAPI it is able to reduce data taken with the Observing Tool (OT) 
defining the required observing blocks (OB), or manually through GEIRS scripts.
PAPI was primarily developed and optimized for reducing broad-band imaging data 
of extragalactic sources (such as imaging data taken for field galaxy surveys and 
galaxy cluster surveys). Other types imaging data have been reduced with PAPI 
but results can not be as good as desired. (See :ref:`troubleshooting` for tips).
PAPI is **not** designed to reduce any kind of field taken with PANIC.  

   
Contents
********
.. toctree::
   :maxdepth: 3
   :numbered:
   
   install
   running
   quick_look
   ref
   data
   photo
   reference
   processing
   acknowledgments
   references
   faq
   troubleshooting
   citation
   license
   glossary


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Documentation last updated on |today|


.. _PANIC: http://www.iaa.es/PANIC
.. _CAHA: http://www.caha.es
.. _Omega2000: http://www.caha.es/CAHA/Instruments/O2000/index.html
.. _HAWK-I: http://www.eso.org/sci/facilities/paranal/instruments/hawki/
.. _sphinx: http://sphinx.pocoo.org
.. _pdf: http://www.iaa.es/~jmiguel/PANIC/PAPI/PAPI.pdf
