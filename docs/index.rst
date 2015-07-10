.. PAPI documentation master file, created by
   sphinx-quickstart on Wed Dec 14 16:53:57 2011.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

############################
PANIC Pipeline Documentation
############################

:Release: |version|

:Date: |today|

.. warning::

   This "Documentation" is still a work in progress; some of the material
   is not organized, and several aspects of PAPI are not yet covered
   with sufficient detail.

Welcome! This is the Documentation for :ref:`papi`, the Data Reduction Pipeline of PANIC
instrument, and for :ref:`pql`, the GUI of PAPI used during the observation;
all processing routines used by PQL are implented in PAPI.


PAPI is the automatic image processing pipeline for data taken with the 
`PAnoramic Near Infrared Camera (PANIC) <http://www.iaa.es/PANIC>`_ for the 2.2m 
and 3.5m Telescopes at `Calar Alto Observatory (CAHA) <http://www.caha.es>`_. 
The pipeline is written in Python_ and developed at the `Institute of Astrophysics 
of Andalusia (CSIC) <http://www.iaa.es/>`_. The automated processing steps 
include basic calibration (removeing instrumental signature), cosmic-ray removal, 
treatment for electronic ghosts (cross-talk), sky subtraction, non-linear 
count-rate correction, robust alignment and registration. 


This manual is a complete description of the data reduction recipes implemented 
by the PAPI pipeline, showing the status of the current pipeline version
and describing data reduction process of the PANIC data using PAPI.

Although PAPI was developed for the PANIC camera, its development was initiated
using data from Omega2000_ camera at the 3.5m CAHA telescope and HAWK-I_ camera at the 
VLT; thus, in principle it should work with data these instruments, although it is not
optimized for them.

In addition to this html version of the manual, there is also a pdf_ version to download.


**Development:** José M. Ibáñez (IAA-CSIC)

**Contribution:** PANIC Team

Caveat
======

Currently PAPI it is able to reduce data taken with the Observing Tool (OT_) 
defining the required observing blocks (OB), or manually through GEIRS scripts.
PAPI was primarily developed and optimized for reducing broad-band imaging data 
of extragalactic sources (such as imaging data taken for field galaxy surveys and 
galaxy cluster surveys). Other types imaging data have been reduced with PAPI 
but results can not be as good as desired. (See :ref:`troubleshooting` for tips).
PAPI is **not** designed to reduce any kind of field taken with PANIC.  


Contents
========

.. toctree::
   :maxdepth: 4
   :numbered:
   
   install
   quick_look
   running_papi
   reference
   data
   processing
   references
   faq
   troubleshooting
   acknowledgments
   license
   glossary

Citation
========

If your research uses PAPI, we'd appreciate it if you could acknowledge the 
fact by including the following citation:

"This research made use of PAPI, the pipeline of PANIC instrument. It is 
funded by the Spanish Ministry of Economy and Competitiveness with funds 
from the European Union (FEDER) and the Spanish national budget, through the grants 
ICTS-2006-15, ICTS-2007-10, ICTS-2008-24, ICTS-2009-32 and the project 
Intramural 200450E458 of the Spanish National Research Council."
    
- Ibáñez Mengual, J.M., Fernández, M., Rodríguez Gómez, J. R., García Segura, A. J., Storz, C., 
  "The PANIC software system", `Proc. SPIE 7740`_, 77402E (2010)

- Ibáñez Mengual,J.M, García A.J, Storz C., Fried J. W., Fernández M., Rodríguez J. F., 
  "Advanced PANIC quick-look tool using Python", `Proc. SPIE 8451`_, (2012)


   
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
.. _Proc. SPIE 7740 : http://proceedings.spiedigitallibrary.org/proceeding.aspx?articleid=751764
.. _Proc. SPIE 8451: http://proceedings.spiedigitallibrary.org/proceeding.aspx?articleid=1363096
.. _OT: http://www.iaa.es/~agsegura/PANIC_OT/PANIC_Observation_Tool.html
.. _Python: http://www.python.org
