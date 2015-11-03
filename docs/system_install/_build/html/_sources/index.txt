.. PANIC System Installation documentation master file, created by
   sphinx-quickstart on Wed Oct 28 08:39:26 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to PANIC Software System Installation's documentation!
==============================================================

.. warning::

   This "Documentation" is still a work in progress; some of the material
   is not organized, and several aspects of PANIC software system installation
   are not yet covered with sufficient detail.

Welcome! PANIC is the `PAnoramic Near Infrared Camera (PANIC) <http://www.iaa.es/PANIC>`_ for the 2.2m 
and 3.5m Telescopes at `Calar Alto Observatory (CAHA) <http://www.caha.es>`_.
This is the documentation for the PANIC software system. We  describe how 
to perform the complete PANIC computer system  installation, from  RAID configuration  
and OS installation on the PANIC computers (PowerEdge R720) to whole PANIC software.

 
The **software** is compound of next main parts:

- GEIRS_: in charge of then instrument control of the wheels, temperature of cryostat and acquistion software for the ROE (Read Out Electronic).

- OT_: the high level software for the desing of observations trough Observation Blocks (OBs) and to make the observations.

- PAPI_: the data reduction software which include pipeline and the Quick-Look (QL)

- LEMON_ (Long-tErm photometric MONitoring): and astronomical pipeline for automated time-series reduction and analysis.



In addition to this html version of the manual, there is also a pdf_ version to download.


**Development:** José-Miguel Ibáñez-Mengual (IAA-CSIC)

**Contribution:** PANIC Team

Contents:

.. toctree::
   :maxdepth: 3
   :numbered:

   system
   os
   geirs
   ot
   papi
   lemon
   
   

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Documentation last updated on |today|


.. _PANIC: http://www.iaa.es/PANIC
.. _CAHA: http://www.caha.es
.. _Omega2000: http://www.caha.es/CAHA/Instruments/O2000/index.html
.. _sphinx: http://sphinx.pocoo.org
.. _pdf: http://www.iaa.es/~jmiguel/PANIC/PAPI/PAPI.pdf
.. _Proc. SPIE 7740 : http://proceedings.spiedigitallibrary.org/proceeding.aspx?articleid=751764
.. _Proc. SPIE 8451: http://proceedings.spiedigitallibrary.org/proceeding.aspx?articleid=1363096
.. _OT: http://www.iaa.es/~agsegura/PANIC_OT/PANIC_Observation_Tool.html
.. _Python: http://www.python.org
.. _GEIRS: http://www2.mpia-hd.mpg.de/~mathar/public/PANIC-SW-DCS-01.pdf
.. _PAPI: http://www.iaa.es/~jmiguel/PANIC/PAPI/html/index.html
.. _LEMON: http://lemon.readthedocs.org/en/latest/