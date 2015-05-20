PANIC Quick-Look Tool (PQL)
===========================

Purpose
*******

The PANIC Quick-Look (hereafter PQL) will  perform some on-line data processing 
for quick-look or quality check of the data being acquired, taking a close look 
at a raw near-infrared image and getting a quick feedback of the running observation.

The PQL is an application with a graphical user interface which monitors the 
GEIRS data output, waiting for new FITS files coming from GEIRS. When a new file 
is detected, it is added to the file list view in the main panel, and then the 
PQL will perform the task previously specified by the user in the setup 
configuration. Some of the available tasks are:

   * Only display the FITS image with no processing
   * Dark subtraction, flat division
   * Sky subtraction (using N-nearest frames or own sky )
   * Field distortion removal
   * Image align and stacking
   * Preliminary astrometric solution
   * Preliminary photometry

In addition, the PQL allows the user to execute manually in an interactive way 
some tasks with the data. For example, he user will be able to select a file, 
compute some statistics values (background, FWHM, min, max, …) or ask for the
sky subtraction looking for the nearest N frames around the selected one. Other 
option available is to select a set of files and request a shift and align of 
them.

The PQL can be operated in both near-real time mode (during the observation) and
offline mode (after the observation, with all data files already stored in the disk);
however, its functionalities have been provided mainly in near-real time to check 
the status and progress of the observation during the night. 



.. index:: quick-look, running


FITS files and headers
**********************
QuickLook only supports FITS_ (Flexible Image Transport System) image formats. Due
PANIC hava an FPA of four detector, the FITS files can be ``Single Extension FITS (SEF)`` 
or ``Multi-Extension FITS (MEF)``, however MEF are prefered.

For general purpose, such as viewing and simple analysis, only minimal headers
keywords are required. However, and in order to group and reduce observing sequences, 
the following header keywords are also required::

    OBS_TOOL= 'OT_V1.1 '           / PANIC Observing Tool Software version          
    PROG_ID = '        '           / PANIC Observing Program ID                     
    OB_ID   = '6       '           / PANIC Observing Block ID                       
    OB_NAME = 'OB CU Cnc Ks 2'     / PANIC Observing Block Name                     
    OB_PAT  = '5-point '           / PANIC Observing Block Pattern Type             
    PAT_NAME= 'OS Ks 2 '           / PANIC Observing Secuence Pattern Name          
    PAT_EXPN=                    1 / PANIC Pattern exposition number                
    PAT_NEXP=                    5 / PANIC Pattern total number of expositions      
    IMAGETYP= 'SCIENCE '           / PANIC Image type                         


These keywords are automatically added to the FITS header by the PANIC Observation Tool,
as each file is created. If these are not saved, the PQL will not work correctly.

Starting the PQL
****************

To start PQL GUI, you can lauch it from the PANIC computer (panic22/panic35) once you are
logged as obs22/obs35 user. Thus, as any one of the workstations of the observing room,
open a X terminal window and log into the PANIC computer as follow:
  
for 2.2m::

    $ ssh -X obs22@panic22 
    (ask Calar Alto staff for password)
   
for 3.5m::

    $ ssh -X obs35@panic35 
    (ask Calar Alto staff for password)
   
Once you are logged into the PANIC computer, to lauch PQL GUI type next command::


    $ start_ql &
    
The next figure shows a snapshot of the main window of the PQL GUI that will bring up the *start_ql* command.
  
.. image:: _static/PQL_GUI_main_window.png
   :align: center
   :scale: 65 %



Configuration files
*******************

The configuration files used by the PQL are located in the $PAPI_HOME/config_files.
The main config file is the same file used by PAPI, ie., $PAPI_CONFIG, and usually
called papi.cfg.

At the end of the $PAPI_CONFIG file, there is section called 'quicklook', where the
user can set next parameters::

    # Next are some configurable options for the PANIC Quick Look tool
    #
    # some important directories
    #
    source = /data1/PANIC/
    output_dir = /data2/out   # the directory to which the resulting images will be saved.
    temp_dir = /data2/tmp    # the directory to which temporal results will be saved
    verbose = True

    # Run parameters
    run_mode = Lazy # default (initial) run mode of the QL; it can be (None, Lazy, Prereduce)


Although the user can edit these values in the config file, they can be set easily
on the PQL GUI. 

PQL's Main Window
*****************

The PQL main window contains a menu bar (1), tool bar (2), four tabbed panels (3) and 
an event log window (4).
Images are displayed in an external well-known application, ds9_. Plots results are displayed in 
the additional windows, usually generated by matplotlib than can be popied to the clipboard, 
printed and saved.

Menu bar
********

The menu bar provides acces to some PQL's capabilities.

1. File
2. View
3. Settings
4. Calibrations
5. Tools
6. Help
  Opens a web browser which shows an on-line HTML version of this user's manual. This will fail 
  if the internet conection and proxy is not correctly configured.
7. Exit
  Quit the PQL application.


Buttons bar
***********

The button bar duplicates some of the options available from the menu bar or the pop-up menu. 
The buttons provide quick access to change the most frecuently-used PQL actions:

- add a file to the current view
- change the source input directory
- display the current selected image 
- open an IRAF_ console
- open Aladin_ tool

.. image:: _static/PQL_GUI_toolbar.png
   :align: center
   :scale: 65 %
   

Main Panel
**********
This tab panel contains the following controls:

- Input directory
- Ouput directory
- Filename filter
- Data list view
- List view filter
- QL mode
- 'Subract last-2' button
- 'START processing' button
- 'Create Calibrations' button


Data Directories
----------------

In the 'Main' tab panel of the PQL main window, the fitst thing to set up are the data directories:

.. image:: _static/PQL_GUI_data_dirs.png
   :align: center
   :scale: 65 %



Input directory
^^^^^^^^^^^^^^^

This is where you tell PQL where the data are or being saved by GEIRS. This directory is specified
at the beggining of the night on the Observation Tool. PQL requieres all data to lie in some main 
directory, not being required to distribute the files in individual sub-directories for darks, flats,
and science images. It is advised that this directory follow the next format::

    /data1/PANIC/YYYYMMDD

To set the value, the user must push the 'Input Dir' button:

.. image:: _static/PQL_GUI_input_dir_but.png
   :align: center
   :scale: 65 %
    
Output directory
^^^^^^^^^^^^^^^^

This is where you tell PQL where the data generated by the PQL, as result of some processing, will be saved.
This directory must also be specified at the begining of the night, and is advised to follow the next format::

   /data2/out_YYYYMMDD
  

To set the value, the user must push the 'Output Dir' button:

.. image:: _static/PQL_GUI_output_dir_but.png
   :align: center
   :scale: 65 %


Temporal directory
^^^^^^^^^^^^^^^^^^

This is where you tell PQL where the temporal files generated by the PQL, as result of some processing, 
will be saved, and probably deleted after at the end of that processing.
This directory must also be specified at the begining of the night, and is advised to follow the next format::

   /data2/tmp_YYYYMMDD

To set the value, the user must push the 'Temporary Dir' button than appears on the 'Setup' tab, 
instead the 'Main' tab used for input and output directory.


.. image:: _static/PQL_GUI_tmp_dir.png
   :align: center
   :scale: 65 %
   

Filename Filter 
---------------

In this box, the user can filter the name of the files should appears on the data list view 
from the input directory (output files are not filtered).
The filter can contains '*' and '?' wildcards. 

For example:

    `*March10_00?1*`

.. image:: _static/PQL_GUI_filter.png
   :align: center
   :scale: 65 %

.. _data_list_view:

Data list view
--------------
Tha data list view control displays all the files found in the input directory, or in the output directory 
if the check box at the right of output directory is checked. Additionaly, the use can add any other FITS file.
The control is a multicolum table with the next fields:

.. image:: _static/PQL_GUI_data_list_view.png
   :align: center
   :scale: 65 %

Filename
  Full path name of the file found in the 
Image type
  The type of the FITS file detected: DARK, DOME_FLAT, SKY_FLAT, FOCUS, SCIENCE 
ExpT
  Exposition time of the file (EXPTIME keyword)
Date-Obs
  Observation data of the file (DATE-OBS keyword)
Object
  Object name (OBJECT keyword)
RA
  Right ascention of center of the image.
Dec
  Declination of the cener of the image.
 
List view filter
----------------
It allows to select the type of files to be shown in the data list view. The options are:


INPUTS
  Files of the input directory
OUTS
  Files of the ouput directory
DARK
  Files marked (IMAGETYP) as DARK images
DOME_FLAT
  Files marked as DOME_FLAT image  
FOUCS
  Files marked as FOCUS image from a focus series
SKY_FLAT
  Files marked as SKY_FLAT images
SCIENCE
  Files marked as SCIENCE image or with unknown type.
MASTERS
  Files marked as MASTER calibration files produced by PAPI
REDUCED
  Files marked as calibrated by PAPI
GROUP
  Special case that show all the files groupped as sequences
ALL
  Show all the files, not matter the type of it
  
 
.. image:: _static/PQL_GUI_listview_filter.png
   :align: center
   :scale: 65 %

   
QuickLook Mode
--------------

The quick look mode combo box allows to select the mode in which the PQL will be run when the **START processing** button is pushed.
The current modes are:

None
  No processing action is done

Lazy
  If the end of a calibration (DARK, FLAT) sequence is detected, the master file is built. Otherwise,
  and the SCIENCE files are processed as specified in the 'Setup->Lazy Mode':
  
  + Apply DARK + FLAT + BPM
  + Subtract Last Frame (Science)
  + Subract Nearest Sky

.. image:: _static/PQL_GUI_qlmodes.png
   :align: center
   :scale: 65 %

  
Pre-Reduction
  If the end of observing sequence is detected, it is processed in a quick mode (single pass for sky subtraction). 
  For calibration sequences, the master file will be built, and for science sequences, a quick 
  reduction will be done using options configured in the 'Setup->Pre-Reduction Mode' and the 
  calibrations found in local database (output directory and external calibration directory).
  Note that the pre-reduction options configured in the config file will be overwritten.
  
.. image:: _static/PQL_GUI_pre-redmode.png
   :align: center
   :scale: 65 %
  
Quick-LEMON
  The same as Pre-reduction, but the processing stops after the 1st sky subtraction, and 
  no final co-added image is produced. It is useful for LEMON_ processing for light curves.

Full-Reduction
  If the end of observing sequence is detected, it is processed in a *science* mode (double pass for sky subtraction). 
  For calibration sequences, the master file will be built, and for science sequences, a *science* 
  reduction will be done using options configured in the 'Setup->Pre-Reduction Mode' and the 
  calibrations found in local database (output directory and external calibration directory).
  Note that the pre-reduction options configured in the config file will be overwritten.

Full-LEMON
  The same as Pre-reduction, but the processing stops after the 2nd sky subtraction, and 
  no final co-added image is produced. It is useful for LEMON_ processing for light curves.



Last file received
------------------
This field shows the last file received (detected) by the PQL.


Buttons
-------

Subract-last2 button
^^^^^^^^^^^^^^^^^^^^
It will produced a new image as result of the subtraction of last two images received.

Create calibrations button
^^^^^^^^^^^^^^^^^^^^^^^^^^

This button will start the processing of all the calibration
sequences received. As result, a list of master calibrations will be generated
in the output directory. 

START button
^^^^^^^^^^^^

This button starts the processing of all the sequences received. You will be 
asked whether to proccess all the current images or only the new ones. 
As result, a list of master calibrations and science calibrated images will be generated
in the output directory. 

Add button
^^^^^^^^^^
This button allows to add manually a single file to the *Data List View* from wherever the
file is.


Remove button
^^^^^^^^^^^^^
This button removes manually from the *Data List View* the currently selected file, but it
does not remove neither from the local database nor the file system.


Clear All button
^^^^^^^^^^^^^^^^
It removes all the current files from the *Data List View*, but they will not be removed
from the file system. As result, it will empty the *Data List View* 
until a new input directory is selected or a new file is detected in the current one.


Setup Panel
***********

Calibrations Panel
******************

Log Panel
*********

Pop-up Menu
***********

It is a context pop-up menu that appears when the user select a file (or
a set of them) in the *Data List View* and click the right mouse button.
Next figure shows the options of that pop-up menu:


.. image:: _static/PQL_GUI_pop_up.png
   :align: center
   :scale: 65 %

Some actions in the menu could be disabled and greyed out if they are not
availabe or applicable to the selected files.
   
Display image
-------------
It displays the currect selected image in the ds9_ display. it will launch the ds9 
application if it is not opened.

Image info
----------
It is a quick way to see some basic information of the selected image. The information
is mainly concerning the FITS structure and exposition times used. The information will
be shown in the *Event Log Console* as follow:

::

  ---------------
  SEF Filename : /data1/PANIC/2015-05-19_SS_zenith_Ks_1_3/SS_Ks_SG1_4_0024.fits
  Image Shape : (32, 32)
  Filter : Ks 
  ITIME : 0.045000 
  NCOADDS : 1 
  EXPTIME : 0.045000 
  TYPE : FOCUS 
  OT keywords : True 
  ---------------

Of course, if you need any other information of the file, you can find it using
the 'ds9->File->Display Header...' option.


Copy files to clipboard
-----------------------
It copies the current selected file to the clipboard. This way you could paste the 
full pathnames to any other file. It is quite useful when using the PAPI command
line to run some operation that is not available on the PQL.
  
Copy files to text file
-----------------------
If copies the current selected files into the specified text file. It is quite useful 
when using the PAPI command line to run some operation that is not available on the PQL.

Show Dither pattern
-------------------
It brings up a plot with the dither offsets obtained from the RA,Dec coordinates 
of the FITS header. You have to select a set of images in the *Data List View* and
then right-button and *Show Dither pattern*.

.. image:: _static/PQL_GUI_dither_pat_ex.png
   :align: center
   :scale: 25 %
   

.. _calibrations:

Calibrations
************
Next options allow to build the master calibration files from a given set of selected files.


Build Master Dark
-----------------
This command is used to produce a master dark file from a set of files currectly selected 
in the *Data List View*. It checks that all the selected files are compliant, ie., 
have the same EXPTIME, NCOADD, ITIME, READMODE and shape. You only have to give the name of 
the master dark file to be created.

The master dark is computed using an average combine with a minmax rejection algorithm.
   

Build Master Dome Flat
----------------------
This command is used to produce a Master DOME FLAT file from a set of files currectly selected 
in the :ref:`Data List View <data_list_view>`. It checks that all the selected files are compliant, ie., 
have the same FILTER, NCOADD, READMODE and shape. You have to select at least one DOME_FLAT_LAMP_OFF 
and one DOME_FLAT_LAMP_ON image, and then provide the name for the master dome flat to create.

The procedure to create the master dome flat is as follow: 

    #. Check the EXPTIME , TYPE(dome) and FILTER of each Flat frame
    #. Separate lamp ON/OFF dome flats
    #. Make the median combine + sigmaclip of Flat LAMP-OFF frames 
    #. Make the median combine + sigmaclip of Flat LAMP-ON frames
    #. Subtract lampON-lampOFF (implicit dark subtraction)
    #. (optionally) Normalize the flat-field with median (robust estimator)
            
    Note that we do **not** need to subtract any MASTER_DARK; it is not required for 
    DOME FLATS (it is done implicitly because both ON/OFF flats are taken 
    with the same Exposition Time).

Build Master Twlight (sky) Flat
-------------------------------
This command is used to produce a Master SKY FLAT file from a set of files currectly selected 
in the :ref:`Data List View <data_list_view>`. It checks that all the selected files are compliant, ie., 
have the same FILTER, NCOADD, READMODE and shape. You have to select at least three SKY_FLAT
images (dusk or dawn). The procedure will look for the required master dark frames to subtract 
in the current output directory and in the external calibration directory. If some of the master dark
are not found, then the procedure will fail.

The procedure to create the master sky flat is as follow:

    #. Check the  TYPE (sky flat) and FILTER of each Flat frame
       If any frame on list missmatch the FILTER, then the master 
       twflat will skip this frame and continue with then next ones.
       EXPTIME do not need be the same, so EXPTIME scaling with 'mode' 
       will be done. 
    
    #. Check either over or under exposed frames ( [10000 < mean_level < 40000] ADUs )
        
    #. We subtract a proper MASTER_DARK, it is required for TWILIGHT FLATS 
       because they might have diff EXPTIMEs.
        
    #. Make the combine (median + sigclip rejection) the dark subtracted Flat 
       frames scaling by 'mode'.
        
    #. Normalize the sky-flat wrt SG1 detector, dividing by its mean value.
    

Build GainMap
-------------
This command is used to produce a Master GainMap file from a set of files currectly selected 
in the :ref:`Data List View <data_list_view>`. It checks that all the selected files are compliant, ie., 
have the same FILTER, NCOADD, READMODE and shape. You have to select at least three
flat frames (dome, dusk or dawn). For sky flats, the procedure will look for the required master dark 
frames to subtract in the current output directory and in the external calibration directory. 
If some of the master dark are not found, then the procedure will fail. Dome flat do not need
dark subtraction.

The procedure to create the master sky flat is as follow:

    #. Check the  TYPE (sky flat) and FILTER of each Flat frame
       If any frame on list missmatch the FILTER, then the master 
       twflat will skip this frame and continue with then next ones.
       EXPTIME do not need be the same, so EXPTIME scaling with 'mode' 
       will be done. 
       
    #. Create the proper master dome/sky flat.
    
    #. Once the master dome flat is created, the procedure will 
    compute the gainmap as follow:
    


Photometric calibration
-----------------------
.. note::

   Your **data is assumed to be calibrated**. Dark subtraction, flat-fielding correction and any 
   other necessary steps should have been performed before any data is fed to the photometric 
   calibration.


    
.. _howto:
   
How to ...?
***********

How to determine the telescope focus ?
--------------------------------------


How to inspect the profile of the stars in an image ?
-----------------------------------------------------
You should follow the next steps:

1. select in the *Data List View* the image to inspect.
2. double-click to display the image into ds9 and zoom to the area you wish to inspect
3. go to the tool bar (or Tool menu) and open an IRAF console
4. type in the iraf console 'imexam'
5. focus the mouse cursor on the ds9_ display and type the *imexam* comand you wish
   for the inspection. For example, type ***r*** to show the *radial profile* of 
   the selected star
6. once you have finished the inspection, type q to exit from *imexam*



How do I make mosaics with PQL? 
-------------------------------
PAPI will automatically warp (using SWARP) your images as thre are located on the sky. 

How do I make use of parallelisation ?
--------------------------------------
Just be sure the number of *parallel* parameter is set to *True* on the $PAPI_CONFIG file.






.. index:: quicklook, off-line, on-line, configuration

.. _FITS: http://fits.gsfc.nasa.gov
.. _IRAF: http://iraf.noao.edu/
.. _ds9: http://ds9.si.edu/site/Home.html
.. _Aladin: http://aladin.u-strasbg.fr
.. _LEMON: https://lemon.readthedocs.org/