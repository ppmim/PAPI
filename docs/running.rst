.. _papi:

Running PAPI
============

.. index:: quickstart, running


Quickstart
**********

Running PAPI can be as simple as executing the following command in a terminal::
	
	$ papi.py -s raw_data -d result 

Where ``raw_data`` is the directory of the raw dataset (uncalibrated) having 
both science or calibration files, and ``result`` is the path to the directory 
where the calibrated data produced by the pìpeline will be saved.

Example::

   $ papi.py -s /my/raw_data/directory -d /my/result/directory

.. index:: uncalibrated, data



Optional Commands
*****************

For most image sets PAPI can be run in the default configuration with no 
additional interaction required. If the default settings are insufficient for 
processing a particular data set, there are a number of run-time options which 
may be applied to help improve the reductions.

The next command will show some of the available options::

   $ papi.py --help


Then, the listing of the PAPI command line options::

    Usage: papi.py [OPTION]... DIRECTORY...
    
    This is the main module of the PANIC data reduction system (PAPI)
    
    Options:
      --version             show program's version number and exit
      -h, --help            show this help message and exit
      -c CONFIG_FILE, --config=CONFIG_FILE
                            Config file for the PANIC Pipeline application.
                            If not specified, './config_files/papi_suse11.cfg' is
                            used.
      -C, --Check           Check if versions of PAPI modules are right.
      -s SOURCE, --source=SOURCE
                            Source file list of data frames. It can be a file
                            or directory name.
      -o OUTPUT_FILE, --output_file=OUTPUT_FILE
                            Final reduced output image
      -t TEMP_DIR, --temp_dir=TEMP_DIR
                            Directory for temporal files
      -d OUTPUT_DIR, --out_dir=OUTPUT_DIR
                            Output dir for product files
      -r ROWS, --rows=ROWS  Use *only* files of the source file-list in the range
                            of rows specified (0 to N, both included)
      -R, --recursive       Does recursive search for files in source directory
      -l, --list            Generate a list with all the source files read from
                            the source only sorted by MJD
      -M REDUCTION_MODE, --red_mode=REDUCTION_MODE
                            Mode of data reduction to do (lemon|quick|science)
      -m OBS_MODE, --obs_mode=OBS_MODE
                            Observing mode (dither|ext_dither|other)
      -p, --print           Print detected sequences in the Data Set
      -S SEQ_TO_REDUCE, --seq_to_reduce=SEQ_TO_REDUCE
                            Sequence number to reduce. By default,
                            all sequences found will be reduced.
      -D MASTER_DARK, --master_dark=MASTER_DARK
                            master dark to subtract
      -F MASTER_FLAT, --master_flat=MASTER_FLAT
                            master flat to divide by
      -b BPM_FILE, --bpm_file=BPM_FILE
                            bad pixel mask file
      -g GROUP_BY, --group_by=GROUP_BY
                            kind of data grouping (based on) to do with the
                            dataset files (ot |filter)
      -k, --check_data      if true, check data properties matching (type, expt,
                            filter, ncoadd, mjd)
      -v, --verbose         Verbose mode [default]

  
	
Configuration files
*******************
PAPI has a set of configuration files required to run properly. They are the next
ones:

   * papi.cfg:  main configuration file

      In addition to the command line options, PAPI has a configuration file in 
      which the user can set both the command line options  and a wider set of 
      additional ones. 
      This config file can be specified with the ``-c`` option, but by default it 
      is looked for it in the ``config_files`` directory defined by PAPI_CONFIG 
      environment variable.

   * scamp.cfg: SCAMP configuration file
   * swarp.conf: SWARP configuration file
   * sextractor.sex : SExtractor configuration file
   * sextractor.conf: 
   * sextractor.cong:
   * sextractor.nnw:
   * sextractor.param:
   
    
.. index:: run, command line, config


Examples
********
Main config file
----------------

File papi.cfg::


    # Default configuration file for PAPI 1.0
    # Updated 27 Feb 2012  
    
    ############################################################################
    [general]
    ############################################################################
    
    #
    # some important directories
    #
    #source = /home/jmiguel/DATA/SIMU_PANIC_3/q1.txt   # it can be a directory or a text file with a list of filenames to be processed
    source = /mnt/GEIRS_DATA
    #source = /data/SIMU # default directory for source raw data files
    output_dir = /data/out   # the directory to which the resulting images will be saved.
    temp_dir = /data/tmp     # the directory to which temporal results will be saved
    output_file = /tmp/reduced.fits  # final reduced produced image
    parallel = True
    ncpus = 2
    verbose = True
    logfile = /tmp/papi.log # to be implemented !!!
    
    #
    #reduction_mode : reduction modo to do with the raw science files
    #
    reduction_mode = quick   # default reduction mode (quick|science|lemon)
    
    #
    obs_mode = dither  #default observing mode of input data files to reduce (dither|ext_dither|other)
    #
    
    # if any, default master calibration files to use
    #master_dark = None
    #master_flat = None
    #master_bpm = None
    
    #
    # External calibration DataBase: directory used as an external calibration database.
    # Then, if during the reduction of a ReductionSet(RS) no calibration (dark, flat) 
    # are found in the current RS, then PAPI will look for them into this directory.
    # If the directory does not exists, or no calibration are found, then no calibrations
    # will be used for the data reduction.
    # Note that the calibrations into the current RS have always higher priority than
    # the ones in the external calibration DB.
    #
    ext_calibration_db = /data2/out/cali_O2K
    
    
    #
    # check data integrity. It consists in checking if TEXP,NCOADD,FILTER and READMODE match properly
    #
    check_data = True
    
    #
    # Remove crosstalk. If True, a procedure to remove the crosstalk will be executed
    # just after the 2nd. sky subtraction. (both O2K or PANIC)
    #
    remove_crosstalk = False
    
    #
    # Cosmic-Ray Removal. If True, a procedure to remove the CR will be executed
    # just after the 2nd. sky subtraction. (both O2K or PANIC)
    #
    remove_cosmic_ray = False
    
    #
    # Purge output. If True, a procedure to remove the temporal or intermediate files
    # (.list, .objs., .ldac, .xml, ...) will be removed from the output directory
    # just after the end of the RS reduction.
    #
    purge_output = False
    
    
    # min_frames : minimun number of frames required to reduce a sequence
    #
    min_frames = 5
    
    #
    # group_by: the pipeline will try to group the data files in two main ways: 
    #           (OT) following the specific keywords provided by the OT as OB_ID, OB_PAT, IMAGETYP, FILTER
    #           and then different observing sequences could be grouped and reduced or
    #           (FILTER) only group by filter band, and then only one observing sequence should be provided
    #           (NONE) No grouping criteria will be taken; force only one group with all the files 
    #
    group_by = filter # (OT or FILTER or NONE)
    
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # The ABOVE option values can be modified at the invokation time of the pipeline in the command line
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    #
    # apply_dark_flat : 0  Neither dark nor flat field will be applied
    #                   1  the pipeline will look for a master dark and master flat 
    #                      field to be applied to the raw science frames.
    #                      Both master DARK and FLAT are optional,i.e., each one can be applied 
    #                      even the other is not present.
    #                   2  master flat will be looked for to be applied AFTER 
    #                      skysubtraction, but no DARK will be subtracted (it is 
    #                      supposed to be done by the skysubtraction) 
    #                   (some people think they are not required !)
    apply_dark_flat = 1 
    
    #
    # some other values (really required ?)
    #
    max_mjd_diff = 6.95e-3 # Maximun seconds (10min=600secs aprox) of temporal distant allowed between two consecutive frames (1/86400.0)*10*60
    
    
    scale = 0.0001333   # scale of the image, in degrees per pixel
                        # (0.0001333 deg/pix = 0.48 arcsec/pixel)
    equinox = 2000      # equinox in years
    radecsys = ICRS     # reference system
    pattern = *.fits    # if specified, only those images that match the pattern (according to the rules used by the Unix shell) will be
                        # considered when autodetecting FITS images in _directories_ no tilde expansion is done, but *, ?, and character
                        # ranges expressed with [] will be correctly matched. NOTE: it is because this feature that images like flatV...
                        # or discarded_.... specify its type at the beginning of they filename (vamos, porque no hay forma de negar un 'match')
    
    filter_name_Z = Z   # the key stored in the FITS header when the filter is Z
    filter_name_Y = Y
    filter_name_J = J
    filter_name_H = H, Filter_H    # admits list of strings if multiple values are possible
    filter_name_K = K
    filter_name_Ks = KS
    
    
    ############################################################################    
    [config_files]
    ############################################################################
    
    terapix_bin = /usr/local/Terapix/bin
    irdr_bin = /home/panic/DEVELOP/PIPELINE/PANIC/trunk/irdr/bin
    #irdr_bin = /home/jmiguel/DEVELOP/PIPELINE/PANIC/trunk/irdr/bin
    sextractor_conf = /home/panic/DEVELOP/PIPELINE/PANIC/trunk/config_files/sextractor.sex     # SExtractor configuration file
    
    sextractor_conf = /home/panic/DEVELOP/PIPELINE/PANIC/trunk/config_files/sextractor.sex     # SExtractor configuration file
    sextractor_param = /home/panic/DEVELOP/PIPELINE/PANIC/trunk/config_files/sextractor.param  # File containing the list of parameters that will be computed and put in the catalog for each object
    sextractor_nnw = /home/panic/DEVELOP/PIPELINE/PANIC/trunk/config_files/sextractor.nnw      # File containing the neutal-network weights for star/galaxy separation
    sextractor_conv = /home/panic/DEVELOP/PIPELINE/PANIC/trunk/config_files/sextractor.conv    # File containing the filter definition
    scamp_conf = /home/panic/DEVELOP/PIPELINE/PANIC/trunk/config_files/scamp.conf              # SCAMP configuration file
    swarp_conf = /home/panic/DEVELOP/PIPELINE/PANIC/trunk/config_files/swarp.conf              # SWarp configuration file
    
    
    ############################################################################        
    [dark]  
    ############################################################################
    
    # object_names: in order to make it possible to work in batch mode, is it
    # possible to run the PANIC dark module in all the images, specifying in
    # this parameter which ones will be considered. That is, only those images 
    # whose object name matchs one of the names listed in this parameter will be
    # considered when generating the master dark.
    #
    # Note that if '*' is contained in the list, _all_ object names will be matched.
    # This symbol, thus, provides a way to easily specify all the images, which is
    # equivalent to saying "do not filter images by their object names".
    # 
    object_names = dark
    
    # check_prop : if true, the dark frames used to build the master dark will be 
    # checkd to have the same acquisition properties (EXPT,NCOADD,ITIME, READMODE)
    #
    check_prop = yes
    
    
    # suffix: the string, if any, to be added to the filename of each resulting
    # image. For example, for suffix = "D" and the imput file /home/images/ferM_0720_o.fits,
    # the resulting image would be saved to /home/images/ferM_0720_o_D.fits.
    # This parameter is optional, as if nothing is specified, nothing will be appended
    #
    suffix = D
    
    
    # min_frames : minimun number of frames required to build a master dark
    #
    min_frames = 5
    
    
    
    ############################################################################
    [dflats] 
    ############################################################################
    
    # object_names: in order to make it possible to work in batch mode, is it
    # possible to run the PANIC flat module in all the images, specifying in
    # this parameter which ones will be considered. That is, only those images 
    # whose object name matchs one of the names listed in this parameter will be
    # considered when generating the master dome flat
    #
    # Note that if '*' is contained in the list, _all_ object names will be matched.
    # This symbol, thus, provides a way to easily specify all the images, which is
    # equivalent to saying "do not filter images by their object names".
    # 
    object_names = DOME_FLAT_LAMP_OFF, DOME_FLAT_LAMP_ON
    
    # check_prop : if true, the frames used to build the master  will be 
    # checkd to have the same acquisition properties (EXPT,NCOADD,ITIME, READMODE, FILTER)
    #
    check_prop = yes
    
    # suffix: the string, if any, to be added to the filename of each resulting
    # image. For example, for suffix = "D" and the imput file /home/images/ferM_0720_o.fits,
    # the resulting image would be saved to /home/images/ferM_0720_o_D.fits.
    # This parameter is optional, as if nothing is specified, nothing will be appended
    #
    suffix = F
    
    
    # min_frames : minimun number of frames required to build a master dome flat
    #
    min_frames = 5
    
    area_width = 1000       # length in pixels of the central area used for normalization
    
    # median_smooth: median filter smooth of combined FF to reduce noise and improve
    # the S/N and preserve the small-scale (high-frequency) features of the flat
    # 
    median_smooth = False
    
    ############################################################################
    [twflats]
    ############################################################################
    
    
    # object_names: in order to make it possible to work in batch mode, is it
    # possible to run the PANIC flat module in all the images, specifying in
    # this parameter which ones will be considered. That is, only those images 
    # whose object name matchs one of the names listed in this parameter will be
    # considered when generating the master twflat
    #
    # Note that if '*' is contained in the list, _all_ object names will be matched.
    # This symbol, thus, provides a way to easily specify all the images, which is
    # equivalent to saying "do not filter images by their object names".
    # 
    object_names = TW_FLAT_DUSK, TW_FLAT_DUSK, SKY_FLAT
    
    # check_prop : if true, the  frames used to build the master will be 
    # checkd to have the same acquisition properties (EXPT,NCOADD,ITIME, READMODE, FILTER)
    #
    check_prop = yes
    
    # suffix: the string, if any, to be added to the filename of each resulting
    # image. For example, for suffix = "D" and the imput file /home/images/ferM_0720_o.fits,
    # the resulting image would be saved to /home/images/ferM_0720_o_D.fits.
    # This parameter is optional, as if nothing is specified, nothing will be appended
    #
    suffix = F
    
    
    # min_frames : minimun number of frames required to build a master twlight flat
    #
    min_frames = 5
    
    area_width = 1000       # length in pixels of the central area used for normalization
    
    # median_smooth: median filter smooth of combined FF to reduce noise and improve
    # the S/N and preserve the large-scale features of the flat
    # 
    median_smooth = False
    
    ############################################################################
    [gainmap] 
    ############################################################################
    
    # object_names: in order to make it possible to work in batch mode, is it
    # possible to run the PANIC gainmap module in all the master flat images, specifying in
    # this parameter which ones will be considered. That is, only those images 
    # whose object name matchs one of the names listed in this parameter will be
    # considered when generating the gain map.
    #
    # Note that if '*' is contained in the list, _all_ object names will be matched.
    # This symbol, thus, provides a way to easily specify all the images, which is
    # equivalent to saying "do not filter images by their object names.
    # 
    object_names = MASTER_SKY_FLAT, MASTER_DOME_FLAT, MASTER_TW_FLAT
    
    mingain = 0.6 # pixels with sensitivity < MINGAIN are assumed bad (0.7) 
    maxgain = 2.0 # pixels with sensitivity > MAXGAIN are assumed bad (1.3)
    nsigma = 5    # badpix if sensitivity > NSIG sigma from local bkg (5.0)
    nxblock = 16  # image size should be multiple of block size (16)
    nyblock = 16  # (16)
    normalize = yes # if 'yes' apply a previous normalization to master flat images 
      
    area_width = 1000   # area to use for normalization (1000) 
    
    ############################################################################
    [skysub] 
    ############################################################################
    
    # object_names: in order to make it possible to work in batch mode, is it
    # possible to run the PANIC skysubtration module in all the images, specifying in
    # this parameter which ones will be considered. That is, only those images 
    # whose object name matchs one of the names listed in this parameter will be
    # considered when generating the master dark.
    #
    # Note that if '*' is contained in the list, _all_ object names will be matched.
    # This symbol, thus, provides a way to easily specify all the images, which is
    # equivalent to saying "do not filter images by their object names".
    #
    object_names = SKY, SKY_FOR
    
    # check_prop : if true, the dark frames used to build the master  will be 
    # checkd to have the same acquisition properties (EXPT,NCOADD,ITIME, READMODE, FILTER)
    #
    check_prop = yes
    
    # suffix: the string, if any, to be added to the filename of each resulting
    # image. For example, for suffix = "D" and the imput file /home/images/ferM_0720_o.fits,
    # the resulting image would be saved to /home/images/ferM_0720_o_D.fits.
    # This parameter is optional, as if nothing is specified, nothing will be appended
    #
    suffix = S
    
    #
    # min_frames : minimun number of frames required to build a master super flat
    #
    min_frames = 5
    
    # half width of sky filter window in frames
    #
    hwidth = 2 
    
    area_width = 1000       # length in pixels of the central area used for normalization
    
    # Object mask
    mask_minarea = 5   # sex:DETECT_MINAREA used for object masking
    mask_thresh = 1.5   # sex:DETECT_THRESH used for object masking
    #expand_mask = 0.5   # amount to expand the object mask regions
    satur_level = 3000000 # sex:SATUR_LEVEL; level (in ADUs) at which arises saturation
    
    # skymodel : sky model used used during the sky subtraction. It will be a 
    #             parameter for the IRDR::skyfilter() executable
    #             (median) the normal way for coarse fields [default]
    #             (min) suitable for crowded fields 
    #
    skymodel = median
    
    ############################################################################
    [offsets] 
    ############################################################################
    
    # single_point: If true, means that the SEextractor objmask will be reduced to a
    # single point (centroid) to run the cross-reference offset algorithm,i.e.,
    # each object is represented by a single, one-valued pixel, located at the
    # coordinates specified by its X_IMAGE and Y_IMAGE parameters in the
    # SExtractor catalog.
    # It is done mainly to avoid problems with large object masks (extended objtects,
    # satured objects, etc ..) that make the cross-reference algorithm too slow 
    # and even might with wrong results.  
    # 
    single_point = True
    
    # Object mask
    mask_minarea = 5 #15   # sex:DDETECT_MINAREA used for object masking
    mask_thresh = 1.5 #5.0   # sex:DDETECT_THRESH used for object masking
    
    satur_level = 3000000 # sex:SATUR_LEVEL; level (in ADUs) at which arises saturation
    
    #
    # Minimun overlap correlation fraction between offset translated images (from irdr::offset.c)
    #
    min_corr_frac = 0.001
    
    ############################################################################
    [astrometry]
    ############################################################################
    
    # Object mask
    mask_minarea = 5   # sex:DETECT_MINAREA used for object masking
    mask_maxarea = 10000 # Not yet implemented in PAPI !!! and not supported by SExtractor
    mask_thresh = 1.5   # sex:DETECT_THRESH used for object masking
    #expand_mask = 0.5   # amount to expand the object mask regions
    satur_level = 3000000 # sex:SATUR_LEVEL; level (in ADUs) at which arises saturation
    catalog = GSC-2.3    # Catalog used in SCAMP configuration (2MASS, USNO-A1, USNO-A2,
                         # USNO-B1,SC-1.3, GSC-2.2, GSC-2.3, UCAC-1, UCAC-2, UCAC-3, 
                         # NOMAD-1, PPMX, DENIS-3, SDSS-R3, SDSS-R5, SDSS-R6 or SDSS-R7)
    
    
    ############################################################################
    [keywords] 
    ############################################################################
    
    # The pipeline is designed for the PANIC data files. You should change
    # this options in case you were going to work with images whose keywords are
    # not the same.
    
    object_name = IMAGETYP    # Target description
    julian_date = MJD-OBS     # Modified Julian date
    x_size = NAXIS1           # Length of x-axis                       
    y_size = NAXIS2           # Length of y-axis
    ra = RA, CRVAL1           # Right ascension, in decimal degrees | The list defines the priority in which the values are read
    dec = DEC, CRVAL2         # Declination, in decimal degrees     | That is, if "DEC" is not found, CRVAL2 will be read, and so on.
    filter = FILTER           # Filter name
    
    
    ############################################################################
    [quicklook] 
    ############################################################################
    
    # Next are some configurable options for the PANIC Quick Look tool
    #
    # some important directories
    #
    #source = /data/O2K/Feb.2012/120213      # it can be a directory or a file (GEIRS datalog file)
    #source = /mnt/GEIRS_DATA
    #source = /home/panic/GEIRS/log/save_CA2.2m.log
    #source = /mnt/tmp/fitsfiles.corrected
    #source = /home/panic/tmp/fitsfiles.corrected
    source = /home/panic/data/i_test/002/
    output_dir = /data/out   # the directory to which the resulting images will be saved.
    temp_dir = /data/tmp     # the directory to which temporal results will be saved
    verbose = True
    
    # Run parameters
    run_mode = None # default (initial) run mode of the QL; it can be (None, Lazy, Prereduce)



Show grouped files in a row directory
-------------------------------------
Command::

    $papi.py -s /my/raw_data/directory -p

Show grouped files per filter and coordinates in a row directory 
----------------------------------------------------------------
Command::

    $papi.py -s /my/raw_data/directory -g filter -p 


Getting PAPI Data
*****************

The PAPI pipeline requires the full set of uncalibrated data products 
and best reference files for each observation in the input image set. These files 
can be readily obtained through the CAHA_ archive. When
requesting data from CAHA you need to specify:
   
   * Instrument : **PANIC**
   * Science Files Requested: **Uncalibrated - Raw** 
   * Reference Files: **Advanced Data Products**

.. image:: _static/caha_archive.jpg
   :align: center
   :height: 300 px
   :width: 565 px

.. _CAHA: http://caha.sdc.cab.inta-csic.es/calto/index.jsp

.. index:: options


Troubleshooting
***************

As we stated previously, PAPI was developed primarily for reducing NIR imaging
data of any kind of sources (galactic, extragalactic, coarse or crowed fields, 
and extended objects). Here are some tips for reducing each types of data:

* Coarse fields:
* Crowded fields:
* Extended objects:

*Add tips here*