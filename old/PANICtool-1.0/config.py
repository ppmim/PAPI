################################################################################
#
# PANIClook
#
# config.py
#
# Last update 27/02/2008
#
################################################################################

"""
   Provide default settings for the reduction process, and methods to
   operate on these settings.
"""

# Import external modules

import os
import cPickle
import configobj
import messageLog

################################################################################

class Setting: 

  """
     Defines one single configuration option, carrying a value, name, type,
     description and affiliated directory (if applicable)
  """

  def __init__(self, value="", type="", name="", desc="", indir=None):
    self.value    = value
    self.type     = type
    self.name     = name
    self.desc     = desc
    self.indir    = indir


  # Provide the value of this setting as it's default representation
  def __repr__(self):
    return repr(self.value)


  # Then the length should also come from the value field
  def __len__(self):
    return len(self.value)


  def nicevalue(self):
  
    "Provide a nice-looking value for this option if it is a list of files"

    # In a GUI you typically don't want all the files listed, but rather see
    # the number of files contained in the list
    
    # Is this a list of (FITS) files?
    if (self.type in ('filelist', 'fitsfilelist')):

      nfiles = len(self)
      if (nfiles  > 1): outvalue = '%i files selected' % nfiles
      # If there is only one file, show its name
      if (nfiles == 1): outvalue = self.value[0]
      # If there are no files, say so
      if (nfiles == 0): outvalue = 'No files selected'

    # Otherwise, the standard representation is good enough
    else: outvalue = self.value
    
    return outvalue


  def set(self, value):

    "Set the value of this option"
  
    self.value = value


  def get(self):

    "Get the value of this option"

    return self.__repr__()



def save_raw(outfile=None, optionlist=[], ignore_options=[]):

  "Save the current values of options to disk (raw format using cPickle)"

  # Only proceed if outfile is defined
  if outfile:

    # Create an object that will only contain the elements that are
    # supposed to be written to disk
    outoptions = dict()

    # If optionlist is defined, only select the listed items
    if optionlist: keys = optionlist
    # Otherwise, select all
    else: keys = options.keys()

    # If ignore_options was defined, modify the keylist
    if ignore_options: 
      for ignore_key in ignore_options:
        if ignore_key in keys: keys.remove(ignore_key)

    # Loop over options (well... the corresponding keys)
    for key in keys:
      # Copy the options into the output object
      try:
        outoptions[key] = options[key]
	messageLog.put("Wrote value for %s to %s" % (key, outfile), 9)
      except KeyError:
        messageLog.put_error("Cannot write configuration option '%s'." % key)

    # Now, open the output file in write mode
    filehandle = open(outfile, mode='w')

    # Do the actual writing to disk. cPickle takes care of creating
    # a textual representation of the contents of outoptions
    try:
      cPickle.dump(outoptions, filehandle)
    except Exception, errstr:
      messageLog.put_error("Could not save config to %s : %s" % (infile, errstr))
      raise

    # Close the output file
    filehandle.close()

    # Inform user
    messageLog.put("Wrote configuration to %s" % outfile, 5)


def load_raw(infile=None, optionlist=[], ignore_options=[]):

  "Restore the values of options from disk (Raw format using cPickle)"

  # Only proceed if infile is defined
  if infile:

    # Test if file exists
    if not os.path.exists(infile):
      messageLog.put_error("Requested config file does not exist: %s" % infile)
      return

    # Open the input file for reading
    filehandle = open(infile, mode='r')

    # Use cPickle to read the newoptions object from disk
    try:
      newoptions = cPickle.load(filehandle)
    except Exception, errstr:
      messageLog.put_error("Could not load config from %s : %s" % (infile, errstr))
      raise

    # Close the input file
    filehandle.close()

    # Get the list of keys to read:
    keys = newoptions.keys()

    # If ignore_options was defined, modify the keylist
    if ignore_options: 
      for ignore_key in ignore_options:
        if ignore_key in keys: keys.remove(ignore_key)

    # Loop over keylist (i.e. the options)
    for key in keys:
      # Copy the read options into the current set of options
      try:
	options[key] = newoptions[key]
	messageLog.put("Restored value for %s from %s" % (key, infile), 9)
      except KeyError:
          messageLog.put_error("Cannot restore configuration option '%s'." % key)

    # Inform user
    messageLog.put("Restored configuration from %s" % infile, 5)


def save(outfile=None, optionlist=[], ignore_options=[]):

  "Save the current values of options to disk"

  # Only proceed if outfile is defined
  if outfile:

    # Create an object that will only contain the elements that are
    # supposed to be written to disk
    outconfig = configobj.ConfigObj(options['default_savefile'].value)

    # If optionlist is defined, only select the listed items
    if optionlist: keys = optionlist
    # Otherwise, select all
    else: keys = outconfig['main'].keys()

    # If ignore_options was defined, modify the keylist
    if ignore_options: 
      for ignore_key in ignore_options:
        if ignore_key in keys:
	  keys.remove(ignore_key)
	  del(outconfig[key])

    # Loop over options
    for key in keys:

      try:
        outvalue = options[key].value
      except KeyError:
	messageLog.put_warning("Could not write value for %s - continuing" % key)
        return

      # Process arrays of files
      if options[key].type in ('filelist', 'fitsfilelist'):

        # Cast to array of strings
        if type(outvalue) != type([]): outvalue = [outvalue]

        # And attach directory name is necessary
	filelist = outvalue
	outvalue = []

	for filename in filelist:
	  if os.path.dirname(filename) == options[key].indir.value :
	    filename = os.path.basename(filename)
          outvalue.append(filename)

      if options[key].type in ('fitsfile', 'outfitsfile'):
	if os.path.dirname(outvalue) == options[key].indir.value :
	  outvalue = os.path.basename(outvalue)


      # Copy the options into the output object
      try:
        outconfig['main'][key] = outvalue
	messageLog.put("Wrote value for %s to %s" % (key, outfile), 9)
      except KeyError:
        messageLog.put_error("Cannot write configuration option '%s'." % key)

    try:
      # Now, open the output file in write mode
      filehandle = open(outfile, mode='w')
      # Write the data
      outconfig.write(filehandle)
      # Close the output file
      filehandle.close()
    except Exception, errstr:
      messageLog.put_error("Could not save config to %s : %s" % (outfile, errstr))

    # Inform user
    messageLog.put("Wrote configuration to %s" % outfile, 5)



def load(infile=None, optionlist=[], ignore_options=[]):

  "Restore the values of from disk"

  
  # Only proceed if infile is defined
  if infile:

    # Test if file exists
    if not os.path.exists(infile):
      raise IOError

    try:
      newoptions = configobj.ConfigObj(infile)
    except Exception, errstr:
      messageLog.put_error("Syntax error in configuration file : %s" % infile)
      messageLog.put_error("Error: %s " % errstr)
      return

    # Get the list of keys to read:
    try:
      keys = newoptions['main'].keys()
    except KeyError:
      messageLog.put_error("Syntax error in config file: no 'main' section defined")
      raise

    # If ignore_options was defined, modify the keylist
    if ignore_options: 
      for ignore_key in ignore_options:
	if ignore_key in keys: keys.remove(ignore_key)


    # Loop over keylist (i.e. the options) and process the different formats
    for key in keys:

      # Convert to integer if possible
      try:
	newvalue = newoptions['main'].as_int(key)
      except:
	newvalue = newoptions['main'][key]

	try:
	  options[key]
	except:
	  messageLog.put_error("Option %s in %s does not exist - ignored" % (key, infile))
	  return
	
	# Process arrays of files
	if options[key].type in ('filelist', 'fitsfilelist'):

          # Cast to array of strings
          if type(newvalue) != type([]): newvalue = [newvalue]

          # And attach directory name is necessary
	  filelist = newvalue
	  newvalue = []

	  for filename in filelist:
	    if not os.path.dirname(filename):
	      filename = (os.path.join(options[key].indir.value, filename))
            newvalue.append(filename)

	# Process arrays of files
	if options[key].type in ('fitsfile', 'outfitsfile'):

          # Cast to string
	  if type(newvalue) == type([]): newvalue = newvalue[0]

	  if newvalue != "":
	    if not os.path.dirname(newvalue):
	      newvalue = os.path.join(options[key].indir.value, newvalue)


      # Now copy the new values into the option
      try:
	options[key].value = newvalue
	messageLog.put("Restored value for %s from %s" % (key, infile), 9)
      except KeyError:
          messageLog.put_error("Cannot restore configuration option '%s'." % key)


    # Now, store the available reduction tasks
    options['reductionmodes'].value = {}
    options['availablemodes'].value = newoptions.keys()
    options['availablemodes'].value.remove('main')

    for mode in options['availablemodes'].value:
      options['reductionmodes'].value[mode] = []
      tasks = newoptions[mode].keys()
      for task in tasks:
	options['reductionmodes'].value[mode].append((task, newoptions[mode].as_bool(task)))
      messageLog.put("Read list of tasks for mode: %s" % mode, 9)

    # Inform user
    messageLog.put("Restored configuration from %s" % infile, 5)


def dump(optionlist=[]):

  """
     Create a quick dump of the current options in the logfile. Mainly for
     debugging purposes
  """

  # Only do this for the options given in optionlist
  if optionlist: keys = optionlist
  # ...or for all options if not list was defined
  else: keys = options.keys()

  # Output header
  messageLog.put("="*80)
  messageLog.put("Dumping current settings to logfile")
  messageLog.put("="*80)

  # Write one by one the name and values to the logfile
  for key in keys:
    try:
      messageLog.put("%s: %s" % (options[key].name, options[key].value))
    except KeyError:
      messageLog.put_error("Cannot find value of '%s'." % key)

  # Output footer
  messageLog.put("="*80)


###################################################################

# The items below should be self-explanatory. You can always define
# more options. Remember that new default option values may not
# necessarily be reflected in the program, because the first thing
# the program does is restore values from the file 'default.cfg',
# which may contain other values than those below.
#
# option.name	is the name that is displayed in the editing menus
# option.type	is the option type, this is used to determine which
#		function to use when editing its value
# option.indir	is the input directory associated with this option
#		(only relevant for files and directories)
# option.desc	is the text that will appear in pop-up help windows


# Create the options dictionary
options = dict()


################## CONSTRUCTORS for input options ###############

options['sextractor_dir']	= Setting('')
options['sextractor_dir'].name 	= 'Sextractor  directory'
options['sextractor_dir'].type	= 'directory'
options['sextractor_dir'].indir	= options['sextractor_dir']
options['sextractor_dir'].desc 	= """Directory where the SExtractor files are located.
                                  """

options['iraf_taskconfig_dir']		= Setting('')
options['iraf_taskconfig_dir'].name 	= 'IRAF tasks config directory'
options['iraf_taskconfig_dir'].type	= 'directory'
options['iraf_taskconfig_dir'].indir	= options['iraf_taskconfig_dir']
options['iraf_taskconfig_dir'].desc 	= """Directory where the IRAF task configuration files are located.
                            	          """

options['config_dir']		   = Setting('')
options['config_dir'].name	   = 'Config file directory'
options['config_dir'].type	   = 'directory'
options['config_dir'].indir	   = options['config_dir']
options['config_dir'].desc	   = """Directory where configuration files and pipeline configuration
                                             files are located.
                            		  """

options['inpath'] 		= Setting('')
options['inpath'].name		= 'Input directory'
options['inpath'].type		= 'directory'
options['inpath'].indir		= options['inpath']
options['inpath'].desc	 	= """Directory from where input data are taken. This is also the directory
				     that is monitored when doing automatic checking for new files.
                                  """

options['filename_filter']	= Setting('*.fit*')
options['filename_filter'].name	= 'Filename filter'
options['filename_filter'].type	= 'other'
options['filename_filter'].desc	= """The default filename pattern used for filtering when
				     autochecking, browsing and selecting files.
                                  """

options['biaslist'] 		= Setting([])
options['biaslist'].name	= 'List of BIAS frames'
options['biaslist'].type	= 'fitsfilelist'
options['biaslist'].indir	= options['inpath']
options['biaslist'].desc	= """A list of BIAS frames. These frames will be combined into the 'combined BIAS.'

                                  """
options['darklist'] 		= Setting([])
options['darklist'].name	= 'List of DARKS frames'
options['darklist'].type	= 'fitsfilelist'
options['darklist'].indir	= options['inpath']
options['darklist'].desc	= """A list of DARKS frames. These frames will be combined into the 'combined DARKS.'
                                  """


options['flatlist'] 		= Setting([])
options['flatlist'].name	= 'List of FLAT frames'
options['flatlist'].type	= 'fitsfilelist'
options['flatlist'].indir	= options['inpath']
options['flatlist'].desc	= """A list of FLAT frames. These frames will be combined into a 'combined FLAT.'
                                  """

options['sciencelist'] 		= Setting([])
options['sciencelist'].name	= 'List of SCIENCE frames'
options['sciencelist'].type	= 'fitsfilelist'
options['sciencelist'].indir	= options['inpath']
options['sciencelist'].desc	= """A list of SCIENCE frames. These frames will be reduced into a 'science ready frames'
                                  """


options['orderdef'] 		= Setting('')
options['orderdef'].name	= 'Order definition frame'
options['orderdef'].type	= 'fitsfile'
options['orderdef'].indir	= options['inpath']
options['orderdef'].desc	= """A well-exposed frame (esp. in the blue) that will be used to determine the
				     location of the spectral orders and the inter-order regions that are used
				     to estimate scattered light.
                                  """

options['interlacedorderdef']	       = Setting('')
options['interlacedorderdef'].name     = 'Interlaced order def. frame'
options['interlacedorderdef'].type     = 'fitsfile'
options['interlacedorderdef'].indir    = options['inpath']
options['interlacedorderdef'].desc     = """A well-exposed frame (esp. in the blue) that will be used to determine the
					     location of the interlaced calibration lamp spectra.
	                                  """

options['wavedef'] 			= Setting('')
options['wavedef'].name			= 'Wavelength definition frame'
options['wavedef'].type			= 'fitsfile'
options['wavedef'].indir		= options['inpath']
options['wavedef'].desc			= """A frame that can be used for determining the wavelength solution, normally
					     a well-exposed ThAr frame.
                        	          """

options['interlacedwavedef']	       = Setting('')
options['interlacedwavedef'].name      = 'Interlaced wavel. def. frame'
options['interlacedwavedef'].type      = 'fitsfile'
options['interlacedwavedef'].indir     = options['inpath']
options['interlacedwavedef'].desc      = """A frame that can be used for determining the interlaced wavelength solution,
					     normally a well-exposed ThAr frame.
	                                  """

################## CONSTRUCTORS for calibration output files ##############


options['refpath'] 		= Setting('')
options['refpath'].name		= 'Reference directory'
options['refpath'].type		= 'directory'
options['refpath'].indir	= options['refpath']
options['refpath'].desc		= """Directory where the calibration frames (combined BIAS, combined FLAT, 2D-normalized
				     combined FLAT, order definition reference and the wavelength reference)
				     are stored.
                                  """

options['pixelmask'] 		= Setting('')
options['pixelmask'].name	= 'Bad pixel mask'
options['pixelmask'].type	= 'fitsfile'
options['pixelmask'].indir	= options['refpath']
options['pixelmask'].desc	= """A frame containing the value of 1 for good pixels and 0 for bad pixels. Bad
				     pixels will be disregarded during the reduction.
                                  """

options['masterbias']		= Setting('')
options['masterbias'].name	= 'Combined BIAS frame'
options['masterbias'].type	= 'outfitsfile'
options['masterbias'].indir	= options['refpath']
options['masterbias'].desc	= """The combined BIAS frame gives the zero-level of the detector readout for each pixel.
                                  """
options['masterdark']		= Setting('')
options['masterdark'].name	= 'Combined DARKS frame'
options['masterdark'].type	= 'outfitsfile'
options['masterdark'].indir	= options['refpath']
options['masterdark'].desc	= """The combined DARKS frame gives the zero-level of the detector readout for each pixel.
                                  """


options['masterflat']		= Setting('')
options['masterflat'].name	= 'Combined FLAT frame'
options['masterflat'].type	= 'outfitsfile'
options['masterflat'].indir	= options['refpath']
options['masterflat'].desc	= """The combined FLAT frame gives the response of the  instrument and detector
				     as a function of pixel.
                                  """

options['masternormflat']	= Setting('')
options['masternormflat'].name	= 'Normalized combined FLAT'
options['masternormflat'].type	= 'outfitsfile'
options['masternormflat'].indir	= options['refpath']
options['masternormflat'].desc	= """The 2D-normalized combined FLAT frame is produced by dividing the combined FLAT frame
				     by a 2D-fit of the order shape. Dividing observed frames by the normalized
				     combined FLAT corrects for local changes in the detector and instrument
				     sensitivity, such as fringes.
                                  """

options['mastersky']		= Setting('')
options['mastersky'].name	= 'Combined science frame'
options['mastersky'].type	= 'outfitsfile'
options['mastersky'].indir	= options['refpath']
options['mastersky'].desc	= """The combined science frame gives an estimation of the background level for each pixel.
                                  """

options['masterdark_sc']	= Setting('')
options['masterdark_sc'].name	= 'Combined and scaled(divided) by EXPTIME'
options['masterdark_sc'].type	= 'outfitsfile'
options['masterdark_sc'].indir	= options['refpath']
options['masterdark_sc'].desc	= """The combined and scaled (divided) by EXPTIME master dark frame.
                                  """

options['blazeshape']		= Setting('')
options['blazeshape'].name	= 'Fitted blaze shape'
options['blazeshape'].type	= 'outfitsfile'
options['blazeshape'].indir	= options['refpath']
options['blazeshape'].desc	= """The is frame gives the blaze shape of orders extracted from the combined FLAT frame.
                                  """

options['masterorderdef']	= Setting('')
options['masterorderdef'].name	= 'Order definition reference'
options['masterorderdef'].type	= 'outfitsfile'
options['masterorderdef'].indir	= options['refpath']
options['masterorderdef'].desc	= """This frame contains the definition of the
				     location of the spectral orders. It is an
				     output frame, similar to the order
				     definition frame, except that overscan
				     regions are removed.
                                  """

options['masterinterlacedorderdef']	  = Setting('')
options['masterinterlacedorderdef'].name  = 'Interl. order def. reference'
options['masterinterlacedorderdef'].type  = 'outfitsfile'
options['masterinterlacedorderdef'].indir = options['refpath']
options['masterinterlacedorderdef'].desc  = """This frame contains the
					       definition of the location of
					       the interlaced spectral
					       orders.It is an output frame,
					       similar to the order definition
					       frame, except that overscan
					       regions are removed.
	                                    """

options['waveref']   = Setting('')
options['waveref'].name   = 'Wavelength reference frame'
options['waveref'].type   = 'outfitsfile'
options['waveref'].indir  = options['refpath']
options['waveref'].desc   = """This frame contains the wavelength
				     solution. It is an output frame, similar
				     to the wavelength definition frame, except
				     that overscan regions are removed.
                                  """

options['interlacedwaveref']	  = Setting('')
options['interlacedwaveref'].name = 'Interlaced wvl. ref. frame'
options['interlacedwaveref'].type = 'outfitsfile'
options['interlacedwaveref'].indir	  = options['refpath']
options['interlacedwaveref'].desc = """This frame contains the interlaced
					     wavelength solution. It is an output
					     frame, similar to the wavelength
					     definition frame, except that overscan
					     regions are removed.
	                                  """

options['masterwaveref']	= Setting('')
options['masterwaveref'].name	= 'Master wavel. ref. frame'
options['masterwaveref'].type	= 'outfitsfile'
options['masterwaveref'].indir	= options['refpath']
options['masterwaveref'].desc	= """This frame contains a standard wavelength
				     solution that can be used as a first guess
				     when determining new wavelength solutions.
                                  """


################## CONSTRUCTORS for output options ###############

options['outpath']		= Setting('')
options['outpath'].name		= 'Output directory'
options['outpath'].type		= 'directory'
options['outpath'].indir	= options['outpath']
options['outpath'].desc	 	= """Directory where reduced frames are stored.
                                  """

################## CONSTRUCTORS for other options ###############

options['fitsheaders']		= Setting('object exptime')
options['fitsheaders'].name	= 'Listed FITS headers'
options['fitsheaders'].type	= 'other'
options['fitsheaders'].desc	= """These FITS headers will be listed in file selection dialogs. Multiple headers
				     may be listed, separated by whitespace.
                                  """

options['mef_dataext']		= Setting(1)
options['mef_dataext'].name	= 'MEF data extension'
options['mef_dataext'].type	= 'integer'
options['mef_dataext'].desc	= """The number of the Multi-Extension FITS Header Data Unit that contains the echelle frame,
                                     or 0 if no FITS extensions are used.
                                  """

options['frameorientation']	       = Setting(7)
options['frameorientation'].name       = 'Frame Orientation'
options['frameorientation'].type       = 'integer'
options['frameorientation'].desc       = """The orientation of the data frame. For a proper reduction the frame should
                                            have increasing wavelength along the X and Y axes. The value here is used
					    to rotate and flip the data frame before processing in order to obtain the
					    correct orientation. A value of 0 means no processsing, 1-3 a rotation of 90,
					    180 or 270 degrees, and 4-7 are similar to 0-3, but with an additional transpose
					    of the data.
	                                 """

options['x_start']		= Setting(51)
options['x_start'].name		= 'x_start'
options['x_start'].type		= 'integer'
options['x_start'].desc		= """The first data pixel along the x-axis of a frame. Used when clipping the pre-
				     and overscan regions.
                                  """

options['x_start']		= Setting(51)
options['x_start'].name		= 'x_start'
options['x_start'].type		= 'integer'
options['x_start'].desc		= """The first data pixel along the x-axis of a frame. Used when clipping the pre-
				     and overscan regions.
                                  """

options['x_end']		= Setting(2098)
options['x_end'].name		= 'x_end'
options['x_end'].type		= 'integer'
options['x_end'].desc		= """The last data pixel along the x-axis of a frame. Used when clipping the pre-
				     and overscan regions.
                                  """

options['y_start']		= Setting(3)
options['y_start'].name		= 'y_start'
options['y_start'].type		= 'integer'
options['y_start'].desc		= """The first data pixel along the y-axis of a frame. Used when clipping the pre-
				     and overscan regions.
                                  """

options['y_end']		= Setting(2048)
options['y_end'].name		= 'y_end'
options['y_end'].type		= 'integer'
options['y_end'].desc		= """The last data pixel along the y-axis of a frame. Used when clipping the pre-
				     and overscan regions.
                                  """


options['plot_startwave']		= Setting(3500)
options['plot_startwave'].name		= 'Starting wavelength in plot'
options['plot_startwave'].type		= 'integer'
options['plot_startwave'].desc		= """The starting wavelength used in plots generated
					     with the Biggles plotting package.
	                                  """

options['plot_endwave']		= Setting(7500)
options['plot_endwave'].name	= 'Ending wavelength in plot'
options['plot_endwave'].type	= 'integer'
options['plot_endwave'].desc	= """The ending wavelength used in plots generat
                                      with the Biggles plotting package.
                                   """

options['plot_defaultorder']		= Setting(1)
options['plot_defaultorder'].name		= 'Default order to plot'
options['plot_defaultorder'].type		= 'integer'
options['plot_defaultorder'].desc		= """Default spectral order to plot if no
						     wavelengths are available.
	        	                          """



################## CONSTRUCTORS for internal options ###############

options['currentmode']	= Setting(())
options['currentmode'].name	= 'Current reduction mode'
options['currentmode'].type	= 'other'
options['currentmode'].desc	= """Internal variable used to hold the current reduction mode.
                                  """

options['autoloadrules']	= Setting(())
options['autoloadrules'].name	= 'Autoloading rules'
options['autoloadrules'].type	= 'other'
options['autoloadrules'].desc	= """Internal variable used to hold the definition of the AutoLoading rules.
                                  """

options['default_savefile']	= Setting('default.cfg')
options['default_savefile'].name	= 'Default configuration file'
options['default_savefile'].type	= 'other'
options['default_savefile'].desc	= """Internal variable used to hold the default configuration file
                                  """

options['default_autoload_savefile']	= Setting('default_autoload.cfg')
options['default_autoload_savefile'].name	= 'Default autoloader configuration file'
options['default_autoload_savefile'].type	= 'other'
options['default_autoload_savefile'].desc	= """Internal variable used to hold the default autoloader configuration file
                                  """

options['availablemodes']	= Setting({})
options['availablemodes'].name	= 'Available reduction modes'
options['availablemodes'].type	= 'other'
options['availablemodes'].desc	= """Internal variable used to hold the available reduction modes.
                                  """

options['reductionmodes']	= Setting({})
options['reductionmodes'].name	= 'Detailed reduction modes and task informations'
options['reductionmodes'].type	= 'other'
options['reductionmodes'].desc	= """Internal variable used to hold which tasks are part of which
				     reduction mode, and if they are selected or not.
                                  """
