#! /usr/bin/env python

# Copyright (c) 2009 Victor Terron. All rights reserved.
# Institute of Astrophysics of Andalusia, IAA-CSIC
#
# This file is part of PAPI.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


# Adapted to PAPI by jmiguel@iaa.es 29-Nov-2010

import sys
import ConfigParser
import os.path
import tempfile

# PAPI modules
import style


def default_config_file():
    """ Return the name of the default configuration file.

    It is __very__ __important__ to keep in mind that if the current value of
    the configuration file ("papi.cfg") is modified, the file
    ./irdr_panic/src/fitsIO/fitsIO.c has to be updated and recompiled.
    The reason for this is that the IRDR offsets code needs to read from it 
    the name of the needed keywords in the FITS images.
    
    """

    return "./config_files/papi_portatil.cfg"


def check_required_option(option_name, option_value):
    """ Abort the execution in the absence of a required option.

    The method prints an error message to standard output and aborts the 
    execution if 'nothing' was read for the required option. That is, the
    error is raised if the option value is None of has length zero.
    Otherwise, the method does nothing.

    Keyword arguments:
    option_name - the name of the required option.
    option_value - the value of the required option.

    """

    if not option_value:
         print style.prefix() + "The option '" + option_name + "' cannot be left empty"
         sys.exit(style.error_exit_message())



def read_parameter(config_object, section, parameter, type, required = False,
                   config_file = default_config_file()):
    """ Read a parameter from the configuration file.

    The method attemps to read the specified parameter in the section in which
    it is included, returning its value in the configuration file. If an error
    is found (such as not being able to find the section), the appropiate (and
    intendendly useful) error message will be thrown and the execution aborted.
    
    Keyword arguments:
    config_object - the object of the configuration file.
    section - name of the section to which the parameter belongs.
    parameter - name of the parameter whose value is to be retrieved.
    type - the data type of the parameter. The expected values are int, float,
           bool or str. The parameter will be assumed to be str if a different
           data type is specified.
    required - if set to True, an error will be thrown if the parameter was
               left empty. Otherwise, None will be returned when this happens.
    config_file - name of the configuration file, used only in order to display
                  a more helpful error message before aborting the execution.

    """

    try:    
        # Return None if nothing was read and the parameter was not mandatory
        if not config_object.get(section, parameter):
            if not required:
                return None
            else:
                print style.prefix() + "[" + config_file + "] " + \
                "The parameter '" + parameter + "' in section '" + section + \
                "' cannot be left empty."
                sys.exit(style.error_exit_message())    


        if type == int:    
            value = config_object.getint(section, parameter)
        elif type == float:
            value = config_object.getfloat(section, parameter)
        elif type == bool:
            value = config_object.getboolean(section, parameter)
        else:
            value = config_object.get(section, parameter)

        return value
        
    # Is an exception is raised, show error message and abort execution
    except ValueError:
        print style.prefix() + "Error while parsing parameter '" + parameter + "' in " \
              "section '" + section + "'."
        print style.prefix() + "Is it a " + str(type) + " data type?"
        sys.exit(style.error_exit_message())
    except ConfigParser.NoSectionError:
        print style.prefix() + "The section '" + section + "' could not be found."
        sys.exit(style.error_exit_message())
    except ConfigParser.NoOptionError:
        if  not required: return None
        else:
            print style.prefix() + "The parameter '" + parameter + "' could not be found in section '" + section + "'."
            sys.exit(style.error_exit_message())
    except ConfigParser.Error: 
        print style.prefix() + "An error occurred while parsing parameter '" + parameter +"."
        sys.exit(style.error_exit_message())    


def read_file_parameter(config_object, section, parameter, \
                        config_file = default_config_file()):
    """ Read a file path, specified in the config file, and check its existence.

    The method works exactly as read_parameter(), with the addition that it
    also checks that the read parameter refers to an existing file. If that is
    not the case, an error message will be displayed and the execution aborted.

    Keyword arguments:
    config_object - the object of the configuration file.
    section - name of the section to which the parameter belongs.
    parameter - name of the parameter whose value is to be retrieved.
    type - the data type of the parameter. The expected values are int, float,
           bool or str. The parameter will be assumed to be str if a different
           data type is specified.
    required - if set to True, an error will be thrown if the parameter was
               left empty. Otherwise, None will be returned when this happens.
    config_file - name of the configuration file, used only in order to display
                  a more helpful error message before aborting the execution.

    """

    file_path = read_parameter(config_object, section, parameter, str, True, config_file)
    if not os.path.exists(file_path):
        print style.prefix() + "The path '" + file_path + "' specified in parameter '" + parameter + "',"
        print style.prefix() + "section '" + section + "' does not exist."
        sys.exit(style.error_exit_message())
    elif not os.path.isfile(file_path):
        print style.prefix() + "The path '" + file_path + "' specified in parameter '" + parameter + "',"
        print style.prefix() + "section '" + section + "' does not refer to an existing file."
        sys.exit(style.error_exit_message())
    else:
        return file_path


def parse_list_of_strings(strings_list):
    """ Extract the substrings from a comma-delimited string. 

    The method receives a string and returns a list of the substrings, 
    separated by commas, that it contains. For exampe, the string "abc,de,fg"
    would return ["abc", "de, "fg"].

    Keyword arguments:
    strings_list - the string representation of a list of strings.

    """
    
    splitted_line = strings_list.split(",")
    # Remove (possible) leading and trailing whitespaces
    for index in range(len(splitted_line)):
        splitted_line[index] = splitted_line[index].strip()

    return splitted_line


def read_list_of_strings(config_object, section, parameter, required = False):
    """ Read a parameter comprising one or more comma-delimited strings.
    
    The method reads from the configuration file a parameter that is expected
    to consist of at least one string. Additional strings, separated by commas,
    may be specified. A list containing the parsed strings is returned. None
    will be returned if for a non-required parameters none was provided.

    Keyword arguments:
    config_object - the object of the configuration file.
    section - name of the section to which the parameter belongs.
    parameter - name of the parameter whose value is to be retrieved.
    required - if set to True, an error will be thrown if for the parameter at
               least one element was specified. Otherwise, None will be returned
               when this happens.
        
    """

    try:
        line = config_object.get(section, parameter)
        
    except ConfigParser.NoSectionError:
        print style.prefix() + "The section '" + section + "' could not be found."
        sys.exit(style.error_exit_message())
    except ConfigParser.NoOptionError:
        print style.prefix() + "The parameter '" + parameter + "' could not be found in " \
              "section '" + section + "'."
        sys.exit(style.error_exit_message())
    except ConfigParser.Error: 
        print style.prefix() + "An error occurred while parsing parameter '" + parameter +"."
        sys.exit(style.error_exit_message())

    if not line:    # if empty list
        if required:
            print style.prefix() + "[" + config_file + "] " + \
                  "The parameter '" + parameter + "' in section '" + section + \
                  "' cannot be left empty."
            sys.exit(style.error_exit_message())    
        else:
            return None
    else:
        return parse_list_of_strings(line)
   


def parse_str_interval(string_interval, section = None, parameter = None, \
                       separator = "-"):
    """ Extract the values from the string representation of an interval.

    The method receives the string representation of an interval, defined as
    two integers separated by the "-" character, and returns a tuple consisting
    of the two values, casted to integers. That is, for the string "1-56", for
    example, the tuple (1, 56) is returned. An error is shown and the execution
    aborted is an invalid interval is provided. If both section and parameter
    are not None, the error message will provide additional information, 
    indicating exactly where the error is in the configuration file.

    Keyword arguments:
    string_interval: the string representation of the interval, expected to
                     consist of two integers separated by the "-" character.
    section - name of the section to which the parameter in which the interval
              was defined belongs. This argument is required only in order to
              display helpful error messages.
    parameter - name of the parameter in which the interval was defined. This
                argument is required only in order to display helpful error 
                messages.
    separator - the delimiter string.

    """

    splitted_interval = string_interval.split(separator)
   
    # (a) The two integers must be separated by the "-" character. Then, after
    # splitting the string representation of the interval, using "-" as the
    # delimiter string, there should be two elements. We check that '' is not
    # in the returned list because "1-" is returned by split() as ['1', '']
    if len(splitted_interval) != 2 or '' in splitted_interval:

        print style.prefix() + "The interval '" + string_interval + "'",

        if section is not None and parameter is not None:
            print "in section '" + section + "',", \
                  "parameter '" + parameter + ","

            # Does not append a new line or space (as print does if , is used)
            sys.stdout.write(style.prefix())  
        print "does not define two integers separated by the '-' character."
        sys.exit(style.error_exit_message())

    
    # (b) Both values must be castable to integer (that is, must _be_ integers)
    casted_values = []
    for value in splitted_interval:
        try:
            casted_values.append(int(value))
        except:
            print style.prefix() + "The value '" + value + "' specified in" \
                  " the interval '" + string_interval + "',",
            if section is not None and parameter is not None:
                print "parameter '" + parameter + "',"
                print style.prefix() + "section '" + section + "',",
            print "is not an integer."
            sys.exit(style.error_exit_message())

    # Return both values as a tuple
    return casted_values[0], casted_values[1]


def parse_list_of_intervals(string_list_of_intervals, section = None, \
                            parameter = None, separator = "-"):
    """ Extract a list of intervals from its string representation.
    
    The method receives a string paramater that is expected to represent,
    comma separated, a series of intervals, which are defined as two integers
    separated bythe "-" character. A list contained the parsed intervals, each
    one now a tuple of the two specified values, casted to integer, is returned.
    For example, for the string "1-3, 4-6", the returned value would be
    [(1,3), (4,6)]. An exception to this rule is the character "*": if present
    in the list of intervals, it will also included in the returned list.
    That is, for "1-3, 4,6", *" the list [(1,3), (4,6), "*"] is returned.

    If both section and parameter are not None, error messages, if any, will
    provide additional information, indicating exactly where the error is
    located in the configuration file.

    None will be returned if the parameters list is empty. 

    Keyword arguments:   
    string_list_of_intervals - the string representation of the intervals
    section - name of the section to which the parameter belongs.
    parameter - name of the parameter whose value is to be retrieved.
    separator - the delimiter string.    

    """

    if not string_list_of_intervals:
        return None
    else:
        parsed_intervals = []
        for str_interval in string_list_of_intervals:
            # The '*' character is not parsed
            if "*" == str_interval:
                parsed_intervals.append("*")
            # Ignore strings of length zero
            elif len(str_interval) != 0:
                parsed_intervals.append(parse_str_interval(str_interval, section, parameter, separator))
        return parsed_intervals
    

def read_list_of_intervals(config_object, section, parameter, separator = "-"):
    """ Read a parameter comprising one or more comma-delimited intervals.

    The method is just a wrapper for the previous method,
    parse_list_of_intervals, so that lists of intervals can be directly read
    from the configuration method. For further information, you should refer to
    the documentation of that method.

    Keyword arguments:
    config_object - the object of the configuration file.
    section - name of the section to which the parameter belongs.
    parameter - name of the parameter whose value is to be retrieved.
    separator - the delimiter string.   

    """

    intervals = read_list_of_strings(config_object, section, parameter)
    return parse_list_of_intervals(intervals, section, parameter)
           

class CommentlessFile(file):
    """ Implements a commentless file subclass.

    ConfigParser forces comments to start on their only line (that is, only
    line beginnings with '#' or ';' are ignored), but in the configuration file
    we would like to be able to add comments on the same line than an option.
    Therefore, we are subclassing the file class, so that everything after the
    first "#" character is discarded.

    """

    def readline(self):
        """ Read a single commentless line from the file. """

        line = super(CommentlessFile, self).readline()
        if line:
            line = line.split('#', 1)[0].strip()
            return line + '\n'
        else:
            return ''


def read_config_file(config_file = default_config_file()):
    """ Read the PAPI configuration file.

    Reads the configuration file for PAPI. The method returns a dictionary,
    each one of whose entries maps the name of a section in the file
    (e.g., "general" or "overscan") to a (second) dictionary which contains
    each one of the parameters in the section and the value that was specified
    in the configuration file. If one or more section or parameters are
    invalid, missing or have been set to an invalid value, the appropiate (and
    intendedly useful) error message will be thrown and the execution aborted.

    Keyword arguments:
    config_file - path to the configuration file.

    """
   
    config = ConfigParser.SafeConfigParser()
    config.optionxform = str # make ConfigParser case sensitive

    try:
        config.readfp(CommentlessFile(config_file))
    except IOError:
        print style.prefix() + "The configuration file '" + config_file + \
              "' could not be read."
        sys.exit(style.error_exit_message())    

    options = {}     # empty dictionaryl

    ##################### "general" section ##################################
    ## ACTUALIZAR LLAMADAS a read_parameter()
    general = {} 
 
    #general["config_dir"] = read_parameter(config, "general", "config_dir", str, True, config_file)
    
    general["source"] = read_parameter(config, "general", "source", str, True, config_file)
    general["output_dir"] = read_parameter(config, "general", "output_dir", str, True, config_file)
    general["temp_dir"] = read_parameter(config, "general", "temp_dir", str, True, config_file)
    general["output_file"] = read_parameter(config, "general", "output_file", str, True, config_file)
    
    general["master_dark"] = read_parameter(config, "general", "master_dark", str, False, config_file)
    general["master_flat"] = read_parameter(config, "general", "master_flat", str, False, config_file)
    general["master_bpm"] = read_parameter(config, "general", "master_bpm", str, False, config_file)
    
    general["obs_mode"] = read_parameter(config, "general", "obs_mode", str, False, config_file)
    general["reduction_mode"] = read_parameter(config, "general", "reduction_mode", str, False, config_file)
    general["check_data"] = read_parameter(config, "general", "check_data", bool, False, config_file)
    general["group_by"] = read_parameter(config, "general", "group_by", str, False, config_file)
    
    general["max_mjd_diff"] = read_parameter(config, "general", "max_mjd_diff", float, False, config_file)
    
    
    general["apply_dark_flat"] = read_parameter(config, "general", "apply_dark_flat", bool, True, config_file)
    general["scale"] = read_parameter(config, "general", "scale", float, True, config_file)
    general["equinox"] = read_parameter(config, "general", "equinox", int, True, config_file)
    general["radecsys"] = read_parameter(config, "general", "radecsys", str, True, config_file)
    general["pattern"] = read_parameter(config, "general", "pattern", str, False, config_file)
    general["parallel"] = read_parameter(config, "general", "parallel", bool, True, config_file)
    general["ncpus"] = read_parameter(config, "general", "ncpus", int, True, config_file)
    general["verbose"] = read_parameter(config, "general", "verbose", bool, False, config_file)
    
    filter_prefix = "filter_name_"
    for filter_letter in ["Z", "Y", "J", "H", "K", "Ks"]:
        general[filter_prefix + filter_letter] = read_list_of_strings(config, "general", filter_prefix + filter_letter)

        # At least one string must be specified, comma-separated
        if not general[filter_prefix + filter_letter]:
            print style.prefix() + "In section 'general', the parameter '" + filter_prefix + filter_letter + "' must list at least one string."
            sys.exit(style.error_exit_message())

    options["general"] = general

    ###################### Astromatic config files #####################################
    config_files = {}
    config_files["terapix_bin"] = read_parameter(config, "config_files", "terapix_bin", str, False, config_file)
    config_files["irdr_bin"] = read_parameter(config, "config_files", "irdr_bin", str, False, config_file)
    config_files["sextractor_conf"] = read_file_parameter(config, "config_files", "sextractor_conf")
    config_files["sextractor_param"] = read_file_parameter(config, "config_files", "sextractor_param")
    config_files["sextractor_nnw"] = read_file_parameter(config, "config_files", "sextractor_nnw")
    config_files["sextractor_conv"] = read_file_parameter(config, "config_files", "sextractor_conv")
    config_files["scamp_conf"] = read_file_parameter(config, "config_files", "scamp_conf")
    config_files["swarp_conf"] = read_file_parameter(config, "config_files", "swarp_conf")
    options["config_files"] = config_files


    ######################### "dark" section ##################################

    dark = {}
    
    dark["object_names"] = read_list_of_strings(config, "dark", "object_names")

    # At least one string must be specified, comma-separated
    if not dark["object_names"]:
        print style.prefix() + "In section 'dark', the parameter 'object_names' must list at least one string."
        sys.exit(style.error_exit_message())


    dark["suffix"] = read_parameter(config, "dark", "suffix", str, False, config_file)

    dark["check_prop"] = read_parameter(config, "dark", "check_prop", bool, False, config_file)

    dark["min_frames"] = read_parameter(config, "dark", "min_frames", int, False, config_file)
    

    options["dark"] = dark


    ########################## "dflats" section ################################
    dflats = {}
    
    dflats["object_names"] = read_list_of_strings(config, "dflats", "object_names")
    
    dflats["suffix"] = read_parameter(config, "dflats", "suffix", str, False, config_file)

    dflats["check_prop"] = read_parameter(config, "dflats", "check_prop", bool, False, config_file)
    
    dflats["min_frames"] = read_parameter(config, "dflats", "min_frames", int, False, config_file)
    
    area_width = read_parameter(config, "dflats", "area_width", int, True, config_file)
    if not area_width > 1:
        print style.prefix() + "[" + config_file + "] The value of 'area_width' in section 'dflats' must be a positive integer."
        sys.exit(style.error_exit_message())    
    else:
        dflats["area_width"] = area_width

    options["dflats"] = dflats

    ########################## "twflats" section ################################
    twflats = {}
    
    twflats["object_names"] = read_list_of_strings(config, "twflats", "object_names")
    
    twflats["suffix"] = read_parameter(config, "twflats", "suffix", str, False, config_file)

    twflats["check_prop"] = read_parameter(config, "twflats", "check_prop", bool, False, config_file)
    
    twflats["min_frames"] = read_parameter(config, "twflats", "min_frames", int, False, config_file)
    
    area_width = read_parameter(config, "twflats", "area_width", int, True, config_file)
    if not area_width > 1:
        print style.prefix() + "[" + config_file + "] The value of 'area_width' in section 'twflats' must be a positive integer."
        sys.exit(style.error_exit_message())    
    else:
        twflats["area_width"] = area_width

    options["twflats"] = twflats  

    
    ########################## "skysub" section ################################
    skysub = {}
    
    skysub["object_names"] = read_list_of_strings(config, "skysub", "object_names")
    skysub["check_prop"] = read_parameter(config, "skysub", "check_prop", bool, False, config_file)
    skysub["suffix"] = read_parameter(config, "skysub", "suffix", str, False, config_file)
    skysub["min_frames"] = read_parameter(config, "skysub", "min_frames", int, False, config_file)
    skysub["hwidth"] = read_parameter(config, "skysub", "hwidth", int, False, config_file)
    skysub["mask_minarea"] = read_parameter(config, "skysub", "mask_minarea", int, False, config_file)
    skysub["mask_thresh"] = read_parameter(config, "skysub", "mask_thresh", float, False, config_file)
    skysub["satur_level"] = read_parameter(config, "skysub", "satur_level", long, False, config_file)
    
    area_width = read_parameter(config, "skysub", "area_width", int, True, config_file)
    if not area_width > 1:
        print style.prefix() + "[" + config_file + "] The value of 'area_width' in section 'skysub' must be a positive integer."
        sys.exit(style.error_exit_message())    
    else:
        skysub["area_width"] = area_width
        
    options["skysub"] = skysub     
    
    ########################## "offsets" section ################################
    offsets = {}
    
    offsets["mask_minarea"] = read_parameter(config, "offsets", "mask_minarea", int, False, config_file)
    offsets["mask_thresh"] = read_parameter(config, "offsets", "mask_thresh", float, False, config_file)
    offsets["min_corr_frac"] = read_parameter(config, "offsets", "min_corr_frac", float, False, config_file)
    offsets["satur_level"] = read_parameter(config, "offsets", "satur_level", long, False, config_file)
    
    options["offsets"] = offsets     
    
    ########################## "skysub" section ################################
    gainmap = {}
    
    gainmap["object_names"] = read_list_of_strings(config, "gainmap", "object_names")
    gainmap["mingain"] = read_parameter(config, "gainmap", "mingain", float, True, config_file)
    gainmap["maxgain"] = read_parameter(config, "gainmap", "maxgain", float, True, config_file)
    gainmap["nxblock"] = read_parameter(config, "gainmap", "nxblock", int, True, config_file)
    gainmap["nyblock"] = read_parameter(config, "gainmap", "nyblock", int, True, config_file)
    gainmap["nsigma"] = read_parameter(config, "gainmap", "nsigma", float, True, config_file)
    
    area_width = read_parameter(config, "gainmap", "area_width", int, True, config_file)
    if not area_width > 1:
        print style.prefix() + "[" + config_file + "] The value of 'area_width' in section 'gainmap' must be a positive integer."
        sys.exit(style.error_exit_message())    
    else:
        gainmap["area_width"] = area_width
        
    options["gainmap"] = gainmap    
    
    ####################### "fits" section #################################
    ## ACTUALIZAR LLAMADAS a read_parameter()
    # This options, although loaded from the configuration file, are not used
    # by our Python code, but instead by the IRDR offsets.c code. However, we
    # are loading them here in order to check that all the required parameters
    # are present, so that errors can be detected early.

    fits = {}

    cannot_be_left_empty = ["right_ascension", "declination", "posang"]
    may_be_left_empty = ["scale_arcsec_per_pix", "pixel_sixe_in_x_microns",
                         "pixel_sixe_in_y_microns", "focus_scale_mm"]

    for option_name in cannot_be_left_empty + may_be_left_empty:
        fits[option_name] = read_parameter(config, "fits", option_name, str)

        # Abort execution if a required option was not defined
        if option_name in cannot_be_left_empty and not fits[option_name]:
            print style.prefix() + "In section 'fits', the value of '" + \
                  option_name + "' cannot be left empty."
            sys.exit(style.error_exit_message())

    # If no value was given for the option 'scale_arcsec_per_pix', the options
    # (1) 'pixel_sixe_in_x_microns', (2) 'pixel_sixe_in_y_microns' and 
    # (3) 'focus_scale_mm' will be used instead, so they must have been
    # defined in the configuration file.
    if not fits["scale_arcsec_per_pix"] and \
       (not fits["pixel_sixe_in_x_microns"] or \
        not fits["pixel_sixe_in_y_microns"] or \
        not fits["focus_scale_mm"]):
        print style.prefix() + "In section 'fits', if the value of " \
              "'scale_arcsec_per_pix' is left empty\n" + style.prefix() + \
              "the options 'pixel_sixe_in_x_microns', " \
              "'pixel_sixe_in_y_microns' and\n" + style.prefix() + \
              "'focus_scale_mm' must be defined."
        sys.exit(style.error_exit_message())

    options["fits"] = fits



    #################### FITS Keywords #####################
    keywords = {}
    keywords["object_name"] = read_parameter(config, "keywords", "object_name", str, True, config_file)
    keywords["julian_date"] = read_parameter(config, "keywords", "julian_date", str, True, config_file)
    keywords["file_creation"] = read_parameter(config, "keywords", "file_creation", str, True, config_file)    
    keywords["x_size"] = read_parameter(config, "keywords", "x_size", str, True, config_file)        
    keywords["y_size"] = read_parameter(config, "keywords", "y_size", str, True, config_file)        
    keywords["ra"] = read_parameter(config, "keywords", "ra", str, True, config_file)      
    # The keywords for alpha and delta are read as a list of strings.
    # There will (should) be at least one element in the list.        
    keywords["ra"] = read_list_of_strings(config, "keywords", "ra", True)
    keywords["dec"] = read_list_of_strings(config, "keywords", "dec", True)            

    #keywords["crpix1"] = read_parameter(config, "keywords", "crpix1", str, True, config_file)
    #keywords["crval1"] = read_parameter(config, "keywords", "crval1", str, True, config_file)        
    #keywords["crpix2"] = read_parameter(config, "keywords", "crpix2", str, True, config_file)
    #keywords["crval2"] = read_parameter(config, "keywords", "crval2", str, True, config_file)        
    keywords["filter"] = read_parameter(config, "keywords", "filter", str, True, config_file) 


    options["keywords"] = keywords

    
    #################### Quick-Look tool #####################
    quicklook = {}
    quicklook["source"] = read_parameter(config, "quicklook", "source", str, True, config_file)
    quicklook["output_dir"] = read_parameter(config, "quicklook", "output_dir", str, True, config_file)
    quicklook["temp_dir"] = read_parameter(config, "quicklook", "temp_dir", str, True, config_file)    
    quicklook["run_mode"] = read_parameter(config, "quicklook", "run_mode", str, True, config_file)        

    options["quicklook"] = quicklook

    #########################################################################
    # TODO: complete as the config file grows!

    return options        
    
    
    

def read_options(options, section, config_file = default_config_file()):
    """ Read the PAPI configuration file, overwriting some of the options.

    Parameters needed by the modules of the pipeline have not only to be 
    checked for existence or type, but also tested in order to ensure that they
    satisfy some criterion. For example, a 'threshold' parameter, which defines
    the maximum allowed relative percent difference between two values, must 
    have a value in the range [0,1]. These are tests performed by the method
    config.read_config_file(), but the problem arises when some of the
    parameters defined in the configuration file are overwritten after being
    passed to a module as an option. Then, until this method was written, these
    values had to be checked _again_, which was nothing but a painful,
    unnecessary duplication of code, and which can be avoided thanks to this
    routine.

    What this method does is to parse the configuration file, as it would 
    be normally done, but _after_ having overwritten the parameters of the 
    'section' section that were passed as an option to the program. Note that
    the configuraton file remains untouched, as what is actually modified and
    parsed is a temporary copy of it, which is automatically deleted when it is
    no longer needed. The process is, in this manner, transparent to the user:
    if the paramater of 'section' was passed as an option, the value specified
    in the configuration file is ignored.

    The method returns a dictionary, each one of whose entries maps the name of
    a section in the file (e.g., 'general' or 'overscan') to a (second)
    dictionary which contains each one of the parameters in the section and the
    corresponding value. If one or more section or parameters are invalid,
    missing or have been set to an invalid value, the appropiate (and intendedly
    useful) error message will be printed and the execution aborted.

    Keyword arguments:
    options - an object containing values for all the options that were passed
              to the program. For example, if '--file' takes a single string
              argument, then options.file would be the filename supplied by the
              user, or None if the user did not supply that option. Note that
              this is the first value of the two-element tuple returned by
              the method optparse.parse_args().
    section - section of the configuration file whose values will be overwritten
              by those options specified in 'options'.
    config_file - path to the configuration file.
    """

    # Make a temporary copy of the configuration file
    temp_config = ConfigParser.SafeConfigParser()
    temp_config.optionxform = str # make ConfigParser case sensitive

    if not temp_config.read(config_file):
        print style.prefix() + "The configuration file '%s' could not be read." % config_file
        sys.exit(style.error_exit_message())    

    # Update in the copy of the configuration file those parameters which are
    # not None or empty (that is, is the option suffix was defined but left
    # empty (i.e., --sufix=''), it is ignored and does not overwrites the value
    # specified in the configuration file.

    for option in options.__dict__:  # get the object's symbol table.
        value = options.__dict__[option]
        if value is not None:
            temp_config.set(section, option, str(value))

    # Write the modified copy of the configuration file to disk        
    temp_fd, temp_filename = tempfile.mkstemp(suffix='.cfg')
    os.close(temp_fd)   # but ConfigParser needs a file object...
    fd = open(temp_filename, 'wt')
    temp_config.write(fd)
    fd.close()    
    # Finally, parse the modified configuration file, as we would normally do
    options = read_config_file(temp_filename)
    print "FILE=",temp_filename
    #os.remove(temp_filename)
    return options

    # TODO: the above method has not still been thoroughly tested!!!!!!!

if __name__ == "__main__":

    print read_config_file()
        
    

