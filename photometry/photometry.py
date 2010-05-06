#! /usr/bin/env python

# Copyright (c) 2009 Victor Terron. All rights reserved.
# Institute of Astrophysics of Andalusia, IAA-CSIC
#
# This file is part of IAAraf.
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


""" This module implements a wrapper for the IRAF task 'qphot' (quick aperture
    photometer) and 'txdump' (print fields from selected records in an
    APPHOT/DAOPHOT text database). The routines implemented in this module
    provide a way to automatically do photometry on an image and save to a file
    the specified fields.

    Functions:

    extract_star_coordinates() -- dump to a file the coordinates of the stars.
    qphot() -- run IRAF's qphot on an image, and extract selected fields.

"""

import pyraf.iraf
from pyraf.iraf import digiphot, apphot     # digiphot is not used, but has to
import tempfile                             # be imported before appphot
import subprocess
import time
import sys
import os

# IAAraf modules
import methods
import fitsimage
import coordinates
import seeing
import style


def extract_star_coordinates(path_to_image):
    """ Extract the coordinates of the stars contained in a SExtractor catalog.

    TODO: update this
 
    The method receives the path to an image file, runs SExtractor on it and
    returns a list that contains an instance of Point for each star that was
    detected. That is, the coordinates of the stars present in the SExtracted
    catalog are returned as Point instances.

    Keyword arguments:
    path_to_image - path to the image the coordinates of whose stars are to be
                    extracted.
    
    """


# TODO: this should be properly commented
""" Returns the string repr of x, y, replacing qphot's INDEF with 0.0000 """
# TODO; borrar???
def parse_qphot_indef(x, y):

    if x == "INDEF":
        parsed_x = str(0.0000)
    else:
        parsed_x = str(x)

    if y == "INDEF":
        parsed_y = str(0.0000)
    else:
        parsed_y = str(y)

    return parsed_x, parsed_y





def qphot(path_to_image,
          list_of_pixels = None,
          fields = None,
          cbox = 6,
          annulus = 13,
          dannulus = 8,
          aperture = 11,
          output_path = None):

    """ Run IRAF's qphot on an image, extracting from its output the specified 
        fields.

    This method is a wrapper, equivalent to (1) running 'qphot' on an image and
    (2) using 'txdump' in order to extract some fields from the text database
    that the former outputs. The goal of this routine is to make it possible
    to easily do photometry on an image, saving to a file the specified fields
    for all the objects (that is, when running txdump the parameter
    'parameters' is set to yes).

    In the first step, photometry will be calculated for all the objects in the
    image. This means that, when calling 'qphot', it is not necessary to pass a
    file to which the x and y coordinates of the objects to be measured have
    been saved. Instead, SExtractor will be run on the image and the catalog
    parsed in order to extract the coordinates of all the objects that were
    detected and save them to a temporary file (a).

    The output of 'qphot' is saved to another temporary file (b), on which
    'txdump' is then run to extract the specified fields from the APPHOT text
    database that was produced by 'qphot'. These fields are saved to a third
    temporary file (c), whose path is returned by the method upon successful
    execution.

    The temporary files (a) and (b) are automatically deleted by the method. 
    Note, however, that the user is responsible for deleting the last one, (c),
    when done with it.
    
    Keyword arguments:
    path_to_image - path to the image containing the objects to be measured.
    fields - list containing the fields to be extracted from each object on
             which photometry is done. If the photometry were to be calculated
             manually, these would be the fields that would be passed to
             txdump when extracting data from the file output by qphot.
    cbox - the width of the centering box, in pixels.
    annulus - the inner radius of the sky annulus, in pixels.
    dannulus - the width of the sky annulus, in pixels.
    aperture - the single aperture radius, in pixels.
    output_path - path to the file to which the calculated photometry will be
                  saved. If not specified, the photometry will be saved to a
                  temporary file.

    TODO: update with explanation of how now the method receives the coordinates
          on which to do photometry. If none, it its done on the stars of the file
        Decir tb que los INDEF pasan a ser 0.0000

    """

    # Mutable objects, such as lists, are dangerous to use as default values.
    # The reason is that default parameter values are always evaluated when,
    # and only when, the 'def' statement they belong to is executed.
    # The workaround is to use a placeholder value instead of modifying the
    # default value. None is a commong value. For details, see:
    # http://effbot.org/zone/default-values.htm
    if fields is None:
        fields = ["xcenter", "ycenter", "mag", "msky", "stdev"]
    
    # At least one field for txdump must have been specified
    if not fields:
        raise  Exception, "The list of specified fields cannot be left empty."

   
    # Temporary file to which the coordinates received by qphot will be saved
    coords_fd, temporary_coords = tempfile.mkstemp(suffix='.qphot_coords', text=True)
    # Temporaty file to which the text database produced by qphot will be saved
    output_fd , temporary_output = tempfile.mkstemp(suffix='.qphot_output', text=True)
    os.close(output_fd)

    # File whose path is to be returned by the method. The fields from the 
    # selected records in the output of qphot will be saved here.
    if output_path is None:
        photom_fd, output_path = tempfile.mkstemp(suffix='.qphot_photom', text=True)
        os.close(photom_fd)  # we need an object-like descriptor, not the OS-level
                             # handle to the open file that mkstemp() returns.

    photom_fd = open(output_path, 'wt')
    photom_fd.write("# Photometry for %s\n" % path_to_image)
    photom_fd.write("# IAAraf | Institute of Astrophysics of Andalusia, IAA-CSIC\n")
    photom_fd.write("# File created on %s UTC\n" % time.asctime(time.gmtime(time.time())))
    photom_fd.write("# Meaning of the columns: %s\n" % ",".join(fields))
    photom_fd.write("#\n")

    # Even if the file is empty, temporary_output must be deleted before
    # calling qphot. Otherwise, an error message, stating that the operation 
    # "would overwrite existing file", will be thrown.
    os.remove(temporary_output)

    # qphot must receive a text file with the initial coordinates for the
    # objects to be measured. We need, then, to create a temporary file, in
    # which the coordinates of the stars, extracted from the SExtractor
    # catalog, will be saved, one per line.

    # No coordinates were specified, so extract the stars of the image
    if list_of_pixels is None:
        image_instance = seeing.FITSeeingImage(path_to_image)
        image_instance.catalog().load_stars()
        for image_star in image_instance.catalog():
            os.write(coords_fd, str(image_star.x()) + ' ' + str(image_star.y()) + '\n')
    else:
        for pixel in list_of_pixels:
            os.write(coords_fd, str(pixel.x()) + ' ' + str(pixel.y()) + '\n')

    os.close(coords_fd)
   
    # Run qphot on the image, save to temporary_coords...
    apphot.qphot(path_to_image, cbox = cbox, annulus = annulus, \
                 dannulus = dannulus, aperture = aperture, \
                 coords = temporary_coords, output = temporary_output, \
                 interactive = 'no')
    
    # ... and extract the specified records from the qphot output. The task 
    # txdump does not seem to include an option to redirect the output to a
    # file, so we will redirect the entire standard output to it.
    sys.stdout = photom_fd 
    pyraf.iraf.txdump(temporary_output, fields = ",".join(fields), expr = 'yes')
    sys.stdout = sys.__stdout__  # restore stdout back to normal
    photom_fd.close()

    # Two of the three temporary files are no longer needed
    os.remove(temporary_coords)
    os.remove(temporary_output)    

    # Replace qphot's "IRAF" values with 0.000, our temporary solution in order
    # to deal with undetermined values when doing differential photometry.
    retcode = subprocess.call("sed -i 's/INDEF/0.0000/g' " + output_path, shell=True)
    if retcode != 0:
        raise RuntimeError("There was an error while replacing qphot's 'INDEF' values with zeroes.")
        
    return output_path


def fotdif_i_mej4_dir():
    """ Return the path to the time-series Fortran 95 code directory. """ 
    return os.path.normpath(os.getcwd() + "/diffphot/")


def fotdif_i_mej4_binaries():
    """ Return the path to the compiled time-series Fortran 95 code. """
    return os.path.normpath(fotdif_i_mej4_dir() + "/fotdif_i_mej4")


def compile_fotdif_i_mej4():
    """ Compile the time-series Fortran 95 code.
    
    This method compiles the Fortran 95 code contained in the diffphot/
    directory, generating an executable file which will be saved to
    fotdif_i_mej4_binaries().

    The return code of the 'make all' command is returned.

    """
    # Compile the Fortran source code and get the return code
    retcode = subprocess.call('cd ' + fotdif_i_mej4_dir() + '; make all', shell=True)

    # If the code was successfullt compiled, make sure that it is executable
    if os.path.exists(fotdif_i_mej4_binaries()):
        os.chmod(fotdif_i_mej4_binaries(), 0755) # u=rwx, og=rx

    return retcode


def build_fotdif_i_mej4():
    """ Make sure that the compiled time-series program is where expected.

    The method tests for the existence of the compiled Fortran 95 code that
    does time-series photometry. In case it can not be found, the code will be
    compiled and the resulting file made executable. If the compilation fails,
    an exception will be raised. A successful invocation of this method, thus,
    means that the time-series program is compiled and ready to go.

    """

    if not os.path.exists(fotdif_i_mej4_binaries()):
        print style.prefix() + "The time-series photometry code has not yet been compiled. Please wait."
        sys.stdout.flush()
        if compile_fotdif_i_mej4() != 0:
            raise RuntimeError("The compilation of the time-series photometry code failed.")
        else:
            print style.prefix() + "The compilation was successful. Hooray! You are ready to go."

if __name__ == "__main__":

    print qphot(sys.argv[1])


    sys.exit()
    build_fotdif_i_mej4()   # make sure that the Fortan 95 is compiled and in place TODO: que tambien detecte el otro, TODO encapsular en un mismo metodo
    working_directory = "qphot_files" #tempfile.mkdtemp()
    print style.prefix() + "Working directory is '" + working_directory + "', cleaning it up."  # TODO: print "cleaning it up" only if NOT EMPTY
    methods.empty_dir(working_directory)

    ref_image = fitsimage.find_fits(sys.argv[1])
    images_set = fitsimage.find_fits(sys.argv[2:], False)

    print ">> Reference image: %s" % ref_image.path()
    print ">> %d FITS images found" % len(images_set)
    for img in images_set.sort():
        qphot_output_path = os.path.normpath(working_directory + '/' + str(img.julian_date()))
        print ">> Saving photometry of %s to %s..." % (img.path(), qphot_output_path),
        sys.stdout.flush()
        qphot(img.path(), output_path = qphot_output_path)
        print "done"
    print "Done! | Hold me closer, tiny dancer... ^_^"
        
