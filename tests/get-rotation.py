#! /usr/bin/env python

""" Run Astrometry.net, get field rotation angle. """

# Author: Victor Terron (c) 2014
# Email: `echo vt2rron1iaa32s | tr 132 @.e`
# License: GNU GPLv3

import os
import os.path
import pyfits
import re
import shutil
import sys
import tempfile
import subprocess

ASTROMETRY_COMMAND = 'solve-field'

class AstrometryNetNotInstalled(StandardError):
    """ Raised if Astrometry.net is not installed on the system """
    pass

class AstrometryNetError(subprocess.CalledProcessError):
    """ Raised if the execution of Astrometry.net fails """
    pass

class AstrometryNetUnsolvedField(subprocess.CalledProcessError):
    """ Raised if Astrometry.net could not solve the field """

    def __init__(self, path):
        self.path = path

    def __str__(self):
        return "%s: could not solve field" % self.path

def astrometry_net(path, ra = None, dec = None, radius = 1):
    """ Do astrometry on a FITS image using Astrometry.net.

    Use a local build of the amazing Astrometry.net software [1] in order to
    compute the astrometric solution of a FITS image. This software has many,
    many advantages over the well-respected SCAMP, but the most important one
    is that it is a blind astrometric calibration service. We do not need to
    know literally anything about the image, including approximate coordinates,
    scale and equinox. It just works, giving us a new FITS file containing the
    WCS header. The second element of the returned tuple is a list with all the
    lines that were written by Astrometry.net to both standard output and error.

    In order for this function to work, you must have built and installed the
    Astrometry.net code in your machine [2]. The main high-level command-line
    user interface, 'solve-field', is expected to be available in your PATH;
    otherwise, the AstrometryNetNotInstalled exception is raised. Note that you
    also need to download the appropriate index files, which are considerably
    heavy. At the time of this writing, the entire set of indexes built from
    the 2MASS catalog [4] has a total size of ~32 gigabytes.

    Raises AstrometryNetError if Astrometry.net exits with a non-zero status
    code and AstrometryNetUnsolvedField if the CPU time limit, set in the
    backend.cfg file (by default located in /usr/local/astrometry/etc/) is hit.

    [1] http://astrometry.net/
    [2] http://astrometry.net/doc/build.html
    [3] http://astrometry.net/doc/readme.html#getting-index-files
    [4] http://data.astrometry.net/4200/
    [5] https://groups.google.com/d/msg/astrometry/ORVkOk0jSZg/PeCMeAJodyAJ

    Keyword arguments:

    ra,
    dec,
    radius - restrict the Astrometry.net search to those indexes within
             'radius' degrees of the field center given by ('ra', 'dec').
             Both the right ascension and declination must be given in order
             for this feature to work. The three arguments must be expressed
             in degrees.

    """

    emsg = "'%s' not found in the current environment"
    if not subprocess.check_output(['which', ASTROMETRY_COMMAND]):
        raise AstrometryNetNotInstalled(emsg % ASTROMETRY_COMMAND)

    root, ext = os.path.splitext(os.path.basename(path))
    tempfile_prefix = '%s_' % root
    # Place all output files in this directory
    kwargs = dict(prefix = tempfile_prefix, suffix = '_astrometry.net')
    output_dir = tempfile.mkdtemp(**kwargs)

    # Path to the temporary FITS file containing the WCS header
    kwargs = dict(prefix = '%s_astrometry_' % root, suffix = ext)
    with tempfile.NamedTemporaryFile(**kwargs) as fd:
        output_path = fd.name

    # If the field solved, Astrometry.net creates a <base>.solved output file
    # that contains (binary) 1. That is: if this file does not exist, we know
    # that an astrometric solution could not be found.
    solved_file = os.path.join(output_dir, root + '.solved')

    # --dir: place all output files in the specified directory.
    # --no-plots: don't create any plots of the results.
    # --new-fits: the new FITS file containing the WCS header.
    # --no-fits2fits: don't sanitize FITS files; assume they're already valid.
    # --overwrite: overwrite output files if they already exist.

    args = [ASTROMETRY_COMMAND, path,
            '--dir', output_dir,
            '--no-plots',
            '--new-fits', output_path,
            '--no-fits2fits',
            '--overwrite']

    # -3 / --ra <degrees or hh:mm:ss>: only search in indexes within 'radius'
    # of the field center given by 'ra' and 'dec'
    # -4 / --dec <degrees or [+-]dd:mm:ss>: only search in indexes within
    # 'radius' of the field center given by 'ra' and 'dec'
    # -5 / --radius <degrees>: only search in indexes within 'radius' of the
    # field center given by ('ra', 'dec')

    if ra is not None:
        args += ['--ra', '%f' % ra]

    if dec is not None:
        args += ['--dec', '%f' % dec]

    if radius is not None:
        args += ['--radius', '%f' % radius]

    try:

        with tempfile.NamedTemporaryFile() as cmd_output:
            p = subprocess.Popen(args, stderr = subprocess.STDOUT, stdout = subprocess.PIPE)
            while True:
                line = p.stdout.readline()
                if not line and p.poll() is not None:
                    break
                if line:
                    print line.strip()
                    cmd_output.write(line)

            cmd_output.seek(0)
            output_lines = cmd_output.readlines()

        sys.stdout.flush()
        sys.stderr.flush()
        
        # .solved file must exist and contain a binary one
        with open(solved_file, 'rb') as fd:
            if ord(fd.read()) != 1:
                raise AstrometryNetUnsolvedField(path)

        return output_path, output_lines

    except subprocess.CalledProcessError, e:
        raise AstrometryNetError(e.returncode, e.cmd)
    # If .solved file doesn't exist or contain one
    except (IOError, AstrometryNetUnsolvedField):
        raise AstrometryNetUnsolvedField(path)
    finally:
        shutil.rmtree(output_dir, ignore_errors=True)

def get_PANIC_rotation(path):

    with pyfits.open(path, readonly = True) as hdu:
        header = hdu[0].header
        ra  = header['RA']
        dec = header['DEC']

    output_path, output_lines = astrometry_net(path, ra=ra, dec=dec)
    os.unlink(output_path)

    # Extract rotation angle from a line like the following one:
    # Field rotation angle: up is 0.214027 degrees E of N
    float_regexp = r"[-+]?\d*\.\d+|\d+"
    regexp = "Field rotation angle: up is ({0}) degrees".format(float_regexp)

    for line in output_lines:
        match = re.search(regexp, line)
        if match:
            return float(match.group(1))

if __name__ == "__main__":

    path = sys.argv[1]
    rotation = get_PANIC_rotation(path)
    print "Rotation (deg): {0}".format(rotation)
    sys.exit()

