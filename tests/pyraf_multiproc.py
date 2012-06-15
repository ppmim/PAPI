#! /usr/env/python

# Victor Terron |vterron@iaa.es
# This is an example of what seems to be a bug in PyRAF.

# It seems that PyRAF is messing with the multiprocessing module. The reason
# why I think so is because two consecutive calls to difiphot.apphot.qphot work
# fine, as well as using it as the function passed to map_async. But calling
# qphot and then map_async, with the instantiation of the pool of workers done
# in-between, makes the scripts fail with the most arcane of errors, such as
# (yes, this is a real example):
#
# "Exception Exception OSErrorOSError: : ((1010, , ''NNoo cchhiilldd
# pprroocceesssseess'')) in in <bound method Subprocess.__del__ of <Subprocess
# '/iraf/iraf/noao/bin.linux/x_apphot.e -c', at 3c8c2d8>><bound method
# Subprocess.__del__ of <Subprocess '/iraf/iraf/noao/bin.linux/x_apphot.e -c',
# at 3c8c2d8>> ignored"
#
# However, everything goes as expected is the pool of workers is created before
# qphot is ever used, directly or indirectly, in the script. It seems the first
# execution of PyRAF's qphot is affecting how the multiprocessing module
# works. That's why we need to create the pool before PyRAF runs, so that we
# use the original, 'unmodified' code.

from __future__ import with_statement

import pyraf.iraf
from pyraf.iraf import digiphot, apphot  # 'digiphot.apphot' package
import multiprocessing
import os
import random
import sys
import tempfile

def get_coordinates(size = 250):
    """ Return the path to a temporary file with the coords of random pixels """
    fd, path = tempfile.mkstemp(text = True)
    for _ in xrange(size):
        random_x = random.uniform(1, 2000)
        random_y = random.uniform(1, 2000)
        os.write(fd, "%f\t%f\n" % (random_x, random_y))
    os.close(fd)
    return path

def photometry(args):
    path, coordinates = args
    with open(os.devnull, 'w') as stdout:
        apphot.qphot(path, interactive = 'no', cbox = 0,
                     aperture = 11, annulus = 13, dannulus = 8,
                     coords = coordinates, Stdout = stdout)

if __name__ == "__main__":

    try:
        path = sys.argv[1]
    except IndexError:
        print "usage: %s PATH_TO_FITS" % sys.argv[0]
        sys.exit(1)

    try:
        coordinates = get_coordinates()

        # [Uncomment this line to make PyRAF fail] The multiprocessing's pool
        # of workers must be created *before* qphot is called the first time.
        photometry((path, coordinates))

        pool = multiprocessing.Pool(multiprocessing.cpu_count())
        queue = multiprocessing.Queue()

        photometry((path, coordinates))
        map_async_args = [(path, coordinates)] * 10
        result = pool.map_async(photometry, map_async_args)
        result.wait()
        result.get()  # re-raise exceptions, if any

    finally:
        try: os.unlink(coordinates)
        except OSError: pass

