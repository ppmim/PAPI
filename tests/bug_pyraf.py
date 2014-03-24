#!/usr/bin/env python

import pyraf.iraf
from pyraf.iraf import noao, imred, ccdred, mscred
from pyraf.iraf import images, imutil

import time
import multiprocessing
import sys

pyraf.iraf.prcacheOff()

def run_iraf():
    print "Start run_iraf"

    pyraf.iraf.unlearn("imstat")
    #import pdb; pdb.set_trace()
    pyraf.iraf.imstat(images = "/tmp/dark.fits", fields="mean")

    print "End run_iraf"

    return "Success!"

if __name__ == '__main__':
    #run_iraf()
    #sys.exit(0)
    pool = multiprocessing.Pool(processes=2)      # start 4 worker processes
    result = pool.apply_async(run_iraf, args=())   
    #result = multiprocessing.Pool(processes=4).apply_async(f, [10])
    print result.get(timeout=5) 
    
    sys.exit(0)
