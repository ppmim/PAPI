#!/usr/bin/env python

import multiprocessing
import os
import time

def f(i):
    print "I am process number",os.getpid(),": i =",i
    time.sleep(1)
    return i*i

pool = multiprocessing.Pool(maxtasksperchild=1)

print pool.map(f, range(10))