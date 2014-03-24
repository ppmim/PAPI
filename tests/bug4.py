#from multiprocessing import Pool
import time
import multiprocessing
import sys

from pyraf import iraf
from iraf import noao
from iraf import imred
from iraf import ccdred
from iraf import mscred

iraf.prcacheOff()

def _pickle_method(method):
    """
    Pickle methods properly, including class methods.
    """
    func_name = method.im_func.__name__
    obj = method.im_self
    cls = method.im_class
    return _unpickle_method, (func_name, obj, cls)

def _unpickle_method(func_name, obj, cls):
    """
    Unpickle methods properly, including class methods.
    """
    for cls in cls.mro():
        try:
            func = cls.__dict__[func_name]
        except KeyError:
            pass
        else:
            break
    return func.__get__(obj, cls)        

import copy_reg 
import types 

copy_reg.pickle(types.MethodType,  
    _pickle_method,  
    _unpickle_method) 


def run_iraf():
    print "Start run_iraf"
    iraf.unlearn("darkcombine")
    #import pdb; pdb.set_trace()
    iraf.mscred.darkcombine(input = "@/tmp/darks.txt",
                        output = "/tmp/dark.fits",
                        combine = "average",
                        ccdtype = '',
                        process = 'no',
                        reject = "minmax",
                        nlow = "0",
                        nhigh = "1",
                        nkeep = "1",
                        scale = "mode")
                        #expname = 'EXPTIME')
                        #ParList = _getparlistname('darkcombine')
                        #)
    print "End run_iraf"

def wrapper():
    c = fooclass()
    c.run_iraf()
    
class fooclass(object):
    def f(self):
        time.sleep(1)
        return "Nothing"

    def run_iraf(self):
        #import pdb; pdb.set_trace()
        #return 0
        iraf.darkcombine(input = "@/tmp/darks.txt",
                        output = "/tmp/dark.fits",
                        combine = 'average',
                        ccdtype = '',
                        process = 'no',
                        reject = 'minmax',
                        nlow = '0',
                        nhigh = '1',
                        nkeep = '1',
                        scale = 'mode',
                        #expname = 'EXPTIME'
                        #ParList = _getparlistname('darkcombine')
                        )
        
def f(x):
    time.sleep(1)
    return x*x

if __name__ == '__main__':
    #run_iraf()
    #sys.exit(0)
    pool = multiprocessing.Pool(processes=4)              # start 4 worker processes
    c = fooclass()
    result = pool.apply_async(c.run_iraf, args=())    # evaluate "f(10)" asynchronously
    #result = multiprocessing.Pool(processes=4).apply_async(f, [10])
    print result.get(timeout=5)           # prints "100" unless your computer is *very* slow
    print pool.map(f, range(10))          # prints "[0, 1, 4,..., 81]"
