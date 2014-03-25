#
# A test of `multiprocessing.Pool` class
#
# Copyright (c) 2006-2008, R Oudkerk
# All rights reserved.
#

import multiprocessing
import time
import random
import sys
import math

      
def _pickle_method(method):
    func_name = method.im_func.__name__
    obj = method.im_self
    cls = method.im_class
    return _unpickle_method, (func_name, obj, cls)

def _unpickle_method(func_name, obj, cls):
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

class MyClass(object):
    def __init__(self):
        self.a = 1.0
    def func1(self,a,b):
        return a*b;

#
# Functions used by test code
#

def calculate(func, args):
    result = func(*args)
    time.sleep(0.4)
    return '%s says that %s%s = %s' % (
        multiprocessing.current_process().name,
        func.__name__, args, result
        )

def calculatestar(args):
    return calculate(*args)

def mul(a, b):
    time.sleep(0.5*random.random())
    return a * b

def plus(a, b):
    time.sleep(0.5*random.random())
    return a + b

def f(x):
    return 1.0 / (x-5.0)

def pow3(x):
    return x**3

def hard(a,b):
    return 0
    r = 0.0
    for i in range(100000000):
        r+=math.sin((math.sin(a)**3*math.sqrt(b)*math.cos(b))**3)
    return r

def noop(x):
    pass

def eval_func_tuple(f_args):
    """Takes a tuple of a function and args, evaluates and returns result"""
    return f_args[0](*f_args[1:])  

def test():
    print 'cpu_count() = %d\n' % multiprocessing.cpu_count()

    #
    # Create pool
    #

    PROCESSES = 4
    print 'Creating pool with %d processes\n' % PROCESSES
    pool = multiprocessing.Pool(PROCESSES)
    print 'pool = %s' % pool
    print
    
    obj = MyClass()

    #
    # Tests
    #

    
    TASKS = [(mul, (i, 7)) for i in range(10)] + \
            [(plus, (i, 8)) for i in range(10)] +\
            [(hard, (i,100)) for i in range(10)]
    TASKS2 = [ (obj.func1, ( i, 5)) for i in range(10)]
    
    results = []
    #results = [pool.apply_async(calculate, t) for t in TASKS]
    results = [pool.map_async(calculatestar, [(mul,(3,5))] )]
    #results+= [pool.apply_async(calculate, t) for t in TASKS2]
    #results+= [pool.map_async(calculate, [(mul,(3,5))])]
    ##results = [pool.map_async(calculatestar, TASKS)]
    #imap_it = pool.imap(calculatestar, TASKS)
    #imap_unordered_it = pool.imap_unordered(calculatestar, TASKS)

    print 'Ordered results using pool.apply_async():'
    for r in results:
        print '\t', r.get()
    

    print "Joinning ..."
    pool.close()
    #pool.terminate()
    pool.join()


if __name__ == '__main__':
    multiprocessing.freeze_support()

    copy_reg.pickle(types.MethodType, _pickle_method, _unpickle_method)

    assert len(sys.argv) in (1, 2)
    print "argv1=", sys.argv[0]

    if len(sys.argv) == 1 or sys.argv[1] == 'processes':
        print ' Using processes '.center(79, '-')
    elif sys.argv[1] == 'threads':
        print ' Using threads '.center(79, '-')
        import multiprocessing.dummy as multiprocessing
    else:
        print 'Usage:\n\t%s [processes | threads]' % sys.argv[0]
        raise SystemExit(2)

    test()
