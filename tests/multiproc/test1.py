#!/usr/bin/python  
import multiprocessing  

#def _pickle_method(method):
#    """
#    Pickle methods properly, including class methods.
#    """
#    func_name = method.im_func.__name__
#    obj = method.im_self
#    cls = method.im_class
#    if isinstance(cls, type):
#        # handle classmethods differently
#        cls = obj
#        obj = None
#    if func_name.startswith('__') and not func_name.endswith('__'):
#        #deal with mangled names
#        cls_name = cls.__name__.lstrip('_')
#        func_name = '_%s%s' % (cls_name, func_name)
#
#    return _unpickle_method, (func_name, obj, cls)
#
#def _unpickle_method(func_name, obj, cls):
#    """
#    Unpickle methods properly, including class methods.
#    """
#    if obj is None:
#        return cls.__dict__[func_name].__get__(obj, cls)
#    for cls in cls.__mro__:
#        try:
#            func = cls.__dict__[func_name]
#        except KeyError:
#            pass
#        else:
#            break
#    return func.__get__(obj, cls)

#
# El codigo de arriba no funciona bien, no se por que ???
#  
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
  
class A(object):  
    def __init__(self):  
        print "A::__init__()"  
        self.weird = "weird"  
  
class B(object):  
    def doAsync(self, lala):  
        print "B::doAsync()"  
        return lala**lala  
  
    def callBack(self, result):  
        print "B::callBack()"  
        self.a.weird="wherio"  
        print result  
  
    def __init__(self, myA):  
        print "B::__init__()"  
        self.a = myA  
  
def callback(result):  
    print "callback result: " + str(result)  
  
def func(x):  
    print "func"  
    return x**x  
  
if __name__ == '__main__':  
        
    for i in range(1):
          
        pool = multiprocessing.Pool(2)  
        a = A()  
        b = B(a)  
        print a.weird  
        print "Starting iteration "  
        result1 = pool.apply_async(func, [4], callback=callback)  
        result2 = pool.apply_async(b.doAsync,
            [8])
    #,  
    #        callback=b.callBack)  
        print a.weird  
        print "result1: " + str(result1.get())  
        print "result2: " + str(result2.get())  
        print a.weird  
        print "End"
        
        del b
        del a
        del pool
        