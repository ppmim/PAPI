#!/usr/bin/env python

__doc__ = """
Check the existence of necessary Python modules for PAPI
"""

import sys
import string

def testmodule(modulename, moduleversion):
    """
    test if a Python module is installed, and
    if yes, if its version is equal or higher than
    a reference version
    """

    bold = "\033[1m"
    probbold = "\033[1;34m"
    reset = "\033[0;0m"
    
    print bold + \
          "Testing Python module installation for module '%s':" % \
          (modulename) + reset
    print "PAPI needs at least version %s" % (moduleversion)
    
    try:
        mod = __import__(modulename)
        refversion = string.split(moduleversion, ".")
        currversion = string.split(mod.__version__, ".")

        if map(int, currversion) < map(int, refversion):            
            print probbold + "PROBLEM: You have it with V%s\n" % \
                  (mod.__version__) + reset
        else:
            print "Your version %s of '%s' is fine!\n" % \
                  (mod.__version__, modulename)
    except:
        print probbold
        print probbold + \
              "PROBLEM: You do not have the Python module '%s' installed!\n" % \
              (modulename) + reset



# define the Python modules, and the versions we need:
PAPImodules = { 'pyraf' : '1.0', 'pyfits' : '1.1', 'numpy' : '1.1', 
               'matplotlib' : '0.98.1', 'scipy': '0.10', 'qt': '0.0'  }

bold = "\033[1m"         # print bold
probbold = "\033[1;34m"  # print bold blue
reset = "\033[0;0m"      # reset special print settings

print bold + "PAPI Python checking tool" + reset
print
print bold + "Checking Python Version:" + reset
print "PAPI needs Python Version 2.Y with Y>=2.7"
pyversion = string.split(string.replace(string.split(sys.version)[0], 
                                        '+', ''), ".")
# well, Python version 3 just gives us a syntax error at the
# first print statement :-)
if map(int, pyversion) >= [3, 0, 0] or map(int, pyversion) < [2, 7, 0]:
    print probbold + "PROBLEM: You have Python V%s.%s.%s\n" \
                      % (pyversion[0], pyversion[1], pyversion[2]) + reset
    print
else:
    print "Your Python version %s is fine!" % (string.split(sys.version)[0])
    print

for modulename in PAPImodules.keys():
    testmodule(modulename, PAPImodules[modulename])
    
