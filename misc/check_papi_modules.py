#!/usr/bin/env python

__doc__ = """
Check the existence of necessary Python modules for PAPI
"""

import sys
import string


bold = "\033[1m"         # print bold
probbold = "\033[1;31m"  # print bold red 
reset = "\033[0;0m"      # reset special print settings


def testmodule(modulename, moduleversion):
    """
    test if a Python module is installed, and
    if yes, if its version is equal or higher than
    a reference version
    """

    print bold + \
          "Testing Python module installation for module '%s':" % \
          (modulename) + reset
    print "PAPI needs at least version %s" % (moduleversion)
    
    try:
        mod = __import__(modulename)
        refversion = string.split(moduleversion, ".")
        if modulename=="PyQt4.QtCore":
            cv = mod.QtCore.PYQT_VERSION_STR
            currversion = string.split(cv , ".")
            currversion = [ a.split('-')[0] for a in currversion ]
        else:
            cv = mod.__version__
            currversion = string.split(cv , ".")
            currversion = [ a.split('-')[0] for a in currversion ]

        if map(int, currversion) < map(int, refversion):            
            print probbold + "PROBLEM: You have it with V%s\n" % \
                  (cv) + reset
        else:
            print "Your version %s of '%s' is fine!\n" % \
                  (cv, modulename)
    except:
        print probbold
        print probbold + \
              "PROBLEM: You do not have the Python module '%s' installed!\n" % \
              (modulename) + reset


def check_modules():
    # --------------------
    # Check Python version
    # --------------------
    print bold + "PAPI Python checking tool" + reset
    print bold + "=========================" + reset
    print
    print bold + "Checking Python Version:" + reset
    print "PAPI needs Python Version 2.Y with Y>=2.7"
    pyversion = string.split(string.replace(string.split(sys.version)[0], 
                                            '+', ''), ".")
    # well, Python version 3 just gives us a syntax error at the
    # first print statement :-)
    if map(int, pyversion) >= [3, 0] or map(int, pyversion) < [2, 7]:
        print probbold + "PROBLEM: You have Python V%s.%s\n" \
                          % (pyversion[0], pyversion[1]) + reset
        print
    else:
        print "Your Python version %s is fine!" % (string.split(sys.version)[0])
        print
    
    # ----------------------------------------------------
    # Define the Python modules, and the versions we need
    # ----------------------------------------------------
    PAPImodules = { 'numpy' : '1.6', 'pyraf' : '1.1', 'pyfits' : '3.0',  
                   'matplotlib' : '0.98.1', 'scipy': '0.10', 'PyQt4.QtCore': '4.8',
                   'pywcs': '1.10', 'astropy': '0.3' }
    
    # -----------------
    # Check the modules
    # -----------------
    for modulename in PAPImodules.keys():
        testmodule(modulename, PAPImodules[modulename])
    
################################################################################
# main
if __name__ == "__main__":
    check_modules()
    sys.exit(0)
    
