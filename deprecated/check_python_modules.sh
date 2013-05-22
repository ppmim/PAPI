#!/bin/sh

# HISTORY INFORMATION
# ===================

# 06.08.2010:
# I modified the check for the Python version so that also
# versions such as '2.6.5+' are treated correctly. The old test
# had problems with the '+' sign.  

#$1: absolute path to the Python Interpreter

"""":
if [ $# -ne 1 ]; then
  echo "$0 python_interpreter_with_absolute_path"  
  exit 1
fi

if [ -x $1 ]; then
  exec $1 "$0" "$@"
else
  tput bold
  tput setf 1
  echo "It seems that you do not have Python installed!" >&2
  echo "THELI needs Python 2.Y with Y >= 5" >& 2
  tput sgr0
  exit 1
fi
"""

__doc__ = """
Check the existence of necessary Python modules for THELI
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
    print "THELI needs at least version %s" % (moduleversion)
    
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
THELImodules = { 'pyfits' : '1.1', 'numpy' : '1.1', 'matplotlib' : '0.98.1' }

bold = "\033[1m"         # print bold
probbold = "\033[1;34m"  # print bold blue
reset = "\033[0;0m"      # reset special print settings

print bold + "THELI Python checking tool" + reset
print
print bold + "Checking Python Version:" + reset
print "THELI needs Python Version 2.Y with Y>=2.5"
pyversion = string.split(string.replace(string.split(sys.version)[0], 
                                        '+', ''), ".")
# well, Python version 3 just gives us a syntax error at the
# first print statement :-)
if map(int, pyversion) >= [3, 0, 0] or map(int, pyversion) < [2, 5, 0]:
    print probbold + "PROBLEM: You have Python V%s.%s.%s\n" \
                      % (pyversion[0], pyversion[1], pyversion[2]) + reset
    print
else:
    print "Your Python version %s is fine!" % (string.split(sys.version)[0])
    print

for modulename in THELImodules.keys():
    testmodule(modulename, THELImodules[modulename])
    
