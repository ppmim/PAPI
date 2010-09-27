import sys

deferror = "FIEStool failed to start!"


def checkImport(package, minversion):

  try:
    thispackage = __import__(package)
  except ImportError:
    print deferror
    print "Required external package '%s' not installed" % package
    raise

  if thispackage.__version__ < minversion:
    print deferror
    print "Version of external package '%s' lower than minimum required version of %s" % (package, minversion)
    print "Your version of '%s' is %s" % (package, thispackage.__version__)
    raise ImportError

  pass





