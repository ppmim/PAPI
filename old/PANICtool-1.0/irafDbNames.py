################################################################################
#
# FIEStool
#
# irafDbNames.py
#
# Last update 29/09/2005
#
################################################################################

"""
   Terrible routine to get the proper path names of IRAF database files.
   ...Sigh...
"""

import os
import pyfits

def apname(filename, extension=0):

  filebase, file    = os.path.split(filename)
  file,     ext     = os.path.splitext(file)

  if extension != 0:

    try:
      data = pyfits.open(filename)
      extname = data[extension].header['EXTNAME']
      extname = "_%s_" % extname
      data.close()
    except: # Make sure to pass if the frame is not defined
      extname = ""

  else :

    extname = ""


  # Create the full path names for the source and target files
  return os.path.join(filebase, 'database', 'ap%s%s' % (file, extname) )


def ecname(filename, extension=0):

  filebase, file    = os.path.split(filename)
  file,     ext     = os.path.splitext(file)

  if extension != 0:

    try:
      data = pyfits.open(filename)
      extname = data[extension].header['EXTNAME']
      extname = "_%s_" % extname
      data.close()
    except: # Make sure to pass if the frame is not defined
      extname = ""

  else :

    extname = ""


  # Create the full path names for the source and target files
  return os.path.join(filebase, 'database', 'ec%s%s' % (file, extname) )

