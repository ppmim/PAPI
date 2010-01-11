import sys
from time import strftime

loglevel = 5

def put(message, level=0, errorflag=0):

  if errorflag == 1:
    sepstring='- ERROR:'
  elif errorflag > 1:
    sepstring='- warning:'
  else:
    sepstring='-'

  # Create a text string to be displayed in the box
  # strftime() gives a nicely formatted time stamp
  textstring = "%s (%i) %s %s\n" % (strftime('%H:%M:%S'), level, sepstring, message)

  if level <= loglevel:
    sys.stdout.write(textstring)


def put_error(message, loglevel=0, errorflag=1):
  put(message, loglevel, errorflag)


def put_warning(message, loglevel=0, errorflag=2):
  put(message, loglevel, errorflag)

