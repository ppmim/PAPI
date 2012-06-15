#!/usr/bin/env python

# Copyright (c) 2009-2012 IAA-CSIC  - All rights reserved. 
# Author: Jose M. Ibanez. 
# Instituto de Astrofisica de Andalucia, IAA-CSIC
#
# This file is part of PAPI (PANIC Pipeline)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
   Main routines for the loggin function.
"""
import logging
import datetime


class ColorFormatter(logging.Formatter):
 
    def color(self, level=None):
        codes = {\
            None:       (0,   0),
	    'DEBUG':    (1,  32), # green
            'INFO':     (0,  33), # yellow
            'WARNING':  (1,  34), # blue
            'ERROR':    (1,  31), # red
            'CRITICAL': (1, 101), # negro, fondo rojo
            }
        return (chr(27)+'[%d;%dm') % codes[level]
 
    def format(self, record):
        retval = logging.Formatter.format(self, record)
        return self.color(record.levelname) + retval + self.color()


### We define two logging handlers (Console and File), each one can have 
### different properties (level, formater, ...)
## Console
console = logging.StreamHandler()
console.setLevel(logging.DEBUG) # here we set the level for console handler
## NOTE: Handler.setLevel() method, just as in logger objects, specifies the 
##lowest severity that will be dispatched to the appropriate destination.
## Why are there two setLevel() methods? The level set in the logger determines 
##which severity of messages it will pass to its handlers. The level set in 
## each handler determines which messages that handler will send on.
console.setFormatter(ColorFormatter('    [%(name)s]: %(asctime)s %(levelname)-8s %(module)s:%(lineno)d: %(message)s'))
logging.getLogger('PAPI').addHandler(console)

## File
datetime_str = str(datetime.datetime.utcnow()).replace(" ","T")

file_hd = logging.FileHandler("/tmp/papi_" + datetime_str + ".log")
file_hd.setLevel(logging.DEBUG) # here we set the level for File handler 
formatter = logging.Formatter('[%(name)s]: %(asctime)s %(levelname)-8s %(module)s:%(lineno)d: %(message)s')
file_hd.setFormatter(formatter)
logging.getLogger('PAPI').addHandler(file_hd)

## define the global log level
logging.getLogger('PAPI').setLevel(logging.DEBUG) # debug is the lowest level
## define the global variable used whole around the PAPI sources
log = logging.getLogger('PAPI')












"""
# STUFF NOT USED 
logging.basicConfig(level=logging.DEBUG,
    format='%(asctime)s %(levelname)-8s %(module)s:%(lineno)d: %(message)s',
    datefmt='%a, %d %b %Y %H:%M:%S',
    filename='/tmp/panicTool.log',
    filemode='w')
#console = logging.StreamHandler()
#console.setLevel(logging.DEBUG)
#formatter = logging.Formatter('%(filename)s %(name)-12s: %(levelname)-8s %(message)s')
# tell the handler to use this format
#console.setFormatter(formatter)
#logging.getLogger('').addHandler(console)
log = logging.getLogger('panicQL.log')
log.setLevel(logging.INFO)


def initLog(log_filename, log_level):
    "" actually not used ...""
    logging.basicConfig(level=log_level,
                        format='%(asctime)s %(levelname)-8s %(module)s:%(lineno)d: %(message)s',
                        datefmt='%a, %d %b %Y %H:%M:%S',
                        filename=log_filename,
                        filemode='w')
    
"""



