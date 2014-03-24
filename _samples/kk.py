#!/usr/bin/env python

################################################################################
#
# symple test programm for PyRAF
#
# symple.py
#
# Last update 14/01/2008
#
################################################################################

#"""
#   Test routines for data reduction with PyRAF.
#"""

################################################################################

# Import necessary modules


import os
import sys
import time
import shutil
import pp ## parallel python

# Interact with FITS files
import pyfits


#Own modules
import mkBadPix


from Tkinter import *
from pyraf import iraf
from iraf import noao
from iraf import imred
from iraf import ccdred

from time import sleep

from pyraf import iraf
from iraf import *
from time import *

import threading
from threading import Thread





class myClass:

    a=1
    b=2
    def __init__(self):
        pass
        #self.a=1
        #self.b=2
        
    def suma(self):
        self.a=self.a+1
        self.b+=2

        
def testf():
    
    print 'Hola'

    
################################################################################      
 
if __name__=="__main__":
    
    clase1=myClass()
    clase1.suma()

    clase2=myClass()
    clase2.suma()

    print "clase1 a=" ,clase1.a
    print "clase2 a=" ,clase2.a
    
    
    print 'Fin'
        
################################################################################

