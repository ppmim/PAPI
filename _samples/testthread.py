import sys
from time import sleep

import threading
from threading import Thread

import pyraf
from pyraf import iraf
from iraf  import noao
from iraf  import imred
from iraf  import ccdred
from iraf  import echelle


class myThread(Thread):

  def __init__(self):
    threading.Thread.__init__(self)

  def run(self):

    # Please substitute below the values for indir, infile, 
    # outfile, reffile and parfile
    #iraf.chdir('indir')
    echelle.apsum(input='infile', output='outfile',
                  references='reffile', ParList='parfile')

    print 'CHILD: apsum finished OK'

if __name__=="__main__":

  aThread = myThread()
  aThread.start()

  while aThread.isAlive():
    print '...parent's heartbeat...'
    sleep(5)
    
  print 'PARENT: Thread finished, but did apsum?'
