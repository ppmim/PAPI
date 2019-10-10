#! /usr/bin/env python 
#
# This program extracts the FWHM measurements from the .mag files using IRAF
#

from os import chdir
import os.path
from pyraf import iraf
from pylab import *
import pyraf.iraf as iraf

import param as par

def fwhm(dummy):
  workdir = sort.get_workdir()
  filt  = par.filter[0]

  iraf.noao(_doprint=0)
  iraf.obsutil(_doprint=0)
  iraf.unlearn("psfmeasure")

  if not os.path.isfile('data/ds92.reg'):
    ds9 = open('data/ds92.reg', 'w+')
    ds9.write(par.pixel_location[0], par.pixel_location[1])
    ds9.close()

  filedir = workdir+"/object/" + filt
  chdir(filedir)

  iraf.noao.obsutil.psfmeasure.coords = "mark1"
  iraf.noao.obsutil.psfmeasure.display = "no"
  iraf.noao.obsutil.psfmeasure.size = "FWHM"
  iraf.noao.obsutil.psfmeasure.radius = "13"
  iraf.noao.obsutil.psfmeasure.iterati = "5"
  iraf.noao.obsutil.psfmeasure.imagecu = str(workdir)+"/data/ds92.reg"
  iraf.noao.obsutil.psfmeasure.logfile = "fwhm.log"

  if not os.path.isfile('fwhm.log'):
    iraf.noao.obsutil.psfmeasure("rfb*.fits")

  F = loadtxt('fwhm.log',skiprows = 3, usecols = (4,5), dtype = str)
  fwhma = F[:-1,0]

  fwhm = []
  for i in range(len(fwhma)):
    fwhm.append(float(fwhma[i]))

  fwhm=array(fwhm)
  fwhm_mean=fwhm.mean()

  print("\nMean FWHM value: "+str(fwhm_mean) + "\n")

  chdir(workdir)

  np.save('fwhm.npy', fwhm)
  np.save('fwhm_mean.npy', fwhm_mean)

  if not os.path.isdir('data'):
    os.mkdir("data")

  os.system('mv fwhm.npy '+workdir+'/data')
  os.system('mv fwhm_mean.npy '+workdir+'/data')
