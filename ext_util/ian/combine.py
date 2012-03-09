#! /usr/bin/env python
# Copyright (c) 2011-2012 IAA-CSIC  - All rights reserved. 
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


################################################################################
#
#
# PAPI (PANIC PIpeline)
#
# dxtalk.py
#
# Created    : 24/02/2012    jmiguel@iaa.es -
# Last update: 
# TODO
################################################################################

# IJMC code
import analysis as an
#import ir
import phot
import numpy as np

#sigma
import sigma_clip as sg


import sys
import fileinput
import pyfits



def red_dark(dark_file, verbose=False, sigma=3, niter=1, clobber=False):  
    """Combine dark frames into a single dark frame:""" 
     
    #_sdark = obs._proc + ('superdark%03g' % num) + obs._suffix
      
    #darklist = obs.rawall[ind]  
    #framecenter = obs.naxis1/2, obs.naxis2/2  
    #framesize = obs.naxis1, obs.naxis2  

    _sdark = "/tmp/prueba_combine.fits"
        
    darklist = [line.replace( "\n", "") for line in fileinput.input(dark_file)]
    framecenter = 1024, 1024  
    framesize = 2048, 2048  

    
    if verbose: print "dark file list>>" , darklist  
    darkstack = phot.subreg2(darklist, framecenter, framesize, verbose=verbose)  
    
    if verbose: print "loaded dark stack" 

    """
    for i in range(0,len(darklist)):
        print "MEAN-r of file %s = %s  %s"%(darklist[i],an.meanr(darkstack[i], 3, 1), np.median(darkstack[i]))
        
 
    return 0
    """
    
    superdark = np.zeros(framesize, np.float32)
    #superdark = np.median(darkstack,0) 
    
    for ii in range(2048):  
        for jj in range(2048):
            ###m_line = darkstack[:,ii,jj]
            ###A = sg.sigmaclip(m_line)[0]
            ###superdark[ii,jj] = np.mean(m_line[A])
            superdark[ii,jj] = np.sqrt(darkstack[:,ii,jj]).sum()#np.median(darkstack[:,ii,jj]) 
            #superdark[ii,jj] = an.meanr(darkstack[:,ii,jj],nsigma=sigma, 
            #                            niter=niter)  
    
    hdr0 = pyfits.getheader(darklist[0])  
    hdr0.update('sup-dark', 'super-dark created by IR.py')  
    pyfits.writeto(_sdark, superdark, header=hdr0, clobber=clobber, output_verify='ignore')  
    
    if verbose: print "Done making superdark!"  
    return



################################################################################
# main
################################################################################
if __name__ == "__main__":
    
    print "Testing my_combine !!"
    
    red_dark( sys.argv[1] )
    
    print "what's up !"
    
    
