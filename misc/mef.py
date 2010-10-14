#!/usr/bin/env python
################################################################################
#
#
# PAPI (PANIC PIpeline)
#
# mef.py
#
# Multi Extension FITS file basic operations
#
# Created    : 27/10/2010    jmiguel@iaa.es -
# Last update: 
# TODO
#       
################################################################################

################################################################################
# Import necessary modules

import getopt
import sys
import os
import logging
import fileinput
import time
from datetime import datetime
from optparse import OptionParser

import misc.fileUtils
import misc.utils as utils


# Interact with FITS files
import pyfits
import numpy as np
import datahandler


# Logging
from misc.paLog import log



class MEF:
    """
    \brief Class used to run some operations with MEF (Multi-Extension FITS) frames
    
    \par Class:
         MEF
    \par Purpose:
         Basic operations with MEF files
    \par Description:
         All PANIC observational data will be recorded in the so-called multi-extension FITS format. A MEF file is comprised of several segments called Header/Data Units (HDUs). Every HDU consists of an Header Unit (the well known FITS headers) in ASCII format followed by an optional Data Unit. The first HDU is called the primary, and any number of additional HDUs may follow. These additional HDUs are referred to as FITS extensions.

         In the PANIC FITS, the primary HDU only contains ASCII header cards describing the observation, but no data. The astronomical data arrays are stored in additional image extensions. There are 4 image extensions , 1 for each detector in the mosaic.
        
    \par Language:
        Python, PyFITS
    \param input_files
        A list FITS files
    \param output_filename_suffix
        Suffix to add to the outfile
    \retval 0
        If no error
    \author
        JMIbannez, IAA-CSIC
        
    """
    def __init__(self,  input_files ):
        """ class initialization """
                 
        self.input_files = input_files
            
    def doJoin(self, output_filename_suffix=".join.fits"):
        """
        \brief Method used to join a MEF into a single FITS frames, coping all the header information required
        
        NOTES:
      
        -An alternative way is to use SWARP to create a single FITS frame ( swarp file_with_ext.fits )
         However, this alternative takes longer.

        -Other alternative is iraf.mscjoin 
        
        TODO:
         - deduce the RA,DEC coordinates of the pointing !!
        """
           
        log.info("Starting JoinMEF")
         
        for file in self.input_files:        
            try:
                hdulist = pyfits.open(file)
            except IOError:
                print 'Error, can not open file %s' %(file)
                return 0
            
            try:
                if hdulist[0].header['EXTEND']!=True:
                    print 'Error, file %s is not a MEF file' %(file)
                    return 0
            except KeyError:
                print 'Error, file %s is not a MEF file' %(file)
                return 0
            
            try:
                next=hdulist[0].header['NEXTEND']
            except KeyError:
                print 'Warning, card NEXTEND not found. Counting number of extensions...'
                next=0
                while (1):
                    try:
                        if (hdulist[next+1].header['XTENSION'] == 'IMAGE'): next += 1
                    except:
                        break
            print "->Found %d extensions" %(next)
                         
            out_filename=file.replace(".fits", output_filename_suffix)
            width=naxis1=2048
            height=naxis2=2048
            temp12=np.zeros((height,width*2), dtype=np.float32)
            temp34=np.zeros((height,width*2), dtype=np.float32)
            for i in range(0, height):
                # Q1 i-row
                temp12[i,0 : width]=hdulist[1].data[i, 0 : width]
                # Q2 i-row
                temp12[i, width: 2*width]=hdulist[2].data[i, 0 : width]
                # Q3 i-row
                temp34[i, 0 : width]=hdulist[3].data[i, 0 : width]
                # Q4 i-row
                temp34[i, width : 2*width]=hdulist[4].data[i, 0 : width]

            joined_data = np.append(temp12, temp34).reshape(4096,4096)
            hdu = pyfits.HDUList([pyfits.PrimaryHDU(header=hdulist[0].header, data=joined_data)])
            #hdu.verify('silentfix')
            # now, copy extra keywords required
            try:
                hdu[0].header.update("BITPIX",-32)
                hdu[0].header.update("NAXIS1",4096)
                hdu[0].header.update("NAXIS2",4096)
                # TODO: deduce RA,DEC pointing coordinates 
            except KeyError:
                print 'Warning, some key cannot not be copied'
                            
            hdu.writeto(out_filename, output_verify='ignore', clobber=True)
            hdu.close(output_verify='ignore')
            del hdu
            print "File %s created " %(out_filename)
        
        log.info("End of JoinMEF. %d files created", next)
        
        return next, out_filename
            
    def doSplit( self , out_filename_suffix=".Q%02d.fits", out_dir=None, copy_keyword=['DATE','OBJECT','DATE-OBS','RA','DEC','EQUINOX','RADECSYS','UTC','LST','UT','ST','AIRMASS','IMAGETYP','EXPTIME','TELESCOP','INSTRUME','MJD-OBS','FILTER2']):
        """ 
        Method used to split a MEF into single FITS frames, coping all the header information required                           
        """
        log.info("Starting SplitMEF")
        
        out_filenames=[]
        n=0 
        for file in self.input_files:        
            try:
                hdulist = pyfits.open(file)
            except IOError:
                print 'Error, can not open file %s' %(file)
                raise Exception("Error, can not open file %s"%file)
            
            try:
                if hdulist[0].header['EXTEND']!=True:
                    print 'Error, file %s is not a MEF file' %(file)
                    raise Exception("Error, file %s is not a MEF file"%(file))
            except KeyError:
                raise Exception("Error, file %s is not a MEF file"%(file))
            
            try:
                next=hdulist[0].header['NEXTEND']
            except KeyError:
                print 'Warning, card NEXTEND not found. Counting number of extensions...'
                next=0
                while (1):
                    try:
                        if (hdulist[next+1].header['XTENSION'] == 'IMAGE'): next += 1
                    except:
                        break
                print "->Found %d extensions" %(next)
                         
            
            for i in range(1,next+1):
                suffix=out_filename_suffix%i
                new_filename = file.replace(".fits", suffix)
                if out_dir!=None: new_filename=new_filename.replace(os.path.dirname(new_filename), out_dir) 
                out_filenames.append(new_filename)
                out_hdulist = pyfits.HDUList([pyfits.PrimaryHDU(header=hdulist[i].header, data=hdulist[i].data)])
                out_hdulist.verify('silentfix')
                # now, copy extra keywords required
                for key in copy_keyword:
                    try:
                        value=hdulist[0].header.ascardlist()[key].value
                        comment=hdulist[0].header.ascardlist()[key].comment
                        out_hdulist[0].header.update(key,value,comment)
                        # We DON'T need to update RA,DEC (pointing coordinates), because each 
                        # extension should have CRVAL/CRPIX values!!
                    except KeyError:
                        print 'Warning, key %s cannot not be copied, is not in the header' %(key)
                # delete some keywords not required anymore
                del out_hdulist[0].header['EXTNAME']                
                out_hdulist.writeto(out_filenames[n], output_verify='ignore', clobber=True)
                out_hdulist.close(output_verify='ignore')
                del out_hdulist
                print "File %s created " %(out_filenames[n])
                n+=1
                
            
        log.info("End of SplitMEF. %d files created", n)
        return next , out_filenames
                    
    def createMEF( self, output_file=os.getcwd()+"/mef.fits" , primaryHeader=None):
        """ 
        Method used to create a MEF from a set of n>0 FITS frames                           
        """
        
        log.info("Starting createMEF")
         
        
        # Add primary header to output file...
        #prihdr = myflat[0].header.copy()
        fo = pyfits.HDUList()
        prihdu=pyfits.PrimaryHDU(data=None, header=primaryHeader)
        # Start by updating PRIMARY header keywords...
        prihdu.header.update('EXTEND', pyfits.TRUE,after='NAXIS')
        prihdu.header.update('NEXTEND', 0)
        prihdu.header.update('FILENAME', output_file)
        
        fo.append(prihdu)
        nExt=0
        for file in self.input_files:        
            try:
                f = pyfits.open(file)
            except IOError:
                print 'Error, can not open file %s' %(file)
                raise Exception("Error, can not open file %s"%file)
            
            # Check if is a MEF, thus
            try:
                if f[0].header['EXTEND']==True:
                    print 'Error, file %s is already a MEF file. Only can create a MEF from single FITS' %(file)
                    f.close()
                    del prihdu
                    del fo
                    raise Exception("Error, file %s is already a MEF file"%(file))
            except KeyError:
                #ok, we suppose file is not a MEF
                pass                        
            hdu = pyfits.ImageHDU(data=f[0].data, header=f[0].header)
            #hdu.header.update('EXTVER',1)
            fo.append(hdu)
            del hdu
            nExt+=1
       
        prihdu.header.update('NEXTEND',nExt)
        misc.fileUtils.removefiles(output_file)
        fo.writeto(output_file, output_verify='ignore')
        fo.close(output_verify='ignore')
        del fo            
        
        log.info("MEF file %s created"%(output_file))
        
        return nExt,output_file
                                  
################################################################################
# main
if __name__ == "__main__":
    
    # Get and check command-line options
    usage = "usage: %prog [options] arg1 arg2"
    parser = OptionParser(usage)
    
    parser.add_option("-f", "--file",
                  action="store", dest="file",
                  help="Input MEF file. It has to be a fullpath file name")
    
    parser.add_option("-l", "--input",
                  action="store", dest="input_file_list",
                  help="Source file list of data frames. It has to be a fullpath file name")
    
    parser.add_option("-s", "--suffix",
                  action="store", dest="out_suffix", help="suffix to out files (default .%02d.fits)")
    
    parser.add_option("-J", "--join",
                  action="store_true", dest="join", help="make a join/stitch of the FITS extensions creating a single FITS file", default=False)
                                 
    parser.add_option("-S", "--split",
                  action="store_true", dest="split", help="make a split of MEF files, adding a suffix for each extension", default=False)
    
    parser.add_option("-C", "--create",
                  action="store_true", dest="create", help="create a MEF (with N extensions) from a set N of single FITS files", default=False)
                  
    
    (options, args) = parser.parse_args()
    
    if options.file:
        filelist=[options.file]
    elif options.input_file_list:
        filelist=[line.replace( "\n", "") for line in fileinput.input(options.input_file_list)]
        print filelist
    else:
        parser.print_help()
        parser.error("incorrect number of arguments " )
        
    myMEF = MEF(filelist)
    if options.join:
        if not options.out_suffix: options.out_suffix=".join.fits"
        myMEF.doJoin( options.out_suffix )
    elif options.split:
        if not options.out_suffix: options.out_suffix=".%02d.fits"
        myMEF.doSplit( options.out_suffix )
    elif options.create:
        myMEF.createMEF(os.getcwd()+"/mef.fits")      
        
        