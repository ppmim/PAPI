#!/usr/bin/env python
"""Module to do some useful operations with Multi-Extension FITS files"""

__version__ = "$Revision: 64bda015861d $"
# $Source$

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

import os
import fileinput
from optparse import OptionParser

import misc.fileUtils
import misc.utils as utils


# Interact with FITS files
import pyfits
import numpy as np
import pywcs
import numpy

# Logging
from misc.paLog import log

class MEF_Exception (Exception): pass


class MEF (object):
    """
    @summary:  Class used to run some operations with MEF (Multi-Extension FITS) 
    frames
    
    
     All PANIC observational data will be recorded in the so-called 
     multi-extension FITS format. A MEF file is comprised of several segments 
     called Header/Data Units (HDUs). 
     Every HDU consists of an Header Unit (the well known FITS headers) in 
     ASCII format followed by an optional Data Unit. The first HDU is called 
     the primary, and any number of additional HDUs may follow. These additional 
     HDUs are referred to as FITS extensions.

     In the PANIC FITS, a MEF file contains 5 FITS HDU's. The first HDU, 
     called primary HDU, only contains ASCII header cards describing the 
     observation and common to all detectors, but no data. The astronomical 
     data arrays are stored in 4 additional HDU's, called image extensions. 
     Thus, there are 4 image extensions, 1 for each detector in the mosaic.
    
    @note:  Language:
        Python, PyFITS
    @author: JMIbannez, IAA-CSIC
        
    """
    def __init__(self,  input_files ):
        """ class initialization """
                 
        self.input_files = input_files
            
    def doJoin(self, output_filename_suffix = ".join.fits", output_dir = None):
        """
        @summary:  Method used to join a MEF into a single FITS frames, 
        coping all the header information required
        
        
        @param param: output_filename_suffix : suffix for the new joined FITS files.
        @param param: output_dir : directory name where the join-fits file must 
                                    be created. If None given, the same directory 
                                    of the source files will be used.
                
        @note: 
      
        - An alternative way is to use SWARP to create a single FITS frame 
        ( swarp file_with_ext.fits )
         However, this alternative takes longer.

        - Other alternative is iraf.mscjoin 
        
        @todo: deduce the RA,DEC coordinates of the pointing !!
        """
           
        log.info("Starting JoinMEF")
        
            
        for file in self.input_files:        
            if output_dir == None:
                output_dir = os.path.dirname(file)
            try:
                hdulist = pyfits.open(file)
            except IOError:
                print 'Error, can not open file %s' % (file)
                return 0
            
            #Check if is a MEF file 
            if len(hdulist)>1:
                n_ext = len(hdulist)-1
                log.debug("Found a MEF file with %d extensions", n_ext)
            else:
                n_ext = 1
                log.error("Found a simple FITS file. Join operation only\
                 supported for MEF files")
                
                         
            out_filename = output_dir + "/" + \
                os.path.basename(file).replace(".fits", output_filename_suffix)
            width = 2048
            height =  2048
            temp12 = np.zeros((height, width*2), dtype = np.float32)
            temp34 = np.zeros((height, width*2), dtype = np.float32)
            for i in range(0, height):
                # Q1 i-row
                temp12[i, 0 : width] = hdulist[1].data[i, 0 : width]
                # Q2 i-row
                temp12[i, width: 2*width] = hdulist[2].data[i, 0 : width]
                # Q3 i-row
                temp34[i, 0 : width] = hdulist[3].data[i, 0 : width]
                # Q4 i-row
                temp34[i, width : 2*width] = hdulist[4].data[i, 0 : width]

            joined_data = np.append(temp12, temp34).reshape(4096, 4096)
            hdu = pyfits.HDUList([pyfits.PrimaryHDU(header=hdulist[0].header, data=joined_data)])
            #hdu.verify('silentfix')
            # now, copy extra keywords required
            try:
                hdu[0].header.update("BITPIX", -32)
                hdu[0].header.update("NAXIS1", 4096)
                hdu[0].header.update("NAXIS2", 4096)
                # TODO: deduce RA,DEC pointing coordinates 
            except KeyError:
                print 'Warning, some key cannot not be copied'
                            
            hdu.writeto(out_filename, output_verify = 'ignore', clobber = True)
            hdu.close(output_verify = 'ignore')
            del hdu
            print "File %s created " % (out_filename)
        
        log.info("End of JoinMEF. %d files created", n_ext)
        
        return n_ext, out_filename
            
    def doSplit( self , out_filename_suffix = ".Q%02d.fits", out_dir = None, 
                 copy_keyword = None):
        """ 
        @summary: Method used to split a MEF into single FITS frames, 
        coping all the header information required                           
        """
        log.info("Starting SplitMEF")
        
        if copy_keyword == None:
            copy_keyword = ['DATE', 'OBJECT', 'DATE-OBS', 'RA', 'DEC', 'EQUINOX', 
                    'RADECSYS', 'UTC', 'LST', 'UT', 'ST', 'AIRMASS', 'IMAGETYP', 
                    'EXPTIME', 'TELESCOP', 'INSTRUME', 'MJD-OBS', 'FILTER2']
                
        out_filenames = []
        n = 0 
        for file in self.input_files:        
            try:
                hdulist = pyfits.open(file)
            except IOError:
                print 'Error, can not open file %s' % (file)
                raise MEF_Exception ("Error, can not open file %s" % file)
            
            #Check if is a MEF file 
            if len(hdulist)>1:
                n_ext = len(hdulist)-1
                log.debug("MEF file with %d extensions", n_ext)
            else:
                n_ext = 1
                log.error("Warning, Found a simple FITS file")
            
            for i in range (1, n_ext+1):
                suffix = out_filename_suffix % i
                new_filename = file.replace(".fits", suffix)
                if out_dir != None: 
                    new_filename = new_filename.replace( 
                                    os.path.dirname(new_filename), out_dir
                                    ) 
                out_filenames.append(new_filename)
                out_hdulist = pyfits.HDUList( [pyfits.PrimaryHDU( 
                                header = hdulist[i].header,
                                data = hdulist[i].data)]
                                             )
                
                out_hdulist.verify('silentfix')
                # now, copy extra keywords required
                for key in copy_keyword:
                    try:
                        value = hdulist[0].header.ascardlist ()[key].value
                        comment = hdulist[0].header.ascardlist ()[key].comment
                        out_hdulist[0].header.update (key, value,comment)
                        # We DON'T need to update RA, DEC (pointing coordinates), because each 
                        # extension should have CRVAL/CRPIX values!!
                    except KeyError:
                        print 'Warning, key %s cannot not be copied, is not in the header' % (key)
                # delete some keywords not required anymore
                del out_hdulist[0].header['EXTNAME']                
                out_hdulist.writeto (out_filenames[n], output_verify = 'ignore', clobber = True)
                out_hdulist.close (output_verify = 'ignore')
                del out_hdulist
                print "File %s created " % (out_filenames[n])
                n += 1
                
            
        log.info("End of SplitMEF. %d files created", n)
        return n_ext , out_filenames
                    
    def createMEF (self, output_file = os.getcwd()+"/mef.fits" ,
                   primaryHeader = None):
        """ 
        @summary: Method used to create a MEF from a set of n>0 FITS frames                           
        """
        
        log.info("Starting createMEF")
         
        
        # Add primary header to output file...
        #prihdr = myflat[0].header.copy()
        fo = pyfits.HDUList()
        prihdu = pyfits.PrimaryHDU (data = None, header = primaryHeader)
        # Start by updating PRIMARY header keywords...
        prihdu.header.update('EXTEND', pyfits.TRUE, after = 'NAXIS')
        prihdu.header.update('NEXTEND', 0)
        prihdu.header.update('FILENAME', output_file)
        
        fo.append(prihdu)
        n_ext = 0
        for file in self.input_files:        
            try:
                f = pyfits.open(file)
            except IOError:
                print 'Error, can not open file %s' %(file)
                raise MEF_Exception ("Error, can not open file %s"%file)
            
            #Check if is a MEF file 
            if len(f)>1:
                mef = True
                n_ext = len(f)-1
                
                log.error ("Found a MEF file with %d extensions.\
                 Operation only supported for simple FITS", n_ext)
                
                raise MEF_Exception ("Found a MEF file with %d extensions.\
                Operation only supported for simple FITS"%n_ext)
            else:
                mef = False
                n_ext = 1
                log.debug("Simple FITS file")
                                        
            hdu = pyfits.ImageHDU (data = f[0].data, header = f[0].header)
            #hdu.header.update('EXTVER',1)
            fo.append(hdu)
            del hdu
            n_ext += 1
       
        prihdu.header.update ('NEXTEND', n_ext)
        misc.fileUtils.removefiles (output_file)
        fo.writeto (output_file, output_verify ='ignore')
        fo.close (output_verify = 'ignore')
        del fo            
        
        log.info ("MEF file %s created" % (output_file))
        
        return n_ext, output_file
    
    def convertGEIRSToMEF( self, out_filename_suffix = ".mef.fits",
                           out_dir = None, copy_keyword = []):
        """ 
        @summary: Method used to convert a single FITS file (PANIC-GEIRS v0) 
        having a 4-detector-frame to 1-MEF with
        a 4-extension FITS, one per each frame.

        @attention: it is NOT valid for cubes of data, by the moment
        
        @param out_filename_suffix: suffix added to the original input filename
        @param out_dir: directory where the new output file will be created
        @param copy_keyword: list of keyword from the PrimaryHDU of the original file to be copied
        to the other header extensions  
        @return: the number of extensions per MEF and the list of output files (MEFs) created                               

        @author: jmiguel@iaa.es
        """
        
        log.info("Starting convertGEIRSToMEF")
        n = 0
        out_filenames = []
        n_ext = 4 # number of extension expected

        # Iterate over the whole input file list
        for file in self.input_files:        
            try:
                in_hdulist = pyfits.open(file)
            except IOError:
                log.error('Error, can not open file %s', file)
                raise MEF_Exception ("Error, can not open file %s" % file)
            
            #Check if is a MEF file 
            if len(in_hdulist)>1:
                log.error("Found a MEF file with %d extensions.\
                 Cannot convert to MEF file", len(in_hdulist)-1)
                raise MEF_Exception("Cannot convert an already MEF to MEF file")
            else:
                log.info("OK, found a single FITS file")
                
                
            # copy primary header from input file
            primaryHeader = in_hdulist[0].header.copy()
            out_hdulist = pyfits.HDUList()
            
            # Create primary HDU (without data, only the common header)
            prihdu = pyfits.PrimaryHDU(data=None, header = primaryHeader)
            # Start by updating PRIMARY header keywords...
            prihdu.header.update ('EXTEND', pyfits.TRUE, after = 'NAXIS')
            prihdu.header.update ('NEXTEND', n_ext)
            out_hdulist.append (prihdu)
            
            new_filename = file.replace(".fits", out_filename_suffix)
            if out_dir != None: 
                new_filename = new_filename.replace (os.path.dirname (new_filename), out_dir) 
            out_filenames.append (new_filename)
            
            # Read all image sections (n_ext frames) and create the associated HDU
            pix_centers = numpy.array ([[1024, 1024], [3072, 1024], 
                                        [3072, 3072], [1024, 3072]], numpy.float_)
            for i in range (0, n_ext/2):
                for j in range (0, n_ext/2):
                    log.debug("Reading %d-quadrant ..." % (i + j))
                    hdu_data_i = in_hdulist[0].data[2048*i:2048*(i+1), 2048*j:2048*(j+1)]
                    hdu_i = pyfits.ImageHDU (data = hdu_data_i)
                    log.debug("Data size of %d-quadrant = %s" % (i, hdu_data_i.shape))
                    
                    # Create the new WCS
                    try:
                        orig_ar = float(primaryHeader['RA'])
                        orig_dec = float(primaryHeader['DEC'])
                    except ValueError:
                        # no ar,dec values in the header, then can't re-compute ra,dec coordinates 
                        # or update the wcs header
                        pass
                    else:
                        #due to PANICv0 header hasn't a proper WCS header, we built a basic one 
                        #in order to computer the new ra, dec coordinates
                        new_wcs = pywcs.WCS(primaryHeader)
                        new_wcs.wcs.crpix = [primaryHeader['NAXIS1']/2, primaryHeader['NAXIS2']/2]  
                        new_wcs.wcs.crval =  [ primaryHeader['RA'], primaryHeader['DEC'] ]
                        new_wcs.wcs.ctype = ['RA-TAN', 'DEC-TAN']
                        new_wcs.wcs.cunit = ['deg', 'deg']
                        pix_scale = 0.45 # arcsec per pixel
                        new_wcs.wcs.cd = [[-pix_scale/3600.0, 0], [pix_scale/3600.0, 0]]

                        new_pix_center = new_wcs.wcs_pix2sky ([pix_centers[i*2+j]], 1)
                        
                        prihdu.header.update ('RA', new_pix_center[0][0])
                        prihdu.header.update ('DEC', new_pix_center[0][1])
                        
                        # Now update the new-wcs for the new subframe header
                        hdu_i.header.update ('CRPIX1', 1024)
                        hdu_i.header.update ('CRPIX2', 1024)
                        hdu_i.header.update ('CRVAL1', new_pix_center[0][0])
                        hdu_i.header.update ('CRVAL2', new_pix_center[0][1])
                        hdu_i.header.update ('CD1_1', -pix_scale/3600.0, "Axis rotation & scaling matrix")
                        hdu_i.header.update ('CD1_2', 0, "Axis rotation & scaling matrix")
                        hdu_i.header.update ('CD2_1', 0, "Axis rotation & scaling matrix")
                        hdu_i.header.update ('CD2_2', pix_scale/3600.0, "Axis rotation & scaling matrix")
                        hdu_i.header.update ('CTYPE1' , 'RA-TAN') 
                        hdu_i.header.update ('CTYPE2' , 'DEC-TAN')
                        hdu_i.header.update ('CUNIT1', 'deg')
                        hdu_i.header.update ('CUNIT2', 'deg')
                        
                    # now, copy extra keywords required
                    for key in copy_keyword:
                        try:
                            value = in_hdulist[0].header.ascardlist()[key].value
                            comment = in_hdulist[0].header.ascardlist()[key].comment
                            hdu_i.header.update (key, value, comment)
                            # We DON'T need to update RA,DEC (pointing coordinates), because each 
                            # extension should have CRVAL/CRPIX values!!
                        except KeyError:
                            print 'Warning, key %s cannot not be copied, is not in the header' % (key)
                    # Append new HDU to MEF
                    out_hdulist.append (hdu_i)
                    out_hdulist.verify ('ignore')
            
            # Now, write the new MEF file
            #out_hdulist[0].header.update('NEXTEND', 4)
            out_hdulist.writeto(new_filename, output_verify = 'ignore', clobber=True)
            out_hdulist.close(output_verify = 'ignore')
            del out_hdulist
            log.info("MEF file %s created" % (out_filenames[n]))
            n += 1
        
        return n_ext, out_filenames
    
    def splitGEIRSToSimple( self, out_filename_suffix = ".Q%02d.fits", out_dir = None):
        """ 
        @summary: Method used to convert a single FITS file (PANIC-GEIRS v0) 
        having a 4-detector-frame to 4 single 
        FITS files with a frame per file. 

        @param out_filename_suffix: suffix added to the original input filename
        @param out_dir: directory where the new output file will be created
        @return: the list of output files (MEFs) created
        
        @attention: it is NOT valid for cubes of data, by the moment                           
        
        @bug: ar,dec need to be recalculated !!!!!
        
        @author: jmiguel@iaa.es
        """
        
        log.info("Starting splitGEIRSToSimple")
        n = 0 
        out_filenames = []

        # Iterate over the whole input file list 
        for file in self.input_files:        
            try:
                in_hdulist = pyfits.open(file)
            except IOError:
                log.error('Error, can not open file %s', file)
                raise MEF_Exception("Error, can not open file %s" % file)
            
            #Check if is a MEF file 
            if len(in_hdulist) > 1:
                log.error("Found a MEF file with %d extensions. Cannot convert", len(in_hdulist)-1)
                raise MEF_Exception("Error, found a MEF file, expected a single FITS ")
            else:
                log.info("OK, found a single FITS file")
 
            # copy primary header from input file
            primaryHeader = in_hdulist[0].header.copy()
            
            # Read all image sections (4 frames) and create the associated 4-single FITS files
            n_ext = 4
            pix_centers = numpy.array ([[1024, 1025], [3072, 1024], 
                                        [3072, 3072], [1024, 3072]], numpy.float_)
            for i in range (0, n_ext/2):
                for j in range (0, n_ext/2):
                    log.debug("Reading %d-quadrant ..." % (i*2 + j))
                    hdu_data_i = in_hdulist[0].data[2048*i:2048*(i+1), 2048*j:2048*(j+1)]
                    log.debug("Data size of %d-quadrant = %s" % (i*2+j, hdu_data_i.shape))    
                    # Create primary HDU (data + header)
                    out_hdulist = pyfits.HDUList()               
                    prihdu = pyfits.PrimaryHDU (data = hdu_data_i, header = primaryHeader)
                    # Start by updating PRIMARY header keywords...
                    prihdu.header.update ('EXTEND', pyfits.FALSE, after = 'NAXIS')
                    # AR,DEC (WCS !!) need to be re-calculated !!!
                    
                    # Create the new WCS
                    try:
                        orig_ar = float(primaryHeader['RA'])
                        orig_dec = float(primaryHeader['DEC'])
                    except ValueError:
                        # no ar,dec values in the header, then can't re-compute ra,dec coordinates 
                        # or update the wcs header
                        pass
                    else:
                        #due to PANICv0 header hasn't a proper WCS header, we built a basic one 
                        #in order to computer the new ra, dec coordinates
                        new_wcs = pywcs.WCS(primaryHeader)
                        new_wcs.wcs.crpix = [primaryHeader['NAXIS1']/2, primaryHeader['NAXIS2']/2]  
                        new_wcs.wcs.crval =  [ primaryHeader['RA'], primaryHeader['DEC'] ]
                        new_wcs.wcs.ctype = ['RA-TAN', 'DEC-TAN']
                        new_wcs.wcs.cunit = ['deg', 'deg']
                        pix_scale = 0.45 # arcsec per pixel
                        new_wcs.wcs.cd = [[-pix_scale/3600.0, 0], [pix_scale/3600.0, 0]]

                        new_pix_center = new_wcs.wcs_pix2sky ([pix_centers[i*2+j]], 1)
                        
                        prihdu.header.update ('RA', new_pix_center[0][0])
                        prihdu.header.update ('DEC', new_pix_center[0][1])
                        
                        # Now update the new-wcs for the new subframe
                        prihdu.header.update ('CRPIX1', 1024)
                        prihdu.header.update ('CRPIX2', 1024)
                        prihdu.header.update ('CRVAL1', new_pix_center[0][0])
                        prihdu.header.update ('CRVAL2', new_pix_center[0][1])
                        prihdu.header.update ('CTYPE1' , 'RA-TAN') 
                        prihdu.header.update ('CTYPE2' , 'DEC-TAN')
                        prihdu.header.update ('CUNIT1', 'deg')
                        prihdu.header.update ('CUNIT2', 'deg')
                        prihdu.header.update ('CD1_1', -pix_scale/3600.0, "Axis rotation & scaling matrix")
                        prihdu.header.update ('CD1_2', 0, "Axis rotation & scaling matrix")
                        prihdu.header.update ('CD2_1', 0,  "Axis rotation & scaling matrix")
                        prihdu.header.update ('CD2_2', pix_scale/3600.0, "Axis rotation & scaling matrix")
                        
                    
                    prihdu.header.update ('PA_QUAD', (i*2)+j, "PANIC Quadrant number [0,1,2,3]")
                    
                    out_hdulist.append (prihdu)    
                    out_hdulist.verify ('ignore')
                
                    new_filename = file.replace(".fits", out_filename_suffix % (i*2+j))
                    if out_dir != None: 
                        new_filename = new_filename.replace (os.path.dirname(new_filename), out_dir) 
                    out_filenames.append (new_filename)
     
                    # Now, write the new MEF file
                    out_hdulist.writeto (new_filename, output_verify = 'ignore', clobber=True)
                    out_hdulist.close (output_verify = 'ignore')
                    del out_hdulist
                    log.info ("FITS file %s created" % (out_filenames[n]))
                    n += 1
        
        return out_filenames
                                  
################################################################################
# main
if __name__ == "__main__":
    
    # Get and check command-line options
    usage = "usage: %prog [options] arg1 arg2"
    parser = OptionParser (usage)
    
    parser.add_option ("-f", "--file",
                  action = "store", dest = "file",
                  help = "Input MEF file. It has to be a fullpath file name")
    
    parser.add_option ("-l", "--input",
                  action = "store", dest = "input_file_list",
                  help = "Source file list of data frames. \
                  It has to be a fullpath file name")
    
    parser.add_option ("-s", "--suffix",
                  action = "store", dest = "out_suffix", \
                  help = "suffix to out files (default .%02d.fits)")
    
    parser.add_option ("-J", "--join",
                  action = "store_true", dest = "join", \
                  help = "make a join/stitch of the FITS extensions \
                  creating a single FITS file", default = False)
                                 
    parser.add_option ("-S", "--split",
                  action = "store_true", dest = "split", \
                  help = "make a split of MEF files, adding a suffix \
                  for each extension", default = False)
    
    parser.add_option ("-G", "--geirs-split",
                  action = "store_true", dest = "geirs_split", \
                  help = "make a split of GEIRS-v0 files, creating \
                  4-single files and adding a suffix for each extension", \
                  default = False)
    
    parser.add_option ("-C", "--create",
                  action = "store_true", dest = "create", \
                  help = "create a MEF (with N extensions) from a \
                  set N single FITS files", default = False)
   
    parser.add_option ("-g", "--geirs-convert",
                  action = "store_true", dest = "geirs_convert", \
                  help = "convert a GEIRS-v0 file (with N extensions) to a \
                  MEF FITS file with 4 extensions", default = False)
                  
    
    (options, args) = parser.parse_args()
    
    if options.file:
        filelist = [options.file]
    elif options.input_file_list:
        filelist = [line.replace( "\n", "") for line in fileinput.input(options.input_file_list)]
        print filelist
    else:
        parser.print_help()
        parser.error("incorrect number of arguments " )
        
    myMEF = MEF(filelist)
    
    if options.join:
        if not options.out_suffix: 
            options.out_suffix = ".join.fits"
        myMEF.doJoin (options.out_suffix)
        
    elif options.split:
        if not options.out_suffix: 
            options.out_suffix = ".%02d.fits"
        myMEF.doSplit( options.out_suffix )
    
    elif options.create:
        myMEF.createMEF(os.getcwd()+"/mef.fits")
    
    elif options.geirs_split:
        myMEF.splitGEIRSToSimple()
    
    elif options.geirs_convert:
        myMEF.convertGEIRSToMEF()
        
              
        
        