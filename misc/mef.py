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
# Last update: 21/Jun/2013   jmiguel@iaa.es - some header improvements in
#                                             copy_keywords.
#              06/Jun/2014   jmiguel@iaa.es - 
# TODO
#    - identificar SGi en MEF y 4kx4k
#      
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
from astropy import wcs
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
    def __init__(self,  input_files , *a, **k):
        """ class initialization """
                 
        super (MEF, self).__init__ (*a,**k)         
        
        self.input_files = input_files
            
    def doJoin(self, output_filename_suffix = ".join.fits", output_dir = None):
        """
        Method used to join a MEF into a single FITS frames, copying all the 
        header information required.
        Basically used to run test with HAWK-I images, converting them to PANIC 
        format.
        
        Parameters
        ----------
        output_filename_suffix : str
            suffix for the new joined FITS files.
        output_dir : str
            directory name where the join-fits file must be created. If None 
            given, the same directory of the source files will be used.
        
        
        Returns
        -------
        Number of extensions and filename of new joined file created.
                
        Notes
        -----
      
        - An alternative way is to use SWARP to create a single FITS frame 
        ( swarp file_with_ext.fits )
         However, this alternative takes longer.

        - Other alternative is iraf.mscjoin 
        
        Todo
        ----
        Deduce the RA,DEC coordinates of the pointing !!
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
            hdu = pyfits.HDUList([pyfits.PrimaryHDU(header=hdulist[0].header, 
                                                    data=joined_data)])
            #hdu.verify('silentfix')
            # now, copy extra keywords required
            try:
                hdu[0].header.update("BITPIX", -32)
                hdu[0].header.update("NAXIS1", 4096)
                hdu[0].header.update("NAXIS2", 4096)
                # TODO: deduce RA,DEC pointing coordinates 
            except KeyError:
                log.warning("Some key cannot be copied into header")
            
            hdu[0].header.add_history("[MEF.doJoin] MEF created with files %s"%str(self.input_files))                            
            hdu.writeto(out_filename, output_verify = 'ignore', clobber = True)
            hdu.close(output_verify = 'ignore')
            del hdu
            print "File %s created " % (out_filename)
        
        log.info("End of JoinMEF. %d files created", n_ext)
        
        return n_ext, out_filename
            
    def doSplit( self , out_filename_suffix = ".Q%02d.fits", out_dir = None, 
                 copy_keyword = None):
        """ 
        Method used to split a MEF into single FITS frames, 
        copying all the header information required.                           
        """
        log.info("Starting SplitMEF")
        
        if copy_keyword == None:
            # Default set of keywords whether None were given.
            copy_keyword = ['DATE', 'OBJECT', 'DATE-OBS', 'RA', 'DEC', 'EQUINOX', 
                    'RADECSYS', 'UTC', 'LST', 'UT', 'ST', 'AIRMASS', 'IMAGETYP', 
                    'EXPTIME', 'TELESCOP', 'INSTRUME', 'MJD-OBS', 'NCOADDS',
                    'FILTER', 'FILTER1','FILTER2', 'HIERARCH ESO DET NDIT','NDIT',
                    'CASSPOS','PIXSCALE']
                
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
                log.error("Found a simple FITS file, not a MEF file")
                raise MEF_Exception("File %s is not a MEF" % file)
            
            for i in range (1, n_ext+1):
                suffix = out_filename_suffix % i # number from 1 to 4
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
                
                #
                # Copy all keywords of header[0] and not included at header[i>0]
                # TODO 
                #keywords  = [if key not in hdulist[1].header for key in hdulist[0].header.cards

                # now, copy extra keywords required
                for key in copy_keyword:
                    try:
                        value = hdulist[0].header.cards[key].value
                        comment = hdulist[0].header.cards[key].comment
                        if key=='HIERARCH ESO DET NDIT':
                            out_hdulist[0].header.set('NDIT', value, comment)
                        else:
                            out_hdulist[0].header.set(key, value,comment)
                        # We DON'T need to update RA, DEC (pointing coordinates), because each 
                        # extension should have CRVAL/CRPIX values!!
                    except KeyError:
                        log.warning("Key %s cannot be copied, is not in the header"%(key))
                
                out_hdulist[0].header.add_history("[MEF.doSplit] Image split from original MEF %s"%file) 
                # delete some keywords not required anymore
                del out_hdulist[0].header['EXTNAME']                
                out_hdulist.writeto(out_filenames[n], 
                        output_verify = 'ignore', clobber = True)
                out_hdulist.close(output_verify = 'ignore')
                del out_hdulist
                log.debug("File %s created"%(out_filenames[n]))
                n += 1
                
            
        log.info("End of SplitMEF. %d files created", n)
        return n_ext , out_filenames
                    
    def createMEF (self, output_file = os.getcwd()+"/mef.fits" ,
                   primaryHeader = None):
        """ 
        Method used to create a MEF from a set of n>0 FITS frames.                          
        """
        
        log.info("Starting createMEF")
         
        
        # Add primary header to output file...
        #prihdr = myflat[0].header.copy()
        fo = pyfits.HDUList()
        prihdu = pyfits.PrimaryHDU (data = None, header = primaryHeader)
        # Start by updating PRIMARY header keywords...
        prihdu.header.set('EXTEND', pyfits.TRUE, after = 'NAXIS')
        prihdu.header.set('NEXTEND', 0)
        prihdu.header.set('FILENAME', output_file)
        
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
                                        
            hdu = pyfits.ImageHDU(data = f[0].data, header = f[0].header)
            #hdu.header.update('EXTVER',1)
            fo.append(hdu)
            del hdu
            n_ext += 1
       
        prihdu.header.set('NEXTEND', n_ext)
        prihdu.header.add_history("[MEF.createMEF] MEF created with files %s"%str(self.input_files))
        misc.fileUtils.removefiles (output_file)
        fo.writeto (output_file, output_verify ='ignore')
        fo.close (output_verify = 'ignore')
        del fo            
        
        log.info ("MEF file %s created" % (output_file))
        
        return n_ext, output_file
    
    def convertGEIRSToMEF( self, out_filename_suffix = ".mef.fits",
                           out_dir = None, copy_keyword = []):
        """ 
        Method used to convert a single FITS file (PANIC-GEIRS v0) 
        having a full 4-detector-frame to a MEF with a 4-extension FITS, 
        one per each frame.

        It is **NOT valid for cubes** of data, by the moment.
        
        Parameters
        ----------

        out_filename_suffix: str 
            suffix added to the original input filename
        
        out_dir: str 
            directory where the new output file will be created
        
        copy_keyword: list 
            list of keyword from the PrimaryHDU of the original file 
            to be copied to the other header extensions.  
        
        Returns
        -------
        n_ext, out_filenames: The number of extensions per MEF and the list of 
                              output files (MEFs) created.                              

        Notes
        -----
        - For an image of the sky, the coordinate system representation 
        might be something like: 

            CTYPE1 = `RA---TAN'   ! r.a. in tangent plane projection
            CTYPE2 = `DEC--TAN'   ! dec. in tangent plane projection
            CRPIX1 = 400.0        ! reference pixel on first axis 
            CRPIX2 = 400.0        ! reference pixel on second axis 
            CRVAL1 = 180.01234    ! r.a. at pixel (400,400) in degrees 
            CRVAL2 = 30.98765     ! dec. at pixel (400,400) in degrees 
            CD1_1  = -2.777778E-4 ! change in RA per pixel along first
                                  ! axis evaluated at reference pixel (pix. scale)
            CD1_2 = 0.0           ! change in RA per pixel along second
                                  ! axis evaluated at reference pixel
            CD2_1 = 0.0           ! change in dec per pixel along first
                                  !axis evaluated at reference pixel
            CD2_2 = 2.777778E-4   ! change in dec per pixel along second
                                  ! axis evaluated at reference pixel
        
        Todo
        ----
        - Header is not well formed.
        - identify of SGi detector properly

        
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
                log.info("Found a single FITS file. Let's convert it to MEF...")
                
                
            # copy primary header from input file
            primaryHeader = in_hdulist[0].header.copy()
            out_hdulist = pyfits.HDUList()
            
            # Create primary HDU (without data, only the common header)
            prihdu = pyfits.PrimaryHDU(data=None, header = primaryHeader)
            # Start by updating PRIMARY header keywords...
            prihdu.header.set('EXTEND', pyfits.TRUE, after = 'NAXIS')
            prihdu.header.set('NEXTEND', n_ext, after = 'EXTEND')
            prihdu.header.add_history("[MEF.convertGEIRSToMEF] MEF created from original filename %s"%file)
            
            # In the Primary Header we do not need the WCS keywords, only RA,DEC
            keys_to_del=['CRPIX1','CRPIX2','CRVAL1','CRVAL2','CDELT1','CDELT2','CTYPE1','CTYPE2']
            for key in keys_to_del: del prihdu.header[key]
            #
            #prihdu.header.update ('RA', new_pix_center[0][0])
            #prihdu.header.update ('DEC', new_pix_center[0][1])
                        
            out_hdulist.append (prihdu)
            
            new_filename = file.replace(".fits", out_filename_suffix)
            if out_dir != None: 
                new_filename = new_filename.replace (os.path.dirname (new_filename), out_dir) 
            out_filenames.append (new_filename)
            
            # Read all image sections (n_ext frames) and create the associated HDU
            pix_centers = numpy.array ([[1024, 1024], [1024, 3072], [3072, 1024], 
                                        [3072, 3072] ], numpy.float_)

            # Taking into account the gap (167pix), we set the new CRPIXi values
            # for each extension, refered to the center of the focal plane.
            new_crpix_center = numpy.array ([[2132, 2132], [2132, -81], [-81, 2132], 
                                        [-81, -81] ], numpy.float_)

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
                        #No RA,DEC values in the header, then can't re-compute 
                        #the coordinates or update the wcs header.
                        pass
                    else:
                        #Due to PANICv0 header hasn't a proper WCS header, 
                        #we built a basic one in order to computer the new 
                        #ra, dec coordinates.
                        
                        # ----
                        # WCS object is not USED !!!
                        # Instead, WCS header is built from primaryHeader. 
                        chip_gap = 167
                        new_wcs = wcs.WCS(primaryHeader)
                        new_wcs.wcs.crpix = [primaryHeader['NAXIS1']/2 + chip_gap/2.0, 
                                             primaryHeader['NAXIS2']/2 + chip_gap/2.0]  
                        new_wcs.wcs.crval =  [ primaryHeader['RA'], 
                                              primaryHeader['DEC'] ]
                        new_wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
                        new_wcs.wcs.cunit = ['deg', 'deg']
                        try:
                            pix_scale = primaryHeader['PIXSCALE']
                        except Exception,e:
                            raise Exception("Cannot find PIXSCALE keyword")

                        # We suppose no rotation, and North is up and East at left.
                        # CD matrix (the CDi_j elements) encode the sky position angle,
                        # the pixel scale, and a possible flipping.
                        # CD1_1 is <0 because East is supposed at Left = flipX
                        # CD2_2 is >0 because North is supposed at Up
                        # In addition, it must be noted that:
                        # CD1_1 = cos(r), CD1_2 = sin(r), CD2_1 = -sin(r), CD2_2 = cos(r)
                        # r = clockwise rotation_angle 
                        new_wcs.wcs.cd = [[-pix_scale/3600.0, 0], 
                                          [0, pix_scale/3600.0]]

                        new_pix_center = new_wcs.wcs_pix2world([pix_centers[i*2+j]], 1)
                        # ----

                        # Now update the new-wcs for the new subframe header.
                        # In principle, this is the same header than is generated
                        # by GEIRS for MEF files (save -M).
                        hdu_i.header.set('CRPIX1', new_crpix_center[i*2+j][0])
                        hdu_i.header.set('CRPIX2', new_crpix_center[i*2+j][1])
                        hdu_i.header.set('CRVAL1', primaryHeader['RA'])
                        hdu_i.header.set('CRVAL2', primaryHeader['DEC'])
                        hdu_i.header.set('CD1_1', -pix_scale/3600.0, 
                                            "Axis rotation & scaling matrix")
                        hdu_i.header.set('CD1_2', 0, 
                                            "Axis rotation & scaling matrix")
                        hdu_i.header.set('CD2_1', 0, 
                                            "Axis rotation & scaling matrix")
                        hdu_i.header.set('CD2_2', pix_scale/3600.0, 
                                            "Axis rotation & scaling matrix")
                        hdu_i.header.set('CTYPE1' , 'RA---TAN') 
                        hdu_i.header.set('CTYPE2' , 'DEC--TAN')
                        hdu_i.header.set('CUNIT1', 'deg')
                        hdu_i.header.set('CUNIT2', 'deg')
                        
                        # Force BITPIX=-32
                        prihdu.scale('float32')

                        if 'DETXYFLIP' in primaryHeader:
                            # When DETXYFLIP=0, no flip; default
                            # When DETXYFLIP=1, Flip Horizontally (interchanges the left and the right) 
                            # When DETXYFLIP=2, Flip Vertically (turns the image upside down) 
                            detflipxy = primaryHeader['DETXYFLIP']
                            log.debug("Found DETXYFLIP in header =%s"%detflipxy)
                        else: 
                            detflipxy = 0

                    
                        if 'DETROT90' in primaryHeader:
                            log.debug("Found DETROT90 in header")
                            # When DETROT90=0, no rotation;
                            # When DETROT90=1, 90deg clockwise rotation
                            # When DETROT90=2, 180deg clockwise rotation (default
                            # for PANIC)
                            if primaryHeader['DETROT90']==0:
                                if (i*2+j)==0: 
                                    if detflipxy==0: det_id = 2
                                    elif detflipxy==1: det_id = 3
                                    elif detflipxy==2: det_id = 1
                                    elif detflipxy==3: det_id = 4
                                elif (i*2+j)==1: 
                                    if detflipxy==0: det_id = 1
                                    elif detflipxy==1: det_id = 4
                                    elif detflipxy==2: det_id = 2
                                    elif detflipxy==3: det_id = 3
                                elif (i*2+j)==2: 
                                    if detflipxy==0: det_id = 3
                                    elif detflipxy==1: det_id = 2
                                    elif detflipxy==2: det_id = 4
                                    elif detflipxy==3: det_id = 1
                                elif (i*2+j)==3:
                                    if detflipxy==0: det_id = 4
                                    elif detflipxy==1: det_id = 1 
                                    elif detflipxy==2: det_id = 3
                                    elif detflipxy==3: det_id = 2
                            elif primaryHeader['DETROT90']==1:
                                if (i*2+j)==0:
                                    if detflipxy==0: det_id = 3
                                    elif detflipxy==1: det_id = 4 
                                    elif detflipxy==2: det_id = 2
                                    elif detflipxy==3: det_id = 1
                                elif (i*2+j)==1:
                                    if detflipxy==0: det_id = 2
                                    elif detflipxy==1: det_id = 1
                                    elif detflipxy==2: det_id = 3
                                    elif detflipxy==3: det_id = 4                                
                                elif (i*2+j)==2: 
                                    if detflipxy==0: det_id = 4
                                    elif detflipxy==1: det_id = 3
                                    elif detflipxy==2: det_id = 1
                                    elif detflipxy==3: det_id = 2                                
                                elif (i*2+j)==3:
                                    if detflipxy==0: det_id = 1
                                    elif detflipxy==1: det_id = 2
                                    elif detflipxy==2: det_id = 4
                                    elif detflipxy==3: det_id = 3                               
                            elif primaryHeader['DETROT90']==2:
                                if (i*2+j)==0:
                                    if detflipxy==0: det_id = 4
                                    elif detflipxy==1: det_id = 1 
                                    elif detflipxy==2: det_id = 3
                                    elif detflipxy==3: det_id = 2                               
                                elif (i*2+j)==1:
                                    if detflipxy==0: det_id = 3
                                    elif detflipxy==1: det_id = 2 
                                    elif detflipxy==2: det_id = 4
                                    elif detflipxy==3: det_id = 1                               
                                elif (i*2+j)==2:
                                    if detflipxy==0: det_id = 1
                                    elif detflipxy==1: det_id = 4 
                                    elif detflipxy==2: det_id = 2
                                    elif detflipxy==3: det_id = 3                              
                                elif (i*2+j)==3: 
                                    if detflipxy==0: det_id = 2
                                    elif detflipxy==1: det_id = 3 
                                    elif detflipxy==2: det_id = 1
                                    elif detflipxy==3: det_id = 4
                            elif primaryHeader['DETROT90']==3:
                                if (i*2+j)==0:
                                    if detflipxy==0: det_id = 1
                                    elif detflipxy==1: det_id = 2 
                                    elif detflipxy==2: det_id = 4
                                    elif detflipxy==3: det_id = 3
                                elif (i*2+j)==1:
                                    if detflipxy==0: det_id = 4
                                    elif detflipxy==1: det_id = 2
                                    elif detflipxy==2: det_id = 1
                                    elif detflipxy==3: det_id = 3
                                elif (i*2+j)==2:
                                    if detflipxy==0: det_id = 2
                                    elif detflipxy==1: det_id = 1
                                    elif detflipxy==2: det_id = 3
                                    elif detflipxy==3: det_id = 4
                                elif (i*2+j)==3:
                                    if detflipxy==0: det_id = 3
                                    elif detflipxy==1: det_id = 4
                                    elif detflipxy==2: det_id = 2
                                    elif detflipxy==3: det_id = 1
                        else:
                            # Then, we suppose DETROT90=2, DETXYFLIP=0, and default for PANIC !
                            log.warning("No DETROT90 found, supposed DETROT90=2 and DETXYFLIP=0")
                            if (i*2+j)==0: det_id = 4
                            elif (i*2+j)==1: det_id = 3
                            elif (i*2+j)==2: det_id = 1
                            elif (i*2+j)==3: det_id = 2                       

                        hdu_i.header.set('DET_ID', "SG%i"%det_id, 
                                            "PANIC Detector id SGi [i=1..4]")
                        hdu_i.header.set('EXTNAME',"SG%i_1"%det_id)
                        hdu_i.header.set('DETSEC',
                            "[%i:%i,%i:%i]"%(2048*i+1, 2048*(i+1), 2048*j+1, 2048*(j+1) ))

                        
                    # now, copy extra keywords required
                    for key in copy_keyword:
                        try:
                            value = in_hdulist[0].header.cards[key].value
                            comment = in_hdulist[0].header.cards[key].comment
                            if key=='HIERARCH ESO DET NDIT':
                                hdu_i.header.set ('NDIT', value, comment)
                            else:
                                hdu_i.header.set (key, value, comment)
                            # We DON'T need to update RA,DEC (pointing coordinates), because each 
                            # extension should have CRVAL/CRPIX values!!
                        except KeyError:
                            log.warning("Key %s cannot be copied, is not in the header"%(key))
                    # Append new HDU to MEF
                    # Force BITPIX=-32
                    hdu_i.scale('float32')
                    out_hdulist.append(hdu_i)
                    out_hdulist.verify('ignore')
            
            # Now, write the new MEF file
            #out_hdulist[0].header.update('NEXTEND', 4)
            out_hdulist.writeto(new_filename, output_verify = 'ignore', 
                                clobber=True)
            out_hdulist.close(output_verify = 'ignore')
            del out_hdulist
            log.info("MEF file %s created" % (out_filenames[n]))
            n += 1
        
        return n_ext, out_filenames
    
    def splitGEIRSToSimple( self, out_filename_suffix = ".Q%02d.fits", 
                            out_dir = None):
        """ 
        Method used to convert a single FITS file (PANIC-GEIRS v0) 
        having a full (4kx4k) 4-detector-frame to 4 single FITS files with a 
        frame per file. 
        Header is fully copied from original file, and added new WCS keywords.
        
        With FITS extensions (MEF) is easy to identify each detector, i.e.,  
        for rotation=2=180deg (DETROT90=2) we have:

        ext1 - [0:2048, 0:2048]      - SG4
        ext2 - [2048:4096, 0:2048]   - SG1
        ext3 - [0:2048, 2048:4096]   - SG3
        ext4 - [2048:4096,2048:4096] - SG2


        But, for single FITS files, I don't know how to identify the detectors 
        on the 4kx4k image; I did not find any keyword (DETROT90) about it. 
        
        DETROT90 = N  (where N=0,1,2,3)

        or much more clearly (for DETROT90=2)

        CHIP_SG1 = '[2048:4096, 0:2048]'
        CHIP_SG2 = '[2048:4096,2048:4096]'
        CHIP_SG3 = '[0:2048, 2048:4096]'
        CHIP_SG4 = '[0:2048,0:2048]' 

        Where SG = science grade detector (Hw id of the detector).

        Parameters
        ----------- 
        out_filename_suffix: suffix added to the original input filename.

        out_dir: directory where the new output file will be created.
        
        Returns
        -------

        n_ext,out_filenames: The number of files and the list of output files 
        created. The numbering may not match with DET_ID numbering (SGi).

                    DETROT90=2     DETROT90=0
        .Q01.fits --> SG4             SG2
        .Q02.fits --> SG3             SG1
        .Q03.fits --> SG1             SG3
        .Q04.fits --> SG2             SG4

        (DETROT90=2=180dec clockwise)
        
            
        Notes
        -----
        - It is NOT valid for cubes of data, by the moment                           
        - The enumeration order of the quadrants read is:
        
        |------------|
        |  Q2  | Q4  |
        | ---- |-----|
        |  Q1  | Q3  |
        |----- |-----|
        
        """
        
        log.info("Starting splitGEIRSToSimple")
        n = 0 
        out_filenames = []

        # Iterate over the whole input file list 
        for file in self.input_files:
            log.debug("Splitting file %s"%file)        
            try:
                in_hdulist = pyfits.open(file)
            except IOError:
                log.error('Error, can not open file %s', file)
                raise MEF_Exception("Error, can not open file %s" % file)
            
            # Check if is a MEF file 
            if len(in_hdulist) > 1:
                log.error("Found a MEF file with %d extensions. Cannot convert", 
                          len(in_hdulist)-1)
                raise MEF_Exception("Error, found a MEF file, expected a single FITS ")
            else:
                log.info("OK, found a single FITS file")
     
            if in_hdulist[0].header['NAXIS1']!=4096 or in_hdulist[0].header['NAXIS2']!=4096:
                log.error('Error, file %s is not a full frame image', file)
                raise MEF_Exception("Error, file %s is not a full frame image" % file)

            # copy primary header from input file
            primaryHeader = in_hdulist[0].header.copy()
            
            # Read all image sections (4 frames) and create the associated 
            # 4-single FITS files.
            n_ext = 4
            pix_centers = numpy.array ([[1024, 1024], [1024, 3072], 
                                        [3072, 1024], [3072, 3072]], 
                                       numpy.float_)

            # Taking into account the gap (167pix), we set the new CRPIXi values
            # for each extension, refered to the center of the focal plane.
            new_crpix_center = numpy.array ([[2132, 2132], [2132, -81], [-81, 2132], 
                                        [-81, -81] ], numpy.float_)
            for i in range (0, n_ext/2):
                for j in range (0, n_ext/2):
                    log.debug("Reading %d-quadrant ..." % (i*2 + j))
                    hdu_data_i = in_hdulist[0].data[2048*i:2048*(i+1), 
                                                    2048*j:2048*(j+1)]
                    log.debug("Data size of %d-quadrant = %s" % (i*2+j, 
                                                                 hdu_data_i.shape))    
                    # Create primary HDU (data + header)
                    out_hdulist = pyfits.HDUList()               
                    prihdu = pyfits.PrimaryHDU(data = hdu_data_i, 
                                               header = primaryHeader)
                    # Start by updating PRIMARY header keywords...
                    prihdu.header.set('EXTEND', 
                                          pyfits.FALSE, after = 'NAXIS')
                    
                    # #############################################
                    # AR,DEC (WCS !!) need to be re-calculated !!!
                    # #############################################
                    
                    # Create the new WCS
                    try:
                        orig_ar = float(primaryHeader['RA'])
                        orig_dec = float(primaryHeader['DEC'])
                    except ValueError:
                        # No RA,DEC values in the header, then can't re-compute 
                        # ra,dec coordinates nor update the wcs header.
                        # Then, we DO NOT create a new WCS header !!!
                        pass
                    else:
                        # Due to PANICv0 header hasn't a proper WCS header, we 
                        # built a basic one in order to computer the new RA and 
                        # DEC coordinates.

                        # ----
                        # Actually, WCS object is not USED !!!
                        # Instead, WCS header is built from primaryHeader. 

                        chip_gap = 167 
                        new_wcs = wcs.WCS(primaryHeader)
                        new_wcs.wcs.crpix = [primaryHeader['NAXIS1']/2 + chip_gap/2.0, 
                                             primaryHeader['NAXIS2']/2 + chip_gap/2.0]  
                        new_wcs.wcs.crval =  [primaryHeader['RA'], 
                                              primaryHeader['DEC'] ]
                        new_wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
                        new_wcs.wcs.cunit = ['deg', 'deg']
                        
                        try:
                            pix_scale = primaryHeader['PIXSCALE']
                        except Exception,e:
                            raise Exception("Cannot find PIXSCALE keyword")
                            #pix_scale = 0.45 # arcsec per pixel

                        # We suppose no rotation, and North is up and East at left.
                        # CD matrix (the CDi_j elements) encode the sky position angle,
                        # the pixel scale, and a possible flipping.
                        # CD1_1 is <0 because East is supposed at Left = flipX
                        # CD2_2 is >0 because North is supposed at Up
                        # In addition, it must be noted that:
                        # CD1_1 = cos(r), CD1_2 = sin(r), CD2_1 = -sin(r), CD2_2 = cos(r)
                        # r = clockwise rotation_angle 
                        new_wcs.wcs.cd = [[-pix_scale/3600.0, 0], 
                                          [0, pix_scale/3600.0]]
                        new_pix_center = new_wcs.wcs_pix2world([pix_centers[i*2+j]], 1)
                        
                        #print "Pix Centers = ",pix_centers[i*2+j]
                        #print "New Pix Centers = ", new_pix_center
                        
                        #prihdu.header.set('RA', new_pix_center[0][0])
                        #prihdu.header.set('DEC', new_pix_center[0][1])
                        prihdu.header.set('RA', primaryHeader['RA'])
                        prihdu.header.set('DEC', primaryHeader['DEC'])
                        
                        # ----

                        # Now update the new-wcs for the new subframe
                        prihdu.header.set('CRPIX1', new_crpix_center[i*2+j][0])
                        prihdu.header.set('CRPIX2', new_crpix_center[i*2+j][1])
                        prihdu.header.set('CRVAL1', primaryHeader['RA'])
                        prihdu.header.set('CRVAL2', primaryHeader['DEC'])
                        prihdu.header.set('CTYPE1' , 'RA---TAN') 
                        prihdu.header.set('CTYPE2' , 'DEC--TAN')
                        prihdu.header.set('CUNIT1', 'deg')
                        prihdu.header.set('CUNIT2', 'deg')
                        prihdu.header.set('CD1_1', -pix_scale/3600.0, 
                                             "Axis rotation & scaling matrix")
                        prihdu.header.set('CD1_2', 0, 
                                             "Axis rotation & scaling matrix")
                        prihdu.header.set('CD2_1', 0,  
                                             "Axis rotation & scaling matrix")
                        prihdu.header.set('CD2_2', pix_scale/3600.0, 
                                             "Axis rotation & scaling matrix")
                        prihdu.header.add_history("[MEF.splitGEIRSToSimple] File created from %s"%file)
                        
                    
                    if 'DETXYFLIP' in primaryHeader:
                        # When DETXYFLIP=0, no flip; default
                        # When DETXYFLIP=1, Flip Horizontally (interchanges the left and the right) 
                        # When DETXYFLIP=2, Flip Vertically (turns the image upside down) 
                        detflipxy = primaryHeader['DETXYFLIP']
                        log.debug("Found DETXYFLIP in header =%s"%detflipxy)
                    else: 
                        detflipxy = 0

                    # Force BITPIX=-32
                    prihdu.scale('float32')
                    if 'DETROT90' in primaryHeader:
                        log.debug("Found DETROT90 in header")
                        # When DETROT90=0, no rotation;
                        # When DETROT90=1, 90deg clockwise rotation
                        # When DETROT90=2, 180deg clockwise rotation (default
                        # for PANIC)
                        if primaryHeader['DETROT90']==0:
                            if (i*2+j)==0: 
                                if detflipxy==0: det_id = 2
                                elif detflipxy==1: det_id = 3
                                elif detflipxy==2: det_id = 1
                                elif detflipxy==3: det_id = 4
                            elif (i*2+j)==1: 
                                if detflipxy==0: det_id = 1
                                elif detflipxy==1: det_id = 4
                                elif detflipxy==2: det_id = 2
                                elif detflipxy==3: det_id = 3
                            elif (i*2+j)==2: 
                                if detflipxy==0: det_id = 3
                                elif detflipxy==1: det_id = 2
                                elif detflipxy==2: det_id = 4
                                elif detflipxy==3: det_id = 1
                            elif (i*2+j)==3:
                                if detflipxy==0: det_id = 4
                                elif detflipxy==1: det_id = 1 
                                elif detflipxy==2: det_id = 3
                                elif detflipxy==3: det_id = 2
                        elif primaryHeader['DETROT90']==1:
                            if (i*2+j)==0:
                                if detflipxy==0: det_id = 3
                                elif detflipxy==1: det_id = 4 
                                elif detflipxy==2: det_id = 2
                                elif detflipxy==3: det_id = 1
                            elif (i*2+j)==1:
                                if detflipxy==0: det_id = 2
                                elif detflipxy==1: det_id = 1
                                elif detflipxy==2: det_id = 3
                                elif detflipxy==3: det_id = 4                                
                            elif (i*2+j)==2: 
                                if detflipxy==0: det_id = 4
                                elif detflipxy==1: det_id = 3
                                elif detflipxy==2: det_id = 1
                                elif detflipxy==3: det_id = 2                                
                            elif (i*2+j)==3:
                                if detflipxy==0: det_id = 1
                                elif detflipxy==1: det_id = 2
                                elif detflipxy==2: det_id = 4
                                elif detflipxy==3: det_id = 3                               
                        elif primaryHeader['DETROT90']==2:
                            if (i*2+j)==0:
                                if detflipxy==0: det_id = 4
                                elif detflipxy==1: det_id = 1 
                                elif detflipxy==2: det_id = 3
                                elif detflipxy==3: det_id = 2                               
                            elif (i*2+j)==1:
                                if detflipxy==0: det_id = 3
                                elif detflipxy==1: det_id = 2 
                                elif detflipxy==2: det_id = 4
                                elif detflipxy==3: det_id = 1                               
                            elif (i*2+j)==2:
                                if detflipxy==0: det_id = 1
                                elif detflipxy==1: det_id = 4 
                                elif detflipxy==2: det_id = 2
                                elif detflipxy==3: det_id = 3                              
                            elif (i*2+j)==3: 
                                if detflipxy==0: det_id = 2
                                elif detflipxy==1: det_id = 3 
                                elif detflipxy==2: det_id = 1
                                elif detflipxy==3: det_id = 4
                        elif primaryHeader['DETROT90']==3:
                            if (i*2+j)==0:
                                if detflipxy==0: det_id = 1
                                elif detflipxy==1: det_id = 2 
                                elif detflipxy==2: det_id = 4
                                elif detflipxy==3: det_id = 3
                            elif (i*2+j)==1:
                                if detflipxy==0: det_id = 4
                                elif detflipxy==1: det_id = 2
                                elif detflipxy==2: det_id = 1
                                elif detflipxy==3: det_id = 3
                            elif (i*2+j)==2:
                                if detflipxy==0: det_id = 2
                                elif detflipxy==1: det_id = 1
                                elif detflipxy==2: det_id = 3
                                elif detflipxy==3: det_id = 4
                            elif (i*2+j)==3:
                                if detflipxy==0: det_id = 3
                                elif detflipxy==1: det_id = 4
                                elif detflipxy==2: det_id = 2
                                elif detflipxy==3: det_id = 1
                    else:
                        # Then, we suppose DETROT90=2, DETXYFLIP=0, and default for PANIC !
                        log.warning("No DETROT90 found, supposed DETROT90=2 and DETXYFLIP=0")
                        if (i*2+j)==0: det_id = 4
                        elif (i*2+j)==1: det_id = 3
                        elif (i*2+j)==2: det_id = 1
                        elif (i*2+j)==3: det_id = 2

                            
                    prihdu.header.set('DET_ID', 'SG%s'%det_id, 
                                          "PANIC Detector id SGi [i=1..4]")
                    
                    out_hdulist.append(prihdu)    
                    out_hdulist.verify('ignore')
                
                    new_filename = file.replace(".fits", 
                                                out_filename_suffix % (i*2+j+1)) # number from 1 to 4
                    if out_dir != None: 
                        new_filename = new_filename.replace(os.path.dirname(new_filename), 
                                                            out_dir) 
                    out_filenames.append (new_filename)
     
                    # Now, write the new MEF file
                    out_hdulist.writeto (new_filename, 
                                         output_verify = 'ignore', clobber=True)
                    out_hdulist.close (output_verify = 'ignore')
                    del out_hdulist
                    log.info ("FITS file %s created" % (out_filenames[n]))
                    n += 1
        
        log.debug("End of [splitGEIRSToSimple]")
        return n_ext,out_filenames
                                  
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
                  help = "convert a GEIRS-v0 file (with 1 extension) to a \
                  MEF FITS file with 4 extensions", default = False)
                  
    
    (options, args) = parser.parse_args()
    
    if options.file:
        filelist = [options.file]
    elif options.input_file_list:
        filelist = [line.replace( "\n", "") 
                    for line in fileinput.input(options.input_file_list)]
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
        
              
        
        
