#!/usr/bin/env python

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
#########################################################################
# PANIC data processing
#########################################################################
# Example of nonlinearity correction for full frame MEF files
# Includes many sanity checks, the actual calculation is very simple
# Incorrectable and saturated pixels are set to NaN. Therefore the ouptut
# data is float.
#
# 1.0 14/07/2014 BD Creation (correct_nonlinearity.py)
#
# 1.1 18/07/2014 JMIM Adaption to PAPI (correctNonLinearity.py)
# 
#

_name = 'correctNonLinearity.py'
_version = '1.1'
################################################################################
# Import necessary modules
import numpy as np
import astropy.io.fits as fits
import dateutil.parser
import sys
import os
import fileinput
from optparse import OptionParser
import multiprocessing


# PAPI modules
import misc.mef
from misc.paLog import log
from misc.version import __version__

# If you want to use the new multiprocessing module in Python 2.6 within a class, 
# you might run into some problems. Here's a trick how to do a work-around. 
def unwrap_self_applyModel(arg, **kwarg):
    #print "ARG=",arg
    return NonLinearityCorrection.applyModel(*arg, **kwarg)

class NonLinearityCorrection(object):
    """
    Class used to correct the Non-linearity of the PANIC detectors based on the 
    algorithm described by Bernhard Dorner at PANIC-DEC-TN-02_0_1.pdf.
    """
    def __init__(self, model, input_files, out_dir='/tmp', 
                suffix='_LC', force=False, coadd_correction=True):
        """
        Init the object.
        
        Parameters
        ----------
                
        model : str 
            FITS filename of the Non-Linearity model, ie., containing polynomial 
            coeffs (4th order) for correction that must has been previously 
            computed. It must be a cube with 4 planes (a plane for each coeff c4
            to c1, c0 is not used in the correction and not stored in the cube), 
            and N extensions (one per detector). Planes definitions:

                plane_0 = coeff_4 
                plane_1 = coeff_3 
                plane_2 = coeff_2 
                plane_3 = coeff_1 
        
        input_files: list 
            A list of FITS files to be corrected (MEF or single FITS).

        out_dir: str
            Directory where new corredted files will be saved

        suffix: str
            Suffix to add to the input filename to generate the output filename.

        force: bool
            If true, no check of input raw header is done (NCOADDS, DETROT90, 
            INSTRUME,...)
            
        coadd_correction: bool
            If true and NCOADDS>1, divide the data by NCOADDS, apply NLC and then 
            multiply by NCOADDS.

        Returns
        -------
        outfitsname: list
            The list of new corrected FITS files created.
        
        """
        self.input_files = input_files
        self.model = model
        self.suffix = suffix
        self.out_dir = out_dir
        self.force = force
        self.coadd_correction = coadd_correction
        
        if not os.access(self.out_dir, os.F_OK):
            try:
                os.mkdir(self.out_dir)
            except Exception,e:
                raise e
        
        if len(self.input_files)<1:
            msg = "Found empty list of input files"
            log.error(msg)
            raise Exception(msg)
        
        if not os.path.exists(self.model):
            msg = "Cannot read non-linearity model file '%s'" % self.model
            log.error(msg)
            raise Exception(msg)

    def checkHeader(self, modelHeader, dataHeader):
        """
        Performs some data checking (readmode, orientation, date-obs, 
        savemode,...).

        Parameters
        ----------
        modelHeader: 
            Header of the NLC model.
        dataHeader:
            Header of an raw data file.

        Returns
        -------
           If a non-compliant header is found, some exception will be raised.

        """


        # First, checks model is a MASTER_LINEARITY
        if modelHeader['PAPITYPE'] != 'MASTER_LINEARITY':
            raise ValueError('Wrong type of nonlinearity correction file')

        # Check input files are non-integrated files (NCOADDS)
        # It is done on applyModel, where can be skipped.
        #if dataHeader['NCOADDS']>1:
        #    log.info("Found a wrong type of source file. Use -F to user ncoadd correction")
        #    raise ValueError('Wrong type of file. Only non-integrated files (NCOADDS=1) allowed.')

        # Check NLC model is used with newer data (USE_AFTER->USE_AFT)
        datadate = dateutil.parser.parse(dataHeader['DATE-OBS'])
        nldate = dateutil.parser.parse(modelHeader['USE_AFT'])
        if datadate < nldate:
            raise ValueError('Nonlinearity calibration too new for input data')           
        
        # Check some other keys related with READOUT configuration 
        keys = ['INSTRUME', 'PREAD', 'PSKIP', 'LSKIP', 'READMODE', 'IDLEMODE', 'IDLETYPE']
        for key in keys:
            if str(dataHeader[key]).lower() != str(modelHeader[key]).lower():
                raise ValueError('Mismatch in header data for keyword \'%s\'' %key)            
        
        # some may not be present in old data
        keys = ['DETROT90', 'DETXYFLI']
        for key in keys:
            if not key in dataHeader:
                print 'Warning: key \'%s\' not in data header' %key
            elif dataHeader[key] != modelHeader[key]:
                raise ValueError('Mismatch in header data for keyword \'%s\'' %key)
        keys = ['B_EXT', 'B_DSUB', 'B_VREST', 'B_VBIAG']
        for key in keys:
            for i in range (1, 5):
                if dataHeader[key + '%i'%i] != modelHeader[key + '%i'%i]:
                    raise ValueError('Mismatch in header data for keyword \'%s%i\'' %(key, i))


    def applyModel(self, data_file):
        """
        Do the Non-linearity correction using the supplied model. In principle,
        it should be applied to all raw images (darks, flats, science, ...).
        
        Parameters
        ----------
        data_file: str
            input data FITS filename to be corrected.

        Returns
        -------
        outfitsname: str
            The list of new corrected files created.
                
        """   
        
        # load raw data file
        hdulist = fits.open(data_file)
        dataheader = hdulist[0].header
        
        # Check if input files are in MEF format or saved as a full-frame with 
        # the 4 detectors 'stitched'.
        to_delete = None
        if len(hdulist)==1:
            # we need to convert to MEF
            log.warning("Mismatch in header data format. Converting to MEF file.")
            hdulist.close()
            
            # Convert single-FITS to MEF
            mef = misc.mef.MEF([data_file])
            mef_suffix = ".mef.fits"
            n_ext, new_mef_files = mef.convertGEIRSToMEF(mef_suffix, self.out_dir)
            if n_ext !=4:
                raise ValueError('Mismatch in header data format. Only MEF files allowed.')
            
            # load new MEF raw data file
            hdulist = fits.open(new_mef_files[0])
            dataheader = hdulist[0].header
            # copy the filename to be deleted after processing (?)
            to_delete = new_mef_files[0]

        # load model
        nlhdulist = fits.open(self.model)
        nlheader = nlhdulist[0].header

        # Check headers
        try:
            if self.force == False:
                self.checkHeader(nlheader, dataheader)
        except Exception,e:
            log.error("Mismatch in header data for input NLC model %s"%str(e))
            raise e

        # Creates output fits HDU
        linhdu = fits.PrimaryHDU()
        linhdu.header = dataheader.copy()
        hdus = []

        # Check which version of MEF we have.
        # Since GEIRS-r731M-18 version, new MEF extension naming:
        #    EXTNAME = 'Qi_j'
        #    DET_ID = 'SGi_j'
        # and the order in the MEF file shall be Q1, Q2, Q3, Q4
        try:
            hdulist['Q1_1'].header
            ext_name = 'Q%i_1'
            ext_order = (1, 2, 3, 4)
        except KeyError:
            ext_name = 'SG%i_1'
            ext_order = (4, 1, 3, 2)
            
        # loop over detectors
        # To avoid the re-arrange of the MEF extensions
        for iSG in ext_order:
            extname = ext_name % iSG
            # check detector sections
            # another way would be to loop until the correct one is found
            datadetsec = hdulist[extname].header['DETSEC']
            nldetsec = nlhdulist['LINMAX%i' %iSG].header['DETSEC']
            if datadetsec != nldetsec:
                raise ValueError('Mismatch of detector sections for SG%i' %iSG)
            # or check SG IDs (as long as they are reliable)
            datadetid = hdulist[extname].header['DET_ID']
            nldetid = nlhdulist['LINMAX%i' %iSG].header['DET_ID']
            if datadetid != nldetid:
                raise ValueError('Mismatch of detector IDs for extension' %extname)

            # Work around to correct data when NCOADDS>1
            if hdulist[0].header['NCOADDS']>1:
                if self.coadd_correction:
                    log.info("NCOADDS>1; Doing ncoadd correction...")
                    n_coadd = hdulist[0].header['NCOADDS']
                else:
                    log.info("Found a wrong type of source file. Use -c to user ncoadd correction")
                    raise ValueError('Cannot apply model, found NCOADDS > 1.')
            else:
                n_coadd = 1
                
            # load file data (and fix coadded images => coadd_correction)
            data = hdulist[extname].data / n_coadd
            nlmaxs = nlhdulist['LINMAX%i' %iSG].data
            nlpolys = np.rollaxis(nlhdulist['LINPOLY%i' %iSG].data, 0, 3)

            # calculate linear corrected data
            lindata = self.polyval_map(nlpolys, data)
            
            # mask saturated inputs - to use nan it has to be a float array
            lindata[data > nlmaxs] = np.nan
            # mask where max range is nan
            # (first, take into account the option of cubes as input images) 
            if len(lindata.shape)==3: 
                # we have a 3D image (cube)
                for i in range(lindata.shape[0]):
                    lindata[i, np.isnan(nlmaxs)] = np.nan
            else:
                # we have a single 2D image
                lindata[np.isnan(nlmaxs)] = np.nan

            # Undo the coadd_correction
            lindata = lindata * n_coadd
            
            exthdu = fits.ImageHDU(lindata.astype('float32'), header=hdulist[extname].header.copy())
            # this may rearrange the MEF extensions, otherwise loop over extensions
            hdus.append(exthdu)

        # add some info in the header
        linhdu.header['HISTORY'] = 'Nonlinearity correction applied'
        linhdu.header['HISTORY'] = 'Nonlinearity data: %s' %nlheader['ID']
        linhdu.header['HISTORY'] = '<-- The German team made this on 2014/07/13'
        linhdu.header.set('PAPIVERS', __version__,'PANIC Pipeline version')
        linhdulist = fits.HDUList([linhdu] + hdus)
        
        # Compose output filename
        mfnp = os.path.basename(data_file).partition('.fits')
        # add suffix before .fits extension, or at the end if no such extension present
        outfitsname = self.out_dir + '/' + mfnp[0] + self.suffix + mfnp[1] + mfnp[2]
        outfitsname = os.path.normpath(outfitsname)

        # overwrite the output file if exists
        linhdulist.writeto(outfitsname, clobber=True)
        
        if to_delete:
            os.unlink(to_delete)

        return outfitsname


    def runMultiNLC(self):
        """
        Run a parallel proceesing of NL-correction for the input files taking
        advantege of multi-core CPUs.

        Returns
        -------
        On succes, a list with the filenames of the corrected files.
        """

        # use all CPUs available in the computer
        n_cpus = multiprocessing.cpu_count()
        log.debug("N_CPUS :" + str(n_cpus))
        pool = multiprocessing.Pool(processes=n_cpus)
        
        results = []
        solved = []
        for i_file in self.input_files:
            red_parameters = [i_file]
            try:
                # Instead of pool.map() that blocks until
                # the result is ready, we use pool.map_async()
                results += [pool.map_async(unwrap_self_applyModel, 
                        zip([self]*len(red_parameters), red_parameters) )]
            except Exception,e:
                log.error("Error processing file: " + i_file)
                log.error(str(e))
                
        for result in results:
            try:
                result.wait()
                # the 0 index is *ONLY* required if map_async is used !!!
                solved.append(result.get()[0])
                log.info("New file created => %s"%solved[-1])
            except Exception,e:
                log.error("Cannot process file \n" + str(e))
                
        

        # Prevents any more tasks from being submitted to the pool. 
        # Once all the tasks have been completed the worker 
        # processes will exit.
        pool.close()

        # Wait for the worker processes to exit. One must call 
        #close() or terminate() before using join().
        pool.join()
        
        log.info("Finished parallel NL-correction")
        
        return solved

    def polyval_map(self, poly, map):
        """
        Evaluate individual polynomials on an array. Looping over each pixel
        is stupid, therefore we loop over the order and calculate the
        polynomial directly.
        Note: The output is a float array!
        
        Input
        -----
        poly : array_like
               Polynomial coefficients without constant offset. The order
               must be along the last axis.
        map : array_like
              Data array, that can be a cube of 2D images.
              
        Returns
        -------
        polymap : array_like
                  Result of evaluation, same shape as map, dtype float
        """

        order = poly.shape[-1]
        polymap = map * 0.
        for io in range(order):
            polymap += poly[Ellipsis, -io-1] * map**(io+1)
        return polymap

################################################################################
# main
if __name__ == "__main__":
    
    
    usage = "usage: %prog [options]"
    desc= """Performs the non-linearity correction of the PANIC raw data files
using the proper NL-Model (FITS file).
"""
    parser = OptionParser(usage, description=desc)

    # Basic inputs
    parser.add_option("-m", "--model",
                  action="store", dest="model",
                  help="FITS MEF-cube file of polinomial coeffs (c4, c3, c2, c1) of the NL model.")
    
    parser.add_option("-i", "--input_file",
                  action="store", dest="input_file",
                  help="FITS file file to be corrected.")
    
    parser.add_option("-s", "--source",
                  action="store", dest="source_file_list",
                  help="Source file list of FITS files to be corrected.")
    
    parser.add_option("-o", "--out_dir", type="str", dest="out_dir",
                  action="store", default="/tmp",
                  help="filename of out data file (default=%default)")
    
    parser.add_option("-S", "--suffix", type="str",
                  action="store", dest="suffix", default="_NLC", 
                  help="Suffix to use for new corrected files (default=%default)")

    parser.add_option("-f", "--force",
                  action="store_true", dest="force", default=False, 
                  help="Force Non-linearity correction with no check of header"
                  "values (NCOADD, DATE-OBS, DETROT90, ...")
    
    parser.add_option("-c", "--coadd_correction",
                  action="store_true", dest="coadd_correction", default=True, 
                  help="Force NCOADDS correction and apply NLC")
    
    (options, args) = parser.parse_args()
    
    if len(sys.argv[1:])<1:
       parser.print_help()
       sys.exit(0)

    # Check required parameters
    if ((not options.source_file_list and not options.input_file) or not options.out_dir 
        or not options.model  or len(args)!=0): # args is the leftover positional arguments after all options have been processed
        parser.print_help()
        parser.error("incorrect number of arguments " )
    
    if options.input_file and os.path.isfile(options.input_file): 
        filelist = [options.input_file]
    elif options.source_file_list and os.path.isfile(options.source_file_list):
        # Read the source file list     
        filelist = [line.replace( "\n", "") for line in fileinput.input(options.source_file_list)]
    else:
        parser.print_help()
        parser.error("incorrect number of arguments " )

    NLC = NonLinearityCorrection(options.model, filelist, options.out_dir, 
                                   options.suffix, options.force,
                                   options.coadd_correction)

    try:
        corr = NLC.runMultiNLC()
    except Exception,e:
        log.error("Error running NLC: %s"%str(e))
        raise e

    print "\nCorrected files: ",corr
    
    """
    # Non parallel processing
    for i_file in filelist:
        try:
            NLC.applyModel(i_file)
        except Exception,e:
            log.error("Error applying NLC model to file '%s': %s"%(i_file, str(e)))
    """ 
    
