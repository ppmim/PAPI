#!/usr/bin/env python
"""
Module to do some photometry functionalities 
"""

################################################################################
#
#
# PAPI (PANIC PIpeline)
#
# photometry.py
#
# Tools used:
#
#    STILTS command line application  (http://www.star.bris.ac.uk/~mbt/stilts)
#
# Created    : 25/05/2011    jmiguel@iaa.es -
# Last update: 
# TODO
#       
################################################################################

################################################################################

# Import necessary modules
from optparse import OptionParser
import sys
import os
import atpy 
import matplotlib.pyplot as plt
import numpy

            
import catalog_query

# Logging
from misc.paLog import log
import misc.utils
import astromatic
import datahandler

class CmdException(Exception): pass

def catalog_xmatch ( cat1, cat2, out_filename, out_format='votable', error=2.0 ):
    """
    @summary: do catalogs cross-match
    
    @param cat1,cat2: catalogs for cross-matching
    @param err: max. error for finding objects within (arcseconds)
    @param out_filename: filename where results will be saved;if absent, 
            the location will be a tempfile with a generated name
    @param out_format: format of the output generated; current options available 
            are:
        - VO Table (XML) (votable) (default)
        - SVC (Software handshaking structure) message (svc)
        - ASCII table (ascii)
    @return: filename where results where saved (VOTABLE, ASCII_TABLE, ...)
    """
    
    in1 = cat1
    in2 = cat2
    out = out_filename
    s_error = str(error)
    ar1 = ''
    dec1 = ''
    ar2 = ''
    dec2 = ''
    
    # del old instances
    if os.path.exists(out): os.remove(out)
    
    command_line = STILTSwrapper._stilts_pathname + " tskymatch2 " + \
         " in1=" + in1 + " in2=" + in2 + " out=" + out + \
         " error=" + s_error
          
    rcode = misc.utils.runCmd(command_line)
    
    if rcode==0 or not os.path.exists(out):
        log.error("Some error while running command: %s", command_line)
        raise CmdException("XMatch failed")
    else:
        return out   

def generate_phot_comp_plot ( input_catalog, expt = 1.0 , 
                              out_filename=None, out_format='pdf'):
    """
    @summary: generate a photometry comparison plot, comparing instrumental magnitude
              versus 2MASS photometry
    
    @param catalog : VOTABLE catalog  having the photometric values (instrumental and 2MASS);
    @param expt : exposure time of original input image; needed to 
                compute the Instrumental Magnitude (Inst_Mag)  
    @param out_filename: filename where results will be saved;if absent, 
            the location will be a tempfile with a generated name
    @param out_format: format of the output generated; current options available are:
        - 'pdf' (default)
        - 'jpg' 
        - 'gif'
        - for more formats, see  http://www.star.bris.ac.uk/~mbt/stilts/sun256/plot2d-usage.html
        
    @return: filename where results where saved (VOTABLE, ASCII_TABLE, ...)
    
    @note: > stilts tpipe  ifmt=votable cmd='addcol MyMag_K "-2.5*log10(FLUX_BEST/46.0)"' omode=out in=/home/panic/SOFTWARE/STILTS/match_double.vot out=match_D_b.vot
           > stilts plot2d in=match_S_b.vot subsetNS='j_k<=1.0 & k_snr>10' lineNS=LinearRegression xdata=k_m ydata=MyMag_k xlabel="2MASS k_m / mag" ylabel="PAPI k_m / mag"


    """
    
    log.debug("entering in <generate_phot_comp_plot>")
    ## 1 - First, we add a column with the Instrumental magnitude, computed as:
    ##    inst_mag = -2.5 log10(FLUX_BEST/TEXP)
    ## where FLUX_BEST is obtained from SExtractor output catalog
    
    input = input_catalog
    output_1 ="/tmp/output_1.xml"
    
    command_line = STILTSwrapper._stilts_pathname + " tpipe " + \
                    " ifmt=votable" + \
                    " cmd='addcol Inst_Mag \"-2.5*log10(FLUX_BEST/%f)\"'" % expt + \
                    " omode=out" + " in=" + input + " out=" + output_1

    rcode = misc.utils.runCmd(command_line)
    
    if rcode==0 or not os.path.exists(output_1):
        log.error("Some error while running command: %s", command_line)
        raise CmdException("STILTS command failed !")
    
        
    ## 2 - Secondly, generate the plot for photometric comparison
    input = output_1
    command_line = STILTSwrapper._stilts_pathname + " plot2d " + \
                    " in=" + input + \
                    " subsetNS='j_k<=1.0 & k_snr>10'  lineNS=LinearRegression" + \
                    " xdata=k_m ydata=Inst_Mag xlabel=\"2MASS K_m / mag\" " + \
                    " ylabel=\"Instrumental_Mag K / mag\" " 
                     
    if out_filename :
        command_line += " ofmt=" + out_format + " out=" + out_filename
                    
    rcode = misc.utils.runCmd(command_line)
    
    if rcode==0:
        log.error("Some error while running command: %s", command_line)
        raise CmdException("STILTS command failed !")
    
    else:
        if out_filename: return out_filename
        else: return 'stdout'


def compute_regresion ( vo_catalog, column_x, column_y ):
    """
    @summary: Compute the linear regression of two columns of the input vo_catalog
    @param column_x: column number for X values of the regression
    @param column_y: column number for Y values of the regression
    
    @return: tuple with linear fit parameters 
             a - intercept 
             b - slope of the linear fit
             r - estimated error
             
             None if happen some error  
    """
    
    try:
        table = atpy.Table(vo_catalog)
    except Exception,e:
        log.error("Canno't read the input table")
        return None
    
    X = table[column_x]
    Y = table[column_y]
   
    #remove the NaN values 
    validdata_X = ~numpy.isnan(X)
    validdata_Y = ~numpy.isnan(Y)
    validdataBoth = validdata_X & validdata_Y
    n_X = X[validdataBoth]
    n_Y = Y[validdataBoth] 

    # Compute the polyfit
    res = numpy.polyfit(n_X, n_Y, 1, None, True)
    a = res[0][1] # intercept
    b = res[0][0] # slope
    r = res[3][1] # regression coeff
    
    # Plot the results
    pol = numpy.poly1d(res[0])
    plt.plot(n_X, n_Y, '.', n_X, pol(n_X), '-')
    plt.title("Poly fit: %f X + %f  r=%f" %(b,a,r))
    plt.xlabel("2MASS Mag K / mag")
    plt.ylabel("Inst_Mag K / mag")
    plt.show()

    
class STILTSwrapper (object):
    """ Make a wrapper to some functionalities of STILTS 
    """
    
    cat_names = {'2MASS':'fp_psc', 
                     'USNOB1': 'usno_b1',
                     'IRAS': 'iraspsc'
                     }
    outfmt = {'votable': 3,
              'ascii': 1
              }
    
    _stilts_pathname = "/home/panic/SOFTWARE/STILTS/stilts"
    
    def __init__(self, *a, **k):
        """ The constructor """
        super (STILTSwrapper, self).__init__ (*a, **k)
        
    
    def runXMatch(self, cat1, cat2, out_filename=None, out_format='votable', error=2.0):    
        """
        @summary: do catalogs cross-match
        
        @param cat1,cat2: catalogs for cross-matching
        @param err: max. error for finding objects within (arcseconds)
        @param out_filename: filename where results will be saved;if absent, 
                the location will be a tempfile with a generated name
        @param out_format: format of the output generated; current options available are:
            - VO Table (XML) (votable) (default)
            - SVC (Software handshaking structure) message (svc)
            - ASCII table (ascii)
        @return: filename where results where saved (VOTABLE, ASCII_TABLE, ...)
        """
        
        in1 = cat1
        in2 = cat2
        out = out_filename
        s_error = str(error)
        ar1 = ''
        dec1 = ''
        ar2 = ''
        dec2 = ''
        
        command_line = STILTSwrapper._stilts_pathname + " tskymatch2 " + \
             " in1=" + in1 + " in2=" + in2 + " out=" + out + \
             " error=" + error
              
        rcode = misc.utils.runCmd(command_line)
        
        if rcode==0:
            log.error("Some error while running command: %s", command_line)
        else:
            return out_filename    
        
        """
        ./stilts tskymatch2 in1=/tmp/alh_single.fits.xml in2=/tmp/prueba.xml out=match.xml error=2
        ./stilts tpipe  ifmt=votable cmd='addcol MyMag_K "-2.5*log10(FLUX_BEST/46.0)"' omode=out in=/home/panic/SOFTWARE/STILTS/match_double.vot out=match_D_b.vot
        ./stilts plot2d in=match_S_b.vot subsetNS='j_k<=1.0 & k_snr>10' lineNS=LinearRegression xdata=k_m ydata=MyMag_k xlabel="2MASS k_m / mag" ylabel="PAPI k_m / mag"
        """

################################################################################
# main
################################################################################
if __name__ == "__main__":

    log.debug( 'Testing Photometric calibration comparison with 2MASS')
    
        # Get and check command-line options
        
    usage = "usage: %prog [options] arg1 arg2 ..."
    parser = OptionParser(usage)
    
    parser.add_option("-i", "--input_image",
                  action="store", dest="input_image", help="input image to calibrate\
                  to do photometric comparison with")
                  
    parser.add_option("-c", "--base_catalog (2MASS, USNO-B)",
                  action="store", dest="base_catalog",
                  help="Name of base catalog to compare with (2MASS, USNO-B) -- not used !!!")
    
    
    parser.add_option("-o", "--output",
                  action="store", dest="output_filename", help="output plot filename")
    
    
    parser.add_option("-v", "--verbose",
                  action="store_true", dest="verbose", default=True,
                  help="verbose mode [default]")
    
                                
    (options, args) = parser.parse_args()
    
    
    if not options.input_image or len(args)!=0: 
    # args is the leftover positional arguments after all options have been processed
        parser.print_help()
        parser.error("wrong number of arguments " )
    if not options.output_filename:
        options.output_filename=None
    

    if not os.path.exists(options.input_image):
        log.error ("Input image %s does not exist", options.input_image)
        sys.exit(0)
        
    ## 0 - Generate image catalog (VOTable) -> SExtractor
    log.debug("*** Creating SExtractor VOTable catalog ....")
    #tmp_fd, tmp_name = tempfile.mkstemp(suffix='.xml', dir=os.getcwd())
    #os.close(tmp_fd)
    image_catalog = os.path.splitext(options.input_image)[0]  + ".xml"
    sex = astromatic.SExtractor()
    sex.ext_config['CHECKIMAGE_TYPE'] = "NONE"
    sex.config['CATALOG_TYPE'] = "ASCII_VOTABLE"
    sex.config['CATALOG_NAME'] = image_catalog
    sex.config['DETECT_THRESH'] = 1.5
    sex.config['DETECT_MINAREA'] = 5
    try:
        sex.run(options.input_image, updateconfig=True, clean=False)
    except Exception,e:
        log.errro("Canno't create SExtractor catalog : %s", str(e)) 
        sys.exit(0)
    
    ## 0.1 - Read the RA,Dec and TEXP values from the input image
    try:
        my_fits = datahandler.ClFits(options.input_image)
        exptime  = my_fits.expTime()
        ra = my_fits.ra
        dec = my_fits.dec
    except Exception,e:
        log.error("Cannot read properly FITS file : %s:", str(e))
        sys.exit(0)
            
    ## 1 - Generate region of base catalog
    icat = catalog_query.ICatalog ()
    out_base_catalog = os.getcwd() + "/catalog_region.xml"
    sr = 500 # arcsec
    try:
        res_file = icat.queryCatalog(ra, dec, sr, 
                                     catalog_query.ICatalog.cat_names['2MASS'], 
                                     out_base_catalog, 'votable')[0]
        log.debug("Output file generated : %s", res_file) 
    except Exception,e:
        log.error("Sorry, cann't solve the query to ICatalog: %s", str(e))
        sys.exit(0)

    ## 2- XMatch the catalogs (image_catalog VS just base_catalog generated)          
    out_xmatch_file = os.getcwd() + "/xmatch.xml"
    try:
        match_cat = catalog_xmatch( image_catalog, res_file, 
                                   out_xmatch_file, out_format='votable', 
                                   error=2.0 ) 
        log.debug("XMatch done !")
    except Exception,e:
        log.error("XMatch failed %s", str(e))
        sys.exit(0)
    
    log.debug("My Own regression !!!")    
    compute_regresion (out_xmatch_file, 'k_m', 'Inst_Mag' )
    
    log.debug("Well DONE !!")    
    
    ## 3- Generate the plot file with the photometric comparison     
    try:
        plot_file = generate_phot_comp_plot ( match_cat, exptime, 
                                              options.output_filename, 
                                              out_format='pdf')
        log.debug("Plot file generated : %s", plot_file) 
    except Exception,e:
        log.error("Sorry, can't generate plot file: %s", str(e))
        sys.exit(0)
        
        
        
