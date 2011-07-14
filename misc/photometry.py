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
import math
import pylab

            
import catalog_query

# Logging
from misc.paLog import log
import misc.utils
import astromatic
import datahandler

class CmdException(Exception): pass


##################################
### Some useful functions ########
##################################       
#from scipy.stats import norm, median
#from scipy.stats.stats import nanmedian,_nanmedian
def MAD(a, c=0.6745, axis=0):
    """
    Median Absolute Deviation along given axis of an array:

    median(abs(a - median(a))) / c

    c = 0.6745 is the constant to convert from MAD to std; it is used by
    default

    """

    good = (a==a)
    a = numpy.asarray(a, numpy.float64)
    if a.ndim == 1:
        d = numpy.median(a[good])
        m = numpy.median(numpy.fabs(a[good] - d) / c)
    else:
        d = numpy.median(a[good], axis=axis)
        # I don't want the array to change so I have to copy it?
        if axis > 0:
            aswp = swapaxes(a[good],0,axis)
        else:
            aswp = a[good]
        m = numpy.median(numpy.fabs(aswp - d) / c, axis=0)

    return m

def nanmedian(arr):
    """
    Returns median ignoring NAN
    """
    return numpy.median(arr[arr==arr])


######### End of useful functions #################

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

def generate_phot_comp_plot ( input_catalog, filter, expt = 1.0 , 
                              out_filename=None, out_format='pdf'):
    """
    @summary: generate a photometry comparison plot, comparing instrumental magnitude
              versus 2MASS photometry
    
    @param catalog : VOTABLE catalog  having the photometric values (instrumental and 2MASS);
    @param filter : NIR wavelength filter (J,H,K,Z)
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
    #command_line = STILTSwrapper._stilts_pathname + " plot2d " + \
    #                " in=" + input + \
    #                " subsetNS='j_k<=1.0 & k_snr>10'  lineNS=LinearRegression" + \
    #                " xdata=k_m ydata=Inst_Mag xlabel=\"2MASS K_m / mag\" " + \
    #                " ylabel=\"Instrumental_Mag K / mag\" " 
                     
    command_line = STILTSwrapper._stilts_pathname + " plot2d " + \
                    " in=" + input + \
                    " subsetNS='h_snr>10 && j_snr>10 && FLAGS==0'  lineNS=LinearRegression" + \
                    " ydata=%s xdata=MAG_AUTO ylabel=\"2MASS %s / mag\" "%(filter,filter) + \
                    " xlabel=\"Instrumental_Mag %s/ mag\" "%filter
                                     
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
    @summary: Compute and Plot the linear regression of two columns of the input vo_catalog
    @param column_x: column number for X values of the regression
    @param column_y: column number for Y values of the regression ( filter )
    
    @return: tuple with linear fit parameters and a Plot showing the fit
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
    
    # ToBeDone
    ## Filter data by FLAGS=0, FLUX_AUTO>0, ...
    ## SNR = FLUX_AUTO / FLUXERR_AUTO
    min_snr = 10
    table_new = table.where( (table.FLAGS==0) & (table.FLUX_BEST > 0) &
                             (table.j_snr>min_snr) & (table.h_snr>min_snr) &
                             (table.k_snr>min_snr) & (table.j_k<1.0) &
                             (table.FLUX_AUTO/table.FLUXERR_AUTO>min_snr))
    
    # If there aren't enough 2MASS objects, don't use color cut
    if len(table_new)<25:
        #table_new = table.where( (table.FLAGS==0) & (table.FLUX_BEST > 0) &
        #                         (table.j_snr>min_snr))
        table_new = table.where( (table.FLAGS==0) & (table.FLUX_BEST > 0) &
                             (table.j_snr>min_snr) & (table.h_snr>min_snr) &
                             (table.k_snr>min_snr) & 
                             (table.FLUX_AUTO/table.FLUXERR_AUTO>min_snr))
        
    #X = -2.5 * numpy.log10(table_new['FLUX_AUTO']/1.0)
    filter = column_y
    X = table_new[column_x]
    Y = table_new[column_y]
   
    #remove the NaN values 
    validdata_X = ~numpy.isnan(X)
    validdata_Y = ~numpy.isnan(Y)
    validdataBoth = validdata_X & validdata_Y
    n_X = X[validdataBoth] #- (0.5*0.05) # a row extinction correction
    n_Y = Y[validdataBoth] 

    # Compute the linear fit
    res = numpy.polyfit(n_X, n_Y, 1, None, True)
    a = res[0][1] # intercept == Zero Point
    b = res[0][0] # slope
    r = res[3][1] # regression coeff ??? not exactly
    
    print "Coeffs =", res
    
    # Plot the results
    pol = numpy.poly1d(res[0])
    plt.plot(n_X, n_Y, '.', n_X, pol(n_X), '-')
    plt.title("Filter %s  -- Poly fit: %f X + %f  r=%f ZP=%f" %(filter, b,a,r,a))
    plt.xlabel("Inst_Mag (MAG_AUTO) / mag")
    plt.ylabel("2MASS Mag  / mag")
    plt.savefig("/tmp/linear_fit.pdf")
    plt.show()
    log.debug("Zero Point (from polyfit) =%f", a)
    
    
    # Compute the ZP as the median  of all per-star ZP=(Mag_2mass-Mag_Inst)
    zps = n_Y - n_X 
    zp = numpy.median(zps)
    #zp = a
    log.debug("Initial ZP(median) = %f"%zp)
    zp_sigma = numpy.std(zps)
    print "ZPS=",zps
    # do a kind of sigma-clipping
    zp_c = numpy.median(zps[numpy.where(numpy.abs(zps-zp)<zp_sigma*2)])
    log.debug("ZP_sigma=%f"%zp_sigma)
    log.debug("Clipped ZP = %f"%zp_c)
    zp = zp_c
    #zp = a
    
    # Now, compute the histogram of errors
    m_err_for_radial_systematic = n_Y - (n_X + zp)
    n_X = n_X[numpy.where(numpy.abs(zps-zp)<zp_sigma*2)]
    n_Y = n_Y[numpy.where(numpy.abs(zps-zp)<zp_sigma*2)]
    log.debug("Number of points = %d"%len(n_X))
    #m_err = n_Y - (n_X*b + zp ) 
    m_err = n_Y - (n_X + zp)
    rms = numpy.sqrt( numpy.mean( (m_err)**2 ) )
    MAD = numpy.median( numpy.abs(m_err-numpy.median(m_err)))
    std = numpy.std(m_err)
    #std2 = numpy.std(m_err[numpy.where(numpy.abs(m_err)<std*2)])
    

    log.debug("ZP = %f"%zp)
    log.debug("MAD(m_err) = %f"%MAD)
    log.debug("MEAN(m_err) = %f"%numpy.mean(m_err))
    log.debug("MEDIAN(m_err) = %f"%numpy.median(m_err))
    log.debug("STD(m_err) = %f"%std)
    log.debug("RMS(m_err) = %f"%rms)
    #log.debug("STD2 = %f"%std2)
    
    #my_mag = n_X*b + zp
    my_mag = n_X + zp
    #plt.plot( my_mag[numpy.where(m_err<std*2)], m_err[numpy.where(m_err<std*2)], '.')
    plt.plot( my_mag, m_err, '.')
    plt.xlabel("Inst_Mag")
    plt.ylabel("2MASS_Mag-Inst_Mag")
    plt.title("(1) Calibration with 2MASS - STD = %f"%std)
    #plt.plot((b * n_X + a), m_err, '.')
    plt.grid(color='r', linestyle='-', linewidth=1)
    plt.savefig("/tmp/phot_errs.pdf")
    plt.show()
    
    # Plot radial distance VS m_err
    radial_distance = numpy.sqrt((table_new['X_IMAGE']-1024)**2 
                                 + (table_new['Y_IMAGE']-1024)**2)*0.45 #arcsecs
    
    plt.plot( radial_distance, m_err_for_radial_systematic, '.')
    plt.xlabel("Radial distance ('arcsec')")
    plt.ylabel("2MASS_Mag-Inst_Mag")
    plt.title("(1) Spatial systematics - STD = %f"%std)
    plt.savefig("/tmp/espatial_systematics_errs.pdf")
    plt.show()
    
    
    # Second ZP
    log.debug("Second iteration")
    temp = n_Y - n_X 
    print "LEN1=",len(temp)
    print "LEN2=",len(temp[numpy.where(numpy.abs(m_err)<std*2)])
    zp2 = numpy.median( temp[numpy.where(numpy.abs(m_err)<std*2)])
    m_err2 = n_Y - (n_X + zp2) 
    rms2 = numpy.sqrt( numpy.mean( (m_err2)**2 ) )
    std2 = numpy.std(m_err2)
    MAD2 = numpy.median( numpy.abs(m_err2-numpy.median(m_err2)))
    MAD2b = numpy.sqrt( numpy.median( (m_err2-numpy.median(m_err2))**2 ) )
    #std3 = numpy.std(m_err[numpy.where(numpy.abs(m_err2)<std2*2)])

    log.debug("ZP2 = %f"%zp2)
    log.debug("MAD2(m_err2) = %f"%MAD2)
    log.debug("MAD2b(m_err2) = %f"%MAD2b)
    log.debug("MEAN2(m_err2) = %f"%numpy.mean(m_err2))
    log.debug("STD(m_err2) = %f"%std2)
    log.debug("RMS(m_err2) = %f"%rms2)

    #log.debug("STD3 = %f"%std3)
    
    
    

    #print m_err
    # Lo normal es que coincida la RMS con la STD, pues la media de m_err en este caso es 0
    my_mag = n_X+zp2
    plt.plot( my_mag[numpy.where(m_err<std*2)], m_err[numpy.where(m_err<std*2)], '.')
    plt.xlabel("Inst_Mag")
    plt.ylabel("2MASS_Mag-Inst_Mag")
    plt.title("(2) Calibration with 2MASS - STD = %f"%std2)
    #plt.plot((b * n_X + a), m_err, '.')
    plt.savefig("/tmp/phot_errs.pdf")
    plt.show()
    
    
    pylab.hist(m_err2, bins=50, normed=0)
    pylab.title("Mag error Histogram - RMS = %f mag STD = %f"%(rms2,std2))
    pylab.xlabel("Mag error")
    pylab.ylabel("Frequency")
    plt.savefig("/tmp/phot_hist.pdf")
    pylab.show()
    

    
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
        filter = my_fits.getFilter()
        log.debug("EXPTIME = %f"%exptime)
        log.debug("RA = %f"%ra)
        log.debug("DEC = %f"%dec)
        log.debug("Filter = %s"%filter)
        if ra<0 or exptime<0:
            log.debug("Found wrong RA,DEC,EXPTIME or FILTER value")
            sys.exit(0)
        
        # Checking Filter    
        if filter=='J':
            two_mass_col_name = 'j_m'
        elif filter=='H':
            two_mass_col_name = 'h_m'
        elif filter=='K' or filter=='K-PRIME' or filter=='KS':
            two_mass_col_name = 'k_m'   
        #elif filter=='OPEN':
        #    two_mass_col_name = 'j_m'
        else:
            log.error("Filter %s not supported" %filter)
            sys.exit(0)
            
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
                                   error=1.0 ) 
        log.debug("XMatch done !")
    except Exception,e:
        log.error("XMatch failed %s", str(e))
        sys.exit(0)
    
    ## 3a- Compute the linear regression (fit) of Inst_Mag VS 2MASS_Mag
    ###### and generate the plot file with the photometric comparison
    ###### 2MASS_Mag = Inst_Mag*b + ZP  
    log.debug("Compute & Plot regression !!!")    
    try:
        compute_regresion (out_xmatch_file, 'MAG_AUTO', two_mass_col_name )
        log.debug("Well DONE !!")
        #sys.exit(0)
    except Exception,e:
        log.error("Sorry, can't some error while computing linear fit or\
        ploting the results: %s", str(e))
        sys.exit(0)
    
    ## 3b- Compute the linear regression (fit) of Inst_Mag VS 2MASS_Mag
    ###### and generate the plot file with the photometric comparison  
    try:
        exptime = 1.0 # SWARP normalize flux to 1 sec
        plot_file = generate_phot_comp_plot ( match_cat, two_mass_col_name, exptime, 
                                              options.output_filename, 
                                              out_format='pdf')
        log.debug("Plot file generated : %s", plot_file) 
    except Exception,e:
        log.error("Sorry, can't some error while computing linear fit or \
        ploting the results: %s", str(e))
        sys.exit(0)
        
        
        
