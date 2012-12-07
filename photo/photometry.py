#!/usr/bin/env python

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
#   - Add Extintion coefficient computation using Airmass
#   
################################################################################

# Import necessary modules
from optparse import OptionParser
import sys
import os

import atpy 
import numpy

#import matplotlib
#matplotlib.use('TkAgg')

import matplotlib.pyplot as plt
#import matplotlib.mlab as mlab
import pylab
            
import catalog_query
import photo.coords as coords


# Logging
from misc.paLog import log
import misc.utils
import astromatic
import datahandler
import misc.robust as robust

class CmdException(Exception): 
    pass


def catalog_xmatch( cat1, cat2, out_filename, out_format='votable', error=2.0 ):
    """
    Takes two input catalogues (VOTables) and performs a cross match to 
    find objects within 'error' arcseconds of each other. 
    The result is a new VOTable (default) containing only rows where a match 
    was found. 
    
    Notes
    -----
    It runs without explicit specification of the sky position columns 
    in either table (OBS_RA,OBS_DEC). It will work only if those columns are 
    identified with appropriate UCDs, for instance pos.eq.ra;meta.main and 
    pos.eq.dec:meta.main. If no suitable UCDs are in place this invocation will
    fail with an error. 
    
    Parameters
    ----------
    
    cat1, cat2: str
        catalogs for cross-matching
    err: float
        max. error for finding objects within (arcseconds)
    out_filename: str
        filename where results will be saved
    out_format: str
        format of the output generated; current options 
            available are:
            - VO Table (XML) (votable) (default)
            - SVC (Software handshaking structure) message (svc)
            - ASCII table (ascii)
    
    Returns
    -------
        Filename where results where saved (VOTABLE, ASCII_TABLE, ...)
    """
    
    in1 = cat1
    in2 = cat2
    out = out_filename
    s_error = str(error)
    
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

def catalog_xmatch2( cat1, cat2, filter_column, error=2.0, min_snr=10):
    """
    Takes two input catalogues (VOTables) and performs a cross match to 
    find objects within 'error' arcseconds of each other. 
    The result is a new VOTable (default) containing only rows where a match 
    was found. 
    
    Notes
    -----
    It runs without explicit specification of the sky position columns 
    in either table (OBS_RA,OBS_DEC). It will work only if those columns are 
    identified with appropriate UCDs, for instance pos.eq.ra;meta.main and 
    pos.eq.dec:meta.main. If no suitable UCDs are in place this invocation will
    fail with an error. 
    
    Parameters
    ----------
    cat1: str
        filename of sextractor image vo-catalog 
    cat2: str
        filename of reference (2MASS, USNO, ...) vo-catalog for cross-matching
    err: float
        max. error for finding objects within tolerance (arcseconds) 
    min_snr: float
        min. SNR value for matched stars
    
    Returns
    -------
        Filename where results where saved (VOTABLE, ASCII_TABLE, ...)
    """
    
    
    # Read vo-catalogs
    #Sextractor catalog
    try:
        table1 = atpy.Table(cat1, type='vo', tid=0)
    except Exception:
        log.error("Canno't read the input catalog %s"%cat1)
        return None
    
    # References (2MASS) catalog
    try:
        table2 = atpy.Table(cat2, type='vo', tid=0)
    except Exception:
        log.error("Canno't read the input catalog %s"%cat2)
        return None
    
    try:
        ind1, ind2 = coords.indmatch(table1['X_WORLD'], table1['Y_WORLD'], 
                                     table2['ra'], table2['dec'],
                                     error)
    except Exception, e:
        log.error("Erron in xmatch-ing tables :%s"%str(e))
        raise e

    table1_xm = table1.where(ind1)
    table2_xm = table2.where(ind2)
    
    
    # Because Atpy does not support table merging, a new xmatch is done with
    # selected rows.
    # Firstly, selection
    table2_tmp = table2_xm.where( (table2_xm.j_snr > min_snr) & 
                                  (table2_xm.h_snr > min_snr) &
                                  (table2_xm.k_snr > min_snr) & 
                                  (table2_xm.j_k < 1.0) )
    
    # if there is no enough stars, then don't use color cut (J-K)<1 
    if len(table2_tmp)<25:
        table2_tmp = table2_xm.where( (table2_xm.j_snr > min_snr) & 
                                      (table2_xm.h_snr > min_snr) &
                                      (table2_xm.k_snr > min_snr))
        
    table1_tmp = table1_xm.where( (table1_xm.FLAGS==0) & 
                                  (table1_xm.FLUX_BEST > 0) &
                                  (table1_xm.FLUX_AUTO/table1_xm.FLUXERR_AUTO > min_snr))
    
    

    #Sencondly, new xmatch
    try:
        ind1_p, ind2_p = coords.indmatch(table1_tmp['X_WORLD'], table1_tmp['Y_WORLD'], 
                                     table2_tmp['ra'], table2_tmp['dec'], 
                                     error)
    except Exception, e:
        log.error("Erron in xmatch-ing tables :%s"%str(e))
        raise e

    #table_res = table1_tmp.where(ind1)

    log.info("Number of matched stars :%s"%(len(ind1_p)))
        
    return table1_tmp.where(ind1_p)['MAG_AUTO'], table2_tmp.where(ind2_p)[filter_column]


    
def generate_phot_comp_plot( input_catalog, filter, expt = 1.0 , 
                              out_filename=None, out_format='pdf'):
    """
    Generate a photometry comparison plot, comparing instrumental magnitude
    versus 2MASS photometry.
    
    Parameters
    ----------
    catalog: str 
        VOTABLE catalog  having the photometric values (instrumental and 2MASS);
    
    filter: str 
        NIR wavelength filter (J,H,K,Z)
    
    expt: float 
        exposure time of original input image; needed to compute the 
        Instrumental Magnitude (Inst_Mag)  
    out_filename: str
        filename where results will be saved;if absent, the location will be a 
        tempfile with a generated name
    out_format: str
        format of the output generated; current options available are:
        - 'pdf' (default)
        - 'jpg' 
        - 'gif'
        - for more formats, see  http://www.star.bris.ac.uk/~mbt/stilts/sun256/plot2d-usage.html
    
    Returns
    -------    
        Filename where results where saved (VOTABLE, ASCII_TABLE, ...)
    
    Example:
      > stilts tpipe  ifmt=votable cmd='addcol MyMag_K "-2.5*log10(FLUX_BEST/46.0)"' omode=out in=/home/panic/SOFTWARE/STILTS/match_double.vot out=match_D_b.vot
      > stilts plot2d in=match_S_b.vot subsetNS='j_k<=1.0 & k_snr>10' lineNS=LinearRegression xdata=k_m ydata=MyMag_k xlabel="2MASS k_m / mag" ylabel="PAPI k_m / mag"

    """
    
    log.debug("entering in <generate_phot_comp_plot>")
    ## 1 - First, we add a column with the Instrumental magnitude, computed as:
    ##    inst_mag = -2.5 log10(FLUX_BEST/TEXP)
    ## where FLUX_BEST is obtained from SExtractor output catalog
    
    input = input_catalog
    output_1 = "/tmp/output_1.xml"
    
    #create a new column
    command_line = STILTSwrapper._stilts_pathname + " tpipe " + \
                    " ifmt=votable" + \
                    " cmd='addcol Inst_Mag \"-2.5*log10(FLUX_BEST/%f)\"'" % expt + \
                    " omode=out" + " in=" + input + " out=" + output_1

    rcode = misc.utils.runCmd(command_line)
    
    if rcode == 0 or not os.path.exists(output_1):
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

def compute_regresion2( column_x, column_y , filter_name,
                        output_filename="/tmp/linear_fit.pdf"):
    """
    Compute and Plot the linear regression of two columns of the input vo_catalog
    
    Parameters
    ----------
    
    column_x: str
        column name/number for X values of the regression (MAG_AUTO)
        
    column_y: str
        column name/number for Y values of the regression (2MASS column name 
        for photometric value )
    
    Returns
    -------
    3-tuple with linear fit parameters and a Plot showing the fit
             a - intercept (ZP)
             b - slope of the linear fit
             r - estimated error
             
    None if happen some error  
    """
    
    
    X = column_x # MAG_AUTO = -2.5 * numpy.log10(table_new['FLUX_AUTO']/1.0)
    Y = column_y # 2MASS photometric value
   
   
    #remove the NaN values 
    validdata_X = ~numpy.isnan(X)
    validdata_Y = ~numpy.isnan(Y)
    validdataBoth = validdata_X & validdata_Y
    n_X = X[validdataBoth] #- (0.5*0.05) # a row extinction correction
    n_Y = Y[validdataBoth] 

    if len(n_X)<3 or len(n_Y)<3:
        raise Exception("Not enough number of data, only %d points found"%len(n_X))

    # Compute the linear fit
    res = numpy.polyfit(n_X, n_Y, 1, None, True)
    a = res[0][1] # intercept == Zero Point
    b = res[0][0] # slope
    r = res[3][1] # regression coeff ??? not exactly
    
    
    # testing robust.polyfit()
    #res = robust.polyfit(n_X, n_Y, 1)
    #a = res[0]
    #b = res[1]
    #r = 0
    # end of test
    print "Coeffs =", res
    
    # Plot the results
    pol = numpy.poly1d(res[0])
    plt.plot(n_X, n_Y, '.', n_X, pol(n_X), '-')
    plt.title("Filter %s  -- Poly fit: %f X + %f  r=%f ZP=%f" %(filter_name, b,a,r,a))
    plt.xlabel("Inst_Mag (MAG_AUTO) / mag")
    plt.ylabel("2MASS Mag  / mag")
    plt.savefig(output_filename)
    plt.show()
    log.debug("Zero Point (from polyfit) =%f", a)
    
    
    # Compute the ZP as the median of all per-star ZP=(Mag_2mass-Mag_Inst)
    zps = n_Y - n_X 
    zp = numpy.median(zps)
    #zp = a
    log.debug("Initial ZP(median) = %f"%zp)
    zp_sigma = numpy.std(zps)
    #print "ZPS=",zps
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
#    radial_distance = numpy.sqrt((table_new['X_IMAGE']-1024)**2 
#                                 + (table_new['Y_IMAGE']-1024)**2)*0.45 #arcsecs
#    
#    plt.plot( radial_distance, m_err_for_radial_systematic, '.')
#    plt.xlabel("Radial distance ('arcsec')")
#    plt.ylabel("2MASS_Mag-Inst_Mag")
#    plt.title("(1) Spatial systematics - STD = %f"%std)
#    plt.savefig("/tmp/espatial_systematics_errs.pdf")
#    plt.show()
    
    
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
    
    return (zp, b, r)


def compute_regresion( vo_catalog, column_x, column_y , 
                        output_filename="/tmp/linear_fit.pdf", min_snr=10.0):
    """
    Compute and Plot the linear regression of two columns of the input vo_catalog
    
    Parameters
    ----------
    column_x: str
        column number for X values of the regression (MAG_AUTO)
    
    column_y: str
        column number for Y values of the regression (2MASS column name for 
        photometric value )
    
    Returns
    -------
    3-tuple with linear fit parameters and a Plot showing the fit
             a - intercept (ZP)
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
    
    table_new = table.where( (table.FLAGS==0) & (table.FLUX_BEST > 0) &
                             (table.j_snr>min_snr) & (table.h_snr>min_snr) &
                             (table.k_snr>min_snr) & (table.j_k<1.0) &
                             (table.FLUX_AUTO/table.FLUXERR_AUTO>min_snr))
    
    """table_new = table.where( (table.FLAGS==0) & (table.FLUX_BEST > 0) &
                             (table.k_snr>min_snr) &
                             (table.FLUX_AUTO/table.FLUXERR_AUTO>min_snr))
    
    """
    print ">> Number of matched points = ", len(table_new)
    
    # If there aren't enough 2MASS objects, don't use color cut (J-K)<1
    if len(table_new)<25:
        #table_new = table.where( (table.FLAGS==0) & (table.FLUX_BEST > 0) &
        #                         (table.j_snr>min_snr))
        table_new = table.where( (table.FLAGS==0) & (table.FLUX_BEST > 0) &
                             (table.j_snr>min_snr) & (table.h_snr>min_snr) &
                             (table.k_snr>min_snr) & 
                             (table.FLUX_AUTO/table.FLUXERR_AUTO>min_snr))
        
        print ">> Number of matched points (no color cut)= ",len(table_new)

        
    #X = -2.5 * numpy.log10(table_new['FLUX_AUTO']/1.0)
    filter = column_y
    X = table_new[column_x] # MAG_AUTO = -2.5 * numpy.log10(table_new['FLUX_AUTO']/1.0)
    Y = table_new[column_y] # 2MASS photometric value
   
   
    #remove the NaN values 
    validdata_X = ~numpy.isnan(X)
    validdata_Y = ~numpy.isnan(Y)
    validdataBoth = validdata_X & validdata_Y
    n_X = X[validdataBoth] #- (0.5*0.05) # a row extinction correction
    n_Y = Y[validdataBoth] 

    if len(n_X)<3 or len(n_Y)<3:
        raise Exception("Not enough number of data, only %d points found"%len(n_X))

    # Compute the linear fit
    res = numpy.polyfit(n_X, n_Y, 1, None, True)
    a = res[0][1] # intercept == Zero Point
    b = res[0][0] # slope
    r = res[3][1] # regression coeff ??? not exactly
    
    # testing robust.polyfit()
    #res = robust.polyfit(n_X, n_Y, 1)
    #a = res[0]
    #b = res[1]
    #r = 0
    
    print "Coeffs =", res
    
    # Plot the results
    pol = numpy.poly1d(res[0])
    plt.plot(n_X, n_Y, '.', n_X, pol(n_X), '-')
    plt.title("Filter %s  -- Poly fit: %f X + %f  r=%f ZP=%f" %(filter, b,a,r,a))
    plt.xlabel("Inst_Mag (MAG_AUTO) / mag")
    plt.ylabel("2MASS Mag  / mag")
    plt.savefig(output_filename)
    plt.show()
    log.debug("Zero Point (from polyfit) =%f", a)
    
    
    # Compute the ZP as the median of all per-star ZP=(Mag_2mass-Mag_Inst)
    zps = n_Y - n_X 
    zp = numpy.median(zps)
    #zp = a
    log.debug("Initial ZP(median) = %f"%zp)
    zp_sigma = numpy.std(zps)
    #print "ZPS=",zps
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
    
    return (zp, b, r)

    
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
    
    _stilts_pathname = "/home/panic/SOFTWARE/PAPI/STILTS/stilts"
    #_stilts_pathname = "/home/panicmgr/PAPI/downloads/STILTS/stilts"
    
    def __init__(self, *a, **k):
        """ The constructor """
        super (STILTSwrapper, self).__init__ (*a, **k)
        
    
    def runXMatch(self, cat1, cat2, out_filename=None, out_format='votable', 
                  error=2.0):    
        """
        Do catalogs cross-match
        
        Parameters
        ----------
        cat1, cat2: str
            catalogs for cross-matching
        err: float
            max. error for finding objects within (arcseconds)
        out_filename: str
            file name where results will be saved; if absent, the location will 
            be a tempfile with a generated name.
        out_format: str
            format of the output generated; current options available are:
            - VO Table (XML) (votable) (default)
            - SVC (Software handshaking structure) message (svc)
            - ASCII table (ascii)
        
        Returns
        -------
        Filename where results where saved (VOTABLE, ASCII_TABLE, ...)
        """
        
        in1 = cat1
        in2 = cat2
        out = out_filename
        
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

def doPhotometry(input_image, catalog, output_filename, snr, zero_point=0.0):
    """
    Run the rough photometric calibraiton
    
    Parameters
    ----------
    input_image: str 
        filename reduced science image
    
    catalog: str
        photometric catalog to compare to (2MASS, USNO-B, ...)
    
    output_filename: str
        filename where output figure will be saved
    
    snr: float
        minimum SNR of stars used of linear fit
    
    zero_point: float 
        initial magnitud zero point for SExtractor (default 0.0)
    
    Returns
    -------
    if all was ok, return ouput_filename
    """
    
    
    ## 0.1 - Read the RA,Dec and TEXP values from the input image
    log.debug("Reading FITS file EXPT, RA, Dec & FILTER values ...")
    try:
        my_fits = datahandler.ClFits(input_image)
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
            raise Exception("Found wrong RA,Dec,EXPTIME or FILTER value.")
        
        # Checking Filter    
        if filter == 'J':
            two_mass_col_name = 'j_m'
        elif filter == 'H':
            two_mass_col_name = 'h_m'
        elif filter == 'K' or filter == 'K-PRIME' or filter == 'KS':
            two_mass_col_name = 'k_m'   
        #elif filter=='OPEN':
        #    two_mass_col_name = 'j_m'
        else:
            log.error("Filter %s not supported" %filter)
            raise Exception("Filter %s not supported"%filter)
            
    except Exception, e:
        log.error("Cannot read properly FITS file : %s:", str(e))
        raise e
    
    ## 0.2 - Generate image catalog (VOTable) -> SExtractor
    log.debug("*** Creating SExtractor VOTable catalog ....")
    
    #tmp_fd, tmp_name = tempfile.mkstemp(suffix='.xml', dir=os.getcwd())
    #os.close(tmp_fd)
    image_catalog = os.path.splitext(input_image)[0]  + ".xml"
    
    sex = astromatic.SExtractor()
    sex.ext_config['CHECKIMAGE_TYPE'] = "NONE"
    sex.config['CATALOG_TYPE'] = "ASCII_VOTABLE"
    sex.config['CATALOG_NAME'] = image_catalog
    sex.config['DETECT_THRESH'] = 1.5
    sex.config['DETECT_MINAREA'] = 5
    sex.config['MAG_ZEROPOINT'] = zero_point
    
    
    try:
        sex.run(input_image, updateconfig=True, clean=False)
    except Exception, e:
        log.error("Canno't create SExtractor catalog : %s", str(e))
        raise e 
    
      
    ## 1 - Generate region of base catalog (2MASS)
    icat = catalog_query.ICatalog ()
    out_base_catalog = os.getcwd() + "/catalog_region.xml"
    pixel_scale = 0.45 # Panic at 2.2m
    sc_area = 0.9 # percentage of area used (default 90%) 
    search_radius = numpy.minimum(my_fits.getNaxis1(), 
                                  my_fits.getNaxis2())*pixel_scale*sc_area/2
    try:
        i_catalog = icat.queryCatalog(ra, dec, search_radius, 
                                     catalog_query.ICatalog.cat_names['2MASS'], 
                                     out_base_catalog, 'votable')[0]
        log.debug("Output file generated : %s", i_catalog) 
    except Exception, e:
        log.error("Sorry, cann't solve the query to ICatalog: %s", str(e))
        raise e

    ## 2- XMatch the catalogs (image_catalog VS just base_catalog generated)          
    out_xmatch_file = os.getcwd() + "/xmatch.xml"
    try:
        #match_cat = catalog_xmatch( image_catalog, i_catalog, 
        #                           out_xmatch_file, out_format='votable',
        #                           error=1.0 ) 
        
        xm, ym = catalog_xmatch2( image_catalog, i_catalog, two_mass_col_name, 
                                   error=1.0 , min_snr=snr)
        log.debug("XMatch done !")
    except Exception, e:
        log.error("XMatch failed %s", str(e))
        raise e
    
    ## 3a- Compute the linear regression (fit) of Inst_Mag VS 2MASS_Mag
    ###### and generate the plot file with the photometric comparison
    ###### using Numpy & Matplotlib
    ###### 2MASS_Mag = Inst_Mag*b + ZP  
    log.debug("Compute & Plot regression !!!")    
    est_zp_err = zero_point
    try:
        #est_zp_err = compute_regresion(match_cat, 'MAG_AUTO', 
        #                   two_mass_col_name, output_filename, snr )[0]
        est_zp_err = compute_regresion2(xm, ym, filter, output_filename)[0] 
        
        log.info("Estimated ZP_err=%s"%est_zp_err)
        #sys.exit(0)
    except Exception, e:
        log.error("Sorry, some error while computing linear fit or\
        ploting the results: %s", str(e))
        raise e
    
    ## 3b- Compute the linear regression (fit) of Inst_Mag VS 2MASS_Mag
    ###### and generate the plot file with the photometric comparison 
    ###### using STILTS
    try:
        exptime = 1.0 # SWARP normalize flux to 1 sec
        #file_ext = os.path.splitext(output_filename)[1]
        #output_filename_2 = output_filename.replace(file_ext, "_b"+file_ext) 
        #plot_file = generate_phot_comp_plot( match_cat, two_mass_col_name, exptime, 
        #                                      output_filename_2, 
        #                                      out_format='pdf')
        #log.debug("Plot file generated : %s", plot_file) 
    except Exception, e:
        log.error("Sorry, some error while computing linear fit or \
        ploting the results: %s", str(e))
        raise e

    ## 4.0 - Finally, generate image catalog (VOTable) -> SExtractor with the 
    ## new estimated ZP
    log.info("*** Creating SExtractor VOTable catalog (ZP=%f)...."%(est_zp_err+zero_point))
    image_catalog = os.path.splitext(input_image)[0]  + ".xml"
    
    sex = astromatic.SExtractor()
    sex.ext_config['CHECKIMAGE_TYPE'] = "NONE"
    sex.config['CATALOG_TYPE'] = "ASCII_VOTABLE"
    sex.config['CATALOG_NAME'] = image_catalog
    sex.config['DETECT_THRESH'] = 1.5
    sex.config['DETECT_MINAREA'] = 5
    sex.config['MAG_ZEROPOINT'] = est_zp_err + zero_point
    
    
    try:
        sex.run(input_image, updateconfig=True, clean=False)
    except Exception, e:
        log.error("Canno't create SExtractor catalog : %s", str(e))
        raise e 

    return output_filename

################################################################################
# main
################################################################################
if __name__ == "__main__":

    log.debug( 'Testing Photometric calibration comparison with 2MASS')
    
    # Get and check command-line options
        
    usage = "usage: %prog [options] arg1 arg2 ..."
    desc = """This module receives a reduced image of any known NIR filter and
match to 2MASS catalog performing a fit in order to get a estimation of the 
Zero Point.""" 
    parser = OptionParser(usage, description=desc)
    
    parser.add_option("-i", "--input_image",
                  action="store", dest="input_image", 
                  help="Input image to calibrate to do photometric comparison with")
                  
    parser.add_option("-c", "--base_catalog (2MASS, USNO-B)",
                  action="store", dest="base_catalog",
                  help="Name of base catalog to compare with (2MASS, USNO-B) -- not used !!! (default = %default)",
                  default="2MASS")
    
    parser.add_option("-S", "--snr",
                  action="store", dest="snr", type=int,
                  help="Min SNR of stars used for linear fit (default = %default)",
                  default=10.0)
    
    parser.add_option("-z", "--zero_point",
                  action="store", dest="zero_point", type=float, default=25.0,
                  help="Initial Magnitude Zero Point estimation [%default]; used for SExtractor")
                  
    parser.add_option("-o", "--output",
                  action="store", dest="output_filename", 
                  help="Output plot filename (default = %default)",
                  default="photometry.pdf")
    
    
                                
    (options, args) = parser.parse_args()
    
    
    if not options.input_image or len(args)!=0: 
    # args is the leftover positional arguments after all options have been processed
        parser.print_help()
        parser.error("wrong number of arguments " )
    if not options.output_filename:
        options.output_filename = None
    

    if not os.path.exists(options.input_image):
        log.error ("Input image %s does not exist", options.input_image)
        sys.exit(0)
        
    try:
        catalog = "2MASS"
        doPhotometry(options.input_image, catalog, options.output_filename, 
                     options.snr, options.zero_point)
    except Exception, e:
        log.info("Some error while running photometric calibration: %s"%str(e))
        sys.exit(0)
        
