""" My module for various data analysis tasks.

:REQUIREMENTS: :doc:`numpy`


2008-07-25 16:20 IJC: Created.

2009-12-08 11:31 IJC: Updated transit flag in planet objects and 
                      :func:`rveph` function.

2010-02-18 14:06 IJC: Added :func:`medianfilter`

2010-08-03 15:38 IJC: Merged versions.

2010-10-28 11:53 IJMC: Updated documentation strings for Sphinx;
                      moved pylab import inside individual
                      function.

2011-04-13 14:29 IJMC: Added keyword functionality to :func:`fmin`
                      (taken from scipy.optimize).

2011-04-14 09:48 IJMC: Added a few fundamental constants.

2011-04-22 14:48 IJMC: Added :func:`trueanomaly` and
                       :func:`eccentricanomaly`.

2011-06-21 13:15 IJMC: Added :func:`get_t23` and :func:`get_t14` to
                       planet objects.
"""


from numpy import ones, std, sum, mean, median, array, linalg, tile, concatenate, floor, Inf, arange, meshgrid, zeros, sin, cos, tan, arctan, sqrt, exp, nan, max
import pdb
import numpy as np

from scipy import optimize


c   = 299792458  # speed of light, m/s
h = 6.626068e-34 # SI units: Planck's constant
k = 1.3806503e-23 # SI units: Boltzmann constant
G = 6.67300e-11 # SI units: Gravitational constant

pi  = 3.14159265358979
AU  = 149597870691.0  # AU in meters
day = 86400.0  # seconds in a Julian Day
rsun = 6.95508e8  # Solar mean radius, in m
msun = 1.9891e30  # solar mass in kg
rearth = 6.378136e6 # Earth's equatorial radius in m; [Allen's]
mearth = 5.9737e24 # in kg; [Allen's]
rjup = 7.1492e7 # Jupiter equatorial radius in m
mjup = 1898.7e24  # Jupiter mass in kg
pc = 3.08568025e16 # parsec in meters


class planet:
    """Very handy planet object.

    Best initialized using :func:`getobj`.

    :REQUIREMENTS: Database file `exoplanets.csv` from http://exoplanets.org/

    """
    # 2010-03-07 22:21 IJC: Created
    # 2011-01-22 17:15 IJMC: Updated format b/c CPS folks keep changing
    #                       their database file format

    # 2011-05-19 16:11 IJMC: Updated format b/c CPS folks keep
    #                       changing their database file format -- but
    #                       it's almost worth it to finally have
    #                       stellar radii.
    def __init__(self, *args):
        keys = ['name','comp','ncomp','mult','discmeth','firstref','firsturl','date','jsname','etdname','per','uper','t0','ut0','ecc','uecc','ueccd','om','uom','k','uk','msini','umsini','a','ua','orbref','orburl','transit','t14','ut14','tt','utt','ar','uar','uard','i','ui','uid','b','ub','ubd','depth','udepth','udepthd','r','ur','density','udensity','gravity','ugravity','transitref','transiturl','trend','dvdt','udvdt','freeze_ecc','rms','chi2','nobs','star','hd','hr','hipp','sao','gl','othername','sptype','binary','v','bmv','j','h','ks','ra','dec','ra_string','dec_string','rstar', 'urstar', 'urstard', 'rstarref', 'rstarurl','mstar','umstar','umstard','teff','uteff','vsini','uvsini','fe','ufe','logg','ulogg','shk','rhk','par','upar','distance','udistance','lambd', 'ulambd', 'massref','massurl','specref','specurl','distref','disturl','simbadname','nstedid','binaryref', 'binaryurl']

        if len(keys)<>len(args):
            print "Incorrect number of input arguments (%i, but should be %i)" % (len(args), len(keys))
            return None

        for key,arg in zip(keys, args):
            try:
                temp = float(arg)+1
                isnumber = True
            except:
                isnumber = False
                

            if isnumber:
                exec('self.%s=%s' % (key,arg) )
            else:
                exec('self.%s="%s"' % (key,arg) )

        return None

    def get_t23(self, *args):
        """Compute full transit duration (in days) for a transiting planet.

        Returns:
           nan if required fields are missing.

        Using Eq. 15 of J. Winn's chapter in S. Seager's book "Exoplanets." 

        :SEE ALSO:
          :func:`get_t14`
          """ 
        # 2011-06-21 12:53 IJMC: Created

        # Get necessary parameters:
        per, ra, k, b, inc, ecc, om = self.per, 1./self.ar, sqrt(self.depth), self.b, self.i, self.ecc, self.om

        inc *= np.pi/180.
        om *= np.pi/180.

        ret = 0.

        # Check that no parameters are empty:
        if per=='' or k=='':
            ret = nan

        if b=='':
            try:
                b = (cos(inc)/ra) * (1. - ecc**2) / (1. + ecc * sin(om))
            except:
                ret = nan
        elif ra=='':
            try:
                ra = (cos(inc)/b) * (1. - ecc**2) / (1. + ecc * sin(om))
            except:
                ret = nan
        elif inc=='':
            try:
                inc = np.arccos((b * ra) / ((1. - ecc**2) / (1. + ecc * sin(om))))
            except:
                ret = nan
        
        if np.isnan(ret):
            print "Could not compute t_23.  Parameters are:"
            print "period>>", per
            print "r_*/a>>", ra
            print "r_p/r_*>>", k
            print "impact parameter (b)>>", b
            print "inclination>>", inc*180./np.pi, " deg"
            print "eccentricity>>", ecc
            print "omega>>", om*180./np.pi, " deg"
        else:
            ret = (per/np.pi) * np.arcsin(ra * np.sqrt((1. - k)**2 - b**2) / np.sin(inc))

        return ret


    def get_t14(self, *args):
        """Compute total transit duration (in days) for a transiting planet.

        Returns:
           nan if required fields are missing.

        Using Eq. 14 of J. Winn's chapter in S. Seager's book "Exoplanets." 

        :SEE ALSO:
          :func:`get_t23`
          """ 
        # 2011-06-21 12:53 IJMC: Created

        # Get necessary parameters:
        per, ra, k, b, inc, ecc, om = self.per, 1./self.ar, sqrt(self.depth), self.b, self.i, self.ecc, self.om

        inc *= np.pi/180.
        om *= np.pi/180.

        ret = 0.

        # Check that no parameters are empty:
        if per=='' or k=='':
            ret = nan

        if b=='':
            try:
                b = (cos(inc)/ra) * (1. - ecc**2) / (1. + ecc * sin(om))
            except:
                ret = nan
        elif ra=='':
            try:
                ra = (cos(inc)/b) * (1. - ecc**2) / (1. + ecc * sin(om))
            except:
                ret = nan
        elif inc=='':
            try:
                inc = np.arccos((b * ra) / ((1. - ecc**2) / (1. + ecc * sin(om))))
            except:
                ret = nan
        
        if np.isnan(ret):
            print "Could not compute t_14.  Parameters are:"
            print "period>>", per
            print "r_*/a>>", ra
            print "r_p/r_*>>", k
            print "impact parameter (b)>>", b
            print "inclination>>", inc*180./np.pi, " deg"
            print "eccentricity>>", ecc
            print "omega>>", om*180./np.pi, " deg"
        else:
            ret = (per/np.pi) * np.arcsin(ra * np.sqrt((1. + k)**2 - b**2) / np.sin(inc))

        return ret


    def rv(self, jd):
        """Compute radial velocity on a planet object for given Julian Date.

        :EXAMPLE:
           ::

               import analysis
               p = analysis.getobj('HD 189733 b')
               jd = 2400000.
               print p.rv(jd)
        
        refer to function `analysis.rv` for full documentation.

        SEE ALSO: :func:`analysis.rv`, :func:`analysis.mjd`
        """
        return rv(self, jd)

    def rveph(self, jd):
        """Compute the most recently elapsed RV emphemeris of a given
        planet at a given JD.  RV ephemeris is defined by the having
        radial velocity equal to zero.
        
        refer to  :func:`analysis.rv` for full documentation.
        
        SEE ALSO: :func:`analysis.getobj`, :func:`analysis.phase`
        """
        return rveph(self, jd)

    def phase(self, hjd):
        """Get phase of an orbiting planet.

        refer to function `analysis.getorbitalphase` for full documentation.

        SEE ALSO: :func:`analysis.getorbitalphase`, :func:`analysis.mjd`
        """
        # 2009-09-28 14:07 IJC: Implemented object-oriented version
        return getorbitalphase(self, hjd)


def loaddata(filelist, path='', band=1):
    """Load a set of reduced spectra into a single data file.

   datalist  = ['file1.fits', 'file2.fits']
   datapath = '~/data/'

   data = loaddata(datalist, path=datapath, band=1)


   The input can also be the name of an IRAF-style file list.
   """
    # 2008-07-25 16:26 IJC: Created
    # 2010-01-20 12:58 IJC: Made function work for IRTF low-res data.
    #                       Replaced 'warn' command with 'print'.
    # 2011-04-08 11:42 IJC: Updated; moved inside analysis.py.
    # 2011-04-12 09:57 IJC: Fixed misnamed imports

    import pyfits

    data = array([])

    if filelist.__class__==str or isinstance(filelist,np.string_):
        filelist = ns.file2list(filelist)
    elif filelist.__class__<>list:
        print('Input to loaddata must be a python list or string')
        return data

    num = len(filelist)

    # Load one file just to get the dimensions right.
    irspec = pyfits.getdata(filelist[0])

    ii = 0
    for element in filelist:
        irspec = pyfits.getdata(element)

        if ii==0:
            irsh = irspec.shape
            data  = zeros((num,)+irsh[1::], float)
            
        if len(irsh)>2:
            for jj in range(irsh[1]):
                data[ii, jj, :] = irspec[band-1, jj, :]
        else:
            data[ii,:] = irspec[band-1,:]

        ii=ii+1

    return data


def getobj(*args, **kw):
    """Get data for a specified planet.  

       :INPUTS:  (str) -- planet name, e.g. "HD 189733 b"

       :OPTIONAL INPUTS:  
           datafile : str 
                datafile name

           verbose : bool
                verbosity flag                          

       :EXAMPLE:
           ::

                   p = getobj('55cnce')
                   p.period   ---> 2.81705

       The attributes of the returned object are many and varied, and
       can be listed using the 'dir' command on the returned object.

       This looks up data from the local datafile, which could be out
       of date.

       SEE ALSO: :func:`rv`"""
    # 2008-07-30 16:56 IJC: Created
    # 2010-03-07 22:24 IJC: Updated w/new exoplanets.org data table! 
    # 2010-03-11 10:01 IJC: If planet name not found, return list of
    #                       planet names.  Restructured input format.
    # 2010-11-01 13:30 IJC: Added "import os"
    # 2011-05-19 15:56 IJC: Modified documentation.

    from pylab import find
    import os

    if kw.has_key('datafile'):
        datafile=kw['datafile']
    else:
        datafile=os.path.expanduser('~/python/exoplanets.csv')
    if kw.has_key('verbose'):
        verbose = kw['verbose']
    else:
        verbose=False

    if len(args)==0:
        inp = 'noplanetname'
    else:
        inp = args[0]

    if verbose:  print "datafile>>" + datafile
    f = open(datafile, 'r')
    data = f.readlines()
    f.close()
    data = data[1::]  # remove header line
    names = array([line.split(',')[0] for line in data])

    foundobj = (names==inp)

    if (not foundobj.any()):
        if verbose: print "could not find desired planet; returning names of known planets"
        return names
    else:
        planetind = int(find(foundobj))
        pinfo  = data[planetind].strip().split(',')
        myplanet = planet(*pinfo)
        
    return myplanet
 

def getorbitalphase(planet, hjd, **kw):
    """Get phase of an orbiting planet.
    
    INPUT: 
       planet  -- a planet from getobj; e.g., getobj('55cnce')
       hjd

    OUTPUT:
       orbital phase -- from 0 to 1.

    NOTES: 
       If planet.transit==True, phase is based on the transit time ephemeris.  
       If planet.transit==False, phase is based on the RV ephemeris as
          computed by function rveph
    

    SEE ALSO:  :func:`getobj`, :func:`mjd`, :func:`rveph`
    """

    hjd = array(hjd).copy()
    if bool(planet.transit)==True:
        thiseph = planet.tt
    else:
        thiseph = planet.rveph(hjd.max())

    orbphase = ((hjd - thiseph) ) / planet.per
    orbphase -= int(orbphase.mean())

    return orbphase

def mjd(date):
    """Convert Julian Date to Modified Julian Date, or vice versa.

    if date>=2400000.5, add 2400000.5
    if date<2400000.5, subtract 2400000.5
    """
    # 2009-09-24 09:54 IJC: Created 

    date = array(date, dtype=float, copy=True, subok=True)
    offset = 2400000.5
    if (date<offset).all():
        date += offset
    elif (date>=offset).all():
        date -= offset
    else:
        print "Input contains date both below and above %f" % offset

    return date

def rveph(p, jd):
    """Compute the most recently elapsed RV emphemeris of a given
    planet at a given JD.  RV ephemeris is defined by the having
    radial velocity equal to zero.



        :EXAMPLE:
            ::

                from analysis import getobj, rveph
                jd = 2454693            # date: 2008/08/14
                p  = getobj('55cnce')   # planet: 55 Cnc e
                t = rveph(p, jd)
            
            returns t ~ 


        SEE ALSO: :func:`getobj`, :func:`phase`
        """
    # 2009-12-08 11:20 IJC: Created.  Ref: Beauge et al. 2008 in
    #    "Extrasolar Planets," R. Dvorak ed.
    # 2010-03-12 09:34 IJC: Updated for new planet-style object.

    from numpy import cos, arccos, arctan, sqrt, tan, pi, sin, int

    if p.__class__<>planet:
        raise Exception, "First input must be a 'planet' object."
    
    omega = p.om*pi/180 
    tau = p.t0
    ecc = p.ecc
    per = p.per
    f = arccos(-ecc * cos(omega)) - omega  # true anomaly
    u = 2.*arctan(sqrt((1-ecc)/(1+ecc))*tan(f/2.)) # eccentric anomaly
    n = 2*pi/per

    time0 = tau+ (u-ecc*sin(u))/n
    norb = int((time0-jd)/per)
    time = time0-norb*per
    return time


def rv(p, jd):
    """ Compute unprojected astrocentric RV of a planet for a given JD in m/s.

        :EXAMPLE:
            ::

                jd = 2454693            # date: 2008/08/14
                p  = getobj('55 Cnc e')   # planet: 55 Cnc e
                vp = rv(p, jd)
            
            returns vp ~ 1.47e5  [m/s]

        The result will need to be multiplied by the sine of the
        inclination angle (i.e., "sin i").  Positive radial velocities
        are directed _AWAY_ from the observer.

        To compute the barycentric radial velocity of the host star,
        scale the returned values by the mass ratio -m/(m+M).

        SEE ALSO: :func:`getobj`
    """
    # 2008-08-13 12:59 IJC: Created with help from Debra Fischer,
    #    Murray & Dermott, and Beauge et al. 2008 in "Extrasolar
    #    Planets," R. Dvorak ed.
    # 2008-09-25 12:55 IJC: Updated documentation to be more clear.
    # 2009-10-01 10:25 IJC: Moved "import" statement within func.
    # 2010-03-12 09:13 IJC: Updated w/new planet-type objects

    if p.__class__<>planet:
        raise Exception, "First input must be a 'planet' object."
    jd = array(jd, copy=True, subok=True)
    if jd.shape==():  
        singlevalueinput = True
        jd = array([jd])
    else:
        singlevalueinput = False
    if ((jd-2454692) > 5000).any():
        raise Exception, "Implausible Julian Date entered."

    try:
        per = p.per
        tau = p.t0
        ecc = p.ecc
        a   = p.a
        omega = p.om * pi/180.0
    except:
        raise Exception, "Could not load all desired planet parameters."

    n = 2.*pi/per
    m = n*(jd - tau)
    e = []
    
    for element in m:
        def kep(e): return element - e + ecc*sin(e)
        e.append(optimize.fsolve(kep, 0))

    e = array(e)
    f = 2. * arctan(  sqrt((1+ecc)/(1.-ecc)) * tan(e/2.)  )
    #r = a*(1-ecc*cos(e))
    #x = r*cos(f)
    #y = r*sin(f)

    K = n * a / sqrt(1-ecc**2)

    vz = -K * ( cos(f+omega) + ecc*cos(omega) )
    vzms = vz * AU/day  # convert to m/s

    if singlevalueinput:
        vzms = vzms[0]

    return vzms



def dopspec(starspec, planetspec, starrv, planetrv, disp, starphase=[], planetphase=[], wlscale=True):
    """ Generate combined time series spectra using planet and star
    models, planet and star RV profiles.

    D = dopspec(sspec, pspec, sRV, pRV, disp, sphase=[], pphase=[])

    :INPUTS:
       sspec, pspec:  star, planet spectra; must be on a common
                         logarithmic wavelength grid 
       sRV, pRV:      star, planet radial velocities in m/s
       disp:          constant logarithmic dispersion of the wavelength 
                         grid: LAMBDA_i/LAMBDA_(i-1)

    :OPTIONAL INPUTS:
       sphase, pphase:  normalized phase functions of star and planet.  
                           The inputs sspec and pspec will be scaled
                           by these values for each observation.
       wlscale:         return relative wavelength scale for new data

    NOTE: Input spectra must be linearly spaced in log wavelength and increasing:
            that is, they must have [lambda_i / lambda_(i-1)] = disp =
            constant > 1
          Positive velocities are directed AWAY from the observer."""

#2008-08-19 16:30 IJC: Created

# Options: 1. chop off ends or not?  2. phase function. 
#   4. Noise level 5. telluric? 6. hold star RV constant


# Initialize:
    starspec   = array(starspec  ).ravel()
    planetspec = array(planetspec).ravel()
    starrv     = array(starrv    ).ravel()
    planetrv   = array(planetrv  ).ravel()

    ns = len(starspec)
    nr = len(starrv)

    starphase   = array(starphase  ).ravel()
    planetphase = array(planetphase).ravel()
    if len(starphase)==0:
        starphase = ones(nr, float)
    if len(planetphase)==0:
        planetphase = ones(nr, float)
    

# Check inputs: 
    if ns<>len(planetspec):
        raise Exception, "Star and planet spectra must be same length."
    if nr<>len(planetrv):
        raise Exception, "Star and planet RV profiles must be same length."

    logdisp = log(disp)

# Calculate wavelength shift limits for each RV point
    sshift = ( log(1.0+starrv  /c) / logdisp ).round() 
    pshift = ( log(1.0+planetrv/c) / logdisp ).round() 

    limlo =  int( concatenate((sshift, pshift)).min() )
    limhi =  int( concatenate((sshift, pshift)).max() )

    ns2 = ns + (limhi - limlo)

    data = zeros((nr, ns2), float)

# Iterate over RV values, constructing spectra
    for ii in range(nr):
        data[ii, (sshift[ii]-limlo):(ns+sshift[ii]-limlo)] = \
            starphase[ii] * starspec
        data[ii, (pshift[ii]-limlo):(ns+pshift[ii]-limlo)] = \
            data[ii, (pshift[ii]-limlo):(ns+pshift[ii]-limlo)] + \
            planetphase[ii] * planetspec

    if wlscale:
        data = (data, disp**(arange(ns2) + limlo))
    return data



def loadatran(filename, wl=True, verbose=False):
    """ Load ATRAN atmospheric transmission data file.

        t = loadatran(filename, wl=True)

        INPUT: 
           filename -- filename of the ATRAN file.  This should be an
                       ASCII array where the second column is
                       wavelength and the third is the atmospheric
                       transmission.
                     (This can also be a list of filenames!)  
           
        :OPTIONAL INPUT: 
           wl -- if True (DEFAULT) also return the wavelength scale.
                 This can take up to twice as long for large files.

        RETURNS:
           if wl==True:   returns a 2D array, with columns [lambda, transmission]
           if wl==False:  returns a 1D Numpy array of the transmission
           
           NOTE: A list of these return values is created if
                 'filename' is actually an input list."""
    # 2008-08-21 09:42 IJC: Created to save myself a bit of time
    # 2008-08-25 10:08 IJC: Read in all lines at once; should go
    #                        faster with sufficient memory
    # 2008-09-09 13:56 IJC: Only convert the wavelength and flux
    #                        columns(#1 & #2) -- speed up slightly.
    if filename.__class__==list:
        returnlist = []
        for element in filename:
            returnlist.append(loadatran(element, wl=wl))
        return returnlist
    
    f = open(filename, 'r')
    dat = f.readlines()
    f.close()
    if verbose:
        print dat[0]
        print dat[0].split()
        print dat[0].split()[1:3]
        print dat[0].split()[2]
    
    if wl:
        data = array([map(float, line.split()[1:3]) for line in dat])
    else:
        data = array([float(line.split()[2]) for line in dat])

    return data

def poly2cheby(cin):
    """Convert straight monomial coefficients to chebychev coefficients.

    INPUT: poly coefficients (e.g., for use w/polyval)
    OUTPUT: chebyt coefficients

    SEE ALSO: :func:`gpolyval`; scipy.special.chebyt
    """
    # 2009-07-07 09:41 IJC: Created
    from scipy.special import poly1d, chebyt
    
    cin = poly1d(cin)
    cout = []

    ord = cin.order
    for ii in range(ord+1):
        chebyii = chebyt(ord-ii)
        cout.append(cin.coeffs[0]/chebyii.coeffs[0])
        cin -= chebyii*cout[ii]

    return cout

def cheby2poly(cin):
    """Convert  chebychev coefficients to 'normal' polyval coefficients .

    :INPUT: chebyt coefficients 

    :OUTPUT: poly coefficients (e.g., for use w/polyval)

    :SEE ALSO: :func:`poly2cheby`, :func:`gpolyval`; scipy.special.chebyt
    """
    # 2009-10-22 22:19 IJC: Created
    from scipy.special import poly1d, chebyt
    
    cin = poly1d(cin)
    cout = poly1d(0)

    ord = cin.order
    for ii in range(ord+1):
        cout += chebyt(ii)*cin[ii]

    return cout


def gpolyval(c,x, mode='poly', retp=False):
    """Generic polynomial evaluator.

    INPUT:
            c (1D array) -- coefficients of polynomial to evaluate,
                          from highest-order to lowest.
            x (1D array) -- pixel values at which to evaluate C

    OPTINAL INPUTS:
            MODE='poly' -- 'poly' uses standard monomial coefficients
                            as accepted by, e.g., polyval.  Other
                            modes -- 'cheby' (1st kind) and 'legendre'
                            (P_n) -- convert input 'x' to a normalized
                            [-1,+1] domain
            RETP=False  -- Return coefficients as well as evaluated poly.

    OUTPUT:
            y -- array of shape X; evaluated polynomial.  
            (y, p)  (if retp=True)

    SEE ALSO:  :func:`poly2cheby`

     """
    # 2009-06-17 15:42 IJC: Created 
    # 2011-12-29 23:11 IJMC: Fixed normalization of the input 'x' array
    from scipy import special
    from numpy import zeros, polyval, concatenate


    c = array(c).copy()
    nc = len(c)

    if mode=='poly':
        totalcoef = c
        ret = polyval(totalcoef, x)
    elif mode=='cheby':
        x = 2. * (x - 0.5*(x.max() + x.min())) / (x.max() - x.min())
        totalcoef = special.poly1d([0])
        for ii in range(nc):
            totalcoef += c[ii]*special.chebyt(nc-ii-1)
        ret = polyval(totalcoef, x)
    elif mode=='legendre':
        x = 2. * (x - 0.5*(x.max() + x.min())) / (x.max() - x.min())
        totalcoef = special.poly1d([0])
        for ii in range(nc):
            totalcoef += c[ii]*special.legendre(nc-ii-1)
        ret = polyval(totalcoef, x)

    if retp:
        return (ret, totalcoef)
    else:
        return ret
    return -1

def stdr(x, nsigma=3, niter=Inf, finite=True, verbose=False, axis=None):
    """Return the standard deviation of an array after removing outliers.
    
    :INPUTS:
      x -- (array) data set to find std of

    :OPTIONAL INPUT:
      nsigma -- (float) number of standard deviations for clipping
      niter -- number of iterations.
      finite -- if True, remove all non-finite elements (e.g. Inf, NaN)
      axis -- (int) axis along which to compute the mean.

    :EXAMPLE:
      ::

          from numpy import *
          from analysis import stdr
          x = concatenate((randn(200),[1000]))
          print std(x), stdr(x, nsigma=3)
          x = concatenate((x,[nan,inf]))
          print std(x), stdr(x, nsigma=3)

    SEE ALSO: :func:`meanr`, :func:`medianr`, :func:`removeoutliers`, 
              :func:`numpy.isfinite`
    """
    # 2010-02-16 14:57 IJC: Created from mear
    # 2010-07-01 14:06 IJC: ADded support for higher dimensions
    from numpy import array, isfinite, zeros, swapaxes

    x = array(x, copy=True)
    xshape = x.shape
    ndim =  len(xshape)
    if ndim==0:
        return x

    if  ndim==1 or axis is None:
        # "1D" array
        x = x.ravel()
        if finite:
            x = x[isfinite(x)]
        x = removeoutliers(x, nsigma, niter=Inf, verbose=verbose)
        return x.std()
    else:
        newshape = list(xshape)
        oldDimension = newshape.pop(axis) 
        ret = zeros(newshape, float)

        # Prevent us from taking the action along the axis of primary incidices:
        if axis==0: 
            x=swapaxes(x,0,1)

        if axis>1:
            nextaxis = axis-1
        else:
            nextaxis = 0

        for ii in range(newshape[0]):
            #print 'x.shape>>',x.shape, 'newshape>>',newshape, 'x[ii].shape>>',x[ii].shape, 'ret[ii].shape>>',ret[ii].shape,'ii>>',ii
            ret[ii] = stdr(x[ii], nsigma=nsigma,niter=niter,finite=finite,\
                                 verbose=verbose, axis=nextaxis)
        return ret


def meanr(x, nsigma=3, niter=Inf, finite=True, verbose=False,axis=None):
    """Return the mean of an array after removing outliers.
    
    :INPUTS:
      x -- (array) data set to find mean of

    :OPTIONAL INPUT:
      nsigma -- (float) number of standard deviations for clipping
      niter -- number of iterations.
      finite -- if True, remove all non-finite elements (e.g. Inf, NaN)
      axis -- (int) axis along which to compute the mean.

    :EXAMPLE:
      ::

          from numpy import *
          from analysis import meanr
          x = concatenate((randn(200),[1000]))
          print mean(x), meanr(x, nsigma=3)
          x = concatenate((x,[nan,inf]))
          print mean(x), meanr(x, nsigma=3)

    SEE ALSO: :func:`medianr`, :func:`stdr`, :func:`removeoutliers`, 
              :func:`numpy.isfinite`
    """
    # 2009-10-01 10:44 IJC: Created
    # 2010-07-01 13:52 IJC: Now handles higher dimensions.
    from numpy import array, isfinite, zeros, swapaxes

    x = array(x, copy=True)
    xshape = x.shape
    ndim =  len(xshape)
    if ndim==0:
        return x

    if  ndim==1 or axis is None:
        # "1D" array
        x = x.ravel()
        if finite:
            x = x[isfinite(x)]
        x = removeoutliers(x, nsigma, niter=Inf, verbose=verbose)
        return x.mean()
    else:
        newshape = list(xshape)
        oldDimension = newshape.pop(axis) 
        ret = zeros(newshape, float)

        # Prevent us from taking the mean along the axis of primary incidices:
        if axis==0: 
            x=swapaxes(x,0,1)

        if axis>1:
            nextaxis = axis-1
        else:
            nextaxis = 0

        for ii in range(newshape[0]):
            #print 'x.shape>>',x.shape, 'newshape>>',newshape, 'x[ii].shape>>',x[ii].shape, 'ret[ii].shape>>',ret[ii].shape,'ii>>',ii
            ret[ii] = meanr(x[ii], nsigma=nsigma,niter=niter,finite=finite,\
                                 verbose=verbose, axis=nextaxis)
        return ret
            
        


def medianr(x, nsigma=3, niter=Inf, finite=True, verbose=False,axis=None):
    """Return the median of an array after removing outliers.
    
    :INPUTS:
      x -- (array) data set to find median of

    :OPTIONAL INPUT:
      nsigma -- (float) number of standard deviations for clipping
      niter -- number of iterations.
      finite -- if True, remove all non-finite elements (e.g. Inf, NaN)
      axis -- (int) axis along which to compute the mean.

    :EXAMPLE:
      ::

          from numpy import *
          from analysis import medianr
          x = concatenate((randn(200),[1000]))
          print median(x), medianr(x, nsigma=3)
          x = concatenate((x,[nan,inf]))
          print median(x), medianr(x, nsigma=3)

    SEE ALSO: :func:`meanr`, :func:`stdr`, :func:`removeoutliers`, 
              :func:`numpy.isfinite`
    """
    # 2009-10-01 10:44 IJC: Created
    #2010-07-01 14:04 IJC: Added support for higher dimensions
    from numpy import array, isfinite, zeros, median, swapaxes

    x = array(x, copy=True)
    xshape = x.shape
    ndim =  len(xshape)
    if ndim==0:
        return x

    if  ndim==1 or axis is None:
        # "1D" array
        x = x.ravel()
        if finite:
            x = x[isfinite(x)]
        x = removeoutliers(x, nsigma, niter=Inf, verbose=verbose)
        return median(x)
    else:
        newshape = list(xshape)
        oldDimension = newshape.pop(axis) 
        ret = zeros(newshape, float)

        # Prevent us from taking the action along the axis of primary incidices:
        if axis==0: 
            x=swapaxes(x,0,1)

        if axis>1:
            nextaxis = axis-1
        else:
            nextaxis = 0

        for ii in range(newshape[0]):
            #print 'x.shape>>',x.shape, 'newshape>>',newshape, 'x[ii].shape>>',x[ii].shape, 'ret[ii].shape>>',ret[ii].shape,'ii>>',ii
            ret[ii] = medianr(x[ii], nsigma=nsigma,niter=niter,finite=finite,\
                                 verbose=verbose, axis=nextaxis)
        return ret

        
    
def amedian(a, axis=None):
    """amedian(a, axis=None)

    Median the array over the given axis.  If the axis is None,
    median over all dimensions of the array.  

    Think of this as normal Numpy median, but preserving dimensionality.

    """
#    2008-07-24 14:12 IJC: Copied from
#         http://projects.scipy.org/scipy/numpy/ticket/558#comment:3,
#         with a few additions of my own (i.e., keeping the same
#         dimensionality).
#  2011-04-08 11:49 IJC: Moved to "analysis.py"

    sorted = array(a, subok=True, copy=True)

    if axis==None:
        return median(sorted.ravel())

    ash  = list(sorted.shape)
    ash[axis] = 1

    sorted = rollaxis(sorted, axis)
    
    sorted.sort(axis=0)
    
    index = int(sorted.shape[0]/2)
    
    if sorted.shape[0] % 2 == 1:
        returnval = sorted[index]
    else:
        returnval = (sorted[index-1]+sorted[index])/2.0
     
    return returnval.reshape(ash)



def polyfitr(x, y, N, s, fev=100, w=None, diag=False, clip='both', \
                 verbose=False, plotfit=False, plotall=False):
    """Matplotlib's polyfit with weights and sigma-clipping rejection.

    :DESCRIPTION:
      Do a best fit polynomial of order N of y to x.  Points whose fit
      residuals exeed s standard deviations are rejected and the fit is
      recalculated.  Return value is a vector of polynomial
      coefficients [pk ... p1 p0].

    :OPTIONS:
        w:   a set of weights for the data; uses CARSMath's weighted polynomial 
             fitting routine instead of numpy's standard polyfit.

        fev:  number of function evaluations to call before stopping

        'diag'nostic flag:  Return the tuple (p, chisq, n_iter)

        clip: 'both' -- remove outliers +/- 's' sigma from fit
              'above' -- remove outliers 's' sigma above fit
              'below' -- remove outliers 's' sigma below fit

    :REQUIREMENTS:
       :doc:`CARSMath`

    :NOTES:
       Iterates so long as n_newrejections>0 AND n_iter<fev. 


     """
    # 2008-10-01 13:01 IJC: Created & completed
    # 2009-10-01 10:23 IJC: 1 year later! Moved "import" statements within func.
    # 2009-10-22 14:01 IJC: Added 'clip' options for continuum fitting
    # 2009-12-08 15:35 IJC: Automatically clip all non-finite points
    # 2010-10-29 09:09 IJC: Moved pylab imports inside this function

    from CARSMath import polyfitw
    from numpy import polyfit, polyval, isfinite, ones
    from pylab import plot, legend, title

    xx = array(x, copy=True)
    yy = array(y, copy=True)
    noweights = (w==None)
    if noweights:
        ww = ones(xx.shape, float)
    else:
        ww = array(w, copy=True)

    ii = 0
    nrej = 1

    if noweights:
        goodind = isfinite(xx)*isfinite(yy)
    else:
        goodind = isfinite(xx)*isfinite(yy)*isfinite(ww)
    
    xx = xx[goodind]
    yy = yy[goodind]
    ww = ww[goodind]

    while (ii<fev and (nrej<>0)):
        if noweights:
            p = polyfit(xx,yy,N)
        else:
            p = polyfitw(xx,yy, ww, N)
            p = p[::-1]  # polyfitw uses reverse coefficient ordering
        residual = yy - polyval(p,xx)
        stdResidual = std(residual)
        if clip=='both':
            ind =  abs(residual) <= (s*stdResidual) 
        elif clip=='above':
            ind = residual < s*stdResidual
        elif clip=='below':
            ind = residual > -s*stdResidual
        else:
            ind = ones(residual.shape, bool)
        xx = xx[ind]
        yy = yy[ind]
        if (not noweights):
            ww = ww[ind]
        ii = ii + 1
        nrej = len(residual) - len(xx)
        if plotall:
            plot(x,y, '.', xx,yy, 'x', x, polyval(p, x), '--')
            legend(['data', 'fit data', 'fit'])
            title('Iter. #' + str(ii) + ' -- Close all windows to continue....')

        if verbose:
            print str(len(x)-len(xx)) + ' points rejected on iteration #' + str(ii)

    if (plotfit or plotall):
        plot(x,y, '.', xx,yy, 'x', x, polyval(p, x), '--')
        legend(['data', 'fit data', 'fit'])
        title('Close window to continue....')

    if diag:
        chisq = ( (residual)**2 / yy ).sum()
        p = (p, chisq, ii)

    return p


def spliner(x, y, k=3, sig=5, s=None, fev=100, w=None, clip='both', \
                 verbose=False, plotfit=False, plotall=False, diag=False):
    """Matplotlib's polyfit with sigma-clipping rejection.

    Do a scipy.interpolate.UnivariateSpline of order k of y to x.
     Points whose fit residuals exeed s standard deviations are
     rejected and the fit is recalculated.  Returned is a spline object.

     Iterates so long as n_newrejections>0 AND n_iter<fev. 

     :OPTIONAL INPUTS:
        err:   a set of errors for the data
        fev:  number of function evaluations to call before stopping
        'diag'nostic flag:  Return the tuple (p, chisq, n_iter)
        clip: 'both' -- remove outliers +/- 's' sigma from fit
              'above' -- remove outliers 's' sigma above fit
              'below' -- remove outliers 's' sigma below fit
     """
    # 2010-07-05 13:51 IJC: Adapted from polyfitr
    from numpy import polyfit, polyval, isfinite, ones
    from scipy import interpolate
    from pylab import plot, legend, title

    xx = array(x, copy=True)
    yy = array(y, copy=True)
    noweights = (w==None)
    if noweights:
        ww = ones(xx.shape, float)
    else:
        ww = array(w, copy=True)

    #ww = 1./err**2

    ii = 0
    nrej = 1

    goodind = isfinite(xx)*isfinite(yy)*isfinite(ww)
    
    #xx = xx[goodind]
    #yy = yy[goodind]
    #ww = ww[goodind]

    while (ii<fev and (nrej<>0)):
        spline = interpolate.UnivariateSpline(xx[goodind],yy[goodind],w=ww[goodind],s=s,k=k)
        residual = yy[goodind] - spline(xx[goodind])
        stdResidual = std(residual)
        #if verbose: print stdResidual
        if clip=='both':
            ind =  abs(residual) <= (sig*stdResidual) 
        elif clip=='above':
            ind = residual < sig*stdResidual
        elif clip=='below':
            ind = residual > -sig*stdResidual
        else:
            ind = ones(residual.shape, bool)

        goodind *= ind
        #xx = xx[ind]
        #yy = yy[ind]
        #ww = ww[ind]

        ii += 1
        nrej = len(residual) - len(xx)
        if plotall:
            plot(x,y, '.', xx[goodind],yy[goodind], 'x', x, spline(x), '--')
            legend(['data', 'fit data', 'fit'])
            title('Iter. #' + str(ii) + ' -- Close all windows to continue....')

        if verbose:
            print str(len(x)-len(xx[goodind])) + ' points rejected on iteration #' + str(ii)

    if (plotfit or plotall):
        plot(x,y, '.', xx[goodind],yy[goodind], 'x', x, spline(x), '--')
        legend(['data', 'fit data', 'fit'])
        title('Close window to continue....')

    if diag:
        chisq = ( (residual)**2 / yy )[goodind].sum()
        spline = (spline, chisq, ii, goodind)

    return spline



def neworder(N):
    """Generate a random order the integers [0, N-1] inclusive.

    """
    #2009-03-03 15:15 IJC: Created
    from numpy import random, arange, int
    from pylab import find

    neworder = arange(N)
    random.shuffle(neworder)

    #print N

    return  neworder
    

def im2(data1, data2, xlab='', ylab='', tit='', bar=False, newfig=True, \
            cl=None, x=[], y=[], fontsize=16):
    """Show two images; title will only be drawn for the top one."""
    from pylab import figure, subplot, colorbar, xlabel, ylabel, title, clim
    from nsdata import imshow

    if newfig:  figure()
    subplot(211)
    imshow(data1,x=x,y=y); 
    if clim<>None:   clim(cl)
    if bar:  colorbar()
    xlabel(xlab, fontsize=fontsize); ylabel(ylab, fontsize=fontsize)
    title(tit)

    subplot(212)
    imshow(data2,x=x,y=y); 
    if clim<>None:   clim(cl)
    if bar:  colorbar()
    xlabel(xlab, fontsize=fontsize); ylabel(ylab, fontsize=fontsize)

    return

def dumbconf(vec, sig, type='central', mid='mean', verbose=False):
    """
    Determine the displacement from the mean which encloses a fraction
    'sig' of the input vector.  This generates two-sided and lower-
    and upper-one-sided confidence limits.

    Try to do it using the bisector technique, until you reach some tolerance.

    :OPTIONAL INPUTS:
       type='both' -- 'upper', 'lower', or 'central' confidence limits
       mid='mean'  -- compute middle with mean or median

    :EXAMPLES:
       ::

           from numpy import random
           from analysis import dumbconf
           x = random.randn(10000)
           dumbconf(x, 0.683)    #  --->   1.0  (one-sigma)
           dumbconf(3*x, 0.954)  #  --->   6.0  (two-sigma)
           dumbconf(x+2, 0.997, type='lower')   #  --->   -0.74
           dumbconf(x+2, 0.997, type='upper')   #  --->    4.7

    
    Some typical values for a Normal distribution:
    one-sigma	0.6826895
    2 sigma	0.9544997
    3 sigma	0.9973002
    4 sigma	0.9999366
    5 sigma	0.9999994
"""
    #2009-03-25 22:47 IJC: A Better way to do it (using bisector technique)
    # 2009-03-26 15:37 IJC: Forget bisecting -- just sort.
    from numpy import sort, array

    vec = array(vec).copy()
    sig = array([sig]).ravel()

    #vec = sort(vec)
    N = len(vec)
    Ngoal = sig*N

    if mid=='mean':
        mid = vec.mean()
    elif mid=='median':
        mid = median(vec)
    else:
        try:
            mid = mid + 0.0
        except MidMethodException:
            print "mid must be median, mean, or numeric in value!"
            return -1
    
    if type =='central':
        vec2 = sort(abs(vec-mid))
    elif type=='upper':
        vec2 = sort(vec)
    elif type=='lower':
        vec2 = -sort(-vec)
    else:
        print "Invalid type -- must be central, upper, or lower"
        return -1
    
    ret = []
    for ii in range(len(sig)):
        ret.append(vec2[int(Ngoal[ii])])
    return ret


def error_dropoff(data):
    """   Calculates the dropoff in measurement uncertainty with increasing
    number of samples (a random and uncorrelated set of data should
    drop of as 1/sqrt[n] ).  

    E(0), the first returned element, always returns the uncertainty
    in a single measurement (i.e., the standard deviation).
    
    :EXAMPLE: Compute the errors on an array of 3 column vectors
       ::

           data = randn(1000,3)
           e = error_dropoff(data)
           plot(e[1], 1./sqrt(e[1]), '-', e[1], e[0], 'x')
           xlim([0,30])
           legend(['Theoretical [1/sqrt(N)]', 'Computed errors'])
    """
    # 2009-05-05 08:58 IJC: Adapted to Python from Matlab.
    # 2006/06/06 IJC: Made it work with arrays of column vectors.
    #                 Added the '--nomean' option.
    
    
# PARSE INPUTS
    data = array(data).copy()
    
#interval = max([1 round(extract_from_options('--interval=', 1, options))]);
    
    
    if len(data)==len(data.ravel()):
        data = data.ravel()
        data = data.reshape(len(data),1)
    
    nsets = data.shape[1]
    npts_vec = arange(data.shape[0]/2)+1.0
    errors = zeros((data.shape[0]/2, nsets))
 
# LOOP AND CALCULATE STUFF
    for ii in range(len(npts_vec)):
        npts = npts_vec[ii]  # number of points we average over
        nsamp = floor(data.shape[0]/npts) # number of subsamples
        dre = reshape(data[0:nsamp*npts,:], (npts, nsamp, nsets))
        error_values = std(dre.mean(1))
        errors[ii,:] = error_values
        
    return (errors, npts_vec)


def binarray(img,ndown, axis=None):
    """downsample a 2D image

    Takes a 1D vector or 2D array and reduce resolution by an integer factor
    "ndown".  This is done by binning the array -- i.e., subaveraging
    over square blocks of pixels of width "ndown"

    If keyword "axis" is None, bin over all axes.  Otherwise, bin over
    the single specified axis.

    :EXAMPLE:
       ::
       
           [img_ds]=binarray(img,ndown)        
    """
    # Renamed (and re-commented) by IJC 2007/01/31 from "downsample.m" to 
    #      "binarray.m"
    # 2009-05-06 14:15 IJC: Converted from Matlab to Python
    # 2009-09-27 14:37 IJC: Added 1D vector support.
    # 2010-09-27 16:08 IJC: ADded axis support.
    # 2010-10-26 10:55 IJC: Much simpler algorithm for 1D case
    # 2011-05-21 19:46 IJMC: Actually got axis support to work (2D
    #                        only). Fixed indexing bug (in x, y meshgrid)
    
    if ndown==1:
        return img

    img = array(img).copy()
    
    ndownx, ndowny = ndown, ndown
    if axis==0:
        ndownx = 1
    elif axis==1:
        ndowny = 1

    #if axis is None:
    if len(img.shape)==2:
        nrows, ncols = img.shape

        x = arange(0,ncols,ndownx, dtype=int)
        y = arange(0,nrows,ndowny, dtype=int)
        x,y = meshgrid(x,y)
        img_ds = 0*img[y,x]

        for kx in range(ndownx):
            for ky in range(ndowny):
                xs = (x+kx);
                xs = (xs*(xs<=ncols) + ncols*(xs>ncols));
                ys = (y+ky);
                ys = (ys*(ys<=nrows) + nrows*(ys>nrows));
                img_ds = img_ds + img[ys, xs];

    elif len(img.shape)==1:
        ncols = len(img)
        n_rep = ncols / ndown
        img_ds = array([img[ii::ndown][0:n_rep] for ii in range(ndown)]).mean(0)


    return img_ds


def fixval(arr, repval, retarr=False):
    """Fix non-finite values in a numpy array, and replace them with repval.

    :INPUT:
       arr -- numpy array, with values to be replaced.

       repval -- value to replace non-finite elements with

    :OPTIONAL INPUT:

       retarr -- if False, changes values in arr directly (more
       efficient).  if True, returns a fixed copy of the input array,
       which is left unchanged.

    :example:
     ::

       fixval(arr, -1)
       """
    # 2009-09-02 14:07 IJC: Created
    from numpy import isfinite, array
    from pylab import find
    if retarr:
        arr2 = arr.ravel().copy()
    else:
        arr2 = arr.ravel()

    badIndex = find((1-isfinite(arr2)))
    arr2[badIndex] = repval

    if retarr:
        return arr2.reshape(arr.shape)
    else:
        return


def removeoutliers(data, nsigma, remove='both', center='mean', niter=Inf, retind=False, verbose=False):
    """Strip outliers from a dataset, iterating until converged.

    :INPUT:
      data -- 1D numpy array.  data from which to remove outliers.

      nsigma -- positive number.  limit defining outliers: number of
                standard deviations from center of data.

    :OPTIONAL INPUTS:               
      remove -- ('min'|'max'|'both') respectively removes outliers
                 below, above, or on both sides of the limits set by
                 nsigma.

      center -- ('mean'|'median'|value) -- set central value, or
                 method to compute it.

      niter -- number of iterations before exit; defaults to Inf,
               which can occasionally result in empty arrays returned

      retind -- (bool) whether to return index of good values as
                second part of a 2-tuple.

    :EXAMPLE: 
       ::

           from numpy import hist, linspace, randn
           from analysis import removeoutliers
           data = randn(1000)
           hbins = linspace(-5,5,50)
           d2 = removeoutliers(data, 1.5, niter=1)
           hist(data, hbins)
           hist(d2, hbins)

       """
    # 2009-09-04 13:24 IJC: Created
    # 2009-09-24 17:34 IJC: Added 'retind' feature.  Tricky, but nice!
    # 2009-10-01 10:40 IJC: Added check for stdev==0
    # 2009-12-08 15:42 IJC: Added check for isfinite

    from numpy import median, ones, isfinite
    from pylab import find

    def getcen(data, method):
        "Get central value of a 1D array (helper function)"
        if method.__class__==str:
            if method=='median':
                cen = median(data)
            else:
                cen = data.mean()
        else:
            cen = method
        return cen

    def getgoodindex(data, nsigma, center, stdev, remove):
        "Get number of outliers (helper function!)"
        if stdev==0:
            distance = data*0.0
        else:
            distance = (data-center)/stdev
        if remove=='min':
            goodind = distance>-nsigma
        elif remove=='max':
            goodind = distance<nsigma
        else:
            goodind = abs(distance)<=nsigma
        return goodind

    data = data.ravel().copy()

    ndat0 = len(data)
    ndat = len(data)
    iter=0
    goodind = ones(data.shape,bool)
    goodind *= isfinite(data)
    while ((ndat0<>ndat) or (iter==0)) and (iter<niter) and (ndat>0) :
        ndat0 = len(data[goodind])
        cen = getcen(data[goodind], center)
        stdev = data[goodind].std()
        thisgoodind = getgoodindex(data[goodind], nsigma, cen, stdev, remove)
        goodind[find(goodind)] = thisgoodind
        if verbose:
            print "cen>>",cen
            print "std>>",stdev
        ndat = len(data[goodind])
        iter +=1
        if verbose:
            print ndat0, ndat
    if retind:
        ret = data[goodind], goodind
    else:
        ret = data[goodind]
    return ret
        

def xcorr2_qwik(img0, img1):
    """
    Perform quick 2D cross-correlation between two images.

    Images must be the same size.

    Computed via squashing the images along each dimension and
    computing 1D cross-correlations.
    """
    # 2009-12-17 10:13 IJC: Created.  Based on idea by J. Johnson.
    from numpy import zeros, max, min, sum
    from pylab import find

    im00 = img0.sum(0)
    im01 = img0.sum(1)
    im10 = img1.sum(0)
    im11 = img1.sum(1)
    n0 = len(im00)
    n1 = len(im01)
    noffsets0 = 2*n0-1
    noffsets1 = 2*n1-1
    corr0 = zeros(noffsets0,float)
    corr1 = zeros(noffsets1,float)

    for ii in range(noffsets0):
        firstind0 = max((ii-n0+1,0))
        lastind0 = min((n0,  ii+1))
        firstind1 = max((n0-ii-1,0))
        lastind1 = min((2*n0-ii-1,n0))
        corr0[ii] = sum(im00[firstind0:lastind0]*im10[firstind1:lastind1])

    for jj in range(noffsets1):
        firstind0 = max((jj-n0+1,0))
        lastind0 = min((n0,  jj+1))
        firstind1 = max((n0-jj-1,0))
        lastind1 = min((2*n0-jj-1,n0))
        corr1[jj] = sum(im01[firstind0:lastind0]*im11[firstind1:lastind1])

    ret = find([corr0==corr0.max()])-n0+1, find([corr1==corr1.max()])-n0+1
    return  ret
        
def lsq(x, z, w=None, retcov=False):
    """Do weighted least-squares fitting.  "x" input should be a tuple of
    1D vectors of equal lengths N; "z" should be a vector of length N,
    and "w" should be an N-vector or NxN array of weights.

    retcov -- bool.  Also return covariance matrix.

    Returns the tuple of (coef, coeferrs, {cov_matrix})"""
    # 2010-01-13 18:36 IJC: Created
    # 2010-02-08 13:04 IJC: Works for lists or tuples of x

    from numpy import vstack, dot, sqrt, isfinite,diag,ones,float, array
    from numpy.linalg import pinv


    if isinstance(x,tuple) or isinstance(x,list):
        Xmat = vstack(x).transpose()
    else:
        Xmat = x.reshape(len(x),1)

    if w==None:
        w = diag(ones(Xmat.shape[0],float))
    else:
        w = array(w,copy=True)
        if len(w.shape)<2:
            w = diag(w)

    goodind = isfinite(Xmat.sum(1))*isfinite(z)*isfinite(diag(w))
    Wmat = w[goodind][:,goodind]

    XtW = dot(Xmat[goodind,:].transpose(),Wmat)
    fitcoef = dot(dot(pinv(dot(XtW,Xmat[goodind,:])),XtW),z[goodind])
    
    covmat = (pinv(dot(XtW,Xmat[goodind,:])))
    efitcoef = sqrt(covmat)

    if retcov:
        return fitcoef, efitcoef, covmat
    else:
        return fitcoef, efitcoef


def medianfilter(data, filtersize, threshold=None,verbose=False):
    """ For now, we assume FILTERSIZE is odd, and that DATA is square!
    
    filt = medianfilter(data, filtersize)
    
    This is about the slowest way to do it, but it was easy to write.
    """
    # 2006/02/01 IJC at the Jet Propulsion Laboratory
    # 2010-02-18 13:52 IJC: Converted to python
    from numpy import zeros, median, abs, std

    print "Just use scipy.signal.medfilt !!!"
    print "Just use scipy.signal.medfilt !!!"
    print "Just use scipy.signal.medfilt !!!"

    if len(filtersize)<1:
        print 'medianfilter2 requires that filtersize be a 1- or 2-element vector'
        return -1
    elif len(filtersize)==1:
        filtersize = [filtersize[0], filtersize[0]]
    else:
        filtersize = filtersize[0:2]

    npix = data.shape[0]
    npiy = data.shape[1]
    bigsize = npix+2*(filtersize[0]-1)
    bigdata = zeros((bigsize,bigsize),float)
    ind = filtersize[0]-1
    if ind==0:
        bigdata = data
    else:
        bigdata[ind:(bigsize-ind), ind:(bigsize-ind)] = data


    # FOR NOW, WE ASSUME FILTERSIZE IS ODD!!
    # AND THAT DATA IS SQUARE!
    niter_x = npix + (filtersize[0]-1)
    niter_y = npiy + (filtersize[1]-1)
    filt = zeros((niter_x,niter_y), float)

    for ii in range(niter_x):
        for jj in range(niter_y):
            if verbose>1:
                print "ii,jj>>",ii,jj
            if filtersize[0]==1:
                indi = 1
            else:
                indi = filtersize[0]-1
            if filtersize[1]==1:
                indj = 1
            else:
                indj = filtersize[1]-1
            select = bigdata[ii:(ii+indi),jj:(jj+indj)].ravel()
            select = select[isfinite(select)]
            #residualSelection = abs(select - median(select))

            if verbose: 
                print "select.shape>>",select.shape
                print "threshold>>",threshold

            if threshold is not None:
                if threshold >= 0: # raw threshold
                    doFilter = abs(bigdata[ii,jj]-median(select))/std(select)>=threshold
                elif threshold<0:  # remove outliers before applying threshold
                    npts_init = len(select)
                    select = removeoutliers(select, abs(threshold), center='median')
                    npts_final = len(select)
                    if verbose>1:
                        print "threshold=",threshold,", removed %i points" % (npts_init-npts_final)
                    
                    doFilter = abs(bigdata[ii,jj]-median(select))/std(select)>=threshold                    
            else:  # filter everything; threshold not set.
                doFilter = True

            if verbose:
                print "doFilter?>>",doFilter
            if verbose>1:
                print "select>>",select

            if doFilter: 
                newval =  median( select )
            else:
                newval = bigdata[ii,jj]

            if verbose>1:
                print "newval>>",newval

            filt[ii,jj] = newval

    print filt.shape, [(filtersize[0]-1)/2,niter_x-(filtersize[0]-1)/2,(filtersize[0]-1)/2,niter_y-(filtersize[0]-1)/2]
    return filt[(filtersize[0]-1)/2:niter_x-(filtersize[0]-1)/2,(filtersize[0]-1)/2:niter_y-(filtersize[0]-1)/2]


def wmean(a, w, axis=None):
    """wmean(a, w, axis=None)

    Perform a weighted mean along the specified axis.

    :INPUTS:
      a : sequence or Numpy array
        data for which weighted mean is computed

      w : sequence or Numpy array
        weights of data -- e.g., 1./sigma^2

    :SEE ALSO:  :func:`wstd`
    """
    # 2008-07-30 12:44 IJC: Created this from ...
    # 2012-02-28 20:31 IJMC: Added a bit of documentation


    newdata    = array(a, subok=True, copy=True)
    newweights = array(w, subok=True, copy=True)

    if axis==None:
        newdata    = newdata.ravel()
        newweights = newweights.ravel()
        axis = 0

    ash  = list(newdata.shape)
    wsh  = list(newweights.shape)

    nsh = list(ash)
    nsh[axis] = 1

    if ash<>wsh:
        warn('Data and weight must be arrays of same shape.')
        return []
    
    wsum = newweights.sum(axis=axis).reshape(nsh) 
    
    weightedmean = (a * newweights).sum(axis=axis).reshape(nsh) / wsum

    return weightedmean


def wstd(a, w, axis=None):
    """wstd(a, w, axis=None)

    Perform a weighted standard deviation along the specified axis.
    If axis=None, then the weighted standard deviation of the entire
    array is computed.

    Note that this computes the _sample_ standard deviation;
    Numpy/Scipy computes the _population_ standard deviation, which is
    greater by a factor sqrt(N/N-1).  This effect is small for large
    datasets.

    :SEE ALSO:  :func:`wmean`

    Taken from http://en.wikipedia.org/wiki/Weighted_standard_deviation
    """
# 2008-07-30 12:44 IJC: Created this from 


    newdata    = array(a, subok=True, copy=True)
    newweights = array(w, subok=True, copy=True)

    if axis==None:
        newdata    = newdata.ravel()
        newweights = newweights.ravel()
        axis = 0

    ash  = list(newdata.shape)
    wsh  = list(newweights.shape)

    nsh = list(ash)
    nsh[axis] = 1

    if ash<>wsh:
        warn('Data and weight must be arrays of same shape.')
        return []
    
    wsum = newweights.sum(axis=axis).reshape(nsh) 
    omega = 1.0 * wsum / (wsum**2 - (newweights**2).sum(axis=axis).reshape(nsh))
    
    weightedstd = omega * (newweights * (newdata-wmean(newdata, newweights, axis=axis))**2 ).sum(axis=axis).reshape(nsh)

    return sqrt(weightedstd)


def fmin(func, x0, args=(), kw=dict(),  xtol=1e-4, ftol=1e-4, maxiter=None, maxfun=None,
         full_output=0, disp=1, retall=0, callback=None, zdelt = 0.00025, nonzdelt = 0.05, 
         holdfixed=None):
    """Minimize a function using the downhill simplex algorithm -- now with KEYWORDS.

    :Parameters:

      func : callable func(x,*args)
          The objective function to be minimized.
      x0 : ndarray
          Initial guess.
      args : tuple
          Extra arguments passed to func, i.e. ``f(x,*args)``.
      callback : callable
          Called after each iteration, as callback(xk), where xk is the
          current parameter vector.

    :Returns: (xopt, {fopt, iter, funcalls, warnflag})

      xopt : ndarray
          Parameter that minimizes function.
      fopt : float
          Value of function at minimum: ``fopt = func(xopt)``.
      iter : int
          Number of iterations performed.
      funcalls : int
          Number of function calls made.
      warnflag : int
          1 : Maximum number of function evaluations made.
          2 : Maximum number of iterations reached.
      allvecs : list
          Solution at each iteration.

    *Other Parameters*:

      xtol : float
          Relative error in xopt acceptable for convergence.
      ftol : number
          Relative error in func(xopt) acceptable for convergence.
      maxiter : int
          Maximum number of iterations to perform.
      maxfun : number
          Maximum number of function evaluations to make.
      full_output : bool
          Set to True if fval and warnflag outputs are desired.
      disp : bool
          Set to True to print convergence messages.
      retall : bool
          Set to True to return list of solutions at each iteration.
      zdelt : number
          Set the initial stepsize for x0[k] equal to zero
      nonzdelt : number
          Set the initial stepsize for x0[k] nonzero
      holdfixed : sequence
          Indices of x0 to hold fixed (e.g., [1, 2, 4])


    :TBD:  gprior : tuple or sequence of tuples
          Set a gaussian prior on the indicated parameter, such that
          chisq += ((x0[p] - val)/unc_val)**2, where the parameters
          are defined by the tuple gprior=(param, val, unc_val)

    :Notes:

        Uses a Nelder-Mead simplex algorithm to find the minimum of
        function of one or more variables.

    """
    # 2011-04-13 14:26 IJMC: Adding Keyword option
    # 2011-05-11 10:48 IJMC: Added the zdelt and nonzdelt options
    # 2011-05-30 15:36 IJMC: Added the holdfixed option

    def wrap_function(function, args, **kw):
        ncalls = [0]
        def function_wrapper(x):
            ncalls[0] += 1
            return function(x, *args, **kw)
        return ncalls, function_wrapper

    # Set up holdfixed arrays
    if holdfixed is not None:
        holdfixed = np.array(holdfixed)
        #x0[holdfixed] = x0[holdfixed]
        holdsome = True
    else:
        holdsome = False
        #holdfixed = np.zeros(params.size, dtype=bool)
    
    #if holdsome:
    #    print "holdfixed>>", holdfixed

    fcalls, func = wrap_function(func, args, **kw)
    x0 = np.asfarray(x0).flatten()
    xoriginal = x0.copy()
    N = len(x0)
    rank = len(x0.shape)
    if not -1 < rank < 2:
        raise ValueError, "Initial guess must be a scalar or rank-1 sequence."
    if maxiter is None:
        maxiter = N * 200
    if maxfun is None:
        maxfun = N * 200

    rho = 1; chi = 2; psi = 0.5; sigma = 0.5;
    one2np1 = range(1,N+1)

    if rank == 0:
        sim = np.zeros((N+1,), dtype=x0.dtype)
    else:
        sim = np.zeros((N+1,N), dtype=x0.dtype)
    fsim = np.zeros((N+1,), float)
    sim[0] = x0
    if retall:
        allvecs = [sim[0]]
    #print func.__name__
    #print x0
    fsim[0] = func(x0)
    for k in range(0,N):
        y = np.array(x0,copy=True)
        if y[k] != 0:
            y[k] = (1+nonzdelt)*y[k]
        else:
            y[k] = zdelt
        if holdsome and k in holdfixed:
            y[k] = xoriginal[k]
        sim[k+1] = y
        f = func(y)
        fsim[k+1] = f

    ind = np.argsort(fsim)
    fsim = np.take(fsim,ind,0)
    # sort so sim[0,:] has the lowest function value
    sim = np.take(sim,ind,0)

    iterations = 1

    while (fcalls[0] < maxfun and iterations < maxiter):
        ### IJC Edit to understand fmin!
        ##print 'xtol>> ' + str(max(np.ravel(abs(sim[1:]-sim[0])))) + ' > ' + str(xtol)
        ##print 'ftol>> ' + str(max(abs(fsim[0]-fsim[1:]))) + ' > ' + str(ftol)
        if (max(np.ravel(abs(sim[1:]-sim[0]))) <= xtol \
            and max(abs(fsim[0]-fsim[1:])) <= ftol):
            break

        xbar = np.add.reduce(sim[:-1],0) / N
        xr = (1+rho)*xbar - rho*sim[-1]
        if holdsome:
            xr[holdfixed] = xoriginal[holdfixed]
        fxr = func(xr)
        doshrink = 0

        if fxr < fsim[0]:
            xe = (1+rho*chi)*xbar - rho*chi*sim[-1]
            if holdsome:
                xe[holdfixed] = xoriginal[holdfixed]
            fxe = func(xe)

            if fxe < fxr:
                sim[-1] = xe
                fsim[-1] = fxe
            else:
                sim[-1] = xr
                fsim[-1] = fxr
        else: # fsim[0] <= fxr
            if fxr < fsim[-2]:
                sim[-1] = xr
                fsim[-1] = fxr
            else: # fxr >= fsim[-2]
                # Perform contraction
                if fxr < fsim[-1]:
                    xc = (1+psi*rho)*xbar - psi*rho*sim[-1]
                    if holdsome:
                        xc[holdfixed] = xoriginal[holdfixed]
                    fxc = func(xc)

                    if fxc <= fxr:
                        sim[-1] = xc
                        fsim[-1] = fxc
                    else:
                        doshrink=1
                else:
                    # Perform an inside contraction
                    xcc = (1-psi)*xbar + psi*sim[-1]
                    if holdsome:
                        xcc[holdfixed] = xoriginal[holdfixed]
                    fxcc = func(xcc)

                    if fxcc < fsim[-1]:
                        sim[-1] = xcc
                        fsim[-1] = fxcc
                    else:
                        doshrink = 1

                if doshrink:
                    for j in one2np1:
                        sim[j] = sim[0] + sigma*(sim[j] - sim[0])
                        if holdsome:
                            sim[j, holdfixed] = xoriginal[holdfixed]
                        fsim[j] = func(sim[j])

        ind = np.argsort(fsim)
        sim = np.take(sim,ind,0)
        fsim = np.take(fsim,ind,0)
        if callback is not None:
            callback(sim[0])
        iterations += 1
        if retall:
            allvecs.append(sim[0])

    x = sim[0]
    fval = min(fsim)
    warnflag = 0

    if fcalls[0] >= maxfun:
        warnflag = 1
        if disp:
            print "Warning: Maximum number of function evaluations has "\
                  "been exceeded."
    elif iterations >= maxiter:
        warnflag = 2
        if disp:
            print "Warning: Maximum number of iterations has been exceeded"
    else:
        if disp:
            print "Optimization terminated successfully."
            print "         Current function value: %f" % fval
            print "         Iterations: %d" % iterations
            print "         Function evaluations: %d" % fcalls[0]


    if full_output:
        retlist = x, fval, iterations, fcalls[0], warnflag
        if retall:
            retlist += (allvecs,)
    else:
        retlist = x
        if retall:
            retlist = (x, allvecs)

    return retlist

def fmin_powell(func, x0, args=(), kw=dict(), xtol=1e-4, ftol=1e-4, maxiter=None,
                maxfun=None, full_output=0, disp=1, retall=0, callback=None,
                direc=None, holdfixed=None):
    """Minimize a function using modified Powell's method -- now with KEYWORDS.

    :Parameters:

      func : callable f(x,*args)
          Objective function to be minimized.
      x0 : ndarray
          Initial guess.
      args : tuple
          Eextra arguments passed to func.
      kw : dict
          Keyword arguments passed to func.
      callback : callable
          An optional user-supplied function, called after each
          iteration.  Called as ``callback(xk)``, where ``xk`` is the
          current parameter vector.
      direc : ndarray
          Initial direction set.

    :Returns: (xopt, {fopt, xi, direc, iter, funcalls, warnflag}, {allvecs})

        xopt : ndarray
            Parameter which minimizes `func`.
        fopt : number
            Value of function at minimum: ``fopt = func(xopt)``.
        direc : ndarray
            Current direction set.
        iter : int
            Number of iterations.
        funcalls : int
            Number of function calls made.
        warnflag : int
            Integer warning flag:
                1 : Maximum number of function evaluations.
                2 : Maximum number of iterations.
        allvecs : list
            List of solutions at each iteration.

    *Other Parameters*:

      xtol : float
          Line-search error tolerance.
      ftol : float
          Relative error in ``func(xopt)`` acceptable for convergence.
      maxiter : int
          Maximum number of iterations to perform.
      maxfun : int
          Maximum number of function evaluations to make.
      full_output : bool
          If True, fopt, xi, direc, iter, funcalls, and
          warnflag are returned.
      disp : bool
          If True, print convergence messages.
      retall : bool
          If True, return a list of the solution at each iteration.


    :Notes:

        Uses a modification of Powell's method to find the minimum of
        a function of N variables.

    """
    # 2010-07-01 11:17 IJC: Added keyword option

    from scipy import optimize
    from numpy import asarray, eye, pi, squeeze

    def wrap_function(function, args, **kw):
        ncalls = [0]
        def function_wrapper(x):
            ncalls[0] += 1
            return function(x, *args, **kw)
        return ncalls, function_wrapper

    def _linesearch_powell(func, p, xi, tol=1e-3):
        """Line-search algorithm using fminbound.

        Find the minimium of the function ``func(x0+ alpha*direc)``.

        """
        def myfunc(alpha):
            return func(p + alpha * xi)
        alpha_min, fret, iter, num = optimize.brent(myfunc, full_output=1, tol=tol)
        xi = alpha_min*xi
        return squeeze(fret), p+xi, xi


    # Set up holdfixed arrays
    if holdfixed is not None:
        holdfixed = np.array(holdfixed)
        #x0[holdfixed] = x0[holdfixed]
        holdsome = True
    else:
        holdsome = False
        #holdfixed = np.zeros(params.size, dtype=bool)

    # we need to use a mutable object here that we can update in the
    # wrapper function
    fcalls, func = wrap_function(func, args, **kw)
    x = asarray(x0).flatten()
    xoriginal = x.copy()
    if retall:
        allvecs = [x]
    N = len(x)
    rank = len(x.shape)
    if not -1 < rank < 2:
        raise ValueError, "Initial guess must be a scalar or rank-1 sequence."
    if maxiter is None:
        maxiter = N * 1000
    if maxfun is None:
        maxfun = N * 1000


    if direc is None:
        direc = eye(N, dtype=float)
    else:
        direc = asarray(direc, dtype=float)

    fval = squeeze(func(x))
    x1 = x.copy()
    iter = 0;
    ilist = range(N)
    while True:
        fx = fval
        bigind = 0
        delta = 0.0
        for i in ilist:
            direc1 = direc[i]
            fx2 = fval
            if (not holdsome) or (i not in holdfixed):
                fval, x, direc1 = _linesearch_powell(func, x, direc1, tol=xtol*100)
                if (fx2 - fval) > delta:
                    delta = fx2 - fval
                    bigind = i
        iter += 1
        if callback is not None:
            callback(x)
        if retall:
            allvecs.append(x)
        if (2.0*(fx - fval) <= ftol*(abs(fx)+abs(fval))+1e-20): break
        if fcalls[0] >= maxfun: break
        if iter >= maxiter: break

        # Construct the extrapolated point
        direc1 = x - x1
        x2 = 2*x - x1
        if holdsome:
            x2[holdfixed] = xoriginal[holdfixed]
        x1 = x.copy()
        fx2 = squeeze(func(x2))

        if (fx > fx2):
            t = 2.0*(fx+fx2-2.0*fval)
            temp = (fx-fval-delta)
            t *= temp*temp
            temp = fx-fx2
            t -= delta*temp*temp
            if t < 0.0:
                fval, x, direc1 = _linesearch_powell(func, x, direc1,
                                                     tol=xtol*100)
                if holdsome:
                    x[holdfixed] = xoriginal[holdfixed]
                direc[bigind] = direc[-1]
                direc[-1] = direc1

    warnflag = 0
    if fcalls[0] >= maxfun:
        warnflag = 1
        if disp:
            print "Warning: Maximum number of function evaluations has "\
                  "been exceeded."
    elif iter >= maxiter:
        warnflag = 2
        if disp:
            print "Warning: Maximum number of iterations has been exceeded"
    else:
        if disp:
            print "Optimization terminated successfully."
            print "         Current function value: %f" % fval
            print "         Iterations: %d" % iter
            print "         Function evaluations: %d" % fcalls[0]

    x = squeeze(x)

    if full_output:
        retlist = x, fval, direc, iter, fcalls[0], warnflag
        if retall:
            retlist += (allvecs,)
    else:
        retlist = x
        if retall:
            retlist = (x, allvecs)

    return retlist


def gaussian2d(p, x, y):
    """ Compute a 2D gaussian distribution at the points x, y.

        p is a four- or five-component array, list, or tuple:

        y =  [p4 +] p0/(2*pi*p1**2) * exp(-((x-p2)**2+(y-p3)) / (2*p1**2))

        p[0] -- Area of the gaussian
        p[1] -- one-sigma dispersion
        p[2] -- x-central offset
        p[3] -- y-central offset
        p[4] -- optional constant, vertical offset

        SEE ALSO:  :func:`egaussian`
        """
    #2010-06-08 20:00 IJC: Created
    
    x = array(x, dtype=float).copy()
    y = array(y, dtype=float).copy()
    p = array(p).copy()

    if len(p)==4:
        p = concatenate((p, [0]))

    z = p[4] + p[0]/(2*pi*p[1]**2) * exp(-((x-p[2])**2 + (y-p[3])**2) / (2*p[1]**2))
    
    return z

def gaussian2d_ellip(p, x, y):
    """ Compute a 2D elliptical gaussian distribution at the points x, y.

        p is a 5-, 6-, or 7-component sequence, defined as:

            p[0] -- Amplitude (Area of the function)

            p[1] -- x-dispersion

            p[2] -- y-dispersion

            p[3] -- x-central offset

            p[4] -- y-central offset

            p[5] -- optional rotation angle (radians)

            p[6] -- optional constant, vertical offset

        X, Y are gridded data from :func:`numpy.meshgrid`

        First define:
          x' = (x - p[3]) cos p[5] - (y - p[4]) sin p[5]

          y' = (x - p[3]) sin p[5] + (y - p[4]) cos p[5]

        Then calculate:
          U = (x' / p[1])**2 + (y' / p[2])**2

          z =  p[6] + p0/(2*pi*p1*p2) * exp(-U / 2)

        SEE ALSO:  :func:`gaussian2d`, :func:`lorentzian2d`
        """
    #2012-02-11 18:06 IJMC: Created from IDL GAUSS2DFIT
    
    x = array(x, dtype=float).copy()
    y = array(y, dtype=float).copy()
    p = array(p).copy()

    if len(p)==5:
        p = concatenate((p, [0, 0]))
    elif len(p)==6:
        p = concatenate((p, [0]))

    xp = (x - p[3]) * np.cos(p[5]) - (y - p[4]) * np.sin(p[5])
    yp = (x - p[3]) * np.sin(p[5]) + (y - p[4]) * np.cos(p[5])

    z = p[6] + p[0]/(2*pi*p[1]*p[2]) * exp(-0.5 * ((xp / p[1])**2 + (yp / p[2])**2))
    
    return z

def lorentzian2d(p, x, y):
    """ Compute a 2D Lorentzian distribution at the points x, y.

        p is a 5-, 6-, or 7--component sequence:

          z = (x-p3) ** 2 / p1 ** 2 + (y-p4) ** 2 / p2 ** 2  [ + (x-p3) * (y-p4) * p5 ]
            lorentz = p0 / (1.0 + z)   [ + p6]

        p[0] -- Amplitude (Area of the function)
        p[1] -- x-dispersion
        p[2] -- y-dispersion
        p[3] -- x-central offset
        p[4] -- y-central offset
        p[5] -- optional ellipticitity parameter
        p[6] -- optional constant, vertical offset

        SEE ALSO:  :func:`gaussian2d`
        """
    #2012-02-04 11:38 IJMC: Created
    
    x = array(x, dtype=float).copy()
    y = array(y, dtype=float).copy()
    p = array(p).copy()

    if len(p)==5:
        p = concatenate((p, [0, 0]))
    elif len(p)==6:
        p = concatenate((p, [0]))

    z = ((x - p[3]) / p[1])**2 + ((y - p[4]) / p[2])**2 + p[5] * (x - p[3]) * (y - p[4])
    
    return p[6] + p[0]/(1. + z)

def egaussian2d(p,x,y,z,w=None):
    """ Return the error associated with a 2D gaussian fit, using gaussian2d.

    w is an array of weights, typically 1./sigma**2"""
    # 2010-06-08 20:02 IJC: Created
    from numpy import ones, array
    x = array(x, dtype=float).copy()
    if w is None:
        w = ones(w.shape,float)
    z0 = gaussian2d(p,x,y)
    return (((z-z0)*w)**2).sum()


def getblocks(vec):
    """Return start and end indices for consecutive sequences of
    integer-valued indices.

    :Example:
        ::

             import analysis as an
             vec = range(5) +range(10,14) + range(22,39)
             starts,ends = an.getblocks(vec)
             print zip(starts,ends)
             """
    # 2010-08-18 17:01 IJC: Created
    # 2010-11-15 23:26 IJC: Added numpy imports

    from numpy import sort, diff, nonzero

    vec = sort(vec)
    starts = [vec[0]]
    ends = []
    dvec = diff(vec)
    start_inds = nonzero(dvec>1)[0]
    for ind in start_inds:
        starts.append(vec[ind+1])
        ends.append(vec[ind])

    ends.append(vec[-1])

    return starts, ends

def snr(data, axis=None, nsigma=None):
    """
    Compute the quantity:
    data.mean(axis=axis)/data.std(axis=axis)
    
    for the specified data array/vector along the specified axis.
    
    Output will be a scalar (axis is None) or numpy array, as
    appropriate.
    """
    # 2010-09-02 08:10 IJC: Created
    # 2011-12-16 14:51 IJMC: Added optional nsigma flag

    data = array(data)
    if nsigma is None:
        ret = data.mean(axis=axis)/data.std(axis=axis)
    else:
        ret = meanr(data, axis=axis, nsigma=nsigma) / \
            stdr(data, axis=axis, nsigma=nsigma)
    return  ret

def pad(inp, npix_rows, npix_cols=None):
    """Pads input matrix to size specified. 
    ::

       out = pad(in, npix)     
       out = pad(in, npix_rows, npix_cols);  # alternate usage
      
    Written by J. Green @ JPL; converted to Python by I. Crossfield"""
    #2008-10-18 12:50 IJC: Converted from Matlab function
    # 2010-10-29 09:35 IJC: Moved from nsdata.py to analysis.py

    from numpy import imag, zeros, complex128

    inp = array(inp, copy=True)
    if len(inp.shape)==0:
        inp = inp.reshape((1,1))
    elif len(inp.shape)==1:
        inp = inp.reshape((1, len(inp)))

    if npix_cols==None:
        npix_cols = npix_rows
        
    if (imag(inp)**2).sum()==0:
        out = zeros((npix_rows, npix_cols))
    else:
        out = zeros((npix_rows, npix_cols), complex128)
    
    nrows, ncols  = inp.shape

    ixc = floor(ncols/2 + 1);
    iyc = floor(nrows/2 + 1);

    oxc = floor(npix_cols/2 + 1);
    oyc = floor(npix_rows/2 + 1);

    dx = npix_cols-ncols;
    dy = npix_rows-nrows;

    if dx<=0:
        ix1 = ixc - floor(npix_cols/2);
        ix2 = ix1 + npix_cols - 1;
        ox1 = 1;
        ox2 = npix_cols;
    else:
        ix1 = 1;
        ix2 = ncols;
        ox1 = oxc - floor(ncols/2);
        ox2 = ox1 + ncols - 1;

   
    if dy<=0:
        iy1 = iyc - floor(npix_rows/2);
        iy2 = iy1 + npix_rows - 1;
        oy1 = 1;
        oy2 = npix_rows;
    else:
        iy1 = 1;
        iy2 = nrows;
        oy1 = oyc - floor(nrows/2);
        oy2 = oy1 + nrows - 1;

    out[ oy1-1:oy2, ox1-1:ox2]  = inp[ iy1-1:iy2, ix1-1:ix2];

# Uncomment for testing
#    print inp
#    print ixc, iyc, iy1, iy2, ix1, ix2 
#    print oxc, oyc, oy1, oy2, ox1, ox2 

    return out


def fftfilter1d(vec, bandwidth, retfilter=False):
    """ Apply a hard-edged low-pass filter to an input vector.

    :INPUTS:

       vec -- sequence -- 1D vector, assumed to be evenly sampled

       bandwidth -- integer -- size of the filter passband in cycles
                              per signal duration.  f <= bandwidth is
                              passed through; f > bandwidth is
                              rejected.

       retfilter -- bool -- if True, return the tuple (filtered_vec, filter)

    :OUTPUT:
       
       Lopass-filtered version of vec

    :NOTE:

       Assumes the input is real-valued.
       """
    # 2011-03-11 14:15 IJC: Created

    from numpy import real, fft, floor, ceil

    vec = array(vec, copy=True)

    # Errorchecking
    if len(vec.shape)<>1:
        print "Input array must be 1D -- try using .ravel()"
        return -1


    npts = vec.size

    filter = concatenate((zeros(floor(npts/2.) - bandwidth), 
                          ones(bandwidth * 2 + 1),
                          zeros(ceil(npts/2.) - bandwidth - 1)))

    ret = real(fft.ifft(fft.ifftshift( fft.fftshift(fft.fft(vec)) * filter )))

    if retfilter==True:
        ret = [ret, filter]

    return ret


def stdres(data, bins=None):
    """Compute the standard deviation in the residuals of a data
    series after average-binning by specified amounts.

    :INPUTS:
      data - 1D numpy array
         Data to analyze.

      bins - sequence
         Factors by which to bin down.  If None, use 1:sqrt(data.size)
         """
    # 2011-06-16 16:50 IJMC: Created

    ndata = data.size
    if bins is None:
        bins = arange(1, sqrt(int(ndata)))

    nout = len(bins)

    ret = zeros(nout, dtype=float)
    for jj, binfactor in enumerate(bins):
        if binfactor>0:
            ret[jj] = binarray(data, binfactor).std()
        else:
            ret[jj] = 0

    return ret
         

def allanvariance(data, dt=1):
    """Compute the Allan variance on a set of regularly-sampled data (1D).

       If the time between samples is dt and there are N total
       samples, the returned variance spectrum will have frequency
       indices from 1/dt to (N-1)/dt."""
    # 2008-07-30 10:20 IJC: Created
    # 2011-04-08 11:48 IJC: Moved to analysis.py

    newdata = array(data, subok=True, copy=True)
    dsh = newdata.shape

    newdata = newdata.ravel()

    nsh = newdata.shape

    alvar = zeros(nsh[0]-1, float)

    for lag in range(1, nsh[0]):
        alvar[lag-1]  = mean( (newdata[0:-lag] - newdata[lag:])**2 )

    return (alvar*0.5)
    
def trueanomaly(ecc, eanom=None, manom=None):
    """Calculate (Keplerian, orbital) true anomaly.

    One optional input must be given.

    :INPUT:
       ecc -- scalar.  orbital eccentricity.

    :OPTIONAL_INPUTS:
       eanom -- scalar or Numpy array.  Eccentric anomaly.  See
               :func:`eccentricanomaly`

       manom -- scalar or sequence.  Mean anomaly, equal to 
                2*pi*(t - t0)/period
    """
    # 2011-04-22 14:35 IJC: Created

    if manom is not None:
        eanom = eccentricanomaly(ecc, manom=manom)

    if eanom is not None:
        ret = 2. * np.arctan(  np.sqrt((1+ecc)/(1.-ecc)) * np.tan(eanom/2.)  )
    else:
        ret = None

    return ret

def eccentricanomaly(ecc, manom=None, tanom=None):
    """Calculate (Keplerian, orbital) eccentric anomaly.

    One optional input must be given.

    :INPUT:
       ecc -- scalar.  orbital eccentricity.

    :OPTIONAL_INPUTS:

       manom -- scalar or sequence.  Mean anomaly, equal to 
                2*pi*(t - t0)/period

       tanom -- scalar or Numpy array.  True anomaly.  See
               :func:`trueanomaly`.
    """
    # 2011-04-22 14:35 IJC: Created

    ret = None
    if manom is not None:
        if not hasattr(manom, '__iter__'):
            mwasscalar = True
            manom = [manom]
        else:
            mwasscalar = False
        # Solve Kepler's equation for each element of mean anomaly:
        e = np.zeros(len(manom))
        for ii,element in enumerate(manom):
            def kep(e): return element - e + ecc*sin(e)
            e[ii] = optimize.fsolve(kep, 0.)

        if mwasscalar:
            e = e[0]
        ret = e
    
    elif tanom is not None:
        ret = 2. * np.arctan(np.tan(0.5 * tanom) / \
                                 np.sqrt((1. + ecc) / (1. - ecc)))
        
    else:
        ret = None

    return ret
    

def gaussian(p, x):
    """ Compute a gaussian distribution at the points x.

        p is a three- or four-component array, list, or tuple:

        y =  [p3 +] p0/(p1*sqrt(2pi)) * exp(-(x-p2)**2 / (2*p1**2))

        p[0] -- Area of the gaussian
        p[1] -- one-sigma dispersion
        p[2] -- central offset (mean location)
        p[3] -- optional constant, vertical offset

        NOTE: FWHM = 2*sqrt(2*ln(2)) * p1  ~ 2.3548*p1

        SEE ALSO:  :func:`egaussian`"""
    #2008-09-11 15:11 IJC: Created for LINEPROFILE
    # 2011-05-18 11:46 IJC: Moved to analysis.
    
    x = array(x, dtype=float).copy()
    p = array(p).copy()

    if len(p)==3:
        p = concatenate((p, [0]))

    y = p[3] + p[0]/(p[1]*sqrt(2*pi)) * exp(-(x-p[2])**2 / (2*p[1]**2))
    
    return y

def egaussian(p, x, y, e=None):
    """ Compute the deviation between the values y and the gaussian defined by p, x:

    p is a three- or four-component array, list, or tuple.

    Returns:   y - p3 - p0/(p1*sqrt(2pi)) * exp(-(x-p2)**2 / (2*p1**2))

    if an error array, e (typ. one-sigma) is entered, the returned value is divided by e.

    SEE ALSO:  :func:`gaussian`"""
    # 2008-09-11 15:19 IJC: Created
    # 2009-09-02 15:20 IJC: Added weighted case
    # 2011-05-18 11:46 IJMC: Moved to analysis.
    from numpy import ones

    if e==None:
        e=ones(x.shape)
    fixval(e,y.max()*1e10)

    z = (y - gaussian(p, x))/e
    fixval(z,0)

    return z


def generic_mcmc(*arg, **kw):
    """Run a Markov Chain Monte Carlo (Metropolis-Hastings algorithm)
    on an arbitrary function.

    :INPUTS:
      EITHER:
        func : function to generate model.  
          First argument must be "params;" subsequent arguments are
          passed in via the "args" keyword

        params : 1D sequence
          parameters to be fit

        stepsize :  1D or 2D array
                If 1D: 1-sigma change in parameter per iteration
                If 2D: covariance matrix for parameter changes.

        z : 1D array
                    Contains dependent data (to be modeled)

        sigma : 1D array
                    Contains standard deviation (errors) of "z" data

        numit : int
                Number of iterations to perform

      OR:
        (allparams, (arg1, arg2, ...), numit)

        where allparams is a concatenated list of parameters for each
        of several functions, and the arg_i are tuples of (func_i,
        stepsize_i, z_i, sigma_i). In this case the keyword 'args'
        must also be a tuple of sequences, one for each function to be
        MCMC'ed.

    :OPTIONAL_INPUTS:    
        args : 1D sequence
                Second, third, etc.... arguments to "func"

        nstep : int
                Saves every "nth" step of the chain

        posdef : None, 'all', or sequences of indices.
                Which elements should be restricted to positive definite?
                If indices, it should be of the form (e.g.): [0, 1, 4]

        holdfixed : None, or sequences of indices.
                    Which elements should be held fixed in the analysis?
                    If indices, it should be of the form (e.g.): [0, 1, 4]

        jointpars : None, or sequence of 2-tuples.
                    Only for simultaneous multi-function fitting.  For
                    each pair of values passed, we set the parameters
                    values so: allparams[pair[1]] = allparams[pair[0]]

    :OUTPUTS:
        allparams : 2D array
                Contains all parameters at each step
        bestp : 1D array
                Contains best paramters as determined by lowest Chi^2
        numaccept: int
                Number of accepted steps
        chisq: 1D array
                Chi-squared value at each step
    
    :REFERENCES:
        Numerical Recipes, 3rd Edition (Section 15.8)
        Wikipedia

    """
    # 2011-06-07 07:50 IJMC: Created from various other MCMC codes,
    #                        eventually deriving from K. Stevenson's
    #                        sample code.
    # 2011-06-14 09:48 IJMC: Allow for covariance matrix pass-in to stepsize
    # 2011-06-27 17:39 IJMC: Now link joint parameters for initial chisq. 
    # 2011-09-16 13:31 IJMC: Fixed bug for nextp when nfits==1
    # 2011-11-02 22:08 IJMC: Now cast numit as an int

    import numpy as np
    
    # Parse keywords/optional inputs:
    defaults = dict(args=(), nstep=1, posdef=None, holdfixed=None, \
                        jointpars=None, verbose=False)
    for key in defaults:
        if (not kw.has_key(key)):
            kw[key] = defaults[key]

    args = kw['args']
    nstep = kw['nstep']
    posdef = kw['posdef']
    holdfixed = kw['holdfixed']
    jointpars = kw['jointpars']
    verbose = kw['verbose']

    # Parse inputs:
    if len(arg)==6:
        func, params, stepsize, z, sigma, numit = arg
        stepsize = np.array(stepsize, copy=True)
        weights = 1./sigma**2
        nfits = 1

    elif len(arg)==3:
        params, allargs, numit = arg[0:3]
        nfits = len(allargs)
        funcs = []
        stepsizes = []
        zs = []
        multiargs = []
        multiweights = []
        npars = []
        for ii, these_args in enumerate(allargs):
            funcs.append(these_args[0])
            stepsizes.append(np.array(these_args[1], copy=True))
            zs.append(these_args[2])
            multiweights.append(1./these_args[3]**2)
            multiargs.append(args[ii])
            npars.append(stepsizes[-1].shape[0])
    else:
        print "Must pass either 3 or 6 parameters as input."
        print "You passed %i." % len(arg)
        return -1
    
    

    #Initial setup
    numaccept  = 0
    numit = int(numit)
    nout = numit/nstep
    bestp      = np.copy(params)
    original_params = np.copy(params)
    allparams  = np.zeros((len(params), nout))
    allchi     = np.zeros(nout,float)

    # Set indicated parameters to be positive definite:
    if posdef=='all':
        params = np.abs(params)
        posdef = np.arange(params.size)
    elif posdef is not None:
        posdef = np.array(posdef)
        params[posdef] = np.abs(params[posdef])
    else:
        posdef = np.zeros(params.size, dtype=bool)

    # Set indicated parameters to be held fixed:
    if holdfixed is not None:
        holdfixed = np.array(holdfixed)
        params[holdfixed] = np.abs(params[holdfixed])
    else:
        holdfixed = np.zeros(params.size, dtype=bool)

    if verbose:
        print params[posdef]

    # Set joint parameters:
    if jointpars is not None:
        for jp in jointpars:
            params[jp[1]] = params[jp[0]]


    #Calc chi-squared for model using current params
    if nfits==1:
        zmodel     = func(params, *args)
        currchisq  = (((zmodel - z)**2)*weights).ravel().sum()
        bestchisq  = currchisq
    else:
        tempchisq = 0
        for ii in range(nfits):
            i0 = sum(npars[0:ii])
            i1 = i0 + npars[ii]
            this_zmodel = funcs[ii](params[i0:i1], *multiargs[ii])
            thischisq = (((this_zmodel - zs[ii])**2) * multiweights[ii]).ravel().sum()
            tempchisq += thischisq
        currchisq = tempchisq
        bestchisq = currchisq

    if verbose:
        print currchisq

    #Run Metropolis-Hastings Monte Carlo algorithm 'numit' times
    for j in range(numit):
        #Take step in random direction for adjustable parameters
        if nfits==1:
            if len(stepsize.shape)==1:
                nextp    = np.array([np.random.normal(params,stepsize)]).ravel()
            else:
                nextp = np.random.multivariate_normal(params, stepsize)
        else:
            nextstep = np.zeros(len(params), dtype=float)
            for ii in range(nfits):
                i0 = sum(npars[0:ii])
                i1 = i0 + npars[ii]
                if len(stepsizes[ii].shape)==1:
                    nextstep[i0:i1] = np.random.normal([0]*npars[ii], stepsizes[ii])
                else:                
                    nextstep[i0:i1] = np.random.multivariate_normal([0]*npars[ii], stepsizes[ii])
            nextp = params + nextstep

        # Constrain the desired parameters:
        nextp[posdef] = np.abs(nextp[posdef])
        nextp[holdfixed] = original_params[holdfixed]
        if jointpars is not None:
            for jp in jointpars:
                nextp[jp[1]] = nextp[jp[0]]
                #print nextp[jp[1]], nextp[jp[0]], jp

        #COMPUTE NEXT CHI SQUARED AND ACCEPTANCE VALUES
        if nfits==1:
            zmodel     = func(nextp, *args)
            nextchisq  = (((zmodel - z)**2)*weights).ravel().sum() 
        else:
            tempchisq = 0
            for ii in range(nfits):
                i0 = sum(npars[0:ii])
                i1 = i0 + npars[ii]
                this_zmodel = funcs[ii](nextp[i0:i1], *multiargs[ii])
                thischisq = (((this_zmodel - zs[ii])**2) * multiweights[ii]).ravel().sum()
                tempchisq += thischisq
            nextchisq = tempchisq

        if verbose:
            print nextchisq
            print nextp == original_params

        accept = np.exp(0.5 * (currchisq - nextchisq))

        if (accept >= 1) or (np.random.uniform(0, 1) <= accept):
                #Accept step
                numaccept += 1
                params  = np.copy(nextp)
                currchisq  = nextchisq
        if (currchisq < bestchisq):
                        #New best fit
                        bestp     = np.copy(params)
                        bestchisq = currchisq

        if (j%nstep)==0:
            allparams[:, j/nstep] = params
            allchi[j/nstep] = currchisq

    return allparams, bestp, numaccept, allchi



def scale_mcmc_stepsize(accept, func, params, stepsize, z, sigma, numit=1000, scales=[0.1, 0.3, 1., 3., 10.], args=(), nstep=1, posdef=None, holdfixed=None, retall=False, jointpars=None):
    """Run :func:`generic_mcmc` and scale the input stepsize to match
    the desired input acceptance rate.

    :INPUTS:
       mostly the same as for :func:`generic_mcmc`, but also with:

       accept : float between 0 and 1
          desired acceptance rate; typically between 0.15 - 0.5.

       scales : sequence of floats
          test scaling parameters; measure for these, then interpolate to get 'accept'

       retall : bool
          if True, return tuple (scalefactor, acceptances, scales).
          Otherwise return only the scaler 'scalefactor.'

    :REQUIREMENTS:
       :doc:`pylab` (for :func:`pylab.interp`)
          """
    # 2011-06-13 16:06 IJMC: Created

    from pylab import interp

    stepsize = array(stepsize, copy=True)

    nfactors = len(scales)
    mcmc_accept = []
    for factor in scales:
        out = generic_mcmc(func, params, stepsize/factor, z, sigma, numit, args=args, nstep=nstep, posdef=posdef, holdfixed=holdfixed, jointpars=jointpars)
        mcmc_accept.append(1.0 * out[2]/numit)

    final_factor = interp(accept, mcmc_accept, scales)

    if retall:
        ret = (final_factor, mcmc_accept, scales)
    else:
        ret = final_factor

    return ret


def unityslope(slope, ttt):
    return 1. + slope*(ttt - ttt.mean())

def travisplanet(p):
    """Generate a line of text for Travis Barman's planet table.

    INPUT: a planet object from :func:`getobj`.
    """
    # 2012-02-14 14:52 IJMC: Created

    vals1 = (p.name.replace(' ','').replace('b', ''), p.msini, p.umsini, p.r, p.r-p.ur, p.r+p.ur)
    vals2 = (p.per, p.a, p.mstar, p.umstar, p.rstar, p.urstar, p.teff, p.uteff, p.fe, p.ufe)
    nvals = len(vals1)+len(vals2)
    fstr = '%s' + ' '*(10-len(vals1)) + '  %1.2f'*5 + '  %1.6f  %1.4f' + '  %1.2f'*4 + '  %i'*2 + '  %1.2f'*2

    return  fstr  % (vals1 + vals2 ) 
