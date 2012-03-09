"""Python script to test photometry tasks.  Maybe I can develop it
into something more reusable someday.
"""

# 2008-12-21 18:27 IJC: Created

from numpy import array, sign, pi, nan, arange
import pdb
#try:
#    import psyco
#    psyco.full()
#except ImportError:
#    print 'Psyco not installed, the program will just run slower'

class phot:
    def __init__(self, **kw):
        """Generate a photometry object w/keywords and kw+'str' descriptive keywords

        Inputs cannot be numpy arrays"""
        # 2009-09-14 10:07 IJC: Created
        keylist = ['time', 'phot', 'ephot', 'bg', 'ebg', 'aper', 'position', 'eposition', \
            'ditherposition', 'object', 'filename', 'exptime', 'ra', 'dec']
        defaults = dict(photunit=None)
        for key in keylist:
            defaults[key]=None
            defaults[key+'str']=None
            keylist = keylist + [key+'str']
        for keyword in defaults:
            if (not kw.has_key(keyword)):
                kw[keyword] = defaults[keyword]

        for k in kw.keys():
            if kw[k].__class__==str:
                exec('self.'+k+'="'+str(kw[k]))+'"'
            else:
                exec('self.'+k+'='+str(kw[k]))
        self.comment = "Created by IJC's spitzer.phot"
        return 

 
def estbg(im, mask=None, bins=None, plotalot=False, rout=(3,200), badval=nan):
    """Estimate the background value of a masked image via histogram fitting.

    INPUTS:
      im -- numpy array.  Input image.

    OPTIONAL INPUTS:
      mask -- numpy array. logical mask, False/0 in regions to ignore
      bins -- sequence.  edges of bins to pass to HIST
      plotalot -- bool.  Plot the histogram and fit.
      rout -- 2-tuple of (nsigma, niter) for analysis.removeoutliers.
              Set to (Inf, 0) to not cut any outliers.
      badval -- value returned when things go wrong.

    OUTPUT:
      b, s_b -- tuple of (background, error on background) from gaussian fit.
                 Note that the error is analagous to the standard deviation on the mean

    COMMENTS:
      The fit parameters appear to be robust across a fairly wide range of bin sizes.  """
    # 2009-09-02 17:13 IJC: Created!
    # 2009-09-04 15:07 IJC: Added RemoveOutliers option. Use only non-empty bins in fit.
    # 2009-09-08 15:32 IJC: Error returned is now divided by sqrt(N) for SDOM
    # 2009-11-03 00:16 IJC: Improved guess for gaussian dispersion
    # 2011-05-18 11:47 IJMC: Moved (e)gaussian imports to analysis.
    # 2012-01-01 21:04 IJMC: Added badval option

    from numpy import histogram, mean, median, sqrt, linspace, isfinite, ones,std
    from pylab import find
    from scipy import optimize
    from analysis import removeoutliers, egaussian, gaussian
    if plotalot:
        from pylab import figure, errorbar, plot, colorbar, title, hist, mean, std
        #from analysis import imshow

    def gaussianChiSquared(guess, x, y, err):
        return (egaussian(guess, x, y, e=err)**2).sum()


    if mask==None:
        mask = ones(im.shape)
    dat = im.ravel()[find(mask<>0)]
    if plotalot:
        figure(); plot(im.ravel()); plot(dat)
        print mean(dat), std(dat), rout[0]*std(dat)
        print len(dat), (abs(dat-mean(dat))<(rout[0]*std(dat))).sum()
        figure(); plot(dat-mean(dat)); 
        plot([0,len(dat)], [rout[0]*std(dat),rout[0]*std(dat)],'--k')
        plot([0,len(dat)], [-rout[0]*std(dat),-rout[0]*std(dat)],'--k')
    dat = removeoutliers(dat, rout[0], remove='both', center='mean', niter=rout[1], verbose=plotalot)
    ndat = len(dat)

    if ndat==0:
        print "No data to work with!"
        return (badval, badval)
    if bins==None:
        if plotalot: print "no bins entered!"
        datmean = dat.mean()
        datstd = dat.std()
        bins = linspace(datmean-3*datstd,datmean+3*datstd, len(dat)/10)
    if plotalot: 
        print "dat.mean, dat.std>>" + str((dat.mean(), dat.std()))

    binwidth = mean(bins[1::]-bins[:-1])
    bincenter = 0.5*(bins[1::]+bins[:-1])
    datIndex = (dat>=bins.min()) * (dat<=bins.max())

    #if plotalot:
    #    figure(); hout = hist(dat[datIndex],bins)
    #else:
    hout = histogram(dat[datIndex], bins,new=True)

    gy = hout[0]
    erry = sqrt(gy)
    usableIndex = gy>0

    guess = [gy.sum()*binwidth, std(dat[datIndex]), median(dat[datIndex])]
    out = optimize.fmin(gaussianChiSquared, guess, \
                            args=(bincenter[usableIndex],gy[usableIndex], erry[usableIndex]), \
                            disp=plotalot)

    if plotalot:
        from pylab import figure, errorbar, plot, colorbar, title
        from nsdata import imshow
        print 'guess>>',guess
        print 'fit>>',out
        figure()
        imshow(im*mask); colorbar()
        figure()
        errorbar(bincenter, gy, erry, fmt='ob')
        plot(bincenter, gaussian(out, bincenter),'-r', linewidth=2)
        title('Mean: %f, Std. Dev.: %f' % (out[2], out[1]))

    return out[2], out[1]/sqrt(ndat)


def makemask(x,y,params, shape='circle'):
    """Generate a binary (logical) mask with given shape and location.

    INPUTS:
         x      = x-coodinate system (made with meshgrid)
         y      = y-coodinate system (made with meshgrid)
      params:
        shape='circle':
           params(1)  = x-offset
           params(2)  = y-offset
           params(3)  = x-diameter
           params(4)  = OPTIONAL y-diameter
        shape='quad':
           params: list of quadrants to include in the mask.  The
              midpoint for determining quadrants will be
              mid = (xmax+xmin)/2.  Quadrants are:
              0:   x<midX  and  y<midY
              1:   x<midX  and  y>=midY
              2:   x>=midX  and  y<midY
              3:   x>=midX  and  y>=midY

    OPTIONAL INPUTS:
        shape=:  desired mask shape.  Currently only 'circle' is valid.

    OUTPUTS:
      mask   = NxM grided representation of specified mask
                  where NxM are the size of the x,y input meshgrids

    EXAMPLE:
        x=arange(32); y=arange(64)
        xx,yy = meshgrid(x,y)
        m = makemask(xx,yy,(24, 53, 10))
        m[53,24]   # ----> True
        """
    # 2009-09-02 10:47 IJC: Created.  Adapted from Joe Green's MakeEllipse.m
    # 2009-09-27 13:37 IJC: Added 'quad' quadrant mask option.
    from numpy import zeros, bool

    if not x.shape==y.shape:
        print "x,y meshgrid coordinates must be the same shape! Exiting."
        return -1
    if shape=='circle':
        if len(params)<3:
            print "Must give at least 3 parameters to mask."
            return -1
        x0 = params[0]
        y0 = params[1]
        xdia =params[2]
        if len(params)==3:
            ydia = xdia
        else:
            ydia = params[3]
        mask = (  (((x-x0)/(xdia/2.))**2 + ((y-y0)/(ydia/2.))**2) < 1 )
    elif shape=='quad':
        midx = (x.max()+x.min())/2.
        midy = (y.max()+y.min())/2.
        mask = zeros(x.shape, bool)
        for ii in range(len(params)):
            if params[ii]==0:
                mask += (x<midx) * (y<midy)
            elif params[ii]==1:
                mask += (x<midx) * (y>=midy)
            elif params[ii]==2:
                mask += (x>=midx) * (y<midy)
            elif params[ii]==3:
                mask += (x>=midx) * (y>=midy)

    return mask

def centroid(im, mask=None, w=None, x=None, y=None):
    """Compute the centroid of an image with a specified binary mask projected upon it.
    
    INPUT:
      im -- image array
      mask -- binary mask, 0 in ignored regions and 1 in desired regions
      w is typically 1.0/u**2, where u is the uncertainty on im
      x,y are those generated by meshgrid.

    OUTPUT:
      (x0,y0) tuple of centroid location"""
    from numpy import ones, arange, meshgrid
    # 2009-09-02 13:35 IJC: Created
    if mask==None:
        mask = ones(im.shape)
    if w==None:
        w = ones(im.shape)
    if not (im.shape==mask.shape and im.shape==w.shape):
        print "Image, mask, and weights must have same shape! Exiting."
        return -1
    if x==None or y==None:
        xx = arange(im.shape[1])
        yy = arange(im.shape[0])
        x,y = meshgrid(xx,yy)
    x0 = (x*im*mask*w).sum()/(im*mask*w).sum()
    y0 = (y*im*mask*w).sum()/(im*mask*w).sum()

    return (x0,y0)


def lmcextinct(ra, dec, **kw):
    """Use the Zaritsky & Harris (ZH) map to get A_V extinction in the LMC.

    INPUT:
       ra  -- decimal degrees of right ascension
       dec -- decimal degrees of declination

    OPTIONAL INPUT:
       method='griddata'  /'nearest'  -- interpolation method
       map='both' /'hot'/'cool'       -- Which ZH map to use
       hot='lmc_hotav.fits'           -- filename of hot map
       cool='lmc_coolav.fits'         -- filename of cool map
       null=0.0                       -- "no data" value in map
       verbose=False / True           -- print out some spam and eggs

    EXAMPLE:
       ra = hms('05:51:56.6')
       dec = dms('-66:25:08.5')
       lmcextinct(ra, dec)
       
    If their map returns null, interpolate from the surroundings.

    Note that these calculations are definitely _not_ optimized.

    SEE ALSO:
       hms, dms"""

    # 2009-02-15 22:11 IJC: Created and tested.

    import pyfits
    from matplotlib.mlab import griddata
    from pylab import meshgrid, arange, array, sqrt, cos, sin, arctan2, arcsin

    ra = array([ra]).copy().ravel()
    dec = array([dec]).copy().ravel()

    defaults = dict(method='griddata', map='hot', hot='lmc_hotav_clean2.fits', 
                    cool='lmc_coolav.fits', null=0.0, verbose=False)

    for key in defaults:
        if (not kw.has_key(key)):
            kw[key] = defaults[key]
    
    if kw['map']=='both':
        kw['map'] = 'hot'
        ret1 = lmcextinct(ra, dec, **kw)
        kw['map'] = 'cool'
        ret2 = lmcextinct(ra, dec, **kw)
        return (ret1, ret2)
    else:
        map = kw[kw['map']]

    verbose = kw['verbose']

    avmap = pyfits.getdata(map)
    avhdr = pyfits.getheader(map)
    
    nx = avhdr['naxis1']
    ny = avhdr['naxis2']
    cx = avhdr['crpix1']
    cy = avhdr['crpix2']
    x0 = avhdr['crval1']
    y0 = avhdr['crval2']
    dx = avhdr['cd1_1']
    dy = avhdr['cd2_2']

    xx,yy = meshgrid(arange(nx),arange(ny))
    
    goodind = (avmap<>kw['null'])

    # I don't know how the following WCS calculation works, but it
    # does -- thank you, Calabretta & Greisen 2002!
    d2r = pi/180.

    phi = arctan2( xx-cx+1, yy-cy+1) + pi
    theta = arctan2(1./(d2r*dy), sqrt((yy-cy)**2 + (xx-cx)**2))
    mapdec = arcsin(sin(theta)*sin(y0*d2r) - cos(theta)*cos(phi)*cos(y0*d2r))/d2r
    mapra = arcsin(cos(theta) * sin(phi) / cos(mapdec*d2r))/d2r + x0

    if kw['method']=='griddata':
        ragood = mapra[goodind]
        decgood = mapdec[goodind]
        avgood = avmap[goodind]

    if verbose:
        print 'ra.shape>>' + str(ra.shape)
    
    # TBD: Vectorize this calculation; an interpolative solution
    # should exist.
    avlist = []
    for ii in range(len(ra)):
        if verbose:
            print 'ra, dec>>' + str((ra[ii], dec[ii]))
        if kw['method']=='nearest':
            distmap = (mapra-ra[ii])**2 + (mapdec-dec[ii])**2
            # If multiple values are equally near, average them:
            val = avmap[distmap==distmap.min()].mean()
            avlist.append(val)
        elif kw['method']=='griddata':
            avlist.append(griddata(ragood, decgood, avgood, array([ra[ii]]), array([dec[ii]])))
        else:
            print "Invalid method specified!"
            avlist.append(-1)
        
 
    return avlist
    


def readogle(filename, **kw):
    """ Read in an OGLE-II photometry map file and output the data.

    Returns a 6-tuple with each element an array:
       0  -- RA.  Strings of Right Ascension values.
       1  -- DEC.  Strings of Declination values
       2  -- x_ref.  OGLE x-coordinate, pixels
       3  -- y_ref.  OGLE y-coordinate, pixels
       4  -- v_mag.  OGLE V magnitude
       5  -- i_mag.  OGLE I magnitude

    If you like, use HMS and DMS to convert the RA and DEC values returned."""
    
    # 2008-12-21 18:53 IJC: Created

    f = open(filename, 'r')
    raw = f.readlines()
    f.close()

    nstars = len(raw)

    raw2 = array([line.split() for line in raw])
    ra = raw2[:,1]
    dec = raw2[:,2]
    xref = raw2[:,3]
    yref = raw2[:,4]
    vmag = raw2[:,5]
    imag = raw2[:,7]
    
    xref = [map(float, [x]) for x in xref]
    yref = [map(float, [y]) for y in yref]
    vmag = [map(float, [v]) for v in vmag]
    imag = [map(float, [i]) for i in imag]

    return (ra, dec, xref, yref, vmag, imag)


def hms(d, delim=':'):
    """Convert hours, minutes, seconds to decimal degrees, and back.

    EXAMPLES:

    hms('15:15:32.8')
    hms([7, 49])
    hms(18.235097)
    
    Also works for negative values.

    SEE ALSO:  dms
    """
    # 2008-12-22 00:40 IJC: Created
    # 2009-02-16 14:07 IJC: Works with spaced or colon-ed delimiters
    from numpy import sign

    if d.__class__==str or hasattr(d, '__iter__'):   # must be HMS
        if d.__class__==str:
            d = d.split(delim)
            if len(d)==1:
                d = d[0].split(' ')
            if (len(d)==1) and (d.find('h')>-1):
                d.replace('h',delim)
                d.replace('m',delim)
                d.replace('s','')
                d = d.split(delim)
        s = sign(float(d[0]))
        if s==0:  s=1
        degval = float(d[0])*15.0
        if len(d)>=2:
            degval = degval + s*float(d[1])/4.0
        if len(d)==3:
            degval = degval + s*float(d[2])/240.0
        return degval
    else:    # must be decimal degrees
        hour = int(d/15.0)
        d = abs(d)
        min = int((d-hour*15.0)*4.0)
        sec = (d-hour*15.0-min/4.0)*240.0
        return [hour, min, sec]


def dms(d, delim=':'):
    """Convert degrees, minutes, seconds to decimal degrees, and back.

    EXAMPLES:

    dms('150:15:32.8')
    dms([7, 49])
    dms(18.235097)
    
    Also works for negative values.

    SEE ALSO:  hms
    """
    # 2008-12-22 00:40 IJC: Created
    # 2009-02-16 14:07 IJC: Works with spaced or colon-ed delimiters
    from numpy import sign

    if d.__class__==str or hasattr(d, '__iter__'):   # must be HMS
        if d.__class__==str:
            d = d.split(delim)
            if len(d)==1:
                d = d[0].split(' ')
        s = sign(float(d[0]))
        if s==0:  s=1
        degval = float(d[0])
        if len(d)>=2:
            degval = degval + s*float(d[1])/60.0
        if len(d)==3:
            degval = degval + s*float(d[2])/3600.0
        return degval
    else:    # must be decimal degrees
        if d<0:  
            sgn = -1
        else:
            sgn = +1
        d = abs(d)
        deg = int(d)
        min = int((d-deg)*60.0)
        sec = (d-deg-min/60.0)*3600.0
        return [sgn*deg, min, sec]


def hess(color, mag, binsize, **kw):
    """Compute a hess diagram (surface-density CMD) on photometry data.

    INPUT:
       color  
       mag
       binsize -- width of bins, in magnitudes

    OPTIONAL INPUT:
       cbin=  -- set the centers of the color bins
       mbin=  -- set the centers of the magnitude bins

    OUTPUT:
       A 3-tuple consisting of:
         Cbin -- the centers of the color bins
         Mbin -- the centers of the magnitude bins
         Hess -- The Hess diagram array"""

    # cbin = out[0]
    # mbin = out[1]
    # imshow(out[2])
    # yticks(range(0, len(mbin), 4), mbin[range(0,len(mbin),4)])
    # xticks(range(0, len(cbin), 4), cbin[range(0,len(cbin),4)])
    # ylim([ylim()[1], ylim()[0]])

    # 2009-02-08 23:01 IJC: Created, on a whim, for LMC data (of course)
    # 2009-02-21 15:45 IJC: Updated with cbin, mbin options

    from numpy import arange, zeros

    defaults = dict(mbin=None, cbin=None, verbose=False)

    for key in defaults:
        if (not kw.has_key(key)):
            kw[key] = defaults[key]

    if kw['mbin']==None:
        mbin = arange(mag.min(), mag.max(), binsize)
    else:
        mbin = array(kw['mbin']).copy()
    if kw['cbin']==None:
        cbin = arange(color.min(), color.max(), binsize)
    else:
        cbin = array(kw['cbin']).copy()

    hess = zeros((len(mbin), len(cbin)), float)
    for ii in range(len(cbin)):
        cindex = (color<(cbin[ii]+binsize/2)) * (color>(cbin[ii]-binsize/2)) 
        for jj in range(len(mbin)):
            index = cindex * (mag<(mbin[jj]+binsize/2)) * (mag>(mbin[jj]-binsize/2)) 
            hess[jj,ii] = index.sum()


    return (cbin, mbin, hess)


def snr(mag=20, itime=1., read=24.5, sky=8.43, npix=24., zero=26.44, dark=0.0):
    """Calculate SNR of a photometric CCD observation of given parameters.

    The default keyword parameters are for the CTIO Blanco MOSAIC
    imager, operating in V band.
    """
    # 2009-02-20 14:40 IJC: Initiated
    
    star = itime * 10**(0.4*(zero-mag))
    noise = npix * (itime*(sky+dark)+read**2)

    return star * (star+noise)**-0.5
    
    

#def sens(texp, mag, npix, tdead=15., rn=8., csky=12., zp=22.8, z=1.0,h=1280., d=1.0):
def sens(texp, mag, lin=3.8e5, tdead=15., rn=8., csky=12., zp=22.8, z=1.0,h=1280., d=1.0):
    """ Calculate sensitiviy for various exposure times.

    texp -- 1D array of times
    mag -- target magnitude
    npix -- number of pixels read out
    
    lin -- linearity limit of the detector (ADU)
    tdead -- detector dead time (reset + readout) between integrations
    rn -- read noise, in ADU
    csky -- sky background per second, in ADU
    zp -- photometric zero point (magnitude for 1 ADU/sec)
    z -- airmass
    h -- observatory altitude, in meters
    d -- telescope mirror diameter, in meters

    TBD:  Flat-field noise!


    For now, I impose an artificial floor to the number of pixels
    needed -- it won't spread it out over less than 8.
    """

    # 2009-03-13 11:48 IJC: Written based on http://arxiv.org/abs/0903.2139

    from numpy import exp, log10, sqrt, floor

    texp = array(texp).copy()

    ctarg = 10**(0.4*(zp-mag))

    Ntarg = ctarg * texp
    npix = floor(Ntarg / (lin-csky))  # spread out the light
    pixind = npix<8
    npix[pixind] = 8

    Nsky = csky * texp * npix
    
    sScint = 0.004 * d**(-2./3.) * z**1.75 * exp(-h/8e3) * (2*texp)**-0.5
    
    sTarg = -2.5*log10((Ntarg - sqrt(Ntarg)) / Ntarg)
    sSky = -2.5*log10((Ntarg - sqrt(Nsky)) / Ntarg)
    sRON = -2.5*log10((Ntarg - rn*sqrt(npix)) / Ntarg)

    sTotal = sqrt(sScint**2 + sSky**2 + sRON**2)

    snrpertime = sTotal * sqrt(texp+tdead)

    return sTotal, snrpertime, npix

def subreg(fn, center=None, dim=None, verbose=False):
    """Load a subsection of a list of 2D FITS files.

    INPUTS:
       fn -- (str) filename or (list) list of filenames to load
       center -- (2-tuple) (x0,y0) center of region to load.
       dim -- (2-tuple) (dx,dy) sizes of region to load. 

    OUTPUT: 
       region -- (array) (N x dx x dy)-sized stack of subregions

    NOTE: center can also be a list of tuples (1 per file), but dim cannot."""
    # 2009-09-24 11:20 IJC: Created
    print "Function subreg is deprecated; please use subreg2 (which is much faster!)"

    import pyfits
    from numpy import zeros, nan, concatenate, int, float, max, min

    if not hasattr(fn, '__iter__'):
        fn = [fn]

    nfiles = len(fn)
    if (center is None) and (dim is None):
        try:
            temp = pyfits.getdata(fn[0], ignore_missing_end=True)
        except:
            temp = pyfits.getdata(fn[0])
        dosubreg = False
        dim = list(temp.shape)
        center = dim[0]/2, dim[1]/2

    if not hasattr(center[0], '__iter__'):
        center = [center]*nfiles
    stack = zeros((0, dim[0], dim[1]), float)
    dim = [dim]*nfiles

    if verbose: print "fn, center, dim",(fn,center,dim)

    iter = 0
    for file, cen, sz in zip(fn, center, dim):
        if verbose: print "file, cen, sz", (file,cen,sz)
        try:
            temp = pyfits.getdata(file, ignore_missing_end=True)
        except:
            temp = pyfits.getdata(file)
        if verbose: print "file "+file+" shape is: "+str(temp.shape)
        if len(temp.shape)<>2:
            print "File %s did not return a 2D FITS file! Using nan." % file
            subreg = zeros((1,dx,dy))+nan
        else:
            temp = temp.reshape((1,)+temp.shape)
            xstartval = max([0,cen[0] - (sz[0]-1)/2.]).astype(int)
            xendval = min([temp.shape[1], 1+cen[0] + (sz[0]-1)/2.]).astype(int)
            ystartval = max([0,cen[1] - (sz[1]-1)/2.]).astype(int)
            yendval = min([temp.shape[2], 1+cen[1] + (sz[1]-1)/2.]).astype(int)
            if verbose: print "subregion limits: ", xstartval,xendval,ystartval,yendval
            subreg = temp[:,xstartval:xendval,ystartval:yendval]
        if file==fn[0] and stack.shape[1::]<>subreg.shape:
            stack = zeros((0,)+subreg.shape[1::],float)
        if verbose: 
            print "stack.shape>>", stack.shape
            print "subreg.shape>>", subreg.shape
        stack=concatenate((stack,subreg),0)
        iter += 1
            
    return stack
    

def subreg2(fn, center=None, dim=None, verbose=False):
    """Load a subsection of a list of 2D FITS files.

    INPUTS:
       fn -- (str) filename or (list) list of filenames to load
       center -- (2-tuple) (x0,y0) center of region to load.
       dim -- (2-tuple) (dx,dy) sizes of region to load. 

    OUTPUT: 
       region -- (array) (N x dx x dy)-sized stack of subregions

    NOTE: center can also be a list of tuples (1 per file), but dim cannot."""
    # 2009-09-24 11:20 IJC: Created
    # 2011-12-15 15:53 IJMC: Updated to use memmap. MUCH faster!
    import pyfits
    from numpy import zeros, nan, concatenate, int, float, max, min

    if not hasattr(fn, '__iter__'):
        fn = [fn]

    nfiles = len(fn)
    if (center is None) and (dim is None):
        try:
            temp = pyfits.getdata(fn[0], ignore_missing_end=True)
        except:
            temp = pyfits.getdata(fn[0])
        dosubreg = False
        dim = list(temp.shape)
        center = dim[0]/2, dim[1]/2

    if not hasattr(center[0], '__iter__'):
        center = [center]*nfiles
    stack = zeros((0, dim[0], dim[1]), float)
    dim = [dim]*nfiles

    if verbose: print "fn, center, dim",(fn,center,dim)

    iter = 0
    for file, cen, sz in zip(fn, center, dim):
        if verbose: print "file, cen, sz", (file,cen,sz)
        try:
            f = pyfits.open(file, ignore_missing_end=True, memmap=True)
        except:
            f = pyfits.open(file, memmap=True)
        #if verbose: print "file "+file+" shape is: "+str(temp.shape)
        #if len(temp.shape)<>2:
        #    print "File %s did not return a 2D FITS file! Using nan." % file
        #    subreg = zeros((1,dx,dy))+nan
        #else:
        if iter==0:
            temp = f[0].data#.reshape([1]+dim[0])
            temp = temp.reshape((1,) + temp.shape)
            xstartval = max([0,cen[0] - (sz[0]-1)/2.]).astype(int)
            xendval = min([temp.shape[1], 1+cen[0] + (sz[0]-1)/2.]).astype(int)
            ystartval = max([0,cen[1] - (sz[1]-1)/2.]).astype(int)
            yendval = min([temp.shape[2], 1+cen[1] + (sz[1]-1)/2.]).astype(int)
            if verbose: print "subregion limits: ", xstartval,xendval,ystartval,yendval
            subreg = temp[0,xstartval:xendval,ystartval:yendval]
            stack = zeros((nfiles,) + subreg.shape, float)
        else:
            this_dy =  center[0][0] - center[iter][0]
            this_dx =  center[0][1] - center[iter][1]
            try:
                subreg = f[0].data[xstartval+this_dx:xendval+this_dx, \
                                       ystartval+this_dy:yendval+this_dy]
            except IndexError:
                print "Likely invalid center position for Frame %i (%s): %s." % \
                    (iter, file, str(center[iter]))
                print "Using offsets for Frame 0 (%s) instead." % str(center[0])
                subreg = f[0].data[xstartval:xendval, ystartval:yendval]
            except:
                print "Some error occurred. Using offsets for Frame 0 (%s) instead." % \
                    str(center[0])
                subreg = f[0].data[xstartval:xendval, ystartval:yendval]                

        stack[iter] = subreg
        iter += 1
            
    return stack
    



def aperphot(fn, timekey=None, pos=[0,0], dap=[2,4,6], mask=None, verbose=False, nanval=999, resamp=None):
    """Do aperture photometry on a specified file.

    pos : 2-sequence
      center of apertures (as if indexing: fn[pos[0], pos[1]])

    dap : 3-sequence
      Photometry aperture DIAMETERS:
         -- target aperture (within which to sum flux)
         -- inner sky aperture
         -- outer sky aperture

    resamp : int
      Factor by which to interpolate frame before measuring photometry
      (in essence, does partial-pixel aperture photometry)

    Aperture masking:
       If no mask is passed in, use the star's input position and
       aperture diameters to create binary pixel masks for stellar and
       background photometry.  If a mask is passed in, stellar and
       background aperture masks are generated from where all input
       mask elements equal 1 and 2, respectively.

    :REQUIREMENTS:
      scipy.interpolate, pyfits, numpy...
           """
    # 2009-09-14 10:49 IJC: Created
    # 2010-01-15 14:20 IJC: Added numpy "_string" check
    # 2011-12-29 12:01 IJMC: Added peak pixel values to photometry report.
    # 2012-01-25 11:26 IJMC: Adding "resamp" option -- thanks to
    #                        K. Stevenson and J. Harrington of UCF for
    #                        the suggestion.
    # 2012-02-26 11:53 IJMC: Now return 'ntarg' and 'nsky' -- number of pixels used.

    from numpy import meshgrid, median,isfinite,sort,ndarray,string_
    import numpy as np
    import pyfits
    from analysis import fixval
    from os import path
    from scipy import interpolate

    thisobs = phot()
    x0, y0 = pos
    dap_targ, dap_skyinner, dap_skyouter = dap
    if resamp is None or resamp<1:
        resamp = 1
    else:
        resamp = float(resamp)
        
    # Determine size:
    if isinstance(fn,str):
        nx = pyfits.getval(fn, 'NAXIS1')
        ny = pyfits.getval(fn, 'NAXIS2')
    elif isinstance(fn,ndarray):
        nx,ny = fn.shape

    nx0, ny0 = nx, ny
    nx = ((nx - 1)*resamp + 1.)  # Avoid resampling at pixel locations
    ny = ((ny - 1)*resamp + 1.)  #   outside the original boundaries.

    # Generate or load masks:
    if mask==None:
        xx,yy = meshgrid(np.arange(nx)/resamp, np.arange(ny)/resamp)
        mask_targ = makemask(xx, yy, (x0, y0, dap_targ))
        mask_s1 = makemask(xx, yy, (x0,y0, dap_skyinner))
        mask_s2 = makemask(xx, yy, (x0,y0, dap_skyouter))
        mask_sky = mask_s2 - mask_s1
    else:
        mask_targ = mask==1
        mask_sky = mask==2
        if resamp>1:
            print "In aperphot, resamp>1 and user-specified mask passed in... beware!"

    # Load data frame:
    thisobs = phot()
    if isinstance(fn,ndarray):
        frame = fn
    elif isinstance(fn, str) or isinstance(fn,string_):
        if not path.isfile(fn):
            print "file %s not found! exiting..." % fn
            return thisobs
        frame = pyfits.getdata(fn)
        fixval(frame, nanval)

    # Resample data frame
    if resamp>1:
        frame0 = frame.copy()
        xx0 = range(nx0)
        yy0 = range(ny0)
        x1,y1 = np.arange(nx)/resamp, np.arange(ny)/resamp
        rectspline = interpolate.fitpack2.RectBivariateSpline(xx0, yy0, frame0, kx=1, ky=1, s=0)
        frame = rectspline(x1, y1)/resamp/resamp

    #from pylab import *
    #pdb.set_trace()
    # Measure background and aperture photometry
    thisbg, thisebg = estbg(frame, mask=mask_sky, plotalot=verbose, rout=[3,99])
    thisphot = (mask_targ*(frame - thisbg)).sum()
    peak = frame.max()
    peak_targ = (mask_targ * frame).max()
    peak_annulus = (mask_sky * frame).max()

    thisobs.bg=thisbg
    thisobs.ebg=thisebg
    thisobs.bgstr='phot.estbg: SDOM on bg histogram mean & dispersion after outlier rejection'
    thisobs.phot=thisphot
    thisobs.photstr='by-hand background-subtracted aperture photometry'
    thisobs.ntarg = mask_targ.sum()/resamp/resamp
    thisobs.nsky = mask_sky.sum()/resamp/resamp

    thisobs.peak = peak
    thisobs.peak_targ = peak_targ
    thisobs.peak_annulus = peak_annulus
    thisobs.peakstr = 'peak pixel value in frame'
    thisobs.peak_targstr = 'peak pixel value in target aperture'
    thisobs.peak_annulusstr = 'peak pixel value in sky annulus'
    thisobs.position = pos
    thisobs.positionstr = 'user-specified, zero-indexed pixel coordinates.'
    if isinstance(fn, str):
        if not timekey==None:
            thisobs.time=pyfits.getval(fn, timekey)
        thisobs.timestr='heliocentric modified julian date'
        thisobs.object=pyfits.getval(fn,'object')
        thisobs.exptime=pyfits.getval(fn,'exptime')
        thisobs.aper = dap
    thisobs.aperstr = 'target, inner, outer aperture diameters, in pixels.'
    thisobs.filename=fn
    if verbose:
        from pylab import figure, colorbar
        from nsdata import imshow
        figure(); imshow(frame*mask_targ); colorbar()

    return thisobs

def psffiterr(xyoffset, psf, frame, w=None, scale=100, dframe=9, verbose=False):
    """ Usage:
          optimize.fmin_powell(psffiterr, [0,0], args=(psf, stack[6],goodind[6],100,13),xtol=0.5,ftol=1)
          """
    # 2009-11-13 11:27 IJC: Created
    from numpy import round
    import pdb
    out = psffit(psf, frame, loc=None, w=w, scale=scale, dframe=dframe, \
                     xoffs=[round(xyoffset[0])], yoffs=[round(xyoffset[1])], verbose=verbose)
    sumsqres = out[2][0,0]
    #pdb.set_trace()
    return sumsqres

def psffit(psf, frame, loc=None, w=None, scale=100, dframe=9, xoffs=range(0,100,10), yoffs=range(0,100,10), verbose=False):
    """
    INPUT:
       psf -- model PSF (supersampled by 'scale')
       frame -- science frame to which psf will be fit. best if it's odd-sized
       
    OPTIONAL INPUTS:
       loc -- (x,y) integers, star location in data frame (e.g., data[x,y])
       w -- [array] weights of pixels in data frame, (typ. 1/sigma^2)
       scale -- supersampling level of PSF
       dframe -- [odd int] diameter of square box around target location
       xoffs -- [int array] subpixel x offsets to test.  
       yoffs -- [int array] subpixel y offsets to test.  
    """
    # 2009-10-07 14:18 IJC: Created.
    # 2011-05-11 22:16 IJC: Added a try/except/pdb debugging step
    from numpy import arange, prod, zeros,vstack,dot, nonzero, diag,floor,abs,max,hstack
    from numpy.linalg import pinv
    from analysis import binarray
    from time import time
    import pdb

    tic = time()

    
    if w==None:
        w = zeros(frame.shape)+1
    if loc==None:
        loc = ((frame.shape[1]-1)/2,(frame.shape[0]-1)/2)
    if xoffs==None:
        xoffsets = arange(scale) #arange(scale)
    if yoffs==None:
        yoffs = arange(scale) #arange(scale)

    ycen = int(loc[0])
    xcen = int(loc[1])
    pycen, pxcen = nonzero(psf==psf.max())
    # Pick a thumbnail frame size to use, and cut it out of the larger frame
    if verbose:
        print "frame.shape>>", frame.shape
        print "(xcen, ycen, dframe)>>", xcen,ycen,dframe
        print "limits>>" , (ycen-(dframe-1)/2.), (ycen+(dframe+1)/2.), (xcen-(dframe-1)/2.), (xcen+(dframe+1)/2.)
        print "xoffs, yoffs>>", xoffs, yoffs
    data = frame[(ycen-(dframe-1)/2.):(ycen+(dframe+1)/2.), 
                 (xcen-(dframe-1)/2.):(xcen+(dframe+1)/2.)]
    weights = w[(ycen-(dframe-1)/2.):(ycen+(dframe+1)/2.), 
                 (xcen-(dframe-1)/2.):(xcen+(dframe+1)/2.)]
    wmat = diag(weights.ravel())

    extrasize = 2*max(abs(floor(1.0*hstack((xoffs,yoffs))/scale)))
    exs = extrasize*scale/2
    if verbose: print "extrasize>> ",extrasize

    # Determine how much of the PSF to cut out, and cut it out
    dpsf0 = (dframe+1)*scale-1
    dpsf1 = dframe*scale-1
    pxmin = int(pxcen-(dpsf0-1)/2-exs)
    pxmax = int(pxcen+(dpsf0+1)/2+exs)
    pymin = int(pycen-(dpsf0-1)/2-exs)
    pymax = int(pycen+(dpsf0+1)/2+exs)
    if verbose: print psf.shape, pymin, pymax, pxmin, pxmax
    smpsf = psf[pymin:pymax, pxmin:pxmax]

    if verbose: print "data.shape>>" , data.shape
    ndata = prod(data.shape)
    if verbose: print "ndata>> %i" % ndata
    const = zeros(ndata,float)+1
    background = zeros((len(xoffs),len(yoffs)), float)
    fluxscale =  zeros((len(xoffs),len(yoffs)), float)
    chisq = zeros((len(xoffs),len(yoffs)), float)
    if verbose: 
        print "wmat.shape>>", wmat.shape
        print "data.ravel().shape>>", data.ravel().shape

    wpmat = dot(wmat, data.ravel()) # outside of loop for efficiency
    dfs = dframe * scale            # outside of loop for efficiency
    initoffset_min = scale - 1 + exs    # outside of loop for efficiency
    initoffset_max = scale - 1 + exs + dfs    # outside of loop for efficiency
    nx = len(xoffs)
    ny = len(yoffs)
    rangeny = range(ny)
    for ii in range(nx):
        xoffset = xoffs[ii]
        xmin, xmax = int(initoffset_min-xoffset), int(initoffset_max-xoffset)
        for jj in rangeny:
            #   Cut out the correct portion of the model PSF
            yoffset = yoffs[jj]
            ymin, ymax = int(initoffset_min-yoffset), int(initoffset_max-yoffset)
            #   Bin down the PSF by the correct factor.  Sizes should now match!
            binpsf = binarray(smpsf[ymin:ymax, xmin:xmax],scale)
            #   Compute the best-fit background & PSF scaling factor
            if verbose:
                print "xmat shapes>> ",const.shape, binpsf.ravel().shape
                print "ymin,ymax,xmin,xmax>> ",ymin,ymax, xmin,xmax
                print "binpsf.shape>> ",binpsf.shape
            try:
                xmat = vstack((const,binpsf.ravel())).transpose()
            except:
                print "vstacking bombed..."
                pdb.set_trace()
                
            wxmat = dot(wmat,xmat)
            fitvec = dot(pinv(wxmat), wpmat) 
            background[ii,jj], fluxscale[ii,jj] = fitvec
            chisq[ii,jj] = ((dot(wxmat, fitvec) - wpmat)**2).sum()

    ii,jj = nonzero(chisq==chisq.min())
    if len(ii)>1: # more than one chisquared minimum found!
        ii = ii[0]
        jj = jj[0]
    xoffset = xoffs[ii]
    xmin, xmax = int(scale-xoffset-1+exs), int(scale-xoffset-1 + dfs+exs)
    yoffset = yoffs[jj]
    ymin, ymax = int(scale-yoffset-1+exs), int(scale-yoffset-1 + dfs+exs)
    #   Bin down the PSF by the correct factor.  Sizes should now match!
    binpsf = binarray(smpsf[ymin:ymax, xmin:xmax],scale)
    modelpsf = background[ii,jj]+fluxscale[ii,jj]*binpsf
    #print "%f seconds for completion!" % (time()-tic)
    return modelpsf, data, chisq, background, fluxscale, xoffset, yoffset, chisq[ii,jj],background[ii,jj], fluxscale[ii,jj]




def prffit(prf, frame, loc=None, w=None, scale=100, dframe=9, xoffs=range(0,100,10), yoffs=range(0,100,10), verbose=False):
    """
    INPUT:
       prf -- model PRF (supersampled by 'scale')
       frame -- science frame to which psf will be fit. best if it's odd-sized
       
    OPTIONAL INPUTS:
       loc -- (x,y) integers, star location in data frame (e.g., data[x,y])
       w -- [array] weights of pixels in data frame, (typ. 1/sigma^2)
       scale -- supersampling level of PSF
       dframe -- [odd int] diameter of square box around target location
       xoffs -- [int array] subpixel x offsets to test.  
       yoffs -- [int array] subpixel y offsets to test.  
    """
    # 2009-10-07 14:18 IJC: Created.
    from numpy import arange, prod, zeros,vstack,dot, nonzero, diag,floor,abs,max,hstack
    from numpy.linalg import pinv
    from analysis import binarray
    from time import time
    tic = time()

    
    if w==None:
        w = zeros(frame.shape)+1
    if loc==None:
        loc = ((frame.shape[1]-1)/2,(frame.shape[0]-1)/2)
    if xoffs==None:
        xoffsets = arange(scale) #arange(scale)
    if yoffs==None:
        yoffs = arange(scale) #arange(scale)

    ycen = int(loc[0])
    xcen = int(loc[1])
    pycen, pxcen = nonzero(prf==prf.max())
    # Pick a thumbnail frame size to use, and cut it out of the larger frame
    if verbose:
        print "frame.shape>>", frame.shape
        print "(xcen, ycen, dframe)>>", xcen,ycen,dframe
        print "limits>>" , (ycen-(dframe-1)/2.), (ycen+(dframe+1)/2.), (xcen-(dframe-1)/2.), (xcen+(dframe+1)/2.)
    data = frame[(ycen-(dframe-1)/2.):(ycen+(dframe+1)/2.), 
                 (xcen-(dframe-1)/2.):(xcen+(dframe+1)/2.)]
    weights = w[(ycen-(dframe-1)/2.):(ycen+(dframe+1)/2.), 
                 (xcen-(dframe-1)/2.):(xcen+(dframe+1)/2.)]
    wmat = diag(weights.ravel())

    extrasize = 2*max(abs(floor(1.0*hstack((xoffs,yoffs))/scale)))
    exs = extrasize*scale/2
    if verbose: print "extrasize>> ",extrasize


    if verbose: print "data.shape>>" , data.shape
    ndata = prod(data.shape)
    if verbose: print "ndata>> %i" % ndata
    const = zeros(ndata,float)+1
    background = zeros((len(xoffs),len(yoffs)), float)
    fluxscale =  zeros((len(xoffs),len(yoffs)), float)
    chisq = zeros((len(xoffs),len(yoffs)), float)
    if verbose: 
        print "wmat.shape>>", wmat.shape
        print "data.ravel().shape>>", data.ravel().shape

    wpmat = dot(wmat, data.ravel()) # outside of loop for efficiency
    dfs = dframe * scale            # outside of loop for efficiency
    initoffset_min = scale - 1 + exs    # outside of loop for efficiency
    initoffset_max = scale - 1 + exs + dfs    # outside of loop for efficiency
    nx = len(xoffs)
    ny = len(yoffs)
    rangeny = range(ny)
    for ii in range(nx):
        xoffset = xoffs[ii]
        xmin, xmax = int(initoffset_min-xoffset), int(initoffset_max-xoffset)
        for jj in rangeny:
            #   Extract the correct indices of the model PRF. Sizes should now match!
            yoffset = yoffs[jj]
            ymin, ymax = int(initoffset_min-yoffset), int(initoffset_max-yoffset)
            binpsf = prf[indy, indx]
            #   Compute the best-fit background & PSF scaling factor
            if verbose and (ii==0) and (jj==0): 
                print "xmat shapes>> ",const.shape, binpsf.ravel().shape
                print "ymin,ymax,xmin,xmax>> ",ymin,ymax, xmin,xmax
                print "binpsf.shape>> ",binpsf.shape
            xmat = vstack((const,binpsf.ravel())).transpose()
            wxmat = dot(wmat,xmat)
            fitvec = dot(pinv(wxmat), wpmat) 
            background[ii,jj], fluxscale[ii,jj] = fitvec
            chisq[ii,jj] = ((dot(wxmat, fitvec) - wpmat)**2).sum()

    ii,jj = nonzero(chisq==chisq.min())
    if len(ii)>1: # more than one chisquared minimum found!
        ii = ii[0]
        jj = jj[0]
    xoffset = xoffs[ii]
    xmin, xmax = int(scale-xoffset-1+exs), int(scale-xoffset-1 + dfs+exs)
    yoffset = yoffs[jj]
    ymin, ymax = int(scale-yoffset-1+exs), int(scale-yoffset-1 + dfs+exs)
    #   Bin down the PSF by the correct factor.  Sizes should now match!
    binpsf = binarray(smpsf[ymin:ymax, xmin:xmax],scale)
    modelpsf = background[ii,jj]+fluxscale[ii,jj]*binpsf
    print "%f seconds for completion!" % (time()-tic)
    return modelpsf, data, chisq, background, fluxscale, xoffset, yoffset, chisq[ii,jj],background[ii,jj], fluxscale[ii,jj]


def gauss2d(param, x, y):
    """ Compute a gaussian distribution at the points x, y.

        p is a three- or four-component array, list, or tuple:

        z =  [p4 +] p0/(p1*sqrt(4pi)) * exp(-r**2 / (2*p1**2))

             where r = sqrt((x-p2)**2 + (y-p3)**2)

        p[0] -- Volume under the gaussian
        p[1] -- one-sigma dispersion
        p[2] -- X- offset of center 
        p[3] -- Y- offset of center
        p[4] -- optional constant, vertical offset

        x & y must be equal-sized arrays from numpy.meshgrid

        NOTE: FWHM = 2*sqrt(2*ln(2)) * p1  ~ 2.3548*p1

        SEE ALSO: egauss2d, numpy.meshgrid"""
    #2010-01-11 22:46 IJC: Created
    from numpy import array, abs, concatenate, exp
    x = array(x, dtype=float).copy()
    y = array(y, dtype=float).copy()
    p = array(param).copy()

    r = abs((x-p[2]) + 1j*(y-p[3]))

    if len(p)==4:
        p = concatenate((p, [0]))

    z = p[4] + p[0]/(p[1]*4*pi) * exp(-r**2 / (2*p[1]**2))
    
    return z
    
def egauss2d(param, x, y, z, ez=None):
    """ Return the chi-squared error on a 2D gaussian fit.  See gauss2d
        for more details on the input parameters.

        z is the array of data to be fit
        ez is an optional array of one-sigma errors to the data in z.
        """
    # 2010-01-11 22:59 IJC: Created
    from numpy import array, float, ones
    x = array(x, dtype=float).copy()
    y = array(y, dtype=float).copy()
    z = array(z, dtype=float).copy()
    p = array(param).copy()
    if ez==None:
        ez = ones(z.shape,float)
    else:
        ez = array(ez, dtype=float).copy()        

    model = gauss2d(param,x,y)
    chisq = (((z-model)/ez)**2).sum()
    return chisq




