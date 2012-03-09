"""
Time-series photometric reduction for ground-based IR photometry.

Assumes a single filter and a regular dither pattern (or 'staring'),
and sky (& dark) removal via subtraction of temporally adjacent frames.
"""
# 2009-12-16 08:36 IJC: Begun.

#from pyraf import iraf
import astrolib
import phot, pyfits, os
import analysis as an
import numpy as np

#_iraf = "/Users/ianc/iraf/"
#iraf.load('fitsutil')
#iraf.load('astutil')
#iraf.unlearn('ccdproc')
#iraf.unlearn('imcombine')


class baseObject:
    """Empty object container.
    """
    def __init__(self):
        return


def initobs(obsname='20091203k', verbose=False, cal='dome'):
    """Initialize variables for ground-based IR data analysis.

    OPTIONAL INPUT:
       obsname: (str) -- e.g., 20091203k, 20091203h

       cal : str
         Valid options are 'dome' and 'sky', in case two sets of
         calibration files are available.  If not, this option will
         have no effect.

    OUTPUT:
       a DICT containing the following values, in order:
         """
    # 2009-12-14 22:48 IJC: Created.
    # 2010-01-20 16:16 IJC: Added IRTF/SpeX Dec2010 observations
    # 2011-11-09 08:50 IJMC: Added some of the 2011 IRTF/SpeX Gj1214b observations
    # 2011-12-16 15:14 IJMC: Added cal keyword, and Subaru/MOIRCS photometry
    from numpy import tile, sort
    from nsdata import gd2jd


    if obsname=='20091203k':
        planet = 'WASP-12 b'
        _raw = '/Users/ianc/proj/transit/data/raw/20091203g/'
        _proc = '/Users/ianc/proj/transit/data/proc/20091203g/'
        _prefix = '03DEB'
        _suffix = '.FTS'
        flats = range(2,43) + range(45,52)
        darks = range(52,92)
        sci = range(156,620)
        system = 'lick/geminib'
        timezone = 'PST'
        ndither = 8
        procsuffix = '_df'
        numdigits = 3
        
    elif obsname=='20091203h':
        planet = 'WASP-12 b'
        _raw = '/Users/ianc/proj/transit/data/raw/20091203g/'
        _proc = '/Users/ianc/proj/transit/data/proc/20091203g/'
        _prefix = '03DEA'
        _suffix = '.FTS'
        flats = range(2,43)
        darks = range(53,113)
        sci = range(177,641)
        system = 'lick/geminia'
        timezone = 'PST'
        ndither = 8
        procsuffix = '_df'
        numdigits = 3

    elif obsname=='20091228a':
        planet = 'WASP-12 b'
        _raw = '/Users/ianc/proj/pcsa/data/raw/2009dec28/'
        _proc = '/Users/ianc/proj/pcsa/data/proc/wasp12_2009dec28/'
        _prefix = 'spectra'
        _suffix = '.fits'
        _rawprefix = 'spc'
        _rawsuffix = '.a.fits'
        flats = [-1]
        darks = [-1]
        #sci = range(176,678)
        sci = list(np.sort(np.hstack(([176],np.hstack((178+4*np.arange(125),181+4*np.arange(125)))))))
        system = 'irtf/spexprism'
        timezone = 'HST'
        ndither = 2
        procsuffix = ''
        numdigits = 4
    
    elif obsname=='20091228b':
        planet = 'WASP-12 b'
        _raw = '/Users/ianc/proj/pcsa/data/raw/2009dec28/'
        _proc = '/Users/ianc/proj/pcsa/data/proc/wasp12_2009dec28/'
        _prefix = 'spectra'
        _suffix = '.fits'
        _rawprefix = 'spc'
        _rawsuffix = '.b.fits'
        flats = [-1]
        darks = [-1]
        #sci = range(176,678)
        sci = list(np.sort(np.hstack(([177],np.hstack((179+4*np.arange(125),180+4*np.arange(125)))))))
        system = 'irtf/spexprism'
        timezone = 'HST'
        ndither = 2
        procsuffix = ''
        numdigits = 4

    elif obsname=='20091230a':
        planet = 'WASP-12 b'
        _raw = '/Users/ianc/proj/pcsa/data/raw/2009dec30/'
        _proc = '/Users/ianc/proj/pcsa/data/proc/wasp12_2009dec30/'
        _prefix = 'spectra'
        _suffix = '.fits'
        _rawprefix = 'spc'
        _rawsuffix = '.a.fits'
        flats = [-1]
        darks = [-1]
        #sci = range(1,389)
        sci = list(np.sort(np.hstack((3+np.arange(30)*4, 6+np.arange(30)*4, 123, 125+np.arange(66)*4,128+np.arange(66)*4))))
        for badframe in [158,159,160,161,162,163,164]:
            try:
                sci.remove(badframe)
            except:
                trash = False
        system = 'irtf/spexprism'
        timezone = 'HST'
        ndither = 2
        procsuffix = ''
        numdigits = 4
    
    elif obsname=='20091230b':
        planet = 'WASP-12 b'
        _raw = '/Users/ianc/proj/pcsa/data/raw/2009dec30/'
        _proc = '/Users/ianc/proj/pcsa/data/proc/wasp12_2009dec30/'
        _prefix = 'spectra'
        _suffix = '.fits'
        _rawprefix = 'spc'
        _rawsuffix = '.b.fits'
        flats = [-1]
        darks = [-1]
        #sci = range(1,389)
        sci = list(np.sort(np.hstack((4+np.arange(30)*4, 5+np.arange(30)*4, 124, 126+np.arange(66)*4,127+np.arange(66)*4))))
        for badframe in [158,159,160,161,162,163,164]:
            try:
                sci.remove(badframe)
            except:
                trash = False
        system = 'irtf/spexprism'
        timezone = 'HST'
        ndither = 2
        procsuffix = ''
        numdigits = 4

    elif obsname=='20091230guide':
        planet = 'WASP-12 b'
        _raw = '/Users/ianc/proj/transit/data/raw/20091230/guidedog/'
        _proc = '<noprocdata>'
        _prefix = 'spc'
        _suffix = '.a.fits'
        flats = [-1]
        darks = [-1]
        sci = range(1,1257)
        for badframe in []:
            try:
                sci.remove(badframe)
            except:
                trash = False
        system = 'irtf/guide'
        timezone = 'HST'
        ndither = 1
        procsuffix = ''
        numdigits = 4

    elif obsname=='2011aug05':
        planet = 'GJ 1214 b'
        _raw = '/Users/ianc/proj/pcsa/data/raw/2011aug05/spec/'
        _proc = '/Users/ianc/proj/pcsa/data/proc/gj1214_2011aug05/'
        _prefix = 'spectra'
        _suffix = '.fits'
        _rawprefix = 'filename'
        _rawsuffix = '.a.fits'
        flats = [-1]
        darks = [-1]
        sci = range(251,462)
        system = 'irtf/spexsxd'
        timezone = 'HST'
        ndither = 1
        procsuffix = ''
        numdigits = 4
    
    elif obsname=='2011aug05cal':
        planet = 'HD 161289'
        _raw = '/Users/ianc/proj/pcsa/data/raw/2011aug05/spec/'
        _proc = '/Users/ianc/proj/pcsa/data/proc/gj1214_2011aug05/'
        _prefix = 'spectra'
        _suffix = '.fits'
        _rawprefix = 'filename'
        _rawsuffix = '.a.fits'
        flats = [-1]
        darks = [-1]
        sci = range(462,502)
        system = 'irtf/spexsxd'
        timezone = 'HST'
        ndither = 1
        procsuffix = ''
        numdigits = 4

    elif obsname=='2011aug16':
        planet = 'GJ 1214 b'
        _raw = '/Users/ianc/proj/pcsa/data/raw/2011aug16/spec/'
        _proc = '/Users/ianc/proj/pcsa/data/proc/gj1214_2011aug16/'
        _prefix = 'spectra01_'
        _suffix = '.fits'
        _rawprefix = 'filename'
        _rawsuffix = '.a.fits'
        flats = [-1]
        darks = [-1]
        sci = range(1, 199)
        system = 'irtf/spexsxd'
        timezone = 'HST'
        ndither = 1
        procsuffix = ''
        numdigits = 4

    elif obsname=='2011sep12':
        planet = 'GJ 1214 b'
        _raw = '/Users/ianc/proj/pcsa/data/raw/2011sep12/spec/'
        _proc = '/Users/ianc/proj/pcsa/data/proc/gj1214_2011sep12/'
        _prefix = 'spectra01_'
        _suffix = '.fits'
        _rawprefix = 'filename'
        _rawsuffix = '.a.fits'
        flats = [-1]
        darks = [-1]
        sci = range(161, 185) + range(190, 353)
        system = 'irtf/spexsxd'
        timezone = 'HST'
        ndither = 1
        procsuffix = ''
        numdigits = 4

    elif obsname=='2011oct01':
        planet = 'GJ 1214 b'
        _raw = '/Users/ianc/proj/pcsa/data/raw/2011oct01/spec/'
        _proc = '/Users/ianc/proj/pcsa/data/proc/gj1214_2011oct01/'
        _prefix = 'spectra01_'
        _suffix = '.fits'
        _rawprefix = 'spc'
        _rawsuffix = '.a.fits'
        flats = [-1]
        darks = [-1]
        sci = range(151, 283)
        system = 'irtf/spexsxd'
        timezone = 'HST'
        ndither = 1
        procsuffix = ''
        numdigits = 4


    elif obsname=='20111214_chip1':
        planet = 'WASP-12 b'
        _raw = '/Users/ianc/proj/transit/data/20111214/'
        _proc = '/Users/ianc/proj/transit/data/20111214_proc/'
        _prefix = 'proc_MCSA0018'
        _suffix = '.fits'
        _rawprefix = 'MCSA0018'
        _rawsuffix = '.fits'
        skyflats = range(7737, 7770, 2)
        skydarks = range(7833, 7853, 2) # Corresponds to flat-field frames
        domeflats = range(7789, 7803, 2)
        domedarks = range(7773, 7787, 2) 
        sci = range(6943, 7728, 2)
        system = 'subaru/moircs'
        timezone = 'UTC'
        ndither = 1
        procsuffix = ''
        numdigits = 4
        badpixelmask = _raw + 'badpixmap_20111214_IC_ch1.fits'
        for badframe in [7271, 7361, 7363]:
            try:
                sci.remove(badframe)
            except:
                pass

    elif obsname=='20111214_chip2':
        planet = 'WASP-12 b'
        _raw = '/Users/ianc/proj/transit/data/20111214/'
        _proc = '/Users/ianc/proj/transit/data/20111214_proc/'
        _prefix = 'proc_MCSA0018'
        _suffix = '.fits'
        _rawprefix = 'MCSA0018'
        _rawsuffix = '.fits'
        skyflats = range(7738, 7770, 2)
        skydarks = range(7834, 7853, 2) # Corresponds to flat-field frames
        domeflats = range(7790, 7803, 2)
        domedarks = range(7774, 7787, 2) 
        sci = range(6944, 7729, 2)
        system = 'subaru/moircs'
        timezone = 'UTC'
        ndither = 1
        procsuffix = ''
        numdigits = 4
        badpixelmask = _raw + 'badpixmap_20111214_IC_ch2.fits'
        for badframe in [7272, 7362, 7364]:
            try:
                sci.remove(badframe)
            except:
                pass

    masterflat = ''
    if cal=='dome': 
        if 'domeflats' in locals():
            flats = domeflats
            masterflat = _proc + 'masterdomeflat' + _suffix
        if 'domedarks' in locals():
            darks = domedarks
    elif cal=='sky': 
        if 'skyflats' in locals():
            flats = skyflats
            masterflat = _proc + 'masterskyflat' + _suffix
        if 'skydarks' in locals():
            darks = skydarks
    else:
        print "'cal' keyword must be either 'dome' or 'sky'!"

    if masterflat=='':
        masterflat = _proc + 'masterflat' + _suffix
    
    if _rawprefix and _rawsuffix:
        pass
    else:
        _rawprefix = _prefix
        _rawsuffix = _suffix

    if timezone == 'PST':
        toffset = 8
    elif timezone=='PDT':
        toffset = 7
    elif timezone=='HST':
        toffset = 10
    elif timezone=='UTC':
        toffset = 0

    fmtstr = '%0'+str(numdigits)+'g'

    allindices = darks+flats+sci
    allindex = np.argsort(allindices)
    rawall = [(_raw+_rawprefix+ (fmtstr % num) +_rawsuffix) for num in np.sort(allindices)]
    types = np.array(['dark']*len(darks) + ['flat']*len(flats) + ['sci']*len(sci))[allindex]

    rawdark = [(_raw + _rawprefix + (fmtstr % num) + _rawsuffix) for num in darks]
    rawflat = [(_raw + _rawprefix + (fmtstr % num) + _rawsuffix) for num in flats]
    rawsci = [(_raw + _rawprefix + (fmtstr % num) + _rawsuffix) for num in sci]

    dsubflat = [(_proc + _prefix +(fmtstr % num) +'_d'+ _suffix) for num in flats]
    dsubsci = [(_proc + _prefix + (fmtstr % num) +'_d'+ _suffix) for num in sci]
    dsuball = [(_proc + _prefix + (fmtstr % num) +'_d'+ _suffix) for num in np.sort(allindices)]
    procsci = [(_proc + _prefix + (fmtstr % num) +procsuffix + _suffix) for num in sci]

    if system=='lick/geminia':
        gain = 0      # photons (i.e., electrons) per data unit  
        readnoise = 0 # photons (i.e., electrons) per individual read               
        naxis1, naxis2 = 256, 256
    elif system=='lick/geminib':
        gain = 0      # photons (i.e., electrons) per data unit  
        readnoise = 0 # photons (i.e., electrons) per individual read               
        naxis1, naxis2 = 256, 256
    elif system.find('irtf/spex')>=0:
        gain = 13      # photons (i.e., electrons) per data unit  
        readnoise = 50 # photons (i.e., electrons) per individual read
        naxis1,naxis2 = 1024,1024
    elif system.find('irtf/guide')>=0:
        gain = 14.7      # photons (i.e., electrons) per data unit  
        readnoise = 65 # photons (i.e., electrons) per individual read
        naxis1,naxis2 = 512,512

    if system.find('lick/gemini')>=0:
        headers = [pyfits.getheader(file) for file in rawall]
        timestrs = ['%s %s' % (h['date_obs'], h['time_obs']) \
                        for h in headers]
        jd = np.array([gd2jd(t) for t in timestrs])+toffset/24.
        coadds = np.array([h['coadds'] for h in headers])

    elif system.find('irtf/spex')>=0:
        headers = [pyfits.getheader(file) for file in procsci]
        timestrs = ['%s %s' % (h['date_obs'], h['time_obs']) \
                        for h in headers]
        jd = np.array([gd2jd(t) for t in timestrs])
        ra = map(phot.hms, [h['RA'] for h in headers])
        dec = map(phot.dms, [h['DEC'] for h in headers])
        hjd = np.array([astrolib.helio_jd(jdi - 2400000., rai, deci) + 2400000. for \
                            (jdi, rai, deci) in zip(jd, ra, dec)])
        coadds = np.array([h['co_adds'] for h in headers])

    elif system.find('irtf/guide')>=0:
        headers = [pyfits.getheader(file) for file in rawsci]
        timestrs = ['%s %s' % (h['date_obs'], h['time_obs']) \
                        for h in headers]
        jd = np.array([gd2jd(t) for t in timestrs])
        coadds = np.array([h['co_adds'] for h in headers])

    elif system.find('subaru/moircs')>=0:
        try:
            headers = [pyfits.getheader(file) for file in procsci]
        except:
            headers = [pyfits.getheader(file) for file in rawsci]
        mjd = [0.5 * (h['mjd-str'] + h['mjd-end']) for h in headers]
        jd = np.array(mjd) + 2400000.5
        planet_obj = an.getobj(planet)
        ra  = planet_obj.ra * 15
        dec = planet_obj.dec 
        hjd = np.array([astrolib.helio_jd(mjdi + 0.5, ra, dec) + \
                            2400000. for mjdi in mjd])
        coadds = np.array([h['coadd'] for h in headers])

        # Note that this gain is the average of the Chip 1 & 2 values:
        gain = 3.4      # photons (i.e., electrons) per data unit  
        readnoise = 30  # photons (i.e., electrons) per individual read
        naxis1,naxis2 = 2048, 2048


    if headers[0].has_key('itime'):
        itimes = np.array([h['itime'] for h in headers])
    elif headers[0].has_key('exp1time'):
        itimes = np.array([h['exp1time'] for h in headers])
    else:
        itimes = None

    obs = baseObject()
    obs._raw = _raw
    obs._proc = _proc
    obs.rawall  = np.array(rawall)
    obs.rawdark  = np.array(rawdark )
    obs.rawflat  = np.array(rawflat )
    obs.rawsci   = np.array(rawsci  )
    obs.dsubflat = np.array(dsubflat)
    obs.dsubsci  = np.array(dsubsci )
    obs.dsuball  = np.array(dsuball )
    obs.procsci  = np.array(procsci )
    obs.target   = planet 
    obs.gain = gain
    obs.readnoise = readnoise
    obs._suffix = _suffix
    obs._prefix = _prefix
    obs.types = types
    obs.jd= jd
    obs.itimes= itimes
    obs.coadds = coadds
    obs.naxis1 = headers[0]['naxis1']
    obs.naxis2 = headers[0]['naxis2']
    obs.ndither = ndither
    obs.system = system
    obs.headers = headers
    obs.masterflat = masterflat
    try:
        obs.hjd = hjd
    except:
        obs.hjd = None

    try:
        obs.badpixelmask = badpixelmask
    except:
        obs.badpixelmask = None

    return obs


def red_dark(obs, ind, num=np.nan, verbose=False, sigma=3, niter=1, clobber=False):
    """Combine dark frames into a single dark frame:"""

    _sdark = obs._proc + ('superdark%03g' % num) + obs._suffix
    darklist = obs.rawall[ind]
    framecenter = obs.naxis1/2, obs.naxis2/2
    framesize = obs.naxis1, obs.naxis2
    
    if verbose: print "dark file list>>" , darklist

    darkstack = phot.subreg2(darklist, framecenter, framesize, verbose=verbose)
    if verbose: print "loaded dark stack"
    superdark = np.zeros(framesize,float)
    for ii in range(obs.naxis1):
        for jj in range(obs.naxis2):
            superdark[ii,jj] = an.meanr(darkstack[:,ii,jj],nsigma=sigma, \
                                            niter=niter)
    hdr0 = pyfits.getheader(darklist[0])
    hdr0.update('sup-dark', 'super-dark created by IR.py')
    pyfits.writeto(_sdark, superdark, header=hdr0, clobber=clobber)

    if verbose: print "Done making superdark!"
    return

def red_sub(ifiles, ofiles, subfile, num=np.nan, clobber=False):
    """Subtract one file from a list of others, and save in a new location.

    First two inputs are lists of strings, 'subfile' is a string.
    """
    subtractor = pyfits.getdata(subfile)
    for infile, outfile in zip(ifiles,ofiles):
        thisin = pyfits.getdata(infile)
        thishdr = pyfits.getheader(infile)
        thishdr.update('subtract', '%s minus %s' % (os.path.split(infile)[-1], os.path.split(subfile)[-1]))
        pyfits.writeto(outfile, thisin-subtractor, header=thishdr, clobber=clobber)
    return
        

def red_flat(obs, ind, num=np.nan, verbose=False, sigma=3, niter=1, clobber=False):
    """Dark-correct flats, make normalized super-flat."""
    _sdark = obs._proc + ('superdark%03g' % num) + obs._suffix
    _sflat = obs._proc + ('superflat%03g' % num) + obs._suffix
    rawflat = obs.rawall[ind]
    dsubflat = obs.dsuball[ind]
    nflat = len(rawflat)
    framecenter = obs.naxis1/2, obs.naxis2/2
    framesize = obs.naxis1, obs.naxis2

    red_sub(rawflat, dsubflat, _sdark, clobber=True)
    if verbose: print "dark-subtracted flat frames"
    
    flatstack = phot.subreg(dsubflat, framecenter, framesize, verbose=verbose)
    if verbose: print "read in flat stack"

    flatmedians = an.amedian(flatstack.reshape((nflat,np.prod(framesize))),1).reshape((nflat,1,1))
    superflat = np.zeros(framesize,float)
    for ii in range(obs.naxis1):
        for jj in range(obs.naxis2):
            superflat[ii,jj] = an.meanr(flatstack[:,ii,jj]/flatmedians,nsigma=sigma, \
                                            niter=niter)
    if verbose: print "computed median stack"

    hdr0 = pyfits.getheader(dsubflat[0])
    hdr0.update('sup-flat', 'super-flat created by IR.py')
    pyfits.writeto(_sflat, superflat, header=hdr0, clobber=clobber)
    if verbose: print "wrote super-flat"
    return

def red_masterflat(obs, nums=[], verbose=False, sigma=3, niter=1, clobber=False):
    """Combine numbered super-flats into a single flat."""
    sflats = [obs._proc + ('superflat%03g' % num) + obs._suffix for num in nums]
    nflat = len(sflats)
    framecenter = obs.naxis1/2, obs.naxis2/2
    framesize = obs.naxis1, obs.naxis2
    
    sflatstack = phot.subreg(sflats, framecenter, framesize, verbose=verbose)
    if verbose: print "read in super-flat stack"

    sflatmedians = an.amedian(sflatstack.reshape((nflat,np.prod(framesize))),1).reshape((nflat,1,1))
    masterflat = np.zeros(framesize,float)
    for ii in range(obs.naxis1):
        for jj in range(obs.naxis2):
            masterflat[ii,jj] = an.meanr(sflatstack[:,ii,jj]/sflatmedians,nsigma=sigma, \
                                            niter=niter, verbose=verbose)
    if verbose: print "computed master-flat"

    hdr0 = pyfits.getheader(sflats[0])
    hdr0.update('mas-flat', 'master-flat created by IR.py')
    pyfits.writeto(obs.masterflat, masterflat, header=hdr0, clobber=clobber)
    if verbose: print "wrote master-flat"
    return


def red_sci(obs, ind, num=np.nan, verbose=False, sigma=3, niter=1):
    """Dark-subtract science frames.

    I don't currently use this (Jan 2010)"""
    red_sub(rawflat, dsubflat, _sdark, clobber=True)

    rawsci = obs.rawall[ind]
    dsubsci = obs.dsuball[ind]
    nsci = len(rawsci
)
    framecenter = obs.naxis1/2, obs.naxis2/2
    framesize = obs.naxis1, obs.naxis2

    red_sub(rawsci, dsubsci, _sdark, clobber=True)
    
    return

def ditherskysub(obs, sci_index, clobber=False, verbose=False, flatten=True, skyramp=False):
    """
    Subtract sky (and dark) levels from a set of dithered exposures
    with equal exposure settings.

    For now, compute a robust mean sky frame from each ndither set.
    Then subtract them, flat-field them, and output. 

    NOTE: unlike other reduction tasks, the input index here refers
    ONLY to science frames.
    """
    from pyraf import iraf

    rawsci = obs.rawsci[sci_index]
    procsci = obs.procsci[sci_index]
    nsci = len(rawsci)
    ndither = obs.ndither
    
    masterflat = pyfits.getdata(obs.masterflat)

    framecenter = obs.naxis1/2, obs.naxis2/2
    framesize = obs.naxis1, obs.naxis2
    rawstack = phot.subreg(rawsci, framecenter, framesize)
    procstack = np.zeros(rawstack.shape,float)
    if verbose: print "Read in raw science stack!"

    if (nsci % ndither)<>0:
        print "WARNING: Input file list length (%i) not evenly"+ \
            "divisible by size of dither set (%i)." \
            % (nsci, ndither)
    nsets = np.ceil(nsci/ndither)
    if skyramp:
        print "not tested yet -- and probably quite slow."
        j0 = int(obs.jd[0])
        for ii in range(nsets):
            i0 = ii*ndither
            i1 = min(ii*ndither+ndither, nsci)
            thesejd = obs.jd[i0:i1]-j0
            for kk in range(obs.naxis1):
                for ll in range(obs.naxis2):
                    thisfit = an.polyfitr(thesejd, rawstack[i0:i1,kk,ll], \
                                              1, sigma)
                    procstack[i0,i1,kk,ll] = rawstack[i0:i1,kk,ll] - \
                        polyval(thisfit, thesejd)
                    
            
            thissky = an.amedian(rawstack[i0:i1,:,:],0)
            procstack[i0:i1,:,:] = rawstack[i0:i1,:,:]-thissky
        if verbose: print "Sky-subtracted raw science data"
    else:
        subtractedSky = np.zeros(nsci,float)
        for ii in range(nsets):
            i0 = ii*ndither
            i1 = min(ii*ndither+ndither, nsci)
            thissky = an.amedian(rawstack[i0:i1,:,:],0)
            procstack[i0:i1,:,:] = rawstack[i0:i1,:,:]-thissky
            subtractedSky[i0:i1] = np.median(thissky.ravel())
        if verbose: print "Sky-subtracted raw science data"

    if flatten:
        procstack = procstack/masterflat
        if verbose: print "Flat-fielded science data"

    iraf.observatory('set','lick')
    for jj in range(nsci):
        hdr = pyfits.getheader(rawsci[jj])
        hdr.update('sky-sub', subtractedSky[jj])
        hdr.update('sky-note','Sky-subtracted using %i-frame dither set' \
                       % ndither)
        hdr.update('flatten','Flattened using %s' \
                       % os.path.split(obs.masterflat)[1] )
        hdr.update('exptime', hdr['coadds']*hdr['itime'])
        pyfits.writeto(procsci[jj], procstack[jj,:,:], hdr, clobber=clobber)
        setairmassG(proscsci[jj])

        


    if verbose: print "Wrote corrected science data to disk."
        
    return
    

def calibrate(obs, verbose=False, clobber=False):
    """
    Take raw NIR frames and output sky-subtracted, flat-fielded frames
    ready for photometry extraction.
    """
    # 2009-12-16 09:02 IJC: Begun.
    
    # initialize and create arrays of filenames, itimes, coadds, all
    # raw filenames (passed in via 'obs')

    # Create indices with which to extract various observation types:
    darkindex = obs.types=='dark'
    flatindex = obs.types=='flat'
    sciindex = obs.types=='sci'

    # From darks and flats, determine a list of overlapping ITIMEs
    darktimes = obs.itimes[darkindex]
    flattimes = obs.itimes[flatindex]

    caltimes = list(set(darktimes).intersection(flattimes))
    ncal = len(caltimes)

    for ii in range(ncal):
        thiscaltime = caltimes[ii]
        # combine darks into super-darks ('dark001', 'dark002', etc.) --
        #   for each dark-frame ITIME
        darkindex = (obs.types=='dark') * (obs.itimes==thiscaltime)
        red_dark(obs, darkindex, num=ii, verbose=verbose, clobber=clobber)

        # dark-subtract flats (_proc+procflat) and combine dark-subtracted
        #   flats into super-flats ('flat001', 'flat002', etc.) -- for
        #   each flat-frame ITIME
        flatindex = (obs.types=='flat') * (obs.itimes==thiscaltime)
        red_flat(obs, flatindex, num=ii, verbose=verbose, clobber=clobber)

    # Combine all super-flats into one master flat.
    red_masterflat(obs, range(ncal), verbose=verbose, clobber=clobber)

    # dark-subtract science frames (_proc+procsci) and flatten
    #   dark-subtracted science frames

    # red_sci()

    #iraf.setjd(output, date='date_obs', time='TIME_OBS', 
    #         epoch='equinox', jd='JD', hjd='hjd')


    return


def setairmassG(filename):
    """Use PyRaf tasks to set airmass for Lick/GEMINI camera observations.
    Assumes PST times."""
    from numpy import pi,sin,cos,arctan2,abs
    from pyraf import iraf
    import matplotlib.dates as dates

    lat =   (37. +20.6/60.)*pi/180.
    long = (121. +38.2/60.)*pi/180.
    epoch = 2000
    #alt = 1290 # meters


    
    hdr = pyfits.open(filename,mode='update')
    
    # Get data and time and put them in the correct format
    datestr = '%s %s PST' % (hdr[0].header['date_obs'],hdr[0].header['time_obs'])
    exptime = hdr[0].header['exptime']
    date = dates.num2date(dates.datestr2num(datestr))  # to UTC
    UTCdate = ('%i-%02g-%02gT%02g:%02g:%02g' % (date.year, date.month,date.day,date.hour,date.minute,date.second))
    
    # Get RA & dec in decimal degrees
    ra = (phot.hms(hdr[0].header['ra']))*pi/180.
    ha = (phot.hms(hdr[0].header['ha']))*pi/180. + 2*pi*exptime/86400.
    dec= (phot.dms(hdr[0].header['dec']))*pi/180.

    # Compute terms for coordinate conversion
    term1 = sin(lat)*sin(dec)+cos(lat)*cos(dec)*cos(ha)
    term2 = cos(lat)*sin(dec)-sin(lat)*cos(dec)*cos(ha)
    term3 = -cos(dec)*sin(ha)
    rad = abs(term3 +1j*term2)
    az = arctan2(term3,term2)
    alt = arctan2(term1, rad)

    # Compute airmass
    z = pi/2. - alt
    airmass = 1./(cos(z) + 0.50572*(96.07995-z*180./pi)**(-1.6364))


    hdr[0].header.update('airmass',airmass,comment='airmass at middle of observation')
    hdr[0].header.update('date-obs',UTCdate,comment='UTC date/time at start of observation')
    hdr[0].header.update('epoch',epoch)
    hdr.flush()
    hdr.close()

    return


def moircs_cal(skydark, skyflat, domedark, domeflat, objdark, rawframes, procframes, badpixelmask=None, clobber=False, retmaster=False, skysub=True, twoatatime=False, verbose=False):
    """Tanaka-san's prescription for MOIRCS flat-fielding.

    :INPUTS:
      skydark : list of str, or 2D Numpy array
        List of raw, long sky darks, or a 2D Numpy array Master Sky Dark

      skyflat : list of str, or 2D Numpy array
        List of raw, long sky flat frames, or a 2D Numpy array Master
        Sky Flat. Exposure time should be the same as for skydark.

      domedark : list of str, or 2D Numpy array
        List of raw, short dome darks, or a 2D Numpy array Master Dome Dark

      domeflat : list of str, or 2D Numpy array
        List of raw, short dome flats, or a 2D Numpy array Master Dome
        Flat. Exposure time per frame should be same as for skydark

      objdark : list of str, or 2D Numpy array
        List of raw, short darks, or a 2D Numpy array Master Dark.
        Per-frame exposure time should be the same as for `rawframes`

      rawframes : list of str
        List of raw, unprocessed science frames

      procframes : list of str
        List of processed science frames, to be written to disk

      clobber : bool
        Whether to overwrite existing FITS files

      retmaster : bool
        Whether to return master dark and (normalized) flat frames

      skysub : bool
        Whether to subtract the sky (not sure why you wouldn't want to...)

      twoatatime: bool
        Whether to correct both chips simultaneously, in 2-frame
        pairs.  If true, all the file inputs above must be len-2 lists
        of filename lists (or 2xN str arrays).

       
    :EXAMPLE:
      ::
          import ir
          obs2dome = ir.initobs('20111214_chip2', cal='dome')
          obs2sky  = ir.initobs('20111214_chip2', cal='sky')
          obj2dark = ['%sMCSA00187%i.fits' % (obs2sky._raw, ii) for ii in range(804,833,2)]
          ir.moircs_cal(obs2sky.rawdark, obs2sky.rawflat, \
                        obs2dome.rawdark, obs2dome.rawflat, obj2dark, \
                        obs2sky.rawsci, obs2sky.procsci, \
                        clobber=False, retmaster=False, verbose=True)

          obs1dome = ir.initobs('20111214_chip1', cal='dome')
          obs1sky  = ir.initobs('20111214_chip1', cal='sky')
          obj1dark = ['%sMCSA00187%i.fits' % (obs1sky._raw, ii) for ii in range(803,833,2)]
          ir.moircs_cal(obs1sky.rawdark, obs1sky.rawflat, \
                        obs1dome.rawdark, obs1dome.rawflat, obj1dark, \
                        obs1sky.rawsci, obs1sky.procsci, \
                        clobber=False, retmaster=False, verbose=True)

      ::
          import ir
          obs1dome  = ir.initobs('20111214_chip1', cal='dome')
          obs1sky  = ir.initobs('20111214_chip1', cal='sky')
          obs2dome  = ir.initobs('20111214_chip2', cal='dome')
          obs2sky  = ir.initobs('20111214_chip2', cal='sky')
          obj1dark = ['%sMCSA00187%i.fits' % (obs1sky._raw, ii) for ii in range(803,833,2)]
          obj2dark = ['%sMCSA00187%i.fits' % (obs2sky._raw, ii) for ii in range(804,833,2)]

          outfns = [[fn1.replace('proc_M', 'proc_simul_M') for fn1 in obs1sky.procsci], 
                    [fn2.replace('proc_M', 'proc_simul_M') for fn2 in obs2sky.procsci]]
          ir.moircs_cal([obs1sky.rawdark, obs2sky.rawdark], \
                              [obs1sky.rawflat, obs2sky.rawflat], \
                              [obs1dome.rawdark, obs2dome.rawdark], \
                              [obs1dome.rawflat, obs2dome.rawflat], \
                              [obj1dark, obj2dark], \
                              [obs1sky.rawsci, obs2sky.rawsci], \
                              outfns, \
                              clobber=False, retmaster=True, verbose=True, twoatatime=True)
                        

 
    :NOTES:
      - Make 300-sec and 21-sec master dark frames.
  
      - Subtract 300sec dark frame from each 300-sec raw data for sky
        flats
  
      - Make object masks for each raw data. [Not implemented yet!]
  
      - normalize each raw data.
  
      - Median-combine these normalized raw sky data with object
        masks. This is master median sky.
  
      - Make dome flat data.
  
      - Subtract 21-sec master dark from each WASP12 raw data.
  
      - Flatfield them by Dome Flats.
  
      - Measure each sky counts. Multiply the master normalized median
        sky frame by the sky counts (i.e. scale the sky for each data)
        then subtract that sky frame from the dome-flatfielded
        data. Apply the procedure to all frames.
  
      - Subtract the residual global pattern by the plane-fitting
        (imsurfit in IRAF, for example) if necessary.

    As you may already know, making object mask is a critical part of the
    good median-sky data.


                    
    """
    # 2011-12-22 13:16 IJMC: Created
    # 2012-01-01 21:38 IJMC: Added Skysub option
    # 2012-01-02 09:59 IJMC: Added twoatatime option so that sky
    #                        subtraction will be common between both
    #                        chips.

    # Define helper functions:
    def makemasterframe(list_or_array):
        """If an array is passed, return it.  Otherwise, return a
        median stack of the input filename list."""
        if hasattr(list_or_array, 'shape') and len(list_or_array.shape)>1:
            masterframe = np.array(list_or_array, copy=False)
        else:
            masterframe = np.median(map(pyfits.getdata, list_or_array), axis=0)

        return masterframe

    def makedark(flats, dark, badpix=None):
        "Helper function; input is list of FITS filenames and dark (to subtract)."
        nflat = len(flats)
        shape = (pyfits.getval(flats[0], 'naxis1'), \
                     pyfits.getval(flats[0], 'naxis2'))
        masterflat = np.zeros((nflat,) + shape, dtype=np.float32)
        medvalues = np.zeros((nflat, 1, 1), dtype=np.float32)
        for ii in range(nflat):
            masterflat[ii] = pyfits.getdata(flats[ii]) - dark
            if badpix is None:
                medvalues[ii] = np.median(masterflat[ii].ravel())
            else:
                medvalues[ii] = np.median(masterflat[ii][badpix].ravel())

        ## Redefine masterflat, from a stack to a single frame:
        masterflat = np.median(masterflat / medvalues, 0)
        return masterflat

    if twoatatime:
        if badpixelmask is None:
            badpixelmask = [None, None]


    # Load in a bad pixel mask array, if one is not provided:
    if twoatatime:
        for ii in range(2):
            if badpixelmask[ii] is not None and \
                    (not hasattr(badpixelmask[ii], 'shape')):
                badpixelmask[ii] = pyfits.getdata(badpixelmask[ii])
    else:
        if badpixelmask is not None and (not hasattr(badpixelmask, 'shape')):
            badpixelmask = pyfits.getdata(badpixelmask)

    #  - Make 300-sec and 21-sec master dark frames.
    if verbose:  print "Creating master dark frames"
    if twoatatime:
        masterskydark  = [makemasterframe(skydark[ii]) for ii in range(2)]
        masterdomedark = [makemasterframe(domedark[ii]) for ii in range(2)]
        masterobjdark  = [makemasterframe(objdark[ii]) for ii in range(2)]
    else:
        masterskydark  = makemasterframe(skydark)
        masterdomedark = makemasterframe(domedark)
        masterobjdark  = makemasterframe(objdark)

    #  - Subtract 300sec dark frame from each 300-sec raw data for sky
    #    flats

    #  - Make object masks for each raw data.
    print "For now, skipping object masking for flat frames..."

    #  - normalize each raw data.

    #  - Median-combine these normalized raw sky data with object
    #    masks. This is master median sky.
    if verbose:  print "Creating and dark-subtracting master flat frames"
    if twoatatime:
        masterskyflat = []
        for ii in range(2):
            if hasattr(skyflat[ii], 'shape') and len(skyflat[ii].shape)>1:
                masterskyflat.append( np.array(skyflat[ii], copy=False) )
            else:
                masterskyflat.append( makedark(skyflat[ii], masterskydark[ii], \
                                                   badpix=badpixelmask[ii]) )
    else:
        if hasattr(skyflat, 'shape') and len(skyflat.shape)>1:
            masterskyflat = np.array(skyflat, copy=False)
        else:
            masterskyflat = makedark(skyflat, masterskydark, \
                                         badpix=badpixelmask)

    #  - Make dome flat data.
    if twoatatime:
        masterdomeflat = []
        for ii in range(2):
            if hasattr(domeflat, 'shape') and len(domeflat.shape)>1:
                masterdomeflat.append( np.array(domeflat[ii], copy=False) )
            else:
                masterdomeflat.append( makedark(domeflat[ii], masterdomedark[ii], \
                                                    badpix=badpixelmask[ii]) )
    else:
        if hasattr(domeflat, 'shape') and len(domeflat.shape)>1:
            masterdomeflat = np.array(domeflat, copy=False)
        else:
            masterdomeflat = makedark(domeflat, masterdomedark, \
                                          badpix=badpixelmask)

    if verbose:  print "Calibrating and writing raw science frames to disk"
    #  - Subtract 21-sec master dark from each WASP12 raw data.
    #  - Flatfield them by Dome Flats.

    #  - Measure each sky counts. Multiply the master normalized median
    #    sky frame by the sky counts (i.e. scale the sky for each data)
    #    then subtract that sky frame from the dome-flatfielded
    #    data. Apply the procedure to all frames.

    #  - Subtract the residual global pattern by the plane-fitting
    #    (imsurfit in IRAF, for example) if necessary.
    ## 2011-12-22 14:26 IJMC: SURFFIT: TBD

    if twoatatime:
        for rawfile1, procfile1, rawfile2, procfile2 in \
                zip(rawframes[0], procframes[0], rawframes[1], procframes[1]):
            frame0_1 = pyfits.getdata(rawfile1)
            frame0_2 = pyfits.getdata(rawfile2)
            hdr_1 = pyfits.getheader(rawfile1)
            hdr_2 = pyfits.getheader(rawfile2)
            frame1_1 = (frame0_1 - masterobjdark[0]) / masterdomeflat[0]
            frame1_2 = (frame0_2 - masterobjdark[1]) / masterdomeflat[1]

            if skysub:
                skyval = np.median(np.concatenate((frame1_1, frame1_2)).ravel())
                frame1_1 -= masterskyflat[0] * skyval
                frame1_2 -= masterskyflat[1] * skyval
                hdr_1.update('N_SKYSUB', skyval)
                hdr_2.update('N_SKYSUB', skyval)
            else:
                hdr_1.update('N_SKYSUB', 0)
                hdr_2.update('N_SKYSUB', 0)

            hdr_1.update('COMMENT1', 'Flat- and dark-corrected by MOIRCS_cal()')
            hdr_1.update('COMMENT2', '   using TWOATATIME mode.')
            pyfits.writeto(procfile1, frame1_1, header=hdr_1, clobber=clobber)
            hdr_2.update('COMMENT1', 'Flat- and dark-corrected by MOIRCS_cal()')
            hdr_2.update('COMMENT2', '   using TWOATATIME mode.')
            pyfits.writeto(procfile2, frame1_2, header=hdr_2, clobber=clobber)

    else:
        for rawfile, procfile in zip(rawframes, procframes):
            frame0 = pyfits.getdata(rawfile)
            hdr0 = pyfits.getheader(rawfile)
            frame1 = (frame0 - masterobjdark) / masterdomeflat

            if skysub:
                skyval = np.median(frame1.ravel())
                frame1 -= masterskyflat * skyval
                hdr0.update('N_SKYSUB', skyval)
            else:
                hdr0.update('N_SKYSUB', 0)

            hdr0.update('COMMENT1', 'Flat- and dark-corrected by MOIRCS_cal()')
            pyfits.writeto(procfile, frame1, header=hdr0, clobber=clobber)

    if retmaster:
        ret = masterskydark, masterskyflat, masterdomedark, masterdomeflat, masterobjdark
    else:
        ret = None

    return ret


