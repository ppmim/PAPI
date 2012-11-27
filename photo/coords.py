"""Astronomical coordinate functions."""

import re,pdb
import numpy as np
from numpy import arccos,sin,cos
from math import pi

# constants
DEG_PER_HR = 360. / 24.             # degrees per hour
DEG_PER_MIN = DEG_PER_HR / 60.      # degrees per min
DEG_PER_S = DEG_PER_MIN / 60.       # degrees per sec
DEG_PER_AMIN = 1./60.               # degrees per arcmin
DEG_PER_ASEC = DEG_PER_AMIN / 60.   # degrees per arcsec
RAD_PER_DEG = pi / 180.             # radians per degree

def ang_sep(ra1, dec1, ra2, dec2):
    """ Returns the angular separation (in units of degrees) on the
    celestial sphere between two ra/dec coordinates given in degrees.
    Accepts numpy arrays.

    Note: only works for separations larger than about 0.1 arcsec.
    Smaller separations are always returned as zero, because of
    floating point effects.

    >>> np.allclose(ang_sep(2, 0, 4, 0), 2.)
    True
    >>> np.allclose(ang_sep(359, 0, 1, 0), 2.)
    True
    >>> np.allclose(ang_sep(0, 20, 0, -10), 30)
    True
    >>> np.allclose(ang_sep(7, 20, 8, 40), 20.018358)
    True
    >>> np.allclose(ang_sep(7, 20, 250, -50.3), 122.388401)
    True
    >>> ang_sep(-1, 10, 240, -10)           # doctest: +ELLIPSIS
    Traceback (most recent call last):
    ...
    ValueError: RA outside sensible limits. -1
    >>> ras = [24.5,23.6]; decs = [66.89,67.01]
    >>> ref_sep = np.array([3.520032,  3.26675])
    >>> np.allclose(ang_sep(20. ,70., ras, decs), ref_sep)
    True
    >>> ras = [24.5, 23.6]; decs = [91.89, 67.01]
    >>> ang_sep(20.0, 70.0, ras, decs)      # doctest: +ELLIPSIS
    Traceback (most recent call last):
    ...
    ValueError: Dec outside sensible limits. [ 91.89  67.01]
    """ 
    ra1 = np.asarray(ra1);  ra2 = np.asarray(ra2)
    dec1 = np.asarray(dec1);  dec2 = np.asarray(dec2)
    # error checking
    # for ra in ra1,ra2:
    #     if not ((0. <= ra)&(ra < 360.)).all():
    #         raise ValueError("RA outside sensible limits. %s" % ra)
    # for dec in dec1,dec2:
    #     if np.any(np.abs(dec) > 90.):
    #         raise ValueError("Dec outside sensible limits. %s" % dec)

    ra1 = ra1 * RAD_PER_DEG           # convert to radians
    ra2 = ra2 * RAD_PER_DEG
    dec1 = dec1 * RAD_PER_DEG
    dec2 = dec2 * RAD_PER_DEG

    sra1 = sin(ra1);  sra2 = sin(ra2)
    cra1 = cos(ra1);  cra2 = cos(ra2)
    sdec1 = sin(dec1);  sdec2 = sin(dec2)
    cdec1 = cos(dec1);  cdec2 = cos(dec2)

    csep = cdec1*cdec2*(cra1*cra2 + sra1*sra2) + sdec1*sdec2

    # An ugly work-around for floating point issues.
    #if np.any(csep > 1):  print csep
    csep = np.where(csep > 1., 1., csep)

    degsep = arccos(csep) / RAD_PER_DEG
    # only works for separations > 0.1 of an arcsec or  >~2.7e-5 dec
    degsep = np.where(degsep < 1e-5, 0, degsep)
    return degsep


def _dec2s(ra, dec,
          raformat='%02.0f %02.0f %06.3f', decformat='%02.0f %02.0f %05.2f'):
    """
    Converts decimal RA and Dec (both in degrees) to sexigesimal RA
    (hours/minutes/seconds) and Dec (degrees/arcmin/arcsec). Returns
    two strings, RA and Dec.

    doctests:

    >>> dec2s(156.1125638,-10.12986)
    ('10 24 27.015', '-10 07 47.50')
    >>> dec2s(0.0,-90.0)
    ('00 00 00.000', '-90 00 00.00')
    >>> dec2s(148.2,95.0)
    Traceback (most recent call last):
    ...
    ValueError: Decimal RA or Dec outside sensible limits.
    >>> dec2s(360.0,-30.1)
    Traceback (most recent call last):
    ...
    ValueError: Decimal RA or Dec outside sensible limits.
    """
    ra = float(ra)
    dec = float(dec)
    if dec < 0.:
        dec *= -1.
        negdec = True
    else:  negdec = False
    # error checking
    if not (0.0 <= ra < 360.) or dec > 90.:
        raise ValueError("Decimal RA or Dec outside sensible limits.")

    rah, temp = divmod(ra, DEG_PER_HR)
    ram, temp = divmod(temp, DEG_PER_MIN)
    ras = temp / DEG_PER_S
    s_ra = raformat % (rah, ram, ras)

    decd, temp = divmod(dec, 1)
    decm, temp = divmod(temp, DEG_PER_AMIN)
    decs = temp / DEG_PER_ASEC
    if negdec:
        s_dec = '-' + decformat % (decd, decm, decs)
    else:  s_dec = '+' + decformat % (decd, decm, decs)

    return s_ra,s_dec

def dec2s(ra,dec):
    """ Convert and ra and dec (or list of ras and decs) to decimal
    degrees.
    """
    try:
        return _dec2s(ra, dec)
    except TypeError:
        pass
    radec = [_dec2s(r,d) for r, d in zip(ra, dec)]
    return zip(*radec)

def _s2dec(ra,dec):
    """ Converts two strings of sexigesimal RA (hms) and Dec (dms) to
    decimal RA and Dec (degrees).  The separators between h/m/s and
    deg/arcmin/arcsec can be whitespace or colons.  Returns a tuple of
    two floats, (ra, dec).

    doctests:

    >>> s2dec('00:00:00', '90:00:00')
    (0.0, 90.0)
    >>> temp = np.array(s2dec ('10 24 27.015', '-10 07 47.50'))
    >>> reference = np.array([156.1125625,-10.129861111111111])
    >>> np.all(temp - reference < 1.e-10)
    True
    >>> s2dec('25:11:19', '-18:4:88')
    Traceback (most recent call last):
    ...
    ValueError: Either RA or Dec is outside sensible limits.
    RA = 25 11 19, Dec = -18 4 88
    """
    # Convert to floats, noting sign of dec
    ra = re.sub('[:hms]', ' ', ra)
    dec = re.sub('[:dms]', ' ', dec)
    rah,ram,ras = [float(item) for item in ra.split()]
    if dec.lstrip()[0] == '-':
        negdec = True
    else:  negdec = False
    decd,decm,decs = [float(item) for item in dec.split()]
    if negdec:  decd *= -1.
    # Error checking
    if (not 0. <= rah < 24. or not 0. <= ram <= 60. or not 0. <= ras <= 60.
        or decd > 90. or decm >= 60. or decs > 60):
        raise ValueError('Either RA or Dec is outside sensible '
                         'limits.\nRA = %s, Dec = %s' % (ra,dec))
    # Conversion
    d_ra = DEG_PER_HR * rah + DEG_PER_MIN * ram + DEG_PER_S * ras
    d_dec = decd + DEG_PER_AMIN * decm + DEG_PER_ASEC * decs
    if negdec:  d_dec *= -1.

    return d_ra, d_dec

def s2dec(ra, dec):
    """ Convert a sexigesimal ra and dec (or list of ras and decs) to
    decimal degrees.
    """
    try:
        return _s2dec(ra, dec)
    except TypeError:
        pass
    radec = [_s2dec(r,d) for r, d in zip(ra, dec)]
    return map(np.array, zip(*radec))

def match(ra1, dec1, ra2, dec2, tol, allmatches=False):
    """
    match(ra1, dec1, ra2, dec2, tol)

    Given two sets of numpy arrays of ra,dec and a tolerance tol
    (float), returns an array of indices and separations with the same
    length as the first input array.  If and index is > 0, it is the
    index of the closest matching second array element within tol
    arcsec.  If it's -1, then there was no matching ra/dec within tol
    arcsec.

    if allmatches = True, then for each object in the first array,
    return the index and separation of everything in the second array
    within the search tolerance, not just the closest match.

    Note to get the indices of objects in ra2, dec2 without a match, use

    imatch = match(ra1, dec1, ra2, dec2, 2.)
    inomatch = numpy.setdiff1d(np.arange(len(ra2)), set(imatch))

    """
    from numpy.core.records import fromarrays
    
    ra1,ra2,dec1,dec2 = map(np.asarray, (ra1, ra2, dec1, dec2))

    abs = np.abs

    isorted = ra2.argsort()
    sdec2 = dec2[isorted]
    sra2 = ra2[isorted]

    LIM = tol * DEG_PER_ASEC

    match = []
    # use mean dec, assumes decs similar
    decav = np.mean(sdec2.mean() + dec1.mean())
    RA_LIM = LIM / cos(decav * RAD_PER_DEG)

    for ra,dec in zip(ra1,dec1):
        i1 = sra2.searchsorted(ra - RA_LIM)
        i2 = i1 + sra2[i1:].searchsorted(ra + RA_LIM)
        #print i1,i2
        close = []
        for j in xrange(i1,i2):
            if abs(dec - sdec2[j]) > LIM:
                continue
            else:
                # if ras and decs are within LIM arcsec, then
                # calculate actual separation:
                disq = ang_sep(ra, dec, sra2[j], sdec2[j])
                close.append((disq, j))

        close.sort()
        if not allmatches:
            # Choose the object with the closest separation inside the
            # requested tolerance, if one was found.
            if len(close) > 0:
                min_dist, jmin = close[0]
                if min_dist < LIM:
                    match.append((isorted[jmin], min_dist))
                    continue
            # otherwise no match
            match.append((-1,-1))
        else:
            # append all the matching objects
            jclose = []
            seps = []
            for dist,j in close:
                if dist < LIM:
                    jclose.append(j)
                    seps.append(dist)
                else:
                    break
            match.append(fromarrays([isorted[jclose], seps],
                                    dtype=[('ind','i8'),('sep','f8')]))

    if not allmatches:
        # return both indices and separations in a recarray
        temp = np.rec.fromrecords(match, names='ind,sep')
        # change to arcseconds
        temp.sep *= 3600.
        temp.sep[temp.sep < 0] = -1.
        return temp
    else:
        return match

def indmatch(ra1, dec1, ra2, dec2, tol):
    """
    Finds objects in ra1, dec1 that have a matching object in ra2,dec2
    within tol arcsec.

    Returns i1, i2 where i1 are indices into ra1,dec1 that have
    matches, and i2 are the indices into ra2, dec2 giving the matching
    objects.
    """
    m = match(ra1, dec1, ra2, dec2, tol)
    c = m.ind > -1
    i1 = c.nonzero()[0]
    i2 = m.ind[c]
    return i1, i2
    
def unique_radec(ra, dec, tol):
    """ Find unique ras and decs in a list of coordinates.

    RA and Dec must be array sof the same length, and in degrees.

    tol is the tolerance for matching in arcsec. Any coord separated by
    less that this amount are assumed to be the same.

    Returns two arrays.  The first is an array of indices giving the
    first occurence of a unique coordinate in the input list.  The
    second is a list giving the indices of all coords that were
    matched to a given unique coord.

    The matching algorithm is confusing, but hopefully correct and not too
    slow. Potential for improvement...

    Example
    -------

    >>> ra,dec = np.loadtxt('radec.txt.gz', unpack=1)
    >>> iunique, iextra = unique_radec(ra,dec,2)
    >>> iknown, extraknown = np.loadtxt('radec_known.txt.gz', unpack=1)
    >>> np.allclose(iunique, iknown)
    >>> np.allclose(iextra, extraknown)
    """
    matches = match(ra, dec, ra, dec, tol, allmatches=True)
    imatchflat = []
    for m in matches:
        imatchflat.extend(m.ind)
    #pdb.set_trace()
    inomatch = np.setdiff1d(np.arange(len(ra)), list(set(imatchflat)))

    assert len(inomatch) == 0
    # Indices giving unique ra, decs
    iunique = []
    # Will be same length as iunique. Gives all indices in original
    # coords that are matched to each unique coord.
    iextras = []
    assigned = set()
    for j,m in enumerate(matches):
        if not (j % 1000):
            print j
        # get the lowest index in this group
        isort = sorted(m.ind)
        ilow = isort[0]
        if ilow not in assigned:
            iunique.append(ilow)
            assigned.add(ilow)
            iextras.append([ilow])
            # assign any extra indices to this unique coord.
            for i in isort[1:]:
                # check not already been assigned to another coord
                if i not in assigned:
                    iextras[-1].append(i)
                    assigned.add(i)

    return np.array(iunique), iextras