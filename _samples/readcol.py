import string
from numpy import asarray
from scipy.stats import mode

#
#
#Reading columns of ASCII text is one of the most basic functions a programmer can need. For example, reading the output of a FORTRAN code.
#
#IDL has the flexible and powerful readcol.pro from Goddard Space Flight Center's astrolib. Python has no equivalent, as far as I can tell. #This code is meant to be the pythong equivalent of readcol.
#
#I don't claim that it's perfectly flexible at the moment, but it's fast, and I'm willing to work to improve it. Originally I had used #readline() but this is terribly slow for large files; now readlines() is used. See below.
#
#

def readcol(filename,skipline=0,names=False,dtype='float',fsep=None,twod=True,comment='#',verbose=True):
    """
    The default return is a two dimensional float array.  You can specify
    the data type (e.g. dtype='S') in the normal python way.  If you want
    a list of columns output instead of a 2D array, pass 'twod=False'.
   
    Example usage:
    CASE 1) a table has the format:
     X    Y    Z
    0.0  2.4  8.2
    1.0  3.4  5.6
    0.7  3.2  2.1
    ...
    names,[x,y,z]=readcol("myfile.tbl",names=True,twod=False)
    or
    x,y,z=readcol("myfile.tbl",skipline=1,twod=False)
    or
    names,xx = readcol("myfile.tbl",names=True)

    CASE 2) no title is contained into the table, then there is
    no need to skipline:
    x,y,z=readcol("myfile.tbl")
   
    CASE 3) there is a names column and then more descriptive text:
     X      Y     Z
    (deg) (deg) (km/s)
    0.0    2.4   8.2
    1.0    3.4.  5.6
    ...
    then use:
    names,x,y,z=readcol("myfile.tbl",names=True,skiplines=1,twod=False)
    or
    x,y,z=readcol("myfile.tbl",skipline=2,twod=False)

    INPUTS:
        fsep - field separator, e.g. for comma separated value (csv) files
        skiplines - number of lines to ignore at the start of the file
        names - read / don't read in the first line as a list of column names
        dtype - datatype of numpy array
        twod - two dimensional or one dimensional output
    """
    f=open(filename,'r').readlines()
    if names is True:
        nms=f.pop(0).split(fsep)

    null=[f.pop(0) for i in range(skipline)]
   
    fstrip = map(string.strip,f)
    fseps = [ fsep for i in range(len(f)) ]
    splitarr = map(string.split,fstrip,fseps)
    for i in xrange(splitarr.count([''])):
        splitarr.remove([''])

    # check to make sure each line has the same number of columns to avoid
    # "ValueError: setting an array element with a sequence."
    nperline = map(len,splitarr)
    ncols,nrows = mode(nperline)
    if nrows != len(splitarr):
        if verbose:
            print "Removing %i rows that don't match most common length.  \
             \n%i rows read into array." % (len(splitarr) - nrows,nrows)
        for i in xrange(len(splitarr)-1,-1,-1):  # need to go backwards
            if nperline[i] != ncols:
                splitarr.pop(i)

    # remove comment lines
    if comment != None:
        splitarr = filter(lambda a: comment.find(a[0][0]),splitarr)

    try:
        x = asarray( splitarr , dtype=dtype)
    except:
        if verbose:
            print "WARNING: reading as string array because %s array failed" % dtype
        x = asarray( splitarr , dtype='S')

    if names is True:
        if twod:
            return nms,x
        else:
            return nms,[ x.T[i] for i in xrange(x.shape[1]) ]
    else:
        if twod:
            return x
        else:
            return [ x.T[i] for i in xrange(x.shape[1]) ]
