#!/bin/env python

"""calibrate IRAF phot output text file (.mag) with standard star table"""

import os, sys
from math import sqrt, cos, log10
from astLib import astWCS

###################

def getirafmags(filename, nsigma=3):

    #nsigma = 3

    f=open(filename, 'r')
    obj = 0
    l = 0
    ras = []
    decs = []
    mags = []
    errs = []
    lims = []
    nowcs = True
    fitsfile = ''
    for line in f:
        if line[0]=='#':
            if line[0:5]=='#K AP':
                naperture = line.count(',')+1
            if line[0:8]=='#K EPADU':
                epadu = float(line[16:21])
            if line[0:7]=='#K ZMAG':
                zmag = float(line[16:21])
            continue
        list = line.split()
        if l==0:
            oldfitsfile = fitsfile
            fitsfile = list[0] #####'../' +
            if (fitsfile != oldfitsfile):
                #print fitsfile
                try:
                    WCS = astWCS.WCS(fitsfile)
                    nowcs = False
                except:
                    nowcs = True
            obj = list[len(list)-2]
        if l==1:
            mag = []
            err = []
            lim = []
            x=list[0]
            y=list[1]
            if nowcs == False:
                [ra, dec] = WCS.pix2wcs(float(x), float(y))
            else:
                [ra, dec] = [x, y]
            ras.append(ra)
            decs.append(dec)
        if l==2:
            try:
                stdev = float(list[1])
                nsky = float(list[3])
            except:
                stdev = -1
                nsky = -1
        if l==3:
            ne = len(list)
            itime=float(list[0])
            air=list[1]
            filt=(list[2])[0]
            #ut=list[3]
        if l>=4 and l<=3+naperture:
            area = float(list[2])
            flux = float(list[3])
            if (flux < 0): flux = 0.0
            skyfluxerr = sqrt(area*stdev**2 + area**2 * stdev**2 / nsky)
            fluxnsigma = flux*0 + nsigma*skyfluxerr #"image" mag limit, ignores the real flux here
            if (flux < nsigma*skyfluxerr):
                magnsigma = zmag - 2.5*log10(fluxnsigma) + 2.5*log10(itime)
                mag.append(magnsigma)
                err.append(0)
                lim.append(True)
            elif (fluxnsigma <= 0.0 or nsky == -1 or stdev == -1):
                mag.append(0)
                err.append(0)
                lim.append(True)
            else:
                mag.append(float(list[4]))
                err.append(float(list[5]))
                lim.append(False)
        if l==3+naperture:
            #print x.ljust(8), y.ljust(8), obj, filt, air.ljust(5),
            mags.append(mag)
            errs.append(err)
            lims.append(lim)
            l = -1
        l=l+1
    f.close()
    
    
    return ras, decs, mags, errs, lims

#####################


def autocal(magfilename,refstarfilename,calcol=2,nsigma=3):

    magfile = open(magfilename,'r')
    calfile = open(refstarfilename,'r')

    cal = calfile.readlines()
    calfile.close()
    
    print '# input file:', magfilename
    print '# calibration file:', refstarfilename
    
    [ras, decs, mags, errs, lims] = getirafmags(magfilename,nsigma=nsigma)
    nap = len(mags[0])
    nobj = len(ras)
    
    outlier = [False] * nobj
    magdiff = []
    magdiffgood = []
    
    #ignore leading comments and blank lines
    head = 0
    while len(cal[head]) <= 1 or cal[head][0] == '#': head = head + 1
    
    cal = cal[head:]
    lc = 0  #line
    obj = 0 #object number - blanks, comments skipped
    refmags = []
    
    while lc < len(cal):
        #in-table comments, blank lines, no magnitude cause this source to be excluded from the median sample
        listcal =  cal[lc].split()
        if (len(cal[lc]) <= 1 or cal[lc][0] == '#'):
            lc += 1
            continue
        
        if (len(listcal) > calcol): 
            refmag = listcal[calcol]
        else:
            refmag = '-'
        
        refmags.append(refmag)
        if (refmag[0] != '-' and refmag[0] != '?' and refmag[0] != 'x' and refmag[0] != 'ERROR'):
            for a in range(len(mags[obj])):
                magdiff.append([])
                magdiff[a].append(float(refmag)-mags[obj][a])
	lc += 1
        obj += 1
    
    if nobj < obj: print '#Warning: #photobjects ('+str(nobj)+') < #refobjects *' +str(obj)+')'
    
    if obj < nobj: 
       while obj < nobj:
          refmags.append(['-']) #Need to append several for different apertures?
          obj += 1
    
    
    med = [0] * nap
    std = [0] * nap
    for a in range(nap):
        med[a] = median(magdiff[a])
        std[a] = stdev(magdiff[a])
    
    #flag outliers
    obj = 0
    refap = 0
    for refmag in refmags:
        if refmag[0] != '-' and refmag[0] != '?' and refmag[0] != 'x' and refmag[0] != 'ERROR':
            isoutlier = abs(float(refmag)-mags[obj][refap] - med[refap]) > min(max(0.3, std[refap]),0.6) or float(refmag) == 0.0
            outlier[obj] = isoutlier
            if not isoutlier: 
                for a in range(len(mags[obj])):
                    magdiffgood.append([])
                    magdiffgood[a].append(float(refmag)-mags[obj][a])
        obj = obj + 1
    
    for a in range(nap):
        med[a] = median(magdiffgood[a])
        std[a] = stdev(magdiffgood[a])
    
    print '# zeropoint adjustment: ', 
    
    for a in range(nap):
        print str(med[a]).ljust(5)[0:5], 
        print '+/-', str(std[a]).ljust(5)[0:5],
        if (a < nap-1): 
            print ', ',
	else:
            print    
    
    zeropt = med 
        
    obj = 0
    
    for refmag in refmags:
        print str(ras[obj]).ljust(10)[0:10], str(decs[obj]).ljust(9)[0:9], 
        if refmag[0] != '-' and refmag[0] != '?' and refmag[0] != 'x':
            print refmag.ljust(6)[0:6], 
        else:
            print ' '*6,
        for a in range(len(mags[obj])):
            obsmag = mags[obj][a]
            calobsmag = obsmag + zeropt[a]
            if not lims[obj][a]:
                print str(calobsmag).ljust(6)[0:6], str(errs[obj][a]).ljust(6)[0:6],
            else:
                print '>'+str(calobsmag).ljust(5)[0:5], ' '*6,
        if outlier[obj]:
            print 'x'
        else:
            print
        obj = obj + 1

    return        

##############

def median(inlist):
    sinlist = sorted(inlist)
    nm = len(sinlist)
    if (nm % 2) == 1:
        med = float(sinlist[nm/2])
    else:
        med = (float(sinlist[nm/2-1]) + float(sinlist[nm/2]))/2.0
    
    return med

def stdev(inlist):
    tot = 0
    var = 0
    for s in inlist:
        tot = tot + s
    mean = tot / len(inlist)
    for s in inlist:
        var = var + (s-mean)**2
    var /= len(inlist)
    return sqrt(var)


######################################################################
def usage():
    (xdir,xname)=os.path.split(sys.argv[0])
    print "Usage:  %s magfilename reffilename calcol" % xname

######################################################################
def main():

    argv = sys.argv
    
    if (len(argv)<3):
        usage()
        sys.exit(1)

    if (len(argv)==3):
        autocal(argv[1], argv[2])
    elif (len(argv)==4):
        autocal(argv[1], argv[2], calcol=int(argv[3]))
    elif (len(argv)>=5):
        autocal(argv[1], argv[2], calcol=int(argv[3]), nsigma=float(argv[4]))
    #elif (len(argv)>=5):
    #    autocal(argv[1], argv[2], calcol=int(argv[3]), ap=int(argv[4]-1))


######################################################################
if __name__=='__main__':
    main()
