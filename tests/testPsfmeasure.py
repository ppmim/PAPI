from pyraf import iraf
#import pyraf.iraf as iraf
#from iraf import obsutil
import re


def getAverageFWHMfromPsfmeasure(image, coord_file):
    """
    Calculate the average Full Width Half Max for the objects in image
    at the coords specified in coord_file calling iraf.obsutil.psfmeasure.

    The coordinates in coord_file should be in the same world coordiantes
    as the WCS applied to the image.
    
    Exam. coord_file:
    1024  1024
    
    """

    iraf.noao(_doprint=0)
    iraf.obsutil(_doprint=0)
    iraf.unlearn("psfmeasure")
  
    psfmeasure = iraf.obsutil.psfmeasure
    # setup all paramaters
    psfmeasure.coords = "mark1"
    psfmeasure.wcs = "world"
    psfmeasure.display = "no"
    psfmeasure.size = "MFWHM"
    psfmeasure.radius = 5
    psfmeasure.iterations = 2
    psfmeasure.imagecur = coord_file
    psfmeasure.graphcur = '/dev/null' #file that is empty by definition
    psfmeasure.logfile = "fwhm.log"
    res = psfmeasure(image, Stdout=1)[-1] # get last linet of output
    numMatch = re.compile(r'(\d+(\.\d+)?)')
    match = numMatch.search(res)

    return float (match.group(1))

   
#myfile = "/home/panic/pruebaF.fits"
myfile = "/home/panic/prueba.fits"
myfile = "/home/panic/DATA/Foco/focus_H0007.fits"
coord_file = "/tmp/coord_file.txt"

fwhm = getAverageFWHMfromPsfmeasure(myfile, coord_file)

print "FWHM = ",fwhm
