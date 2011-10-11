################################################################################
#
# PANICtool
#
# fileUtils.py
#
# Last update 30/10/2008
#
################################################################################

"""
   File-related utility routine(s) that make life easier.
"""

# Import external modules

import os
import fnmatch

#import misc.paLog
from misc.paLog import log
import pyfits

################################################################################

#log=logging.getLogger('panic.fileutils')

def removefiles(*patterns):

    """
        Remove files with UNIX-like filename pattern matching.
        No questions asked, so be careful!
    """
    
    # 'patterns' is an array or tuple of filename patterns. Can, of course,
    # also be inidividual file names.
    
    # Loop over patterns
    for pattern in patterns:
        # Split off the directory path
        (dirname, filepattern) = os.path.split(pattern)
        #print "dirname= %(dirname)s, file= %(filepattern)s" %vars()
        
        # Get a list of all files in the directory
        try:
            filelist = os.listdir(dirname)
        except OSError, errstr:
            # Succeed even if the directory was not there (and put warning in log)
            log.warning(errstr)
            raise Exception('Cannot list dir %s' % dirname)
        
        # Check each file in the directory list
        for file in filelist:
            # And see if it's name matches the pattern
            if fnmatch.fnmatch(file, filepattern):
                # If yes, (try to) remove it from the system
                #log.debug('Removing file : "%s"' % file)
                try:
                    #print "not removing...debug...."
                    os.remove(os.path.join(dirname, file))
                    #print "removefiles[DEBUG]: file %s removed"%os.path.join(dirname, file)
                except OSError, errstr:
                    # Succeed even if there were no files (and put warning in log)
                    log.error(errstr)
                    raise

################################################################################

def splitMEF_deprecated(fnameMEF, out_filenames):

 """
 
   Split a MEF file into individual files (one per extension).
   The name of the output files will be given in 'output_filenames'
   The new output filename will be 'filenameorig_N.fits' where N goes 1 to NEXTENSIONS

   Input 
      fnameMEF:  filename of MEF
   Output
      out_filenames : list of output filenames 
      Return: the number of extension extracted to outputfiles
 
 
    @deprecated: actually not used anymore ! now we use mef.py module
 """

 next = 0
 
 try:
   hdulist = pyfits.open(fnameMEF)
 except IOError:
   print 'Error, can not open file %s' %(fnameMEF)
   return 0

 try:
   if hdulist[0].header['EXTEND']!=True:
     print 'Error, file %s is not a MEF file' %(fnameMEF)
     return 0
 except KeyError:
   print 'Error, file %s is not a MEF file' %(fnameMEF)
   return 0

 try:
   next=hdulist[0].header['NEXTEND']
 except KeyError:
   print 'Error, card NEXTEND not found'

 for i in range(1,next+1):
   sufix="_%d.fits" %i
   out_filenames.append(fnameMEF.replace('.fits', sufix))
   out_hdulist = pyfits.HDUList([pyfits.PrimaryHDU(header=hdulist[i].header, data=hdulist[i].data)])
   out_hdulist.verify('silentfix')
   out_hdulist.writeto(out_filenames[i-1], output_verify='ignore', clobber=True)
   out_hdulist.close(output_verify='ignore')
   del out_hdulist
   print "File %s created " %(out_filenames[i-1])

 
 return next 
 print "End of createMEF"
  
################################################################################  
def linkSourceFiles( source, dest ):
    """Create a symbolic link to all the sources specified in 'source', which can be a file a dir"""

    if (type(source)==type(list())):
        #we suppose is a list
        for file in source:
            #print "FILE=", file
            if not os.path.exists(dest+"/"+os.path.basename(file)):
                os.symlink(file, dest+"/"+os.path.basename(file))
    elif (os.path.isfile(source)):
        #We have a source-file with absolute path for the data files
        file=open(source, 'r')
        for line in file:
            #print "FILE=", line
            if line.endswith('\n'):
                line=line.replace('\n','')
            print
            #print "DEST=", dest+"/"+os.path.basename(line)
            if not os.path.exists(dest+"/"+os.path.basename(line)):
                os.symlink(line, dest+"/"+os.path.basename(line))
    elif (os.path.isdir(source)):
        #We have a source-directory as input data
        for file in glob.glob(source+"*.fits"):
            print "FILE=", file
            if not os.path.exists(dest+"/"+os.path.basename(file)):
                os.symlink(file, dest+"/"+os.path.basename(file))
    else:
        print "Error reading source; cannot indentify the type of input"
