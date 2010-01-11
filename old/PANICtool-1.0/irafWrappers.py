################################################################################
#
# FIEStool
#
# irafWrappers.py
#
# Last update 29/09/2005
#
################################################################################

"""
   A collection of short routines that wrap around the PyRAF calls to IRAF
   routines. Apart from making the PyRAF call, these routines take care of
   everything surrounding the call, such as copying the calibration files,
   preparing the task using parameter definition files, as well as
   the final cleaning of temporary files.
"""

# Import external modules

# System modules

import os
import config
import string
import irafDbNames


# PyRAF modules that relate to echelle data reduction
from pyraf import *
from iraf import noao
from iraf import imred
from iraf import ccdred
from iraf import echelle
from iraf import astutil

# Additional module for workaround with PyRAF on slow displays
from time import sleep

################################################################################

# Upon loading this module, set the default imtype to 'fits' (IMPORTANT!)
iraf.set(imtype='fits')

################################################################################

def _tempname(filename):

  """
     Creates valid temporary file name, even if filename contains full path
  """

  newname = os.path.join(os.path.dirname(filename),
                        'TMP_'+os.path.basename(filename))

  return newname


################################################################################

def _getparlistname(taskname):

  """
     Determine the path to the IRAF parameter definition files for the
     current reduction mode. (There is one set for each mode).
  """

  parlistname = os.path.join(config.options['iraf_taskconfig_dir'].value,
                             config.options['currentmode'].value,
			     taskname+'.par')

  if not os.path.exists(parlistname) :
    raise IOError, 'IRAF task configuration file not found: %s\n(Check your taskconf directory setting)' % parlistname
    

  return parlistname


################################################################################

def _copy_apref(file, refframe):

  """
     Copy the aperture definition from the reference database to the
     outdata working directory
  """

  import shutil

  # Get the name of the aperture definition file in the IRAF database
  apfile = irafDbNames.apname(refframe)

  # And separate out the source directory
  refbase, apfile = os.path.split(apfile)
 
  # Separate out the target directory
  filebase, file = os.path.split(file)
  filebase       = os.path.join(filebase, 'database')

  # Get the full path names for the source and target files
  origfile = os.path.join(refbase,  apfile)
  targfile = os.path.join(filebase, apfile)

  # Make sure the target directory exists; if not - create
  if not os.path.isdir(os.path.dirname(targfile)):
    os.makedirs(os.path.dirname(targfile))

  # Copy the file (as long as source and target are different)
  if (refbase != filebase):
    # Always succeed with this action!
    try:
      shutil.copy(origfile, targfile)
    except IOError: pass


################################################################################

def _copy_ecref(file, refframe):

  """
     Copy the wavelength definition AND the master wavelength definition frame
     from the reference database to the outdata working directory
  """

  import shutil

  # file is actually the name of the file being reduced, and therefore
  # the directory of this file is the only really interesting part

  # First step : separate out the source and target directories
  filebase, dummy  = os.path.split(file)
  refbase, refname = os.path.split(refframe)
  
  # Second step : copy the master wavelength definition file
  targfile = os.path.join(filebase, refname)
  if (refbase != filebase):
    shutil.copy(refframe, targfile)

  # Third step : copy the corresponding wavelength definition file

  # Get the name of the wavelength definition file in the IRAF database
  ecfile = irafDbNames.ecname(refframe)
  
  # And separate out the source directory
  filebase, dummy  = os.path.split(file)
  ecbase, ecname   = os.path.split(ecfile)

  # Append IRAFs database directory to target dir
  filebase        = os.path.join(filebase, 'database')

  # Get the full path names for the target file
  targfile = os.path.join(filebase, ecname)

  if not os.path.isdir(os.path.dirname(targfile)):
    os.makedirs(os.path.dirname(targfile))

  # Now do it: copy the wavelength definition file
  if (ecbase != filebase):
    # Always succeed with this action!
    try:
      shutil.copy(ecfile, targfile)
    except IOError: pass



################################################################################

def _dbclean(file):

  """
     Remove all existing files in the (IRAF-)database related to a single
     frame, in order to clean the table before a new reduction. This is
     especially useful for the wavelength definition files, because IRAF only
     appends new definitions to the existing ones, generating very large files
     that eventually leads to non-functional wavelength definition files.
  """

  import fileUtils

  # Separate out the directory and basename of the frame to clean up for
  filedir, file = os.path.split(file)
  filebase, ext = os.path.splitext(file)

  # Create pattern of files to remove
  zapfiles = os.path.join(filedir, 'database', '??'+filebase+'*')

  # Get rid of them!  
  fileUtils.removefiles(zapfiles)


################################################################################

def scattering(infile, outfile, reffile):

  """
     Wrapper for calling the task 'echelle.apscatter'
  """

  # Copy the aperture reference from the database
  _copy_apref(infile, reffile)

  # Fill IRAF parameters with the settings from 'apscat1' and 'apscat2', which
  # are tasks dedicated to setting parameters
  iraf.apscat1.setParList(ParList = _getparlistname('apscat1'))
  iraf.apscat2.setParList(ParList = _getparlistname('apscat2'))

  # Separate out the source and target directory and file name, as well as
  # the name of the (copied) aperture reference file
  base, infile   = os.path.split(infile)
  dummy, reffile = os.path.split(reffile)

  # Change to the output directory
  iraf.chdir(base)

  # Call the echelle.apscatter task through PyRAF
  iraf.apscatter(input=infile, output=outfile,
                 references=reffile,
		 ParList = _getparlistname('apscatter')
		 )

  # Change back to the original working directory
  iraf.chdir()


################################################################################

def normalize(infile, outfile, reffile, threshold):

  """
     Wrapper for calling the task 'echelle.apflatten'
  """

  # Copy the aperture reference from the database
  _copy_apref(infile, reffile)

  # Following line is necessary to avoid bug (?) in PyRAF
  # Problem : 'apflat' is referenced instead of 'apflat1'
  # Solution: Generate 'apflat', identical to 'apflat1'
  _dummy = iraf.IrafTaskFactory(taskname='apflat',
                value='apflat1.par', function=iraf.echelle.apflat1)

  # Now set the fitting parameters
  iraf.apflat1.setParList(ParList = _getparlistname('apflat1'))  

  # Separate out the source and target directory and file name, as well as
  # the name of the (copied) aperture reference file
  base, infile   = os.path.split(infile)
  dummy, reffile = os.path.split(reffile)

  # Change to the output directory
  iraf.chdir(base)

  # Call the echelle.apflatten task through PyRAF
  iraf.apflatten(input=infile, output=outfile,
                 references=reffile,
		 threshold=threshold,
		 ParList = _getparlistname('apflatten')
		 )

  # Change back to the original working directory
  iraf.chdir()


################################################################################

def fit1d(infile, outfile):

  """
     Wrapper for calling the task 'echelle.fit1d'
  """

  # Separate out the directory part of the target directory
  base, outfile = os.path.split(outfile)

  # Change to the output directory
  iraf.chdir(base)

  # Call the echelle.fit1d task through PyRAF
  iraf.fit1d(input=infile, output=outfile,
             ParList = _getparlistname('fit1d')
	     )

  # Change back to the original working directory
  iraf.chdir()


################################################################################

def findtraceord(infile, trace='yes', reference=""):


  """
     Wrapper for calling the task 'echelle.apfind' and 'echelle.aptrace'
  """

  # Separate out the reference directory
  base, infile = os.path.split(infile)

  # Change to the reference directory
  iraf.chdir(base)

  # Set the default aperture parameters    
  iraf.apdefault.setParList(ParList = _getparlistname('apdefault'))

  # Configure the default answers to questions posed in the IRAF task
  # (Unfortunately, it is not possible to do this with a parameter file
  #  because it won't remember the capitalized answers!)
  iraf.apparams.initialize='no'
  iraf.apparams.ansfind='YES'
  iraf.apparams.ansrecenter='NO'
  iraf.apparams.ansresize='NO'
  iraf.apparams.ansedit='YES'
  iraf.apparams.ansdbwrite='YES'

  # Call the echelle.apfind task through PyRAF

  # Modified code (workaround) - thanks John!
  gwm.window("The incredible FIES-tool !      ApFind")
  sleep(1.5)
  iraf.apfind.setParList(ParList = _getparlistname('apfind'))
  sleep(1.5)
  iraf.apfind(input=infile, references=reference)

  # Old code - may lock up PyRAF window on slow workstations
  # iraf.apfind(input=infile, references=reference,
  #	      ParList = _getparlistname('apfind')
  #	      )

  # It is possible to set 'trace' to 'no' and skip this part, but this
  # is currently not used
  if trace == 'yes':

    # Configure the default answers to questions posed in the IRAF task
    # (see remarks above)
    iraf.apparams.ansfind='NO'
    iraf.apparams.ansedit='NO'
    iraf.apparams.anstrace='YES'
    iraf.apparams.ansfittrace='NO'
    iraf.apparams.ansfittrace1='NO'

    # Call the echelle.apfind aptrace through PyRAF
    iraf.aptrace(input=infile,
        	 ParList = _getparlistname('aptrace')
		 )

  # Reset the apparms.initialize parameter to its default (maybe not
  # necessary, but try to be nice anyway)
  iraf.apparams.initialize='yes'

  # Change back to the original working directory
  iraf.chdir()

  # Delete the current graphics window
  gwm.delete(gwm.getActiveWindowName())

################################################################################

def findwavesol(infile, keep_old_solution=0):

  """
     Wrapper for calling the task 'echelle.ecidentify'
  """

  # At this point, there is a 'bug' in IRAF that needs to be circumvented. The
  # problem is that'ecidentify'  uses '/' as a directory separator when
  # writing to the database directory, instead of '_' (as is the case for
  # aperture definitions). This makes it impossible to use 'ecidentify' on
  # frames outside the current directory. As this is hard-coded in the IRAF
  # source, the only way to avoid this is to perform the line identification
  # of the wavelength definition in its own directory, and then copy the
  # wavelength definition from the source to the output/reference directory
  # when applying the wavelength definition on data frames

  # Clean the database for entries related to this frame
  if not keep_old_solution: _dbclean(infile)

  # Separate out the source directory
  base, infile = os.path.split(infile)

  # Change to the source directory
  iraf.chdir(base)

  # Call the echelle.ecidentify task through PyRAF

  # Modified code (workaround) - thanks John!
  gwm.window("The incredible FIES-tool !    EcIdentify")
  sleep(2.5)
  iraf.ecidentify.setParList(ParList = _getparlistname('ecidentify'))
  sleep(2.5)
  iraf.ecidentify(images=infile)

  # Old code - may lock up PyRAF window on slow workstations
  # iraf.ecidentify(images=infile,
  #               ParList = _getparlistname('ecidentify')
  #	          )

  # Change back to the original working directory
  iraf.chdir()

  # Delete the current graphics window
  gwm.delete(gwm.getActiveWindowName())

################################################################################

def refindwavesol(infile, reffile):

  """
     Wrapper for calling the task 'echelle.ecreidentify'
  """
  
  # (the comment mentioned in 'findwavesol' is also relevant here)

  # Clean the database for entries related to this frame
  _dbclean(infile)

  # Copy the wavelength definition from the source to the reference directory
  _copy_ecref(infile, reffile)

  # Separate out the source and reference directory
  base, infile   = os.path.split(infile)
  dummy, reffile = os.path.split(reffile)
  reffile, dummy = os.path.splitext(reffile)

  # Change to the source directory
  iraf.chdir(base)

  # Call the echelle.ecreidentify task through PyRAF
  iraf.ecreidentify(images=infile,
                    reference=reffile,
		    ParList = _getparlistname('ecreidentify')
	            )

  # Change back to the original working directory
  iraf.chdir()
 
################################################################################

def findwaveshift(infile, reffile):

  """
     Wrapper for calling the task 'echelle.ecreidentify' and return (!) the
     reported pixel and wavelength shift. Useful for determining and correcting
     the instrumental drift.
  """

  # Copy the wavelength definition from the reference to the input directory
  _copy_ecref(infile, reffile)

  # Separate out the source and reference directory
  base, infile   = os.path.split(infile)
  dummy, reffile = os.path.split(reffile)
  reffile, dummy = os.path.splitext(reffile)

  # Change to the source directory
  iraf.chdir(base)

  # Call the echelle.ecreidentify task through PyRAF, and store the
  # stdout stream in 'output'
  output = iraf.ecreidentify(images=infile,
                    reference=reffile,
		    ParList = _getparlistname('ecreidentify'),
		    Stdout=1
	            )

  # Print the generated output to the screen
  print output

  # Get the forth line, containing the values of the shifts
  shiftline = string.split(output[4])

  # Determine the ratio of the number of reidentified lines and the
  # number of defined lines
  n1, n2   = string.split(shiftline[1], '/')
  ratio   = float(n1) / float(n2)

  # Extract the shift in pixel and in wavelength
  pixshift = float(shiftline[3])
  wvlshift = float(shiftline[4])

  # Change back to the original working directory
  iraf.chdir()

  # Return a tuple with three output values
  return (pixshift, wvlshift, ratio)
 
################################################################################

def extractspec(infile, outfile, reffile):

  """
  Wrapper for calling the task 'echelle.apsum'
  """

  # Clean the database for entries related to this frame
  _dbclean(outfile)

  # Copy the aperture definition from the reference to the input directory
  _copy_apref(infile, reffile)

  # Separate out the source and reference directory
  base, infile   = os.path.split(infile)
  dummy, reffile = os.path.split(reffile)

  # Change to the source directory
  iraf.chdir(base)

  # Call the echelle.apsum task through PyRAF
  iraf.apsum(input=infile,
             output=outfile,
 	     references=reffile,
	     ParList = _getparlistname('apsum')
	     )

  # Change back to the original working directory
  iraf.chdir()


################################################################################

def extractwave(infile, outfile, reffile):

  """
  Wrapper for calling the task 'echelle.apsum', optimized for ThAr frames
  (which means that spectrum extraction is kept as simple as possible)
  """

  # Clean the database for entries related to this frame
  _dbclean(outfile)

  # Copy the aperture definition from the reference to the input directory
  _copy_apref(infile, reffile)

  # Separate out the source and reference directory
  base, infile   = os.path.split(infile)
  dummy, reffile = os.path.split(reffile)

  # Change to the source directory
  iraf.chdir(base)

  # Call the echelle.apsum task through PyRAF
  iraf.apsum(input=infile,
             output=outfile,
 	     references=reffile,
	     ParList = _getparlistname('apsum_wave')
	     )

  # Change back to the original working directory
  iraf.chdir()


################################################################################

def specplot(infile, order=0):

  """
     Wrapper for calling the task 'echelle.splot' - the plotting window!
  """

  # order gives the number of the order that is displayed initially

  # Separate out the source directory
  base, infile = os.path.split(infile)

  # Change to the source directory
  iraf.chdir(base)

### added 20/12/2006 by JHT
  gwm.window("FIES-tool !      IRAF splot")
  sleep(1.5)
###

  # Call the echelle.splot task through PyRAF
  iraf.splot(images=infile, line=order, mode='h')

  # Change back to the original working directory
  iraf.chdir()

  # Delete the current graphics window
  gwm.delete(gwm.getActiveWindowName())


################################################################################

def dispcor(infile, outfile, reffile):

  """
     Wrapper for calling the task 'echelle.refspec' and 'echelle.dispcor'
  """

  # Just ignore and exit if no wavelength reference exists
  if not os.path.exists(irafDbNames.ecname(reffile)):
    return

  # Copy the wavelength definition from the reference to the input directory
  _copy_ecref(infile, reffile)

  # Separate out the source and reference directory
  base, infile   = os.path.split(infile)
  refpath, reffile = os.path.split(reffile)
  refroot, dummy = os.path.splitext(reffile)

  # Change to the source directory
  iraf.chdir(base)

  # Call the echelle.refspec task through PyRAF
  iraf.refspec(input=infile,
               references=reffile,
               ParList = _getparlistname('refspec')
	       )

  # Call the echelle.dispcor task through PyRAF
  iraf.dispcor(input=infile, output=outfile,
               ParList = _getparlistname('dispcor')
	       )

  # Change back to the original working directory
  iraf.chdir()



################################################################################

def shiftorder(file, order, newoffset):

  """
     Wrapper for calling the task 'echelle.sapertures' to change (offset)
     the base order number
  """

  # Separate out the source directory
  base, infile = os.path.split(file)

  # Change to the source directory
  iraf.chdir(base)

  # Call the echelle.sapertures task through PyRAF
  iraf.sapertures(input=file,
		  apertures=order,
		  w1=newoffset,
                  ParList = _getparlistname('sapertures')
	         )

  # Change back to the original working directory
  iraf.chdir()


################################################################################

def mergeorders(infile, outfile, reffile):

  """
     Wrapper for calling the task 'echelle.scombine'
  """

  # Copy the wavelength definition from the reference to the input directory
  _copy_ecref(infile, reffile)

  # Separate out the source directory
  base, infile   = os.path.split(infile)

  # Change to the source directory
  iraf.chdir(base)

  # Call the echelle.scombine task through PyRAF
  iraf.scombine(input=infile,
		output=outfile,
		ParList = _getparlistname('scombine')
		)

  # Change back to the original working directory
  iraf.chdir()


################################################################################

def rvcorrect(year, month, day, ut, ra, dec, epoch):

  """
     Wrapper for calling the task 'astutil.rvcorrect'
  """

  # Call the echelle.scombine task through PyRAF
  iraf.rvcorrect(year=year,
                 month=month,
		 day=day,
		 ut=ut,
		 ra=ra,
		 dec=dec,
		 epoch=epoch,
		 observatory='lapalma'
#		 ParList = _getparlistname('rvcorrect')
		)

  vhelio = iraf.rvcorrect.vhel

  # Change back to the original working directory
  iraf.chdir()

  return  vhelio


################################################################################

def darkcombine(infile, outfile ):

  """
     Wrapper for calling the task 'noao.imred.ccdred.darkcombine'
  """

  # Separate out the source directory
  #base, infile   = os.path.split(infile)

  print infile

  for nfile in infile:
    list_of_files.append (nfile)

  print list_of_files
     
  
  
  # Change to the source directory
  iraf.chdir('/disk-a/caha/panic/TMP/data')

  # Call the noao.imred.ccdred task through PyRAF
  iraf.darkcombine(input='/disk-a/caha/panic/TMP/data/A0408060036.fits, /disk-a/caha/panic/TMP/data/A0408060037.fits',
		output=outfile,
		ParList = _getparlistname('darkcombine')
		)

  # Change back to the original working directory
  iraf.chdir()


################################################################################

