#! /usr/bin/env python

""" 
    Run Get the unique values of [read_mode, itime, ncoadd, save_mode]
    for the files of a given directory.
"""

# Author: Jose M. Ibanez (c) 2014
# Email: jmiguel@iaa.es
# License: GNU GPLv3

from optparse import OptionParser
import sys
import os
import shutil
import glob
import astropy.io.fits as fits


  
def getItimesNcoadds(path, output_file, recursive=False):
    """
    Read all the FITS file of the given path and create a 
    output_file with a table with the unique combinations
    of [read_mode, itime, ncoadds, save_mode].
    
    - DARK and DOME_LAMP_[ON|OFF] images are ignored.
    - Only lir and rrr-mpia images are recognized.
    - Cubes (SEF or MEF) are detected and the number of layers is used as NCOADDS.
    
    Example of output file:
    
    lir       2.7   2  mef
    lir       5.0  10  mef
    rrr-mpia  1.3   5  mef
    rrr-mpia  1.3   5  sef
    
    
    """
    
    # Read the files
    if os.path.exists(path):
        if os.path.isfile(path):
            try:
                hdulist = fits.open(path)
                filelist = [path]
            except:
                filelist = [line.replace( "\n", "") 
                            for line in fileinput.input(path)]
        elif os.path.isdir(path):
            filelist = glob.glob(path + "/*.fit*")
            # Look for subdirectories
            if recursive:
                subdirectories = [ name for name in os.listdir(path) \
                    if os.path.isdir(os.path.join(path, name)) ]
                for subdir in subdirectories:
                    filelist += glob.glob(os.path.join(path, subdir)+"/*.fit*")
    else:
        print("Error, file %s does not exits" % path)
        return 0
     
    file_types = []
    fd = open(output_file + "_tmp_", "w+")
    fd.write("# READMODE\tITIME\tNCOADDS\tSAVEMODE\n")
     
    for my_file in filelist:
        try:
            my_fits = fits.open(my_file)
            # Image type
            if 'IMAGETYP' in  my_fits[0].header:
                papitype = my_fits[0].header['IMAGETYP']
            else:
                papitype = 'unknown'
                
            # Ignore DARKs and DOME_LAMP_[ON|OFF] images
            if "DARK" in papitype  or "LAMP" in papitype:
                print("Image is a DARK or DOME_LAMP_. Skipping file %s" % my_file)
                continue
            
            # Read-Mode
            read_mode = my_fits[0].header['READMODE']
            if read_mode == 'line.interlaced.read': 
                read_mode = 'lir'
            elif read_mode == 'fast-reset-read.read':
                read_mode = 'rrr-mpia'
            else:
                print("Read mode not recognized. Skipping file %s" % my_file)
                continue
            
            # ITIME
            itime = my_fits[0].header['ITIME']
            
            # Save mode
            if len(my_fits) > 1:
                save_mode = 'mef'
            else:
                save_mode = 'sef'
            
            # NCOADDS: if we have a cube, we take into account the number of layers
            if save_mode == 'mef' and len(my_fits[1].data.shape) > 2:
                ncoadds = my_fits[1].data.shape[0]
            elif save_mode == 'sef' and len(my_fits[0].data.shape) > 2:
                ncoadds = my_fits[0].data.shape[0]
            else:
                ncoadds = my_fits[0].header['NCOADDS']
                
            # Finaly, write out the parameters to text file
            if [read_mode, itime, ncoadds, save_mode] not in file_types:
                file_types.append([read_mode, itime, ncoadds, save_mode])
                # Insert into output file
                fd.write("%s\t%3.03f\t%d\t%s\n"%(read_mode, float(itime), int(ncoadds), save_mode))
        except Exception as e:
            print("Error while reading file: %s\n %s"%(my_file,str(e)))

    fd.close()
    shutil.move(output_file + "_tmp_", output_file)
    
    return len(file_types)
    

if __name__ == "__main__":

    
    usage = "usage: %prog [options] "
    desc = """Get the unique values of [read_mode, itime, ncoadd, save_mode]
for the files of a given directory to know the DARKs required for them."""

    parser = OptionParser(usage, description=desc)
    
                  
    parser.add_option("-s", "--source",
                  action="store", dest="source_file",
                  help="Source directory or txt file listing the filenames of input images to read.")
    
    parser.add_option("-o", "--output_file",
                  action="store", dest="output_file",default="filetypes.txt", 
                  help="Output file to be generated [default: %default]")
    
    parser.add_option("-r", "--recursive",
                  action="store_true", dest="recursive", default=False,
                  help="Recursive subdirectories if source is a directory name (only first level)")
    
    (options, args) = parser.parse_args()
    
    if len(sys.argv[1:])<1:
       parser.print_help()
       sys.exit(0)

    if not options.source_file or not options.output_file:
        parser.print_help()
        parser.error("incorrent number of arguments")
        
    n = getItimesNcoadds(options.source_file, options.output_file)
    print("%d types found" % n)
    
    sys.exit()
