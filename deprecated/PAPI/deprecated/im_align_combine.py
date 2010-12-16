#!/usr/bin/env python

################################################################################
#
#
# im_align_combine.py
#
# Last update 26/05/2009
#
################################################################################

import os
import glob
import shutil

# Pyraf modules
from pyraf import iraf
from iraf import noao
from iraf import imred
from iraf import ccdred

# Interact with FITS files
import pyfits

# Look for reference image
images_list='im.list'
coords_list='images.coords'
shifts_list='off.list'
ref_image=''
output_filename='im_align_combined.fits'


infile = open(images_list,"r")
ref_image = os.getcwd()+"/"+infile.readline()[:-1]
infile.close()

# Generation of images.coords file
print "REF_IMAGE", ref_image
os.system("/disk-a/caha/panic/DEVELOP/PIPELINE/PANIC/PAPI/irdr/makemask.pl 15 7.5 %s" %ref_image)
cat_file = open("/disk-a/caha/panic/DEVELOP/PIPELINE/PANIC/PAPI/test.cat","r")
coord_file = open(coords_list, "w")
n=0 # counter of number of objects included in coord_file
for line in cat_file.readlines():
    m_line=line.split()
    print "M_LINE=%s        %s" %(m_line[1], m_line[2])
    #if m_line[1]>1:
    if float(m_line[1])>200.0 and float(m_line[1])<1900.0 and float(m_line[2])>200.0 and float(m_line[2])<1900.0:
        coord_file.write( "%s   %s\n" %(m_line[1], m_line[2]) )
        print("INSERTED %s   %s\n" %(m_line[1], m_line[2]) )
        if n>5:
            break
        else:
            n=n+1
    else:
        print "Out of limits"

cat_file.close()
coord_file.close()
          
if not os.path.exists( ref_image ):
    print('No reference image "%s" found.' %ref_image)  
else:
    #Clean for old files
    for filename in glob.glob(os.getcwd()+"/"+output_filename) : os.remove(filename)
    
    iraf.imalign(input='@'+images_list,
                output='temp//@'+images_list,
                referenc=ref_image,
                coords= coords_list,
                shifts= shifts_list,
                boxsize= 7,
                bigbox=11,
                negative='no',
                backgro='INDEF',
                lower='INDEF',
                upper='INDEF',
                niterat=3,
                toleran=0,
                interp= 'poly3',
                boundar='constant',
                constan=0,
                trimima='no',
                verbose='yes')
             

    iraf.imcombine(input='temp//@'+images_list,
                output=output_filename,
                combine='aver',
                offset='none',
                reject='sigclip',
                lsigma=2.5,
                hsigma=2.5,
                scale='none',
                zero='median',
                masktype='goodv',
                maskvalue=1
                #masktype='none'
                #scale='exposure',
                #expname='EXPTIME'
                #ParList = _getparlistname ('flatcombine')
                )
                
    # Clean temp files
    #infile = open(images_list,"r")
    #for line in infile.readlines():
    #    os.remove(os.getcwd()+"/temp"+line.split('\n')[0])
    #infile.close()
            
                