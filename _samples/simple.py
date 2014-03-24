#!/usr/bin/env python

################################################################################
#
# symple test programm for PyRAF
#
# symple.py
#
# Last update 14/01/2008
#
################################################################################

#"""
#   Test routines for data reduction with PyRAF.
#"""

################################################################################

# Import necessary modules


import os
import sys
import shutil
import time

# Interact with FITS files
import pyfits


#Own modules
import mkBadPix


from pyraf import iraf
from iraf import noao
from iraf import imred
from iraf import ccdred


import threading
from threading import Thread


dark_framelist=['A0408060036.fits', 'A0408060037.fits', 'A0408160001.fits','A0408160002.fits','A0408160003.fits','A0408160004.fits','A0408160005.fits',
                'A0408160006.fits','A0408160007.fits','A0408160008.fits','A0408160009.fits','A0408160010.fits','A0408160011.fits','A0408160012.fits',
                'A0408160013.fits','A0408160014.fits','A0408160015.fits','A0408160016.fits','A0408160017.fits','A0408160018.fits','A0408160019.fits',
               'A0408160020.fits']

#dark_framelist=['A0408060036.fits', 'A0408060037.fits', 'A0408160001.fits']

flatK_begin_framelist=['A0408160021.fits','A0408160022.fits','A0408160023.fits','A0408160024.fits','A0408160025.fits','A0408160026.fits',
                       'A0408160327.fits','A0408160328.fits','A0408160329.fits','A0408160330.fits','A0408160331.fits']

flatK_begin_framelist_D=['A0408160021_D.fits','A0408160022_D.fits','A0408160023_D.fits','A0408160024_D.fits','A0408160025_D.fits','A0408160026_D.fits',
                         'A0408160327_D.fits','A0408160328_D.fits','A0408160329_D.fits','A0408160330_D.fits','A0408160331_D.fits']

science_framelist=['A0408160047.fits','A0408160048.fits','A0408160049.fits','A0408160050.fits','A0408160051.fits','A0408160052.fits','A0408160053.fits',
                   'A0408160054.fits','A0408160055.fits','A0408160056.fits','A0408160057.fits','A0408160058.fits','A0408160059.fits','A0408160060.fits']

science_framelist_D=['A0408160047_D.fits','A0408160048_D.fits','A0408160049_D.fits','A0408160051_D.fits','A0408160052_D.fits','A0408160053_D.fits',
                     'A0408160054_D.fits','A0408160055_D.fits','A0408160056_D.fits','A0408160057_D.fits','A0408160058_D.fits','A0408160059_D.fits','A0408160060_D.fits']

science_framelist_D_F=['A0408160047_D_F.fits','A0408160048_D_F.fits','A0408160049_D_F.fits','A0408160051_D_F.fits','A0408160052_D_F.fits','A0408160053_D_F.fits',
                     'A0408160054_D_F.fits','A0408160055_D_F.fits','A0408160056_D_F.fits','A0408160057_D_F.fits','A0408160058_D_F.fits','A0408160059_D_F.fits','A0408160060_D_F.fits']

data_dir='/disk-a/caha/panic/DATA/test2/'
out_dir=data_dir+'//out/'


dark_framelist_b=['/disk-a/caha/panic/TMP/data/A0408060036.fits', '/disk-a/caha/panic/TMP/data/A0408060037.fits']
flat_framelist=['/disk-a/caha/panic/TMP/data/A0408060036.fits', '/disk-a/caha/panic/TMP/data/A0408060037.fits']
sky_framelist=['/disk-a/caha/panic/TMP/data/A0408060036.fits', '/disk-a/caha/panic/TMP/data/A0408060037.fits']
sky_frame='/disk-a/caha/panic/TMP/data/A0408060036.fits'
new_pixmask='/disk-a/caha/panic/TMP/data/badpixmask'

global wk_dir
wk_dir='/disk-a/caha/panic/DATA/test1/'


################################################################################

def pipe( param ):

    try:
        flatframe = pyfits.open('/disk-a/caha/panic/TMP/data/A0408060036.fits','update')
        flatframe[0].header.add_history('Combined images by averaging (%i files) ' % 1)
        #flatframe.flush()
        #print flatframe[0].header
        #flatframe.writeto('/disk-a/caha/panic/TMP/data/A0408060036kk.fits')
        flatframe.close(output_verify='fix')
        #flatframe.close()
    except:
        pyfits.info(flatframe)
        print "Error !!!!"

    return

    checkIntegrity ( '/disk-a/caha/panic/TMP/data/A0408060037.fits')

    return 

    apply_dark ( science_framelist, 'out/masterDARK.fits')

    apply_flatfield (  science_framelist, 'out/masterFLAT_K_N.fits' )
    
    run_compute_sky( science_framelist_D_F )

    return
                 
    #if check_dirs()==-1 :
    #    print 'Error, pipe cannot continue.'
    #    return
    
    run_darkcombine ( dark_framelist, data_dir  )

    apply_dark ( flatK_begin_framelist, 'out/masterDARK.fits')

    run_flatcombine ( flatK_begin_framelist_D, data_dir )
    
    nomalize_flatfield ( out_dir+'masterFLAT_K.fits', out_dir)
    
    return 

    #Parallel Python
    #job_server = pp.Server() 
    #dark_framelistA=['A0408060036.fits', 'A0408060037.fits', 'A0408160001.fits']
    #f1 = job_server.submit(run_darkcombine, (dark_framelistA, data_dir,), (), ("os","pyraf","pyraf.iraf","time","sys",))
    #r1 = f1()

    mkBadPix.compute_badPixMask ( flatK_begin_framelist, 19000 )

    print '****** PIPE initiation ******'

    
    
    
    #run_compute_sky ( science_framelist_1 )
        
    print '->Computing master DARK frame ........' 
    ## Compute DARK of a frame list
    run_darkcombine( dark_framelist, data_dir )

    print '->Computing master FLAT frame (K filter) .......'
    ## Compute FLAT of a frame list 
    run_flatcombine ( flatK_begin_framelist, data_dir )

    print '->Normalizing master FLAT frame (K filter) .......'
    nomalize_flatfield ( out_dir+'masterFLAT_K.fits', out_dir)

    return

    print '->Applying dark frame......'
    apply_dark ( science_framelist, 'masterDARK.fits')

    print '->Applying flat-field frame .....'
    apply_flatfield ( science_framelist, 'masterFLAT_K.fits')

    
    return

    ## Compute SKY of a frame list of science frames
    run_compute_sky ( sky_framelist )
    
    subtract_sky ( framelist, 'masterDARK.fits' )
    
    
    # Apply flatfield to a framelist
    
    apply_flatfield ( dark_framelist , '/disk-a/caha/panic/TMP/data/A0408060036.fits' )
    
    create_pixmask ( sky_frame , new_pixmask)

    apply_pixmask ( sky_framelist , new_pixmask )

    ## Shift and Align the frame set
    # TODO !!!!

    ## Add all the frames 
    add_reduced_frameset ( sky_framelist,'reduced_frame.fits')
    
    print 'PIPE finished'
    
################################################################################
    
def run_darkcombine(framelist, dir):

    print 'Step 1 of Darkcombine'
    start_time = time.clock()
    
    st_frames=''

    for nframe in framelist:
        st_frames+=nframe+ ' , '

    print st_frames
        
    # Change to the source directory
    iraf.chdir(dir)

    print 'Step 2'
    
    # Call the noao.imred.ccdred task through PyRAF
    iraf.darkcombine(input=st_frames,
                     output=out_dir+'masterDARK.fits',
                     combine='median',
                     ccdtype='none',
                     process='no',
                     reject='sigclip',
                     #nlow='3',
                     #nhigh='3',
                     scale='none',
                     #expname='EXPTIME'
                     #ParList = _getparlistname('darkcombine')
                     )
    print 'Step 3'
    
    # Change back to the original working directory
    iraf.chdir()
    
    print "Time elapsed: ", time.time() - start_time, "s"

    
################################################################################

def run_flatcombine ( framelist , dir ):

    print 'Initiation of Flatcomine'
    start_time = time.time()
    st_frames=''

    for nframe in framelist:
        st_frames+=nframe+' , '

    # Change to the source directory
    iraf.chdir(dir)
    
    print 'Combining images : ' + st_frames
    iraf.flatcombine(input=st_frames,
                     output=out_dir+'masterFLAT_K.fits',
                     combine='median',
                     ccdtype='none',
                     process='no',
                     reject='sigclip',
                     scale='none',
                     #scale='exposure',
                     #expname='EXPTIME'
                     #ParList = _getparlistname ('flatcombine')
                     )
    # Change back to the original working directory
    iraf.chdir ()
    print "Time elapsed: ", time.time() - start_time, "s"
    
################################################################################

def nomalize_flatfield ( flat_frame, dir ):

    mean = 0.0
    normalized_flat_frame='masterFLAT_K_N.fits'
    
    print ' ******  Normalizing Flat-Field *********'

    #Change to the source directory
    iraf.chdir(dir)

    #Compute the mean of the image
    mean=float(iraf.imstat (
        images=flat_frame,
        fields='mean',Stdout=1)[1])

    print 'The mean value of image is: %(mean)f' %vars() 
    
    if os.path.exists(data_dir+normalized_flat_frame):
        print 'A Normalized flatfield already exist: '+data_dir+normalized_flat_frame+'...it will be removed..'
        os.remove(data_dir+normalized_flat_frame)
        
    iraf.imarith(operand1=flat_frame,
                 operand2=mean,
                 op='/',
                 result=normalized_flat_frame,
                 )

    print 'Flat normalized ---->' + normalized_flat_frame
    
    # Change back to the original working directory
    iraf.chdir ()

################################################################################
def apply_dark ( framelist, darkframe ):

    print ' ******** Applying dark_frame to images **********'

    st_frames=''
    st_frames_res=''

    # Build the frame list for IRAF 
    for nframe in framelist:
        st_frames+=nframe+' , '
        st_frames_res+=nframe.replace(".fits","_D.fits")+' , '

    print 'Ficheros entrada: ' + st_frames   
    print 'Ficheros salida: ' + st_frames_res
    
    #Change to the source directory
    iraf.chdir(data_dir)
    
    iraf.imarith(operand1=st_frames,
                 operand2=darkframe,
                 op='-',
                 result=st_frames_res,
                 )
    
    # Change back to the original working directory
    iraf.chdir ()  

################################################################################  
def run_compute_sky ( framelist ):

     print 'Computing sky ... '
     st_frames=''

     for nframe in framelist:
         st_frames+=nframe+' , '

     # Change to the source directory
     iraf.chdir(data_dir)
    
     print 'Images to be combined : ' + st_frames
     iraf.imcombine(input=st_frames,
                    output=out_dir+'sky_K.fits',
                    combine='median',
                    reject='minmax',
                    nlow=2,
                    nhigh=2
                     #ccdtype='none',
                     #process='no',
                     #ParList = _getparlistname ('flatcombine')
                     )
     # Change back to the original working directory
     iraf.chdir ()
     
################################################################################
# DO NOT WORK !!!
#################
def substrac_sky ( framelist, sky_frame ):

    print 'Subtracting SKY'

    st_frames=''
    st_frames_res=''

    for nframe in framelist:
        st_frames+=nframe+' , '
        st_frames_res+=nframe.replace(".fits","_nsky.fits")+' , '

    print 'Ficheros entrada: ' + st_frames   
    print 'Ficheros salida: ' + st_frames_res 
        
    #Change to the source directory
    iraf.chdir('/disk-a/caha/panic/TMP/data')
    
    print 'Images to be subtracted : ' + st_frames
    iraf.imarith(operand1=st_frames,
                 operand2='/disk-a/caha/panic/TMP/data/sky_K.fits',
                 #operand2=st_frames,
                 op='-',
                 result=st_frames_res
                 )
    
    # Change back to the original working directory
    iraf.chdir ()

    

    
################################################################################

def apply_flatfield ( framelist , flatframe ):

    print ' ******** Applying flat field to images **********'

    st_frames=''
    st_frames_res=''

    # Build the frame list for IRAF 
    for nframe in framelist:
        st_frames+=nframe.replace(".fits","_D.fits") +' , '
        st_frames_res+=nframe.replace(".fits","_D_F.fits")+' , '

    print 'Ficheros entrada: ' + st_frames   
    print 'Ficheros salida: ' + st_frames_res
    
    #Change to the source directory
    iraf.chdir(data_dir)
    
    iraf.imarith(operand1=st_frames,
                 operand2=flatframe,
                 op='/',
                 result=st_frames_res,
                 )
    
    # Change back to the original working directory
    iraf.chdir ()


################################################################################
def create_pixmask ( frame , new_pixmask):

    print ' ******** Creating pixel mask from image **********'
    
    #Change to the source directory
    iraf.chdir('/disk-a/caha/panic/TMP/data')
    
    iraf.ccdmask(image=frame,
                 mask=new_pixmask,
                 )
    
    # Change back to the original working directory
    iraf.chdir ()

################################################################################
def apply_pixmask ( framelist , pixmask ):

    print ' ******** Applying pixel mask from image **********'
    
    #Change to the source directory
    iraf.chdir('/disk-a/caha/panic/TMP/data')

    st_frames=''

    # Build the frame list for IRAF 
    for nframe in framelist:
        st_frames+=nframe+' , '

    print 'Ficheros entrada: ' + st_frames   
    
    iraf.fixpix(images=st_frames,
                 masks=pixmask,
                 )
    
    # Change back to the original working directory
    iraf.chdir ()


################################################################################
def add_reduced_frameset ( framelist, result_frame  ):

    print ' ******** Adding the reduced frameset  **********'
    
    #Change to the source directory
    iraf.chdir('/disk-a/caha/panic/TMP/data')

    st_frames=''

    # Build the frame list for IRAF 
    for nframe in framelist:
        st_frames+=nframe+' , '

    print 'Ficheros entrada: ' + st_frames   

    #iraf.clobber='yes'
    iraf.imcombine(input=st_frames,
                   output=result_frame,
                   combine='median',
                   reject='avsigclip',
                   scale='none',
                   zero='mode',
                   statsec='[500:600,500:600]',
                   )
    
    # Change back to the original working directory
    iraf.chdir ()    

################################################################################
def check_dirs():

    print '********* Previus check before running pipe *********'

    if not os.path.exists(data_dir):
        print 'Data directory : '+data_dir +'not available'
        return -1
    
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    else:
        for f in os.listdir(out_dir):
            if f.endswith('.fits'):
                print 'File -->' + f
                os.remove(out_dir+f)
                

################################################################################
def checkIntegrity(image):

    print '********* Checking Image integrity  *********'

    mean = 0.0
    mode = 0.0

    #Change to the source directory
    iraf.chdir('/disk-a/caha/panic/TMP/data')
    
    #Compute the mean of the image
    out=iraf.imstat (
        images=image,
        fields='mean,mode',Stdout=1)[1]

    values = out.strip().split()
    mean   = float(values[0])
    mode   = float(values[1])

    if (mode > 10.0 ):
        shutil.mv (image, "/tmp/")
    
    print 'The MEAN value of image is: %(mean)f' %vars()
    print 'The MODE value of image is: %(mode)f' %vars()
    
    
    # Change back to the original working directory
    iraf.chdir ()
                        

#################################################################################        


class myThread(Thread):
    
    dark_framelist=['A0408060036.fits', 'A0408060037.fits', 'A0408160001.fits']
    data_dir='/disk-a/caha/panic/DATA/data_alh'

    global flatK_begin_framelist_D
    def __init__(self):
        threading.Thread.__init__(self)
        
    def run(self):
        
        #The task to do
        #run_darkcombine( dark_framelist, data_dir )
        #run_flatcombine ( flatK_begin_framelist_D, data_dir ) 
        #out=iraf.imstat (
        #    images="/disk-a/caha/panic/DATA/data_alh/A0408060036.fits",
        #    fields='mean,mode',Stdout=1)[1]
        
        print "mean,mode=",out
        
        print 'CHILD: run_darkcombine  finished OK'

def testf():
    dark_framelist=['A0408060036.fits', 'A0408060037.fits', 'A0408160001.fits']
    data_dir='/disk-a/caha/panic/DATA/test1/'
    
    
    print "Start testf"
    
    #The task to do
    #run_darkcombine( dark_framelist, data_dir )
    global flatK_begin_framelist_D
    run_flatcombine ( flatK_begin_framelist_D, data_dir ) 
    
    print 'run_darkcombine  finished OK'

    
################################################################################      
 
if __name__=="__main__":
    
    #pipe(sys.argv[1:])
    #return

    testf()

    #return

    """aThread1 = myThread()
    aThread1.start()
    
    aThread2 = myThread()
    aThread2.start()
    
    while aThread1.isAlive() or aThread1.isAlive():
        print 'Threads...parent\'s heartbeat...'
        sleep(1)
    """    
    #while aThread2.isAlive():
    #    print 'Thread2...parent\'s heartbeat...'
    #    sleep(2)
    
    
#    print 'PARENT: Thread finished, but did apsum?'
        
################################################################################

