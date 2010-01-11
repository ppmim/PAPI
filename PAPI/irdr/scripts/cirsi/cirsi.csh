#!/bin/csh -f
#
# run the CIRSI data reduction pipeline.  edit the parameters below then:
# ./cirsi.csh >& cirsi.log &
#
# set on_runs = (0 0) to make/use skyflats, otherwise domeflats are used
# or gain.c*.fits are taken from FLAT_DIR if you set makeflat to 0
#
# the script uses dithers.pl to determine the begin and end run numbers
# of the dither sets but this can fail, eg, if header RA/DEC are wrong.
# in this case, create a "dithers" file yourself and set dithers = 0 below.
# similarly, you can set makeruns = 0 and create the "runs" file yourself,
# eg, if you wanted to not use a particular run from a dither set
#
# If you want clipping during coaddition of fewer than 5 image planes,
# change mean.o to meangain.o in Makefile OBJ_MATH definition before compiling
# irdr and setenv GAIN to the actual detector gain.
#
# to turn off object mask dilation, set expandmask = 0.  the normal value
# of expandmask = 0.5 will grow the object mask regions by 50%
#
# if you chose 2mass for the astrometry, retreive a catalogue corresponding to
# your field on http://irsa.ipac.caltech.edu/applications/Gator/,
# with the first fields columns being :
# ra  dec  j_m j_msig  h_m  h_msig  k_m  k_msig
#
# May 2002
#

#---- Edit these parameters

setenv IRDR_BASEDIR ~/IRDR                      # path to IRDR location
setenv IRDR_DATADIR ~/CIRSI/16AUG00_1           # path to raw data
set DARK_DIR = "~/CIRSI/AUG00/DARK"	        # path to dark.c*.fits
set FLAT_DIR = "~/CIRSI/AUG00/FLAT_16AUG00_H"	# path to gain.c*.fits
set NLIN_DIR = "~/CIRSI/NLIN"                   # path to nlincoeff.dat
setenv CAT_2MASS ~/2MASS/cat_0129003.tbl        # path to 2MASS catalogue
setenv GAIN 2.0                                 # used in clipping bad pts

set runs = (18 53)            # run begin and run end of data set to process
set on_runs = (0 0)           # run begin and end of domeflat lamp on frames
set off_runs = (0 0)          # run begin and end of domeflat lamp off frames
set astrocat = "2mass"        # catalogue used for astrometry (usno,apm,2mass)
set expandmask = 0.5          # amount to expand the object mask regions

#---- Edit these parameters if you want

set darkcorr = 1           # 1 to substract dark using DARK_DIR files
set lincorr = 1            # 1 to correct for non-linearity using NLIN_DIR file
set makeflat = 0           # 0 to use FLAT_DIR files
set normalize = 0          # 0 to turn off final image normalization
set makedithers = 1        # 0 to create "dithers" file yourself
set makeruns = 1           # 0 to create "runs" file yourself
set deshadow = 1           # 0 to turn off deshadowing
set bpm_sigma = 5.0        # sigma threshold for identifying bad pixels
set bpm_nxblock = 16       # local bkg block size in x
set bpm_nyblock = 16       # local bkg block size in y
set bpm_mingain = 0.7      # min gain for good pixel
set bpm_maxgain = 1.3      # max gain for good pixel
set sky_halfwidth = 3      # half-width for running sky frame creation
set position_err = 15      # FITS hdr RA/DEC relative uncertainties [arcsec]
set offsets_minarea = 5    # SExtractor DETECT_MINAREA used in offsets step
set offsets_thresh = 4     # SExtractor DETECT_THRESH used in offsets step
set mask_minarea = 10      # DETECT_MINAREA used in object masking step
set mask_thresh = 2.7      # DETECT_THRESH used in object masking step

#---- Don't touch below here

if (-e $IRDR_BASEDIR/README) cp -f $IRDR_BASEDIR/README ./README.irdrversion

set status = 0

echo ""
echo ">>>> Running CIRSI Pipeline -- `date`"
echo ""

echo Basedir: $IRDR_BASEDIR
echo Datadir: $IRDR_DATADIR
echo Data runs: $runs
echo Dark Correction: $darkcorr
echo Non-linearity Correction: $lincorr
echo Make Flat: $makeflat
if ($makeflat) then
  echo Domeflat lamp on runs: $on_runs
  echo Domeflat lamp off runs: $off_runs
endif
echo Astrometry Catalogue: $astrocat
echo Make Dithers: $makedithers
echo Make Runs: $makeruns
echo Deshadow: $deshadow
echo Expand mask: $expandmask
echo BPM sigma: $bpm_sigma
echo BPM nxblock: $bpm_nxblock
echo BPM nyblcok: $bpm_nyblock
echo BPM mingain: $bpm_mingain
echo BPM maxgain: $bpm_maxgain
echo Sky Halfwidth: $sky_halfwidth
echo Position_Err: $position_err
echo Offsets Min Area: $offsets_minarea
echo Offsets Thresh: $offsets_thresh
echo Mask Min Area: $mask_minarea
echo Mask Thresh: $mask_thresh

echo ""
echo ">>>> Determining dither sets -- `date`"
echo ""

if ($makedithers) then
    $IRDR_BASEDIR/scripts/cirsi/dithers.pl $runs

    if ($status) then
        echo "** dithers.pl failed, aborting"
        exit
    endif
endif

if (! -e "dithers") then
    echo "Unable to find 'dithers' file"
    exit
endif

echo Dither sets are:
cat dithers

if ($status) then
    echo "Problem with 'dithers' file?"
    exit
endif

if ($makeruns) then
    $IRDR_BASEDIR/scripts/cirsi/runs.pl

    if ($status) then
        echo "** runs.pl failed, aborting"
        exit
    endif
endif

if (! -e "runs") then
    echo "Unable to find 'runs' file"
    exit
endif

echo Runs are:
cat runs

if ($status) then
    echo "Problem with 'runs' file?"
    exit
endif

if($makeflat) then
   echo ""
   echo ">>>> Producing flatfield and gain/bpm maps -- `date`"
   echo ""
   if ($on_runs[1] != 0) then
      $IRDR_BASEDIR/scripts/cirsi/domeflat.pl $on_runs $off_runs $bpm_sigma $bpm_nxblock $bpm_nyblock $bpm_mingain $bpm_maxgain
      if ($status) then
        echo "** domeflat.pl failed, aborting"
        exit
      endif
   else
      $IRDR_BASEDIR/scripts/cirsi/skyflat.pl
      if ($status) then
        echo "** skyflat.pl failed, aborting"
        exit
      endif
   endif
   $IRDR_BASEDIR/scripts/cirsi/normalize_gain.pl
   if ($status) then
     echo "** normalize_gain.pl failed, aborting"
     exit
   endif    
else
   echo ""
   echo ">> Copy flat.c*fits and gain.c*fits from $FLAT_DIR"
   echo ""
   cp $FLAT_DIR/flat.c*.fits .
   cp $FLAT_DIR/gain.c*.fits .
endif   

echo ""
echo ">>>> Applying flatfield to data -- `date`"
echo ""

if ($darkcorr) then
  cp $DARK_DIR/dark.c*.fits .
endif
if ($lincorr) then
  cp $NLIN_DIR/nlincoeff.dat .
endif

if ($darkcorr) then
  if ($lincorr) then
    $IRDR_BASEDIR/scripts/cirsi/lindarkflatten.pl
    if ($status) then
      echo "** lindarkflatten.pl failed, aborting"
      exit
    endif
  else
    $IRDR_BASEDIR/scripts/cirsi/darkflatten.pl
    if ($status) then
      echo "** darkflatten.pl failed, aborting"
      exit
    endif
  endif
else
  if ($lincorr) then
    $IRDR_BASEDIR/scripts/cirsi/linflatten.pl
    if ($status) then
      echo "** linflatten.pl failed, aborting"
      exit
    endif
  else
    $IRDR_BASEDIR/scripts/cirsi/flatten.pl
    if ($status) then
      echo "** flatten.pl failed, aborting"
      exit
    endif
  endif
endif  


echo ""
echo ">>>> Combining loops, firstpass -- `date`"
echo ""

$IRDR_BASEDIR/scripts/cirsi/combineloops.pl firstpass

if ($status) then
    echo "** combineloops.pl failed, aborting"
    exit
endif

echo ""
echo ">>>> Sky subtraction, firstpass -- `date`"
echo ""

$IRDR_BASEDIR/scripts/cirsi/skypass1.pl 2

if ($status) then
    echo "** skypass1.pl failed, aborting"
    exit
endif

echo ""
echo ">>>> Determining dither offsets -- `date`"
echo ""

$IRDR_BASEDIR/scripts/cirsi/runcmd.pl $IRDR_BASEDIR/scripts/cirsi/offsets.pl \
    $offsets_minarea $offsets_thresh $position_err

if ($status) then
    echo "** offsets.pl failed, aborting"
    exit
endif

echo ""
echo ">>>> Dither coaddition, first pass -- `date`"
echo ""

$IRDR_BASEDIR/scripts/cirsi/runcmd.pl $IRDR_BASEDIR/scripts/cirsi/coadd.pl

if ($status) then
    echo "** coadd.pl failed, aborting"
    exit
endif

echo ""
echo ">>>> Produce object masks -- `date`"
echo ""

$IRDR_BASEDIR/scripts/cirsi/runcmd.pl $IRDR_BASEDIR/scripts/cirsi/mask.pl \
    $expandmask $mask_minarea $mask_thresh

if ($status) then
    echo "** mask.pl failed, aborting"
    exit
endif

echo ""
echo ">>>> Sky subtraction, second pass -- `date`"
echo ""

$IRDR_BASEDIR/scripts/cirsi/skypass2.pl $sky_halfwidth

if ($status) then
    echo "** skypass2.pl failed, aborting"
    exit
endif

echo ""
echo ">>>> Combine loops, second pass -- `date`"
echo ""

$IRDR_BASEDIR/scripts/cirsi/combineloops.pl secondpass

if ($status) then
    echo "** combineloops.pl failed, aborting"
    exit
endif

if ($deshadow) then
    echo ""
    echo ">>>> Deshadow -- `date`"
    echo ""

    $IRDR_BASEDIR/scripts/cirsi/deshadow.pl

    if ($status) then
        echo "** deshadow.pl failed, aborting"
        exit
    endif
endif

echo ""
echo ">>>> Second pass coadd -- `date`"
echo ""

$IRDR_BASEDIR/scripts/cirsi/runcmd.pl $IRDR_BASEDIR/scripts/cirsi/coadd.pl

if ($status) then
    echo "** coadd.pl failed, aborting"
    exit
endif

echo ""
echo ">>>> Astrometry -- `date`"
echo ""

if !($astrocat == "") then
  mkdir noastro
  cp -f irx.r*c?.fits noastro
  $IRDR_BASEDIR/scripts/cirsi/astrometry.pl $astrocat irx.r*c?.fits
  $IRDR_BASEDIR/scripts/cirsi/avgwcs.pl irx.r*c?.fits
#  $IRDR_BASEDIR/scripts/cirsi/astrometry_refine.pl $astrocat irx.r*c?.fits
endif

if ($normalize) then
  mkdir nonorm
  cp -f irx.r*c?.fits nonorm
  $IRDR_BASEDIR/scripts/cirsi/normalize.pl irx.r*c?.fits
endif
  
echo End: `date`
