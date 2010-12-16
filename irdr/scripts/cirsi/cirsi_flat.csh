#!/bin/csh -f
#
# run the creation of the CIRSI flatfield maps.
#
# edit the parameters below then:
# ./cirsi_flat.csh >& cirsi_flat.log &
#
# May 2002
#

#---- Edit these parameters

setenv IRDR_BASEDIR ~/IRDR      # path to IRDR location
setenv IRDR_DATADIR ~/CIRSI/16AUG00_Dome_H   # path to raw data

set on_runs = (206 206)        # run begin and end of domeflat lamp on frames
set off_runs = (207 207)       # run begin and end of domeflat lamp off frames

#---- Edit these parameters if you want

set bpm_sigma = 5.0           # sigma threshold for identifying bad pixels
set bpm_nxblock = 16          # local bkg block size in x
set bpm_nyblock = 16          # local bkg block size in y
set bpm_mingain = 0.7         # min gain for good pixel
set bpm_maxgain = 1.3         # max gain for good pixel

#---- Don't touch below here

set status = 0

echo ""
echo ">>>> Running CIRSI Flatfield creation -- `date`"
echo ""

echo Basedir: $IRDR_BASEDIR
echo Datadir: $IRDR_DATADIR
echo Domeflat lamp on runs: $on_runs
echo Domeflat lamp off runs: $off_runs
echo BPM sigma: $bpm_sigma
echo BPM nxblock: $bpm_nxblock
echo BPM nyblcok: $bpm_nyblock
echo BPM mingain: $bpm_mingain
echo BPM maxgain: $bpm_maxgain

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

echo End: `date`
