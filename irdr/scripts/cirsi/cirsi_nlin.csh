#!/bin/csh -f
#
# create linearity coefficient images for CIRSI data reduction pipeline.
# edit the parameters below then:
# ./cirsi_nlin.csh >& cirsi_nlin.log &
#
# runs are supposed to be taken with a exptimeRef exposure between other
# exposures for calibration
#
# May 2002
#

#---- Edit these parameters

setenv IRDR_BASEDIR ~/IRDR                  # path to IRDR location
setenv IRDR_DATADIR ~/CIRSI/APR01/20010415  # path to linearity raw data

set runs = (3064 3103)        # run begin and run end of data set to process
set exptimeRef = 4            # exptime used as reference
set satlim = 40000            # don't use points over this high saturation limit
set makeruns = 1              # 0 to create "runs" file yourself

#---- Don't touch below here

set status = 0

echo ""
echo ">>>> Running CIRSI non-Linearity Mapping -- `date`"
echo ""

echo Basedir: $IRDR_BASEDIR
echo Datadir: $IRDR_DATADIR
echo Data runs: $runs
echo Make Runs: $makeruns
echo Reference ExpTime: $exptimeRef
echo Saturation Limit: $satlim

if ($makeruns) then
    $IRDR_BASEDIR/scripts/cirsi/runs.pl $runs

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

echo ""
echo ">>>> Combining loops -- `date`"
echo ""

$IRDR_BASEDIR/scripts/cirsi/combineloops.pl standalone

if ($status) then
    echo "** combineloops.pl failed, aborting"
    exit
endif

echo ""
echo ">>>> Non-Linearity Mapping -- `date`"
echo ""

$IRDR_BASEDIR/scripts/cirsi/nlinmap.pl $satlim $exptimeRef

if ($status) then
    echo "** nlinmap.pl failed, aborting"
    exit
endif

echo End: `date`
