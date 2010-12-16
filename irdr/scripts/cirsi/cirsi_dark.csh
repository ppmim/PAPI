#!/bin/csh -f
#
# run the creation of the CIRSI dark maps.
#
# edit the parameters below then:
# ./cirsi_dark.csh >& cirsi_dark.log &
#
# May 2002
#

#---- Edit these parameters

setenv IRDR_DATADIR ~/CIRSI/16AUG00   # path to raw data
setenv IRDR_BASEDIR ~/IRDR            # path to IRDR location

set dark_runs = (206 206)           # run begin and end of dark frames

#---- Don't touch below here

set status = 0

echo ""
echo ">>>> Running CIRSI Dark creation -- `date`"
echo ""

echo Basedir: $IRDR_BASEDIR
echo Datadir: $IRDR_DATADIR
echo Dark runs: $dark_runs

echo ""
echo ">>>> Producing dark maps -- `date`"
echo ""

$IRDR_BASEDIR/scripts/cirsi/dark.pl $dark_runs

if ($status) then
    echo "** dark.pl failed, aborting"
    exit
endif

echo End: `date`
