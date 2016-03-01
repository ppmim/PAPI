#!/bin/bash

#------------------------------------------------------------------------------
# User Configurable Settings
#------------------------------------------------------------------------------
# path to PAPI source directory
PAPI_HOME=${HOME}/papi
PAPI_BIN=${HOME}/bin
PAPI_CONFIG=${PAPI_HOME}/config_files/papi.cfg


#--------------------------------
# IRAF settings
#--------------------------------
mkdir $HOME/.iraf
ln -s $HOME/.iraf $HOME/iraf
cd $HOME/.iraf
mkiraf 
# Copy custom iraf scripts
cp ${PAPI_HOME}/scripts/login.cl ${HOME}/iraf
cp ${PAPI_HOME}/scripts/papi_ql_user.cl ${HOME}/iraf
user=$(whoami)
sed -i "s/panic/$user/g" ${HOME}/iraf/login.cl


#-------------------------------
# IRDR build
#------------------------------
cd $PAPI_HOME/irdr
make clean; make all


#-------------------------------
# Create symlinks to PAPI_BIN
#------------------------------
mkdir $PAPI_BIN

cp $PAPI_HOME/scripts/start_ql.sh $PAPI_BIN/start_ql
chmod a+x $PAPI_BIN/start_ql
cp $PAPI_HOME/scripts/start_iraf.sh $PAPI_BIN/start_iraf
chmod a+x $PAPI_BIN/start_iraf

# To check the PANIC temperatures and press
cp -av $PAPI_HOME/scripts/panic_status $PAPI_BIN/


chmod a+x $PAPI_HOME/papi.py
ln -s $PAPI_HOME/papi.py $PAPI_BIN/papi

chmod a+x $PAPI_HOME/reduce/*.py
ln -s $PAPI_HOME/reduce/applyDarkFlat.py $PAPI_BIN/applyDarkFlat
ln -s $PAPI_HOME/reduce/astrowarp.py $PAPI_BIN/astrowarp
ln -s $PAPI_HOME/reduce/calBPM.py $PAPI_BIN/calBPM
ln -s $PAPI_HOME/reduce/calCombineFF.py $PAPI_BIN/calCombineFF
ln -s $PAPI_HOME/reduce/calDark.py $PAPI_BIN/calDark
ln -s $PAPI_HOME/reduce/calDarkModel.py $PAPI_BIN/calDarkModel
ln -s $PAPI_HOME/reduce/calDomeFlat.py $PAPI_BIN/calDomeFlat
ln -s $PAPI_HOME/reduce/calGainMap.py $PAPI_BIN/calGainMap
ln -s $PAPI_HOME/reduce/calSuperFlat.py $PAPI_BIN/calSuperFlat
ln -s $PAPI_HOME/reduce/checkQuality.py $PAPI_BIN/checkQuality
ln -s $PAPI_HOME/reduce/correctNonLinearity.py $PAPI_BIN/correctNonLinearity
ln -s $PAPI_HOME/reduce/dxtalk.py $PAPI_BIN/dxtalk
ln -s $PAPI_HOME/reduce/eval_focus_serie.py $PAPI_BIN/eval_focus
ln -s $PAPI_HOME/reduce/makeobjmask.py $PAPI_BIN/makaobjmask
ln -s $PAPI_HOME/reduce/remove_cosmics.py $PAPI_BIN/remove_cosmics
ln -s $PAPI_HOME/reduce/solveAstrometry.py $PAPI_BIN/solveAstrometry

chmod a+x $PAPI_HOME/misc/*.py
ln -s $PAPI_HOME/misc/check_papi_modules.py $PAPI_BIN/check_papi_modules
ln -s $PAPI_HOME/misc/collapse.py $PAPI_BIN/collapse
ln -s $PAPI_HOME/misc/genLogsheet.py $PAPI_BIN/genLogsheet
ln -s $PAPI_HOME/misc/health.py $PAPI_BIN/health
ln -s $PAPI_HOME/misc/mef.py $PAPI_BIN/mef
ln -s $PAPI_HOME/misc/modFITS.py $PAPI_BIN/modFITS

chmod a+x $PAPI_HOME/photo/photometry.py
ln -s $PAPI_HOME/photo/photometry.py $PAPI_BIN/photometry

# Some tools for commissioning
ln -s $PAPI_HOME/commissioning/runStarfocus.py $PAPI_BIN/runStarfocus.py
ln -s $PAPI_HOME/commissioning/p_50_tiltcheck.py $PAPI_BIN/p_50_tiltcheck.py
ln -s $PAPI_HOME/commissioning/getImageOffsets.py $PAPI_BIN/getImageOffsets.py
# getDarks: Used by the OT to genereate all the darks for a given night directory.
ln -s $PAPI_HOME/commissioning/getDarks.py $PAPI_BIN/getDarks.py
chmod a+x $PAPI_HOME/commissioning/getDarks.py
# tool to edit FITS headers
ln -s $PAPI_HOME/irdr/extern/wcstools/bin/edhead $PAPI_BIN/edhead
chmod a+x $PAPI_BIN/edhead

# ---------------------------------------
# Add Environment Variables to bash shell
# ---------------------------------------
cp ${PAPI_HOME}/scripts/bashrc.papi ${HOME}/.papirc

# Create a backup of current files
if [ -f ~/.bashrc ]; then
    cp -av ~/.bashrc ~/.bashrc.bak
else
    echo "***Error***: no ~/.bashrc found"
fi

if [ -f ~/.papirc ]; then
    cp -av ~/.papirc ~/.papirc.bak
fi

echo "if [ -f ~/.papirc ]; then . ~/.papirc; fi">> ~/.bashrc
source ~/.bashrc


echo "--------------------"
echo "PAPI was installed !"
echo "Run 'papi -e' to check all is succesful installed."
echo "--------------------------------------------------"
