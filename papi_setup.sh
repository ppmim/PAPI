#------------------------------------------------------------------------------
# User Configurable Settings
#------------------------------------------------------------------------------
# path to PAPI source directory
export PAPI_HOME=${HOME}/DEVELOP/papi
export PAPI_BIN=${HOME}/bin

# path to PAPI output data products
export PAPI_PROD=${HOME}/DataProd

#------------------------------------------------------------------------------
# Settings
#------------------------------------------------------------------------------
# path to PAPI reference files
export PAPI_CONFIG=${PAPI_HOME}/config_files/papi_panic2_PANIC.cfg
export PATH=${PATH}:${PAPI_BIN}
export PYTHONPATH=${PYTHONPATH}:${PAPI_HOME}

#--------------------------------
# IRAF settings
#--------------------------------
mkdir $HOME/.iraf
ln -s $HOME/.iraf $HOME/iraf
cd $HOME/.iraf
mkiraf 

# comment-out chkupdate on login.cl
# and add the next:
# set stdimage        = imt4096
# set imtype          = "fits"
#if (access ("home$papi_ql_user.cl"))
#   if (access ("/tmp/focus_seq.txt"))
#      cl < "home$papi_ql_user.cl"
#;
 
## To be done

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
chmod $PAPI_BIN/start_iraf

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
ln -s $PAPI_HOME/reduce/refPixelCorrection.py $PAPI_BIN/refPixelCorrection
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


# ---------------------------------------
# Add Environment Variables to bash shell
# ---------------------------------------
cp ${PAPI_HOME}/scripts/bashrc ${HOME}/.bashrc

#echo "export PAPI_HOME=${HOME}/papi" >> ~/.bashrc
#echo "export PAPI_BIN=${HOME}/bin" >> ~/.bashrc 
#echo "export PAPI_PROD=${HOME}/DataProd" >> ~/.bashrc
#echo "export PAPI_CONFIG=${PAPI_HOME}/config_files/papi_panic2_PANIC.cfg" >> ~/.bashrc
#echo "export PYTHONPATH=${PYTHONPATH}:${PAPI_HOME}" >> ~/.bashrc
#echo "export PATH=\$PATH:${PAPI_BIN}" >> ~/.bashrc

source ~/.bashrc

echo "--------------------"
echo "PAPI was installed !"
echo "Run 'papi --test' to check all is succesful installed."
echo "------------------------------------------------------"
