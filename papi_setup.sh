#------------------------------------------------------------------------------
# User Configurable Settings
#------------------------------------------------------------------------------
# path to PAPI directory
export PAPI_HOME=${HOME}/DEVELOP/papi

# path to PAPI output data products
export PAPI_PROD=${HOME}/DataProd

#------------------------------------------------------------------------------
# Fixed Settings
#------------------------------------------------------------------------------
# path to PAPI reference files
export PAPI_CONFIG=${PAPI_HOME}/config_files/papi_panic2_o2k.cfg
export PATH=${PATH}:${PAPI_HOME}/bin
export PYTHONPATH=${PYTHONPATH}:${PAPI_HOME}


#--------------------------------
# IRAF settings
#--------------------------------
mkdir $HOME/.iraf
ln -s $HOME/.iraf $HOME/iraf
cd $HOME/.iraf
mkiraf 

# comment-out chkupdate on login.cl
## To be done

#-------------------------------
# IRDR build
#------------------------------
cd irdr
make all

