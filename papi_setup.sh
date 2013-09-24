#------------------------------------------------------------------------------
# User Configurable Settings
#------------------------------------------------------------------------------
# path to PAPI directory
export PAPI_HOME=${HOME}/papi

# path to PAPI output data products
export PAPI_PROD=${HOME}/DataProd

#------------------------------------------------------------------------------
# Fixed Settings
#------------------------------------------------------------------------------
# path to PAPI reference files
export PAPI_CONFIG=${PAPI_HOME}/config_files
export PATH=${PATH}:${PAPI_HOME}/bin
export PYTHONPATH=${PYTHONPATH}:${PAPI_HOME}
