#!/bin/bash

#/home/panic/SOFTWARE/PAPI/PyQt-x11-gpl-4.9.4/pyuic/pyuic4 panicQL.ui > panicQL.py

export QT_GRAPHICSSYSTEM=native
QL_DIR=/home/panic/DEVELOP/papi/QL4

cd $QL_DIR
pyuic4 $QL_DIR/panicQL.ui > $QL_DIR/panicQL.py

if [[ `hostname -s` = udit22* ]]; then
	python $QL_DIR/runQL.py  -c $QL_DIR/../config_files/papi.cfg
else
	python $QL_DIR/runQL.py  -c $QL_DIR/../config_files/papi_panic2_hawki.cfg
fi
