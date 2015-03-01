#!/bin/bash

#/home/panic/SOFTWARE/PAPI/PyQt-x11-gpl-4.9.4/pyuic/pyuic4 panicQL.ui > panicQL.py

export QT_GRAPHICSSYSTEM=native
QL_DIR=$PAPI_HOME/QL4

cd $QL_DIR
pyuic4 $QL_DIR/panicQL.ui > $QL_DIR/panicQL.py

$QL_DIR/runQL.py  -c $QL_DIR/../config_files/papi.cfg
