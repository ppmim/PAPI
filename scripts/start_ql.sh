#!/bin/bash


export QT_GRAPHICSSYSTEM=native
QL_DIR=$PAPI_HOME/QL4

cd $QL_DIR
pyuic4 $QL_DIR/panicQL.ui > $QL_DIR/panicQL.py
pyrcc4 -o panicQL_resources_rc.py panicQL_resources.qrc

$QL_DIR/runQL.py  -c $QL_DIR/../config_files/papi.cfg
