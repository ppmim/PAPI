#!/bin/bash

#/home/panic/SOFTWARE/PAPI/PyQt-x11-gpl-4.9.4/pyuic/pyuic4 panicQL.ui > panicQL.py

export QT_GRAPHICSSYSTEM=native

pyuic4 panicQL.ui > panicQL.py
pyrcc4 -o panicQL_resources_rc.py panicQL_resources.qrc

python runQL.py  -c ../config_files/papi.cfg
