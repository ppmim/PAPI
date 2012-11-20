#!/bin/bash
#cp /disk-a/caha/cafos/TMP/panicQL.* .

/home/panic/SOFTWARE/PAPI/PyQt-x11-gpl-4.9.4/pyuic/pyuic4 panicQL.ui > panicQL.py

if [[ `hostname -s` = udit22* ]]; then
	python runQL.py  -c ../config_files/papi_suse11.cfg
else
	#python runQL.py  -c ../config_files/papi_irws2.cfg
	python runQL.py  -c ../config_files/papi_portatil.cfg
fi