#!/bin/bash
#cp /disk-a/caha/cafos/TMP/panicQL.* .

pyuic panicQL.ui > panicQL.py

if [[ `hostname -s` = udit43* ]]; then
	python runQL.py  -c ../config_files/papi_suse11.cfg
else
	python runQL.py  -c ../config_files/papi_irws2.cfg
fi
