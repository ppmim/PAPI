#!/bin/bash
#cp /disk-a/caha/cafos/TMP/panicQL.* .
pyuic panicQL.ui > panicQL.py
python runQL.py  -c ../config_files/papi_suse11.cfg
