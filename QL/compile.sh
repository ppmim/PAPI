#!/bin/bash
#cp /disk-a/caha/cafos/TMP/panicQL.* .
pyuic panicQL.ui > panicQL.py
python runGUI.py
#echo "una prueba"
