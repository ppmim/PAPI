#!/bin/bash

# After the complete installation, put this file in your home directory and make it executable (chmod u=rwx iraf). The command "./iraf" will then launch a complete IRAF session containing DS9, xgterm and ecl, based in ~/iraf

cd ~/iraf

ds9 &
sleep 2

xgterm -sb -title "Image Reduction and Analysis Facility (IRAF)" -fn 12x24 -bg grey -fg black -e "ecl" &
