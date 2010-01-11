#!/bin/csh
#
# Collect statistics on FITS files in the current working directory and
# make sky brightness plot.
#
# Uses IRDR and ECLIPSE, which should be in your path, on APM2:
# set path = ($path /home/sabbey/irdr/bin /home/sabbey/src/eclipse/bin)
#
# imagetype can be gif or postscript or ...
#

set imagetype=gif

stats updatehdr irx*fits

dfits irx*fits | \
    fitsort OBJECT FILTER NRUN CHIP DATAMODE DATASIG >! stats.txt

dfits irx*c1*fits | \
    fitsort NRUN DATAMODE DATASIG | tail +2 | awk '{print $2,$3,$4}' >! \
        stats.c1.txt

dfits irx*c2*fits | \
    fitsort NRUN DATAMODE DATASIG | tail +2 | awk '{print $2,$3,$4}' >! \
        stats.c2.txt

dfits irx*c3*fits | \
    fitsort NRUN DATAMODE DATASIG | tail +2 | awk '{print $2,$3,$4}' >! \
        stats.c3.txt

dfits irx*c4*fits | \
    fitsort NRUN DATAMODE DATASIG | tail +2 | awk '{print $2,$3,$4}' >! \
        stats.c4.txt

gnuplot <<EOT
set xlabel "Run Number"
set ylabel "Background Level"
set title "Sky Brightness Verse Run Number"
set term $imagetype
set output "stats.gif"
plot "stats.c1.txt" title "Chip 1" with errorbars, \
     "stats.c2.txt" title "Chip 2" with errorbars, \
     "stats.c3.txt" title "Chip 3" with errorbars, \
     "stats.c4.txt" title "Chip 4" with errorbars
EOT

exit 0;
