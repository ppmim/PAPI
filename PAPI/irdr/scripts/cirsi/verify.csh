
#!/bin/csh
#
# Produce SExtractor catalogs and update FITS headers with data statistics
# (image mode, sigma, etc) for all sky-subtracted data frames and coadded
# dither sets
#
# Feb 2001
#

setenv IRDR_BASEDIR /home/sabbey/irdr

set area = 10
set thresh = 3
set files = (*.skysub irx.r*c?.fits)

if (-e gifs) rm -rf gifs

mkdir gifs

$IRDR_BASEDIR/bin/stats updatehdr $files
$IRDR_BASEDIR/scripts/cirsi/sextractor.pl $area $thresh $files

exit 0;
