#!/usr/bin/env perl
#
# Create object masks (SExtractor OBJECTS images) for a list of FITS images.
# Expects the command "sex" (SExtractor Version 2+) in path.  If weight maps
# exist they will be used (assume weight map filename given by replacing .fits
# with .weight.fits).
#
# January 2001
#

die "Usage: makemask.pl detect_minarea detect_thresh *.fits\n" 
    unless ($#ARGV > 1);

$basedir = $ENV{'IRDR_BASEDIR'};

print $basedir

die "Please setenv IRDR_BASEDIR /path/irdr\n"
    unless ($basedir);

$area = shift(@ARGV);
$thresh = shift(@ARGV);
$config = $basedir . "/src/config/default.sex";

die "Sorry, unable to find SExtractor default.sex file: $config\n"
    unless (-e $config);
  
for $fn (@ARGV) {
    $wfn = $fn;
    $wfn =~ s|\.fits$|\.weight\.fits|;

    die "Unable to find file: $fn\n"
        unless (-e $fn);

    $wstr = (-e $wfn) ? " -WEIGHT_IMAGE $wfn -WEIGHT_TYPE MAP_WEIGHT " : "";

    $cmd = "sex $fn -c $config -FITS_UNSIGNED Y " . $wstr .
               " -DETECT_MINAREA $area -DETECT_THRESH $thresh -PIXEL_SCALE 0.45 -GAIN 4.15" .
               " -CHECKIMAGE_TYPE OBJECTS -CHECKIMAGE_NAME $fn.objs";

    die "Sorry, command failed: $cmd\n"              # run SExtractor
        unless (system($cmd) == 0);
}
