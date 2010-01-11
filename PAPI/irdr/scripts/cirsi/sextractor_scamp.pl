#!/usr/bin/env perl
#
# Run SExtractor on a list of FITS image (with optional weight images) to 
# produce catalog files.  If weight maps exist they will be used (assume 
# weight map filename given by replacing .fits with .weight.fits).
#
# Jan 2001
#

die "Usage: sextractor_scamp.pl detect_minarea detect_thresh *.fits\n" 
    unless ($#ARGV >= 2);

$basedir = $ENV{'IRDR_BASEDIR'};

die "Please setenv IRDR_BASEDIR /path/irdr\n"
    unless ($basedir);

$area   = shift(@ARGV);
$thresh = shift(@ARGV);
$config = $basedir . "/src/config/default_scamp.sex";

die "Sorry, unable to find SExtractor default.sex file: $config\n"
    unless (-e $config);

for $fn (@ARGV) {
    $wfn = $fn;
    $wfn =~ s|\.fits|\.weight\.fits|;

    die "Unable to find file: $fn\n"
        unless (-e $fn);

    $wstr = (-e $wfn) ? " -WEIGHT_IMAGE $wfn -WEIGHT_TYPE MAP_WEIGHT " : "";

    $cmd = "sex $fn -c $config " . $wstr .
              " -PARAMETERS_NAME $basedir/src/config/verify_scamp.param " .
              " -FILTER N " .
              " -FITS_UNSIGNED Y " .
              " -CATALOG_TYPE FITS_LDAC ".
              " -CATALOG_NAME $fn.ldac " .
              " -DETECT_MINAREA $area " .
              " -DETECT_THRESH $thresh " .
              " -SATUR_LEVEL 300000.0".
              " -CHECKIMAGE_TYPE NONE " ;

    die "Sorry, command failed: $cmd\n"              # run SExtractor
        unless (system($cmd) == 0);
}
