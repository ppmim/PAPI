#!/usr/bin/env perl
#
# use SExtractor to produce a master object mask for each coadded dither 
# set.  the expandmask parameter is passed to dilate.c for doing object mask 
# dilation.  if expandmask is <= 0 then no dilation is done.  for the typical
# value of expandmask = 0.5, the mask regions are made 50% larger.
#
# Jan 2001
#

die "Usage: mask.pl run_begin run_end expandmask [detect_minarea " .
    "detect_thresh]\n" unless ($#ARGV == 2 or $#ARGV == 4);

@chips   = (1, 2, 3, 4);
$infn    = "./irx.r$ARGV[0]to$ARGV[1]";
$area    = ($#ARGV == 4) ? $ARGV[3] : 12;
$thresh  = ($#ARGV == 4) ? $ARGV[4] : 2.5;

$basedir = $ENV{'IRDR_BASEDIR'};

die "Please setenv IRDR_BASEDIR /path/irdr\n"
    unless ($basedir);

$dilateprog = "$basedir/bin/dilate";

die "Expected executable: $dilateprog\n"
    unless (-x $dilateprog);

for ($i = 0; $i <= $#chips; $i++) {                  # foreach chip
    $fn  = "$infn.c$chips[$i].fits";
    $wfn = "$infn.c$chips[$i].weight.fits";

    if (-e $fn) {
        die "Sorry, unable to find weight map: $wfn\n"
            unless (-e $wfn);

        $cmd = "sex $fn " .
               " -c $basedir/src/config/default.sex " .
               " -FITS_UNSIGNED Y " .
               " -DETECT_MINAREA $area " .
               " -DETECT_THRESH $thresh " .
               " -CATALOG_NAME $fn.cat " .
               " -WEIGHT_IMAGE $wfn " .
               " -WEIGHT_TYPE MAP_WEIGHT " .             # run SExtractor
               " -CHECKIMAGE_TYPE OBJECTS " .
               " -CHECKIMAGE_NAME $fn.objs ";

        print "Running SExtractor command:\n$cmd\n";

        die "Sorry, SExtractor failed on file: $fn\n"
            unless(system($cmd) == 0);

        $cmd = "$dilateprog $fn.objs $ARGV[2]";          # dilate object mask

        print "Dilating: $fn.objs\n";
        print "Running dilate command:\n$cmd\n";

        die "Sorry, program $dilateprog failed\n"
            unless (system($cmd) == 0);
    }
}
