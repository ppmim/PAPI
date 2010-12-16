#!/usr/bin/env perl
#
# use SExtractor to produce a background map for a coadded dither set, then
# subtract that bkgmap to defringe the image.
#
# Jan 2001
#

die "Usage: defringe.pl irx.r*c?.fits\n"
    unless ($#ARGV >= 0);

$basedir = $ENV{'IRDR_BASEDIR'};

die "Please setenv IRDR_BASEDIR /path/irdr\n"
    unless ($basedir);

$prog = "$basedir/bin/defringe";

die "Expected executable: $prog\n"
    unless (-x $prog);

while ($fn = shift(@ARGV)) {
    $wfn = $fn;
    $wfn =~ s/\.fits/\.weight\.fits/;

    die "Sorry, expected to find $wfn and $fn\n"
        unless (-e $wfn and -e $fn);

    $cmd = "sex $fn -c $basedir/src/config/default.sex -WEIGHT_IMAGE $wfn " .
           " -WEIGHT_TYPE MAP_WEIGHT -CHECKIMAGE_TYPE BACKGROUND " .
           " -CHECKIMAGE_NAME $fn.bkg -FITS_UNSIGNED Y";

    die "Sorry, SExtractor failed on file: $fn\n"
        unless(system($cmd) == 0);

    $cmd = "$prog $fn $fn.bkg";

    print "Running defringing command:\n$cmd\n";

    die "Sorry, command $cmd failed\n"
        unless (system($cmd) == 0);
}
