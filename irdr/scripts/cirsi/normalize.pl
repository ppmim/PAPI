#!/usr/bin/env perl
#
# run statsmask.c on input FITS files to obtain background levels, then run
# normalize.c to offset each to the average background level.
#
# Jan 2001
#

die "Usage: normalize.pl irx.r*c?.fits\n" 
    unless ($#ARGV >= 0);

$basedir = $ENV{'IRDR_BASEDIR'};

die "Please setenv IRDR_BASEDIR /path/irdr\n"
    unless ($basedir);

$statsprog = "$basedir/bin/statsmask";
$normprog = "$basedir/bin/normalize";
$headprog = "$basedir/extern/wcstools/bin/gethead";

die "Expected executables: $statsprog, $normprog\n"
    unless (-x $statsprog and -x $normprog and -x $headprog);

$avgbkg = 0.0;

for ($i = 0; $i <= $#ARGV; $i++) {
    $wfn = $fn = $ARGV[$i];                    # image filename
    $wfn =~ s/\.fits/\.weight\.fits/;          # weight map filename
    $mfn = $fn . ".objs";                      # object mask filename
   
    $cmd = "$statsprog $fn $wfn $mfn";

    print "Running stats command:\n$cmd\n";

    die "Sorry, command failed: $cmd\n"
        unless (system($cmd) == 0);

    chomp($bkg = `$headprog $fn DATAMODE`);

    print "$fn, BKG: $bkg\n";

    $avgbkg += $bkg;
}

$avgbkg = $avgbkg / ($#ARGV + 1);

print "Normalizing images to background level: $avgbkg\n";

for ($i = 0; $i <= $#ARGV; $i++) {
    $cmd = "$normprog $ARGV[$i] $avgbkg offset";

    print "Running normalize command:\n$cmd\n";

    die "Sorry, command $cmd failed\n"
        unless (system($cmd) == 0);
}
