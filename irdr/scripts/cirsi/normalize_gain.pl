#!/usr/bin/env perl
#
# Use stats.c to calculate the signal level in the flatfield image for
# each chip (flat.c?.fits), then use normalize.c to scale the gainmap per 
# chip (gain.c?.fits) by the chip 1 flatfield signal level.
#
# Jan 2001
#

die "Usage: normalize_gain.pl\n" 
    unless ($#ARGV == -1);

$basedir = $ENV{'IRDR_BASEDIR'};

die "Please setenv IRDR_BASEDIR /path/irdr\n"
    unless ($basedir);

$statsprog = "$basedir/bin/stats";
$normprog = "$basedir/bin/normalize";
$headprog = "$basedir/extern/wcstools/bin/gethead";

die "Expected to find programs: $statsprog, $normprog, and $headprog\n"
    unless (-x $statsprog and -x $normprog and -x $headprog);

@signal = ();                      # signal level per chip in flatfield images
$nchips = 4;

for ($i = 1; $i <= $nchips; $i++) {
    $flat = "./flat.c$i.fits";

    die "Sorry, unable to find chip 1 flatfield: $flat\n"
        unless (-e $flat or $i != 1);                # require chip 1 flatfield

    if (-e $flat) {
        $cmd = "$statsprog updatehdr $flat > /dev/null";

        print "Running stats command:\n$cmd\n";

        die "Sorry, command failed: $cmd\n"
            unless (system($cmd) == 0);

        chomp($signal[$i-1] = `$headprog $flat DATAMODE`);

        print "$flat, BKG: $signal[$i-1]\n";
    }
}

print "Normalizing gain maps wrt Chip 1\n";

for ($i = 2; $i <= $nchips; $i++) {
    $norm = $signal[$i-1] / $signal[0];

    $gain = "./gain.c$i.fits";

    if (-e $gain) {
        $cmd = "$normprog $gain $norm scale";

        print "Running normalize command:\n$cmd\n";

        die "Sorry, command $cmd failed\n"
            unless (system($cmd) == 0);
    }
}
