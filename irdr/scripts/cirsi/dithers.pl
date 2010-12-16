#!/usr/bin/env perl
#
# run dithers.c to identify the dither sets within a data set.  use
# loop 1 raw data.
#
# Jan 2001
#

die "Usage: dithers.pl run_begin run_end\n"
    unless ($#ARGV == 1);

@chips = (1, 2, 3, 4);
@runs  = ($ARGV[0] .. $ARGV[1]);

$basedir = $ENV{'IRDR_BASEDIR'};

die "Please setenv IRDR_BASEDIR /path/irdr\n"
    unless ($basedir);

$datadir = $ENV{'IRDR_DATADIR'};

die "Please setenv IRDR_DATADIR /path/irdr\n"
    unless ($datadir);

$ditherprog = "$basedir/bin/dithers";

die "Expected executable: $ditherprog\n"
    unless (-x "$ditherprog");

@dithers = ();

for ($i = 0; $i <= $#runs; $i++) {                # foreach run
    $run = sprintf("%05d", $runs[$i]);

    for ($j = 0; $j <= $#chips; $j++) {           # foreach chip
        $chip = "c$chips[$j]";

        $fn = "$datadir/irx_${run}_${chip}_001.fits";

        if (-e $fn) { 
            push(@dithers, $fn);
        }
    }
}

die "No files found\n"
    unless ($#dithers >= 0);

$cmd = "$ditherprog @dithers";

print "Running command:\n$cmd\n";

die "Sorry, command $cmd failed\n"
    unless (system($cmd) == 0);
