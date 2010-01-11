#!/usr/bin/env perl
#
# Usage: avgwcs.pl irx.r*c?.fits
#
# Use avgkey.c to average the WCS RA/DEC scale factors (CD1_1, CD2_2) over 
# all chips, then average the cross terms (CD1_2, CD2_1) over each chip
# separately.
#
# Jan 2001
#

die "Usage: avgwcs.pl irx.r*c?.fits\n" 
    unless ($#ARGV >= 0);

$nchips = 4;

$basedir = $ENV{'IRDR_BASEDIR'};

die "Please setenv IRDR_BASEDIR /path/irdr\n"
    unless ($basedir);

$avgprog = "$basedir/bin/avgkey";

die "Expected executable $avgprog\n"
    unless (-x "$avgprog");

$cmd = "$avgprog CD1_1 @ARGV";

die "Sorry, command failed: $cmd\n"
    unless (system($cmd) == 0);

$cmd = "$avgprog CD2_2 @ARGV";

die "Sorry, command failed: $cmd\n"
    unless (system($cmd) == 0);

for ($i = 1; $i <= $nchips; $i++) {
    $chip = "c${i}.fits";
    @args = ();

    for ($j = 0; $j <= $#ARGV; $j++) {
        if ($ARGV[$j] =~ /${chip}$/) {
            push(@args, $ARGV[$j]);
        }
    }

    if ($#args > 0) {
        $cmd = "$avgprog CD1_2 @args";

        die "Sorry, command failed: $cmd\n"
            unless (system($cmd) == 0);

        $cmd = "$avgprog CD2_1 @args";

        die "Sorry, command failed: $cmd\n"
            unless (system($cmd) == 0);
    }
}
