#!/usr/bin/env perl
#
# read beginning and end run numbers of dither sets from ASCII list
# "dithers" and run the specified command on each dither set.  assumes
# that the first two arguments of the specified command are runbegin
# and runend and just passes on the additional arguments.
#
# Eg: runcmd.pl offsets.pl 15 2.5
#
# where "dithers" contains:
# 5792 5800
# 5801 5809
# 5810 5818
# 5819 5827
#
# Jan 2001
#

die "Usage: runcmd.pl cmd arg1 arg2\n"
    unless ($#ARGV >= 0);

$| = 1;                             # turn off output buffering

$cmdtorun = shift(@ARGV);

open(DITHERS, "dithers") or
    die "Unable to open ASCII dither list: dithers\n";

while (<DITHERS>) {
    chomp;

    s/^\s*//;                 # remove leading spaces

    if (! /^\d/) {            # skip line unless starts with number
        next; 
    }

    ($runbegin, $runend) = split(/\s+/, $_);

    $cmd = "$cmdtorun $runbegin $runend @ARGV";

    print "Running command: $cmd\n";

    die "Sorry, command failed\n"
        unless (system($cmd) == 0);
}
