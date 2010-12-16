#!/usr/bin/env perl
#
# create runs file
#
# May 2002
#

die "Usage: runs.pl [run_begin run_end]\n"
    unless ($#ARGV == -1 or $#ARGV == 1);

open(RUNS, "> runs") or
    die "Unable to open file: `runs'\n";

if ($#ARGV == 1) { 
    for $run ($ARGV[0] .. $ARGV[1]) { print RUNS "$run\n"; }
}
else {
    open(DITHERS, "dithers") or 
    die "Unable to open file: `dithers'\n";

    while (<DITHERS>) {
    	($runbeg, $runend) = split;
    	for $run ($runbeg .. $runend) { print RUNS "$run\n"; }
    }
    close(DITHERS);
}
close(RUNS);
