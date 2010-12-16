#!/usr/bin/env perl
#
# use apmcat.c to extract objects in region specified on command line.
# eg: apm.pl "10 00 00.0 -01 00 00.0" 4 ukst
#
# field size is in arcmin, assumes J2000 coordinates
#

die "Usage: apm.pl radec box poss|ukst [wcstools]\n" 
    unless ($#ARGV == 2 or $#ARGV == 3);

$survey = ($ARGV[2] eq "poss") ? "poss1" : "ukst";

$tmpfile = "./apm.tmp";

$basedir = $ENV{'IRDR_BASEDIR'};

die "Please setenv IRDR_BASEDIR /path/irdr\n"
    unless ($basedir);

$apmprog = "$basedir/extern/apmcat/apmcat";

$catprog = ($#ARGV == 3) ? "$basedir/scripts/cirsi/apmwcstools.pl" :
                           "$basedir/scripts/cirsi/apmgaia.pl";

die "Expected executable: $apmprog\n"
    unless (-x $apmprog);

die "Expected executable: $catprog\n"
    unless (-x $catprog);

$radecstr = shift;

$radecstr =~ s/:/ /g;                # eg, 10:00:00 -> 10 00 00

$box = shift;

$cmd = "$apmprog \"$radecstr\" box=$box equinox=j2000 survey=$survey " .
       " list=$tmpfile";

die "Sorry, apmcat failed\n"
    unless (system($cmd) == 0);

print `$catprog < $tmpfile`;

if (-e $tmpfile) {
    `rm -f $tmpfile`;
}
