#!/usr/bin/perl -w
#
# Use apmcat.c to extract objects in sky region specified on command line,
# and use readcat.pl to convert object list to format expected by ESO SkyCat.
#
# Usage: apm.pl "RA DEC" BOX
# where BOX size is in arcmin, RA and DEC are hms dms J2000
#
# For example: apm.pl "10 00 00.0 -01 00 00.0" 4
#
# CNS, 27 Jan 2000
#

die "Usage: apm.pl RADEC BOX\n"
    unless ($#ARGV == 1);

$survey = "ukst";

$basedir = $ENV{'APM_BASEDIR'};

die "Please: setenv APM_BASEDIR /path/apmdir\n"
    unless ($basedir);

$tmpfile = "$basedir/apm.tmp";
$apmprog = "$basedir/apmcat";
$catprog = "$basedir/readcat.pl";

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
