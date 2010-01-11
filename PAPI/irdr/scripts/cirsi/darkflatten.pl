#!/usr/bin/env perl
#
# run darkflat.c to do dark subtraction and flatfield correction to data.
# 
#
#
# 3-April-2009 : PANIC version (jmiguel@iaa.es) 
#
#
# Notes: 
#   - For the moment, no chip treatment will be done; I hope to run 4 parallel reductions
#
#
#


die "Usage: darkflatten.pl files.list dark.fits flat.fits\n" 
    unless ($#ARGV == 2);

$listfn = $ARGV[0];
$dark = $ARGV[1];
$flat = $ARGV[2];
  
$|      = 1;                             # turn off output buffering

                            

$basedir = $ENV{'IRDR_BASEDIR'};

die "Please setenv IRDR_BASEDIR /path/irdr\n"
    unless ($basedir);

$darkflatprog = "$basedir/bin/darkflat";

die "Expected executable $darkflatprog\n"
    unless (-x "$darkflatprog");

die "Unable to open dark image: $dark\n"
        unless (-e $dark);
          
die "Unable to open flat image: $flat\n"
        unless (-e $flat);
 

print "\nApplying dark $dark and flatfield $flat to data...\n";

$cmd = "$darkflatprog $listfn $dark $flat";

print "Running darkflat command:\n$cmd\n";

die "Sorry, command $cmd failed on: $listfn\n"
    unless (system($cmd) == 0);


