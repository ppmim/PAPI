#!/usr/bin/env perl
#
# median combine for each chip to create dark.c?.fits
#
# 2-April-2009 : PANIC version (jmiguel@iaa.es) 
#
#
# Notes: 
#   - For the moment, no chip treatment will be done; I hope to run 4 parallel reductions
#
#
#
# usage: dark.pl dark_files.list
# 
# example: dark.pl my_dark_files.list
# 
# input: raw dark fits files
#
# output: master dark file per chip 'i'
#
# method: median combine the dark images
#
# notes:

die "Usage: dark.pl files.list\n" 
    unless ($#ARGV == 0);

#@chips   = (1, 2, 3, 4);
@chips   = (1);
$list  = ($ARGV[0]);
$|       = 1;                             # turn off output buffering

$basedir = $ENV{'IRDR_BASEDIR'};

die "Please setenv IRDR_BASEDIR /path/irdr\n"
    unless ($basedir);

$cubeprog = "$basedir/bin/cubemean";

die "Expected executable $cubeprog\n"
    unless (-x $cubeprog);

for ($i = 0; $i <= $#chips; $i++) {                   # foreach chip
    $chip = "c$chips[$i]";


    chomp($nr = `wc $list | awk '{print \$1}'`);
      
    die "Warning: Found no dark run, chip $chip\n"
        unless ($nr > 0);

    $darkfn = "./dark.${chip}.fits";
    print "\nMedian combine dark images to make dark.c${chip}.fits\n";

    $cmd = "$cubeprog $list $darkfn NULL none median noweight float";

    print "Running combine command:\n$cmd\n";

    die "Sorry, command $cmd failed on: $list\n"
         unless (system($cmd) == 0);

}

