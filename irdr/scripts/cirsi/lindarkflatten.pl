#!/usr/bin/env perl
#
# run lindarkflat.c to apply
# - linearity correction (using file nlincoeff.dat)
# - dark subtraction (using dark.c?.fits)
# - flatfield correction (using gain.c?.fits)
#
# 3-April-2009 : PANIC version (jmiguel@iaa.es) 
#
#
# Notes: 
#   - For the moment, no chip treatment will be done; I hope to run 4 parallel reductions
#
#
#

die "Usage: lindarkflatten.pl datafiles.list\n" 
    unless ($#ARGV == -1 or $#ARGV == 1);

#@chips  = (1, 2, 3, 4);
@chips = (1);
$listfn = $ARGV[0];
$|      = 1;                             # turn off output buffering


$basedir = $ENV{'IRDR_BASEDIR'};

die "Please setenv IRDR_BASEDIR /path/irdr\n"
    unless ($basedir);

$lindarkflatprog = "$basedir/bin/lindarkflat";

die "Expected executable $lindarkflatprog\n"
    unless (-x "$lindarkflatprog");

$lincoeff = "./nlincoeff.dat";
open(LINCOEFF,$lincoeff) or
    die "Unable to open file $lincoeff\n";

for ($i = 0; $i <= $#chips; $i++) {                  # foreach chip
    $chip = "c$chips[$i]";

    $_ = <LINCOEFF>; chomp;
    ($a1,$a2) = split;

    die "check $lincoeff file format\n"
        unless ( ($a1!=0)&&($a2!=0) );

    $flat = "./gain.${chip}.fits";
    die "Unable to open flatfield: $flat\n"
        unless (-e $flat);

    $dark = "./dark.${chip}.fits";
    die "Unable to open dark image: $dark\n"
        unless (-e $dark);

    open(LIST, "> $listfn");

    for ($j = 0; $j <= $#runs; $j++) {                # foreach run
        $run = sprintf("%05d", $runs[$j]);

        for ($l = 0; $l <= $#loops; $l++) {             # foreach loop
            $loopno = sprintf("%03d", $loops[$l]);
            $outfn = "irx_${run}_${chip}_${loopno}.fits";
            $infn  = "$datadir/$outfn";

            if (-e $infn) { 
                print LIST "$infn $outfn\n";
            }
        }
    }

    close(LIST);

    print "\nApplying linearity coeffs $a1 $a2, dark $dark and flatfield $flat to data...\n";

    $cmd = "$lindarkflatprog $listfn $a1 $a2 $dark $flat";

    print "Running lindarkflat command:\n$cmd\n";

    die "Sorry, command $cmd failed on: $listfn\n"
        unless (system($cmd) == 0);
}

if (-e $listfn) { `rm $listfn`; }
