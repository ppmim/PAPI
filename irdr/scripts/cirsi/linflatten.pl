#!/usr/bin/env perl
#
# run linflat.c to apply
# - linearity correction (using file nlincoeff.dat)
# - flatfield correction (using gain.c?.fits)
#
# May 2002
#

die "Usage: linflatten.pl [run_begin run_end]\n" 
    unless ($#ARGV == -1 or $#ARGV == 1);

@chips  = (1, 2, 3, 4);
@loops  = (1..99);
$listfn = "./linflatten.list";
$|      = 1;                             # turn off output buffering

if ($#ARGV == 1) {
    @runs = ($ARGV[0] .. $ARGV[1]);
} else {
    die "Unable to find `runs' file\n"
        unless (-e "runs");
    @runs = `cat runs`;
}                            

$datadir = $ENV{'IRDR_DATADIR'};

die "Please setenv IRDR_DATADIR /path/datadir\n"
    unless ($datadir);

$basedir = $ENV{'IRDR_BASEDIR'};

die "Please setenv IRDR_BASEDIR /path/irdr\n"
    unless ($basedir);

$linflatprog = "$basedir/bin/linflat";

die "Expected executable $linflatprog\n"
    unless (-x "$linflatprog");

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
    die "Unable to open file $flat\n"
        unless (-e $flat);

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

    print "\nApplying linearity coeff $a1 $a2 and flatfield $flat to data...\n";

    $cmd = "$linflatprog $listfn $a1 $a2 $flat";

    print "Running linflat command:\n$cmd\n";

    die "Sorry, command $cmd failed on: $listfn\n"
        unless (system($cmd) == 0);
}

if (-e $listfn) { `rm $listfn`; }
