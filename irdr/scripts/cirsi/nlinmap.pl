#!/usr/bin/env perl
#
# run nlinmap.c to compute linearity coefficient images for each chip.  
#
# May 2002
#

die "Usage: nlinmap.pl satlim exptimeRef [run_begin run_end]\n"
    unless ($#ARGV == 1 or $#ARGV == 3);

@chips   = (1, 2, 3, 4);
$nlinfn  = "./nlincoeff.dat";
$listfn  = "./nlinmap.list";
$|       = 1;

$basedir = $ENV{'IRDR_BASEDIR'};

die "Please setenv IRDR_BASEDIR /path/irdr\n"
    unless ($basedir);

$datadir = $ENV{'IRDR_DATADIR'};

die "Please setenv IRDR_DATADIR /path/irdr\n"
    unless ($datadir);

if ($#ARGV == 3) {
    @runs = ($ARGV[2] .. $ARGV[3]);
} else {
    die "Unable to find `runs' file\n" 
        unless (-e "runs");
    @runs = `cat runs`;
}

$prog = "$basedir/bin/nlinmap";

die "Expected executable: $prog\n"
    unless (-x "$prog");

$satlim = $ARGV[0];
$exptimeRef = $ARGV[1];

for ($j = 0; $j <= $#chips; $j++) {           # foreach chip
     $chip = "c$chips[$j]";
     $ncombine = 0;

     open(COMBINE, "> $listfn");

     for ($i = 0; $i <= $#runs; $i++) {                # foreach run
         $run = sprintf("%05d", $runs[$i]);

         $infn = "irx_${run}_${chip}.fits";

            if (-e $infn) { 
                print COMBINE "$infn\n";
                $ncombine += 1; 
            }
     }
     close(COMBINE);
 
     die "Can't find files for non-linearity calibration\n"
         unless($ncombine>0);

     print "\nNon-Linearity mapping of Chip $chip\n";

     $feOfmfn = "./feOfm.${chip}.dat";
     $cmd = "$prog $listfn $satlim $exptimeRef $feOfmfn $nlinfn";

     print "Running combine command:\n$cmd\n";
     die "Sorry, command $cmd failed on list $listfn\n"
        unless (system($cmd) == 0);
}

print "Coefficients written in file $nlinfn\n";

if (-e $listfn) { `rm $listfn`; }

0;                 # Successful return
