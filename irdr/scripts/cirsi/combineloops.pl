#!/usr/bin/env perl
#
# run cubemean.c to combine loops for each run of each chip.  second pass 
# option uses sky-subtracted frames (extension .skysub) and scalar weights.
# standalone option read raw images in DATADIR
#
# Jan 2001
# May 2002 : standalone option added
#

die "Usage: combineloops.pl firstpass|secondpass|standalone [run_begin run_end]\n"
    unless ($#ARGV == 0 or $#ARGV == 2);

@chips   = (1, 2, 3, 4);
$loopbeg = 1;
$loopend = 20;
$listfn  = "./combineloops.list";
$scale   = "offset";
$method  = "mean";
$bitpix  = "short";
$|       = 1;

$basedir = $ENV{'IRDR_BASEDIR'};

die "Please setenv IRDR_BASEDIR /path/irdr\n"
    unless ($basedir);

$datadir = $ENV{'IRDR_DATADIR'};

die "Please setenv IRDR_DATADIR /path/irdr\n"
    unless ($datadir);

if ($#ARGV == 2) {
    @runs = ($ARGV[1] .. $ARGV[2]);
} else {
    die "Unable to find `runs' file\n" 
        unless (-e "runs");
    @runs = `cat runs`;
}

if ($ARGV[0] eq "firstpass") {               # first pass loop combine
    $ext = "";
    $weight = "noweight";
    $dir = "";
} elsif ($ARGV[0] eq "secondpass") {         # second pass
    $ext = ".skysub";
    $weight = "scalarweight";
    $dir = "";
} elsif ($ARGV[0] eq "standalone") {         # standalone loop combine
    $ext = "";
    $weight = "noweight";
    $scale = "none";    # no scale
    $dir = $datadir . "/";
    $bitpix = "float";  # float fits out
} else {
     die "combineloops.pl: incorrect usage\n";
}

$cubeprog = "$basedir/bin/cubemean";

die "Expected executable: $cubeprog\n"
    unless (-x "$cubeprog");

for ($i = 0; $i <= $#runs; $i++) {                # foreach run
    $run = sprintf("%05d", $runs[$i]);

    for ($j = 0; $j <= $#chips; $j++) {           # foreach chip
        $chip = "c$chips[$j]";
        $outfn = "irx_${run}_${chip}.fits" . "$ext";
        $ncombine = 0;

        open(COMBINE, "> $listfn");

        for ($k = $loopbeg; $k <= $loopend; $k++) {          # foreach loop
            $loop = sprintf("%03d", $k);
            $loopfn = "$dir" . "irx_${run}_${chip}_${loop}.fits" . "$ext";

            if (-e $loopfn) { 
                print COMBINE "$loopfn\n";
                $ncombine += 1; 
            }
        }

        close(COMBINE);
 
        if ($ncombine > 0) {
            print "\nCombining Run $run, Chip $chip, $ncombine loops\n";

            $cmd = "$cubeprog $listfn $outfn NULL $scale $method $weight $bitpix";

            print "Running combine command:\n$cmd\n";

            die "Sorry, command $cmd failed on list $listfn\n"
                unless (system($cmd) == 0);
        }
    }
}

if (-e $listfn) { `rm $listfn`; }

0;                 # Successful return
