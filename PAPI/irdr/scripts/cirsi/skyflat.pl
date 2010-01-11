#!/usr/bin/env perl
#
# For each loop number and each chip, calculate the median image over all
# runs numbers.  Take the average of the median images to make a flatfield 
# per chip, then call gainmap.c to convert the flatfields to gain maps with 
# bad pixels automatically identified and set to 0. 
# Should use scale not offset when median combine?
#
# Jan 2001
# May 2002 : assume new controler in place, no loop1 reset problem

die "Usage: skyflat.pl [run_begin run_end]\n" 
    unless ($#ARGV == -1 or $#ARGV == 1);

@chips   = (1, 2, 3, 4);
@loops   = (1 .. 20);
$listfn  = "./skyflat.loops";
$listfn2 = "./skyflat.runs";
$|       = 1;                                  # turn off output buffering

if ($#ARGV == 1) {
    @runs = ($ARGV[0] .. $ARGV[1]);
} else {
    die "Unable to find `runs' file\n"
        unless (-e "runs");
    @runs = `cat runs`;
}

$basedir = $ENV{'IRDR_BASEDIR'};

die "Please setenv IRDR_BASEDIR /path/irdr\n"
    unless ($basedir);

$datadir = $ENV{'IRDR_DATADIR'};

die "Please setenv IRDR_DATADIR /path/irdr\n"
    unless ($datadir);

$cubeprog = "$basedir/bin/cubemean";

die "Expected executable $cubeprog\n"
    unless (-x "$cubeprog");

$bpmprog = "$basedir/bin/gainmap";

die "Expected executable $bpmprog\n"
    unless (-x "$bpmprog");

foreach $chip (@chips) {
    $flat = "./flat.c$chip.fits";
    $gain = "./gain.c$chip.fits";

    open(FLATLIST, "> $listfn2")
        or die "Unable to open $listfn2\n";

    foreach $loop (@loops) {
        $loopmed = sprintf("./flat.c%01d.%03d.fits", $chip, $loop);
        $ncombine = 0;

        open(COMBINE, "> $listfn")
            or die "Unable to open $listfn\n";

        foreach $run (@runs) {
            $infn = sprintf("irx_%05d_c%01d_%03d.fits", $run, $chip, $loop);
            $infn = "$datadir/$infn";

            if (-e $infn) { 
                print COMBINE "$infn\n";
                $ncombine += 1; 
            }
        }

        close(COMBINE);

        if ($ncombine > 0) {
            print "\nMedian combine to make: $loopmed\n";
            print FLATLIST "$loopmed\n";

            $cmd="$cubeprog $listfn $loopmed NULL offset median noweight short";

            print "Running combine command:\n$cmd\n";

            die "Sorry, command $cmd failed on: $listfn\n"
                unless (system($cmd) == 0);
        }
    }

    close(FLATLIST);

    $cmd = "$cubeprog $listfn2 $flat NULL offset median noweight short";

    print "Running combine command:\n$cmd\n";

    die "Sorry, command $cmd failed on: $listfn2\n"
        unless (system($cmd) == 0);

    print "\nScale $flat and set bad pixels to 0.0 to make: $gain\n";

    $cmd = "$bpmprog $flat $gain";

    print "Running gainmap command:\n$cmd\n";

    die "Sorry, command $cmd failed\n"
        unless (system($cmd) == 0);
}

if (-e $listfn) { `rm $listfn`; }
if (-e $listfn2) { `rm $listfn2`; }
