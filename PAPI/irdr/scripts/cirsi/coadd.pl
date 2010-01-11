#!/usr/bin/env perl
#
# coadd a set of dithered observations using weight maps and fractional
# pixel image alignment.  call dithercubemean.c to calculate the weighted, 
# clipped, registered mean plane of the dither set.  dither offsets are read 
# from offsets.r{run_begin}to{run_end}
#
# Jan 2001
#

die "Usage: coadd.pl run_begin run_end\n" 
    unless ($#ARGV == 1);

@chips   = (1, 2, 3, 4);
$runbeg  = $ARGV[0];
$runend  = $ARGV[1];
$listfn  = "./coadd.list";
$infn    = "./offsets.r${runbeg}to${runend}";
$outfn   = "./irx.r${runbeg}to${runend}";

$basedir = $ENV{'IRDR_BASEDIR'};

die "Please setenv IRDR_BASEDIR /path/irdr\n"
    unless ($basedir);

$prog = "$basedir/bin/dithercubemean";

die "Expected executable $prog\n"
    unless (-x "$prog");

for ($i = 0; $i <= $#chips; $i++) {                  # foreach chip
    $chip = "c$chips[$i]";
    $gain = "gain.$chip.fits";
    $nlist = 0;

    open(OFFSET, "$infn") or
        die "Unable to open offsets file: $infn\n";

    open(FILELIST, ">$listfn") or
        die "Unable to open: $listfn\n";

    while (<OFFSET>) {                             # foreach dither frame
        chomp;

        ($runno, $xshift, $yshift, $match) = split;

        $fn = sprintf("./irx_%05d_$chip.fits.skysub", $runno);

        if (-e $fn && $match > 0.0) {
            print FILELIST "$fn $xshift $yshift\n";
            $nlist += 1;
        }
    }

    close(OFFSET);
    close(FILELIST);

    if ($nlist > 1) {                              # coadd dither set
        $cmd = "$prog $listfn $gain $outfn.$chip.fits $outfn.$chip.weight.fits";

        print "Running coadd command:\n$cmd\n";

        die "Sorry, command $cmd failed on list: $listfn\n"
            unless (system($cmd) == 0);
    }
}

if (-e $listfn) { `rm $listfn`; }
