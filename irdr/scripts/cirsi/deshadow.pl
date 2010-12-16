#!/usr/bin/env perl
#
# run reset_correct.c on each sky-subtracted, loop-combined data frame.
# this step is done right before second pass dither set coaddition, and
# removes the shadows seen near extended objects due to reset correction.
# reads "dithers" file.
#
# Jan 2001
#

die "Usage: deshadow.pl\n"
    unless ($#ARGV == -1);

@chips = (1, 2, 3, 4);
$listfn = "./deshadow.list";

$basedir = $ENV{'IRDR_BASEDIR'};

die "Please setenv IRDR_BASEDIR /path/irdr\n"
    unless ($basedir);

die "Unable to read `dithers' file\n"
    unless (-r "dithers");

@dithers = `cat dithers`;

$prog = "$basedir/bin/deshadow";

die "Expected executable: $prog\n"
    unless (-x "$prog");

foreach $chip (@chips) {
    open(FILELIST, "> $listfn")
        or die "Unable to open: $listfn\n";

    $nlist = 0;

    foreach $ditherset (@dithers) {
        ($runbeg, $runend) = split(/\s+/, $ditherset);

        $offsetsfn = "./offsets.r${runbeg}to${runend}";
        $objsfn = "./irx.r${runbeg}to${runend}.c$chip.fits.objs";

        open(OFFSET, "$offsetsfn") or
            die "Unable to open offsets file: $offsetsfn\n";

        while (<OFFSET>) {
            chomp;
            ($run, $xshift, $yshift, $match) = split;
            $fn = sprintf("irx_%05d_c%01d.fits.skysub", $run, $chip);

            if (-e $fn) {
                print FILELIST "$fn $objsfn $xshift $yshift\n";
                $nlist += 1;
            }
        }

        close(OFFSET);
    }

    close(FILELIST);

    if ($nlist > 0) {
        $cmd = "$prog $listfn gain.c$chip.fits rowcol";

        print "Running deshadow command:\n$cmd\n";

        die "Sorry, command $cmd failed on list: $listfn\n"
            unless (system($cmd) == 0);
    }
}
