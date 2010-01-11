#!/usr/bin/env perl
#
# Find the coordinate offset between each frame and the first frame in the
# list.  Use SExtractor to produce images with non-object pixels set to 0, 
# then run offsets.c to find the offsets by cross-correlation.  The offsets 
# are printed to stdout.
#
# Jan 2001
#

die "Usage: offsets.pl run_begin run_end [detect_minarea detect_thresh] " . 
    "[pos_err]\n" unless ($#ARGV > 0);

@chips  = (1, 2, 3, 4);
@runs   = ($ARGV[0] .. $ARGV[1]);
$listfn = "./offsets.list";
$outfn  = "./offsets.r$ARGV[0]to$ARGV[1]";
$poserr = ($#ARGV == 2 or $#ARGV == 4) ? ($poserr = $ARGV[$#ARGV]) : 10;
$area   = ($#ARGV >= 3) ? $ARGV[2] : 15;
$thresh = ($#ARGV >= 3) ? $ARGV[3] : 3.0;
$|      = 1;                                      # turn off output buffering

$basedir = $ENV{'IRDR_BASEDIR'};

die "Please setenv IRDR_BASEDIR /path/irdr\n"
    unless ($basedir);

$offsetprog = "$basedir/bin/offsets";

die "Expected to find program $offsetprog\n"
    unless (-e "$offsetprog");

$mergeprog = "$basedir/bin/avgoffsets";

die "Expected to find program $mergeprog\n"
    unless (-e "$mergeprog");

print STDERR "\nDetection Min Area: $area pixels\n";
print STDERR "Threshold: $thresh sigma\n";
print STDERR "RA/DEC keyword uncertainty: $poserr arcsec\n";

for ($i = 0; $i <= $#chips; $i++) {                  # foreach chip
    $nlist = 0;

    open(OUTLIST, ">$listfn")
        or die "Unable to open output file: $listfn\n";

    $chip = "c$chips[$i]";

    print STDERR "\nRunning SExtractor to create OBJECTS frames...\n";

    for ($j = 0; $j <= $#runs; $j++) {                  # foreach run
        $run = sprintf("%05d", $runs[$j]);
        $fn = "./irx_${run}_${chip}.fits.skysub";

        if (-e $fn) {
            $cmd = "sex $fn " .
                   " -c $basedir/src/config/default.sex " .
                   " -FITS_UNSIGNED Y " .
                   " -DETECT_MINAREA $area " .
                   " -DETECT_THRESH $thresh " .
                   " -CATALOG_NAME $fn.cat " .
                   " -CHECKIMAGE_TYPE OBJECTS " .
                   " -CHECKIMAGE_NAME $fn.objs ";

            print "Running SExtractor command:\n$cmd\n";

            die "Sorry, command $cmd failed\n"          # run SExtractor
                unless (system($cmd) == 0);

            print OUTLIST "$fn.objs\n";

            $nlist++;
        }
    }

    close(OUTLIST);

    if ($nlist > 1) {                       # call offset calculation program
        $cmd = "$offsetprog $listfn $poserr | tee $outfn.${chip}";

        print "Running offsets command:\n$cmd\n";

        die "Sorry, offsets failed on list: $listfn\n"
            unless (system($cmd) == 0);
    }
}

$cmd = "$mergeprog $ARGV[0] $ARGV[1]";

print "Running avg offsets command:\n$cmd\n";

die "Sorry, command $cmd failed on runs $ARGV[0] to $ARGV[1]\n"
    unless (system($cmd) == 0);
