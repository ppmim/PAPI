#!/usr/bin/env perl
#
# call skyfilter.c to do 2nd pass sky subtraction.  reads "dithers" file.
# processes each loop separately.
#
# Jan 2001
#

die "Usage: skypass2.pl half-width\n" 
    unless ($#ARGV == 0);

@chips  = (1, 2, 3, 4);
@loops  = (1..9);
@ramp   = ("rowcol", "colrow", "rowcol", "colrow");
$hwidth = $ARGV[0];
$listfn = "./skypass2.list";

$basedir = $ENV{'IRDR_BASEDIR'};

die "Please setenv IRDR_BASEDIR /path/irdr\n"
    unless ($basedir);

$skyprog = "$basedir/bin/skyfilter";

die "Expected executable $skyprog\n"
    unless (-x "$skyprog");

die "Unable to read `dithers' file\n" 
    unless (-r "dithers");

@dithers = `cat dithers`;

print STDERR "Second pass sky subtraction and reset correction\n\n";

foreach $loop (@loops) {                      # foreach loop
    foreach $chip (@chips) {                        # foreach chip
        open(FILELIST, "> $listfn") or
            die "Unable to open file: $listfn\n";

        $nlist = 0;

        foreach $ditherset (@dithers) {                    # foreach dither set
            ($runbeg, $runend) = split(/\s+/, $ditherset);

            $offsetsfn = "./offsets.r${runbeg}to${runend}";
            $objsfn = "./irx.r${runbeg}to${runend}.c$chip.fits.objs";

            open(OFFSET, "$offsetsfn") or 
                die "Unable to open offsets file: $offsetsfn\n";

            while (<OFFSET>) {
                chomp;
                ($runno, $xshift, $yshift, $match) = split;
                $fn = sprintf("./irx_%05d_c%01d_%03d.fits", 
                                                $runno, $chip, $loop);
                if (-e $fn) {
                    print FILELIST "$fn $objsfn $xshift $yshift\n";
                    $nlist += 1;
                }
            }

            close(OFFSET);
        }

        close(FILELIST);

        if ($nlist > 2 * $hwidth) {
            $gain = "./gain.c$chip.fits";
            $cmd = "$skyprog $listfn $gain $hwidth mask $ramp[$chip-1]";

            print "Running sky filter command:\n$cmd\n";

            die "Sorry, command $cmd failed\n"
                unless (system($cmd) == 0);
        } elsif ($nlist > 0) {
            print "WARN: skypass2.pl too few frames, loop $loop, chip $chip\n";
        }
    }
}

if (-e $listfn) { `rm $listfn`; }
