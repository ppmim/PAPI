#!/usr/bin/env perl
#
# call skyfilter.c to do 1st pass sky subtraction.  reads "dithers" file.
# uses loop-combined data.
#
# Jan 2001
#

die "Usage: skypass1.pl half-width\n" 
    unless ($#ARGV == 0);

@chips  = (1, 2, 3, 4);
@ramp   = ("rowcol", "colrow", "rowcol", "colrow");
$hwidth = $ARGV[0];
$listfn = "./skypass1.list";

$basedir = $ENV{'IRDR_BASEDIR'};

die "Please setenv IRDR_BASEDIR /path/irdr\n"
    unless ($basedir);

$skyprog = "$basedir/bin/skyfilter";

die "Expected executable $skyprog\n"
    unless (-x "$skyprog");

die "Unable to read `runs' file\n" 
    unless (-r "runs");

@runs = `cat runs`;

print STDERR "First pass sky subtraction and reset correction\n\n";

foreach $chip (@chips) {                                   # foreach chip
    open(FILELIST, "> $listfn") or
        die "Unable to open file: $listfn\n";

    $nlist = 0;

    foreach $run (@runs) {                                  # foreach run
        $fn = sprintf("./irx_%05d_c%01d.fits", $run, $chip);

        if (-e $fn) {
            print FILELIST "$fn\n";                    # print filename to list
            $nlist += 1;
        }
    }

    close(FILELIST);

    if ($nlist > 0) {
        $gain = "./gain.c$chip.fits";
        $cmd = "$skyprog $listfn $gain $hwidth nomask $ramp[$chip-1]";

        print "Running sky filter command:\n$cmd\n";

        die "Sorry, command $cmd failed\n"
            unless (system($cmd) == 0);               # call skyfilter on list
    }
}

if (-e $listfn) { `rm $listfn`; }
