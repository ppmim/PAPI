#!/usr/bin/env perl
#
# median combine all loops for each chip with dome lamp off, and with lamp
# on, and take the difference to form flatfield per chip.  (skip loop 1 due
# to reset anomaly).  call gainmap.c to convert the flatfield to a bad 
# pixel / gain map.
#
# 2-April-2009 : PANIC version (jmiguel@iaa.es) 
#
#
# Notes: 
#   - For the moment, no chip treatment will be done; I hope to run 4 parallel reductions
#
#
#



die "Usage: domeflat.pl on_list off_list" .
    " [nsig nxblock nyblock mingain maxgain]\n" 
    unless ($#ARGV == 1 or $#ARGV == 6);

#@chips   = (1, 2, 3, 4);
@chips   = (1);
$onlist  = $ARGV[0];
$offlist = $ARGV[1];
$|       = 1;                             # turn off output buffering

$basedir = $ENV{'IRDR_BASEDIR'};

die "Please setenv IRDR_BASEDIR /path/irdr\n"
    unless ($basedir);

$datadir = $ENV{'IRDR_DATADIR'};

die "Please setenv IRDR_DATADIR /path/irdr\n"
    unless ($datadir);

$bpmprog = "$basedir/bin/gainmap";
$calcprog = "$basedir/bin/imcalc";
$cubeprog = "$basedir/bin/cubemean";

die "Expected executables: $bpmprog, $calcprog, and $cubeprog\n"
    unless (-x $bpmprog and -x $calcprog and -x $cubeprog);

for ($i = 0; $i <= $#chips; $i++) {                   # foreach chip
    $chip = "c$chips[$i]";

    $non = 0;   $noff = 0;

    chomp($non = `wc $onlist | awk '{print \$1}'`);
  
    chomp($noff = `wc $offlist | awk '{print \$1}'`);
      
    print "Warning: Found no lamp off flats, chip $chip, not making gainmap\n"
        unless ($noff > 0);

    print "Warning: Found no lamp on flats, chip $chip, not making gainmap\n"
        unless ($non > 0);

    if ($noff > 0 and $non > 0) {
        print "\nMedian combine dome lamp on images to make on.fits\n";

        $cmd = "$cubeprog $onlist on.fits NULL offset median noweight float";

        print "Running combine command:\n$cmd\n";

        die "Sorry, command $cmd failed on: $onlist\n"
            unless (system($cmd) == 0);

        print "\nMedian combine dome lamp off images to make off.fits\n";

        $cmd = "$cubeprog $offlist off.fits NULL offset median noweight float";

        print "Running combine command:\n$cmd\n";

        die "Sorry, command $cmd failed on: $offlist\n"
            unless (system($cmd) == 0);

        $flat = "./flat.${chip}.fits";
        $gain = "./gain.${chip}.fits";

        print "\nSubtract lamp off from lamp on to make: $flat\n";

        $cmd  = "$calcprog on.fits '-' off.fits $flat";

        print "Running subtract command:\n$cmd\n";

        die "Sorry, command $cmd failed\n"
            unless (system($cmd) == 0);

        print "\nNormalize $flat and set bad pixels to 0.0 to make: $gain\n";

        if ($#ARGV == 1) {
            $cmd = "$bpmprog $flat $gain";
        } elsif ($#ARGV == 6) {
            $cmd = "$bpmprog $flat $gain $ARGV[2] $ARGV[3] $ARGV[4] $ARGV[5]" .
                   " $ARGV[6]";
        }

        print "Running normalize command:\n$cmd\n";

        die "Sorry, command $cmd failed on: $flat\n"
            unless (system($cmd) == 0);
    }
}

