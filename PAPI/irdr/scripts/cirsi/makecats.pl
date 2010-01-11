#!/usr/bin/env perl
#
# Produce SExtractor x,y,ra,dec,mag catalog and APM ra,dec,mag catalog
# for that FOV.  Assumes input FITS files have accurate WCS set.
#
# Jan 2001
#

die "Usage: makecats.pl irx.r*c?.fits\n" 
    unless ($#ARGV >= 0);

sub decdeg;

$area = 10.0;                            # SExtractor DETECT_MINAREA
$thresh = 4.0;                           # SExtractor DETECT_THRESH
$cattype = "apm";

$basedir = $ENV{'IRDR_BASEDIR'};

die "Please setenv IRDR_BASEDIR /path/irdr\n"
    unless ($basedir);

$apmprog = "$basedir/scripts/cirsi/apm.pl";

die "Expected executable: $apmprog\n"
    unless (-x $apmprog);

$usnoprog = "$basedir/scripts/usno.pl";

die "Expected executable: $usnoprog\n"
    unless (-x $usnoprog);

$printprog = "$basedir/bin/printwcs";

die "Expected executable: $printprog\n"
     unless (-x $printprog);

$fitskeyprog = "$basedir/extern/wcstools/bin/gethead";

die "Expected executable: $fitskeyprog\n"
     unless (-x $fitskeyprog);

for ($i = 0; $i <= $#ARGV; $i++) {       # loop over files on command line
    $fn = $ARGV[$i];                          # image filename

    # read field center and size from WCS header

    ($ra, $dec, $scale, $nx, $ny) = split(/\s+/, `$printprog $fn`);

    $ny = 0;                             # unused

    $box = 1.1 * $nx * $scale / 60.;     # assume square image

    print "File: $fn\n";
    print "Scale: $scale, Center: $ra $dec, FOV: $box arcmin\n";

    $radecstr = decdeg($ra, $dec);

    # run apm.pl or usno.pl to extract object list

    if ($cattype eq "usno") {                  # Use USNO A Catalog
        $cmd = "$usnoprog $ra $dec $box > sky.cat";

        die "Sorry, $usnoprog failed\n"
            unless (system($cmd) == 0);
    } else {                                   # Use APM Catalog
        $apmcatalog = ($dec < 0.0) ? "ukst" : "poss";

        $cmd = "$apmprog \"$radecstr\" $box $apmcatalog wcstools > $fn.apmcat";

        die "Sorry, $apmprog failed\n"
            unless (system($cmd) == 0);
    }

    # run SExtractor to create pixel object list

    print `sextractor.pl $area $thresh $fn`;
}

# subroutine to convert decimal degrees to hmsdms
sub decdeg {
    my ($ra, $dec) = @_;

    $hour = int( $ra / 15.0 );
    $min  = int( ($ra - 15.0 * $hour) * 4.0 );
    $sec  = ($ra - 15.0 * $hour - $min / 4.0) * 240.0;

    $deg     = int( $dec );
    $farcmin = abs ( $deg - $dec ) * 60.0;
    $arcmin  = int( $farcmin );
    $arcsec  = ($farcmin - $arcmin) * 60.0;

    if ($deg == 0 && $dec < 0) {
        $radec = sprintf "%02d %02d %4.2f -%02d %02d %3.1f", 
            $hour, $min, $sec, $deg, $arcmin, $arcsec;
    } else {
        $radec = sprintf "%02d %02d %4.2f %02d %02d %3.1f", 
            $hour, $min, $sec, $deg, $arcmin, $arcsec;
    }

    return($radec);
}
