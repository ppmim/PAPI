#!/usr/bin/perl -w
#
# run this script after astrometry.pl. this script uses the first pass 
# WCS parameters to run WCSTools again with a smaller clipping radius.
#
# March 2001
# Mai 2003: debug. add 2mass option

die "Usage: astrometry_refine.pl apm|usno|2mass irx.r*c?.fits\n" 
    unless ($#ARGV >= 1);

$area   = 10.0;                           # SExtractor DETECT_MINAREA
$thresh = 4.0;                           # SExtractor DETECT_THRESH
$tol    = 2.5;                           # offset tolerance in arcsec
#@rots   = (-1.0,-1.0,-1.0,-1.1);

$cattype = shift;

die "Catalog type should be `apm' or `usno' or `2mass'\n"
    unless ($cattype eq "apm" or $cattype eq "usno" or $cattype eq "2mass");

sub decdeg;

$basedir = $ENV{'IRDR_BASEDIR'};

die "Please setenv IRDR_BASEDIR /path/irdr\n"
    unless ($basedir);

$apmprog = "$basedir/scripts/cirsi/apm.pl";
$usnoprog = "$basedir/scripts/cirsi/usno.pl";
$twomassprog = "$basedir/scripts/cirsi/2mass.pl";

die "Expected executables: $apmprog, $usnoprog, $twomassprog\n"
    unless (-x $apmprog and -x $usnoprog and -x $twomassprog);

$sexprog = "$basedir/scripts/cirsi/sextractor.pl";
$initprog = "$basedir/bin/initwcs";
$matchprog = "$basedir/extern/wcstools/bin/imwcs";
$printprog = "$basedir/bin/printwcs";
$fitskeyprog = "$basedir/extern/wcstools/bin/gethead";

die "Expected executables: $initprog, $matchprog, $printprog, $fitskeyprog, $sexprog\n"
 unless (-x $initprog and -x $matchprog and -x $printprog and -x $fitskeyprog and -x $sexprog);

for ($i = 0; $i <= $#ARGV; $i++) {       # loop over files on command line

    $fn = $ARGV[$i];                          # image filename

    # read field center and size from rough WCS header

    ($ra, $dec, $scale, $nx, $ny) = split(/\s+/, `$printprog $fn`);

    $ny = 0;                              # unused
    $box = 1.2 * $nx * $scale / 60.;     # assume square image

    chomp($radecstr = `$fitskeyprog $fn RA DEC`);

    # read expected chip rotation from rough WCS header

#    $fn =~ /c(.)\.fits/;
#    $chipno = $1;
#    $rotation = $rots[$chipno-1];

    chomp($rotation = `$fitskeyprog $fn CHIPROT`);
    if (!$rotation or $rotation < -5.0 or $rotation > 5.0) {
        $rotation = 0.0;
    }

    print "File: $fn\n";
    print "Scale: $scale, Center: $radecstr, FOV: $box arcmin\n";
    print "Rotation: $rotation\n\n";

    # run apm.pl or usno.pl or 2mass.pl to extract object list

    if ($cattype eq "usno") {                  # Use USNO A Catalog
        $cmd = "$usnoprog $ra $dec $box > sky.cat";

        die "Sorry, $usnoprog failed\n"
            unless (system($cmd) == 0);

    } elsif ( $cattype eq "apm") {             # Use APM Catalog
        $apmcatalog = ($dec < 0.0) ? "ukst" : "poss";

        $cmd = "$apmprog \"$radecstr\" $box $apmcatalog wcstools > sky.cat";

        die "Sorry, $apmprog failed\n"
            unless (system($cmd) == 0);
    } else {
                chomp($filter = `$fitskeyprog $fn FILTER`);
        $cmd = "$twomassprog $ra $dec $box $filter > sky.cat";

        die "Sorry, $twomassprog failed\n"
            unless (system($cmd) == 0);
    }

    # run SExtractor to create pixel object list

    print `$sexprog $area $thresh $fn`;

    # check for empty object lists, then match object lists

    die "Sorry, unable to find catalog: $fn.cat\n"
        unless (-e "$fn.cat");

    # sort SExtractor catalog by magnitude

    print `sort +2n $fn.cat > $fn.cat.sort`;
    print `mv -f $fn.cat.sort $fn.cat`;

    die "Sorry, unable to find catalog: sky.cat\n"
        unless (-e "sky.cat");

    chomp($nobjs1 = `wc $fn.cat | awk '{print \$1}'`);
    chomp($nobjs2 = `wc sky.cat | awk '{print \$1}'`);

    die "SExtractor catalog output has too few objects\n"
        unless ($nobjs1 > 5);

    die "Empty sky catalog output\n"
        unless ($nobjs2 > 2);

    $nref = ($nobjs2 > 100) ? 100 : $nobjs2;
    $tolpix = $tol / $scale;

    $matchcmd = "$matchprog -v -c sky.cat -d $fn.cat "
              . " -o $fn.astrom -t $tolpix -h $nref -q irps "
              . " $fn > $fn.imwcs 2>&1";

    if (system($matchcmd) != 0) {
        print "Warning: Command failed: $matchcmd\n";
        next;
    }

    print `grep "^#" $fn.imwcs`;
    print `grep resid $fn.imwcs`;
    print `grep "No reference stars" $fn.imwcs`;
    print `grep Numerical $fn.imwcs`;

    $_ = `grep points $fn.imwcs`;
    ($nbpoints)=split();
    if($nbpoints<10) {
        print "WARN: only $nbpoints found for the fit of $fn\n";
    }


    print `mv $fn $fn.initwcs`;             # save backup copy
    print `mv $fn.astrom $fn`;

    $ra_orig = $ra;
    $dec_orig = $dec;

    ($ra, $dec, $scale, $nx, $ny) = split(/\s+/, `$printprog $fn`);

    $delra = $ra - $ra_orig;
    $deldec = $dec - $dec_orig;

    print "DRA $delra, DDEC $deldec\n";
}
