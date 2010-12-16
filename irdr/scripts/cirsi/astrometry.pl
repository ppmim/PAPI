#!/usr/bin/env perl
#
# Add accurate FITS WCS header information to FITS images with rough WCS
# using the WCSTools imwcs.c program, SExtractor (to produce an x,y,mag list),
# and apmcat.c (to extract an ra,dec,mag list from the APM online catalog)
# or USNO online
# or 2MASS catalogue previously download from Gator
# (http://irsa.ipac.caltech.edu/applications/Gator/) 
# See 2mass.pl for more details
#
# initwcs.c is run first on the image to produce a rough WCS for which RA,DEC
# keywords refer to the chip/image center and not telescope pointing position.
#
# Jan 2001
# May 2002 : add 2mass option
#

die "Usage: astrometry.pl apm|usno|2mass irx.r*c?.fits\n" 
    unless ($#ARGV >= 1);

$area = 10.0;                            # SExtractor DETECT_MINAREA
$thresh = 4.0;                           # SExtractor DETECT_THRESH
$tol = 2.5;                              # offset tolerance in arcsec

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
$sexprog = "$basedir/scripts/cirsi/sextractor.pl";

die "Expected executables: $apmprog, $usnoprog, $twomassprog\n"
    unless (-x $apmprog and -x $usnoprog and -x $twomassprog);

$initprog = "$basedir/bin/initwcs";
$matchprog = "$basedir/extern/wcstools/bin/imwcs";
$printprog = "$basedir/bin/printwcs";
$fitskeyprog = "$basedir/extern/wcstools/bin/gethead";

die "Expected executables: $initprog, $matchprog, $printprog, $fitskeyprog\n"
 unless (-x $initprog and -x $matchprog and -x $printprog and -x $fitskeyprog);

for ($i = 0; $i <= $#ARGV; $i++) {       # loop over files on command line

    $fn = $ARGV[$i];                          # image filename

    # run initwcs.c to initialize rough WCS header

    $initcmd = "$initprog $fn";

    die "Command $initcmd failed\n"
        unless (system($initcmd) == 0);

    # make backup copy of coadded dither set with rough WCS header

#    print `cp -f $fn $fn.initwcs`;

    # read field center and size from rough WCS header

    ($ra, $dec, $scale, $nx, $ny) = split(/\s+/, `$printprog $fn`);

    $ny = 0;                             # unused

    $box = 1.2 * $nx * $scale / 60.;     # assume square image

    print "File: $fn\n";

    print "Scale: $scale, Center: $ra $dec, FOV: $box arcmin Nx:$nx Ny:$ny \n";

    $radecstr = decdeg($ra, $dec);

    # read expected chip rotation from rough WCS header

    chomp($rotation = `$fitskeyprog $fn CHIPROT`);

    if (!$rotation or $rotation < -5.0 or $rotation > 5.0) {
        $rotation = 0.0;
    }

    print "Rotation: $rotation\n\n";

    # run apm.pl or usno.pl or 2mass.pl to extract object list

    if ($cattype eq "usno") {                  # Use USNO A Catalog
        $cmd = "$usnoprog $ra $dec $box > sky.cat";

        die "Sorry, $usnoprog failed\n"
            unless (system($cmd) == 0);

    } elsif ($cattype eq "apm") {              # Use APM Catalog
        $apmcatalog = ($dec < 0.0) ? "ukst" : "poss";

        $cmd = "$apmprog \"$radecstr\" $box $apmcatalog wcstools > sky.cat";

        die "Sorry, $apmprog failed\n"
            unless (system($cmd) == 0);
	
    } else {					# Use 2MASS Catalog
	chomp($filter = `$fitskeyprog $fn FILTER`);
	$cmd = "$twomassprog $ra $dec $box $filter > sky.cat";

        die "Sorry, $twomassprog failed\n"
            unless (system($cmd) == 0);
   }


    # run SExtractor to create pixel object list

    print `$sexprog $area $thresh $fn`;
    print `sort +2n $fn.cat > $fn.cat.sort` ;
    print `mv -f $fn.cat.sort $fn.cat` ;

    # check for empty object lists, then match object lists

    die "Sorry, unable to find catalog: $fn.cat\n"
        unless (-e "$fn.cat");

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

    $matchcmd = "$matchprog -v -p $scale -c sky.cat -d $fn.cat -o $fn.astrom "
              . "-t $tolpix -h $nref -a $rotation -q irps -n -8 "
              . "$fn > $fn.imwcs 2>&1";

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
