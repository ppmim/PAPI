#!/usr/bin/env perl
#
# Add accurate FITS WCS header information to FITS images with rough WCS
# using the SCAMP and SExtractor (to produce an x,y,mag list),
# or USNO online
# or 2MASS catalogue previously download from Gator
# (http://irsa.ipac.caltech.edu/applications/Gator/) 
#
# initwcs.c is run first on the image to produce a rough WCS for which RA,DEC
# keywords refer to the chip/image center and not telescope pointing position.
#
# 24-March 2009
# 

die "Usage: astrometry_scamp.pl apm|usno|2mass file*.fits [-u]\n" 
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
$sexprog = "$basedir/scripts/cirsi/sextractor_scamp.pl";
$scamp_cfg = "/disk-a/caha/panic/DEVELOP/PIPELINE/PAPI/scamp.cfg";

die "Expected executables: $apmprog, $usnoprog, $twomassprog\n"
    unless (-x $apmprog and -x $usnoprog and -x $twomassprog);

$initprog = "$basedir/bin/initwcs";
$matchprog = "/usr/local/Terapix/bin/scamp";
$printprog = "$basedir/bin/printwcs";
$fitskeyprog = "$basedir/extern/wcstools/bin/gethead";
$updateheaderprog = "/usr/local/Terapix/bin/missfits";

die "Expected executables: $initprog, $matchprog, $printprog, $fitskeyprog\n"
 unless (-x $initprog and -x $matchprog and -x $printprog and -x $fitskeyprog);

for ($i = 0; $i <= $#ARGV; $i++) {       # loop over files on command line

    $fn = $ARGV[$i];                          # image filename

  
    $rotation=0.0;
    print "Rotation: $rotation\n\n";


    # run SExtractor to create pixel object list

    print `$sexprog $area $thresh $fn`;
    
    
    # check for empty object lists, then match object lists

    die "Sorry, unable to find catalog: $fn.ldac\n"
        unless (-e "$fn.ldac");



    $matchcmd = "$matchprog $fn.ldac -c $scamp_cfg -POSANGLE_MAXERR $rotation -ASTREF_CATALOG 2MASS -SOLVE_PHOTOM N";
    
    if (system($matchcmd) != 0) {
        print "Warning: Command failed: $matchcmd\n";
        next;
    }
    
    # Save backup copy and Update header with new WCS parameters
    print `cp $fn $fn.orig`;             # save backup copy
    print `$updateheaderprog $fn`;
    
    
    exit;

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
