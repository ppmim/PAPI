#!/usr/bin/env perl
#
# PAPI Astrometry with SCAMP 
#
# Add accurate FITS WCS header information to FITS images with rough WCS
# using the SCAMP and SExtractor (to produce catalog .ldac)
#
# initwcs.c is run first on the image to produce a rough WCS for which RA,DEC
# keywords refer to the chip/image center and not telescope pointing position.
#
# Created: 24-March 2009
# Last Update: 12-Feb-2010

die "Usage: astrometry_scamp.pl apm|usno|2mass regrid|noregrid file*.fits [-u]\n" 
    unless ($#ARGV >= 1);

$area = 10.0;                            # SExtractor DETECT_MINAREA
$thresh = 4.0;                           # SExtractor DETECT_THRESH
$tol = 2.5;                              # offset tolerance in arcsec

$cattype = shift;

die "Catalog type should be `usno' or `2mass' or `ucac' or `sdss'\n"
    unless ($cattype eq "usno" or $cattype eq "ucac" or $cattype eq "2mass" or $cattype eq "sdss");

$regridtype = shift;

die "Regrid type should be `regrid' or `noregrid' \n"
    unless ($regridtype eq "regrid" or $regridtype eq "noregrid" );


sub decdeg;

$basedir = $ENV{'IRDR_BASEDIR'};
$papi_home = $ENV{'PAPI_HOME'};
$terapix_home = $ENV{'TERAPIX'};

die "Please setenv PAPI_HOME /path/papi\n"
    unless ($papi_home);

$sexprog = "$papi_home/irdr/scripts/cirsi/sextractor_scamp.pl";
$scamp_cfg = "$papi_home/config/scamp.cfg";
$initprog = "$papi_home/irdr/bin/initwcs";
$matchprog = "$terapix_home/scamp";
$printprog = "$papi_home/irdr/bin/printwcs";
$fitskeyprog = "$papi_home/irdr/extern/wcstools/bin/gethead";
$updateheaderprog = "$terapix_home/missfits ";
$regridprog = "$terapix_home/swarp ";

die "Expected executables: $initprog, $sexprog, $matchprog, $printprog, $fitskeyprog\n"
 unless (-x $initprog and -x $matchprog and -x $printprog and -x $fitskeyprog and -x $sexprog);

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


    # run SExtractor to create pixel object list

    print `$sexprog $area $thresh $fn`;
    
    
    # check for empty object lists, then match object lists

    die "Sorry, unable to find catalog: $fn.ldac\n"
        unless (-e "$fn.ldac");


#    chomp($nobjs1 = `wc $fn.ldac | awk '{print \$1}'`);


#    die "SExtractor catalog output has too few objects\n"
#        unless ($nobjs1 > 5);

#    die "Empty sky catalog output\n"
#        unless ($nobjs2 > 2);

#    $nref = ($nobjs2 > 100) ? 100 : $nobjs2;
    


    if ($cattype eq "2mass"){ 
        $CAT="2MASS";
    }elsif ($cattype eq "usno"){
        $CAT="USNO-B1";
    }elsif ($cattype eq "ucac"){
        $CAT="UCAC-2";
    }elsif ($cattype eq "sdss"){
        $CAT="SDSS-R5";
    } else {
        $CAT="2MASS";
    }

#    $matchcmd = "$matchprog $fn.ldac -c $scamp_cfg -POSANGLE_MAXERR $rotation -ASTREF_CATALOG $CAT -SOLVE_PHOTOM N -CHECKPLOT_TYPE NONE";
    $matchcmd = "$matchprog $fn.ldac -c $scamp_cfg -POSANGLE_MAXERR $rotation -ASTREF_CATALOG $CAT -SOLVE_PHOTOM N -WRITE_XML N";    
    if (system($matchcmd) != 0) {
        print "Warning: Command failed: $matchcmd\n";
        next;
    }
    
    
    
    # Save backup copy and Update header with new WCS parameters
    print `cp $fn $fn.orig`;             # save backup copy
    $fn_= $fn;
    $fn_=~ s/.fits/.head/;
    print `cp $fn.head $fn_`;
    print `$updateheaderprog $fn -WRITE_XML N`;
    
    # Re-grid the image with the new WCS
    if ($regridtype eq "regrid"){ 
        print `$regridprog $fn -WRITE_XML N`;
    }

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
