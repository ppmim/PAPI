#!/usr/bin/env perl
#
# Use WCSTools xy2sky program to determine tile coordinates from FITS WCS of
# coadded dither sets.
#
# Jan 2001
#

die "Usage: tile.pl irx.r*c?.fits\n" 
    unless ($#ARGV >= 1);

$basedir = $ENV{'IRDR_BASEDIR'};

die "Please setenv IRDR_BASEDIR /path/irdr\n"
    unless ($basedir);

$wcsprog = "$basedir/extern/wcstools/bin/xy2sky";

die "Expected executable: $wcsprog\n"
     unless (-x $wcsprog);

$headprog = "$basedir/extern/wcstools/bin/gethead";

die "Expected executable: $headprog\n"
     unless (-x $headprog);

for ($i = 0; $i <= $#ARGV; $i++) {       # loop over files on command line
    $fn = $ARGV[$i];                          # image filename
    
    ($nx, $ny) = split(/\s+/, `$headprog $fn NAXIS1 NAXIS2`);

    chomp($scale = `$headprog $fn SECPIX1`);

    chomp($scale);

    ($j, $j, $j, $ra1, $dec0) = split(/\s+/, `$wcsprog -j -d $fn 0 0`);
    ($j, $j, $j, $ra0, $dec1) = split(/\s+/, `$wcsprog -j -d $fn $nx $ny`);

    if ($i == 0) {
        $ramin = $ra0;
        $decmin = $dec0;
        $ramax = $ra1;
        $decmax = $dec1;
    } else {
        $ramin = ($ramin < $ra0) ? $ramin : $ra0;
        $decmin = ($decmin < $dec0) ? $decmin : $dec0;
        $ramax = ($ramax > $ra1) ? $ramax : $ra1;
        $decmax = ($decmax > $dec1) ? $decmax : $dec1;
    }

    print "$fn $scale $ra0 $ra1 $dec0 $dec1\n";
}

# print "$ramin $ramax $decmin $decmax\n";

$xcen = $ramin + 0.5 * ($ramax - $ramin);
$ycen = $decmin + 0.5 * ($decmax - $decmin);
print "RA center: $xcen\n";
print "DEC center: $ycen\n";

print "scale is $scale\n";
$xsize = int(($ramax - $ramin) * 3600.0 / $scale + 0.5);
$ysize = int(($decmax - $decmin) * 3600.0 / $scale + 0.5);
$xysize = ($xsize + $ysize) / 2.;
$xysize=0;
# print "$xsize $ysize $xysize\n";
