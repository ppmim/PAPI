#!/usr/bin/env perl
#
# translate apmcat.c output to ra,dec,mag list in decimal degrees
# with format as expected by WCSTools
#

print "APM/j/d/n\n";
print "APM Catalog, J2000, Red magnitudes\n";

while (<>) {
  chomp;
  s/^\s*//;                 # remove leading spaces

  if (! /^\d/) {            # skip line unless starts with number
    next; 
  }

  ($hour, $min, $sec, $deg, $arcmin, $arcsec, $rmag, $class, $sigma,
          $a, $ecc, $pa, $bmag) = split;

  $a     = 0;
  $pa    = 0;
  $bmag  = 0;
  $class = 0;                       # not used
  $sigma = 0;
  $ecc   = 0;

  if ($deg =~ /-/) {                # watch out for -0
    $arcmin = -1 * $arcmin;
    $arcsec = -1 * $arcsec;
  }

  $ra = ($hour / 24. + $min / (60.*24.) + $sec / (3600.*24.)) * 360.;

  $dec = $deg + $arcmin / 60. + $arcsec / 3600.;

  if ($rmag > 5) {                       # print if detected in red
      $rmag = sprintf("%5.1f", $rmag);

      print "$ra $dec $rmag\n";
  }
}
