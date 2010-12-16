#!/usr/bin/perl -w
#
# translate apmcat.c output to ra,dec list in decimal degrees.
# use standard catalog format as expected by ESO SkyCat
#
# CNS, 27 Jan 2000
#

print "id  \t  ra  \t  dec  \t  mag  \t  a  \t  a/b  \t pa\n";
print "--  \t  --  \t  ---  \t  ---  \t  -  \t  ---  \t --\n";

$count = 0;

while (<>) {
  chomp;
  s/^\s*//;                 # remove leading spaces

  if (! /^\d/) {            # skip line unless starts with number
    next; 
  }

  ($hour, $min, $sec, $deg, $arcmin, $arcsec, $rmag, $class, $sigma,
          $a, $ecc, $pa, $bmag) = split;

  $pa    = $pa - 90.;               # position angle from North
  $b     = $a * (1.0 - $ecc);       # minor axis
  $bmag  = 0;
  $class = 0;                       # not used
  $sigma = 0;

  if ($deg =~ /-/) {                # watch out for -0
    $arcmin = -1 * $arcmin;
    $arcsec = -1 * $arcsec;
  }

  $ra = ($hour / 24. + $min / (60.*24.) + $sec / (3600.*24.)) * 360.;

  $dec = $deg + $arcmin / 60. + $arcsec / 3600.;

  if ($rmag > 5) {                  # require detection in red plate
      $ratio = $a / $b;

      $rmag = sprintf("%5.1f", $rmag);

      print "$count  \t  $ra  \t  $dec  \t $rmag  \t $a \t $ratio \t $pa\n";

      $count += 1;
  }

}
