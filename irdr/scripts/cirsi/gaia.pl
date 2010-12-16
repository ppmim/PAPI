#!/usr/bin/perl -w
#
# translate sextract output to ra,dec list in decimal degrees
# use standard catalog format as expected by ESO SkyCatGaia
#

print "id  \t  ra  \t  dec  \t  mag  \t  a  \t  a/b  \t pa\n";
print "--  \t  --  \t  ---  \t  ---  \t  -  \t  ---  \t --\n";

while (<>) {
  chomp;
  s/^\s*//;                 # remove leading spaces

  if (! /^\d/) {            # skip line unless starts with number
    next; 
  }

  ($id, $ra, $dec, $mag, $a, $ellip, $pa) = split;

  $b = $a * (1.0 - $ellip);

  $ratio = $a / $b;

  $mag = $mag + 30.;

  print "$id \t  $ra  \t  $dec  \t $mag  \t $a \t $ratio \t $pa\n";
}
