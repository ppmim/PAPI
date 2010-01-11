#!/usr/bin/env perl
#
# Get the required region of USNO A from remote server as plain text
# and print it to standard output.  Input field center coordinates 
# are in degrees, and field diameter (box size) is in arcmin.  The
# executable lynx (ascii web browser) must be in your path.
#
# Jan 2001
#

die "Usage: usno.pl RAcenter DECcenter Box\n"
    unless ($#ARGV == 2);

sub decdeg;

$url = "http://archive.eso.org/skycat/servers/usnoa-server?";

$lynx = "lynx -dump ";

$ra = $ARGV[0];

$dec = $ARGV[1];

$radius = $ARGV[2] / 2.;

$arg1 = decdeg($ra, $dec);

$arg2 = sprintf("&radius=0.0,%7.2f&sort=mr", $radius);

$req = $arg1 . $arg2;

$req =~ s/ //g;

$cmd = $lynx . "'" . $url . $req . "'" . " | awk '{print \$2, \$3, \$4}'";

open(CMD, "$cmd |") or die "unable to run cmd: $cmd\n";

print "USNO/j/d/n\n";
print "USNO Catalog, J2000, Red magnitudes\n";

while (<CMD>) {
  ($radeg, $decdeg, $rmag) = split;

  if ($radeg and $decdeg and $rmag and $radeg =~ /^\d/) {
      print "$radeg $decdeg $rmag\n";
  }
}

close(CMD);

# convert decimal degrees to hmsdms

sub decdeg {
    my ($ra, $dec) = @_;

    $hour = int( $ra / 15.0 );
    $min  = int( ($ra - 15.0 * $hour) * 4.0 );
    $sec  = ($ra - 15.0 * $hour - $min / 4.0) * 240.0;

    $deg     = int( $dec );
    $farcmin = abs ( $deg - $dec ) * 60.0;
    $arcmin  = int( $farcmin );
    $arcsec  = ($farcmin - $arcmin) * 60.0;

    if ($dec < 0) {
        if ($deg == 0) {
            $radec = sprintf "%02d:%02d:%4.2f -%02d:%02d:%3.1f", 
                $hour, $min, $sec, $deg, $arcmin, $arcsec;
        } else {
            $radec = sprintf "%02d:%02d:%4.2f %02d:%02d:%3.1f", 
                $hour, $min, $sec, $deg, $arcmin, $arcsec;
        }
    } else {
        $radec = sprintf "%02d:%02d:%4.2f +%02d:%02d:%3.1f",
            $hour, $min, $sec, $deg, $arcmin, $arcsec;
    }

    return($radec);
}
