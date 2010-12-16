#!/usr/local/bin/perl
# Get the required region of 2MASS catalogue from catalogue 
# $ENV{'2MASS_CAT'}
# To retreive the catalogue, use Gator 
# http://irsa.ipac.caltech.edu/applications/Gator/
# select the "2MASS Second Incremental Release Point Source Catalog (PSC)"
# Select "Cone Search" Within 1100 arcsec radius (for ex.), centered on the mosaic center"
# Select the Output Column Selection so that the first columns are:
# ra dec j_m j_msig h_m h_msig k_m k_msig
#
# November 2001

die "Usage: 2mass.pl RAcenter DECcenter Box Filter\n"
    unless ($#ARGV == 3);

$cat2mass = $ENV{'CAT_2MASS'};
$basedir = $ENV{'IRDR_BASEDIR'};

die "Please setenv CAT_2MASS /2mass/file.cat\n"
    unless ($cat2mass);
die "Please setenv IRDR_BASEDIR /path/irdr\n"
    unless ($basedir);

$raX = $ARGV[0];
$decX = $ARGV[1];
$radius = $ARGV[2] / 2. / 60. ;
$filter = $ARGV[3];

$tmpfile = "./2mass.tmp";

if ($filter=~/(J|j)/) { $f=2; }
elsif ($filter=~/(H|h)/) { $f=4; }
elsif ($filter=~/(K|k)/) { $f=6; }

die "Filter $filter unknown\n"
    unless ($f);

open(TMP,">$tmpfile") || die "unable to open $tmpfile\n" ;
open(CAT,"<$cat2mass") || die "unable to open $cat2mass\n" ;
while (<CAT>) {
  if ((! /^\\/)&&(! /^\|/)) {
	chomp;
	@tab = split;
 	$mag = $tab[$f] ; 
	$ra = $tab[0]; $dec = $tab[1];
 	if (($mag!=-99.999)&&($mag ne "null")&&(abs($ra-$raX)<$radius)&&(abs($dec-$decX)<$radius)) {
		print TMP "$ra $dec $mag\n";
	}
  }
}
close (CAT);
close (TMP);

# sort / magnitude and print on STDOUT
print "2MASS/j/d/n\n";
print "2MASS Catalog, J2000, $filter magnitudes\n";
print `sort +2n $tmpfile` ;

if (-e $tmpfile) {
   `rm -f $tmpfile`;
}

exit 0;

