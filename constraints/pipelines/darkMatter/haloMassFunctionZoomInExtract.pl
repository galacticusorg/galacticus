#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use XML::Simple;
use PDL;
use PDL::NiceSlice;

# Locate and extract properties of the primary halo in the zoom-in simulation.
# Andrew Benson (15-March-2022)

# Get arguments.
die("Usgae: haloMassFunctionZoomInExtract.pl <pathName> <primaryHaloFileName> <expansionFactor> <hubbleConstant> <massParticle>")
    unless ( scalar(@ARGV) == 5 );
my $pathName            =     $ARGV[0];
my $primaryHaloFileName =     $ARGV[1];
my $expansionFactor     =     $ARGV[2];
my $hubbleConstant      = pdl $ARGV[3];
my $massParticle        = pdl $ARGV[4];

# Extract box and region size from the Music config file.
print "\tExtracting high-resolution region\n";
my $musicFileName = $pathName."music.conf";
my $boxSizeOriginal;
my $regionExtent;
open(my $musicFile,$musicFileName);
while ( my $line = <$musicFile> ) {
    if ( $line =~ m/^boxlength\s*=\s*([\d\.]+)/ ) {
	$boxSizeOriginal = pdl $1;
    }
    if ( $line =~ m/^ref_extent\s*=\s*([\d\.\,\s]+)/ ) {
	my @extent = split(/\s*,\s*/,$1);
	$regionExtent = pdl @extent;
    }
}
close($musicFile);
$boxSizeOriginal /= $hubbleConstant;
my $boxSize       = $boxSizeOriginal*$regionExtent->prodover()->pow(1.0/3.0);
# Extract the name of the tree file being processed.
my $treeFileName = $pathName."tree_0_0_0.dat";
# Read the merger tree file looking for the position and virial radius of the most massive halo.
print "\tReading tree file\n";
system("awk '{if (substr(\$1,1,1) != \"#\" && NF > 1) print \$1,\$2,\$4,\$11,\$12,\$15,\$18,\$19,\$20}' ".$treeFileName." > ".$pathName."treeReduced.dat");
(my $haloExpansionFactor, my $haloID, my $descID, my $massHalo, my $radiusHalo, my $mostMassiveProgenitor, my $x, my $y, my $z) = rcols($pathName."treeReduced.dat",0,1,2,3,4,5,6,7,8, {CHUNKSIZE => 100000000});
# Convert units to Galacticus standards.
print "\tConverting units\n";
$massHalo       /=        $hubbleConstant; # Convert Msun/h to Msun.
$radiusHalo     *= 1.0e-3/$hubbleConstant; # Convert kpc/h to Mpc.
$x              /=        $hubbleConstant; # Convert Mpc/h to Mpc.
$y              /=        $hubbleConstant; # Convert Mpc/h to Mpc.
$z              /=        $hubbleConstant; # Convert Mpc/h to Mpc.
# Find region center.
my $selection = which($haloExpansionFactor > 0.9999);
my $xMinimum = $x->($selection)->minimum();
my $xMaximum = $x->($selection)->maximum();
my $yMinimum = $y->($selection)->minimum();
my $yMaximum = $y->($selection)->maximum();
my $zMinimum = $z->($selection)->minimum();
my $zMaximum = $z->($selection)->maximum();
my $xCenter  = 0.5*($xMaximum+$xMinimum);
my $yCenter  = 0.5*($yMaximum+$yMinimum);
my $zCenter  = 0.5*($zMaximum+$zMinimum);
# Find focus halo.
print "\tLocating focus halo\n";
my $halosCurrent           = which(($haloExpansionFactor > 0.9999) & (log10($massHalo) > 12.0) & (log10($massHalo) < 12.3));
die("ERROR: unable to locate any viable halos")
    if ( nelem($halosCurrent) <= 0 );
my $distanceFromCenter     = sqrt(+($x-$xCenter)**2+($y-$yCenter)**2+($z-$zCenter)**2);
my $indexHalosCurrentFocus = $distanceFromCenter->($halosCurrent)->minimum_ind();
my $indexHaloFocus         = $halosCurrent->(($indexHalosCurrentFocus));
# Determine the upper mass limit.
my $massCentral            = $massHalo->($indexHaloFocus);
# Find the primary progenitor at this expansion factor.
print "\tLocating primary progenitor\n";
my $indexPrimaryProgenitor = $indexHaloFocus;
while ( abs($haloExpansionFactor->($indexPrimaryProgenitor)-$expansionFactor) > 1.0e-3 ) {
    my $selection =
	which
	(
	 ($descID                == $haloID->($indexPrimaryProgenitor))
	 &
	 ($mostMassiveProgenitor == 1                                 )
	);
    if      ( nelem($selection) < 1 ) {
	die('no progenitor found' );
    } elsif ( nelem($selection) > 1 ) {
	die('ambiguous progenitor');
    } else {
	$indexPrimaryProgenitor = $selection->((0));
    }
}
# Store the relevant data.
my $primaryHaloData;
$primaryHaloData->{'x'   } = $x         ->($indexPrimaryProgenitor)->sclr();
$primaryHaloData->{'y'   } = $y         ->($indexPrimaryProgenitor)->sclr();
$primaryHaloData->{'z'   } = $z         ->($indexPrimaryProgenitor)->sclr();
$primaryHaloData->{'r'   } = $radiusHalo->($indexPrimaryProgenitor)->sclr();
$primaryHaloData->{'m'   } = $massHalo  ->($indexPrimaryProgenitor)->sclr();
$primaryHaloData->{'l'   } = $boxSize                              ->sclr();
$primaryHaloData->{'xc'  } = $xCenter                              ->sclr();
$primaryHaloData->{'yc'  } = $yCenter                              ->sclr();
$primaryHaloData->{'zc'  } = $zCenter                              ->sclr();
$primaryHaloData->{'mc'  } = $massCentral                          ->sclr();
my $xml = new XML::Simple();
open(my $primaryHaloFile,">",$primaryHaloFileName);
print $primaryHaloFile $xml->XMLout($primaryHaloData, RootName => "primaryHalo");
close($primaryHaloFile);

exit 0;
