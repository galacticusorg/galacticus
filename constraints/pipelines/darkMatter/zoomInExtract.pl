#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use XML::Simple;
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;

# Locate and extract properties of the primary halo in the zoom-in simulation.
# Andrew Benson (15-March-2022)

# Get arguments.
die("Usage: zoomInExtract.pl <pathName> <primaryHaloFileName> <expansionFactor> <hubbleConstant> <massParticle> <massHostLogMin> <massHostLogMax> <hostHaloID>")
    unless ( scalar(@ARGV) == 8 );
my $pathName            =     $ARGV[0];
my $primaryHaloFileName =     $ARGV[1];
my $expansionFactor     =     $ARGV[2];
my $hubbleConstant      = pdl $ARGV[3];
my $massParticle        = pdl $ARGV[4];
my $massHostLogMin      = pdl $ARGV[5];
my $massHostLogMax      = pdl $ARGV[6];
my $hostHaloID          =     $ARGV[7];

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
my $treeFileName = $pathName."tree_?_?_?.dat";
# Read the merger tree file looking for the position and virial radius of the most massive halo.
print "\tPre-processing tree file\n";
system("awk '{if (substr(\$1,1,1) != \"#\" && NF > 1 && \$15 == 1) print \$1,\$2,\$4,\$11,\$12,\$18,\$19,\$20}' ".$treeFileName." > ".$pathName."treeReduced.dat")
    unless ( -e $pathName."treeReduced.dat" );
system("wc -l ".$pathName."treeReduced.dat > ".$pathName."treeCount.dat")
    unless ( -e $pathName."treeCount.dat" );
open(my $treeMetaData,$pathName."treeCount.dat");
my $wc = <$treeMetaData>;
my @wcs = split(" ",$wc);
my $treeCount = $wcs[0];
close($treeMetaData);
print "\t\tPre-processed tree contains ".$treeCount." halos\n";
print "\tReading tree file\n";
my $haloExpansionFactor;
my $haloID             ;
my $descID             ;
my $massHalo           ;
my $radiusHalo         ;
my $x                  ;
my $y                  ;
my $z                  ; 
if ( -e $pathName."treeReduced.hdf5" ) {
    my $treeRaw = new PDL::IO::HDF5($pathName."treeReduced.hdf5");
    $haloExpansionFactor = $treeRaw->dataset('haloExpansionFactor')->get();
    $haloID              = $treeRaw->dataset('haloID'             )->get();
    $descID              = $treeRaw->dataset('descID'             )->get();
    $massHalo            = $treeRaw->dataset('massHalo'           )->get();
    $radiusHalo          = $treeRaw->dataset('radiusHalo'         )->get();
    $x                   = $treeRaw->dataset('x'                  )->get();
    $y                   = $treeRaw->dataset('y'                  )->get();
    $z                   = $treeRaw->dataset('z'                  )->get();
} else {
    $haloExpansionFactor = pdl zeroes($treeCount);
    $haloID              = pdl zeroes($treeCount);
    $descID              = pdl zeroes($treeCount);
    $massHalo            = pdl zeroes($treeCount);
    $radiusHalo          = pdl zeroes($treeCount);
    $x                   = pdl zeroes($treeCount);
    $y                   = pdl zeroes($treeCount);
    $z                   = pdl zeroes($treeCount);
    open(my $tree,$pathName."treeReduced.dat");
    for(my $i=0;$i<$treeCount;++$i) {
	my $line = <$tree>;
	my @columns = split(" ",$line);
	$haloExpansionFactor->(($i)) .= $columns[0];
	$haloID             ->(($i)) .= $columns[1];
	$descID             ->(($i)) .= $columns[2];
	$massHalo           ->(($i)) .= $columns[3];
	$radiusHalo         ->(($i)) .= $columns[4];
	$x                  ->(($i)) .= $columns[5];
	$y                  ->(($i)) .= $columns[6];
	$z                  ->(($i)) .= $columns[7];
    }
    close($tree);
    my $treeRaw = new PDL::IO::HDF5(">".$pathName."treeReduced.hdf5");
    $treeRaw->dataset('haloExpansionFactor')->set($haloExpansionFactor);
    $treeRaw->dataset('haloID'             )->set($haloID             );
    $treeRaw->dataset('descID'             )->set($descID             );
    $treeRaw->dataset('massHalo'           )->set($massHalo           );
    $treeRaw->dataset('radiusHalo'         )->set($radiusHalo         );
    $treeRaw->dataset('x'                  )->set($x                  );
    $treeRaw->dataset('y'                  )->set($y                  );
    $treeRaw->dataset('z'                  )->set($z                  );
}

# Convert units to Galacticus standards.
print "\tConverting units\n";
$massHalo       /=        $hubbleConstant; # Convert Msun/h to Msun.
$radiusHalo     *= 1.0e-3/$hubbleConstant; # Convert kpc/h to Mpc.
$x              /=        $hubbleConstant; # Convert Mpc/h to Mpc.
$y              /=        $hubbleConstant; # Convert Mpc/h to Mpc.
$z              /=        $hubbleConstant; # Convert Mpc/h to Mpc.
# Find region center.
my $selection = which($haloExpansionFactor > 0.9999);
my $xCenter   = sum($x->($selection)*$massHalo->($selection))/sum($massHalo->($selection));
my $yCenter   = sum($y->($selection)*$massHalo->($selection))/sum($massHalo->($selection));
my $zCenter   = sum($z->($selection)*$massHalo->($selection))/sum($massHalo->($selection));

# Find focus halo unless we have already been given a host halo ID.
print "\tLocating focus halo\n";
my $indexHaloFocus;
if ( $hostHaloID < 0 ) {
    my $halosCurrent = which(($haloExpansionFactor > 0.9999) & (log10($massHalo) > $massHostLogMin) & (log10($massHalo) < $massHostLogMax));
    if ( nelem($halosCurrent) <= 0 ) {
	# No halo meets the criterion. Simply find the most massive halo at z=0.
	undef($halosCurrent);
	my $halosPresent = which($haloExpansionFactor > 0.9999);
	die("ERROR: unable to locate any viable halos")
	    if ( nelem($halosPresent) <= 0 );
	my $orderMass = $massHalo->($halosPresent)->qsorti();
	$halosCurrent = $halosPresent->($orderMass)->(-1);
    }
    my $distanceFromCenter     = sqrt(+($x-$xCenter)**2+($y-$yCenter)**2+($z-$zCenter)**2);
    my $indexHalosCurrentFocus = $distanceFromCenter->($halosCurrent)->minimum_ind();
    $indexHaloFocus            = $halosCurrent->(($indexHalosCurrentFocus));
} else {
    my $halosCurrent = which($haloID == $hostHaloID);
    die('unable to locate host halo')
	if ( nelem($halosCurrent) < 1 );
    die('ambiguous host halo')
	if ( nelem($halosCurrent) > 1 );
    $indexHaloFocus = $halosCurrent->((0));
}
# Determine the upper mass limit.
my $massCentral            = $massHalo->($indexHaloFocus);
# Find the primary progenitor at this expansion factor.
print "\tLocating primary progenitor\n";
my $indexPrimaryProgenitor = $indexHaloFocus;
while ( abs($haloExpansionFactor->(($indexPrimaryProgenitor))-$expansionFactor) > 1.0e-3 ) {
    my $selection =
	which
	(
	 $descID == $haloID->(($indexPrimaryProgenitor))
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
$primaryHaloData->{'i'   } = $haloID    ->(($indexPrimaryProgenitor))->sclr();
$primaryHaloData->{'x'   } = $x         ->(($indexPrimaryProgenitor))->sclr();
$primaryHaloData->{'x'   } = $x         ->(($indexPrimaryProgenitor))->sclr();
$primaryHaloData->{'y'   } = $y         ->(($indexPrimaryProgenitor))->sclr();
$primaryHaloData->{'z'   } = $z         ->(($indexPrimaryProgenitor))->sclr();
$primaryHaloData->{'r'   } = $radiusHalo->(($indexPrimaryProgenitor))->sclr();
$primaryHaloData->{'m'   } = $massHalo  ->(($indexPrimaryProgenitor))->sclr();
$primaryHaloData->{'l'   } = $boxSize                                ->sclr();
$primaryHaloData->{'xc'  } = $xCenter                                ->sclr();
$primaryHaloData->{'yc'  } = $yCenter                                ->sclr();
$primaryHaloData->{'zc'  } = $zCenter                                ->sclr();
$primaryHaloData->{'mc'  } = $massCentral                            ->sclr();
my $xml = new XML::Simple();
open(my $primaryHaloFile,">",$primaryHaloFileName);
print $primaryHaloFile $xml->XMLout($primaryHaloData, RootName => "primaryHalo");
close($primaryHaloFile);

exit 0;
