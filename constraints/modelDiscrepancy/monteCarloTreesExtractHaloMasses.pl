#!/usr/bin/env perl
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V094"}) ) {
    $galacticusPath = $ENV{"GALACTICUS_ROOT_V094"};
    $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
    $galacticusPath = "./";
}
unshift(@INC, $galacticusPath."perl"); 
use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
require Galacticus::HDF5;

# Extract final time halo masses from a Galacticus model and write them to file to allow future runs to use the exact same set of
# halo masses.
# Andrew Benson (10-October-2013)

# Get arugments.
die("Usage: monteCarloTreesExtractHaloMasses.pl <modelDirectory>")
    unless ( scalar(@ARGV) == 1 );
my $modelDirectory = $ARGV[0];

# Construct the Galaticus model file name.
my $galacticusFileName = $modelDirectory."/galacticus.hdf5";
# Read the required data.
my $galacticus;
$galacticus->{'file' } = $galacticusFileName;
$galacticus->{'store'} = 0;
$galacticus->{'tree' } = "all";
&HDF5::Select_Output($galacticus,0.0);
&HDF5::Get_Dataset  ($galacticus,['mergerTreeWeight','basicMass','nodeIsIsolated','mergingStatisticsNodeHierarchyLevel']);
my $dataSets = $galacticus->{'dataSets'};
# Find isolated halos which have never been subhalos.
my $isolated =
    which(
	($dataSets->{'nodeIsIsolated'                     } == 1)
	&
	($dataSets->{'mergingStatisticsNodeHierarchyLevel'} == 0)
    );
# Output to file the masses and weights.
my $treeMassFile = new PDL::IO::HDF5(">".$modelDirectory."/treeMasses.hdf5");
$treeMassFile->dataset('treeRootMass')->set($dataSets->{'basicMass'       }->($isolated));
$treeMassFile->dataset('treeWeight'  )->set($dataSets->{'mergerTreeWeight'}->($isolated));

exit;
