#!/usr/bin/env perl
use strict;
use warnings;
use PDL;
use PDL::IO::HDF5;

# Run a test case which checks that volume weights of merger trees read from file are assigned correctly when they are computed
# from the simulation box size.
# Andrew Benson (16-January-2014)

# Run the model and check for successful completion.
system("./Galacticus.exe testSuite/regressions/mergerAtFinalTimeInTree.xml");
die("FAILED: mergerTreeBoxSizeWeight.pl model failed to complete") 
    unless ( $? == 0 );

# Extract required information from files.
my $trees                       = new PDL::IO::HDF5("testSuite/data/mergerTrees/mergerAtFinalTimeInTree.hdf5"   );
my $model                       = new PDL::IO::HDF5("testSuite/outputs/regressions/mergerAtFinalTimeInTree.hdf5");
(
 my $boxSize
)                               = $trees->group("simulation")
    ->attrGet(
    'boxSize'
    );
(
 my $lengthHubbleExponent     ,
 my $lengthScaleFactorExponent,
 my $lengthUnitsInSI
)                               = $trees->group("units"     )
    ->attrGet(
    'lengthHubbleExponent'     ,
    'lengthScaleFactorExponent',
    'lengthUnitsInSI'
    );
(
 my $hubble
)                               = $trees->group("cosmology" )
    ->attrGet(
    'HubbleParam'
    );
my $treeWeight                  = $model->group("Outputs"   )->group("Output1")->dataset("mergerTreeWeight")->get();

# Compute expected tree weight.
my $scaleFactor    = pdl 1.0;
my $megaParsec     = pdl 3.08567758135e+22;
$boxSize           = $boxSize*($hubble**$lengthHubbleExponent)*($scaleFactor**$lengthScaleFactorExponent)*($lengthUnitsInSI/$megaParsec);
my $weightExpected = 1.0/$boxSize**3;

# Compute fractional difference.
my $fractionalError = abs($treeWeight-$weightExpected)/$weightExpected;

# Check for consistency.
die("FAILED: mergerTreeBoxSizeWeight.pl tree weights do not equal expected values")
    unless ( all($fractionalError < 1.0e-6) );

exit;
