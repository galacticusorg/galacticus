#!/usr/bin/env perl
use strict;
use warnings;
use PDL;
use PDL::IO::HDF5;

# Run a test case which checks that satelliteBoundMass is equal to basicMass in the preset satellite class for isolated halos and
# when no bound mass history is set.
# Andrew Benson (20-September-2016)

# Run the model and check for successful completion.
system("./Galacticus.exe testSuite/regressions/satellitePresetBoundMassNonZero.xml");
die("FAILED: satellitePresetBoundMassNonZero.pl model failed to complete") 
    unless ( $? == 0 );

# Extract required information.
my $model              = new PDL::IO::HDF5("testSuite/outputs/regressions/satellitePresetBoundMassNonZero.hdf5");
my $nodeData           = $model   ->group  ("Outputs"           )->group("Output1")->group("nodeData");
my $satelliteBoundMass = $nodeData->dataset("satelliteBoundMass")->get();
my $basicMass          = $nodeData->dataset("basicMass"         )->get();

# Check for consistency.
die("FAILED: satellitePresetBoundMassNonZero.pl: satelliteBoundMass does not always equal basicMass")
    unless ( all($satelliteBoundMass == $basicMass) );

exit;
