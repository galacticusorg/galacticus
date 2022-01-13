#!/usr/bin/env perl
use strict;
use warnings;
use PDL;
use PDL::IO::HDF5;

# Run a test case for bug #1066052: "Crash due to zero-size spheroids".
# Andrew Benson (12-October-2012)

# Run the model and check for successful completion.
system("./Galacticus.exe testSuite/parameters/bug1066052.xml");
die("FAILED: bug1066052.xml model failed to complete") 
    unless ( $? == 0 );

# Check that the spheroid size is non-zero whenever the spheroid mass is non-zero.
my $model       = new PDL::IO::HDF5("testSuite/outputs/bug1066052.hdf5");
my $nodeData    = $model->group("Outputs")->group("Output1")->group("nodeData");
my $scaleLength = $nodeData->dataset("spheroidRadius"     )->get();
my $stellarMass = $nodeData->dataset("spheroidMassStellar")->get();
my $gasMass     = $nodeData->dataset("spheroidMassGas"    )->get();
my $mass        = $stellarMass+$gasMass;
my $badNodes    = which(($mass > 0.0) & ($scaleLength <= 0.0));
die("FAILED: bug1066052.xml model contains zero-sized spheroids")
    if ( nelem($badNodes) > 0 );

exit;
