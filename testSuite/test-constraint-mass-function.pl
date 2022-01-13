#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;

# Test that the dark matter halo mass function matches the constraint.
# Andrew Benson (15-January-2020)

# Run model.
system("cd ..; mkdir -p testSuite/outputs; ./Galacticus.exe testSuite/parameters/constrainHaloMassFunction.xml");
unless ( $? == 0 ) {
    print "FAILED: Galacticus model failed to run\n";
    exit 0;
}

# Read model results.
my $model             = new PDL::IO::HDF5('outputs/constrainHaloMassFunction.hdf5');
my $outputs           = $model  ->group  ('Outputs'                       )       ;
my $output            = $outputs->group  ('Output1'                       )       ;
my $massHaloModel     = $output ->dataset('haloMass'                      )->get();
my $massFunctionModel = $output ->dataset('haloMassFunctionLnMBinAveraged')->get();

# Read target data.
my $target                  = new PDL::IO::HDF5('data/darkMatterHaloMassFunctionMDPL2.hdf5');
my $massHaloTarget          = $target->dataset('massHalo'         )->get();
my $massFunctionTarget      = $target->dataset('massFunction'     )->get();
my $massFunctionErrorTarget = $target->dataset('massFunctionError')->get();
my $emptyBins               = which($massFunctionTarget <= 0.0);
my $haloCount               = 1.0/($massFunctionErrorTarget/$massFunctionTarget)**2;
$haloCount->($emptyBins)   .= 0.0;

# Check that halo masses agree.
my $massDifferenceFractional = abs($massHaloModel-$massHaloTarget)/$massHaloTarget;
unless ( all($massDifferenceFractional < 1.0e-6) ) {
    print "FAILED: halo masses differ from target dataset\n";
    exit 0;
}

# Specify ranges.
my $massRangeMinimum = pdl  1.0e12;
my $haloCountMinimum =     30     ;

# Select bins within range.
my $inRange = which(($massHaloTarget > $massRangeMinimum) & ($haloCount >= $haloCountMinimum));

# Compute normalized offsets.
my $offsets = abs($massFunctionModel->($inRange)-$massFunctionTarget->($inRange))/$massFunctionErrorTarget->($inRange);

# Check offsets are sufficiently small.
if ( all($offsets < 4.0) ) {
    print "SUCCESS: halo mass function matches target\n";
} else {
    print "FAILED: halo mass function does not match target\n";
}

exit 0;
