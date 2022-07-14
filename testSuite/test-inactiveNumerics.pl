#!/usr/bin/env perl
use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;

# Run models that test that inactive properties are evolved when using an ODE solver that does not support inactive property evaluation.
# Andrew Benson (04-May-2022)

# Make output directory.
system("mkdir -p outputs/");

# Run the model.
system("cd ..; ./Galacticus.exe testSuite/parameters/inactiveNumerics.xml");
unless ( $? == 0 ) {
    print "FAIL: inactiveNumerics model failed to run\n";
}

# Check that luminosities are non-zero.
my $nonZeroLuminosities = 0;
my $model               = new PDL::IO::HDF5("outputs/inactiveNumerics.hdf5");
my $nodes               = $model->group('Outputs/Output1/nodeData');
foreach my $datasetName ( $nodes->datasets() ) {
    next
	unless ( $datasetName =~ m/Luminosities/ );
    my $dataset = $nodes->dataset($datasetName)->get();
    $nonZeroLuminosities = 1
	if ( any($dataset > 0.0) );
}
if ( $nonZeroLuminosities ) {
    print "SUCCESS: non-zero luminosities are present\n";
} else {
    print "FAIL: all luminosities are zero\n";
}

exit;
