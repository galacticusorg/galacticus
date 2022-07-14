#!/usr/bin/env perl
use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;

# Run models that test that the bar instability mechanism is able to create spheroids.
# Andrew Benson (25-July-2018)

# Make output directory.
system("mkdir -p outputs/");

# Run the model.
system("cd ..; ./Galacticus.exe testSuite/parameters/barInstability.xml");
if ( $? == 0 ) {
    my $model               = new PDL::IO::HDF5("outputs/barInstability.hdf5");
    my $outputs             = $model  ->group  ('Outputs'            )       ;
    my $output              = $outputs->group  ('Output1'            )       ;
    my $nodes               = $output ->group  ('nodeData'           )       ;
    my $spheroidMassStellar = $nodes  ->dataset('spheroidMassStellar')->get();
    if ( $spheroidMassStellar->((0)) <= 0.0 ) {
	print "FAIL: spheroid has non-positive mass\n";
    } else {
	print "SUCCESS\n";
    }
} else {
    print "FAIL: bar instability test model failed to run\n";
}

exit;
