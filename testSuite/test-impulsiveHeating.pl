#!/usr/bin/env perl
use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use Data::Dumper;

# Run models that test that the impulsive heating model works as expected.
# Andrew Benson (02-June-2022)

# Make output directory.
system("mkdir -p outputs/");

# Run the impulsive heating model.
system("export OMP_NUM_THREADS=1; cd ..; ./Galacticus.exe testSuite/parameters/impulsiveHeating.xml");
if ( $? == 0 ) {
    print "success: impulsive heating model ran successfully\n";
} else {
    print "FAIL: impulsive heating model failed to run\n";
}

exit;
