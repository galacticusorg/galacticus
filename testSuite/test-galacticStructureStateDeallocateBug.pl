#!/usr/bin/env perl
use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;

# Run models that test for regression of the deallocate bug in galacticStructureState
# Andrew Benson (02-May-2022)

# Make output directory.
system("mkdir -p outputs/");

# Run the model.
system("cd ..; ./Galacticus.exe testSuite/parameters/galacticStructureStateDeallocateBug.xml");
if ( $? == 0 ) {
    print "SUCCESS:  galacticStructureState deallocate bug\n";
} else {
    print "FAIL: galacticStructureState deallocate bug model failed to run\n";
}

exit;
