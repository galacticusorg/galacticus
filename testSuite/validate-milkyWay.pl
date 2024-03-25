#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use Galacticus::Validation;

# Run models to validate a Milky Way model.
# Andrew Benson (10-August-2022)

# Make output directory.
system("mkdir -p outputs/");

# Run the validate model.
system("cd ..; ./Galacticus.exe testSuite/parameters/validate_milkyWay.xml");
unless ( $? == 0 ) {
    print "FAIL: Milky Way validation model failed to run\n";
    exit;
}

# Extract and validate the likelihoods.
&Galacticus::Validation::extract("outputs/validate_milkyWay.hdf5","Milky Way model","milkyWayModel","testSuite/parameters/validate_milkyWay.xml");

print "SUCCESS: Milky Way validation model\n";

exit;
