#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use Galacticus::Validation;

# Run models to validate a dark matter only subhalo evolution model.
# Andrew Benson (05-August-2022)

# Make output directory.
system("mkdir -p outputs/");

# Run the validate model.
system("cd ..; ./Galacticus.exe testSuite/parameters/validate_darkMatterOnlySubHalos.xml");
unless ( $? == 0 ) {
    print "FAIL: dark matter-only subhalos validation model failed to run\n";
    exit;
}

# Extract and validate the likelihoods.
&Galacticus::Validation::extract("outputs/validate_darkMatterOnlySubHalos.hdf5","Dark Matter Only Subhalos","darkMatterOnlySubhalos");

print "SUCCESS: dark matter-only subhalos validation model\n";

exit;
