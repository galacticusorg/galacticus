#!/usr/bin/env perl
use strict;
use warnings;
use PDL;
use PDL::IO::HDF5;

# Run a single Galacticus model and ensure that outputs are created.
# Andrew Benson (04-April-2014)

# Simply run the models.
system("cd ..; mkdir -p testSuite/outputs; ./Galacticus.exe testSuite/parameters/test-output.xml");

# Check for outputs.
die("test-output.pl: FAILED to run Galacticus model")
    unless ( -e "outputs/test-output.hdf5" );
my $file = new PDL::IO::HDF5("outputs/test-output.hdf5");
my @groups1 = $file->groups();
die("test-output.pl FAIL - no Outputs group exists")
    unless ( grep {$_ eq "Outputs"} @groups1 );
my @groups2 = $file->group("Outputs")->groups();
die("test-output.pl FAIL - no Output1 group exists")
    unless ( grep {$_ eq "Output1"} @groups2 );
my @groups3 = $file->group("Outputs")->group("Output1")->groups();
die("test-output.pl FAIL - no nodeData group exists")
    unless ( grep {$_ eq "nodeData"} @groups3 );
my @datasets = $file->group("Outputs")->group("Output1")->group("nodeData")->datasets();
die("test-output.pl FAIL - no nodeIndex dataset exists")
    unless ( grep {$_ eq "nodeIndex"} @datasets );
print "test-output.pl: SUCCESS\n";

exit;
