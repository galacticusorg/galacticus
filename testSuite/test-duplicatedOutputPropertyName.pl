#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use File::Slurp qw(slurp);
use System::Redirect;

# Check that duplicated output property names are caught.
# Andrew Benson (28-February-2024)

# Run the model.
&System::Redirect::tofile("cd ..; export OMP_NUM_THREADS=1; ./Galacticus.exe testSuite/parameters/duplicatedOutputPropertyName.xml","outputs/duplicatedOutputPropertyName.log");
if ( $? == 0 ) {
    print "FAILED: duplicated parameter name was not detected\n";
} else {
    system("grep -q 'duplicate property name' outputs/duplicatedOutputPropertyName.log");
    if ( $? == 0 ) {
	print "SUCCESS: duplicated parameter name was detected\n";
    } else {
	print "FAILED: duplicated parameter name error message was not given\n";
    }
}

exit 0;
