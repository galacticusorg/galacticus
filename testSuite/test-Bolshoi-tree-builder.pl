#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use File::Which;

# Take trees from the Bolshoi simulation and process into both Galacticus and IRATE formats.
# Andrew Benson (12-October-2012)
# Contributions to this file made by: Stephanie Dörschner.

# Create output folder.
system("cd ..; mkdir -p testSuite/outputs");

# Create a file in Galacticus format.
system("cd ..; ./Galacticus.exe testSuite/parameters/mergerTreeFileBuildBolshoi.xml");
die("FAILED: failed to make Galacticus-format merger tree file from Bolshoi merger tree")
    unless ( $? == 0 );

# Create a file in IRATE format.
system("cd ..; ./Galacticus.exe testSuite/parameters/mergerTreeFileBuildBolshoiIRATE.xml");
die("FAILED: failed to make IRATE-format merger tree file from Bolshoi merger tree")
    unless ( $? == 0 );

# Validate the IRATE format file.
my $validator = &File::Which::which('iratevalidate');
if ( $validator ) {
    system("cd ..; iratevalidate testSuite/outputs/bolshoiTestTreesIRATE_in.hdf5");
    die("FAILED: IRATE-format file generated from Bolshoi merger tree did not validate")
	unless ( $? == 0 );
} else {
    print "SKIP: iratevalidate is not installed - validation of IRATE-format file will be skipped\n";
}

# Run Galacticus on file in Galacticus format.
system("cd ..; ./Galacticus.exe testSuite/parameters/bolshoiTestTreesGLC.xml");

exit;
