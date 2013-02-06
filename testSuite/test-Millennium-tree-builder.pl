#!/usr/bin/env perl
use strict;
use warnings;
use File::Which;

# Take trees from the Millennium Database and process into both Galacticus and IRATE formats.
# Andrew Benson (12-October-2012)

# Build the file maker code.
system("cd ..; make Millennium_Merger_Tree_File_Maker.exe");

# Create a file in Galacticus format.
system("cd ..; Millennium_Merger_Tree_File_Maker.exe testSuite/data/mergerTrees/millenniumTestTrees.csv testSuite/data/mergerTrees/millenniumTestTreesParticles.csv testSuite/outputs/millenniumTestTreesGLC.hdf5 galacticus 1");
die("FAILED: failed to make Galacticus-format merger tree file from Millennium database output")
    unless ( $? == 0 );

# Create a file in IRATE format.
system("cd ..; Millennium_Merger_Tree_File_Maker.exe testSuite/data/mergerTrees/millenniumTestTrees.csv testSuite/data/mergerTrees/millenniumTestTreesParticles.csv testSuite/outputs/millenniumTestTreesIRATE.hdf5 irate 1");
die("FAILED: failed to make IRATE-format merger tree file from Millennium database output")
    unless ( $? == 0 );

# Validate the IRATE format file.
my $validator = which('iratevalidate');
if ( $validator ) {
    system("cd ..; iratevalidate testSuite/outputs/millenniumTestTreesIRATE.hdf5");
    die("FAILED: IRATE-format file generated from Millennium database output did not validate")
	unless ( $? == 0 );
} else {
    print "SKIP: iratevalidate is not installed - validation of IRATE-format file will be skipped\n";
}

exit;
