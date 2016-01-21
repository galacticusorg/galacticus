#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V093"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V093"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
 $ENV{"GALACTICUS_ROOT_V093"} = "./";
}
unshift(@INC,$galacticusPath."perl"); 
require File::Which;

# Take trees from the Bolshoi simulation and process into both Galacticus and IRATE formats.
# Andrew Benson (12-October-2012)
# Contributions to this file made by: Stephanie Dörschner.

# Build the file maker code.
system("cd ..; make Bolshoi_Merger_Tree_File_Maker.exe Galacticus.exe");

# Create output folder.
system("cd ..; mkdir -p testSuite/outputs");

# Create a file in Galacticus format.
system("cd ..; ./Bolshoi_Merger_Tree_File_Maker.exe testSuite/data/mergerTrees/bolshoiTestTrees.dat testSuite/outputs/bolshoiTestTreesGLC_in.hdf5 galacticus 1");
die("FAILED: failed to make Galacticus-format merger tree file from Bolshoi merger tree")
    unless ( $? == 0 );

# Create a file in IRATE format.
system("cd ..; ./Bolshoi_Merger_Tree_File_Maker.exe testSuite/data/mergerTrees/bolshoiTestTrees.dat testSuite/outputs/bolshoiTestTreesIRATE_in.hdf5 irate 1");
die("FAILED: failed to make IRATE-format merger tree file from Bolshoi merger tree")
    unless ( $? == 0 );

# Validate the IRATE format file.
my $validator = which('iratevalidate');
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
