#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V094"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V094"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
}
unshift(@INC,$galacticusPath."perl");
require Galacticus::Build::SourceTree;

# Preprocess a Galacticus Fortran source file.
# Andrew Benson (17-April-2015)

# Get arguments.
die("Usage: preprocess.pl <infile> <outfile>")
    unless ( scalar(@ARGV) == 2 );
my $inputFileName  = $ARGV[0];
my $outputFileName = $ARGV[1];

# Parse the file to build a tree.
my $tree = &SourceTree::ParseFile($inputFileName);

# Process the tree.
&SourceTree::ProcessTree($tree);

# Serialize back to source code.
open(my $outputFile,">",$outputFileName);
print $outputFile &SourceTree::Serialize($tree);
close($outputFile);

exit;
