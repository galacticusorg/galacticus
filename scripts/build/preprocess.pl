#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use Galacticus::Build::SourceTree;

# Preprocess a Galacticus Fortran source file.
# Andrew Benson (17-April-2015)

# Get arguments.
die("Usage: preprocess.pl <infile> <outfile>")
    unless ( scalar(@ARGV) == 2 );
my $inputFileName  = $ARGV[0];
my $outputFileName = $ARGV[1];

# Parse the file to build a tree.
my $tree = &Galacticus::Build::SourceTree::ParseFile($inputFileName);

# Process the tree.
&Galacticus::Build::SourceTree::ProcessTree($tree);

# Serialize back to source code.
open(my $outputFile,">",$outputFileName);
print $outputFile &Galacticus::Build::SourceTree::Serialize($tree);
close($outputFile);

exit;
