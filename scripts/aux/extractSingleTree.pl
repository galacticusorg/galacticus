#!/usr/bin/env perl
use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use Data::Dumper;

# Extract a single merger tree from a merger tree file. Useful for creation of test cases.
# Andrew Benson (19-September-2012)

# Get arguments.
die("Usage: extractSingleTree.pl <fromFile> <toFile> <treeIndex>")
    unless ( scalar(@ARGV) == 3 );
my $inFileName  = $ARGV[0];
my $outFileName = $ARGV[1];
my $tree        = $ARGV[2];

# Remove the output file.
unlink($outFileName);

# Open the files.
my $inFile  = new PDL::IO::HDF5(     $inFileName);
my $outFile = new PDL::IO::HDF5(">".$outFileName);

# Find the tree to extract.
my $treeIndex = $inFile->group("treeIndex")->dataset("treeIndex"    )->get();
my $firstNode = $inFile->group("treeIndex")->dataset("firstNode"    )->get();
my $nodeCount = $inFile->group("treeIndex")->dataset("numberOfNodes")->get();
my $selected  = which($treeIndex == $tree);
my $start     = $firstNode->index($selected)->sclr();
my $count     = $nodeCount->index($selected)->sclr();
my $end       = $start+$count-1;

# Read all haloTrees datasets.
foreach my $datasetName ( $inFile->group("haloTrees")->datasets() ) {
    my $dataset    = $inFile->group("haloTrees")->dataset($datasetName)->get();
    my $dimensions = $dataset->ndims();
     if      ( $dimensions == 1 ) {
     	$outFile->group("haloTrees")->dataset($datasetName)->set($dataset->(  $start:$end));
     } elsif ( $dimensions == 2 ) {
     	$outFile->group("haloTrees")->dataset($datasetName)->set($dataset->(:,$start:$end));
     } else {
     	die("??");
     }
}

# Create the treeIndex group.
$outFile->group("treeIndex")->dataset("treeIndex"    )->set(pdl longlong([$tree ]));
$outFile->group("treeIndex")->dataset("firstNode"    )->set(pdl longlong([     0]));
$outFile->group("treeIndex")->dataset("numberOfNodes")->set(pdl longlong([$count]));

# Copy all attributes.
foreach my $groupName ( $inFile->groups() ) {
    my $fromGroup = $inFile ->group($groupName);
    my $toGroup   = $outFile->group($groupName);
    &Copy_Attributes($fromGroup,$toGroup);
}

exit;

sub Copy_Attributes {
    # Copy all attributes from an object in the input files to the output file.
    # Get the objects from and to.
    my $objectFrom = shift;
    my $objectTo   = shift;
  
    # Copy all attributes.
    my @attributes = $objectFrom->attrs();
    foreach my $attribute ( @attributes ) {
	my @attrValue = $objectFrom->attrGet($attribute);
	$objectTo->attrSet($attribute => $attrValue[0]);
    }
}
