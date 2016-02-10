#!/usr/bin/env perl
use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use Data::Dumper;
$PDL::BIGPDL = 1;

# Extract a single merger tree from a merger tree file. Useful for creation of test cases.
# Andrew Benson (19-September-2012)

# Get arguments.
die("Usage: extractSingleTree.pl <fromFile> <toFile> <forestIndex>")
    unless ( scalar(@ARGV) == 3 );
my $inFileName  = $ARGV[0];
my $outFileName = $ARGV[1];
my $tree        = $ARGV[2];

# Remove the output file.
unlink($outFileName);

# Open the files.
my $inFile  = new PDL::IO::HDF5(     $inFileName);
my $outFile = new PDL::IO::HDF5(">".$outFileName);

# Detect file version and set group names appropriately.
(my $inFormatVersion) = (grep {$_ eq "formatVersion"} $inFile->attrs()) ? $inFile->attrGet("formatVersion") : (1);
my $forestName      = $inFormatVersion == 1 ? "treeIndex" : "forestIndex";
my $halosName       = $inFormatVersion == 1 ? "haloTrees" : "forestHalos";

# Find the tree to extract.
my $forestIndex = $inFile->group($forestName)->dataset($forestName    )->get();
my $firstNode   = $inFile->group($forestName)->dataset("firstNode"    )->get();
my $nodeCount   = $inFile->group($forestName)->dataset("numberOfNodes")->get();
my $selected    = which($forestIndex == $tree);
my $start       = $firstNode->index($selected)->sclr();
my $count       = $nodeCount->index($selected)->sclr();
my $end         = $start+$count-1;

# Read all halo datasets.
foreach my $datasetName ( $inFile->group($halosName)->datasets() ) {
    my $dataset    = $inFile->group($halosName)->dataset($datasetName)->get();
    next
	if ( $datasetName =~ m/^particleIndex(Start|Count)$/ );
    my $dataset    = $inFile->group($halosName)->dataset($datasetName)->get();
    my $dimensions = $dataset->ndims();
     if      ( $dimensions == 1 ) {
     	$outFile->group("forestHalos")->dataset($datasetName)->set($dataset->(  $start:$end));
     } elsif ( $dimensions == 2 ) {
     	$outFile->group("forestHalos")->dataset($datasetName)->set($dataset->(:,$start:$end));
     } else {
     	die("extractSingleTree.pl: unable to handle rank > 2 datasets");
     }
}

# Create the forestIndex group.
$outFile->group("forestIndex")->dataset("forestIndex"  )->set(pdl longlong([$tree ]));
$outFile->group("forestIndex")->dataset("firstNode"    )->set(pdl longlong([     0]));
$outFile->group("forestIndex")->dataset("numberOfNodes")->set(pdl longlong([$count]));

# Set format version.
my $formatVersion;
if ( grep {$_ eq "formatVersion"} $outFile->attrs() ) {
    $formatVersion = $inFile->attrGet("formatVersion");
} else {
    $formatVersion = pdl long 2;
}
$outFile->attrSet("formatVersion" => $formatVersion);

# Copy all attributes.
foreach my $groupName ( $inFile->groups() ) {
    my $outGroupName = $groupName;
    $outGroupName = "forestIndex"
	if ( $groupName eq "treeIndex" );
    $outGroupName = "forestHalos"
	if ( $groupName eq "haloTrees" );
    next 
	if ( $groupName eq "particles" );
    my $fromGroup = $inFile ->group(   $groupName);
    my $toGroup   = $outFile->group($outGroupName);
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
