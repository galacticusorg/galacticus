#!/usr/bin/env perl
use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;

# Run models that test that merger trees are output in contiguous blocks.
# Andrew Benson (14-July-2021)

# Make output directory.
system("mkdir -p outputs/");

# Run the model.
system("cd ..; ./Galacticus.exe testSuite/parameters/outputTreeContiguosity.xml");
unless ( $? == 0 ) {
    print "FAIL: output tree contiguosity model failed to run\n";
    exit;
}

# Check for contiguous output.
my $model   = new PDL::IO::HDF5("outputs/outputTreeContiguosity.hdf5");
my $outputs = $model->group('Outputs');
my $success = 1;
foreach my $outputName ( $outputs->groups() ) {
    my $output          = $outputs ->group  ($outputName           )       ;
    my $nodeData        = $output  ->group  ('nodeData'            )       ;
    my $treeIndex       = $output  ->dataset('mergerTreeIndex'     )->get();
    my $treeStart       = $output  ->dataset('mergerTreeStartIndex')->get();
    my $treeCount       = $output  ->dataset('mergerTreeCount'     )->get();
    my $treeIndexByNode = $nodeData->dataset('mergerTreeIndex'     )->get();
    for(my $i=0;$i<nelem($treeIndex);++$i) {
	my $indexTarget =               $treeIndex->(($i));
	my $indexStart  =               $treeStart->(($i));
	my $indexEnd    = $indexStart-1+$treeCount->(($i));
	$success = 0
	    unless ( all($treeIndexByNode->($indexStart:$indexEnd) == $indexTarget) );
    }
}
my $status = $success ? "SUCCESS" : "FAILED";
print $status.": merger tree output contiguosity\n";
exit;
