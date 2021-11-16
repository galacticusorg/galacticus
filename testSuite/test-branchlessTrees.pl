#!/usr/bin/env perl
use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;

# Run models that test that the branchless merger tree algorithm.
# Andrew Benson (15-July-2021)

# Make output directory.
system("mkdir -p outputs/");

# Run the branchless model.
system("cd ..; ./Galacticus.exe testSuite/parameters/mergerTreeBranchless.xml");
unless ( $? == 0 ) {
    print "FAIL: merger tree branchless model failed to run\n";
    exit;
}

my $success   = 1;
my $model     = new PDL::IO::HDF5("outputs/mergerTreeBranchless.hdf5");
my $structure = $model->group('mergerTreeStructures');
foreach my $treeName ( $structure->groups() ) {
    my $tree = $structure->group($treeName);
    my $properties;
    $properties->{$_} = $tree->dataset($_)->get()
	foreach ( 'nodeIndex', 'parentIndex', 'nodeIsOnMainBranch' );
    # If this tree is truly branchless then any node which is not on the main branch should have its parent in the subsequent
    # entry in the arrays, and that parent should be on the main branch.
    my $sideBranch = which($properties->{'nodeIsOnMainBranch'} == 0);
    my $mainBranch = $sideBranch+1;
    $success = 0
	unless ( all($properties->{'parentIndex'}->($sideBranch) == $properties->{'nodeIndex'}->($mainBranch)) && all($properties->{'nodeIsOnMainBranch'}->($mainBranch) == 1) );
}

my $status = $success ? "SUCCESS" : "FAILED";
print $status.": branchless merger tree construction\n";

exit;
