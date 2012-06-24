#!/usr/bin/env perl
use strict;
use warnings;
use PDL;
use PDL::IO::HDF5;
use PDL::NiceSlice;

# Run a set of Galacticus models to test the state store/retrieve functionality.
# Andrew Benson (23-Jun-2012)

# Run full model and a restored state model.
system("export OMP_NUM_THREADS=1; Galacticus.exe testSuite/parameters/state/store.xml"   );
system("export OMP_NUM_THREADS=1; Galacticus.exe testSuite/parameters/state/retrieve.xml");

# Open both output files.
my $store    = new PDL::IO::HDF5("testSuite/outputs/stateStore.hdf5"   );
my $retrieve = new PDL::IO::HDF5("testSuite/outputs/stateRetrieve.hdf5");

# Get data groups.
my $storeData    = $store   ->group('Outputs')->group('Output1')->group('nodeData');
my $retrieveData = $retrieve->group('Outputs')->group('Output1')->group('nodeData');

# Get number of nodes in final tree.
my $storeTreeSize    = $store   ->group('Outputs')->group('Output1')->dataset('mergerTreeCount')->get()->((-1));
my $retrieveTreeSize = $retrieve->group('Outputs')->group('Output1')->dataset('mergerTreeCount')->get()->((-1));

# Check that the number of nodes is the same.
unless ( $storeTreeSize == $retrieveTreeSize ) {
    print "FAIL: number of nodes in output changed after state retrieve\n";
    exit;
}

# Get all available datasets.
my @datasets = $storeData->datasets();

# Check that each dataset is unchanged.
foreach my $dataset ( @datasets ) {
    my $storeDataset    = $storeData   ->dataset($dataset)->get()->(-$storeTreeSize:-1);
    my $retrieveDataset = $retrieveData->dataset($dataset)->get();
    my $equal = all($storeDataset == $retrieveDataset);
    print "FAIL: dataset '".$dataset."' chanegd after state retrieve\n"
	unless ( $equal == 1 );
}

exit;
