#!/usr/bin/env perl
use strict;
use warnings;
use PDL;
use PDL::IO::HDF5;
use PDL::NiceSlice;

# Run a set of Galacticus models to test the state store/retrieve functionality under MPI.
# Andrew Benson (15-Jun-2018)

# Run full store model.
system("export OMP_NUM_THREADS=1; rm -f outputs/state.state* outputs/state.fgsl.state*; cd ..; mpirun -np 8 Galacticus.exe testSuite/parameters/state/store.xml"  );
die("FAILED: failed to run store model")
    unless ( $? == 0 );
# Find which threads ran the final tree.
my $finalTreeThreadOpenMP;
my $finalTreeProcessMPI;
opendir(my $stateDirectory,"outputs");
while ( my $fileName = readdir($stateDirectory) ) {
    if  ( $fileName =~ m/state\.state\.log\:openMP(\d+):MPI(\d+)/ ) {
	my $threadOpenMP = $1;
	my $processMPI   = $2;
	open(my $stateLogFile,"outputs/".$fileName);
	while (my $line = <$stateLogFile> ) {
	    if ( $line =~ m/^\s*Storing state for tree #(\d+)/ ) {
		if ( $1 == 15 ) {
		    $finalTreeThreadOpenMP = $threadOpenMP;
		    $finalTreeProcessMPI   = $processMPI  ;
		}
	    }
	}
	close($stateLogFile);
    }    
}
closedir($stateDirectory);
if ( defined($finalTreeThreadOpenMP) && defined($finalTreeProcessMPI) ) {
    print "Final tree was run by OpenMP thread ".$finalTreeThreadOpenMP." on MPI process ".$finalTreeProcessMPI."\n";
    system("cp -f outputs/state.state:openMP"     .$finalTreeThreadOpenMP.":MPI".$finalTreeProcessMPI." outputs/state.state:MPI0000"     );
    system("cp -f outputs/state.fgsl.state:openMP".$finalTreeThreadOpenMP.":MPI".$finalTreeProcessMPI." outputs/state.fgsl.state:MPI0000");
} else {
    die("FAILED: failed to identify which thread/process ran final tree");
}

# Run the restore model.
system("export OMP_NUM_THREADS=1; cd ..; mpirun -np 1 Galacticus.exe testSuite/parameters/state/retrieve.xml");
die("FAILED: failed to run retrieve model")
    unless ( $? == 0 );

# Open both output files.
die("FAILED: stateStore:MPI".$finalTreeProcessMPI.".hdf5 file is missing")
    unless ( -e "outputs/stateStore:MPI".$finalTreeProcessMPI.".hdf5" );
die("FAILED: stateRetrieve:MPI0000.hdf5 file is missing")
    unless ( -e "outputs/stateRetrieve:MPI0000.hdf5" );
my $store    = new PDL::IO::HDF5("outputs/stateStore:MPI".$finalTreeProcessMPI.".hdf5");
my $retrieve = new PDL::IO::HDF5("outputs/stateRetrieve:MPI0000.hdf5"                 );

# Get data groups.
my $storeData    = $store   ->group('Outputs')->group('Output1')->group('nodeData');
my $retrieveData = $retrieve->group('Outputs')->group('Output1')->group('nodeData');

# Find the tree in the store model.
my $storeTreeIndex   = $store->group('Outputs')->group('Output1')->dataset('mergerTreeIndex')->get();
my $treeFinal        = which($storeTreeIndex == 15);
unless ( nelem($treeFinal) == 1 ) {
    print "FAILED: unable to (uniquely) identify final tree in stored model output\n";	
    exit;
}
my $treeFinalIndex = $treeFinal->((0))->sclr();

# Get number of nodes in final tree.
my $storeTreeStart   = $store   ->group('Outputs')->group('Output1')->dataset('mergerTreeStartIndex')->get()->(($treeFinalIndex));
my $storeTreeSize    = $store   ->group('Outputs')->group('Output1')->dataset('mergerTreeCount'     )->get()->(($treeFinalIndex));
my $retrieveTreeSize = $retrieve->group('Outputs')->group('Output1')->dataset('mergerTreeCount'     )->get()->((             -1));

# Check that the number of nodes is the same.
unless ( $storeTreeSize == $retrieveTreeSize ) {
    print "FAILED: number of nodes in output changed after state retrieve\n";
    exit;
}

# Get all available datasets.
my @datasets = $storeData->datasets();

# Check that each dataset is unchanged.
my $failed = 0;
foreach my $dataset ( @datasets ) {
    my $storeDataset    = $storeData   ->dataset($dataset)->get()->($storeTreeStart:$storeTreeStart+$storeTreeSize-1);
    my $retrieveDataset = $retrieveData->dataset($dataset)->get();
    my $equal = all($storeDataset == $retrieveDataset);
    unless ( $equal == 1 ) {
	print "FAILED: dataset '".$dataset."' changed after state retrieve\n";
	$failed = 1;
    } else {
 	print "SUCCESS: dataset '".$dataset."'\n";
   }
}
print "SUCCESS!\n"
    unless ( $failed );
exit;
