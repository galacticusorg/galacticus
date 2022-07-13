#!/usr/bin/env perl
use strict;
use warnings;
use PDL;
use PDL::IO::HDF5;
use PDL::NiceSlice;

# Run a set of Galacticus models to test the state store/retrieve functionality.
# Andrew Benson (23-Jun-2012)

# Allow two run-throughs. If tabulation files are updated during the "store" run it can result in divergent behavior in the
# "restore" run. This is a limitation of our current state store/restore infrastructure. So, if state store/restore tests fail on
# the first run through we allow a second attempt - for which tabulations should already be sufficient so we would expect this
# issue to not occur.
for(my $pass=0;$pass<2;++$pass) {
    # Report on pass.
    print "State store/restore test pass #".$pass."\n";

    # Run full store model.
    system("export OMP_NUM_THREADS=12; rm -f outputs/state.state*:openMP* outputs/state.gsl.state*:openMP*; cd ..; ./Galacticus.exe testSuite/parameters/state/store.xml"  );
    die("FAILED: failed to run store model")
	unless ( $? == 0 );
    # Find which threads ran the final tree.
    my $finalTreeThread;
    opendir(my $stateDirectory,"outputs");
    while ( my $fileName = readdir($stateDirectory) ) {
	if  ( $fileName =~ m/state\.state\.log:openMP(\d+)/ ) {
	    my $thread = $1;
	    open(my $stateLogFile,"outputs/".$fileName);
	    while (my $line = <$stateLogFile> ) {
		if ( $line =~ m/^\s*Storing state for tree #(\d+)/ ) {
		    $finalTreeThread = $thread
			if ( $1 == 15 );
		}
	    }
	    close($stateLogFile);
	}    
    }
    closedir($stateDirectory);
    if ( defined($finalTreeThread) ) {
	print "Final tree was run by thread ".$finalTreeThread."\n";
	system("cp -f outputs/state.state:openMP"    .$finalTreeThread." outputs/state.state"    );
	system("cp -f outputs/state.gsl.state:openMP".$finalTreeThread." outputs/state.gsl.state");
    } else {
	die("FAILED: failed to identify which thread ran final tree");
    }

    # Run the restore model.
    system("export OMP_NUM_THREADS=1; cd ..; ./Galacticus.exe testSuite/parameters/state/retrieve.xml");
    die("FAILED: failed to run retrieve model")
	unless ( $? == 0 );

    # Open both output files.
    die("FAILED: stateStore.hdf5 file is missing")
	unless ( -e "outputs/stateRetrieve.hdf5" );
    die("FAILED: stateRetrieve.hdf5 file is missing")
	unless ( -e "outputs/stateStore.hdf5" );
    my $store    = new PDL::IO::HDF5("outputs/stateStore.hdf5"   );
    my $retrieve = new PDL::IO::HDF5("outputs/stateRetrieve.hdf5");

    # Get data groups.
    my $storeData    = $store   ->group('Outputs')->group('Output1')->group('nodeData');
    my $retrieveData = $retrieve->group('Outputs')->group('Output1')->group('nodeData');

    # Find the tree in the store model.
    my $storeTreeIndex   = $store->group('Outputs')->group('Output1')->dataset('mergerTreeIndex')->get();
    my $treeFinal        = which($storeTreeIndex == 15);
    unless ( nelem($treeFinal) == 1 ) {
	die("FAILED: unable to (uniquely) identify final tree in stored model output\n");
    }
    my $treeFinalIndex = $treeFinal->((0))->sclr();

    # Get number of nodes in final tree.
    my $storeTreeStart   = $store   ->group('Outputs')->group('Output1')->dataset('mergerTreeStartIndex')->get()->(($treeFinalIndex));
    my $storeTreeSize    = $store   ->group('Outputs')->group('Output1')->dataset('mergerTreeCount'     )->get()->(($treeFinalIndex));
    my $retrieveTreeSize = $retrieve->group('Outputs')->group('Output1')->dataset('mergerTreeCount'     )->get()->((             -1));

    # Set initial failed state to "not failed".
    my $failed        = 0;
    my $statusMessage = "";

    # Check that the number of nodes is the same.
    unless ( $storeTreeSize == $retrieveTreeSize ) {
	$statusMessage .= "FAILED: number of nodes in output changed after state retrieve\n";
	$failed         = 1;
    }

    # Get all available datasets.
    my @datasets = $storeData->datasets();

    # Check that each dataset is unchanged.
    foreach my $dataset ( @datasets ) {
	my $storeDataset    = $storeData   ->dataset($dataset)->get()->($storeTreeStart:$storeTreeStart+$storeTreeSize-1);
	my $retrieveDataset = $retrieveData->dataset($dataset)->get();
	my $equal = all($storeDataset == $retrieveDataset);
	unless ( $equal == 1 ) {
	    $statusMessage .= "FAILED: dataset '".$dataset."' changed after state retrieve\n";
	    $statusMessage .= "   before --> ".$storeDataset   ."\n";
	    $statusMessage .= "   after  --> ".$retrieveDataset."\n";
	    $failed         = 1;
	}
    }

    if ( $failed ) {
	# Test failed. If this is not the first pass, report failure - otherwise allow a second attempt.
	print $statusMessage
	    if ( $pass == 1 );
    } else {
	# Test succeeded - report this and finish.
	print "SUCCESS!\n";
	last;
    }

}

exit;
