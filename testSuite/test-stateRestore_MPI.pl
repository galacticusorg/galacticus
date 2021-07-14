#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use PDL;
use PDL::IO::HDF5;
use PDL::NiceSlice;
use Galacticus::Options;
use Data::Dumper;

# Run a set of Galacticus models to test the state store/retrieve functionality under MPI.
# Andrew Benson (15-Jun-2018)

# Read in any configuration options.
my $config;
if ( -e $ENV{'GALACTICUS_EXEC_PATH'}."/galacticusConfig.xml" ) {
    my $xml = new XML::Simple;
    $config = $xml->XMLin($ENV{'GALACTICUS_EXEC_PATH'}."/galacticusConfig.xml");
}

# Parse config options.
my $queueManager = &Galacticus::Options::Config(                'queueManager' );
my $queueConfig  = &Galacticus::Options::Config($queueManager->{'manager'     });

# Get any command line options.
my %options =
    (
     'processesPerNode' => exists($queueConfig->{'ppn'}) ? $queueConfig->{'ppn'} : 1,
    );
&Galacticus::Options::Parse_Options(\@ARGV,\%options);

# We need at least 8 processes to run this test.
if ( $options{'processesPerNode'} < 8 ) {
    print "SKIPPED: at least 8 processes per node are required for this test\n";
    exit;
}

# Run full store model.
system("export OMP_NUM_THREADS=1; rm -f outputs/state.state* outputs/state.gsl.state*; cd ..; mpirun -np 8 Galacticus.exe_MPI testSuite/parameters/state/store.xml"  );
die("FAILED: failed to run store model")
    unless ( $? == 0 );
# Find which threads ran the final tree.
my $finalTreeProcessMPI;
opendir(my $stateDirectory,"outputs");
while ( my $fileName = readdir($stateDirectory) ) {
    if  ( $fileName =~ m/state\.state\.log:MPI(\d+)/ ) {
	my $processMPI = $1;
	open(my $stateLogFile,"outputs/".$fileName);
	while (my $line = <$stateLogFile> ) {
	    if ( $line =~ m/^\s*Storing state for tree #(\d+)/ ) {
		if ( $1 == 15 ) {
		    $finalTreeProcessMPI = $processMPI;
		}
	    }
	}
	close($stateLogFile);
    }    
}
closedir($stateDirectory);
if ( defined($finalTreeProcessMPI) ) {
    print "Final tree was run on MPI process ".$finalTreeProcessMPI."\n";
    unless ( $finalTreeProcessMPI eq "0000" ) {
	system("cp -f outputs/state.state:MPI"    .$finalTreeProcessMPI." outputs/state.state:MPI0000"     );
	system("cp -f outputs/state.gsl.state:MPI".$finalTreeProcessMPI." outputs/state.gsl.state:MPI0000");
    }
} else {
    die("FAILED: failed to identify which thread/process ran final tree");
}

# Run the restore model.
system("export OMP_NUM_THREADS=1; cd ..; mpirun -np 1 Galacticus.exe_MPI testSuite/parameters/state/retrieve.xml");
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
	print "   before --> ".$storeDataset   ."\n";
	print "   after  --> ".$retrieveDataset."\n";
	$failed = 1;
    } else {
 	print "SUCCESS: dataset '".$dataset."'\n";
   }
}
print "SUCCESS!\n"
    unless ( $failed );
exit;
