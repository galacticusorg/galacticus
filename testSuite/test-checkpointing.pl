#!/usr/bin/env perl
use strict;
use warnings;
use PDL;
use PDL::IO::HDF5;
use PDL::NiceSlice;

# Run a set of Galacticus models to test checkpointing functionality.
# Andrew Benson (05-July-2023)

# Ensure output directory exists.
system("mkdir -p outputs");

# # Run model without checkpointing.
system("export OMP_NUM_THREADS=1; cd ..; ./Galacticus.exe testSuite/parameters/checkpointingNoCheckpoints.xml"  );
die("FAILED: failed to run model without checkpointing")
    unless ( $? == 0 );

# # Run model with checkpointing - interrupted.
system("export OMP_NUM_THREADS=1; cd ..; rm -f checkpoint.chk; ./Galacticus.exe testSuite/parameters/checkpointingCheckpoints.xml");
die("FAILED: failed to produce a checkpoint file")
    unless ( -e "../checkpoint.chk" );

# # Run model with checkpointing - resuming.
system("export OMP_NUM_THREADS=1; cd ..; ./Galacticus.exe testSuite/parameters/checkpointingResume.xml");
die("FAILED: failed to run model resuming from checkpoint")
    unless ( $? == 0 );

# Read model data.
my $data;
foreach my $modelName ( "chechkpointingNoCheckpoints", "chechkpointingCheckpoints" ) {
    my $model = new PDL::IO::HDF5("outputs/".$modelName.".hdf5");
    my $nodes = $model->group("Outputs/Output1/nodeData");
    $data->{$modelName}->{$_} = $nodes->dataset($_)->get()
	foreach ( "nodeIndex", "basicMass" );
}

# Compare results.
my $status = "SUCCESS";
for(my $i=0;$i<nelem($data->{'chechkpointingNoCheckpoints'}->{'nodeIndex'});++$i) {
    my $match = which($data->{'chechkpointingCheckpoints'}->{'nodeIndex'} == $data->{'chechkpointingNoCheckpoints'}->{'nodeIndex'}->(($i)));
    next
	unless ( nelem($match) == 1 );
    my $massNoCheckpoints = $data->{'chechkpointingNoCheckpoints'}->{'basicMass'}          ->(($i));
    my $massCheckpoints   = $data->{'chechkpointingCheckpoints'  }->{'basicMass'}->($match)->(( 0));
    $status = "FAIL"
	unless ( $massNoCheckpoints == $massCheckpoints );
}
print $status.": resume from chekpoint file\n";

exit;
