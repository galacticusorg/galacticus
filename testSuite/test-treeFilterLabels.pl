#!/usr/bin/env perl
use strict;
use warnings;
use PDL;
use PDL::IO::HDF5;
use PDL::NiceSlice;

# Run a Galacticus model to test node labeling functionality.
# Andrew Benson (07-July-2023)

# Ensure output directory exists.
system("mkdir -p outputs");

# Run model without checkpointing.
system("cd ..; ./Galacticus.exe testSuite/parameters/treeFilterLabels.xml"  );
die("FAILED: failed to run model")
    unless ( $? == 0 );

# Read model data.
my $data;
my $model = new PDL::IO::HDF5("outputs/treeFilterLabels.hdf5");
my $nodes = $model->group("Outputs/Output1/nodeData");
die("FAILED: label not output")
    unless ( grep {$_ eq "nodeLabelLMC"} $nodes->datasets() );
$data->{$_} = $nodes->dataset($_)->get()
	foreach ( "nodeIsIsolated", "nodeLabelLMC", "darkMatterProfileDMOVelocityMaximum", "basicTimeLastIsolated" );
my $hosts = which($data->{'nodeIsIsolated'} == 1);
my $lmcs  = which($data->{'nodeLabelLMC'  } == 1);
print "Found ".nelem($lmcs)." LMCs in ".nelem($hosts)." trees\n";
die("FAILED: LMCs not found in all trees")
    unless ( nelem($lmcs) >= nelem($hosts) );
die("FAILED: LMCs labelled prior to infall time")
    if ( any($data->{'basicTimeLastIsolated'              }->($lmcs) < 11.8) );
die("FAILED: LMCs labelled at low Vₘₐₓ")
    if ( any($data->{'darkMatterProfileDMOVelocityMaximum'}->($lmcs) < 55.0) );
print "SUCCESS: filtered tree labeling\n";

exit;
