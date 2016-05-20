#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V094"}) ) {
    $galacticusPath = $ENV{"GALACTICUS_ROOT_V094"};
    $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
    $galacticusPath = getcwd()."/../";
    $ENV{'GALACTICUS_ROOT_V094'} = $galacticusPath;
}
unshift(@INC,$galacticusPath."perl"); 
use PDL;
use PDL::IO::HDF5;
use PDL::NiceSlice;
require Galacticus::Options;

# Run a set of merger trees with forests split and not split. Compare the results which should be identical.
# Andrew Benson (19-May-2016)

# Find the location of the Millennium Database data.
my $millenniumDatabaseConfig = &Options::Config('millenniumDB');
if ( defined($millenniumDatabaseConfig) ) {
    if ( exists($millenniumDatabaseConfig->{'path'}) ) {
	print $millenniumDatabaseConfig->{'path'}."\n";
	unless ( -e $millenniumDatabaseConfig->{'path'}."/milliMillennium/milliMillennium.hdf5" ) {
	    print "SKIPPED: milli-Millennium data not available\n";
	    exit;
	}	    
    } else {
	print "SKIPPED: Millennium database data not available\n";
	exit;
    }
} else {
    print "SKIPPED: Millennium database location undefined\n";
    exit;
}

# Define the set of properties that we want to compare.
my @properties = ( "nodeIndex", "positionPositionX", "positionPositionY", "positionPositionZ", "positionVelocityX", "positionVelocityY", "positionVelocityZ", "satelliteNodeIndex", "satelliteMergeTime", "satelliteBoundMass", "nodeIsIsolated", "basicTimeLastIsolated", "parentIndex", "basicMass" );

# Define the two types of model we want to run.
my @types = ( 'split', 'unsplit' );

# Create directories needed.
system("mkdir -p outputs/test-splitForests");

# Locate a scratch directory.
my $scratchConfig = &Options::Config('scratch');

# Iterate over models, running them, reading in their data, and generating sort indices into the node indices.
my $data;
foreach ( @types ) {
    unless ( -e $_.'.hdf5' ) {
	print "Running ".$_." model...\n";
	# Read and modify parameter file.
	my $xml        = new XML::Simple(RootName => "parameters");
	my $parameters = $xml->XMLin("parameters/test-splitForests-".$_.".xml");
	$parameters->{'mergerTreeReadFileName'}->{'value'} = $millenniumDatabaseConfig->{'path'}."/milliMillennium/milliMillennium.hdf5";
	$parameters->{'treeEvolveSuspendPath' }->{'value'} = defined($scratchConfig) ? $scratchConfig->{'path'} : ".";
	my $outputFileName = "outputs/test-splitForests/".$_.".xml";
	open(my $outputFile,">".$outputFileName);
	print $outputFile $xml->XMLout($parameters);
	close($outputFile);
	# Run the model.
	system("cd ..; ./Galacticus.exe testSuite/".$outputFileName);
	unless ( $? == 0 ) { 
	    print "FAILED: model '".$_."' did not complete\n";
	    exit;
	}
    }
    $data->{$_}->{'file' } = new PDL::IO::HDF5("outputs/test-splitForests/".$_.'.hdf5');
    $data->{$_}->{'nodes'} = $data->{$_}->{'file'}->group("Outputs/Output1/nodeData");
    foreach my $property ( @properties ) {
	$data->{$_}->{'properties'}->{$property} = $data->{$_}->{'nodes'}->dataset($property)->get();
    }
    $data->{$_}->{'rank'} = $data->{$_}->{'properties'}->{'nodeIndex'}->qsorti();
}

# Test for equal numbers of nodes.
print "Testing for equal numbers of nodes...\n";
unless ( nelem($data->{'split'}->{'rank'}) == nelem($data->{'unsplit'}->{'rank'}) ) {
    for(my $i=0;$i<nelem($data->{'split'}->{'properties'}->{'nodeIndex'});++$i) {
	print "In split but not unsplit: ".$data->{'split'}->{'properties'}->{'nodeIndex'}->(($i))."\n"
	    unless ( any($data->{'unsplit'}->{'properties'}->{'nodeIndex'} == $data->{'split'}->{'properties'}->{'nodeIndex'}->(($i))) );
    }
    for(my $i=0;$i<nelem($data->{'unsplit'}->{'properties'}->{'nodeIndex'});++$i) {
	print "In unsplit but not split: ".$data->{'unsplit'}->{'properties'}->{'nodeIndex'}->(($i))."\n"
	    unless ( any($data->{'split'}->{'properties'}->{'nodeIndex'} == $data->{'unsplit'}->{'properties'}->{'nodeIndex'}->(($i))) );
    }
    print "FAILED: number of nodes differ\n";
    exit;
}
print "...done\n";

# Test for equal properties.
foreach my $property ( @properties ) {
    print "Testing for equal node ".$property."...\n";
    for(my $i=0;$i<nelem($data->{'split'}->{'rank'});++$i) {
	my $js = $data->{  'split'}->{'rank'}->(($i));
	my $ju = $data->{'unsplit'}->{'rank'}->(($i));
	unless ( $data->{'split'}->{'properties'}->{$property}->(($js)) == $data->{'unsplit'}->{'properties'}->{$property}->(($ju)) ) {
	    print $i."\t".$data->{'split'}->{'properties'}->{'nodeIndex'}->(($js))."\t".$data->{'unsplit'}->{'properties'}->{'nodeIndex'}->(($ju))."\t:\t".$data->{'split'}->{'properties'}->{$property}->(($js))."\t".$data->{'unsplit'}->{'properties'}->{$property}->(($ju))."\t:\t".$data->{'split'}->{'properties'}->{'nodeIsIsolated'}->(($js))."\t".$data->{'unsplit'}->{'properties'}->{'nodeIsIsolated'}->(($ju))."\n";
	    print "FAILED: '".$property."' mismatch\n";
	    exit;
	}
    }
    print "...done\n";
}
print "SUCCESS: split forests match unsplit forests\n";
exit;
