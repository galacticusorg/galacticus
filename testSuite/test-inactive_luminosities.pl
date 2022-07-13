#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;

# Test inactive luminosity calculations give results identical to active luminosity calculations.
# Andrew Benson (14-May-2019)

# Specify types.
my @types = ( 'active', 'inactive' );

# Run the models.
system("mkdir -p outputs/inactiveLuminosities");
foreach my $type ( @types ) {
    system("cd ..; ./Galacticus.exe testSuite/parameters/".$type."Luminosities.xml");
    die("FAILED: model '".$type."' failed to run")
	unless ( $? == 0 );
}

# Extraxt datasets and compare.
my $models;
foreach my $type ( @types ) {
    $models->{$type}->{'file'   } = new PDL::IO::HDF5("outputs/inactiveLuminosities/".$type."Luminosities.hdf5");
    $models->{$type}->{'outputs'} = $models->{$type}->{'file'}->group('Outputs');
}
my $success = 1;
print "Status\tDataset\tFractional error\n";
for(my $i=1;$i<=4;++$i) {
    foreach my $type ( @types ) {
	$models->{$type}->{'nodes'.$i} = $models->{$type}->{'outputs'}->group('Output'.$i)->group('nodeData');
    }
    my @datasetNames = $models->{'active'}->{'nodes'.$i}->datasets();
    foreach my $datasetName ( @datasetNames ) {
	next
	    unless
	    (
	     $datasetName =~ m/^(disk|spheroid)MassStellar$/
	     ||
	     $datasetName =~ m/^(disk|spheroid)LuminositiesStellar/
	    );
	my $property;
	foreach my $type ( @types ) {
	    $property->{$type} = $models->{$type}->{'nodes'.$i}->dataset($datasetName)->get();
	}
	my $errorFractional = abs($property->{'inactive'}-$property->{'active'})/$property->{'active'};
	my $status = "SUCCESS";
	if ( $errorFractional->((0)) > 1.5e-3 ) {
	    $status  = "FAILED";
	    $success = 0;
	}
	print $status."\t".$datasetName."\t".$errorFractional->((0))."\n";
    }
}
my $status = $success ? "SUCCESS" : "FAILED: some datasets differed";
print "\n".$status."\n";
exit;
