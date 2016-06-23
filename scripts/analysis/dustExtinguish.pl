#!/usr/bin/env perl
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V094"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V094"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
}
unshift(@INC,$galacticusPath."perl"); 
use PDL;
use PDL::NiceSlice;
use Data::Dumper;
require Galacticus::HDF5;
require Galacticus::DustAttenuation;
require Galacticus::Options;
 
# Compute dust-extinguished luminosities for all luminosities in a Galacticus output.
# Andrew Benson (10-March-2016)

# Get name of Galacticus file.
die("dustExtinguish.pl <galacticusFile>")
    unless ( scalar(@ARGV) >= 1 );
my $galacticusFileName = $ARGV[0];

# Parse any command line options.
my %options =
    (
     overwrite => "false"
    );
&Options::Parse_Options(\@ARGV,\%options);

# Create data structure to read the results.
my $galacticus;
$galacticus->{'file' } = $galacticusFileName;
$galacticus->{'store'} = 1;
$galacticus->{'tree' } = "all";

# Iterate over redshifts.
&HDF5::Get_Times($galacticus);
foreach my $redshift ( $galacticus->{'outputs'}->{'redshift'}->list() ) {
    print "Processing z=".$redshift."\n";
    &HDF5::Select_Output         ($galacticus,$redshift);
    &HDF5::Get_Datasets_Available($galacticus          );
    # Iterate over datasets.
    my @properties;
    my $outputGroup = $galacticus->{'hdf5File'}->group('Outputs/Output'.$galacticus->{'output'}.'/nodeData');
    foreach my $datasetName ( keys(%{$galacticus->{'dataSetsAvailable'}}) ) {
	# Identify luminosity datasets.
	if ( $datasetName =~ m/^(disk|spheroid)LuminositiesStellar:([^:]+):(rest|observed):z[\d\.]+$/ ) {
	    # Generate corresponding dusty property name.
	    my $dustyDatasetName = $datasetName.":dustAtlas";
	    push(
		@properties,
		$dustyDatasetName
		);
	    # If overwriting, remove the dataset.
	    if ( $options{'overwrite'} eq "true" && exists($galacticus->{'dataSetsAvailable'}->{$dustyDatasetName}) ) {
		$outputGroup->unlink($dustyDatasetName);
		delete($galacticus->{'dataSetsAvailable'}->{$dustyDatasetName});
	    }
	}
    }
    # Retrieve (and store) dust luminosities.
    &HDF5::Get_Dataset($galacticus,\@properties);
    # Drop all properties.
    delete($galacticus->{'dataSets'});
}

exit;
