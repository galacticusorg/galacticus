#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use Data::Dumper;
use Galacticus::HDF5;
use Galacticus::AGNLuminosities;
use Galacticus::Options;
use Galacticus::Path;
 
# Compute AGN luminosities for redshifts in a Galacticus output.
# Andrew Benson (15-February-2017)

# Get name of Galacticus file.
die("agnFilters.pl <galacticusFile>")
    unless ( scalar(@ARGV) >= 1 );
my $galacticusFileName = $ARGV[0];

# Parse any command line options.
my %options =
    (
     overwrite    => "false",
     noAbsorption => "false"
    );
&Galacticus::Options::Parse_Options(\@ARGV,\%options);

# Create data structure to read the results.
my $galacticus;
$galacticus->{'file' } = $galacticusFileName;
$galacticus->{'store'} = 1;
$galacticus->{'tree' } = "all";

# Iterate over redshifts.
&Galacticus::HDF5::Get_Times($galacticus);
foreach my $redshift ( $galacticus->{'outputs'}->{'redshift'}->list() ) {
    print "Processing z=".$redshift."\n";
    &Galacticus::HDF5::Select_Output         ($galacticus,$redshift);
    &Galacticus::HDF5::Get_Datasets_Available($galacticus          );
    # Iterate over datasets.
    my %properties;
    my $outputGroup = $galacticus->{'hdf5File'}->group('Outputs/Output'.$galacticus->{'output'}.'/nodeData');
    foreach my $datasetName ( keys(%{$galacticus->{'dataSetsAvailable'}}) ) {
	# Identify luminosity datasets.
	if ( $datasetName =~ m/^(disk|spheroid)LuminositiesStellar:([^:]+):(rest|observed):z([\d\.]+)$/ ) {
	    # Generate corresponding dusty property name.
	    my $component      = $1;
	    my $filter         = $2;
	    my $frame          = $3;
	    my $redshift       = $4;
	    my $agnDatasetName = "agnLuminosity:".$filter.":".$frame.":z".$redshift;
	    $agnDatasetName   .= ":noAbsorption"
		if ( $options{'noAbsorption'} eq "true" );
	    $properties{$agnDatasetName} = 1;
	}
    }
    # If overwriting, remove the dataset.
    if ( $options{'overwrite'} eq "true" ) {
	foreach my $datasetName ( keys(%properties) ) {
	    if ( exists($galacticus->{'dataSetsAvailable'}->{$datasetName}) ) {
		$outputGroup->unlink($datasetName);
		delete($galacticus->{'dataSetsAvailable'}->{$datasetName});
	    }
	}
    }
    # Retrieve (and store) AGN luminosities.
    my @propertyList = keys(%properties);
    &Galacticus::HDF5::Get_Dataset($galacticus,\@propertyList);
    # Drop all properties.
    delete($galacticus->{'dataSets'});
}

exit;
