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
use Galacticus::EmissionLines;
use Galacticus::Options;
use Galacticus::Path;
 
# Compute emission line luminosities for redshifts in a Galacticus output.
# Andrew Benson (21-September-2016)

# Get name of Galacticus file.
die("emissionLinesFilters.pl <galacticusFile>")
    unless ( scalar(@ARGV) >= 1 );
my $galacticusFileName = $ARGV[0];

# Parse any command line options.
my %options =
    (
     overwrite => "false"
    );
&Galacticus::Options::Parse_Options(\@ARGV,\%options);

# Get the names of all available lines.
my $emissionLineTable = new PDL::IO::HDF5(&galacticusPath()."data/hiiRegions/emissionLines.hdf5");
my $linesGroup        = $emissionLineTable->group   ('lines');
my @lineNames         = $linesGroup       ->datasets(       );

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
    my @properties;
    my $outputGroup = $galacticus->{'hdf5File'}->group('Outputs/Output'.$galacticus->{'output'}.'/nodeData');
    foreach my $datasetName ( keys(%{$galacticus->{'dataSetsAvailable'}}) ) {
	# Identify luminosity datasets.
	if ( $datasetName =~ m/^(disk|spheroid)LuminositiesStellar:([^:]+):(rest|observed):z([\d\.]+)$/ ) {
	    # Generate corresponding dusty property name.
	    my $component = $1;
	    my $filter    = $2;
	    my $frame     = $3;
	    my $redshift  = $4;
	    unless ( $filter eq "Lyc" || $filter =~ m/Continuum$/ ) {
		foreach my $lineName ( @lineNames ) {
		    my $emissionLineDatasetName = $1."LineLuminosity:".$lineName.":".$filter.":".$frame.":z".$redshift;
		    push(
			@properties,
			$emissionLineDatasetName
			);
		    # If overwriting, remove the dataset.
		    if ( $options{'overwrite'} eq "true" && exists($galacticus->{'dataSetsAvailable'}->{$emissionLineDatasetName}) ) {
			$outputGroup->unlink($emissionLineDatasetName);
			delete($galacticus->{'dataSetsAvailable'}->{$emissionLineDatasetName});
		    }
		}
	    }
	}
    }
    # Retrieve (and store) emission line luminosities.
    &Galacticus::HDF5::Get_Dataset($galacticus,\@properties);
    # Drop all properties.
    delete($galacticus->{'dataSets'});
}

exit;
