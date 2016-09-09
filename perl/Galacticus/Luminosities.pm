# Contains a Perl module which implements total luminosity calculations for Galacticus.

package Galacticus::Luminosities;
use strict;
use warnings;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use PDL;
use Data::Dumper;
use XML::Simple;
use Galacticus::HDF5;
use Galacticus::DustAttenuation;

%Galacticus::HDF5::galacticusFunctions = ( %Galacticus::HDF5::galacticusFunctions,
    "^totalLuminositiesStellar:([^:]+):([^:]+):z([\\d\\.]+)(:dust[^:]+)?" => \&Galacticus::Luminosities::Get_Luminosity,
    "^bulgeToTotalLuminosities:([^:]+):([^:]+):z([\\d\\.]+)(:dust[^:]+)?" => \&Galacticus::Luminosities::Get_BulgeToTotal
    );

sub Get_Luminosity {
    my $dataSet     = shift;
    my $dataSetName = $_[0];

    # Check that the dataset name matches the expected regular expression.
    if ( $dataSetName =~ m/^totalLuminositiesStellar:([^:]+):([^:]+):z([\d\.]+)(:dust[^:]+)?/ ) {
	# Extract the dataset name information.
	my $filter        = $1;
	my $frame         = $2;
	my $redshift      = $3;
	my $dustExtension = $4;
	$dustExtension = ""
	    unless ( defined($dustExtension) );
	# Construct the name of the corresponding luminosity properties.
	my @luminosityDataset;
	$luminosityDataset[0] = "diskLuminositiesStellar:"    .$filter.":".$frame.":z".$redshift.$dustExtension;
	$luminosityDataset[1] = "spheroidLuminositiesStellar:".$filter.":".$frame.":z".$redshift.$dustExtension;
	&Galacticus::HDF5::Get_Dataset($dataSet,\@luminosityDataset);
	my $dataSets = $dataSet->{'dataSets'};
	$dataSets->{$dataSetName} = $dataSets->{$luminosityDataset[0]}+$dataSets->{$luminosityDataset[1]};
    } else {
	die("Get_Luminosity(): unable to parse data set: ".$dataSetName);
    }
}

sub Get_BulgeToTotal {
    my $dataSet     = shift;
    my $dataSetName = $_[0];
  
    # Check that the dataset name matches the expected regular expression.
    if ( $dataSetName =~ m/^bulgeToTotalLuminosities:([^:]+):([^:]+):z([\d\.]+)(:dust[^:]+)?/ ) {
	# Extract the dataset descriptor.
	my $filter        = $1;
	my $frame         = $2;
	my $redshift      = $3;
	my $dustExtension = $4;
	# Construct the name of the corresponding luminosity properties.
	my @luminosityDataset;
	$luminosityDataset[0] = "diskLuminositiesStellar:"    .$filter.":".$frame.":z".$redshift.$dustExtension;
	$luminosityDataset[1] = "spheroidLuminositiesStellar:".$filter.":".$frame.":z".$redshift.$dustExtension;
	&Galacticus::HDF5::Get_Dataset($dataSet,\@luminosityDataset);
	my $dataSets = $dataSet->{'dataSets'};
	$dataSets->{$dataSetName} = $dataSets->{$luminosityDataset[1]}/($dataSets->{$luminosityDataset[0]}+$dataSets->{$luminosityDataset[1]});
	my $nonluminous                 = which($dataSets->{$luminosityDataset[0]}+$dataSets->{$luminosityDataset[1]} <= 0.0);
	$dataSets->{$dataSetName}->index($nonluminous) .= 0.0;
    } else {
	die("Get_Luminosity(): unable to parse data set: ".$dataSetName);
    }
}

1;
