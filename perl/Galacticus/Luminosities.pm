# Contains a Perl module which implements total luminosity calculations for Galacticus.

package Luminosities;
use PDL;
require Galacticus::HDF5;
require Galacticus::DustAttenuation;
use Data::Dumper;
use XML::Simple;

%HDF5::galacticusFunctions = ( %HDF5::galacticusFunctions,
    "^totalLuminositiesStellar:([^:]+):([^:]+):z([\\d\\.]+)(:dust[^:]+)?" => \&Luminosities::Get_Luminosity,
    "^bulgeToTotalLuminosities:([^:]+):([^:]+):z([\\d\\.]+)(:dust[^:]+)?" => \&Luminosities::Get_BulgeToTotal
    );

sub Get_Luminosity {
    $dataSet = shift;
    $dataSetName = $_[0];

    # Check that the dataset name matches the expected regular expression.
    if ( $dataSetName =~ m/^totalLuminositiesStellar:([^:]+):([^:]+):z([\d\.]+)(:dust[^:]+)?/ ) {
	# Extract the dataset name information.
	$filter        = $1;
	$frame         = $2;
	$redshift      = $3;
	$dustExtension = $4;
	# Construct the name of the corresponding luminosity properties.
	$luminosityDataset[0] = "diskLuminositiesStellar:".$filter.":".$frame.":z".$redshift.$dustExtension;
	$luminosityDataset[1] = "spheroidLuminositiesStellar:".$filter.":".$frame.":z".$redshift.$dustExtension;
	&HDF5::Get_Dataset($dataSet,\@luminosityDataset);
	$dataSets = $dataSet->{'dataSets'};
	$dataSets->{$dataSetName} = $dataSets->{$luminosityDataset[0]}+$dataSets->{$luminosityDataset[1]};
    } else {
	die("Get_Luminosity(): unable to parse data set: ".$dataSetName);
    }
}

sub Get_BulgeToTotal {
    $dataSet = shift;
    $dataSetName = $_[0];
  
    # Check that the dataset name matches the expected regular expression.
    if ( $dataSetName =~ m/^bulgeToTotalLuminosities:([^:]+):([^:]+):z([\d\.]+)(:dust[^:]+)?/ ) {
	# Extract the dataset descriptor.
	$filter        = $1;
	$frame         = $2;
	$redshift      = $3;
	$dustExtension = $4;
	# Construct the name of the corresponding luminosity properties.
	$luminosityDataset[0] = "diskLuminositiesStellar:"    .$filter.":".$frame.":z".$redshift.$dustExtension;
	$luminosityDataset[1] = "spheroidLuminositiesStellar:".$filter.":".$frame.":z".$redshift.$dustExtension;
	&HDF5::Get_Dataset($dataSet,\@luminosityDataset);
	$dataSets = $dataSet->{'dataSets'};
	$dataSets->{$dataSetName} = $dataSets->{$luminosityDataset[1]}/($dataSets->{$luminosityDataset[0]}+$dataSets->{$luminosityDataset[1]});
	$nonluminous                 = which($dataSets->{$luminosityDataset[0]}+$dataSets->{$luminosityDataset[1]} <= 0.0);
	$dataSets->{$dataSetName}->index($nonluminous) .= 0.0;
    } else {
	die("Get_Luminosity(): unable to parse data set: ".$dataSetName);
    }
}

1;
