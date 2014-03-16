# Contains a Perl module which implements magnitude calculations for Galacticus.

package Magnitudes;
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V093"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V093"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
}
unshift(@INC,$galacticusPath."perl"); 
use PDL;
use Data::Dumper;
use XML::Simple;
require Galacticus::HDF5;
require Galacticus::DustAttenuation;
require Galacticus::Luminosities;
require Galacticus::Survey;

%HDF5::galacticusFunctions = ( %HDF5::galacticusFunctions,
    "^magnitude([^:]+):([^:]+):([^:]+):z([\\d\\.]+)(:dust[^:]+)?(:vega|:AB)?" => \&Magnitudes::Get_Magnitude                 ,
    "^magnitude:.*(:vega|:AB)?"                                               => \&Magnitudes::Get_Generic_Magnitude         ,
    "^apparentMagnitude:.*"                                                   => \&Magnitudes::Get_Generic_Apparent_Magnitude
    );

my %vegaOffsets;
my $vegaMagnitude;
my $filter;

sub Get_Magnitude {
    my $dataSet = shift;
    my $dataSetName = $_[0];
    # Check that the dataset name matches the expected regular expression.
    if ( $dataSetName =~ m/^magnitude([^:]+):([^:]+):([^:]+):z([\d\.]+)(:dust[^:]+)?(:vega|:AB)?/ ) {
	# Extract the dataset name information.
	my $component     = $1;
	$filter           = $2;
	my $frame         = $3;
	my $redshift      = $4;
	my $dustExtension = $5;
	$dustExtension = ""
	    unless ( defined($dustExtension) );
	if ( defined($6) && $6 eq ":vega" ) {
	    $vegaMagnitude = 1;
	} else {
	    $vegaMagnitude = 0;
	}
	# Construct the name of the corresponding luminosity property.
	my $luminosityDataset = lc($component)."LuminositiesStellar:".$filter.":".$frame.":z".$redshift.$dustExtension;
	&HDF5::Get_Dataset($dataSet,[$luminosityDataset]);
	my $dataSets = $dataSet->{'dataSets'};
	$dataSets->{$dataSetName} = -2.5*log10($dataSets->{$luminosityDataset}+1.0e-40);
	# If a Vega magnitude was requested, add the appropriate offset.
	if ( $vegaMagnitude == 1 ) {
	    unless ( exists($vegaOffsets{$filter}) ) {
		my $filterPath = $galacticusPath."data/filters/".$filter.".xml";
		die("Get_Magnitudes(): can not find filter file for: ".$filter) unless ( -e $filterPath );
		my $xml = new XML::Simple;
		my $filterData = $xml->XMLin($filterPath);
		unless ( exists($filterData->{'vegaOffset'}) ) {
		    # No Vega offset data available for filter - run the script that computes it.
		    system("scripts/filters/vega_offset_effective_lambda.pl");
		    $filterData = $xml->XMLin($filterPath);
		    die ("Get_Magnitudes(): failed to compute Vega offsets for filters") unless ( exists($filterData->{'vegaOffset'}) );
		}
		$vegaOffsets{$filter} = pdl $filterData->{'vegaOffset'};
	    }
	    $dataSets->{$dataSetName} += $vegaOffsets{$filter};
	}
    } else {
	die("Get_Magnitude(): unable to parse data set: ".$dataSetName);
    }
}

sub Get_Generic_Magnitude {
    my $dataSet     = shift;
    my $dataSetName = $_[0];
    # Construct the name of the corresponding luminosity property.
    (my $luminosityDataset = $dataSetName) =~ s/^magnitude:/luminosity:/;
    $luminosityDataset =~ s/(:vega|:AB)$//;
    &HDF5::Get_Dataset($dataSet,[$luminosityDataset]);
    my $dataSets = $dataSet->{'dataSets'};
    $dataSets->{$dataSetName} = -2.5*log10($dataSets->{$luminosityDataset}+1.0e-40);
    # If a Vega magnitude was requested, add the appropriate offset.
    if ( $vegaMagnitude == 1 ) {
	unless ( exists($vegaOffsets{$filter}) ) {
	    my $filterPath = $galacticusPath."data/filters/".$filter.".xml";
	    die("Get_Magnitudes(): can not find filter file for: ".$filter) unless ( -e $filterPath );
	    my $xml = new XML::Simple;
	    my $filterData = $xml->XMLin($filterPath);
	    unless ( exists($filterData->{'vegaOffset'}) ) {
		# No Vega offset data available for filter - run the script that computes it.
		system("scripts/filters/vega_offset_effective_lambda.pl");
		$filterData = $xml->XMLin($filterPath);
		die ("Get_Magnitudes(): failed to compute Vega offsets for filters") unless ( exists($filterData->{'vegaOffset'}) );
	    }
	    $vegaOffsets{$filter} = pdl $filterData->{'vegaOffset'};
	}
	$dataSets->{$dataSetName} += $vegaOffsets{$filter};
    }
}

sub Get_Generic_Apparent_Magnitude {
    my $dataSet     = shift;
    my $dataSetName = $_[0];
    # Construct the name of the corresponding absolute magnitude property.
    (my $absoluteMagnitudeDataset = $dataSetName) =~ s/^apparentMagnitude:/magnitude:/;
    &HDF5::Get_Dataset($dataSet,[$absoluteMagnitudeDataset,"distanceModulus"]);
    my $dataSets = $dataSet->{'dataSets'};
    $dataSets->{$dataSetName} = $dataSets->{$absoluteMagnitudeDataset}+$dataSets->{'distanceModulus'};
}

1;
