# Contains a Perl module which implements calculation of Lyman continuum luminosity in units of 10⁵⁰ photons/s.

package Lyc;
use strict;
use warnings;
use PDL;
use XML::Simple;
require Galacticus::HDF5;

%HDF5::galacticusFunctions = ( %HDF5::galacticusFunctions,
    "^(disk|spheroid)LymanContinuumLuminosity:z[\\d\\.]+\$" => \&Lyc::Get_Lyc_Luminosity
    );

sub Get_Lyc_Luminosity {
    my $model       = shift;
    my $dataSetName = $_[0];

    # Check to see if a particular postprocessing chain was specified for Lyman-continuum flux.
    my $postprocessingChain = "";
    $postprocessingChain = ":".$model->{'lymanContinuum'}->{'postProcessingChain'} 
      if ( exists($model->{'lymanContinuum'}->{'postProcessingChain'}) );

    # Required constants.
    my $plancksConstant = pdl 6.6260680000000000e-34; # J s
    my $luminosityAB    = pdl 4.4659201576470211e+13; # W/Hz
    my $lycUnits        = pdl 1.0000000000000000e+50;

    # Read the Lyman continuum filter.
    my $xml       = new XML::Simple;
    my $lycFilter = $xml->XMLin('./data/filters/Lyc.xml');
    my $wavelengthMaximum = pdl 0.0;
    my $wavelengthMinimum = pdl 1.0e30;
    foreach my $datum ( @{$lycFilter->{'response'}->{'datum'}} ) {
	$datum =~ s/^\s*//;
	my @columns = split(/\s+/,$datum);
	my $wavelength   = pdl $columns[0];
	my $transmission = pdl $columns[1];
	if ( $transmission > 0.0 ) {
	    $wavelengthMaximum = $wavelength
		if ( $wavelength > $wavelengthMaximum );
	    $wavelengthMinimum = $wavelength
		if ( $wavelength < $wavelengthMinimum );
	}
    }

    # Check that the dataset name matches the expected regular expression.
    if ( $dataSetName =~ m/^(disk|spheroid)LymanContinuumLuminosity:z([\d\.]+)$/ ) {
	# Extract the name of the component and redshift.
	my $component = $1;
	my $redshift  = $2;
	# Construct the name of the corresponding luminosity property.
	my $luminosityDataset = $component."StellarLuminosity:Lyc:rest:z".$redshift.$postprocessingChain;
	&HDF5::Get_Dataset($model,[$luminosityDataset]);
	my $dataSets = $model->{'dataSets'};
	$dataSets->{$dataSetName} = $dataSets->{$luminosityDataset}*($luminosityAB/$plancksConstant/$lycUnits)*log($wavelengthMaximum/$wavelengthMinimum);
    } else {
	die("Get_Luminosity(): unable to parse data set: ".$dataSetName);
    }
}

1;
