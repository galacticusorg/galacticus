# Contains a Perl module which implements calculation of Lyman continuum luminosity in units of
# 10⁵⁰ photons/s.

# Contributions to this file from: Andrew Benson; Christoph Behrens.

package IonizingContinuua;
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V094"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V094"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
}
use PDL;
use XML::Simple;
require Galacticus::HDF5;

%HDF5::galacticusFunctions = ( %HDF5::galacticusFunctions,
    "^(disk|spheroid|total)(Lyman|Helium|Oxygen)ContinuumLuminosity:z[\\d\\.]+\$" => \&IonizingContinuua::Get_Ionizing_Luminosity
    );

sub Get_Ionizing_Luminosity {
    my $model       = shift;
    my $dataSetName = $_[0];
    # Check to see if a particular postprocessing chain was specified for Lyman-continuum flux.
    my $postprocessingChain = "";
    $postprocessingChain = ":".$model->{'ionizingContinuum'}->{'postProcessingChain'} 
        if ( exists($model->{'ionizingContinuum'}->{'postProcessingChain'}) );
    # Required constants.
    my $plancksConstant = pdl 6.6260680000000000e-34; # J s
    my $luminosityAB    = pdl 4.4659201576470211e+13; # W/Hz
    my $continuumUnits  = pdl 1.0000000000000000e+50; # photons/s
    # Continuum filter names.
    my %filterName =
	(
	 Lyman  => "Lyc",
	 Oxygen => "OxygenContinuum",
	 Helium => "HeliumContinuum"
	);
    # Find the continuum requested.
    if ( $dataSetName =~ m/^(disk|spheroid|total)(Lyman|Helium|Oxygen)ContinuumLuminosity:z([\d\.]+)$/ ) {
	# Extract the name of the component and redshift.
	my $component     = $1;
	my $continuumName = $2;
	my $redshift      = $3;
	# Read the relevant continuum filter.
	my $xml       = new XML::Simple;
	my $continuumFilter = $xml->XMLin($galacticusPath."data/filters/".$filterName{$continuumName}.".xml");
	my $wavelengthMaximum = pdl 0.0;
	my $wavelengthMinimum = pdl 1.0e30;
	foreach my $datum ( @{$continuumFilter->{'response'}->{'datum'}} ) {
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
	# Construct the name of the corresponding luminosity properties.
	my @luminosityDatasets;
	if ( $component eq "total" ) {
	    push(
		@luminosityDatasets,
		"diskLuminositiesStellar:".$filterName{$continuumName}.":rest:z".$redshift.$postprocessingChain,
		"spheroidLuminositiesStellar:".$filterName{$continuumName}.":rest:z".$redshift.$postprocessingChain
		);
	} else {
	    push(
		@luminosityDatasets,
		$component."LuminositiesStellar:".$filterName{$continuumName}.":rest:z".$redshift.$postprocessingChain
		);
	}
	&HDF5::Get_Dataset($model,\@luminosityDatasets);
	my $dataSets = $model->{'dataSets'};
	foreach ( @luminosityDatasets ) {
	    $dataSets->{$dataSetName} += $dataSets->{$_}*($luminosityAB/$plancksConstant/$continuumUnits)*log($wavelengthMaximum/$wavelengthMinimum);
	}
    } else {
	die("Galacticus::IonizingContinuua::Get_Ionizing_Luminosity: unrecognized continuum name ".$dataSetName);
    }
}

1;
