# Contains a Perl module which implements dust attenuated luminosity calculations for Galacticus using the model of Charlot & Fall
# (2000).

package DustCharlotFall2000;
use strict;
use warnings;
use XML::Simple;
use PDL;
use PDL::NiceSlice;
require Galacticus::HDF5;

%HDF5::galacticusFunctions = ( %HDF5::galacticusFunctions,
    "^(disk|spheroid)(StellarLuminosity|CentralOpticalDepthISM|CentralOpticalDepthClouds):.*:dustCharlotFall2000\$" => \&DustCharlotFall2000::Get_Dust_Charlot_Fall_2000_Luminosity
    );

sub Get_Dust_Charlot_Fall_2000_Luminosity {
    my $model       = shift;
    my $dataSetName = $_[0];

    # Specify constants in the extinction model.
    #  * "wavelengthExponent" is set to the value of 0.7 found by Charlot & Fall (2000).
    #  * "opticalDepthCloudsFactor" is set to unity, such that in gas with Solar metallicity the cloud optical depth will be 1.
    #  * "opticalDepthISMFactor" is set to 1.0 such that we reproduce the standard (Bohlin et al 1978) relation between visual
    #      extinction and column density in the local ISM (essentially solar metallicity).
    my $opticalDepthISMFactor    = pdl    1.0;
    my $opticalDepthCloudsFactor = pdl    1.0;
    my $wavelengthZeroPoint      = pdl 5500.0; # Angstroms
    my $wavelengthExponent       = pdl    0.7;
    
    # Extract any user-specified values of the model parameters.
    if ( exists($model->{'dustCharlotFall2000'}) ) {
	$opticalDepthISMFactor    = $model->{'dustCharlotFall2000'}->{'opticalDepthISMFactor'   }
	  if ( exists($model->{'dustCharlotFall2000'}->{'opticalDepthISMFactor'   }) );
	$opticalDepthCloudsFactor = $model->{'dustCharlotFall2000'}->{'opticalDepthCloudsFactor'}
	  if ( exists($model->{'dustCharlotFall2000'}->{'opticalDepthCloudsFactor'}) );
	$wavelengthExponent       = $model->{'dustCharlotFall2000'}->{'wavelengthExponent'      }
	  if ( exists($model->{'dustCharlotFall2000'}->{'wavelengthExponent'      }) );
    }

    # Get the name of the unattenuated luminosity dataset.
    (my $luminosityDataSet = $dataSetName) =~ s/:dustCharlotFall2000//;

    # Get the name of the corresponding recent luminosity dataset.
    my $recentLuminosityDataSet = $luminosityDataSet.":recent";

    # Determine if this is a disk or a spheroid luminosity.
    my $component;
    if ( $luminosityDataSet =~ m/^(disk|spheroid)/ ) {
	$component = $1;
    } else {
	die("DustCharlotFall2000::Get_Dust_Charlot_Fall_2000_Luminosity(): can not determine if a disk or spheroid luminosity was requested");
    }

    # Extract filter data.
    my $frame;
    my $filter;
    my $redshift;
    my $opticalDepthOrLuminosity;
    if ( $dataSetName =~ m/^(.*?)(StellarLuminosity|CentralOpticalDepthISM|CentralOpticalDepthClouds):([^:]+):([^:]+):z([\d\.]+)/ ) {
	# Extract the dataset name information.
	$component                = $1;
	$opticalDepthOrLuminosity = $2;
	$filter                   = $3;
	$frame                    = $4;
	$redshift                 = $5;
    } else {
	die ("DustCharlotFall2000::Get_Dust_Charlot_Fall_2000_Luminosity(): unable to parse dataset name");
    }

    # List of properties to read.
    my @propertyList = ($component."GasMetals",$component."GasMass",$component."ScaleLength");
    push(@propertyList,$luminosityDataSet,$recentLuminosityDataSet)
	unless ( $dataSetName =~ m/CentralOpticalDepth/ );

    # Get the datasets needed for our calculation.
    &HDF5::Get_Dataset($model,\@propertyList);
    my $dataSets = $model->{'dataSets'};

    # Ensure we have an effective wavelength and a wavelength index for this filter.
    my $filterPath = "./data/filters/".$filter.".xml";
    die("DustCharlotFall2000::Get_Dust_Charlot_Fall_2000_Luminosity(): can not find filter file for: ".$filter) unless ( -e $filterPath );
    my $xml = new XML::Simple;
    my $filterData = $xml->XMLin($filterPath);
    unless ( exists($filterData->{'effectiveWavelength'}) ) {
	# No effective wavelength data available for filter - run the script that computes it.
	system("scripts/filters/vega_offset_effective_lambda.pl");
	$filterData = $xml->XMLin($filterPath);
	die ("DustCharlotFall2000::Get_Dust_Charlot_Fall_2000_Luminosity(): failed to compute effective wavelengths for filters") unless ( exists($filterData->{'effectiveWavelength'}) );
    }
    my $effectiveWavelength;
    if ( $frame eq "rest" ) {
	$effectiveWavelength = pdl $filterData->{'effectiveWavelength'};
    } else {
	$effectiveWavelength = pdl $filterData->{'effectiveWavelength'}/(1.0+$redshift);
    }

    # Specify constants used in calculation of central optical depths.
    my $localISMMetallicity       = pdl 0.02;                # Metallicity in the local ISM.
    my $Pi                        = pdl 3.1415926536;        # Pi.
    my $mega                      = pdl 1.0e6;               # SI prefix.
    my $hecto                     = pdl 1.0e2;               # SI prefix.
    my $parsec                    = pdl 3.0856775807e16;     # Parsec in m.
    my $atomicMass                = pdl 1.66053873e-27;      # Atomic mass in kg.
    my $solarMass                 = pdl 1.9891e30;           # Solar mass in kg.
    my $AV_to_EBV                 = pdl 3.10;                # (A_V/E(B-V); Savage & Mathis 1979)
    my $NH_to_EBV                 = pdl 5.8e21;              # (N_H/E(B-V); atoms/cm^2/mag; Savage & Mathis 1979)
    my $hydrogenMassFractionSolar = pdl 0.707;               # Hydrogen mass fraction in the Sun.
    my $opticalDepthToMagnitudes  = pdl 2.5*log10(exp(1.0)); # Conversion factor from optical depth to magnitudes of extinction.
    my $opticalDepthNormalization = (1.0/$opticalDepthToMagnitudes)*($AV_to_EBV/$NH_to_EBV)*($hydrogenMassFractionSolar/$atomicMass)*($solarMass/($parsec*$hecto)**2)/$localISMMetallicity;

    # Central surface density in M_Solar/pc^2.
    my $gasMetalMass                   = pdl $dataSets->{$component."GasMetals"};
    my $scaleLength                    = pdl $dataSets->{$component."ScaleLength"};
    my $noGasMetals                    = which($gasMetalMass <= 0.0);
    my $notPresent                     = which($scaleLength  <= 0.0);
    my $gasMetallicity                 = pdl $dataSets->{$component."GasMetals"}/$dataSets->{$component."GasMass"};
    $gasMetallicity->index($noGasMetals) .= 0.0;            
    my $gasMetalsSurfaceDensityCentral = $gasMetalMass/(2.0*$Pi*($mega*$scaleLength)**2);
    $gasMetalsSurfaceDensityCentral->index($notPresent) .= 0.0;

    # Compute central optical depths.
    my $opticalDepth                   = $opticalDepthNormalization*$gasMetalsSurfaceDensityCentral/($effectiveWavelength/$wavelengthZeroPoint)**$wavelengthExponent;
    my $opticalDepthISM                = $opticalDepthISMFactor*$opticalDepth;
    my $opticalDepthClouds             = $opticalDepthCloudsFactor*$gasMetallicity/$localISMMetallicity/($effectiveWavelength/$wavelengthZeroPoint)**$wavelengthExponent;

    # Compute the attenutations
    my $attenuationsISM                = exp(-$opticalDepthISM);
    my $attenuationsClouds             = exp(-$opticalDepthClouds);

    # For luminosities, multiply luminosities by attenuations.
    $dataSets->{$dataSetName} = (
	                          $dataSets->{$luminosityDataSet}-$dataSets->{$recentLuminosityDataSet}
	                                                         +$dataSets->{$recentLuminosityDataSet}
	                                                         *$attenuationsClouds
	                        )
	                        *$attenuationsISM
				if (  $opticalDepthOrLuminosity eq "StellarLuminosity" );

    # For optical depths, return the appropriate dataset.
    $dataSets->{$dataSetName} = $opticalDepthISM
	if ( $opticalDepthOrLuminosity eq "CentralOpticalDepthISM" );
    $dataSets->{$dataSetName} = $opticalDepthClouds
	if ( $opticalDepthOrLuminosity eq "CentralOpticalDepthClouds" );

}

1;
