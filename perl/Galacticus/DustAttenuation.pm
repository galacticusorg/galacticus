# Contains a Perl module which implements dust attenuated luminosity calculations for Galacticus.

package DustAttenuation;
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
use PDL::NiceSlice;
use XML::Simple;
require Galacticus::HDF5;
require Galacticus::Inclination;

%HDF5::galacticusFunctions = ( %HDF5::galacticusFunctions,
    "^(disk|spheroid)LuminositiesStellar:.*:dustAtlas(\\[faceOn\\])?\$" => \&DustAttenuation::Get_Dust_Attenuated_Luminosity
    );

# Flag indicating whether dust data is loaded yet.
my $dustDataLoaded = 0;
my %effectiveWavelength;
my %wavelengthIndex;
my $wavelengths;
my $wavelengthCount;
my $wavelengthIndices;
my $diskInclinationsCount;
my $diskInclinations;
my $diskInclinationIndices;
my $spheroidInclinationsCount;
my $spheroidInclinations;
my $spheroidInclinationIndices;
my $diskOpticalDepths;
my $diskOpticalDepthsCount;
my $diskOpticalDepthsIndices;
my $spheroidOpticalDepths;
my $spheroidOpticalDepthsCount;
my $spheroidOpticalDepthsIndices;
my $spheroidSizes;
my $spheroidSizesCount;
my $spheroidSizeIndices;
my $diskAttenuations;
my $spheroidAttenuations;

sub Get_Dust_Attenuated_Luminosity {
    my $dataSet     = shift;
    my $dataSetName = $_[0];

    # Check parameters.
    &HDF5::Get_Parameters($dataSet);
    die ("Get_Dust_Attenuated_Luminosity(): routine assumes exponential disks and Hernquist or Sersic spheroids")
	unless
	(
	 $dataSet ->{'parameters'}->{"treeNodeMethodDisk"      } eq "exponential" &&
	 (
	  $dataSet->{'parameters'}->{"spheroidMassDistribution"} eq "hernquist" ||
	  $dataSet->{'parameters'}->{"spheroidMassDistribution"} eq "sersic" 
	 )
	);

    # Get the name of the unattenuated luminosity dataset.
    (my $luminosityDataSet = $dataSetName) =~ s/:dustAtlas(\[faceOn\])?//;

    # List of properties to read.
    my @propertyList = ("diskAbundancesGasMetals","diskRadius","spheroidRadius",$luminosityDataSet);

    # Check if a face-on magnitude is required.
    my $faceOn;
    if ( $dataSetName =~ m/\[faceOn\]/ ) {
	$faceOn = 1;
    } else {
	$faceOn = 0;
	push(@propertyList,"inclination");
    }

    # Get the datasets needed for our calculation.
    &HDF5::Get_Dataset($dataSet,\@propertyList);
    my $dataSets = $dataSet->{'dataSets'};

    # Load the dust atlas data.
    &Load_Dust_Atlas($dataSet);

    # Determine types of extrapolation to be used.
    my $extrapolateInSize;
    if ( exists($dataSet->{'dustAtlasExtrapolateInSize'}) ) {
	$extrapolateInSize = $dataSet->{'dustAtlasExtrapolateInSize'};
	if ( lc($extrapolateInSize) eq "yes" ) {$extrapolateInSize = 1};
	if ( lc($extrapolateInSize) eq "no"  ) {$extrapolateInSize = 0};
    } else {
	# Assume that we will extrapolate in spheroid sizes by default.
	$extrapolateInSize = 1;
    }
    my $extrapolateInTau;
    if ( exists($dataSet->{'dustAtlasExtrapolateInTau'}) ) {
	$extrapolateInTau = $dataSet->{'dustAtlasExtrapolateInTau'};
	if ( lc($extrapolateInTau) eq "yes" ) {$extrapolateInTau = 1};
	if ( lc($extrapolateInTau) eq "no"  ) {$extrapolateInTau = 0};
    } else {
	# Assume that we will extrapolate in optical depth by default.
	$extrapolateInTau = 1;
    }

    # Extract filter data.
    my $component;
    my $filter;
    my $frame;
    my $redshift;
    if ( $dataSetName =~ m/^(.*?)LuminositiesStellar:([^:]+):([^:]+):z([\d\.]+)/ ) {
	# Extract the dataset name information.
	$component = $1;
	$filter    = $2;
	$frame     = $3;
	$redshift  = $4;
    } else {
	die ("Get_Dust_Attenuated_Luminosity(): unable to parse dataset name");
    }

    # Construct a label for the filter.
    my $filterLabel;
    if ( $frame eq "rest" ) {
	$filterLabel = $filter.":rest";
    } else {
 	$filterLabel = $filter.":observed:z".$redshift;
    }

    # Ensure we have an effective wavelength and a wavelength index for this filter.
    unless ( exists($effectiveWavelength{$filterLabel}) ) {
	my $filterPath = "./data/filters/".$filter.".xml";
	die("Get_Dust_Attenuated_Luminosity(): can not find filter file for: ".$filter) unless ( -e $filterPath );
	my $xml = new XML::Simple;
	my $filterData = $xml->XMLin($filterPath);
	unless ( exists($filterData->{'effectiveWavelength'}) ) {
	    # No effective wavelength data available for filter - run the script that computes it.
	    system("scripts/filters/vega_offset_effective_lambda.pl");
	    $filterData = $xml->XMLin($filterPath);
	    die ("Get_Dust_Attenuated_Luminosity(): failed to compute effective wavelengths for filters") unless ( exists($filterData->{'effectiveWavelength'}) );
	}
	if ( $frame eq "rest" ) {
	    $effectiveWavelength{$filterLabel} = pdl $filterData->{'effectiveWavelength'};
	} else {
	    $effectiveWavelength{$filterLabel} = pdl $filterData->{'effectiveWavelength'}/(1.0+$redshift);
	}
	($wavelengthIndex{$filterLabel}, my $error) = interpolate($effectiveWavelength{$filterLabel},$wavelengths,$wavelengthIndices);
    }

    # Get interpolations of inclinations.
    my $galacticusInclination;
    if ( $faceOn == 1 ) {
	$galacticusInclination = zeroes(nelem($dataSets->{"diskScaleLength"}));
    } else {
	$galacticusInclination = $dataSets->{"inclination"};
    }
    my $inclinationIndex;
    if ( $component      eq "disk"     ) {
	($inclinationIndex, my $error) = interpolate($galacticusInclination,$diskInclinations    ,$diskInclinationIndices    );
    } elsif ( $component eq "spheroid" ) {
	($inclinationIndex, my $error) = interpolate($galacticusInclination,$spheroidInclinations,$spheroidInclinationIndices);
    } else {
	die("Get_Dust_Attenuated_Luminosity(): unknown component");
    }

    # Get interpolations of bulge sizes.
    my $sizeIndex;
    if ( $component eq "spheroid" ) {
	# Compute size as spheroid (assumed to be Hernquist profile) half-mass radius in units of disk scale length.
	my $sizes;
	if ( $dataSet->{'parameters'}->{"spheroidMassDistribution"} eq "Hernquist" ) {
	    $sizes = (1.0+sqrt(2.0))*$dataSets->{"spheroidScaleLength"}/$dataSets->{"diskScaleLength"};
	}
	if ( $dataSet->{'parameters'}->{"spheroidMassDistribution"} eq "Sersic" ) {
	    $sizes =                 $dataSets->{"spheroidScaleLength"}/$dataSets->{"diskScaleLength"};
	}
	my $sizesLimited;
	if ( $extrapolateInSize == 1 ) {
	    my $sizeMinimum  = $spheroidSizes->index(0);
	    my $sizeMaximum  = $spheroidSizes->index($spheroidSizesCount-1);
	    $sizesLimited = $sizeMaximum*($sizes>$sizeMaximum) + $sizeMinimum*($sizes<$sizeMinimum) + $sizes*(($sizes>=$sizeMinimum) & ($sizes<=$sizeMaximum));
	} else {
	    $sizesLimited = $sizes;
	}
	($sizeIndex, my $error) = interpolate($sizesLimited,$spheroidSizes,$spheroidSizeIndices);
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
    my $gasMetalMass                   = pdl $dataSets->{"diskGasMetals"};
    my $diskScaleLength                = pdl $dataSets->{"diskScaleLength"};
    my $noDisks                        = which($diskScaleLength <= 0.0);
    my $gasMetalsSurfaceDensityCentral = $gasMetalMass/(2.0*$Pi*($mega*$diskScaleLength)**2);
    $gasMetalsSurfaceDensityCentral->index($noDisks) .= 0.0;

    # Compute central optical depths.
    my $opticalDepthCentral            = $opticalDepthNormalization*$gasMetalsSurfaceDensityCentral;
    my $opticalDepthIndex;
    if ( $component      eq "disk"     ) {
	my $opticalDepthCentralLimited;
	if ( $extrapolateInTau == 1 ) {
	    my $tauMinimum                 = $diskOpticalDepths->index(0);
	    my $tauMaximum                 = $diskOpticalDepths->index($diskOpticalDepthsCount-1);
	    $opticalDepthCentralLimited = $tauMaximum*($opticalDepthCentral>$tauMaximum)
		+ $tauMinimum*($opticalDepthCentral<$tauMinimum)
		+ $opticalDepthCentral*(($opticalDepthCentral>=$tauMinimum) & ($opticalDepthCentral<=$tauMaximum));
	} else {
	    $opticalDepthCentralLimited = $opticalDepthCentral;
	}
	($opticalDepthIndex,my $error) = interpolate($opticalDepthCentralLimited,$diskOpticalDepths    ,$diskOpticalDepthsIndices    );
    } elsif ( $component eq "spheroid" ) {
	my $opticalDepthCentralLimited;
	if ( $extrapolateInTau == 1 ) {
	    my $tauMinimum                 = $spheroidOpticalDepths->index(0);
	    my $tauMaximum                 = $spheroidOpticalDepths->index($spheroidOpticalDepthsCount-1);
	    $opticalDepthCentralLimited = $tauMaximum*($opticalDepthCentral>$tauMaximum)
		+ $tauMinimum*($opticalDepthCentral<$tauMinimum)
		+ $opticalDepthCentral*(($opticalDepthCentral>=$tauMinimum) & ($opticalDepthCentral<=$tauMaximum));
	} else {
	    $opticalDepthCentralLimited = $opticalDepthCentral;
	}
	($opticalDepthIndex,my $error) = interpolate($opticalDepthCentralLimited,$spheroidOpticalDepths,$spheroidOpticalDepthsIndices);
    } else{
 	die("Get_Dust_Attenuated_Luminosity(): unknown component");
    }

    # Interpolate in the attenuation table.
    $PDL::BIGPDL = 1;
    my $attenuations;
    if ( $component eq "disk" ) {
	my $indices = zeroes(3,nelem($inclinationIndex));
	$indices->((0),0:nelem($inclinationIndex)-1) .= $inclinationIndex;
	$indices->((1),0:nelem($inclinationIndex)-1) .= $opticalDepthIndex;
	$indices->((2),0:nelem($inclinationIndex)-1) .= $wavelengthIndex{$filterLabel};
	$attenuations = $diskAttenuations->interpND($indices);
    } elsif ( $component eq "spheroid" ) {
	my $indices = zeroes(4,nelem($inclinationIndex));
	$indices->((0),0:nelem($inclinationIndex)-1) .= $sizeIndex;
	$indices->((1),0:nelem($inclinationIndex)-1) .= $inclinationIndex;
	$indices->((2),0:nelem($inclinationIndex)-1) .= $opticalDepthIndex;
	$indices->((3),0:nelem($inclinationIndex)-1) .= $wavelengthIndex{$filterLabel};
	$attenuations = $spheroidAttenuations->interpND($indices);
    } else{
 	die("Get_Dust_Attenuated_Luminosity(): unknown component");
    }
    $PDL::BIGPDL = 0;
    
    # Multiply luminosities by attenuations.
    $dataSets->{$dataSetName} = $dataSets->{$luminosityDataSet}*$attenuations;
}

sub Load_Dust_Atlas {
    my $dataSet = shift;
    # Ensure that dust attenuation data is loaded.
    unless ( $dustDataLoaded == 1 ) {
	# Specify the dust file.
	my $dustFile;
	if ( exists($dataSet->{'dustAtlasFile'}) ) {
	    $dustFile = $dataSet->{'dustAtlasFile'};
	} else {
	    $dustFile = "data/dust/atlasFerrara2000/attenuations_MilkyWay_dustHeightRatio1.0.xml";
	}
	
        # Read the dust file.
	my $xml = new XML::Simple;
	my $dustData = $xml->XMLin($dustFile);

        # Process the dust data into a PDL array.

        # Get wavelengths.
	$wavelengths = pdl @{$dustData->{'wavelengths'}->{'lambda'}};
	$wavelengthCount = nelem($wavelengths);
	$wavelengthIndices = pdl [0..$wavelengthCount-1];
	
        # Get attenuations.
	my $diskPdlCreated = 0;
	my $spheroidPdlCreated = 0;
	foreach my $component ( keys(%{$dustData->{'components'}}) ) {
	    if ( $component eq "disk" ) {
		#print "Processing dust data for: ".$component."\n";
		my $componentStruct = $dustData->{'components'}->{$component};
		$diskInclinationsCount = $#{$componentStruct->{'inclination'}}+1;
		$diskInclinations = zeroes($diskInclinationsCount);
		$diskInclinationIndices = pdl [0..$diskInclinationsCount-1];
		for(my $iInclination=0;$iInclination<=$#{$componentStruct->{'inclination'}};++$iInclination) {
		    my $inclinationStruct = ${$dustData->{'components'}->{$component}->{'inclination'}}[$iInclination];
		    $diskInclinations->(($iInclination)) += $inclinationStruct->{'angle'};
		    #print "  Processing inclination ".$inclinationStruct->{'angle'}."\n";
		    $diskOpticalDepthsCount = $#{$inclinationStruct->{'opticalDepth'}}+1;
		    $diskOpticalDepths = zeroes($diskOpticalDepthsCount);
		    for(my $iOpticalDepth=0;$iOpticalDepth<=$#{$inclinationStruct->{'opticalDepth'}};++$iOpticalDepth) {
			my $opticalDepthStruct = ${$inclinationStruct->{'opticalDepth'}}[$iOpticalDepth];
			$diskOpticalDepths->(($iOpticalDepth)) += $opticalDepthStruct->{'tau'};
			$diskOpticalDepthsIndices = pdl [0..$diskOpticalDepthsCount-1];
			#print "    Processing optical depth ".$opticalDepthStruct->{'tau'}."\n";
			unless ( $diskPdlCreated == 1 ) {
			    $diskAttenuations = zeroes($diskInclinationsCount,$diskOpticalDepthsCount,$wavelengthCount);
			    $diskPdlCreated = 1;
			}
			for(my $iAttenuation=0;$iAttenuation<=$#{$opticalDepthStruct->{'attenuation'}};++$iAttenuation) {
			    my $attenuation = ${$opticalDepthStruct->{'attenuation'}}[$iAttenuation];
			    #print "      Got attenuation of ".$attenuation." at ".$wavelengths->index($iAttenuation)." Angstroms\n";
			    $diskAttenuations->(($iInclination),($iOpticalDepth),($iAttenuation)) += $attenuation;
			}
		    }
		}
	    }
	    if ( $component eq "bulge" ) {
		#print "Processing dust data for: ".$component."\n";
		my $componentStruct = $dustData->{'components'}->{$component};
		$spheroidSizesCount = $#{$componentStruct->{'bulgeSize'}}+1;
		$spheroidSizes = zeroes($spheroidSizesCount);
		$spheroidSizeIndices = pdl [0..$spheroidSizesCount-1];
		for(my $iSize=0;$iSize<=$#{$componentStruct->{'bulgeSize'}};++$iSize) {
		    my $sizeStruct = ${$dustData->{'components'}->{$component}->{'bulgeSize'}}[$iSize];
		    $spheroidSizes->(($iSize)) += $sizeStruct->{'size'};
		    #print "  Processing size ".$sizeStruct->{'size'}."\n";
		    $spheroidInclinationsCount = $#{$sizeStruct->{'inclination'}}+1;
		    $spheroidInclinations = zeroes($spheroidInclinationsCount);
		    $spheroidInclinationIndices = pdl [0..$spheroidInclinationsCount-1];
		    for(my $iInclination=0;$iInclination<=$#{$sizeStruct->{'inclination'}};++$iInclination) {
			my $inclinationStruct = ${$sizeStruct->{'inclination'}}[$iInclination];
			$spheroidInclinations->(($iInclination)) += $inclinationStruct->{'angle'};
			#print "  Processing inclination ".$inclinationStruct->{'angle'}."\n";
			$spheroidOpticalDepthsCount = $#{$inclinationStruct->{'opticalDepth'}}+1;
			$spheroidOpticalDepths = zeroes($spheroidOpticalDepthsCount);
			for(my $iOpticalDepth=0;$iOpticalDepth<=$#{$inclinationStruct->{'opticalDepth'}};++$iOpticalDepth) {
			    my $opticalDepthStruct = ${$inclinationStruct->{'opticalDepth'}}[$iOpticalDepth];
			    $spheroidOpticalDepths->(($iOpticalDepth)) += $opticalDepthStruct->{'tau'};
			    $spheroidOpticalDepthsIndices = pdl [0..$spheroidOpticalDepthsCount-1];
			    #print "    Processing optical depth ".$opticalDepthStruct->{'tau'}."\n";
			    unless ( $spheroidPdlCreated == 1 ) {
				$spheroidAttenuations = zeroes($spheroidSizesCount,$spheroidInclinationsCount,$spheroidOpticalDepthsCount,$wavelengthCount);
				$spheroidPdlCreated = 1;
			    }
			    for(my $iAttenuation=0;$iAttenuation<=$#{$opticalDepthStruct->{'attenuation'}};++$iAttenuation) {
				my $attenuation = ${$opticalDepthStruct->{'attenuation'}}[$iAttenuation];
				#print "      Got attenuation of ".$attenuation." at ".$wavelengths->index($iAttenuation)." Angstroms\n";
				$spheroidAttenuations->(($iSize),($iInclination),($iOpticalDepth),($iAttenuation)) += $attenuation;
			    }
			}
		    }
		}
	    }
	}
	$dustDataLoaded = 1;
    }
}

1;
