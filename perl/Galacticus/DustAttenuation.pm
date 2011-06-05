# Contains a Perl module which implements dust attenuated luminosity calculations for Galacticus.

package DustAttenuation;
use PDL;
use PDL::NiceSlice;
use XML::Simple;
use Galacticus::HDF5;
use Galacticus::Inclination;

%HDF5::galacticusFunctions = ( %HDF5::galacticusFunctions,
    "^(disk|spheroid)StellarLuminosity:.*:dustAtlas(\\[faceOn\\])?\$" => "DustAttenuation::Get_Dust_Attenuated_Luminosity"
    );

# Flag indicating whether dust data is loaded yet.
$dustDataLoaded = 0;

my $status = 1;
$status;

sub Get_Dust_Attenuated_Luminosity {
    $dataSet = shift;
    $dataSetName = $_[0];

    # Check parameters.
    &HDF5::Get_Parameters($dataSet);
    die ("Get_Dust_Attenuated_Luminosity(): routine assumes exponential disks and Hernquist spheroids")
	unless ( ${$dataSet}{'parameters'}->{"treeNodeMethodDisk"} eq "exponential"
		 &&  ${$dataSet}{'parameters'}->{"treeNodeMethodSpheroid"} eq "Hernquist" );

    # Get the name of the unattenuated luminosity dataset.
    ($luminosityDataSet = $dataSetName) =~ s/:dustAtlas(\[faceOn\])?//;

    # List of properties to read.
    @propertyList = ("diskGasMetals","diskScaleLength","spheroidScaleLength",$luminosityDataSet);

    # Check if a face-on magnitude is required.
    if ( $dataSetName =~ m/\[faceOn\]/ ) {
	$faceOn = 1;
    } else {
	$faceOn = 0;
	push(@propertyList,"inclination");
    }

    # Get the datasets needed for our calculation.
    &HDF5::Get_Dataset($dataSet,\@propertyList);
    $dataSets = \%{${$dataSet}{'dataSets'}};

    # Load the dust atlas data.
    &Load_Dust_Atlas;

    # Determine types of extrapolation to be used.
    if ( exists(${$dataSet}{'dustAtlasExtrapolateInSize'}) ) {
	$extrapolateInSize = ${$dataSet}{'dustAtlasExtrapolateInSize'};
	if ( lc($extrapolateInSize) eq "yes" ) {$extrapolateInSize = 1};
	if ( lc($extrapolateInSize) eq "no"  ) {$extrapolateInSize = 0};
    } else {
	# Assume that we will extrapolate in spheroid sizes by default.
	$extrapolateInSize = 1;
    }
    if ( exists(${$dataSet}{'dustAtlasExtrapolateInTau'}) ) {
	$extrapolateInTau = ${$dataSet}{'dustAtlasExtrapolateInTau'};
	if ( lc($extrapolateInTau) eq "yes" ) {$extrapolateInTau = 1};
	if ( lc($extrapolateInTau) eq "no"  ) {$extrapolateInTau = 0};
    } else {
	# Assume that we will extrapolate in optical depth by default.
	$extrapolateInTau = 1;
    }

    # Extract filter data.
    if ( $dataSetName =~ m/^(.*?)StellarLuminosity:([^:]+):([^:]+):z([\d\.]+)/ ) {
	# Extract the dataset name information.
	$component = $1;
	$filter    = $2;
	$frame     = $3;
	$redshift  = $4;
    } else {
	die ("Get_Dust_Attenuated_Luminosity(): unable to parse dataset name");
    }

    # Construct a label for the filter.
    if ( $frame eq "rest" ) {
	$filterLabel = $filter.":rest";
    } else {
 	$filterLabel = $filter.":observed:z".$redshift;
    }

    # Ensure we have an effective wavelength and a wavelength index for this filter.
    unless ( exists($effectiveWavelength{$filterLabel}) ) {
	$filterPath = "./data/filters/".$filter.".xml";
	die("Get_Dust_Attenuated_Luminosity(): can not find filter file for: ".$filter) unless ( -e $filterPath );
	$xml = new XML::Simple;
	$filterData = $xml->XMLin($filterPath);
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
	($wavelengthIndex{$filterLabel},$error) = interpolate($effectiveWavelength{$filterLabel},$wavelengths,$wavelengthIndices);
    }

    # Get interpolations of inclinations.
    if ( $faceOn == 1 ) {
	$galacticusInclination = zeroes(nelem(${$dataSets->{"diskScaleLength"}}));
    } else {
	$galacticusInclination = ${$dataSets->{"inclination"}};
    }
    if ( $component      eq "disk"     ) {
	($inclinationIndex,$error) = interpolate($galacticusInclination,$diskInclinations    ,$diskInclinationIndices    );
    } elsif ( $component eq "spheroid" ) {
	($inclinationIndex,$error) = interpolate($galacticusInclination,$spheroidInclinations,$spheroidInclinationIndices);
    } else {
	die("Get_Dust_Attenuated_Luminosity(): unknown component");
    }

    # Get interpolations of bulge sizes.
    if ( $component eq "spheroid" ) {
	# Compute size as spheroid (assumed to be Hernquist profile) half-mass radius in units of disk scale length.
	$sizes = (1.0+sqrt(2.0))*${$dataSets->{"spheroidScaleLength"}}/${$dataSets->{"diskScaleLength"}};
	if ( $extrapolateInSizes == 1 ) {
	    $sizeMinimum  = $spheroidSizes->index(0);
	    $sizeMaximum  = $spheroidSizes->index($spheroidSizesCount-1);
	    $sizesLimited = $sizeMaximum*($sizes>$sizeMaximum) + $sizeMinimum*($sizes<$sizeMinimum) + $sizes*($sizes>=$sizeMinimum & $sizes<=$sizeMaximum);
	} else {
	    $sizesLimited = $sizes;
	}
	($sizeIndex,$error) = interpolate($sizesLimited,$spheroidSizes,$spheroidSizeIndices);
    }

    # Specify constants used in calculation of central optical depths.
    $localISMMetallicity       = pdl 0.02;                # Metallicity in the local ISM.
    $Pi                        = pdl 3.1415926536;        # Pi.
    $mega                      = pdl 1.0e6;               # SI prefix.
    $hecto                     = pdl 1.0e2;               # SI prefix.
    $parsec                    = pdl 3.0856775807e16;     # Parsec in m.
    $atomicMass                = pdl 1.66053873e-27;      # Atomic mass in kg.
    $solarMass                 = pdl 1.9891e30;           # Solar mass in kg.
    $AV_to_EBV                 = pdl 3.10;                # (A_V/E(B-V); Savage & Mathis 1979)
    $NH_to_EBV                 = pdl 5.8e21;              # (N_H/E(B-V); atoms/cm^2/mag; Savage & Mathis 1979)
    $hydrogenMassFractionSolar = pdl 0.707;               # Hydrogen mass fraction in the Sun.
    $opticalDepthToMagnitudes  = pdl 2.5*log10(exp(1.0)); # Conversion factor from optical depth to magnitudes of extinction.
    $opticalDepthNormalization = (1.0/$opticalDepthToMagnitudes)*($AV_to_EBV/$NH_to_EBV)*($hydrogenMassFractionSolar/$atomicMass)*($solarMass/($parsec*$hecto)**2)/$localISMMetallicity;

    # Central surface density in M_Solar/pc^2.
    $gasMetalMass                   = pdl ${$dataSets->{"diskGasMetals"}};
    $diskScaleLength                = pdl ${$dataSets->{"diskScaleLength"}};
    $noDisks                        = which($diskScaleLength <= 0.0);
    $gasMetalsSurfaceDensityCentral = $gasMetalMass/(2.0*$Pi*($mega*$diskScaleLength)**2);
    $gasMetalsSurfaceDensityCentral->index($noDisks) .= 0.0;

    # Compute central optical depths.
    $opticalDepthCentral            = $opticalDepthNormalization*$gasMetalsSurfaceDensityCentral;
    if ( $component      eq "disk"     ) {
	if ( $extrapolateInTau == 1 ) {
	    $tauMinimum                 = $diskOpticalDepths->index(0);
	    $tauMaximum                 = $diskOpticalDepths->index($diskOpticalDepthsCount-1);
	    $opticalDepthCentralLimited = $tauMaximum*($opticalDepthCentral>$tauMaximum)
		+ $tauMinimum*($opticalDepthCentral<$tauMinimum)
		+ $opticalDepthCentral*($opticalDepthCentral>=$tauMinimum & $opticalDepthCentral<=$tauMaximum);
	} else {
	    $opticalDepthCentralLimited = $opticalDepthCentral;
	}
	($opticalDepthIndex,$error) = interpolate($opticalDepthCentralLimited,$diskOpticalDepths    ,$diskOpticalDepthsIndices    );
    } elsif ( $component eq "spheroid" ) {
	if ( $extrapolateInTau == 1 ) {
	    $tauMinimum                 = $spheroidOpticalDepths->index(0);
	    $tauMaximum                 = $spheroidOpticalDepths->index($spheroidOpticalDepthsCount-1);
	    $opticalDepthCentralLimited = $tauMaximum*($opticalDepthCentral>$tauMaximum)
		+ $tauMinimum*($opticalDepthCentral<$tauMinimum)
		+ $opticalDepthCentral*($opticalDepthCentral>=$tauMinimum & $opticalDepthCentral<=$tauMaximum);
	} else {
	    $opticalDepthCentralLimited = $opticalDepthCentral;
	}
	($opticalDepthIndex,$error) = interpolate($opticalDepthCentralLimited,$spheroidOpticalDepths,$spheroidOpticalDepthsIndices);
    } else{
 	die("Get_Dust_Attenuated_Luminosity(): unknown component");
    }

    # Interpolate in the attenuation table.
    if ( $component eq "disk" ) {
	$indices = zeroes(3,nelem($inclinationIndex));
	$indices->((0),0:nelem($inclinationIndex)-1) .= $inclinationIndex;
	$indices->((1),0:nelem($inclinationIndex)-1) .= $opticalDepthIndex;
	$indices->((2),0:nelem($inclinationIndex)-1) .= $wavelengthIndex{$filterLabel};
	$attenuations = $diskAttenuations->interpND($indices);
    } elsif ( $component eq "spheroid" ) {
	$indices = zeroes(4,nelem($inclinationIndex));
	$indices->((0),0:nelem($inclinationIndex)-1) .= $sizeIndex;
	$indices->((1),0:nelem($inclinationIndex)-1) .= $inclinationIndex;
	$indices->((2),0:nelem($inclinationIndex)-1) .= $opticalDepthIndex;
	$indices->((3),0:nelem($inclinationIndex)-1) .= $wavelengthIndex{$filterLabel};
	$attenuations = $spheroidAttenuations->interpND($indices);
    } else{
 	die("Get_Dust_Attenuated_Luminosity(): unknown component");
    }
    
    # Multiply luminosities by attenuations.
    ${$dataSets->{$dataSetName}} = ${$dataSets->{$luminosityDataSet}}*$attenuations;

}

sub Load_Dust_Atlas {
    # Ensure that dust attenuation data is loaded.
    unless ( $dustDataLoaded == 1 ) {
	# Specify the dust file.
	if ( exists(${$dataSet}{'dustAtlasFile'}) ) {
	    $dustFile = ${$dataSet}{'dustAtlasFile'};
	} else {
	    $dustFile = "data/dust_atlas/attenuations_MilkyWay_dustHeightRatio1.0.xml";
	}
	
        # Read the dust file.
	$xml = new XML::Simple;
	$dustData = $xml->XMLin($dustFile);

        # Process the dust data into a PDL array.

        # Get wavelengths.
	$wavelengths = pdl @{$dustData->{'wavelengths'}->{'lambda'}};
	$wavelengthCount = nelem($wavelengths);
	$wavelengthIndices = pdl [0..$wavelengthCount-1];
	
        # Get attenuations.
	$diskPdlCreated = 0;
	$spheroidPdlCreated = 0;
	foreach $component ( keys(%{$dustData->{'components'}}) ) {
	    if ( $component eq "disk" ) {
		#print "Processing dust data for: ".$component."\n";
		$componentStruct = $dustData->{'components'}->{$component};
		$diskInclinationsCount = $#{$componentStruct->{'inclination'}}+1;
		$diskInclinations = zeroes($diskInclinationsCount);
		$diskInclinationIndices = pdl [0..$diskInclinationsCount-1];
		for($iInclination=0;$iInclination<=$#{$componentStruct->{'inclination'}};++$iInclination) {
		    $inclinationStruct = ${$dustData->{'components'}->{$component}->{'inclination'}}[$iInclination];
		    $diskInclinations->(($iInclination)) += $inclinationStruct->{'angle'};
		    #print "  Processing inclination ".$inclinationStruct->{'angle'}."\n";
		    $diskOpticalDepthsCount = $#{$inclinationStruct->{'opticalDepth'}}+1;
		    $diskOpticalDepths = zeroes($diskOpticalDepthsCount);
		    for($iOpticalDepth=0;$iOpticalDepth<=$#{$inclinationStruct->{'opticalDepth'}};++$iOpticalDepth) {
			$opticalDepthStruct = ${$inclinationStruct->{'opticalDepth'}}[$iOpticalDepth];
			$diskOpticalDepths->(($iOpticalDepth)) += $opticalDepthStruct->{'tau'};
			$diskOpticalDepthsIndices = pdl [0..$diskOpticalDepthsCount-1];
			#print "    Processing optical depth ".$opticalDepthStruct->{'tau'}."\n";
			unless ( $diskPdlCreated == 1 ) {
			    $diskAttenuations = zeroes($diskInclinationsCount,$diskOpticalDepthsCount,$wavelengthCount);
			    $diskPdlCreated = 1;
			}
			for($iAttenuation=0;$iAttenuation<=$#{$opticalDepthStruct->{'attenuation'}};++$iAttenuation) {
			    $attenuation = ${$opticalDepthStruct->{'attenuation'}}[$iAttenuation];
			    #print "      Got attenuation of ".$attenuation." at ".$wavelengths->index($iAttenuation)." Angstroms\n";
			    $diskAttenuations->(($iInclination),($iOpticalDepth),($iAttenuation)) += $attenuation;
			}
		    }
		}
	    }
	    if ( $component eq "bulge" ) {
		#print "Processing dust data for: ".$component."\n";
		$componentStruct = $dustData->{'components'}->{$component};
		$spheroidSizesCount = $#{$componentStruct->{'bulgeSize'}}+1;
		$spheroidSizes = zeroes($spheroidSizesCount);
		$spheroidSizeIndices = pdl [0..$spheroidSizesCount-1];
		for($iSize=0;$iSize<=$#{$componentStruct->{'bulgeSize'}};++$iSize) {
		    $sizeStruct = ${$dustData->{'components'}->{$component}->{'bulgeSize'}}[$iSize];
		    $spheroidSizes->(($iSize)) += $sizeStruct->{'size'};
		    #print "  Processing size ".$sizeStruct->{'size'}."\n";
		    $spheroidInclinationsCount = $#{$sizeStruct->{'inclination'}}+1;
		    $spheroidInclinations = zeroes($spheroidInclinationsCount);
		    $spheroidInclinationIndices = pdl [0..$spheroidInclinationsCount-1];
		    for($iInclination=0;$iInclination<=$#{$sizeStruct->{'inclination'}};++$iInclination) {
			$inclinationStruct = ${$sizeStruct->{'inclination'}}[$iInclination];
			$spheroidInclinations->(($iInclination)) += $inclinationStruct->{'angle'};
			#print "  Processing inclination ".$inclinationStruct->{'angle'}."\n";
			$spheroidOpticalDepthsCount = $#{$inclinationStruct->{'opticalDepth'}}+1;
			$spheroidOpticalDepths = zeroes($spheroidOpticalDepthsCount);
			for($iOpticalDepth=0;$iOpticalDepth<=$#{$inclinationStruct->{'opticalDepth'}};++$iOpticalDepth) {
			    $opticalDepthStruct = ${$inclinationStruct->{'opticalDepth'}}[$iOpticalDepth];
			    $spheroidOpticalDepths->(($iOpticalDepth)) += $opticalDepthStruct->{'tau'};
			    $spheroidOpticalDepthsIndices = pdl [0..$spheroidOpticalDepthsCount-1];
			    #print "    Processing optical depth ".$opticalDepthStruct->{'tau'}."\n";
			    unless ( $spheroidPdlCreated == 1 ) {
				$spheroidAttenuations = zeroes($spheroidSizesCount,$spheroidInclinationsCount,$spheroidOpticalDepthsCount,$wavelengthCount);
				$spheroidPdlCreated = 1;
			    }
			    for($iAttenuation=0;$iAttenuation<=$#{$opticalDepthStruct->{'attenuation'}};++$iAttenuation) {
				$attenuation = ${$opticalDepthStruct->{'attenuation'}}[$iAttenuation];
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
