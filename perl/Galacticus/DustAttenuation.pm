# Contains a Perl module which implements dust attenuated luminosity calculations for Galacticus.

package DustAttenuation;
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V094"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V094"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
}
unshift(@INC,$galacticusPath."perl"); 
use PDL;
use PDL::NiceSlice;
use XML::Simple;
use Data::Dumper;
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
	 (
	  ! exists($dataSet ->{'parameters'}->{"diskMassDistribution"    })
	  ||
	           $dataSet ->{'parameters'}->{"diskMassDistribution"    }->{'value'} eq "exponentialDisk"
    	 )
	 &&
    	 (
	  ! exists($dataSet->{'parameters'}->{"spheroidMassDistribution"})
	  ||
	           $dataSet->{'parameters'}->{"spheroidMassDistribution"}->{'value'} eq "hernquist"
	  ||
	           $dataSet->{'parameters'}->{"spheroidMassDistribution"}->{'value'} eq "sersic" 
    	 )
    	);

    # Get the name of the unattenuated luminosity dataset.
    (my $luminosityDataSet = $dataSetName) =~ s/:dustAtlas(\[faceOn\])?//;

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

    # List of properties to read.
    my @propertyList = ("diskAbundancesGasMetals","diskRadius",$luminosityDataSet);
    push(@propertyList,"spheroidRadius")
	if ( $component eq "spheroid" );
    
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

    # Ensure we have an effective wavelength and a wavelength index for this filter.
    unless ( exists($effectiveWavelength{$filterLabel}) ) {
	my $filterPath = $galacticusPath."data/filters/".$filter.".xml";
	die("Get_Dust_Attenuated_Luminosity(): can not find filter file for: ".$filter) unless ( -e $filterPath );
	my $xml = new XML::Simple;
	my $filterData = $xml->XMLin($filterPath);
	unless ( exists($filterData->{'effectiveWavelength'}) ) {
	    # No effective wavelength data available for filter - run the script that computes it.
	    system($galacticusPath."scripts/filters/vega_offset_effective_lambda.pl");
	    $filterData = $xml->XMLin($filterPath);
	    die ("Get_Dust_Attenuated_Luminosity(): failed to compute effective wavelengths for filters") unless ( exists($filterData->{'effectiveWavelength'}) );
	}
	if ( $frame eq "rest" ) {
	    $effectiveWavelength{$filterLabel} = pdl $filterData->{'effectiveWavelength'};
	} else {
	    $effectiveWavelength{$filterLabel} = pdl $filterData->{'effectiveWavelength'}/(1.0+$redshift);
	}
	($wavelengthIndex{$filterLabel}, my $error) = interpolate(log($effectiveWavelength{$filterLabel}),log($wavelengths),$wavelengthIndices);    
    }
    # Get interpolations of inclinations.
    my $galacticusInclination;
    if ( $faceOn == 1 ) {
	$galacticusInclination = zeroes(nelem($dataSets->{"diskRadius"}));
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
    my $diskless;
    if ( $component eq "spheroid" ) {
	# Compute size as spheroid (assumed to be Hernquist profile) half-mass radius in units of disk scale length.
	my $sizes;
	if ( $dataSet->{'parameters'}->{"spheroidMassDistribution"}->{'value'} eq "hernquist" ) {
	    $sizes = (1.0+sqrt(2.0))*$dataSets->{"spheroidRadius"}/$dataSets->{"diskRadius"};
	}
	if ( $dataSet->{'parameters'}->{"spheroidMassDistribution"}->{'value'} eq "sersic" ) {
	    $sizes =                 $dataSets->{"spheroidRadius"}/$dataSets->{"diskRadius"};
	}
	# Identify diskless galaxies and assign an arbitrary size. These will later have attenuation set to unity anyway, so the
	# value here does not matter.
	$diskless = which($dataSets->{'diskRadius'} <= 0.0);
	$sizes->($diskless) .= 1.0
	    if ( nelem($diskless) > 0 );
	my $sizesLimited;
	if ( $extrapolateInSize == 0 ) {
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
    my $gasMetalMass                   = pdl $dataSets->{"diskAbundancesGasMetals"};
    my $diskScaleLength                = pdl $dataSets->{"diskRadius"             };
    my $noDisks                        = which($diskScaleLength <= 0.0);
    my $gasMetalsSurfaceDensityCentral = $gasMetalMass/(2.0*$Pi*($mega*$diskScaleLength)**2);
    $gasMetalsSurfaceDensityCentral->index($noDisks) .= 0.0;
    # Compute central optical depths.
    my $opticalDepthCentral            = $opticalDepthNormalization*$gasMetalsSurfaceDensityCentral;
    my $opticalDepthIndex;
    if ( $component      eq "disk"     ) {
	my $opticalDepthCentralLimited;
	if ( $extrapolateInTau == 0 ) {
	    my $tauMinimum              = $diskOpticalDepths->index(0);
	    my $tauMaximum              = $diskOpticalDepths->index($diskOpticalDepthsCount-1);
	    $opticalDepthCentralLimited = $tauMaximum*($opticalDepthCentral>$tauMaximum)
		+ $tauMinimum*($opticalDepthCentral<$tauMinimum)
		+ $opticalDepthCentral*(($opticalDepthCentral>=$tauMinimum) & ($opticalDepthCentral<=$tauMaximum));
	} else {
	    $opticalDepthCentralLimited = $opticalDepthCentral;
	}
	($opticalDepthIndex,my $error) = interpolate($opticalDepthCentralLimited,$diskOpticalDepths    ,$diskOpticalDepthsIndices    );
    } elsif ( $component eq "spheroid" ) {
	my $opticalDepthCentralLimited;
	if ( $extrapolateInTau == 0 ) {
	    my $tauMinimum              = $spheroidOpticalDepths->index(0);
	    my $tauMaximum              = $spheroidOpticalDepths->index($spheroidOpticalDepthsCount-1);
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
	$attenuations = &dustTableInterpolation($diskAttenuations,$indices);
    } elsif ( $component eq "spheroid" ) {
	my $indices = zeroes(4,nelem($inclinationIndex));
	$indices->((0),0:nelem($inclinationIndex)-1) .= $sizeIndex;
	$indices->((1),0:nelem($inclinationIndex)-1) .= $inclinationIndex;
	$indices->((2),0:nelem($inclinationIndex)-1) .= $opticalDepthIndex;
	$indices->((3),0:nelem($inclinationIndex)-1) .= $wavelengthIndex{$filterLabel};
	$attenuations = &dustTableInterpolation($spheroidAttenuations,$indices);
    } else{
 	die("Get_Dust_Attenuated_Luminosity(): unknown component");
    }
    $attenuations->($diskless) .= 1.0
	if ( nelem($diskless) > 0 );
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
	    $dustFile = $galacticusPath."data/dust/atlasFerrara2000/attenuations_MilkyWay_dustHeightRatio1.0.xml";
	}
	
        # Read the dust file.
	my $xml = new XML::Simple;
	my $dustData = $xml->XMLin($dustFile);

        # Process the dust data into a PDL array.

        # Get wavelengths.
	$wavelengths       = pdl @{$dustData->{'wavelengths'}->{'lambda'}};
	$wavelengthCount   = nelem($wavelengths);
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
		    $diskOpticalDepthsCount   = $#{$inclinationStruct->{'opticalDepth'}}+2; # Add an extra zero optical depth.
		    $diskOpticalDepths        = zeroes($diskOpticalDepthsCount);
		    $diskOpticalDepthsIndices = pdl [0..$diskOpticalDepthsCount-1];
		    $diskOpticalDepths->((0)) .= 0.0;
		    for(my $iOpticalDepth=1;$iOpticalDepth<$diskOpticalDepthsCount;++$iOpticalDepth) {
			my $opticalDepthStruct = ${$inclinationStruct->{'opticalDepth'}}[$iOpticalDepth-1];
			$diskOpticalDepths->(($iOpticalDepth)) += $opticalDepthStruct->{'tau'};
			#print "    Processing optical depth ".$opticalDepthStruct->{'tau'}."\n";
			unless ( $diskPdlCreated == 1 ) {
			    $diskAttenuations = zeroes($diskInclinationsCount,$diskOpticalDepthsCount,$wavelengthCount);
			    $diskAttenuations->(:,(0),:) .= 1.0;
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
			$spheroidOpticalDepthsCount   = $#{$inclinationStruct->{'opticalDepth'}}+2; # Add an extra zero optical depth.
			$spheroidOpticalDepths        = zeroes($spheroidOpticalDepthsCount);
			$spheroidOpticalDepthsIndices = pdl [0..$spheroidOpticalDepthsCount-1];
			for(my $iOpticalDepth=1;$iOpticalDepth<$spheroidOpticalDepthsCount;++$iOpticalDepth) {
			    my $opticalDepthStruct = ${$inclinationStruct->{'opticalDepth'}}[$iOpticalDepth-1];
			    $spheroidOpticalDepths->(($iOpticalDepth)) += $opticalDepthStruct->{'tau'};
			    #print "    Processing optical depth ".$opticalDepthStruct->{'tau'}."\n";
			    unless ( $spheroidPdlCreated == 1 ) {
				$spheroidAttenuations = zeroes($spheroidSizesCount,$spheroidInclinationsCount,$spheroidOpticalDepthsCount,$wavelengthCount);
				$spheroidAttenuations->(:,:,(0),:) .= 1.0;
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

sub dustTableInterpolation {
    # Interpolate (with extrapolation) in dust tables,
    my $dustTable            = shift();
    my $interpolationIndices = shift();
    # Iterate over dimensions, computing interpolation/extrapolation factors.
    my $factors              = pdl zeroes(2,$interpolationIndices->dim(0),$interpolationIndices->dim(1));
    my $indexInRange         = pdl zeroes(  $interpolationIndices->dim(0),$interpolationIndices->dim(1));
    for(my $i=0;$i<$interpolationIndices->dim(0);++$i) {
	# Find allowed range of indices in tabulation.
	my $indexMinimum =                     0;
	my $indexMaximum = $dustTable->dim($i)-1;
	# Construct a set of in-range indices.
	$indexInRange->(($i),:) .= floor($interpolationIndices->(($i),:));
	my $belowRange   = which($indexInRange->(($i),:) <  $indexMinimum);
	my $aboveRange   = which($indexInRange->(($i),:) >= $indexMaximum);
	$indexInRange->(($i),$belowRange) .= $indexMinimum  
	    if ( nelem($belowRange) > 0 );
	$indexInRange->(($i),$aboveRange) .= $indexMaximum-1
	    if ( nelem($aboveRange) > 0 );
	# Compute interpolation factors.
	$factors->((1),($i),:) .= $interpolationIndices->(($i),:)-$indexInRange->(($i),:);
	$factors->((0),($i),:) .= 1.0-$factors->((1),($i),:);
    }
    # Perform the interpolation.
    my $attenuationSmall   = pdl 1.0e-30;
    my $attenuations       = pdl zeroes($interpolationIndices->dim(1))                       ;
    my $attenuationsMaxima = pdl   ones($interpolationIndices->dim(1))*log($attenuationSmall);
    if ( $interpolationIndices->dim(0) == 3 ) {
	for(my $j=0;$j<2;++$j) {
	    for(my $k=0;$k<2;++$k) {
		for(my $l=0;$l<2;++$l) {
		    # Shift the indices to use in the dust tables according to out interpolation hyper-cube.
		    my $indexInRangeShift    = $indexInRange+[$j,$k,$l];
		    # Compute the attenuation at these indices. We use the logarithmic attenuation as that is what we want to interpolate in.
		    my $attenuationTabulated = log($dustTable->indexND($indexInRangeShift)+$attenuationSmall);
		    $attenuations += 
			+$attenuationTabulated
			*$factors->(($j),(0),:)
			*$factors->(($k),(1),:)
			*$factors->(($l),(2),:);
		    # Find the maximum attenuation.
		    my $attenuationsMaximaExceeded = which($attenuationTabulated > $attenuationsMaxima);
		    $attenuationsMaxima->($attenuationsMaximaExceeded) .= 
			+$attenuationTabulated->($attenuationsMaximaExceeded)
			    if ( nelem($attenuationsMaximaExceeded) > 0 );
		}
	    }
	}
    } else {
	for(my $j=0;$j<2;++$j) {
	    for(my $k=0;$k<2;++$k) {
		for(my $l=0;$l<2;++$l) {
		    for(my $m=0;$m<2;++$m) {			
			my $indexInRangeShift    = $indexInRange+[$j,$k,$l,$m];
			my $attenuationTabulated = log($dustTable->indexND($indexInRangeShift)+$attenuationSmall);
			$attenuations += 
			    +$attenuationTabulated
			    *$factors->(($j),(0),:)
			    *$factors->(($k),(1),:)
			    *$factors->(($l),(2),:)
			    *$factors->(($m),(3),:);		    
			# Find the maximum attenuation.
			my $attenuationsMaximaExceeded = which($attenuationTabulated > $attenuationsMaxima);
			$attenuationsMaxima->($attenuationsMaximaExceeded) .= 
			    +$attenuationTabulated->($attenuationsMaximaExceeded)
			    if ( nelem($attenuationsMaximaExceeded) > 0 );
		    }
		}
	    }
	}
    }
    # Limit attenuations to be physical. Attenuations can not be negative - this is automatically handled by the fact that we
    # interpolate in the logarithm of the attenuation. They are also prevented from exceeding the maximum of 1 and the largest
    # tabulated attenuation used in their interpolation.
    my $attenuationAbove = which($attenuations > $attenuationsMaxima);
    $attenuations->flat()->($attenuationAbove) .= $attenuationsMaxima->($attenuationAbove)
	if ( nelem($attenuationAbove) > 0 );
    # Convert attenuations back from logarithmic form.
    $attenuations .= exp($attenuations);
    # Return the computed attenuations.
    return $attenuations;
}

1;
