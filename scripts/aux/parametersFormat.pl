#!/usr/bin/env perl
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
use Scalar::Util 'reftype';
use XML::LibXML qw(:libxml);
use XML::LibXML::PrettyPrint;
use Data::Dumper;
require Galacticus::Options;
require List::ExtraUtils;

# Format a Galacticus parameter file to make it more easily comprehansible.
# Andrew Benson (08-March-2016)

# Get arguments.
die("Usage: parametersFormat.pl <inputFile> <outputFile> [options]")
    unless ( scalar(@ARGV) >= 2 );
my $inputFileName  = $ARGV[0];
my $outputFileName = $ARGV[1];
my %options =
    (
    );
# Parse options.
%options = &Options::Parse_Options(\@ARGV,\%options);

# Parse the input file.
my $parser = XML::LibXML->new();
my $input  = $parser->parse_file($inputFileName);
my @parameterSets;
if ( $input->findnodes('parameters') ) {
    @parameterSets = $input->findnodes('parameters');
} else {
    die('parametersFormat.pl: can not find parameters')
}

# Define groups.
my @groups =
    (
     {
	 name        => "logging",
	 description => "Logging",
	 members     =>
	     [
	      "verbosityLevel"
	     ]
     },
     {
	 name        => "cosmology"         ,
	 description => "Cosmological model",
	 members     =>
	     [
	      "cosmologyParametersMethod",
	      "cosmologyFunctionsMethod",
	      qr/cosmology\d/
	     ]
     },
     {
	 name        => "powerSpectrum" ,
	 description => "Power spectrum",
	 members     =>
	     [
	      "powerSpectrumMethod",
	      "powerSpectrumPrimordialMethod",
	      "transferFunctionMethod",
	      "powerSpectrumPrimordialTransferredMethod",
	      "cosmologicalMassVarianceMethod",
	      "powerSpectrumNonlinearMethod"
	     ]
     },
     {
	 name        => "structureFormation" ,
	 description => "Structure formation",
	 members     =>
	     [
	      "linearGrowthMethod",
	      "criticalOverdensityMethod",
	      "virialDensityContrastMethod",
	      qr/virialDensityContrast.*/,
	      "haloMassFunctionMethod",
	      qr/haloMassFunctionSystematic\d/,
	      "gravitationalLensingMethod"
	     ]
     },
     {
	 name        => "igmReionization" ,
	 description => "Intergalactic medium and reionization",
	 members     =>
	     [
	      "reionizationSuppressionOpticalDepth",
	      "reionizationSuppressionVelocity",
	      "reionizationSuppressionRedshift",
	      "intergalacticMediumStateMethod",
	      qr/intergalaticMediumState.*/
	     ]
     },
     {
	 name        => "treeNodeMethods" ,
	 description => "Tree node component selection",
	 members     =>
	     [
	      qr/^treeNodeMethod/,
	      "nodeComponentBasicExtendedSphericalCollapseType"
	     ]
     },
     {
	 name        => "mergerTreeConstruction" ,
	 description => "Merger tree construction",
	 members     =>
	     [
	      "mergerTreeConstructMethod",
	      "mergerTreeBuilderMethod",
	      "mergerTreeMassResolutionMethod",
	      qr/^mergerTreeBuild.*/,
	      qr/^mergerTreesBuild.*/,
	      "treeBranchingMethod",
	      qr/^modifiedPressSchechter.*/,
	      "haloMassFunctionSamplingMethod",
	      qr/^haloMassFunctionSampling.*/,
	      qr/^mergerTreeRead.*/,
	      qr/mergerTreeImport.*/
	     ]
     },
     {
	 name        => "mergerTreeOperators" ,
	 description => "Merger tree operators",
	 members     =>
	     [
	      "mergerTreePruneBaryons"
	     ]
     },
     {
	 name        => "treeEvolution" ,
	 description => "Merger tree evolution flow control",
	 members     =>
	     [
	      "mergerTreesBuildFixedThreadAssignment",
	      "treeEvolveThreadLock",
	      "treeEvolveThreadLockName",
	      "treeEvolveThreadsMaximum"
	     ]
     },
     {
	 name        => "haloHierarchy" ,
	 description => "Merger tree halo hierarchy",
	 members     =>
	     [
	      "nodeMergersMethod",
	      "nodePromotionIndexShift"
	     ]
     },
     {
	 name        => "darkMatterHalo" ,
	 description => "Dark matter halo structure",
	 members     =>
	     [
	      "darkMatterProfileMethod",
	      "darkMatterProfileConcentrationMethod",
	      "darkMatterProfileMinimumConcentration",
	      "darkMatterProfileScaleCorrectForConcentrationDefinition"
	     ]
     },
     {
	 name        => "haloSpin" ,
	 description => "Dark matter halo spin",
	 members     =>
	     [
	      "haloSpinDistributionMethod",
	      "randomSpinResetMassFactor"
	     ]
     },
     {
	 name        => "satelliteOrbits" ,
	 description => "Satellite orbits",
	 members     =>
	     [
	      "virialOrbitMethod",
	      qr/^virialOrbits.*/,
	      "satelliteMergingTimescalesMethod",
	      "mergingTimescaleMultiplier",
	      "satellitesTidalFieldMethod",
	      "satelliteOrbitStoreOrbitalParameters"
	     ]
     },
     {
	 name        => "hotAtmosphere" ,
	 description => "Hot atmosphere",
	 members     =>
	     [
	      "hotHaloOutflowReincorporationMethod",
	      qr/^hotHaloOutflowReincorporation.*/,
	      "hotHaloTemperatureProfileMethod",
	      "hotHaloTrackStrippedGas",
	      "hotHaloOutflowReturnRate",
	      "hotHaloOutflowToColdMode",
	      "hotHaloAngularMomentumAlwaysGrows",
	      "hotHaloAngularMomentumLossFraction",
	      "hotHaloColdModeCoredIsothermalCoreRadiiMethod",
	      "hotHaloCoreRadiusOverVirialRadius",
	      "hotHaloExcessHeatDrivesOutflow",
	      "hotHaloMassDistributionMethod",
	      "hotHaloNodeMergerLimitBaryonFraction"
	     ]
     },
     {
	 name        => "coldAtmosphere" ,
	 description => "Cold atmosphere",
	 members     =>
	     [
	      "coldModeIsothermalCoreRadiusOverVirialRadius",
	      "coldModeMassDistribution"
	     ]
     },
     {
	 name        => "accretionFromIGM" ,
	 description => "Accretion from the IGM",
	 members     =>
	     [
	      "accretionHaloMethod",
	      "accretionHalosSimpleNegativeAccretionAllowed",
	      "accretionColdModeShockStabilityThreshold",
	      "accretionColdModeShockStabilityTransitionWidth"
	     ]
     },
     {
	 name        => "hotHaloRamPressure" ,
	 description => "Ram pressure stripping of hot atmosphere",
	 members     =>
	     [
	      "starveSatellites",
	      "hotHaloRamPressureForceMethod",
	      "hotHaloRamPressureStrippingMethod",
	      "hotHaloRamPressureStrippingTimescaleMethod",
	      "hotHaloOutflowStrippingEfficiency",
	      "ramPressureStrippingFormFactor"
	     ]
     },
     {
	 name        => "coolingAndInfall" ,
	 description => "Cooling and infall",
	 members     =>
	     [
	      "coolingRateMethod",
	      qr/^coolingRate.*/,
	      "coolingFunctionMethod",
	      "coolingMeanAngularMomentumFrom",
	      "coolingRadiusMethod",
	      "coolingRotationVelocityFrom",
	      "coolingSpecificAngularMomentumMethod",
	      "coolingTimeMethod",
	      qr/^coolingTime.*/,
	      "coldModeInfallRateMethod",
	      qr/^coldModeInfallRate.*/,
	      "zeroCoolingRateAboveVelocity"
	     ]
     },
     {
	 name        => "galacticStructure" ,
	 description => "Galactic structure",
	 members     =>
	     [
	      "diskMassDistribution",
	      "spheroidMassDistribution",
	      "galacticStructureRadiusSolverMethod",
	      qr/^galacticStructureRadii.*/,
	      "adiabaticContractionGnedinA",
	      "adiabaticContractionGnedinOmega",
	      "spheroidAngularMomentumAtScaleRadius"
	     ]
     },
     {
   	 name        => "barInstability" ,
	 description => "Galactic disk bar instability",
	 members     =>
	     [
	      "barInstabilityMethod",
	      "stabilityThresholdGaseous",
	      "stabilityThresholdStellar"
	     ]
     },
     {
	 name        => "diskRamPressure" ,
	 description => "Ram pressure stripping of galactic disk",
	 members     =>
	     [
	      "ramPressureStrippingMassLossRateDisksMethod"
	     ]
     },
     {
	 name        => "spheroidRamPressure" ,
	 description => "Ram pressure stripping of galactic spheroid",
	 members     =>
	     [
	      "ramPressureStrippingMassLossRateSpheroidsMethod"
	     ]
     },
     {
	 name        => "diskTidal" ,
	 description => "Tidal stripping of galactic disk",
	 members     =>
	     [
	      "tidalStrippingMassLossRateDisksMethod"
	     ]
     },
     {
	 name        => "spheroidTidal" ,
	 description => "Tidal stripping of galactic spheroid",
	 members     =>
	     [
	      "tidalStrippingMassLossRateSpheroidsMethod"
	     ]
     },
     {
	 name        => "diskStarFormation" ,
	 description => "Star formation in disks",
	 members     =>
	     [
	      "starFormationTimescaleDisksMethod",
	      qr/^starFormationTimescaleDisks.*/,
	      qr/^diskStarFormation.*/,
	      qr/^diskVerySimpleSurfaceDensity.*/,
	      "starFormationDiskMinimumTimescale",
	      "starFormationFrequencyKMT09",
	      "molecularComplexClumpingFactorKMT09",
	      "molecularFractionFastKMT09",
	      "starFormationRateSurfaceDensityDisksMethod",
	      "starFormationTimescaleIntegratedSurfaceDensityTolerance"
	     ]
     },
     {
	 name        => "spheroidStarFormation" ,
	 description => "Star formation in spheroids",
	 members     =>
	     [
	      "starFormationTimescaleSpheroidsMethod",
	      qr/^starFormationTimescaleSpheroids.*/,
	      qr/^spheroidStarFormation.*/,
	      "starFormationSpheroidEfficiency",
	      "starFormationSpheroidMinimumTimescale",
	      "starFormationSpheroidVelocityExponent",
	     ]
     },
     {
	 name        => "diskStellarFeedback" ,
	 description => "Stellar feedback in disks",
	 members     =>
	     [
	      "starFormationFeedbackDisksMethod",
	      qr/^diskOutflow.*/
	     ]
     },
     {
	 name        => "spheroidStellarFeedback" ,
	 description => "Stellar feedback in spheroids",
	 members     =>
	     [
	      "starFormationFeedbackSpheroidsMethod",
	      qr/^spheroidOutflow.*/
	     ]
     },
     {
	 name        => "initialMassFunction" ,
	 description => "Stellar initial mass function",
	 members     =>
	     [
	      "imfSelectionMethod",
	      "imfSelectionFixed",
	      qr/^imf.*RecycledInstantaneous$/,
	      qr/^imf.*YieldInstantaneous$/,
	     ]
     },
     {
	 name        => "stellarEvolution" ,
	 description => "Stellar evolution",
	 members     =>
	     [
	      "stellarPopulationPropertiesMethod",
	      "stellarPopulationSpectraMethod",
	      "stellarPopulationLuminosityIntegrationToleranceRelative"
	     ]
     },
     {
	 name        => "smbhAccretionDisks",
	 description => "Super-massive black hole accretion disks",
	 members     =>
	     [
	      "accretionDisksMethod",
	      "accretionDiskSwitchedScaleAdafRadiativeEfficiency",
	      "accretionRateThinDiskMaximum",
	      "accretionRateThinDiskMinimum",
	      "adafAdiabaticIndex",
	      "adafEnergyOption",
	      "adafRadiativeEfficiency",
	      "adafRadiativeEfficiencyType",
	      "adafViscosityOption"
	     ]
     },
     {
	 name        => "superMassiveBlackHoles" ,
	 description => "Super massive black holes",
	 members     =>
	     [
	      "blackHoleBinaryMergersMethod",
	      "blackHoleHeatsHotHalo",
	      "blackHoleSeedMass",
	      "blackHoleWindEfficiency",
	      "blackHoleWindEfficiencyScalesWithRadiativeEfficiency",
	      "bondiHoyleAccretionEnhancementHotHalo",
	      "bondiHoyleAccretionEnhancementSpheroid",
	      "bondiHoyleAccretionHotModeOnly",
	      "bondiHoyleAccretionTemperatureSpheroid",
	      "spheroidEnergeticOutflowMassRate"
	     ]
     },
     {
	 name        => "galaxyMerging" ,
	 description => "Galaxy merging",
	 members     =>
	     [
	      "satelliteMergingMassMovementsMethod",
	      "satelliteMergingRemnantSizeMethod",
	      "majorMergerMassRatio",
	      "mergerRemnantSizeOrbitalEnergy",
	      "minorMergerGasMovesTo"
	     ]
     },
     {
	 name        => "solvers" ,
	 description => "Solvers and time-stepping",
	 members     =>
	     [
	      "odeAlgorithm",
	      "odeToleranceAbsolute",
	      "odeToleranceRelative",
	      "diskMassToleranceAbsolute",
	      "spheroidMassToleranceAbsolute",
	      "diskVerySimpleMassScaleAbsolute",
	      "hotHaloVerySimpleDelayedMassScaleRelative",
	      qr/^timestep.*/,
	      qr/.*AnalyticSolver.*/
	     ]
     },
     {
	 name        => "output" ,
	 description => "Output",
	 members     =>
	     [
	      "galacticusOutputFileName",
	      "hdf5CompressionLevel",
	      "mergerTreeOutput",
	      "mergerTreeOutputReferences",
	      "outputRedshifts",
	      "outputSatelliteStatus"
	     ]
     },
     {
	 name        => "analysis" ,
	 description => "Analysis",
	 members     =>
	     [
	      "mergerTreeAnalyses",
	      qr/^analysisMassFunction.*/,
	      qr/^analysisMassFunctions.*/,
	     ]
     },
     {
	 name        => "haloModel" ,
	 description => "Halo model [clustering]",
	 members     =>
	     [
	      "haloModelPowerSpectrumModifierMethod",
	     ]
     },
     {
	 name        => "stellarMassSystematics" ,
	 description => "Stellar mass random and systematic error models [constraints]",
	 members     =>
	     [
	      qr/.*StellarMassFunctionZ[\d\.]+MassSystematic\d+$/,
	      qr/^stellarMassSystematics\d+$/,
	      qr/.*StellarMassFunctionZ[\d\.]+MassRandom\d+$/,
	      qr/.*StellarMassFunctionZ[\d\.]+MassRandomMinimum$/,
	      qr/^sdssStellarMassFunctionHighMass.*/,
	      qr/^gamaStellarMassFunctionZ0\.03/
	     ]
     },
     {
	 name        => "hiMassSystematics" ,
	 description => "HI mass random and systematic error models [constraints]",
	 members     =>
	     [
	      qr/alfalfaHiMassFunctionZ0\.00.*$/,
	     ]
     },
     {
	 name        => "sizeSystematics" ,
	 description => "Galaxy size random and systematic error models [constraints]",
	 members     =>
	     [
	      qr/^.*SizeFunctionZ[\d\.]+MassSystematic\d+$/,
	      qr/^.*SizeFunctionZ[\d\.]+MassRandom\d+$/,
	      qr/^.*SizeFunctionZ[\d\.]+MassRandomMinimum$/,
	      qr/^.*SizeFunctionZ[\d\.]+RadiusSystematic\d+$/,
	     ]
     },
     {
	 name        => "clusteringSystematics" ,
	 description => "Clustering random and systematic error models [constraints]",
	 members     =>
	     [
	      qr/^.*ClusteringZ[\d\.]+MassSystematic\d+$/,
	      qr/^.*ClusteringZ[\d\.]+MassRandom\d+$/,
	      qr/^.*ClusteringZ[\d\.]+MassRandomMinimum$/,
	     ]
     }
    );

# Iterate over parameter sets.
foreach my $parameters ( @parameterSets ) {
    # Create groups in this parameter set.
    my $newGroups;
    # Find parameter nodes.
    my @parameterNodes = map {($_->exists('value') || $_->hasAttribute('value') ) ? $_ : ()} $parameters->findnodes('*');
    # Iterate over parameters.
    for my $parameter ( @parameterNodes ) {
	# Check if this node belongs to a group.
	foreach my $group ( @groups ) {
	    if ( 
		grep {
		      (ref($_) eq "Regexp" && $parameter->nodeName() =~ $_)
		      ||
		      (ref($_) eq ""       && $parameter->nodeName() eq $_)
		} @{$group->{'members'}} 
		) {
		# Remove and add to the relevant group.
		push(
		    @{$newGroups->{$group->{'name'}}},
		    $parameter->parentNode()->removeChild($parameter)
		    );
	    }
	}
    }
    # Insert a final description.
    if ( grep {$_->exists('value') || $_->hasAttribute('value')} $parameters->findnodes('*') ) {
	my $comment = $input->createComment(" Miscellaneous parameters ");
	$parameters->insertBefore($comment,$parameters->firstChild());
    }
    # Insert removed parameters.
    foreach my $group ( reverse(@groups) ) {
	# Skip empty groups.
	next
	    unless ( exists($newGroups->{$group->{'name'}}) );
	# Get a sorted list of group members.
	my @sorted  = sort {&sortParameters($group,$a,$b)} @{$newGroups->{$group->{'name'}}};
	# Iterate over group members.
	foreach my $parameter ( @sorted ) {
	    $parameters->insertBefore($parameter,$parameters->firstChild());
	}
	# Insert a description.
	my $comment = $input->createComment(" ".$group->{'description'}." ");
	$parameters->insertBefore($comment,$parameters->firstChild());
    }
}

# Output the resulting file.
my $pp = XML::LibXML::PrettyPrint->new(indent_string => "  ");
$pp->pretty_print($input);
foreach my $parameters ( @parameterSets ) {
    foreach ( $parameters->childNodes() ) {
	if ( $_->nodeType() == XML_COMMENT_NODE ) {
	    my $text = $input->createTextNode("\n  ");
	    $parameters->insertBefore($text,$_);
	}
    }
}
$input->toFile($outputFileName);

exit;

sub sortParameters {
    my $group    = shift();
    my @elements = @_;
    die("sortParameters(): expect two parameters")
	unless ( scalar(@elements) == 2 );
    my %orderBy = map { ${$group->{'members'}}[$_] => $_ } 0..$#{$group->{'members'}};
    my @rank;
    foreach ( @elements ) {
	if ( exists($orderBy{$_->nodeName()}) ) {
	    # We have a direct name match.
	    push(@rank,{rank => $orderBy{$_->nodeName()}, name => $_->nodeName()});
	} else {
	    # No direct name match, we should have a regex match.
	    foreach my $member ( @{$group->{'members'}} ) {
		if ( ref($member) eq "Regexp" && $_->nodeName() =~ $member ) {
		    push(@rank,{rank => $orderBy{$member}, name => $_->nodeName()});
		    last;
		}
	    }
	}
    }
    die("sortParameters(): matching failed")
	unless ( scalar(@rank) == 2 );
    # Do the sort.
    my $order = $rank[1]->{'rank'} <=> $rank[0]->{'rank'};
    $order = $rank[1]->{'name'} cmp $rank[0]->{'name'}
	if ( $order == 0 );
    return $order;
}
