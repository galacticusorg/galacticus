#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Scalar::Util 'reftype';
use XML::LibXML qw(:libxml);
use XML::LibXML::PrettyPrint;
use Data::Dumper;
use Galacticus::Options;
use List::ExtraUtils;

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
%options = &Galacticus::Options::Parse_Options(\@ARGV,\%options);

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
	      "cosmologyParameters",
	      "cosmologyFunctions",
	      qr/cosmology\d/
	     ]
     },
     {
	 name        => "powerSpectrum" ,
	 description => "Power spectrum",
	 members     =>
	     [
	      "powerSpectrum",
	      "powerSpectrumPrimordial",
	      "transferFunction",
	      "powerSpectrumPrimordialTransferred",
	      "cosmologicalMassVariance",
	      "powerSpectrumNonlinear"
	     ]
     },
     {
	 name        => "structureFormation" ,
	 description => "Structure formation",
	 members     =>
	     [
	      "linearGrowth",
	      "criticalOverdensity",
	      "virialDensityContrast",
	      qr/virialDensityContrast.*/,
	      "haloMassFunction",
	      qr/haloMassFunctionSystematic\d/,
	      "gravitationalLensing"
	     ]
     },
     {
	 name        => "igmReionization" ,
	 description => "Intergalactic medium and reionization",
	 members     =>
	     [
	      "intergalacticMediumState",
	      qr/intergalaticMediumState.*/
	     ]
     },
     {
	 name        => "components" ,
	 description => "Tree node component selection",
	 members     =>
	     [
	      qr/^component/,
	      "nodeComponentBasicExtendedSphericalCollapseType"
	     ]
     },
     {
	 name        => "mergerTreeConstruction" ,
	 description => "Merger tree construction",
	 members     =>
	     [
	      "mergerTreeConstruct",
	      "mergerTreeBuilder",
	      "mergerTreeMassResolution",
	      qr/^mergerTreeBuild.*/,
	      qr/^mergerTreesBuild.*/,
	      "treeBranching",
	      qr/^modifiedPressSchechter.*/,
	      "haloMassFunctionSampling",
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
	 name        => "haloHierarchy" ,
	 description => "Merger tree halo hierarchy",
	 members     =>
	     [
	      "nodeMergers",
	      "nodePromotionIndexShift"
	     ]
     },
     {
	 name        => "darkMatterHalo" ,
	 description => "Dark matter halo structure",
	 members     =>
	     [
	      "darkMatterProfile",
	      "darkMatterProfileConcentration",
	      "darkMatterProfileMinimumConcentration",
	      "darkMatterProfileScaleCorrectForConcentrationDefinition"
	     ]
     },
     {
	 name        => "haloSpin" ,
	 description => "Dark matter halo spin",
	 members     =>
	     [
	      "haloSpinDistribution",
	      "randomSpinResetMassFactor"
	     ]
     },
     {
	 name        => "satelliteOrbits" ,
	 description => "Satellite orbits",
	 members     =>
	     [
	      "virialOrbit",
	      qr/^virialOrbits.*/,
	      "satelliteMergingTimescales",
	      "mergingTimescaleMultiplier",
	      "satellitesTidalField",
	      "satelliteOrbitStoreOrbitalParameters"
	     ]
     },
     {
	 name        => "hotAtmosphere" ,
	 description => "Hot atmosphere",
	 members     =>
	     [
	      "hotHaloOutflowReincorporation",
	      qr/^hotHaloOutflowReincorporation.*/,
	      "hotHaloTemperatureProfile",
	      "hotHaloTrackStrippedGas",
	      "hotHaloOutflowReturnRate",
	      "hotHaloOutflowToColdMode",
	      "hotHaloAngularMomentumAlwaysGrows",
	      "hotHaloAngularMomentumLossFraction",
	      "hotHaloColdModeCoredIsothermalCoreRadii",
	      "hotHaloCoreRadiusOverVirialRadius",
	      "hotHaloExcessHeatDrivesOutflow",
	      "hotHaloMassDistribution",
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
	 name        => "accretionOntoHalos" ,
	 description => "Accretion onto halos",
	 members     =>
	     [
	      "accretionHaloTotal"
	     ]
     },
    {
	 name        => "accretionFromIGM" ,
	 description => "Accretion from the IGM",
	 members     =>
	     [
	      "accretionHalo",
	     ]
     },
     {
	 name        => "hotHaloRamPressure" ,
	 description => "Ram pressure stripping of hot atmosphere",
	 members     =>
	     [
	      "starveSatellites",
	      "hotHaloRamPressureForce",
	      "hotHaloRamPressureStripping",
	      "hotHaloRamPressureStrippingTimescale",
	      "hotHaloOutflowStrippingEfficiency",
	      "ramPressureStrippingFormFactor"
	     ]
     },
     {
	 name        => "coolingAndInfall" ,
	 description => "Cooling and infall",
	 members     =>
	     [
	      "coolingRate",
	      qr/^coolingRate.*/,
	      "coolingFunction",
	      "coolingMeanAngularMomentumFrom",
	      "coolingRadius",
	      "infallRadius",
	      "coolingRotationVelocityFrom",
	      "coolingSpecificAngularMomentum",
	      "coolingTime",
	      qr/^coolingTime.*/,
	      "coldModeInfallRate",
	      qr/^coldModeInfallRate.*/,
	      "zeroCoolingRateAboveVelocity",
	      "chemicalState"
	     ]
     },
     {
	 name        => "galacticStructure" ,
	 description => "Galactic structure",
	 members     =>
	     [
	      "diskMassDistribution",
	      "spheroidMassDistribution",
	      "galacticStructureRadiusSolver",
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
	      "barInstability",
	      "stabilityThresholdGaseous",
	      "stabilityThresholdStellar"
	     ]
     },
     {
	 name        => "diskRamPressure" ,
	 description => "Ram pressure stripping of galactic disk",
	 members     =>
	     [
	      "ramPressureStrippingMassLossRateDisks"
	     ]
     },
     {
	 name        => "spheroidRamPressure" ,
	 description => "Ram pressure stripping of galactic spheroid",
	 members     =>
	     [
	      "ramPressureStrippingMassLossRateSpheroids"
	     ]
     },
     {
	 name        => "diskTidal" ,
	 description => "Tidal stripping of galactic disk",
	 members     =>
	     [
	      "tidalStrippingMassLossRateDisks"
	     ]
     },
     {
	 name        => "spheroidTidal" ,
	 description => "Tidal stripping of galactic spheroid",
	 members     =>
	     [
	      "tidalStrippingMassLossRateSpheroids"
	     ]
     },
     {
	 name        => "diskStarFormation" ,
	 description => "Star formation in disks",
	 members     =>
	     [
	      "starFormationTimescaleDisks",
	      qr/^starFormationTimescaleDisks.*/,
	      qr/^diskStarFormation.*/,
	      qr/^diskVerySimpleSurfaceDensity.*/,
	      "starFormationDiskMinimumTimescale",
	      "starFormationFrequencyKMT09",
	      "molecularComplexClumpingFactorKMT09",
	      "molecularFractionFastKMT09",
	      "starFormationRateSurfaceDensityDisks",
	      "starFormationTimescaleIntegratedSurfaceDensityTolerance"
	     ]
     },
     {
	 name        => "spheroidStarFormation" ,
	 description => "Star formation in spheroids",
	 members     =>
	     [
	      "starFormationTimescaleSpheroids",
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
	      "starFormationFeedbackDisks",
	      qr/^diskOutflow.*/
	     ]
     },
     {
	 name        => "spheroidStellarFeedback" ,
	 description => "Stellar feedback in spheroids",
	 members     =>
	     [
	      "starFormationFeedbackSpheroids",
	      qr/^spheroidOutflow.*/
	     ]
     },
     {
	 name        => "initialMassFunction" ,
	 description => "Stellar initial mass function",
	 members     =>
	     [
	      "imfSelection",
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
	      "stellarPopulationProperties",
	      "stellarPopulationSpectra",
	      "stellarPopulationLuminosityIntegrationToleranceRelative"
	     ]
     },
     {
	 name        => "smbhAccretionDisks",
	 description => "Super-massive black hole accretion disks",
	 members     =>
	     [
	      "accretionDisks",
	      "accretionDiskSwitchedScaleAdafRadiativeEfficiency",
	      "accretionRateThinDiskMaximum",
	      "accretionRateThinDiskMinimum",
	      "accretionRateTransitionWidth",
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
	      "blackHoleBinaryMergers",
	      "blackHoleHeatsHotHalo",
	      "blackHoleSeedMass",
	      "blackHoleWindEfficiency",
	      "blackHoleWindEfficiencyScalesWithRadiativeEfficiency",
	      "blackHoleRadioModeFeedbackEfficiency",
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
	      "satelliteMergingMassMovements",
	      "satelliteMergingRemnantSize",
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
	      "outputSatelliteStatus",
	      "outputColdModeInfallRate"
	     ]
     },
     {
	 name        => "luminosities" ,
	 description => "Stellar luminosities",
	 members     =>
	     [
	      "luminosityFilter",
	      "luminosityRedshift",
	      "luminosityType",
	     ]
     },
     {
	 name        => "numerics" ,
	 description => "Numerics",
	 members     =>
	     [
	      "randomSeed",
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
	      "haloModelPowerSpectrumModifier",
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
	      qr/.*StellarMassFunctionZ[\d\.]+MassRandomMaximum$/,
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
	      qr/^.*SizeFunctionZ[\d\.]+MassRandomMaximum$/,
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
