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
use System::Redirect;

# Update a Galacticus parameter file from its last modified revision to the current revision.
# Andrew Benson (13-September-2014)

# Get arguments.
die("Usage: parametersMigrate.pl <inputFile> <outputFile> [options]")
    unless ( scalar(@ARGV) >= 2 );
my $inputFileName  = $ARGV[0];
my $outputFileName = $ARGV[1];
my %options =
    (
     validate            => "yes",
     prettyify           => "no" ,
     outputFormatVersion => 2
    );
# Parse options.
my %optionsDefined = &Galacticus::Options::Parse_Options(\@ARGV,\%options);

# Parse the input file.
my $parser     = XML::LibXML->new();
my $input      = $parser->parse_file($inputFileName);
my $root;
my @parameterSets;
if ( $input->findnodes('parameters') ) {
    @parameterSets = $input->findnodes('parameters');
    $root          = $parameterSets[0];
} elsif ( $input->findnodes('parameterGrid') ) {
    my $parameterGrid = $input->findnodes('parameterGrid')->[0];
    $root             = $parameterGrid;
    if ( $parameterGrid->findnodes('parameters') ) {
	@parameterSets = $parameterGrid->findnodes('parameters');
    } else {
	die('can not find parameters')
    }
} else {
    die('can not find parameters')
}

# Find the ancestry of the given file. We find the hash at which is was last modified, the current HEAD hash, and then find the
# ancestry between them. We will then check each hash in the ancestry for parameter migrations.
## Check if the file is under version control.
&System::Redirect::tofile("git ls-files --error-unmatch ".$inputFileName,"/dev/null");
my $isInGit = $? == 0;
## Determine ancestry.
## Find the current hash.
my $hashHead;
{
    open(my $git,"git rev-parse HEAD|");
    $hashHead = <$git>;
    chomp($hashHead);
}
## Find last modified hash.
my $hashLastModified;
if ( $isInGit ) {
    # File is in git index, use git to determine the last revision at which it was modified.
    ## Find the hash at which the file was last modified.
    {
	open(my $git,"git log -n 1 --pretty=format:\%H -- ".$inputFileName."|");
	$hashLastModified = <$git>;
	chomp($hashLastModified);
    }
} else {
    # Look for a last modification hash.
    my $elementLastModified = $input->findnodes('lastModified')->[0];
    if ( defined($elementLastModified) ) {
	# A last modified element exists, extract the hash, and update.
	$hashLastModified = $elementLastModified->textContent();
	$elementLastModified->firstChild()->setData($hashHead);
    } else {
	# No last modified element exists. Set the last modified hash that immediately prior to the earliest one that this script
	# is aware of, and add a last modified element.
	$hashLastModified = "6eab8997cd73cb0a474228ade542d133890ad138^";
	my $elementLastModifiedNode = $input->createElement("lastModified");
	my $newLastModified         = $input->createTextNode($hashHead);
	my $newBreak                = $input->createTextNode("\n  "   );
	$elementLastModifiedNode->addChild($newLastModified);
	$root->insertBefore($elementLastModifiedNode,$root->firstChild());
	$root->insertBefore($newBreak               ,$root->firstChild());
    }
}
## Find the ancestry.
my @ancestry;
{
    open(my $git,"git rev-list --ancestry-path ".$hashLastModified."..".$hashHead."|");
    while ( my $hash = <$git> ) {
	chomp($hash);
	push(@ancestry,$hash);
    }
}

# Define translations.
## These enumerate the parameter changes that occurred at each revision hash.
my %translations =
    (
     "0ab22d1b3959fc9242ee6314f3e199cd99600406" => {
	 names => {
	     "diskMassDistribution"                                 => "componentDisk--massDistributionDisk"                            ,
	     "spheroidMassDistribution"                             => "componentSpheroid--massDistributionSpheroid"                    ,
	 }
     },
     "6eab8997cd73cb0a474228ade542d133890ad138" => {
	 names => {
	     "bondiHoyleAccretionEnhancementSpheroid"               => "componentBlackHole--bondiHoyleAccretionEnhancementSpheroid"     ,
	     "bondiHoyleAccretionEnhancementHotHalo"                => "componentBlackHole--bondiHoyleAccretionEnhancementHotHalo"      ,
	     "bondiHoyleAccretionHotModeOnly"                       => "componentBlackHole--bondiHoyleAccretionHotModeOnly"             ,
	     "bondiHoyleAccretionTemperatureSpheroid"               => "componentBlackHole--bondiHoyleAccretionTemperatureSpheroid"     ,
	     "blackHoleSeedMass"                                    => "componentBlackHole--massSeed"                                   ,
             "tripleBlackHoleInteraction"                           => "componentBlackHole--tripleInteraction"                          ,
	     "blackHoleToSpheroidStellarGrowthRatio"                => "componentBlackHole--growthRatioToStellarSpheroid"               ,
	     "blackHoleHeatsHotHalo"                                => "componentBlackHole--heatsHotHalo"                               ,
	     "blackHoleHeatingEfficiency"                           => "componentBlackHole--efficiencyHeating"                          ,
	     "blackHoleRadioModeFeedbackEfficiency"                 => "componentBlackHole--efficiencyRadioMode"                        ,
	     "blackHoleWindEfficiency"                              => "componentBlackHole--efficiencyWind"                             ,
	     "blackHoleWindEfficiencyScalesWithRadiativeEfficiency" => "componentBlackHole--efficiencyWindScalesWithEfficiencyRadiative",
	     "blackHoleOutputAccretion"                             => "componentBlackHole--outputAccretion"                            ,
	     "blackHoleOutputData"                                  => "componentBlackHole--outputData"                                 ,
	     "blackHoleOutputMergers"                               => "componentBlackHole--outputMergers"                              ,
	     "spheroidEnergeticOutflowMassRate"                     => "componentSpheroid--efficiencyEnergeticOutflow"                  ,
	     "spheroidMetallicityTolerance"                         => "componentSpheroid--toleranceRelativeMetallicity"                ,
	     "spheroidMassToleranceAbsolute"                        => "componentSpheroid--toleranceAbsoluteMass"                       ,
	     "spheroidLuminositiesStellarInactive"                  => "componentSpheroid--inactiveLuminositiesStellar"                 ,
	     "spheroidAngularMomentumAtScaleRadius"                 => "componentSpheroid--ratioAngularMomentumScaleRadius"             ,
	     "spheroidVerySimpleMassScaleAbsolute"                  => "componentSpheroid--scaleAbsoluteMass"                           ,
	     "diskMassToleranceAbsolute"                            => "componentDisk--toleranceAbsoluteMass"                           ,
	     "diskMetallicityTolerance"                             => "componentDisk--toleranceRelativeMetallicity"                    ,
	     "diskStructureSolverRadius"                            => "componentDisk--radiusStructureSolver"                           ,
	     "diskRadiusSolverCole2000Method"                       => "componentDisk--structureSolverUseCole2000Method"                ,
	     "diskLuminositiesStellarInactive"                      => "componentDisk--inactiveLuminositiesStellar"                     ,
	     "diskVerySimpleMassScaleAbsolute"                      => "componentDisk--scaleAbsoluteMass"                               ,
	     "diskVerySimpleTrackAbundances"                        => "componentDisk--trackAbundances"                                 ,
	     "diskVerySimpleTrackLuminosities"                      => "componentDisk--trackLuminosities"                               ,
	     "diskVerySimpleUseAnalyticSolver"                      => "componentDisk--useAnalyticSolver"                               ,
	     "diskVerySimpleAnalyticSolverPruneMassGas"             => "componentDisk--pruneMassGas"                                    ,
	     "diskVerySimpleAnalyticSolverPruneMassStars"           => "componentDisk--pruneMassStars"                                  ,
	     "hotHaloOutflowToColdMode"                             => "componentHotHalo--outflowToColdMode"                            ,
	     "hotHaloVerySimpleDelayedMassScaleRelative"            => "componentHotHalo--scaleRelativeMass"                            ,
	     "starveSatellites"                                     => "componentHotHalo--starveSatellites"                             ,
	     "starveSatellitesOutflowed"                            => "componentHotHalo--starveSatellitesOutflowed"                    ,
	     "hotHaloTrackStrippedGas"                              => "componentHotHalo--trackStrippedGas"                             ,
	     "hotHaloOutflowReturnOnFormation"                      => "componentHotHalo--outflowReturnOnFormation"                     ,
	     "hotHaloAngularMomentumAlwaysGrows"                    => "componentHotHalo--angularMomentumAlwaysGrows"                   ,
	     "hotHaloCoolingFromNode"                               => "componentHotHalo--coolingFromNode"                              ,
	     "hotHaloOutflowStrippingEfficiency"                    => "componentHotHalo--efficiencyStrippingOutflow"                   ,
	     "hotHaloExpulsionRateMaximum"                          => "componentHotHalo--rateMaximumExpulsion"                         ,
	     "hotHaloAngularMomentumLossFraction"                   => "componentHotHalo--fractionLossAngularMomentum"                  ,
	     "hotHaloNodeMergerLimitBaryonFraction"                 => "componentHotHalo--fractionBaryonLimitInNodeMerger"              ,
	     "satelliteBoundMassIsInactive"                         => "componentSatellite--inactiveBoundMass"                          ,
	     "satelliteBoundMassInitializeType"                     => "componentSatellite--initializationTypeMassBound"                ,
	     "satelliteMaximumRadiusOverVirialRadius"               => "componentSatellite--radiusMaximumOverRadiusVirial"              ,
	     "satelliteDensityContrast"                             => "componentSatellite--densityContrastMassBound"                   ,
	 },
	 values => {
	     "spheroidVerySimpleTrackAbundances"                    => {}                                                               ,
	     "spheroidVerySimpleTrackLuminosities"                  => {}                                                               ,
	 }
     }
    );
     
# Define known defaults.
my %knownDefaults =
    (
     "componentBlackHole"             => "standard"        ,
     "componentDisk"                  => "standard"        ,
     "componentHotHalo"               => "standard"        ,
     "componentSpheroid"              => "standard"        ,
     "cosmologyParameters"            => "simple"          ,
     "mergerTreeBuilder"              => "cole2000"        ,
     "cosmologicalMassVariance"       => "filteredPower"   ,
     "cosmologyParameters"            => "simple"          ,
     "darkMatterProfileConcentration" => "gao2008"         ,
     "powerSpectrumPrimordial"        => "powerLaw"        ,
     "haloMassFunction"               => "tinker2008"      ,
     "satelliteMergingTimescales"     => "jiang2008"       ,
     "starFormationHistory"           => "metallicitySplit"
  );

# Write starting message.
print "Translating file: ".$inputFileName."\n";

# Iterate over parameter sets.
foreach my $parameters ( @parameterSets ) {
    &Translate($parameters,1,$inputFileName);
}

# Output the resulting file.
my $pp = XML::LibXML::PrettyPrint->new(indent_string => "  ");
$pp->pretty_print($input)
    if ( $options{'prettyify'} eq "yes" );
my $serialized = $input->toString();
$serialized =~ s/><!\-\-/>\n\n  <!\-\-/gm;
$serialized =~ s/></>\n\n  </gm;
open(my $outputFile,">",$outputFileName);
print $outputFile $serialized;
close($outputFile);

exit;

sub Translate {
    my $parameters    = shift();
    my $rootLevel     = shift();
    my $inputFileName = shift();

    # Check for a version element.
    my $version    = $parameters->findnodes('version')->[0];
    if ( defined($version) ) {
	$version->firstChild()->setData("0.9.4");
    } elsif ( $rootLevel ) {
	my $versionNode = $input->createElement ("version");
	my $newVersion  = $input->createTextNode("0.9.4"  );
	my $newBreak    = $input->createTextNode("\n  "   );
	$versionNode->addChild($newVersion);
	$parameters->insertBefore($versionNode,$parameters->firstChild());
	$parameters->insertBefore($newBreak   ,$parameters->firstChild());
    }
    
    # Determine output format version.
    $options{'outputFormatVersion'} = 2
	unless ( exists($options{'outputFormatVersion'}) );
    
    # Check for a format version element.
    my $formatVersion = $parameters->findnodes('formatVersion')->[0];
    if ( defined($formatVersion) ) {
	$formatVersion->firstChild()->setData($options{'outputFormatVersion'});
    } elsif ( $rootLevel ) {
	my $formatVersionNode = $input->createElement ("formatVersion"                );
	my $newFormatVersion  = $input->createTextNode($options{'outputFormatVersion'});
	my $newBreak          = $input->createTextNode("\n  "                         );
	$formatVersionNode->addChild($newFormatVersion);
	$parameters->insertBefore($formatVersionNode,$parameters->firstChild());
	$parameters->insertBefore($newBreak         ,$parameters->firstChild());
    }

    # Validate the parameter file.
    if ( $options{'validate'} eq "yes" ) {
	system($ENV{'GALACTICUS_EXEC_PATH'}."/scripts/aux/validateParameters.pl ".$inputFileName);
	die('input file "'.$inputFileName.'"is not a valid Galacticus parameter file')
	    unless ( $? == 0 );
    }
    
    # Iterate over the revision ancestry.
    foreach my $hashAncestor ( reverse(@ancestry) ) {
	next
	    unless ( exists($translations{$hashAncestor}) );
	my $translation = $translations{$hashAncestor};
	# Report.
	print "Updating to revision ".$hashAncestor."\n";
	# Find parameter nodes.
	my @parameterNodes = map {( $_->exists('value') || $_->hasAttribute('value') ) ? $_ : ()} $parameters->findnodes('*');
	# Apply translations.
	for my $parameter ( @parameterNodes ) {
	    # Get name and value text elements.
	    my $name     = $parameter;
	    my $nameText = $name->nodeName();
	    my @allValues;
	    if ( $parameter->exists('value') ) {
		@allValues = $parameter->findnodes('value');
	    } else {
		@allValues = $parameter;
	    }
	    # Ignore parameters with no value - sometimes this is used as a placeholder.
	    next
		if ( ! grep {my $valueText = $_->isSameNode($name) ? $_->getAttribute('value') : $_->firstChild()->textContent(); $valueText ne ""} @allValues );	    
	    # Translate names.
	    if ( exists($translation->{'names'}->{$nameText}) ) {
		print "   translate parameter name: ".$nameText." --> ".$translation->{'names'}->{$nameText}."\n";
		my $leafName = $translation->{'names'}->{$nameText};
		my $valueTo;
		unless ( $translation->{'names'}->{$nameText} =~ m/\-\-/ ) {
		    if ( $leafName =~ m/(.*)\.(.*)\./ ) {
			$valueTo    = $2;
		    }
		}
		$name->setNodeName($leafName);
		$name->setAttribute('value',$valueTo)
		    if ( $valueTo );
	    }
	    # Translate values.
	    foreach my $value ( @allValues ) {
		if ( exists($translation->{'values'}->{$nameText}) ) {
		    # Split values.
		    my $valuesText;
		    if ( $value->isSameNode($name) ) {
			$valuesText = $value->getAttribute('value');
		    } else {
			$valuesText = $value->firstChild()->textContent();
		    }
		    $valuesText =~ s/^\s*//;
		    $valuesText =~ s/\s*$//;
		    my @values = split(/\s+/,$valuesText);
		    foreach my $thisValue ( @values ) {
			if ( exists($translation->{'values'}->{$nameText}->{$thisValue}) ) {
			    print "   translate parameter value: ".$nameText."\n";
			    if ( ref($translation->{'values'}->{$nameText}->{$thisValue}) ) {
				my $newValue = $translation->{'values'}->{$nameText}->{$thisValue};
				if ( exists($newValue->{'value'}) ) {
				    print "                                 ".$thisValue." --> ".$newValue->{'value'}."\n";
				    $thisValue = $newValue->{'value'};
				} else {
				    print "                                 ".$thisValue." --> {remove}\n";
				    $thisValue = undef();
				}
				if ( exists($newValue->{'new'}) ) {
				    foreach my $newParameter ( &List::ExtraUtils::as_array($newValue->{'new'}) ) {
					print "      add parameter: ".$newParameter->{'name'}." = ".$newParameter->{'value'}."\n";
					my $parameterNode = $input->createElement($newParameter->{'name'});
					$parameterNode->setAttribute('value',$newParameter->{'value'});
					$parameters->insertAfter($parameterNode,$parameter);
					$parameters->insertAfter($input->createTextNode("\n  "),$parameter);
				    }
				}
			    } else {
				print "                                 ".$thisValue." --> ".$translation->{'values'}->{$nameText}->{$thisValue}."\n";
				$thisValue = $translation->{'values'}->{$nameText}->{$thisValue};
			    }
			}
			if ( scalar(@values) == 1 && ! defined($values[0]) ) {
			    if ( $value->isSameNode($name) ) {
				$value->parentNode->removeChild($value);
			    } else {
				$value->removeChild($value->firstChild());
			    }
			} else {
			    if ( $value->isSameNode($name) ) {
				$value->setAttribute('value',join(" ",@values));
			    } else {
				$value->firstChild()->setData(join(" ",@values));
			    }
			}
		    }
		}
	    }
	}
    }

    # Put subparameters into their host parameter.
    for my $parameter ( $parameters->getChildrenByTagName("*") ) {
	if ( $parameter->nodeName() =~ m/(.+)\-\-(.+)/ ) {
	    my $hostName = $1;
	    my $subName  = $2;
	    (my $hostNameTrimmed = $hostName) =~ s/\..+\.//;
	    $parameters->removeChild($parameter->nextSibling());
	    my $sibling = $parameter->nextSibling();
	    my $hostFound;
	    $parameter->setNodeName($subName);
	    for my $hostParameter ( $parameters->getChildrenByTagName("*") ) {
		if ( $hostParameter->nodeName() eq $hostNameTrimmed ) {
		    my $hostChildren = $hostParameter->getChildrenByTagName("*");
		    $hostParameter->addChild($input     ->createTextNode("\n"      ))
			if ( $hostChildren->size() == 0 );
		    $hostParameter->addChild($input     ->createTextNode("  "      ));
		    $hostParameter->addChild($parameters->removeChild   ($parameter));
		    $hostParameter->addChild($input     ->createTextNode("\n  "    ));
		    $hostFound = 1;
		}
	    }
	    unless ( $hostFound ) {
		# Create the new node.
		(my $hostLeafName = $hostName) =~ s/\..+\.//;
		my $defaultValue;
		if ( $hostName =~ m/\..+\./ ) {
		    ($defaultValue = $hostName) =~ s/.*\.(.+)\./$1/;
		} else {
		    die('parametersMigrate.pl: attempting to insert a "'.$hostLeafName.'" element, but no default value is known')
			unless ( exists($knownDefaults{$hostLeafName}) );
		    $defaultValue = $knownDefaults{$hostLeafName};
		}
		my $parameterNode = $input->createElement($hostLeafName);
		$parameterNode->setAttribute('value',$defaultValue);
		$parameterNode->addChild($input     ->createTextNode("\n    "  ));
		$parameterNode->addChild($parameters->removeChild   ($parameter));
		$parameterNode->addChild($input     ->createTextNode("\n  "    ));
		if ( defined($sibling) ) {
		    $parameters->insertBefore($parameterNode,$sibling);
		} else {
		    $parameters->insertBefore($parameterNode,undef());
		}
	    }
	}
    }
    # Strip out value extensions.
    for my $parameter ( $parameters->getChildrenByTagName("*") ) {
	if ( $parameter->nodeName() =~ m/([^\.]+)\./ ) {
	    my $nodeName = $1;
	    $parameter->setNodeName($nodeName);
	}
    }
    # Search for duplicated parameters.
    my %duplicates =
	(
	 mergerTreeOperatorMethod => "sequence"
	);
   for my $parameter ( $parameters->getChildrenByTagName("*") ) {
	my $nodeName = $parameter->nodeName();
	my $duplicates = $parameters->getChildrenByTagName($nodeName);
	if ( $duplicates->size() > 1 && exists($duplicates{$nodeName}) ) {
	    my $sibling = $parameter->nextSibling();
	    my $wrapperNode = $input->createElement($nodeName);
	    $wrapperNode->setAttribute('value',$duplicates{$nodeName});
	    foreach my $duplicate ( $duplicates->get_nodelist() ) {
		$wrapperNode->addChild($input     ->createTextNode("\n    "  ));
		$wrapperNode->addChild($parameters->removeChild   ($duplicate));
		$wrapperNode->addChild($input     ->createTextNode("\n  "    ));
	    }
	    $parameters->insertBefore($wrapperNode,$sibling);
	}
    }
}
