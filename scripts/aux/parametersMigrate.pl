#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Scalar::Util 'reftype';
use XML::LibXML qw(:libxml);
use XML::LibXML::PrettyPrint;
use XML::Simple;
use Data::Dumper;
use Galacticus::Options;
use List::ExtraUtils;
use System::Redirect;
use DateTime;

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
my $elementLastModified = $root->findnodes('//lastModified')->[0];
unless ( defined($elementLastModified) ) {
    my $elementLastModifiedNode = $input->createElement("lastModified");
    my $newBreak                = $input->createTextNode("\n  "   );
    $root->insertBefore($elementLastModifiedNode,$root->firstChild());
    $root->insertBefore($newBreak               ,$root->firstChild());
    $elementLastModified = $root->findnodes('lastModified')->[0];
}
if ( exists($options{'lastModifiedRevision'}) ) {
    $hashLastModified = $options{'lastModifiedRevision'};
} else {
    # Look for a last modification hash.
    if ( defined($elementLastModified) && grep {$_->name() eq 'revision'} $elementLastModified->attributes() ) {
	# A last modified element exists, extract the hash, and update.
	$hashLastModified = $elementLastModified->getAttribute('revision');
    } elsif ( $isInGit ) {
	# File is in git index, use git to determine the last revision at which it was modified.
	## Find the hash at which the file was last modified.
	{
	    open(my $git,"git log -n 1 --pretty=format:\%H -- ".$inputFileName."|");
	    $hashLastModified = <$git>;
	    chomp($hashLastModified);
	}
    } else {
	$hashLastModified           = "6eab8997cd73cb0a474228ade542d133890ad138^";
	my $elementLastModifiedNode = $input->createElement("lastModified");
	my $newBreak                = $input->createTextNode("\n  "   );
	$root->insertBefore($elementLastModifiedNode,$root->firstChild());
	$root->insertBefore($newBreak               ,$root->firstChild());
	$elementLastModified = $root->findnodes('lastModified')->[0];
    }
}
## Update the last modified metadata.
$elementLastModified->setAttribute('revision',$hashHead);
$elementLastModified->setAttribute('time'    ,DateTime->now());

## Find the ancestry.
my @ancestry;
{
    open(my $git,"git rev-list --ancestry-path ".$hashLastModified."..".$hashHead."|");
    while ( my $hash = <$git> ) {
	chomp($hash);
	push(@ancestry,$hash);
    }
}

# Read migration rules.
my $xml        = new XML::Simple();
my $migrations = $xml->XMLin("scripts/aux/migrations.xml");

# Write starting message.
print "Translating file: ".$inputFileName."\n";

# Iterate over parameter sets.
foreach my $parameters ( @parameterSets ) {
    &Migrate($parameters,1,$inputFileName);
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

sub Migrate {
    my $parameters    = shift();
    my $rootLevel     = shift();
    my $inputFileName = shift();

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
	my @matchedMigrations = grep {$_->{'commit'} eq $hashAncestor} &List::ExtraUtils::as_array($migrations->{'migration'});
	next
	    unless ( scalar(@matchedMigrations) == 1 );
	my $migration = $matchedMigrations[0];
	# Report.
	print "Updating to revision ".$hashAncestor."\n";

	# Iterate over translation.
	foreach my $translation ( &List::ExtraUtils::as_array($migration->{'translation'}) ) {
	    
	    # Handle special cases.
	    if ( exists($translation->{'function'}) ) {
		&{\&{$translation->{'function'}}}($input,$parameters);
	    }
	    
	    # Handle removals.
	    if ( exists($translation->{'remove'}) ) {
		foreach my $parameter ( $parameters->findnodes($translation->{'xpath'})->get_nodelist() ) {
		    print "   remove parameter: ".$parameter->nodeName()."\n";
		    $parameter->parentNode->removeChild($parameter);
		}
	    }
	
	    # Translate names.
	    if ( exists($translation->{'name'}) ) {
		foreach my $parameter ( $parameters->findnodes($translation->{'xpath'})->get_nodelist() ) {
		    print "   translate parameter name: ".$parameter->nodeName()." --> ".$translation->{'name'}->{'new'}."\n";
		    my $leafName = $translation->{'name'}->{'new'};
		    my $valueTo;
		    unless ( $translation->{'name'}->{'new'} =~ m/\-\-/ ) {
			if ( $leafName =~ m/(.*)\.(.*)\./ ) {
			    $valueTo    = $2;
			}
		    }
		    $parameter->setNodeName($leafName);
		    $parameter->setAttribute('value',$valueTo)
			if ( $valueTo );
		}
	    }
	
	    # Translate values.
	    if ( exists($translation->{'value'}) ) {
		foreach my $parameter ( $parameters->findnodes($translation->{'xpath'})->get_nodelist() ) {
		    my @allValues;
		    if ( $parameter->exists('value') ) {
			@allValues = $parameter->findnodes('value');
		    } else {
			@allValues = $parameter;
		    }
		    foreach my $value ( @allValues ) {
			# Split values.
			my $valuesText;
			if ( $value->isSameNode($parameter) ) {
			    $valuesText = $value->getAttribute('value');
			} else {
			    $valuesText = $value->firstChild()->textContent();
			}
			$valuesText =~ s/^\s*//;
			$valuesText =~ s/\s*$//;
			my @values = split(/\s+/,$valuesText);
			foreach my $thisValue ( @values ) {
			    if ( $thisValue eq $translation->{'value'}->{'old'} ) {
				print "   translate parameter value: ".$translation->{'xpath'}."\n";
				print "                                 ".$thisValue." --> ".$translation->{'value'}->{'new'}."\n";
				$thisValue = $translation->{'value'}->{'new'};
			    }
			}
			# Update the values in the parameter.
			if ( $value->isSameNode($parameter) ) {
			    $value->setAttribute('value',join(" ",@values));
			} else {
			    $value->firstChild()->setData(join(" ",@values));
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
		    my @defaults = grep {$_->{'parameter'} eq $hostLeafName} @{$migrations->{'default'}};
		    die('parametersMigrate.pl: attempting to insert a "'.$hostLeafName.'" element, but no default value is known')
		    	unless ( scalar(@defaults) == 1 );
		    $defaultValue = $defaults[0]->{'value'};
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

sub radiationFieldIntergalacticBackgroundCMB {
    # Special handling to add CMB radiation into the intergalactic background radiation.
    my $input      = shift();
    my $parameters = shift();
    # Look for "radiationFieldIntergalacticBackground" parameters.
    foreach my $node ( $parameters->findnodes("//radiationFieldIntergalacticBackground[\@value]")->get_nodelist() ) {
	print "   translate special '//radiationFieldIntergalacticBackground[\@value]'\n";
	# Construct new summation and CMB nodes, and a copy of the original node.
	my $summationNode = $input->createElement("radiationFieldIntergalacticBackground");
	my $cmbNode	  = $input->createElement("radiationField"                       );
	my $originalNode  = $input->createElement("radiationField"                       );
	$summationNode->setAttribute('value','summation'                 );
	$cmbNode      ->setAttribute('value','cosmicMicrowaveBackground' );
	$originalNode ->setAttribute('value',$node->getAttribute('value'));
	# Move any children of the original node to the copy.
	my $childNode = $node->firstChild();
	while ( $childNode ) {
	    my $siblingNode = $childNode->nextSibling();
	    $childNode   ->unbindNode(          );
	    $originalNode->addChild  ($childNode);
	    $childNode      = $siblingNode;
	}
	# Assemble our new nodes.
	$summationNode            ->addChild   ($cmbNode            );
	$summationNode            ->addChild   ($originalNode       );
	# Insert the new parameters and remove the original.
	$node         ->parentNode->insertAfter($summationNode,$node);
	$node         ->parentNode->removeChild(               $node);
    }
}

sub blackHoleSeedMass {
    # Special handling to add black hole seed mass models.
    my $input      = shift();
    my $parameters = shift();
    # Look for "componentBlackHole" parameters.
    my $massSeed = 100.0; # This was the default value, to be used if no value was explicitly set.
    my $componentBlackHole;
    foreach my $node ( $parameters->findnodes("//componentBlackHole[\@value='simple' or \@value='standard' or \@value='nonCentral']/massSeed[\@value]")->get_nodelist() ) {
	print "   translate special '//componentBlackHole[\@value]/massSeed[\@value]'\n";
	$componentBlackHole = $node->parentNode;
	# Extract the seed mass.
	$massSeed = $node->getAttribute('value');
	# Delete this node.
	$node->parentNode->removeChild($node);
    }
    # Find nodeOperators.
    my @nodeOperators = $parameters->findnodes("//nodeOperator[\@value='multi']")->get_nodelist();
    die("can not find any `nodeOperator[\@value='multi']` into which to insert a black hole seed operator")
	if ( scalar(@nodeOperators) == 0 );
    die("found multiple `nodeOperator[\@value='multi']` nodes - unknown into which to insert a black hole seed operator")
	if ( scalar(@nodeOperators) >  1 );
    # Create a new node operator and insert into the list.
    my $operatorNode = $input->createElement("nodeOperator"  );
    my $seedNode     = $input->createElement("blackHoleSeeds");
    my $massNode     = $input->createElement("mass"          );
    my $spinNode     = $input->createElement("spin"          );
    $operatorNode->setAttribute('value','blackHolesSeed');
    $seedNode    ->setAttribute('value','fixed'         );
    $massNode    ->setAttribute('value',$massSeed       );
    $spinNode    ->setAttribute('value',0.0             );
    # Assemble our new nodes.
    $seedNode    ->addChild($massNode);
    $seedNode    ->addChild($spinNode);
    # Insert the new parameters.
    $nodeOperators[0]->insertAfter($operatorNode,$nodeOperators[0]->lastChild);
    if ( defined($componentBlackHole) ) {
	$componentBlackHole->parentNode->insertAfter($seedNode,$componentBlackHole);
    } else {
	$parameters->addChild($seedNode);
    }
}

sub blackHolePhysics {
    # Special handling to move black hole physics from components to operators.
    my $input      = shift();
    my $parameters = shift();
    # Set defaults for parameters.
    my $defaults =
    {
	"simple" =>
	{
	    heatsHotHalo                 => "false"    ,
	    efficiencyHeating            => "1.0e-3"   ,
	    efficiencyWind               => "2.2157e-3",
	    growthRatioToStellarSpheroid => "1.0e-3"
	},
	"standard" =>
	{
	    bondiHoyleAccretionEnhancementSpheroid      => "5.0e0" ,
	    bondiHoyleAccretionEnhancementHotHalo       => "6.0e0" ,
	    bondiHoyleAccretionHotModeOnly              => "true"  ,
	    bondiHoyleAccretionTemperatureSpheroid      => "1.0e2" ,
	    efficiencyWind                              => "2.4e-3",
	    efficiencyWindScalesWithEfficiencyRadiative => "false" ,
	    efficiencyRadioMode                         => "1.0"   ,
	    heatsHotHalo                                => "true"
	}
    };
    # Look for "componentBlackHole" parameters.
    my $componentProperties;
    my $componentNode;
    foreach my $node ( $parameters->findnodes("//componentBlackHole[\@value='simple' or \@value='standard' or \@value='nonCentral']")->get_nodelist() ) {
	print "   translate special '//componentBlackHole[\@value]'\n";
	$componentNode = $node;
	# Extract the type of component.
	my $componentType = $node->getAttribute('value');
	if ( $componentType eq "simple" ) {
	    $componentProperties->{'type'} = "simple";
	} elsif ( $componentType eq "standard" || $componentType eq "nonCentral" ) {
	    $componentProperties->{'type'} = "standard";
	}
	# Extract all sub-parameters.
	foreach my $nodeChild ( $node->findnodes("*[\@value]")->get_nodelist() ) {
	    $componentProperties->{$nodeChild->nodeName()} = $nodeChild->getAttribute('value');
	    $node->removeChild($nodeChild);
	}
	# Insert default parameters.
	if ( exists($defaults->{$componentProperties->{'type'}}) ) {
	    foreach my $parameterName ( keys(%{$defaults->{$componentProperties->{'type'}}}) ) {
		$componentProperties->{$parameterName} = $defaults->{$componentProperties->{'type'}}->{$parameterName}
		    unless ( exists($componentProperties->{$parameterName}) );
	    }
	}
    }
    # Find nodeOperators.
    my @nodeOperators = $parameters->findnodes("//nodeOperator[\@value='multi']")->get_nodelist();
    die("can not find any `nodeOperator[\@value='multi']` into which to insert a black hole seed operator")
     	if ( scalar(@nodeOperators) == 0 );
    die("found multiple `nodeOperator[\@value='multi']` nodes - unknown into which to insert a black hole seed operator")
     	if ( scalar(@nodeOperators) >  1 );
    # Build node operators.
    my $operatorAccretionNode = $input->createElement("nodeOperator");
    my $operatorWindsNode     = $input->createElement("nodeOperator");
    my $operatorCGMHeatNode   = $input->createElement("nodeOperator");
    $operatorAccretionNode->setAttribute('value','blackHolesAccretion' );
    $operatorWindsNode    ->setAttribute('value','blackHolesWinds'     );
    $operatorCGMHeatNode  ->setAttribute('value','blackHolesCGMHeating');
    $nodeOperators[0]->insertAfter($operatorAccretionNode,$nodeOperators[0]->lastChild);
    $nodeOperators[0]->insertAfter($operatorWindsNode    ,$nodeOperators[0]->lastChild);
    $nodeOperators[0]->insertAfter($operatorCGMHeatNode  ,$nodeOperators[0]->lastChild)
	if ( $componentProperties->{'heatsHotHalo'} eq "true" );
    # Handle the "simple" black hole component.
    if ( $componentProperties->{'type'} eq "simple" ) {
	{
	    # Accretion.
	    my $accretionNode  = $input->createElement("blackHoleAccretionRate"      );
	    my $growthRateNode = $input->createElement("growthRatioToStellarSpheroid");
	    $accretionNode             ->setAttribute('value'        ,'spheroidTracking'                                    );
	    $growthRateNode            ->setAttribute('value'        ,$componentProperties->{'growthRatioToStellarSpheroid'});
	    $accretionNode             ->addChild    ($growthRateNode                                                       );
	    $componentNode ->parentNode->insertAfter ($accretionNode ,$componentNode                                        );
	}
	{
	    # Winds.
	    my $windNode       = $input->createElement("blackHoleWind" );
	    my $efficiencyNode = $input->createElement("efficiencyWind");
	    $windNode                  ->setAttribute('value'        ,'simple'                                );
	    $efficiencyNode            ->setAttribute('value'        ,$componentProperties->{'efficiencyWind'});
	    $windNode                  ->addChild    ($efficiencyNode                                         );
	    $componentNode ->parentNode->insertAfter ($windNode      ,$componentNode                          );
	}
	if ( $componentProperties->{'heatsHotHalo'} eq "true" ) {
	    # CGM heating.
	    my $heatingNode    = $input->createElement("blackHoleCGMHeating");
	    my $efficiencyNode = $input->createElement("efficiencyHeating"  );
	    $heatingNode               ->setAttribute('value'        ,'quasistatic'                               );
	    $efficiencyNode            ->setAttribute('value'        ,$componentProperties->{'efficiencyHeating'});
	    $heatingNode               ->addChild    ($efficiencyNode                                            );
	    $componentNode ->parentNode->insertAfter ($heatingNode   ,$componentNode                             );
	}
    }
    # Handle the "standard" black hole component.
    if ( $componentProperties->{'type'} eq "standard" ) {
	{
	    # Accretion.
	    my $accretionNode   = $input->createElement("blackHoleAccretionRate"                );
	    my $spheroidNode    = $input->createElement("bondiHoyleAccretionEnhancementSpheroid");
	    my $hotHaloNode     = $input->createElement("bondiHoyleAccretionEnhancementHotHalo" );
	    my $hotModeNode     = $input->createElement("bondiHoyleAccretionHotModeOnly"        );
	    my $temperatureNode = $input->createElement("bondiHoyleAccretionTemperatureSpheroid");
	    $accretionNode             ->setAttribute('value'        ,'standard'                                                     );
	    $spheroidNode              ->setAttribute('value'        ,$componentProperties->{'bondiHoyleAccretionEnhancementSpheroid'});
	    $hotHaloNode               ->setAttribute('value'        ,$componentProperties->{'bondiHoyleAccretionEnhancementHotHalo' });
	    $hotModeNode               ->setAttribute('value'        ,$componentProperties->{'bondiHoyleAccretionHotModeOnly'        });
	    $temperatureNode           ->setAttribute('value'        ,$componentProperties->{'bondiHoyleAccretionTemperatureSpheroid'});
	    $accretionNode             ->addChild    ($spheroidNode                                                                   );
	    $accretionNode             ->addChild    ($hotHaloNode                                                                    );
	    $accretionNode             ->addChild    ($hotModeNode                                                                    );
	    $accretionNode             ->addChild    ($temperatureNode                                                                );
	    $componentNode ->parentNode->insertAfter ($accretionNode ,$componentNode                                                  );
	}
	{
	    # Winds.
	    my $windNode       = $input->createElement("blackHoleWind"                              );
	    my $efficiencyNode = $input->createElement("efficiencyWind"                             );
	    my $scaleNode      = $input->createElement("efficiencyWindScalesWithEfficiencyRadiative");
	    $windNode                  ->setAttribute('value'        ,'ciotti2009'                                                         );
	    $efficiencyNode            ->setAttribute('value'        ,$componentProperties->{'efficiencyWind'                             });
	    $scaleNode                 ->setAttribute('value'        ,$componentProperties->{'efficiencyWindScalesWithEfficiencyRadiative'});
	    $windNode                  ->addChild    ($efficiencyNode                                                                      );
	    $windNode                  ->addChild    ($scaleNode                                                                           );
	    $componentNode ->parentNode->insertAfter ($windNode      ,$componentNode                                                       );
	}
	if ( $componentProperties->{'heatsHotHalo'} eq "true" ) {
	    # CGM heating.
	    my $heatingNode    = $input->createElement("blackHoleCGMHeating");
	    my $efficiencyNode = $input->createElement("efficiencyRadioMode");
	    $heatingNode               ->setAttribute('value'        ,'jetPower'                                   );
	    $efficiencyNode            ->setAttribute('value'        ,$componentProperties->{'efficiencyRadioMode'});
	    $heatingNode               ->addChild    ($efficiencyNode                                              );
	    $componentNode ->parentNode->insertAfter ($heatingNode   ,$componentNode                               );
	}
    }
}

sub modelParameterXPath {
    # Special handling to switch `modelParameter` names to use XPath syntax.
    my $input      = shift();
    my $parameters = shift();
    # Look for "modelParameter" parameters.
    foreach my $modelParameter ( $parameters->findnodes("//modelParameter[\@value='active' or \@value='inactive']")->get_nodelist() ) {
	print "   translate special '//modelParameter[\@value]'\n";
	# Process `name` nodes.
	foreach my $nameNode ( $modelParameter->findnodes("name") ) {
	    my $value = $nameNode->getAttribute('value');	    
	    $value =~ s/::/\//g;                               # Translate the old `::` separator to XPath standard `/`.
	    $value =~ s/([\[\{])(\d+)([\]\}])/$1.($2+1).$3/ge; # Increment indices to XPath standard 1-indexing.
	    $nameNode->setAttribute('value',$value);	    
	}
    }
}
