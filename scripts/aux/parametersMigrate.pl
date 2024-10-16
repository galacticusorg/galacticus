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
if ( exists($options{'lastModifiedRevision'}) ) {
    $hashLastModified = $options{'lastModifiedRevision'};
} elsif ( $isInGit ) {
    # File is in git index, use git to determine the last revision at which it was modified.
    ## Find the hash at which the file was last modified.
    {
	open(my $git,"git log -n 1 --pretty=format:\%H -- ".$inputFileName."|");
	$hashLastModified = <$git>;
	chomp($hashLastModified);
    }
} else {
    # Look for a last modification hash.
    my $elementLastModified = $root->findnodes('lastModified')->[0];
    if ( defined($elementLastModified) ) {
	# A last modified element exists, extract the hash, and update.
	$hashLastModified = $elementLastModified->getAttribute('revision');
   } else {
	# No last modified element exists. Set the last modified hash that immediately prior to the earliest one that this script
	# is aware of, and add a last modified element.
	$hashLastModified           = "6eab8997cd73cb0a474228ade542d133890ad138^";
	my $elementLastModifiedNode = $input->createElement("lastModified");
	my $newBreak                = $input->createTextNode("\n  "   );
	$root->insertBefore($elementLastModifiedNode,$root->firstChild());
	$root->insertBefore($newBreak               ,$root->firstChild());
    }
}
## Update the last modified metadata.
my $elementLastModified = $root->findnodes('lastModified')->[0];
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
