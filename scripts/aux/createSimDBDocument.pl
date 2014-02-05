#!/usr/bin/env perl
use strict;
use warnings;
use XML::Writer;
use XML::Simple;
use IO::File;
use DateTime;
use Data::Dumper;
use LaTeX::Decode;
use LaTeX::BibTeX;
use File::Which;
$LaTeX::Decode::DefaultScheme = 'full';

# Construct a SimDB XML document for Galacticus.
# Andrew Benson (13-August-2011)

# Known parameter group properties.
my $knownGroups;
$knownGroups->{'cosmology'      }->{'name'       } = "Cosmological parameters";
$knownGroups->{'cosmology'      }->{'description'} = "Cosmological parameters";
$knownGroups->{'treeNodeMethods'}->{'name'       } = "Tree node component methods";
$knownGroups->{'treeNodeMethods'}->{'description'} = "Tree node component methods";
$knownGroups->{'starFormation'  }->{'name'       } = "Parameters controlling star formation";
$knownGroups->{'starFormation'  }->{'description'} = "Parameters controlling star formation";
$knownGroups->{'initialMassFunction'  }->{'name'       } = "Parameters controlling the initial mass function";
$knownGroups->{'initialMassFunction'  }->{'description'} = "Parameters controlling the initial mass function";
$knownGroups->{'blackHoles'     }->{'name'       } = "Parameters controlling black holes";
$knownGroups->{'blackHoles'     }->{'description'} = "Parameters controlling black holes";
$knownGroups->{'outputs'        }->{'name'       } = "Parameters controlling outputs";
$knownGroups->{'outputs'        }->{'description'} = "Parameters controlling outputs";
$knownGroups->{'timeStepping'   }->{'name'       } = "Parameters controlling time stepping";
$knownGroups->{'timeStepping'   }->{'description'} = "Parameters controlling time stepping";

# Parse the bibliography.
system("cd doc; Bibliography_Demangle.pl");
my $citations;
my $bibfile = new LaTeX::BibTeX::File "doc/GalacticusAccented.bib";
while (my $entry = new LaTeX::BibTeX::Entry $bibfile)
{
    next unless $entry->parse_ok;

    my @names = $entry->names('author');
    my $authors;
    if ( scalar(@names) == 1 ) {
	my @last = $names[0]->part('last');
	$authors = join(" ",$last[0]);
    } elsif ( scalar(@names) == 2 ) {
	my @last = $names[0]->part('last');
	$authors = join(" ",$last[0]);
	@last = $names[1]->part('last');
	$authors .= " & ".join(" ",$last[0]);
    } elsif ( scalar(@names) == 3 ) {
	my @last = $names[0]->part('last');
	$authors = join(" ",$last[0]);
	@last = $names[1]->part('last');
	$authors .= ", ".join(" ",$last[0]);
	@last = $names[2]->part('last');
	$authors .= " & ".join(" ",$last[0]);
    } else {
	my @last = $names[0]->part('last');
	$authors = join(" ",$last[0])." et al.";
    }
    $citations->{$entry->key()}->{'inline'} = $authors." (".$entry->get('year');
    $citations->{$entry->key()}->{'inline'} .= "; ".$entry->get('url') if ( $entry->exists('url') );
    $citations->{$entry->key()}->{'inline'} .= ")";
    $citations->{$entry->key()}->{'enclosed'} = "(".$authors." ".$entry->get('year');
    $citations->{$entry->key()}->{'enclosed'} .= "; ".$entry->get('url') if ( $entry->exists('url') );
    $citations->{$entry->key()}->{'enclosed'} .= ")";
    $citations->{$entry->key()}->{'alt'} = $authors." ".$entry->get('year');
    $citations->{$entry->key()}->{'alt'} .= "; ".$entry->get('url') if ( $entry->exists('url') );
    $citations->{$entry->key()}->{'author'} = $authors;
    $citations->{$entry->key()}->{'year'} = $entry->get('year');
    $citations->{$entry->key()}->{'year'} .= $entry->get('url') if ( $entry->exists('url') );
}

# Citation markup.
my %cites = (
    "cite" => "inline",
    "citep" => "enclosed",
    "citealt" => "alt",
    "citeauthor" => "author",
    "citeyear" => "year",
    );

# Initalize variable to store parameter groups.
my $groups;

# Initialize variables to store output properties.
my $outputTypes;

# Get code directive locations.
system("make work/build/Code_Directive_Locations.xml &> /dev/null");
my $xml = new XML::Simple;
my $directives = $xml->XMLin("work/build/Code_Directive_Locations.xml");

# Open an output file.
my $output = new IO::File(">galacticus.xml");

# Create an XML writer.
my $writer = new XML::Writer(OUTPUT => $output, ENCODING => "utf-8", NEWLINES => 0, DATA_MODE => 1, DATA_INDENT => 1);

# Write opening.
print $output "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n";

# Write the namespace element.
$writer->startTag("ns2:aSimulator",  "xmlns:ns2" => "http://www.ivoa.net/xml/SimDM/v1.0");

# Write basic properties.
$writer->emptyTag('identity','publisherDID' => "ivo://www.mpa-garching.mpg.de/galacticus");
$writer->dataElement('name','Galacticus');
$writer->dataElement('description','The Galacticus semi-analytic galaxy formation code.');
$writer->dataElement('referenceURL','https://sites.google.com/site/galacticusmodel/');
my $dt = DateTime->now->set_time_zone('local');
(my $tz = $dt->format_cldr("ZZZ")) =~ s/(\d{2})(\d{2})/$1:$2/;
my $now = $dt->ymd."T".$dt->hms.".".$dt->format_cldr("SSS").$tz;
$writer->dataElement('created',$now);
$writer->dataElement('updated',$now);
$writer->dataElement('status','published');
$writer->startTag("contact");
$writer->dataElement('role','owner');
$writer->emptyTag('party','publisherDID' => "ivo://www.mpa-garching.mpg.de/andrewbenson");
$writer->endTag("contact");
$writer->dataElement('code','http://www.ctcp.caltech.edu/galacticus/versions/galacticus_v0.9.0.tar.bz2');
$writer->dataElement('version','0.9.0');

# Write input parameters.    
$writer->comment("InputParameter");

# Open the source diretory, finding F90, cpp, and Inc files.
my %parametersListed;
my @filesToProcess;
opendir(sDir,"source");
while ( my $fileName = readdir(sDir) ) {
    push(@filesToProcess,"source/".$fileName)
	if ( $fileName =~ m/\.F90$/ || $fileName =~ m/\.cpp$/ || $fileName =~ m/\.Inc/ );
}
closedir(sDir);
# Add the component objects file.
push(@filesToProcess,"work/build/objects.nodes.components.Inc");
system("make work/build/objects.nodes.components.Inc");
# Process files.
foreach my $fileName ( @filesToProcess ) {
    # Open the file and scan for parameters.
    my $xmlBuffer;
    open(sFile,$fileName);
    while ( my $line = <sFile> ) {
	if ( $line =~ m/^\s*(\!|\/\/)@(.*)/ ) {
	    my $xmlLine = $2;
	    $xmlBuffer = "" if ( $xmlLine =~ m/^\s*<(inputParameter|outputProperty|outputType|outputPropertyGroup)>\s*$/ );
	    $xmlBuffer .= $xmlLine;
	    if ( $xmlLine =~ m/^\s*<\/inputParameter>\s*$/ ) {
		# Parse the XML.
		my $xml = new XML::Simple;
		my $inputParameter = $xml->XMLin($xmlBuffer);
		$inputParameter->{'description'} =~ s/^\s*//;
		$inputParameter->{'description'} =~ s/\s*$//;
		foreach my $cite ( keys(%cites) ) {
		    while ( $inputParameter->{'description'} =~ m/\\$cite\{([a-zA-Z0-9_\-]+)\}/ ) {
			my $key = $1;
			if ( exists($citations->{$key}) ) {
			    $inputParameter->{'description'} =~ s/\\$cite\{[a-zA-Z0-9_\-]+\}/$citations->{$key}->{$cites{$cite}}/g;
			} else {
			    $inputParameter->{'description'} =~ s/\\$cite\{[a-zA-Z0-9_\-]+\}/???/g;
			}
		    }
		}
		$inputParameter->{'name'} = $inputParameter->{'regEx'}
		    if ( exists($inputParameter->{'regEx'}) );
		unless ( exists($inputParameter->{'name'}) ) {
		    print Dumper($inputParameter);
		    die("createSimDBDocument.pl: inputParameter has no name");
		}
		unless ( exists($parametersListed{$inputParameter->{'name'}}) ) {
		    $parametersListed{$inputParameter->{'name'}} = 1;
		    $writer->startTag("parameter");
		    $writer->emptyTag('identity','publisherDID' => "ivo://www.mpa-garching.mpg.de/galacticus#param/".$inputParameter->{'name'},'xmlId' => 'PAR_'.$inputParameter->{'name'});
		    $writer->dataElement('name'       ,$inputParameter->{'name'       });
		    $writer->dataElement('datatype'   ,$inputParameter->{'type'       });
		    $writer->dataElement('cardinality',$inputParameter->{'cardinality'});
		    $writer->dataElement('description',latex_decode($inputParameter->{'description'}));
		    if ( exists($directives->{$inputParameter->{'name'}}) ) {
			$writer->dataElement('isEnumerated','true');
			my @files;
			if ( UNIVERSAL::isa($directives->{$inputParameter->{'name'}}->{'file'},"ARRAY") ) {
			    @files = @{$directives->{$inputParameter->{'name'}}->{'file'}};
			} else {
			    @files = ( $directives->{$inputParameter->{'name'}}->{'file'} );
			}
			foreach my $file ( @files ) {
			    my $xmlBuffer = "";
			    open(iHndl,$file);
			    my $inModule = 0;
			    my $description = "";
			    while ( my $line = <iHndl> ) {
				$inModule = 1 if ( $line =~ m/^\s*module\s/ );
				$inModule = 0 if ( $line =~ m/^\s*end\+module\s/ || $line =~ m/^\s*contains\s*$/ );
				if ( $inModule == 1 && $line =~ m/^\s*\!\%/ ) {
				    my $text = $line;
				    $text =~ s/^\s*\!\%\s*//;
				    $text =~ s/\s*$//;
				    chomp($text);
				    $description .= $text." ";
				}
				if ( $line =~ m/$inputParameter->{'name'}\s*==\s*['"]([a-zA-Z0-9\_\-\+\s]+)['"]/ ) {
				    my $option = $1;
				    $writer->startTag("validValue");
				    $writer->dataElement('value',$option);
				    foreach my $cite ( keys(%cites) ) {
					while ( $description =~ m/\\$cite\{([a-zA-Z0-9_\-]+)\}/ ) {
					    my $key = $1;
					    if ( exists($citations->{$key}) ) {
						$description =~ s/\\$cite\{[a-zA-Z0-9_\-]+\}/$citations->{$key}->{$cites{$cite}}/g;
					    } else {
						$description =~ s/\\$cite\{[a-zA-Z0-9_\-]+\}/???/g;
					    }
					}
				    }
				    $writer->dataElement('description',$description);
				    $writer->dataElement('title',$option);
				    $writer->endTag("validValue");
				}
			    }
			    close(iHndl);			   
			}
		    }
		    $writer->endTag("parameter");
		    # If the parameter belongs to a group, record this fact.
		    push(@{$groups->{$inputParameter->{'group'}}},$inputParameter->{'name'}) if ( exists($inputParameter->{'group'}) );
		}
	    } elsif ( $xmlLine =~ m/^\s*<\/outputProperty>\s*$/ ) {
		# Parse the XML.
		my $xml = new XML::Simple;
		my $outputProperty = $xml->XMLin($xmlBuffer);
		$outputTypes->{$outputProperty->{'outputType'}}->{'properties'}->{$outputProperty->{'name'}} = $outputProperty;
		$outputTypes->{$outputProperty->{'outputType'}}->{'groups'}->{$outputProperty->{'group'}}->{'members'}->{$outputProperty->{'name'}}= 1 if ( exists($outputProperty->{'group'}) );
	    } elsif ( $xmlLine =~ m/^\s*<\/outputType>\s*$/ ) {
		# Parse the XML.
		my $xml = new XML::Simple;		  
		my $outputType = $xml->XMLin($xmlBuffer);
		$outputTypes->{$outputType->{'name'}}->{'description'} = $outputType->{'description'};
	    } elsif ( $xmlLine =~ m/^\s*<\/outputPropertyGroup>\s*$/ ) {
		# Parse the XML.
		my $xml = new XML::Simple;		  
		my $outputPropertyGroup = $xml->XMLin($xmlBuffer);
		$outputTypes->{$outputPropertyGroup->{'outputType'}}->{'groups'}->{$outputPropertyGroup->{'name'}}->{'description'} = $outputPropertyGroup->{'description'};
	    }
	}
    }
    close(sFile);
}

# Scan for component options.
system("make work/build/objects.tree_node.create.definitions.Inc &> /dev/null");
open(sFile,"work/build/objects.tree_node.create.definitions.Inc");
while ( my $line = <sFile> ) {
    if ( $line =~ m/^\s*type\(varying_string\)\s*::\s*treeNodeMethod(\S*)/ ) {
	my $component = $1;
	$writer->startTag("parameter");
	$writer->emptyTag('identity','publisherDID' => "ivo://www.mpa-garching.mpg.de/galacticus#param/treeNodeMethod".$component,'xmlId' => 'PAR_treeNodeMethod'.$component);
	$writer->dataElement('name'       ,'treeNodeMethod'.$component);
	$writer->dataElement('datatype'   ,"string");
	$writer->dataElement('cardinality',"1");
	$writer->dataElement('description',"Specifies the method to be used for ".$component." components.");
	$writer->dataElement('isEnumerated','true');
	foreach my $file ( @{$directives->{'treeNodeCreateInitialize'}->{'file'}} ) {
	    my $xmlBuffer = "";
	    open(iHndl,$file);
	    my $inModule = 0;
	    my $optionMatches = 0;
	    my $description = "";
	    while ( my $line = <iHndl> ) {
		$optionMatches = 1 if ( $line =~ m/^\s*\!\#\s*<optionName[^>]*>treeNodeMethod$component<\/optionName>/ );
		$inModule = 1 if ( $line =~ m/^\s*module\s/ );
		$inModule = 0 if ( $line =~ m/^\s*end\+module\s/ || $line =~ m/^\s*contains\s*$/ );
		if ( $inModule == 1 && $line =~ m/^\s*\!\%/ ) {
		    my $text = $line;
		    $text =~ s/^\s*\!\%\s*//;
		    chomp($text);
		    $description .= $text;
		}
		if ( $optionMatches == 1 && $line =~ m/componentOption\s*==\s*['"]([a-zA-Z0-9\_\-\+\s]+)['"]/ ) {
		    my $option = $1;
		    $writer->startTag("validValue");
		    $writer->dataElement('value',$option);
		    while ( $description =~ m/\\cite\{([a-zA-Z0-9_\-]+)\}/ ) {
			my $key = $1;
			if ( exists($citations->{$key}) ) {
			    $description =~ s/\\cite\{[a-zA-Z0-9_\-]+\}/$citations->{$key}->{'inline'}/g;
			} else {
			    $description =~ s/\\cite\{[a-zA-Z0-9_\-]+\}/???/g;
			}
		    }
		    while ( $description =~ m/\\citep\{([a-zA-Z0-9_\-]+)\}/ ) {
			my $key = $1;
			if ( exists($citations->{$key}) ) {
			    $description =~ s/\\citep\{[a-zA-Z0-9_\-]+\}/$citations->{$key}->{'enclosed'}/g;
			} else {
			    $description =~ s/\\citep\{[a-zA-Z0-9_\-]+\}/(???)/g;
			}
		    }
		    $writer->dataElement('description',latex_decode($description));
		    $writer->dataElement('title',$option);
		    $writer->endTag("validValue");
		}
	    }
	    close(iHndl);			   
	}
	$writer->endTag("parameter");
	# If the parameter belongs to a group, record this fact.
	push(@{$groups->{'treeNodeMethods'}},"treeNodeMethod".$component);
    }
}
close(sFile);

# Write parameter groups.
$writer->comment("Parameter groups");

foreach my $group ( keys(%{$groups}) ) {
    $writer->startTag("parameterGroup");
    $writer->dataElement('name',$knownGroups->{$group}->{'name'});
    $writer->dataElement('description',$knownGroups->{$group}->{'description'});
    foreach my $parameter ( @{$groups->{$group}} ) {
	$writer->startTag("member");
	$writer->emptyTag('parameter','xmlId' => 'PAR_'.$parameter);
	$writer->endTag("member");
    }
    $writer->endTag("parameterGroup");
}

# Write output types.
$writer->comment("Outputs");

foreach my $outputType ( keys(%{$outputTypes}) ) {
    $writer->startTag("outputType");
    $writer->emptyTag('identity','publisherDID' => 'ivo://www.mpa-garching.mpg.de/galacticus#'.$outputType);
    $writer->dataElement('name'       ,               $outputType                  );
    $writer->dataElement('description',$outputTypes->{$outputType}->{'description'});
    foreach my $outputProperty ( keys(%{$outputTypes->{$outputType}->{'properties'}}) ) {
	$writer->startTag("property");
	$writer->emptyTag('identity','publisherDID' => 'ivo://www.mpa-garching.mpg.de/galacticus#'.$outputType.'/'.$outputProperty, 'xmlId' => $outputType.'_'.$outputProperty);
	$writer->dataElement('name'       ,$outputTypes->{$outputType}->{'properties'}->{$outputProperty}->{'name'       });
	$writer->dataElement('datatype'   ,$outputTypes->{$outputType}->{'properties'}->{$outputProperty}->{'datatype'   });
	$writer->dataElement('cardinality',$outputTypes->{$outputType}->{'properties'}->{$outputProperty}->{'cardinality'});
	$writer->dataElement('description',$outputTypes->{$outputType}->{'properties'}->{$outputProperty}->{'description'});
	$writer->dataElement('label'      ,$outputTypes->{$outputType}->{'properties'}->{$outputProperty}->{'label'      });
	$writer->endTag("property");
    }
    foreach my $outputPropertyGroup ( keys(%{$outputTypes->{$outputType}->{'groups'}}) ) {
	$writer->startTag("propertyGroup");
	$writer->emptyTag('identity');
	$writer->dataElement('name'       ,                                          $outputPropertyGroup                  );
	$writer->dataElement('description',$outputTypes->{$outputType}->{'groups'}->{$outputPropertyGroup}->{'description'});
	foreach my $member ( keys(%{$outputTypes->{$outputType}->{'groups'}->{$outputPropertyGroup}->{'members'}}) ){
	    $writer->startTag("member");
	    $writer->emptyTag('property','xmlId' => $outputType.'_'.$member);
	    $writer->endTag("member");
	}
	$writer->endTag("propertyGroup");
    }
    $writer->endTag("outputType");
}


# Write physics.
my @physicalProcesses = (
    {
	name => "starFormation",
	description => "The process of forming stars from interstellar medium gas"
    },
    {
	name => "dynamicalFriction",
	description => "A dissipative force acting on a gravitating body moving through a background of particles"
    },
    {
	name => "supernovaFeedback",
	description => "The effects of supernova explosions on their host galaxy"
    },
    {
	name => "agnFeedback",
	description => "The effects of active galactic nucleus activity on the host galaxy"
    },
    {
	name => "blackHoleAccretion",
	description => "The accretion of mass onto a black hole from a smooth background"
    },
    {
	name => "blackHoleMerging",
	description => "Coalescence of two black holes"
    },
    {
	name => "blackHoleBinaryHardening",
	description => "Various processes which remove orbital energy from black hole binaries"
    },
    {
	name => "blackHoleMergingRecoil",
	description => "Recoil of merging black hole binaries due to asymmetric emission of gravitational waves"
    },
    {
	name => "gasCooling",
	description => "Any process by which gas is able to radiate energy"
    },
    {
	name => "gravitationalCollapse",
	description => "Runaway collapse of self-gravitating systems"
    },
    {
	name => "hierarchicalStructureFormation",
	description => "The growth of cosmological structures via merging of smaller predecessors"
    },
    {
	name => "galaxyMerging",
	description => "The process of galaxies merging together to form a large galaxy"
    },
    {
	name => "gravitationalInstability",
	description => "The instability (typically of disks) due to self-gravity"
    },
    );

$writer->comment("Physical processes");
foreach my $physicalProcess ( @physicalProcesses ) {
    $writer->startTag("physicalProcess");
    $writer->dataElement('name',$physicalProcess->{'name'});
    $writer->dataElement('description',$physicalProcess->{'description'});
    $writer->dataElement('label',$physicalProcess->{'name'});
    $writer->endTag("physicalProcess");
}

# Close the namespace.
$writer->endTag("ns2:aSimulator");

# Write the document.
$writer->end();

# Close the output file.
$output->close();

# Valiate. DISABLED AS THE WEB SERVICE IS NO LONGER AVAILABLE.
# system("curl --form doc=@`pwd`/galacticus.xml --form action=validate -o validate.html http://galformod.mpa-garching.mpg.de/dev/SimDM-browser/Validate.do");
# my $isOK = 0;
# open(vHndl,"validate.html");
# while ( my $line = <vHndl> ) {
#     $isOK = 1 if ( $line =~ m/Your uploaded document was valid/ );
# }
# close(vHndl);
# print "\n\n\n";
# if ( $isOK == 1 ) {
#     print "SUCCESS: The generated document is valid.\n";
# } else {
#     print "FAIL: The generated document is NOT valid.\n";
# }
# system("firefox file://`pwd`/validate.html") 
#     unless ( $isOK == 1 || ! which('firefox') || ! -t STDIN || ! -t STDOUT );

exit;
