# Contains a Perl module which analyzes pre-processed source code to extract docuemntation on classes.

package Galacticus::Build::SourceTree::Process::ClassDocumentation;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Data::Dumper;
use List::ExtraUtils;
use List::Uniq qw(uniq);
use Fortran::Utils;
use XML::Simple;
use Storable qw(dclone);

# Insert hooks for our functions.
$Galacticus::Build::SourceTree::Hooks::processHooks       {'classDocumentation'} = \&Process_ClassDocumentation;
$Galacticus::Build::SourceTree::Hooks::processDependencies{'classDocumentation'} = [ "functionClass", "eventHooks", "stateStorable", "deepCopyActions", "generics" ];

# List of classes already output.
our @outputPrevious;

sub Process_ClassDocumentation {
    # Skip this if not a documentation build.
    return
	unless ( exists($ENV{'GALACTICUS_BUILD_DOCS'}) && $ENV{'GALACTICUS_BUILD_DOCS'} eq "yes" );
    # Get the tree.
    my $tree  = shift();
    my @trees = ( $tree );
    # For the "objects.nodes.F90" the type-bound functions are currently accessed via an include file. We need to parse that file
    # (and a few others) here so that we can find those functions. This is an ugly hack - ideally the type-bound functions would
    # be pre-processed into the objects.nodes.F90 file instead.
    if ( $tree->{'name'} eq "objects.nodes.F90" ) {
	my @includeFiles = ( "objects.nodes.components.Inc" );
	my $xml = new XML::Simple();
	my $directiveLocations = $xml->XMLin($ENV{'BUILDPATH'}."/directiveLocations.xml");
	foreach my $componentFile ( @{$directiveLocations->{'component'}->{'file'}} ) {
	    push(
		@includeFiles,
		 map
		{exists($_->{'functions'}) && $_->{'functions'} =~ s/\.inc/.Inc/ ? $_->{'functions'} : ()}
		&Galacticus::Build::Directives::Extract_Directives($componentFile,"component")
		);
	}
	foreach my $includeFile ( @includeFiles ) {
	    my $includeTree = &Galacticus::Build::SourceTree::ParseFile($ENV{'BUILDPATH'}."/".$includeFile);
	    push(@trees,$includeTree);
	}
    }
    # Walk the trees.
    my $classes;
    my @functionList;
    foreach my $tree_ ( @trees ) { 
	my $node  = $tree_;
	my $depth = 0;
	while ( $node ) {
	    if ( $node->{'type'} eq "type" ) {
		my $className = $node->{'name'      };
		my $child     = $node->{'firstChild'};
		my $contained = 0;
		while ( $child ) {
		    if ( $child->{'type'} eq "contains" ) {
			$child = $child->{'firstChild'};
			last
			    unless ( $child );
			$contained = 1;
		    }
		    if ( $contained ) {
			if ( $child->{'type'} eq "methods" ) {
			    # Method directive - store a pointer to the directive.
			    push(@{$classes->{$className}->{'descriptions'}},&List::ExtraUtils::as_array($child->{'directive'}->{'method'}));
			} elsif ( $child->{'type'} eq "declaration" ) {
			    # Type-bound function declarations - store a pointer to these.
			    push(@{$classes->{$className}->{'functions'}},$child->{'declarations'});
			}
		    }
		    $child = $child->{'sibling'};
		}
		# Store extension information.
		$classes->{$className}->{'name'} = $className;
		if ( $node->{'opener'} =~ m/,\s*extends\s*\(\s*([a-zA-Z0-9_]+)\s*\)/ ) {
		    $classes->{$className}->{'extends'} = $1;
		}
		# Build a list of methods for this class, along with their associated procedures.
		if ( exists($classes->{$className}->{'descriptions'}) ) {
		    foreach my $method ( @{$classes->{$className}->{'descriptions'}} ) {
			foreach my $functions ( @{$classes->{$className}->{'functions'}} ) {
			    foreach my $function ( &List::ExtraUtils::as_array($functions) ) {
				(my $variableName = $function->{'variables'}->[0]) =~ s/\s*=>.*//;
				if ( $variableName eq lc($method->{'method'}) ) {
				    # This is the matching function. Look for the corresponding function or interface name.
				    if ( grep {$_ eq "deferred"} @{$function->{'attributes'}} ) {
					# A deferred function - there should be an interface defined.
					$method->{'interface'} = $function->{'type'};
				    } else {
					if ( $function->{'variables'}->[0] =~ m/\s*=>\s*([a-zA-Z0-9_,]+)/ ) {
					    foreach my $variable ( @{$function->{'variables'}} ) {
						(my $functionName = $variable) =~ s/\s*([a-zA-Z0-9_\.\(\)\+\-\*\/=<>]+)=>//;					  
						# For generic bindings, check if the bound object is another method. If so, replace
						# the bound name with that of the associated method.
						if ( $function->{'intrinsic'} eq "generic" ) {
						    foreach my $siblingFunctions ( @{$classes->{$className}->{'functions'}} ) {
							foreach my $siblingFunction ( &List::ExtraUtils::as_array($siblingFunctions) ) {
							    if ( $siblingFunction->{'variables'}->[0] =~ m/(\S+)\s*=>\s*(\S+)/ ) {
								my $methodName    = $1;
								my $boundFunction = $2;
								if ( lc($methodName) eq lc($functionName) ) {
								    $functionName = $boundFunction;
								}
							    }
							}
						    }
						}
						push(@{$method->{'boundFunctions'}},$functionName);
					    }
					}
				    }
				}
			    }
			}
		    }
		}
		# Build a list of methods with no matching description.
		my @missingMethods;
		my @genericUses;
		foreach my $functions ( @{$classes->{$className}->{'functions'}} ) {
		    foreach my $function ( &List::ExtraUtils::as_array($functions) ) {
			(my $functionName = $function->{'variables'}->[0]) =~ s/=>.*//;
			unless ( $function->{'intrinsic'} eq "final" || $function->{'variables'}->[0] !~ m/=>/ || grep {lc($_->{'method'}) eq lc($functionName)} @{$classes->{$className}->{'descriptions'}} ) {
			    # Method has no description. Check if it is used by a generic.
			    my $genericUse = 0;
			    foreach my $functionsGeneric ( @{$classes->{$className}->{'functions'}} ) {
				foreach my $functionGeneric ( &List::ExtraUtils::as_array($functions) ) {
				    next
					unless ( $functionGeneric->{'intrinsic'} eq "generic" );
				    my @genericNames = map {$_ =~ s/.*=>\s*//; $_} @{$functionGeneric->{'variables'}};
				    $genericUse = 1
					if ( grep {$_ eq lc($functionName)} @genericNames );
				}
			    }
			    if ( $genericUse ) {
				push(@genericUses   ,$functionName);
			    } else {
				push(@missingMethods,$functionName);
			    }
			}
		    }
		}
		$classes->{$className}->{'missingMethods'} = join(" ",@missingMethods);
		$classes->{$className}->{'genericUses'   } = join(" ",@genericUses   );
	    } elsif ( $node->{'type'} eq "subroutine" || $node->{'type'} eq "function" ) {
		# A function - push it onto the list to be processed later.
		push(@functionList,$node);
	    }
	    # Walk to the next node in the tree.
	    $node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
	}
    }
    # Process all functions.
    foreach my $function ( @functionList ) {
	&processFunction($function,$classes);
    }
    # Remove content not for output.
    foreach my $class ( &List::ExtraUtils::hashList($classes) ) {
	delete($class->{'functions'})
	    if ( exists($class->{'functions'}) );
	delete($classes->{$class->{'name'}})
	    if ( @outputPrevious && grep {$_ eq $class->{'name'}} @outputPrevious );
    }
    # Output the gathered information.
    if ( $classes && $tree->{'name'} =~ m/(.+)\.(F90|Inc)/ && $tree->{'name'} ne "objects.nodes.components.Inc" ) {
	my $classFileName = $ENV{'BUILDPATH'}."/".$1.".classes.xml";
	my $xml = new XML::Simple();
	open(my $classFile,">",$classFileName);
	print $classFile $xml->XMLout($classes,RootName => "classes",NoAttr => 1);
	close($classFile);
	push(@outputPrevious,keys(%{$classes}));
    }
}

sub processFunction {
    # Process a function, extracting argument declarations.
    my $function = shift();
    my $classes  = shift();
    my $functionName = $function->{'name'};
    my @functionArguments;
    if ( my @matches = ($function->{'opener'} =~ $Fortran::Utils::unitOpeners{$function->{'type'}}->{'regEx'}) ) {
	@functionArguments = map {$_ eq "self" ? () : $_} split(/\s*,\s*/,$matches[$Fortran::Utils::unitOpeners{$function->{'type'}}->{'arguments'}])
	    if ( defined($matches[$Fortran::Utils::unitOpeners{$function->{'type'}}->{'arguments'}]) );
    }
    ## Iterate over discovered classes.
    foreach my $class ( &List::ExtraUtils::hashList($classes) ) {
	## Iterate over discovered methods.
	foreach my $method ( @{$class->{'descriptions'}} ) {
	    if ( ( exists($method->{'interface'}) && lc($method->{'interface'}) eq lc($functionName) ) || ( exists($method->{'boundFunctions'}) && grep {lc($_) eq lc($functionName)} @{$method->{'boundFunctions'}} ) ) {
		$method->{'type'} = $function->{'type'};
		push(@{$method->{'argumentList'}},join(", ",@functionArguments));
		my $functionReturnName;
		if ( $method->{'type'} eq "function" ) {
		    if ( $function->{'opener'} =~ m/\sresult\s*\(\s*([a-zA-Z0-9_]+)\s*\)/ ) {
			$functionReturnName = $1;
		    } else {
			$functionReturnName = $function->{'name'};
		    }
		}
		my $functionContent = $function->{'firstChild'};
		my @argumentDeclarations;
		while ( $functionContent ) {
		    if ( $functionContent->{'type'} eq "declaration" ) {
			# Look for arguments.
			foreach my $declaration ( @{$functionContent->{'declarations'}} ) {
			    my @matchedNames;
			    foreach my $argumentName ( @functionArguments ) {
				if ( grep {$_ eq $argumentName} @{$declaration->{'variableNames'}} ) {
				    push(@matchedNames,$argumentName);
				}
			    }
			    if ( @matchedNames ) {
				my $declarationCopy = dclone($declaration);
				@{$declarationCopy->{'variableNames'}} = @matchedNames;
				@{$declarationCopy->{'variables'    }} = @matchedNames;
				push(@argumentDeclarations,$declarationCopy);
			    }
			    # Look for a match to the function return type.
			    if ( $functionReturnName && grep {lc($_) eq lc($functionReturnName)} @{$declaration->{'variableNames'}} ) {
				$method->{'type'} = dclone($declaration);
				$method->{'type'}->{'variableNames'} = [ $functionReturnName ];
				$method->{'type'}->{'variables'    } = [ $functionReturnName ];
			    }
			}			
		    }
		    $functionContent = $functionContent->{'sibling'};
		}
		my $argument;
		$argument->{'argument'} = \@argumentDeclarations;
		push(@{$method->{'arguments'}},$argument);
		# If this is a function, and no return type has been found, determine it from the opener.
		if ( $method->{'type'} eq "function" ) {
		    if ( my @matches = ( $function->{'opener'} =~ $Fortran::Utils::unitOpeners{'function'}->{'regEx'}) ) {
			my $intrinsic = $matches[$Fortran::Utils::unitOpeners{'function'}->{'intrinsic'}];
			my $kind      = $matches[$Fortran::Utils::unitOpeners{'function'}->{'kind'     }];
			delete($method->{'type'});
			$method->{'type'}->{'intrinsic'} = $intrinsic;
			$method->{'type'}->{'type'     } = $kind     ;
		    } else {
			die(" Galacticus::Build::SourceTree::Process::ClassDocumentation::processFunction(): can not determine return type for function");
		    }
		}
	    }
	}
    }
}

1;
