# Contains a Perl module which implements processing of functionClass directives.

package Galacticus::Build::SourceTree::Process::FunctionClass;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Data::Dumper;
use XML::Simple;
use Sort::Topological qw(toposort);
use LaTeX::Encode;
use Scalar::Util qw(reftype);
use List::ExtraUtils;
use Fortran::Utils;
use Text::Template 'fill_in_string';
use Galacticus::Build::SourceTree::Process::SourceIntrospection;

# Insert hooks for our functions.
$Galacticus::Build::SourceTree::Hooks::processHooks{'functionClass'} = \&Process_FunctionClass;

sub Process_FunctionClass {
    # Get the tree.
    my $tree = shift();
    # Get an XML parser.
    my $xml = new XML::Simple();
    # Initialize code directive locations.
    my $directiveLocations;
    # Initialize state storables database.
    my $stateStorables;
    # Initialize deep copy actions database.
    my $deepCopyActions;
    # Determine if debugging output is required.
    my $debugging = exists($ENV{'GALACTICUS_OBJECTS_DEBUG'}) && $ENV{'GALACTICUS_OBJECTS_DEBUG'} eq "yes";
    # Walk the tree, looking for code blocks.
    my $node  = $tree;
    my $depth = 0;
    while ( $node ) {
	if ( $node->{'type'} eq "functionClass" ) {
	    # Assert that our parent is a module.
	    die("Process_FunctionClass: parent node must be a module")
		unless ( $node->{'parent'}->{'type'} eq "module" );
	    my $lineNumber = $node->{'line'};
	    # Extract the directive.
	    my $directive = $node->{'directive'};
	    # Get code directive locations if we do not have them.
	    $directiveLocations = $xml->XMLin($ENV{'BUILDPATH'}."/directiveLocations.xml")
		unless ( $directiveLocations );	    
	    # Get state storables database if we do not have it.
	    $stateStorables = $xml->XMLin($ENV{'BUILDPATH'}."/stateStorables.xml")
		unless ( $stateStorables );
	    # Get state storables database if we do not have it.
	    $deepCopyActions = $xml->XMLin($ENV{'BUILDPATH'}."/deepCopyActions.xml")
		unless ( $deepCopyActions );
	    # Find methods.
	    my %methods;
	    if ( exists($directive->{'method'}) ) {
		if ( exists($directive->{'method'}->{'name'}) && ! reftype($directive->{'method'}->{'name'}) ) {
		    $methods{$directive->{'method'}->{'name'}} = $directive->{'method'};
		} else {
		    %methods = %{$directive->{'method'}};
		}
	    }
	    # Load any functionClassType that the base class extends.
	    my $functionClassType;
	    if ( exists($directive->{'extends'}) ) {
		(my $functionClassTypeFileName) = map {$_->{'name'} eq $directive->{'extends'} ? $_->{'file'} : ()} &List::ExtraUtils::as_array($stateStorables->{'functionClassTypes'});
		die('failed to find file containing functionClassType "'.$directive->{'extends'}.'"')
		    unless ( defined($functionClassTypeFileName) );
		$functionClassType->{'tree'} = &Galacticus::Build::SourceTree::ParseFile($functionClassTypeFileName);
		my $classNode  = $functionClassType->{'tree'};
		my $classDepth = 0;
		while ( $classNode ) {
		    if ( $classNode->{'type'} eq "type" ) {
			if ( $classNode->{'name'} eq $directive->{'extends'} ) {
			    $functionClassType->{'node'} = $classNode;
			}
		    }
		    $classNode = &Galacticus::Build::SourceTree::Walk_Tree($classNode,\$classDepth);
		}
	    }
	    # Find class locations.
	    my @classLocations = &List::ExtraUtils::as_array($directiveLocations->{$directive->{'name'}}->{'file'})
		if ( exists($directiveLocations->{$directive->{'name'}}) );
	    # Parse classes.
	    my %dependencies;
	    my %classes;
	    foreach my $classLocation ( @classLocations ) {
		my $classTree  = &Galacticus::Build::SourceTree::ParseFile($classLocation);
		&Galacticus::Build::SourceTree::ProcessTree($classTree, errorTolerant => 1);
		my $classNode  = $classTree;
		my $classDepth = 0;
		my $className;
		my $class;
		$class->{'file'} = $classLocation;
		while ( $classNode ) {
		    # Collect class directives.
		    if ( $classNode->{'type'} eq $directive->{'name'} ) {
			$class->{'node'} = $classNode;
			$class->{$_    } = $classNode->{'directive'}->{$_}
			    foreach ( keys(%{$classNode->{'directive'}}) );
		    }
		    if ( $classNode->{'type'} eq "type" ) {
			# Parse class openers to find dependencies.
			if (
			    $classNode->{'opener'} =~ m/^\s*type\s*(,\s*abstract\s*|,\s*public\s*|,\s*private\s*|,\s*extends\s*\(([a-zA-Z0-9_]+)\)\s*)*(::)??\s*$directive->{'name'}([a-z0-9_]+)\s*$/i 
			    &&
			    defined($2)
			    ) {
			    my $extends = $2;
			    my $type    = $directive->{'name'}.$4;
			    $class->{'type'   } = $type;
			    $class->{'extends'} = $extends;
			    push(@{$dependencies{$extends}},$type);
			    # Also determine if any other members of this class are used in this type definition, and add suitable dependencies.
			    my $childNode = $classNode->{'firstChild'};
			    while ( $childNode ) {
				if ( $childNode->{'type'} eq "declaration" ) {
				    foreach my $declaration ( @{$childNode->{'declarations'}} ) {
					push(@{$dependencies{$declaration->{'type'}}},$type)
					    if ( ( $declaration->{'intrinsic'} eq "class" || $declaration->{'intrinsic'} eq "type" ) && $declaration->{'type'} =~ m/^$directive->{'name'}/ );
				    }
				}
				$childNode = $childNode->{'sibling'};
			    }
			}
		    }
		    $classNode = &Galacticus::Build::SourceTree::Walk_Tree($classNode,\$classDepth);
		}
		# Store tree.
		$class->{'tree'} = $classTree;
		# Set defaults.
		$class->{'abstract'} = "no"
		    unless ( exists($class->{'abstract'}) );
	        # Store to set of all classes.
		die('Galacticus::Build::SourceTree::Process::FunctionClass::Process_FunctionClass: class is undefined in file "'.$classLocation.'"')
		    unless ( defined($class->{'type'}) );
		$classes{$class->{'type'}} = $class;
	    }
	    # Sort classes.
	    my @unsortedClasses = keys(%classes);
	    my @sortedClasses   = toposort(sub { @{$dependencies{$_[0]} || []}; }, \@unsortedClasses );
	    my @classes         = map($classes{$_},@sortedClasses);
	    # Create a set of non-abstract classes.
	    my @nonAbstractClasses;
	    foreach ( @classes ) {
		push(
		    @nonAbstractClasses,
		    $_
		    )
		    unless ( exists($_->{'abstract'}) && $_->{'abstract'} eq "yes" );
	    }
	    # Add methods to store and retrieve state.
	    $methods{'stateStore'} =
	    {
		description => "Store the state of the object to file.",
		type        => "void",
		pass        => "yes",
		modules     =>
		    [
		     {
			 name => "FGSL",
			 only => [ "fgsl_file" ]
		     }
		    ],
		argument    => [ "integer, intent(in   ) :: stateFile", "type(fgsl_file), intent(in   ) :: fgslStateFile" ],
	    };
	    $methods{'stateRestore'} =
	    {
		description => "Restore the state of the object to file.",
		type        => "void",
		pass        => "yes",
		modules     => 
		    [
		     {
			 name => "FGSL",
			 only => [ "fgsl_file" ]
		     }
		    ],
		argument    => [ "integer, intent(in   ) :: stateFile", "type(fgsl_file), intent(in   ) :: fgslStateFile" ],
	    };
	    # If the function requires calculation reset, add method to do so.
	    if ( exists($directive->{'calculationReset'}) && $directive->{'calculationReset'} eq "yes" ) {
		$methods{'calculationReset'} =
		{
		    description => "Reset the calculation state of the object.",
		    type        => "void",
		    pass        => "yes",
		    argument    => [ "type(treeNode), intent(inout) :: node" ],
		    code        => join("",map {"if (sizeof(".$_.")<0.and.sizeof(".$_.")>0) then\nend if\n"} ('self','node') )
		};
	    }
	    # Add auto-hook function.
	    $methods{'autoHook'} = 
	    {
		description => "Insert any event hooks required by this object.",
		type        => "void",
		pass        => "yes",
		code        => "!GCC\$ attributes unused :: self\n\n! Nothing to do by default.\n"
	    };
	    # Add "descriptor" method.
	    my $descriptorCode;
	    my %descriptorModules = ( "Input_Parameters" => 1 );
	    my %addSubParameters;
	    my $addLabel         = 0;
	    my $descriptorUsed   = 0;
	    $descriptorCode .= "select type (self)\n";
	    foreach my $nonAbstractClass ( @nonAbstractClasses ) {
		(my $label = $nonAbstractClass->{'name'}) =~ s/^$directive->{'name'}//;
		$label = lcfirst($label)
		    unless ( $label =~ m/^[A-Z]{2,}/ );
		my $hasCustomDescriptor = 0;
		my $extensionOf;
		# Build lists of all potential parameter and object names for this class, including any from parent classes.
		my $potentialNames;
		my $class = $nonAbstractClass;
		while ( $class ) {
		    my $node = $class->{'tree'}->{'firstChild'};
		    $node = $node->{'sibling'}
		        while ( $node && ( $node->{'type'} ne "type" || ( ! exists($node->{'name'}) || $node->{'name'} ne $class->{'name'} ) ) );
		    last
			unless ( $node );
		    # Find the parent class.
		    if ( $class == $nonAbstractClass && $node->{'opener'} =~ m/,\s*extends\s*\(\s*([a-zA-Z0-9_]+)\s*\)/ ) {
			$extensionOf = $1;
		    }
		    # Search the node for declarations.
		    $node = $node->{'firstChild'};
		    while ( $node ) {
			if ( $node->{'type'} eq "declaration" ) {
			    foreach my $declaration ( @{$node->{'declarations'}} ) {
				# Identify object pointers.
				push(@{$potentialNames->{'objects'}},map {$_ =~ s/\s*([a-zA-Z0-9_]+).*/$1/; $_} @{$declaration->{'variables'}})
				    if
				    (
				     $declaration->{'intrinsic'} eq "class"
				     &&
				     $declaration->{'type'     } =~ m/Class\s*$/
				     &&
				     grep {$_ eq "pointer"} @{$declaration->{'attributes'}}
				    );
				push(@{$potentialNames->{'parameters'}},$declaration)
				    if
				    (
				     (grep {$_ eq $declaration->{'intrinsic'}} ( "integer", "logical", "double precision" ))
				     ||
				     (
				             $declaration->{'intrinsic'}  eq "type"
				      &&
				      trimlc($declaration->{'type'     }) eq "varying_string"				      
				     )
				    );
				$hasCustomDescriptor = 1
				    if
				    (
				     $declaration->{'intrinsic'} eq "procedure"
				     &&
				     $declaration->{'variables'}->[0] =~ m/^descriptor=>/
				    );
			    }
			}
			$node = $node->{'type'} eq "contains" ? $node->{'firstChild'} : $node->{'sibling'};
		    }
		    # Move to the parent class.
		    $class = ($class->{'extends'} eq $directive->{'name'}) ? undef() : $classes{$class->{'extends'}};
		}		
		# Add any names declared in the base class.
		foreach my $data ( &List::ExtraUtils::as_array($directive->{'data'}) ) {
		    my $declarationSource;
		    if ( reftype($data) ) {
			$declarationSource = $data->{'content'}
			    if ( $data->{'scope'} eq "self" );
		    } else {
			$declarationSource = $data;
		    }
		    next
			unless ( defined($declarationSource) );
		    my $declaration = &Fortran::Utils::Unformat_Variables($declarationSource);
		    die("Galacticus::Build::SourceTree::Process::FunctionClass::Process_FunctionClass(): unable to parse variable declaration")
			unless ( defined($declaration) );
		    push(@{$potentialNames->{'objects'}},map {$_ =~ s/\s*([a-zA-Z0-9_]+).*/$1/; $_} @{$declaration->{'variables'}})
			if
			(
			 $declaration->{'intrinsic'} eq "class"
			 &&
			 $declaration->{'type'     } =~ m/Class\s*$/
			 &&
			 grep {$_ eq "pointer"} @{$declaration->{'attributes'}}
			);
		    push(@{$potentialNames->{'parameters'}},$declaration)
			if
			(
			 (grep {$_ eq $declaration->{'intrinsic'}} ( "integer", "logical", "double precision" ))
			 ||
			 (
			         $declaration->{'intrinsic'}  eq "type"
			  &&
			  trimlc($declaration->{'type'     }) eq "varying_string"				      
			 )
			);
		}
		# Add any names declared in the functionClassType.
		if ( defined($functionClassType) ) {
		    # Search the node for declarations.
		    my $node = $functionClassType->{'node'}->{'firstChild'};
		    while ( $node ) {
			if ( $node->{'type'} eq "declaration" ) {
			    foreach my $declaration ( @{$node->{'declarations'}} ) {
				# Identify object pointers.
				push(@{$potentialNames->{'objects'}},map {$_ =~ s/\s*([a-zA-Z0-9_]+).*/$1/; $_} @{$declaration->{'variables'}})
				    if
				    (
				     $declaration->{'intrinsic'} eq "class"
				     &&
				     $declaration->{'type'     } =~ m/Class\s*$/
				     &&
				     grep {$_ eq "pointer"} @{$declaration->{'attributes'}}
				    );
				push(@{$potentialNames->{'parameters'}},$declaration)
				    if
				    (
				     (grep {$_ eq $declaration->{'intrinsic'}} ( "integer", "logical", "double precision" ))
				     ||
				     (
				             $declaration->{'intrinsic'}  eq "type"
				      &&
				      trimlc($declaration->{'type'     }) eq "varying_string"				      
				     )
				    );
			    }
			}
			$node = $node->{'type'} eq "contains" ? $node->{'firstChild'} : $node->{'sibling'};
		    }
		}		
		# Search the tree for this class to find the interface to the parameters constructor.
		my $node = $nonAbstractClass->{'tree'}->{'firstChild'};
		$node = $node->{'sibling'}
		    while ( $node && ( $node->{'type'} ne "interface" || ( ! exists($node->{'name'}) || $node->{'name'} ne $nonAbstractClass->{'name'} ) ) );
		next
		    unless ( $node );
		# Find all constructor names.
		$node = $node->{'firstChild'};		
		my @constructors;
		while ( $node ) {
		    push(@constructors,@{$node->{'names'}})
			if ( $node->{'type'} eq "moduleProcedure" );
		    $node = $node->{'sibling'};
		}
		# Search for constructors.
		$node = $nonAbstractClass->{'tree'}->{'firstChild'};
		my $descriptorParameters;
		my %subParameters;
		my $declarationMatches    = 0;
		my $supported             = 1;
		my $parentConstructorUsed = 0;
		my @failureMessage;
		while ( $node ) {
		    if ( $node->{'type'} eq "function" && (grep {$_ eq $node->{'name'}} @constructors) && $node->{'opener'} =~ m/^\s*function\s+$node->{'name'}\s*\(\s*parameters\s*\)/ ) {
			# Extract the name of the return variable in this function.
			my $result = ($node->{'opener'} =~ m/result\s*\(\s*([a-zA-Z0-9_]+)\s*\)\s*$/) ? $1 : $node->{'name'};
			# Check if this is the parameters constructor.
			my $constructorNode    = $node->{'firstChild'};
			my $depth = 0;
			while ( $constructorNode ) {
			    # Process node.
			    if ( $constructorNode->{'type'} eq "declaration" ) {
				# Declaration node found - check if we have a parameters argument of the correct type.
				foreach my $declaration ( @{$constructorNode->{'declarations'}} ) {
				    $declarationMatches = 1
					if ( 
					    $declaration->{'intrinsic'}  eq "type"
					    &&
					    trimlc($declaration->{'type'     }) eq "inputparameters"
					    &&
					    grep {$_ eq "parameters"} @{$declaration->{'variables'}}
					);	     
				}
			    }
			    if ( $constructorNode->{'type'} eq "code" ) {
				# Locate any use of sub-parameters and of the parent class constructor.
				open(my $code,"<",\$constructorNode->{'content'});
				do {
				    # Get a line.
				    &Fortran::Utils::Get_Fortran_Line($code,my $rawLine, my $processedLine, my $bufferedComments); 
				    # Identify subparameter usages.
				    if ( $processedLine =~ m/^\s*([a-zA-Z0-9_]+)\s*=\s*([a-zA-Z0-9_]+)\s*\%\s*subParameters\s*\(/ ) {
					$subParameters{$1} =
					{
					    parent => $2,
					    source => $processedLine
					};
				    }
				    # Identify use of parent constructor.
				    $parentConstructorUsed = 1
					if ( $processedLine =~ m/^\s*$result\s*\%\s*$extensionOf\s*=/ );
				} until ( eof($code) );
				close($code);
			    }
			    if ( $constructorNode->{'type'} eq "inputParameter" ) {
				if ( exists($constructorNode->{'directive'}->{'source'}) ) {
				    if      ( exists($constructorNode->{'directive'}->{'name'    }) ) {
					# A regular parameter, defined by its name.
					my $name;
					if ( exists($constructorNode->{'directive'}->{'variable'}) ) {
					    ($name = $constructorNode->{'directive'}->{'variable'}) =~ s/.*\%(.*)/$1/;
					} else {
					    $name = $constructorNode->{'directive'}->{'name'};
					}
					if ( grep {$_ eq lc($name)} (map {@{$_->{'variables'}}} @{$potentialNames->{'parameters'}}) ) {
					    push(@{$descriptorParameters->{'parameters'}},{name => $name, inputName => $constructorNode->{'directive'}->{'name'}, source => $constructorNode->{'directive'}->{'source'}});
					    # Find the matched variable.
					    my $descriptor;
					    foreach my $potentialDescriptor ( @{$potentialNames->{'parameters'}} ) {
						$descriptor = $potentialDescriptor
						    if ( grep {$_ eq lc($name)} @{$potentialDescriptor->{'variables'}} );
					    }					    
					    if ( grep {$_ =~ m/^dimension\s*\([a-z0-9_:,\s]+\)/} @{$descriptor->{'attributes'}} ) {
						# A non-scalar parameter - currently not supported.
						$supported = -8;
						push(@failureMessage,"non-scalar parameters not supported");
					    }
					} else {
					    $supported = -1;
					    push(@failureMessage,"could not find a matching internal variable for parameter [".$name."]");
					}
				    } elsif ( exists($constructorNode->{'directive'}->{'regEx'   }) ) {
					# A regular expression parameter. Currently not supported.
					$supported = -2;
					push(@failureMessage,"regular expression parameter [".$constructorNode->{'directive'}->{'regEx'}."] not supported");
				    } elsif ( exists($constructorNode->{'directive'}->{'iterator'}) ) {
					# A parameter whose name iterates over a set of possible names. Currently not supported.
					$supported = -3;
					push(@failureMessage,"iterator parameter [".$constructorNode->{'directive'}->{'iterator'}."] not supported");
				    }
				} else {
				    $supported = -4;
				    push(@failureMessage,"unsourced parameters not supported");
				}
			    }
			    if ( $constructorNode->{'type'} eq "objectBuilder"  ) {		    
				if ( exists($constructorNode->{'directive'}->{'source'}) ) {				    
				    (my $name = $constructorNode->{'directive'}->{'name'}) =~ s/([a-zA-Z0-9_]+\s*\%\s*)?([a-zA-Z0-9_]+).*/$2/;
				    $name =~ s/\s//g;
				    if ( grep {$_ eq lc($name)} @{$potentialNames->{'objects'}} ) { 
					push(@{$descriptorParameters->{'objects'}},{name => $name, source => $constructorNode->{'directive'}->{'source'}});
				    } else {
					$supported = -5;
					push(@failureMessage,"could not find a matching internal object for object [".$name."]");
				    }
				} else {
				    $supported = -6;
				    push(@failureMessage,"unsourced objects not supported");
				}
			    }
			    $constructorNode = &Galacticus::Build::SourceTree::Walk_Tree($constructorNode,\$depth);
			    last
				if ( $depth < 0 );
			}
		    }
		    $node = $node->{'type'} eq "contains" ? $node->{'firstChild'} : $node->{'sibling'};
		}
		# Validate sub-parameters.
		foreach my $subParameterName ( keys(%subParameters) ) {
		    unless ( exists($subParameters{$subParameters{$subParameterName}->{'parent'}}) || $subParameters{$subParameterName}->{'parent'} eq "parameters" ) {
			$supported = -7;
			push(@failureMessage,"subparameter hierarchy failure");
		    }
		}
		# Build the code.
		$descriptorCode .= "type is (".$nonAbstractClass->{'name'}.")\n";
		if ( $hasCustomDescriptor ) {
		    # The class has its own descriptor function, so we should never arrive at this point in the code.
		    $descriptorCode .= " call Galacticus_Error_Report('custom descriptor exists - this should not happen'//".&Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($nonAbstractClass->{'node'},$nonAbstractClass->{'node'}->{'line'}).")\n";
		    $descriptorModules{'Galacticus_Error'} = 1;
		} else{
		    # Build an auto-descriptor function.
		    if ( $declarationMatches && $supported == 1 ) {
			$descriptorUsed = 1;
			$descriptorCode .= " if (.not.present(includeMethod).or.includeMethod) call descriptor%addParameter('".$directive->{'name'}."Method','".$label."')\n";
			if ( defined($descriptorParameters) ) {			    
			    # Get subparameters.
			    $addSubParameters{'parameters'} = 1;
			    $descriptorCode   .= "parameters=descriptor%subparameters('".$directive->{'name'}."Method')\n";
			    foreach my $subParameterName (keys(%subParameters) ) {
				$addSubParameters{$subParameterName} = 1;
				$descriptorCode .= $subParameters{$subParameterName}->{'source'};
			    }
			    # Handle parameters set via inputParameter directives.
			    if ( defined($descriptorParameters->{'parameters'}) ) {
				foreach my $parameter ( @{$descriptorParameters->{'parameters'}} ) {
				    foreach my $declaration ( @{$potentialNames->{'parameters'}} ) {
					if ( grep {$_ eq lc($parameter->{'name'})} @{$declaration->{'variables'}} ) {
					    if      ( $declaration->{'intrinsic'} eq "type" ) {
						$descriptorCode .= "call ".$parameter->{'source'}."%addParameter('".$parameter->{'inputName'}."',char(self%".$parameter->{'name'}."))\n";
					    } elsif ( $declaration->{'intrinsic'} eq "logical" ) {
						$descriptorCode .= "if (self%".$parameter->{'name'}.") then\n";
						$descriptorCode .= "  call ".$parameter->{'source'}."%addParameter('".$parameter->{'inputName'}."','true' )\n";
						$descriptorCode .= "else\n";
						$descriptorCode .= "  call ".$parameter->{'source'}."%addParameter('".$parameter->{'inputName'}."','false')\n";
						$descriptorCode .= "end if\n";
					    } else {
						$addLabel = 1;
						my $format;
						$format = "e17.10"
						    if ( $declaration->{'intrinsic'} eq "double precision" );
						$format = "i17"
						    if ( $declaration->{'intrinsic'} eq "integer"          );
						$descriptorCode .= "write (parameterLabel,'(".$format.")') self%".$parameter->{'name'}."\n";
						$descriptorCode .= "call ".$parameter->{'source'}."%addParameter('".$parameter->{'inputName'}."',trim(adjustl(parameterLabel)))\n";
					    }
					}
				    }
				}
			    }
			    # Handle objects built via objectBuilder directives.
			    if ( defined($descriptorParameters->{'objects'}) ) {
				foreach ( @{$descriptorParameters->{'objects'}} ) {
				    $descriptorCode .= "call self%".$_->{'name'}."%descriptor(".$_->{'source'}.")\n";
				}
			    }
			}
			# If the parent constructor was used, call its descriptor method.
			if ( $parentConstructorUsed ) {
			    $descriptorCode .= "call self%".$extensionOf."%descriptor(descriptor,includeMethod=.false.)\n";
			}
		    } elsif ( ! $declarationMatches     ) {		    
			$descriptorCode .= " call Galacticus_Error_Report('auto-descriptor not supported for this class: parameter-based constructor not found'//".&Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($nonAbstractClass->{'node'},$nonAbstractClass->{'node'}->{'line'}).")\n";
			$descriptorModules{'Galacticus_Error'} = 1;
		    } elsif (   $supported         != 1 ) {
			$descriptorCode .= " call Galacticus_Error_Report('auto-descriptor not supported for this class because:'//char(10)//".join("//char(10)// &\n & ",map {"'  --> ".$_."'"} @failureMessage)."//".&Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($nonAbstractClass->{'node'},$nonAbstractClass->{'node'}->{'line'}).")\n";
			$descriptorModules{'Galacticus_Error'} = 1;
		    }
		}
	    }	    
	    $descriptorCode .= "end select\n";
	    $descriptorCode  = " !GCC\$ attributes unused :: descriptor, includeMethod\n".$descriptorCode
		unless ( $descriptorUsed );
 	    $descriptorCode  = "type(inputParameters) :: ".join(",",keys(%addSubParameters))."\n".$descriptorCode
		if ( %addSubParameters );
 	    $descriptorCode  = "character(len=18) :: parameterLabel\n".$descriptorCode
		if ( $addLabel );
	    $methods{'descriptor'} = 
	    {
		description => "Return an input parameter list descriptor which could be used to recreate this object.",
		type        => "void",
		pass        => "yes",
		modules     => join(" ",keys(%descriptorModules)),
		argument    => [ "type(inputParameters), intent(inout) :: descriptor", "logical, intent(in   ), optional :: includeMethod" ],
		code        => $descriptorCode
	    };
	    # Add a "hashedDescriptor" method.
	    $code::directiveName = $directive->{'name'};
	    my $hashedDescriptorCode = fill_in_string(<<'CODE', PACKAGE => 'code');
type(inputParameters)       :: descriptor
type(varying_string )       :: descriptorString
type(varying_string ), save :: descriptorStringPrevious, hashedDescriptorPrevious
!$omp threadprivate(descriptorStringPrevious,hashedDescriptorPrevious)
descriptor=inputParameters()
! Disable live nodeLists in FoX as updating these nodeLists leads to memory leaks.
call setLiveNodeLists(descriptor%document,.false.)
call self%descriptor(descriptor)
descriptorString=descriptor%serializeToString()
call descriptor%destroy()
if (present(includeSourceDigest).and.includeSourceDigest) then
select type (self)
CODE
	    foreach my $nonAbstractClass ( @nonAbstractClasses ) {
		$code::type = $nonAbstractClass->{'name'};
		(my $classFile = $tree->{'source'}) =~ s/^.*\/([^\/]+)$/$1/;
		my @sourceFiles = ( $classFile );
		my $class = $nonAbstractClass;
		while ( $class ) {
		    (my $sourceFile = $class->{'file'}) =~ s/^.*\/([^\/]+)$/$1/;
		    push(@sourceFiles,$sourceFile);
		    $class = ($class->{'extends'} eq $directive->{'name'}) ? undef() : $classes{$class->{'extends'}};
		}
		$code::digest = &Galacticus::Build::SourceTree::Process::SourceDigest::Find_Hash(@sourceFiles);
		$hashedDescriptorCode .= fill_in_string(<<'CODE', PACKAGE => 'code');
type is ({$type})
descriptorString=descriptorString//":sourceDigest\{{$digest}\}"
CODE
	    }
	    $hashedDescriptorCode .= fill_in_string(<<'CODE', PACKAGE => 'code');
end select
end if
if (descriptorString /= descriptorStringPrevious) then
   descriptorStringPrevious=         descriptorString
   hashedDescriptorPrevious=Hash_MD5(descriptorString)
end if
{$directiveName}HashedDescriptor=hashedDescriptorPrevious
CODE
	    $methods{'hashedDescriptor'} = 
	    {
		description => "Return a hash of the descriptor for this object, optionally include the source code digest in the hash.",
		type        => "type(varying_string)",
		pass        => "yes",
		modules     => "ISO_Varying_String Input_Parameters Hashes_Cryptographic FoX_DOM",
		argument    => [ "logical, intent(in   ), optional :: includeSourceDigest" ],
		code        => $hashedDescriptorCode
	    };
	    # Add "allowedParameters" method.
	    my $allowedParametersCode;
	    my $parametersPresent = 0;
	    foreach my $class ( @classes ) {
		(my $label = $class->{'name'}) =~ s/^$directive->{'name'}//;
		$label = lcfirst($label)
		    unless ( $label =~ m/^[A-Z]{2,}/ );
		# Search the tree for this class to find the interface to the parameters constructor.
		my $node = $class->{'tree'}->{'firstChild'};
		$node = $node->{'sibling'}
		    while ( $node && ( $node->{'type'} ne "interface" || ( ! exists($node->{'name'}) || $node->{'name'} ne $class->{'name'} ) ) );
		next
		    unless ( $node );
		# Find all constructor names.
		$node = $node->{'firstChild'};		
		my @constructors;
		while ( $node ) {
		    push(@constructors,@{$node->{'names'}})
			if ( $node->{'type'} eq "moduleProcedure" );
		    $node = $node->{'sibling'};
		}
		# Search for constructors.
		$node = $class->{'tree'}->{'firstChild'};
		my $allowedParameters;
		my $declarationMatches = 0;
		while ( $node ) {
		    if ( $node->{'type'} eq "function" && (grep {$_ eq $node->{'name'}} @constructors) && $node->{'opener'} =~ m/^\s*function\s+$node->{'name'}\s*\(\s*parameters\s*\)/ ) {
			# Extract the name of the return variable in this function.
			my $result = ($node->{'opener'} =~ m/result\s*\(\s*([a-zA-Z0-9_]+)\s*\)\s*$/) ? $1 : $node->{'name'};
			# Check if this is the parameters constructor.
			my $constructorNode    = $node->{'firstChild'};
			my $depth = 0;
			while ( $constructorNode ) {
			    # Process node.
			    if ( $constructorNode->{'type'} eq "declaration" ) {
				# Declaration node found - check if we have a parameters argument of the correct type.
				foreach my $declaration ( @{$constructorNode->{'declarations'}} ) {
				    $declarationMatches = 1
					if ( 
					           $declaration->{'intrinsic'}  eq "type"
					    &&
					    trimlc($declaration->{'type'     }) eq "inputparameters"
					    &&
					    grep {$_ eq "parameters"} @{$declaration->{'variables'}}
					);	     
				}
			    }
			    if ( $constructorNode->{'type'} eq "inputParameter" ) {
				my $source = exists($constructorNode->{'directive'}->{'source'}) ? $constructorNode->{'directive'}->{'source'} : "globalParameters";
				if      ( exists($constructorNode->{'directive'}->{'name'    }) ) {
				    # A regular parameter, defined by its name.
				    push(@{$allowedParameters->{$source}->{'all'}},         $constructorNode->{'directive'}->{'name' });
				} elsif ( exists($constructorNode->{'directive'}->{'regEx'   }) ) {
				    # A regular expression parameter.
				    push(@{$allowedParameters->{$source}->{'all'}},"regEx:".$constructorNode->{'directive'}->{'regEx'});
				} elsif ( exists($constructorNode->{'directive'}->{'iterator'}) ) {
				    # A parameter whose name iterates over a set of possible names.
				    if ( $constructorNode->{'directive'}->{'iterator'} =~ m/\(\#([a-zA-Z0-9]+)\-\>([a-zA-Z0-9]+)\)/ ) {
					my $directiveName = $1;
					my $attributeName = $2;
					die('Process_FunctionClass(): locations not found for directives')
					    unless ( exists($directiveLocations->{$directiveName}) );
					foreach my $fileName ( &List::ExtraUtils::as_array($directiveLocations->{$directiveName}->{'file'}) ) {
					    foreach ( &Galacticus::Build::Directives::Extract_Directives($fileName,$directiveName) ) {
						(my $parameterName = $constructorNode->{'directive'}->{'iterator'}) =~ s/\(\#$directiveName\-\>$attributeName\)/$_->{$attributeName}/;
						push(@{$allowedParameters->{$source}->{'all'}},$parameterName);
					    }
					}
				    } else {
					die('Process_FunctionClass(): nothing to iterate over');
				    }
				}
			    }
			    if ( $constructorNode->{'type'} eq "objectBuilder"  ) {		    
				my $source = exists($constructorNode->{'directive'}->{'source'}) ? $constructorNode->{'directive'}->{'source'} : "globalParameters";
				push(@{$allowedParameters->{$source}->{'all'}},exists($constructorNode->{'directive'}->{'parameterName'}) ? $constructorNode->{'directive'}->{'parameterName'} : $constructorNode->{'directive'}->{'class'}."Method");
				# Check if the class contains a pointer of the expected type and name for this object.
				my $typeNode = $class->{'tree'}->{'firstChild'};
				while ( $typeNode ) {
				    if ( $typeNode->{'type'} eq "type" && lc($typeNode->{'name'}) eq lc($class->{'name'}) ) {
					$typeNode = $typeNode->{'firstChild'};
					while ( $typeNode ) {
					    if ( $typeNode->{'type'} eq "declaration" ) {
						foreach my $declaration ( @{$typeNode->{'declarations'}} ) {
						    if ( 
							$declaration->{'intrinsic'} eq "class"
							&&
							trimlc($declaration->{'type'}) eq trimlc($constructorNode->{'directive'}->{'class'})."class"
							) {
							push(
							     @{$allowedParameters->{$source}->{'objects'}},
							     map {
								  (
								   lc(            $_) eq striplc($constructorNode->{'directive'}->{'name'})
								   ||
								   lc($result."%".$_) eq striplc($constructorNode->{'directive'}->{'name'})
								  )
								  ?
								  $_
								  :
								  ()
							    } @{$declaration->{'variables'}}
							    );
						    }
						}
					    }
					    $typeNode = $typeNode->{'sibling'};  
					}
					last;
				    }
				    $typeNode = $typeNode->{'sibling'};
				}
			    }
			    $constructorNode = &Galacticus::Build::SourceTree::Walk_Tree($constructorNode,\$depth);
			    last
				if ( $depth < 0 );
			}
		    }
		    $node = $node->{'type'} eq "contains" ? $node->{'firstChild'} : $node->{'sibling'};
		}		
		if ( $declarationMatches && defined($allowedParameters) ) {
		    $parametersPresent      = 1;
		    $allowedParametersCode .= "select type (self)\n";		    
		    # Include the class and all parent classes here - in the parent class constructor we want to accept parameters
		    # that are valid in child classes.
		    my $className = $class->{'name'};
		    while ( defined($className) ) {		    
			$allowedParametersCode .= "class is (".$className.")\n";
			foreach my $source ( keys(%{$allowedParameters}) ) {
			    my $parameterCount = scalar(@{$allowedParameters->{$source}->{'all'}});
			    if ( $parameterCount > 0 ) {
				$allowedParametersCode .= "  if (sourceName == '".$source."') then\n";
				$allowedParametersCode .= "    if (allocated(allowedParameters)) then\n";
				$allowedParametersCode .= "      call move_alloc(allowedParameters,allowedParametersTmp)\n";
				$allowedParametersCode .= "      allocate(allowedParameters(size(allowedParametersTmp)+".$parameterCount."))\n";
				$allowedParametersCode .= "      allowedParameters(1:size(allowedParametersTmp))=allowedParametersTmp\n";
				$allowedParametersCode .= "      deallocate(allowedParametersTmp)\n";
				$allowedParametersCode .= "    else\n";
				$allowedParametersCode .= "      allocate(allowedParameters(".$parameterCount."))\n";
				$allowedParametersCode .= "    end if\n";
				# The following is done as a sequence of scalar assignments, instead of assigning a single array
				# using an array constructor, as that approach lead to a memory leak.
				for(my $i=0;$i<$parameterCount;++$i) {
				    $allowedParametersCode .= "    allowedParameters(size(allowedParameters)-".($parameterCount-1-$i).")='".$allowedParameters->{$source}->{'all'}->[$i]."'\n";
				}
				$allowedParametersCode .= "  end if\n";
			    }
			    # Call the allowedParameters() method of any stored obejcts.
			    if ( $className eq $class->{'name'} ) {
				foreach ( @{$allowedParameters->{$source}->{'objects'}} ) {
				    $allowedParametersCode .= "  if (associated(self%".$_.")) call self%".$_."%allowedParameters(allowedParameters,'".$source."')\n";
				}
			    }
			}
			if ( $classes{$className}->{'extends'} eq $directive->{'name'}."Class" ) {
			    undef($className);
			} else {
			    $className = $classes{$className}->{'extends'};
			}
		    }
		    $allowedParametersCode .= "end select\n";		    
		}
	    }
	    if ( $parametersPresent ) {		
		$allowedParametersCode = "type(varying_string), allocatable, dimension(:) :: allowedParametersTmp\n".$allowedParametersCode;
	    } else {
		$allowedParametersCode = "!GCC\$ attributes unused :: self, allowedParameters, sourceName\n";
	    }
	    $methods{'allowedParameters'} = 
	    {
		description => "Return a list of parameter names allowed for this object.",
		type        => "void",
		recursive   => "yes",
		pass        => "yes",
		modules     => "ISO_Varying_String",
		argument    => [ "type(varying_string), dimension(:), allocatable, intent(inout) :: allowedParameters", "character(len=*), intent(in   ) :: sourceName" ],
		code        => $allowedParametersCode
	    };
	    # Add "deepCopy" method.
            my %deepCopyModules;
            if ( $debugging ) {
		$deepCopyModules{'MPI_Utilities'     } = 1;
		$deepCopyModules{'ISO_Varying_String'} = 1;
		$deepCopyModules{'String_Handling'   } = 1;
		$deepCopyModules{'Galacticus_Display'} = 1;
            }
            my $rankMaximum = 0;
	    my $deepCopyCode;
	    $deepCopyCode .= "select type (self)\n";
	    foreach my $nonAbstractClass ( @nonAbstractClasses ) {
		# Search the tree for this class.
		my $class = $nonAbstractClass;
		my $assignments;
		while ( $class ) {
		    my $node = $class->{'tree'}->{'firstChild'};
		    $node = $node->{'sibling'}
		        while ( $node && ( $node->{'type'} ne "type" || ( ! exists($node->{'name'}) || $node->{'name'} ne $class->{'name'} ) ) );
		    last
			unless ( $node );
		    # Search the node for declarations.
		    my @ignore = exists($class->{'deepCopy'}->{'ignore'}) ? split(/\s*,\s*/,$class->{'deepCopy'}->{'ignore'}->{'variables'}) : ();
		    $node = $node->{'firstChild'};
		    while ( $node ) {
			if ( $node->{'type'} eq "declaration" ) {
			    foreach my $declaration ( @{$node->{'declarations'}} ) {
				# Deep copy of functionClass objects.
				(my $type = $declaration->{'type'}) =~ s/(^\s*|\s*$)//g
				    if ( $declaration->{'intrinsic'} eq "class" || $declaration->{'intrinsic'} eq "type" );
				if
				    (
				     $declaration->{'intrinsic'} eq "class"
				     &&
				     (grep {$_ eq $type    } (@{$stateStorables->{'functionClasses'}},@{$stateStorables->{'functionClassInstances'}}))
				     &&
				     grep {$_ eq "pointer"}  @{$declaration   ->{'attributes'     }}
				    ) 
				{
				    foreach my $object ( @{$declaration->{'variables'}} ) {
					(my $name = $object) =~ s/^([a-zA-Z0-9_]+).*/$1/; # Strip away anything (e.g. assignment operators) after the variable name.
					next
					    if ( grep {lc($_) eq lc($name)} @ignore );
					$assignments .= "nullify(destination%".$name.")\n";
					$assignments .= "allocate(destination%".$name.",mold=self%".$name.")\n";
					$assignments .= "call self%".$name."%deepCopy(destination%".$name.")\n";
					$assignments .= "if (mpiSelf\%isMaster()) call Galacticus_Display_Message(var_str('functionClass[own] (class : ownerName : ownerLoc : objectLoc : sourceLoc): ".$name." : [destination] : ')//loc(destination)//' : '//loc(destination%".$name.")//' : '//".&Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($node,$lineNumber,compact => 1).",verbositySilent)\n"
					    if ( $debugging );
					$assignments .= "call destination%".$name."%autoHook()\n";
				    }
				};
				# Deep copy of objects with explicit deep copy actions.
				if
				    (
				     (
				      $declaration->{'intrinsic'} eq "class"
				      ||
				      $declaration->{'intrinsic'} eq "type"
				     )
				     &&
				     (grep {$_->{'type'} eq $type} @{$deepCopyActions->{'deepCopyActions'}})
				    ) {
					my $rank = 0;
					if ( grep {$_ =~ m/^dimension\s*\(/} @{$declaration->{'attributes'}} ) {
					    my $dimensionDeclarator = join(",",map {/^dimension\s*\(([a-zA-Z0-9_,]+)\)/} @{$declaration->{'attributes'}});
					    $rank        = ($dimensionDeclarator =~ tr/,//)+1;
					    $rankMaximum = $rank
						if ( $rank > $rankMaximum );
					}
					foreach my $variableName ( @{$declaration->{'variables'}} ) {
					    for(my $i=1;$i<=$rank;++$i) {
						$assignments .= (" " x $i)."do i".$i."=1,size(self%".$variableName.",dim=".$i.")\n";
					    }
					    my $arrayElement = $rank > 0 ? "(".join(",",map {"i".$_} 1..$rank).")" : "";
					    $assignments .= (" " x $rank)."call destination%".$variableName.$arrayElement."%deepCopyActions()\n";
					    for(my $i=1;$i<=$rank;++$i) {
						    $assignments .= (" " x ($rank+1-$i))."end do\n";
					    }					    
					}
				}
				# Deep copy of HDF5 objects.
				if
				    (
				     $declaration->{'intrinsic'} eq "type"
				     &&
				     $declaration->{'type'     } =~ m/^\s*hdf5object\s*$/i
				    ) {
					$deepCopyModules{'IO_HDF5'} = 1;
					$assignments .= "!\$ call hdf5Access%set  ()\n";
					$assignments .= "call self%".$_."%deepCopy(destination%".$_.")\n"
					    foreach ( @{$declaration->{'variables'}} );
					$assignments .= "!\$ call hdf5Access%unset()\n";
				}
				# Deep copy of non-(class,pointer) functionClass objects.
				if ( exists($class->{'deepCopy'}->{'functionClass'}) ) {
				    foreach my $object ( @{$declaration->{'variables'}} ) {
					(my $name = $object) =~ s/^([a-zA-Z0-9_]+).*/$1/; # Strip away anything (e.g. assignment operators) after the variable name.
					if ( grep {lc($_) eq lc($name)} split(/\s*,\s*/,$class->{'deepCopy'}->{'functionClass'}->{'variables'}) ) {
					    if ( grep {$_ eq "pointer"}  @{$declaration->{'attributes'}} ) {
						$assignments .= "nullify(destination%".$name.")\n";
						$assignments .= "if (associated(self%".$name.")) then\n";
						$assignments .= "allocate(destination%".$name.",mold=self%".$name.")\n";
					    }
					    $assignments .= "call self%".$name."%deepCopy(destination%".$name.")\n";
					    $assignments .= "if (mpiSelf\%isMaster()) call Galacticus_Display_Message(var_str('functionClass[own] (class : ownerName : ownerLoc : objectLoc : sourceLoc): ".$name." : [destination] : ')//loc(destination)//' : '//loc(destination%".$name.")//' : '//".&Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($node,$lineNumber,compact => 1).",verbositySilent)\n"
					    if ( $debugging );
					    $assignments .= "call destination%".$name."%autoHook()\n";
					    if ( grep {$_ eq "pointer"}  @{$declaration->{'attributes'}} ) {
						$assignments .= "end if\n";
					    }
					}
				    }
				}
				# Perform any increments.
				if ( exists($class->{'deepCopy'}->{'increment'}) ) {
				    my @increments = map {{variable => $_}} split(/\s*,\s*/,$class->{'deepCopy'}->{'increment'}->{'variables'});
				    foreach ( @increments ) {
					($_->{'host'} = $_->{'variable'}) =~ s/^([^%]+)%.+/$1/;
				    }
				    foreach my $object ( @{$declaration->{'variables'}} ) {
					(my $name = $object) =~ s/^([a-zA-Z0-9_]+).*/$1/; # Strip away anything (e.g. assignment operators) after the variable name.
					foreach my $increment ( @increments ) {
					    if ( lc($increment->{'host'}) eq lc($name) ) {
						$assignments .= "!\$omp atomic\n"
						    if ( exists($class->{'deepCopy'}->{'increment'}->{'atomic'}) && $class->{'deepCopy'}->{'increment'}->{'atomic'} eq "yes" );
						$assignments .= "destination\%".$increment->{'variable'}."=destination\%".$increment->{'variable'}."+1\n";
					    }
					}
				    }
				}
				# Deallocate FGSL interpolators.
				if
				    (
				     $declaration->{'intrinsic'} eq "type"
				     &&
				     $declaration->{'type'     } =~ m/^\s*fgsl_interp\s*$/i
				    ) {
					$assignments .= "destination%".$_."=fgsl_interp()\n"
					    foreach ( @{$declaration->{'variables'}} );
				}
				if
				    (
				     $declaration->{'intrinsic'} eq "type"
				     &&
				     $declaration->{'type'     } =~ m/^\s*fgsl_interp_accel\s*$/i
				    ) {
					$assignments .= "destination%".$_."=fgsl_interp_accel()\n"
					    foreach ( @{$declaration->{'variables'}} );
				}
				# Reinitialize OpenMP locks.
				if
				    (
				     $declaration->{'intrinsic'} eq "integer"
				     &&
				     exists ($declaration->{'type'})
				     &&
				     defined($declaration->{'type'})
				     &&
				     $declaration->{'type'     } =~ m/^\s*omp_lock_kind\s*$/i
				    ) {
					$assignments .= "!\$ call OMP_Init_Lock(destination\%".$_.")\n"
					    foreach ( @{$declaration->{'variables'}} );
				}
				# Reinitialize OpenMP read/write locks.
				if
				    (
				     $declaration->{'intrinsic'} eq "type"
				     &&
				     exists ($declaration->{'type'})
				     &&
				     defined($declaration->{'type'})
				     &&
				     $declaration->{'type'     } =~ m/^\s*ompReadWriteLock\s*$/i
				    ) {
					foreach ( @{$declaration->{'variables'}} ) {
					    my @dimensions =
						exists($declaration->{'attributes'}) 
						?
						map {/^dimension\s*\(([:,]+)\)/} @{$declaration->{'attributes'}} 
					        :
						undef();	    
					    if ( @dimensions ) {
						my @rank = split(",",$dimensions[0]);
						# Add loop index variables.
						$rankMaximum = scalar(@rank)
						    if ( scalar(@rank) > $rankMaximum );
						for(my $i=1;$i<=scalar(@rank);++$i) {
						    $assignments .= "!\$ do i".$i."=lbound(destination\%".$_.",dim=".$i."),ubound(destination\%".$_.",dim=".$i.")\n";
						}
						$assignments .= "!\$    call destination\%".$_."(".join(",",map {"i".$_} 1..scalar(@rank)).")%initialize()\n";
						for(my $i=1;$i<=scalar(@rank);++$i) {
						    $assignments .= "!\$ end do\n";
						}
					    } else {
						# Scalar lock.
						$assignments .= "!\$ call destination\%".$_."%initialize()\n";
					    }
					}
				}		    
			    }
			}
			$node = $node->{'sibling'};
		    }
		    # Move to the parent class.
		    $class = ($class->{'extends'} eq $directive->{'name'}) ? undef() : $classes{$class->{'extends'}};
		}
		# Add any objects declared in the base class.
		foreach my $data ( &List::ExtraUtils::as_array($directive->{'data'}) ) {
		    my $declarationSource;
		    if ( reftype($data) ) {
			$declarationSource = $data->{'content'}
			    if ( $data->{'scope'} eq "self" );
		    } else {
			$declarationSource = $data;
		    }
		    next
			unless ( defined($declarationSource) );
		    my $declaration = &Fortran::Utils::Unformat_Variables($declarationSource);
		    (my $type = $declaration->{'type'}) =~ s/(^\s*|\s*$)//g
			if ( $declaration->{'intrinsic'} eq "class" || $declaration->{'intrinsic'} eq "type" );
		    if
			(
			 $declaration->{'intrinsic'} eq "class"
			 &&
			 $declaration->{'type'     } =~ m/Class\s*$/
			 &&
			 grep {$_ eq "pointer"} @{$declaration->{'attributes'}}
			) {
			    foreach my $object ( @{$declaration->{'variables'}} ) {
				(my $name = $object) =~ s/^([a-zA-Z0-9_]+).*/$1/; # Strip away anything (e.g. assignment operators) after the variable name.
				$assignments .= "nullify(destination%".$name.")\n";
				$assignments .= "allocate(destination%".$name.",mold=self%".$name.")\n";
				$assignments .= "call self%".$name."%deepCopy(destination%".$name.")\n";
				$assignments .= "if (mpiSelf\%isMaster()) call Galacticus_Display_Message(var_str('functionClass[own] (class : ownerName : ownerLoc : objectLoc : sourceLoc): ".$name." : [destination] : ')//loc(destination)//' : '//loc(destination%".$name.")//' : '//".&Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($node,$lineNumber,compact => 1).",verbositySilent)\n"
				    if ( $debugging );
			       	$assignments .= "call destination%".$name."%autoHook()\n";
			    }
		    }
		    # Deep copy of objects with explicit deep copy actions.
		    if
			(
			 (
			  $declaration->{'intrinsic'} eq "class"
			  ||
			  $declaration->{'intrinsic'} eq "type"
			 )
			 &&
			 (grep {$_->{'type'} eq $type} @{$deepCopyActions->{'deepCopyActions'}})
			) {
			    my $rank = 0;
			    if ( grep {$_ =~ m/^dimension\s*\(/} @{$declaration->{'attributes'}} ) {
				my $dimensionDeclarator = join(",",map {/^dimension\s*\(([a-zA-Z0-9_,]+)\)/} @{$declaration->{'attributes'}});
				$rank        = ($dimensionDeclarator =~ tr/,//)+1;
				$rankMaximum = $rank
				    if ( $rank > $rankMaximum );
			    }
			    foreach my $variableName ( @{$declaration->{'variables'}} ) {
				for(my $i=1;$i<=$rank;++$i) {
				    $assignments .= (" " x $i)."do i".$i."=1,size(self%".$variableName.",dim=".$i.")\n";
				}
				my $arrayElement = $rank > 0 ? "(".join(",",map {"i".$_} 1..$rank).")" : "";
				$assignments .= (" " x $rank)."call destination%".$variableName.$arrayElement."%deepCopyActions()\n";
				for(my $i=1;$i<=$rank;++$i) {
					$assignments .= (" " x ($rank+1-$i))."end do\n";
				}					    
			    }
		    }
		    # Deep copy of HDF5 objects.
		    if
			(
			 $declaration->{'intrinsic'} eq "type"
			 &&
			 $declaration->{'type'     } =~ m/^\s*hdf5object\s*$/i
			) {
			    $deepCopyModules{'IO_HDF5'} = 1;
			    $assignments .= "!\$ call hdf5Access%set  ()\n";
			    $assignments .= "call self%".$_."%deepCopy(destination%".$_.")\n"
				foreach ( @{$declaration->{'variables'}} );
			    $assignments .= "!\$ call hdf5Access%unset()\n";
		    }
		    # Deep copy of non-(class,pointer) functionClass objects.
		    if ( exists($class->{'deepCopy'}->{'functionClass'}) ) {
			foreach my $name ( split(/\s*,\s*/,$class->{'deepCopy'}->{'functionClass'}->{'variables'}) ) {
			    if ( grep {$_ eq "pointer"}  @{$declaration->{'attributes'}} ) {
				$assignments .= "nullify(destination%".$name.")\n";
				$assignments .= "if (associated(self%".$name.")) then\n";
				$assignments .= "allocate(destination%".$name.",mold=self%".$name.")\n";
			    }
			    $assignments .= "call self%".$name."%deepCopy(destination%".$name.")\n";
			    $assignments .= "if (mpiSelf\%isMaster()) call Galacticus_Display_Message(var_str('functionClass[own] (class : ownerName : ownerLoc : objectLoc : sourceLoc): ".$name." : [destination] : ')//loc(destination)//' : '//loc(destination%".$name.")//' : '//".&Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($node,$lineNumber,compact => 1).",verbositySilent)\n"
				if ( $debugging );
			    $assignments .= "call destination%".$name."%autoHook()\n";
			    if ( grep {$_ eq "pointer"}  @{$declaration->{'attributes'}} ) {
				$assignments .= "end if\n";
			    }
			}
		    }
		    # Perform any increments.
		    if ( exists($class->{'deepCopy'}->{'increment'}) ) {
			my @increments = map {{variable => $_}} split(/\s*,\s*/,$class->{'deepCopy'}->{'increment'}->{'variables'});
			foreach ( @increments ) {
			    ($_->{'host'} = $_->{'variable'}) =~ s/^([^%]+)%.+/$1/;
			}
			foreach my $object ( @{$declaration->{'variables'}} ) {
			    (my $name = $object) =~ s/^([a-zA-Z0-9_]+).*/$1/; # Strip away anything (e.g. assignment operators) after the variable name.
			    foreach my $increment ( @increments ) {
				if ( lc($increment->{'host'}) eq lc($name) ) {
				    $assignments .= "!\$omp atomic\n"
					if ( exists($class->{'deepCopy'}->{'increment'}->{'atomic'}) && $class->{'deepCopy'}->{'increment'}->{'atomic'} eq "yes" );
				    $assignments .= "destination\%".$increment->{'variable'}."=destination\%".$increment->{'variable'}."+1\n";
				}
			    }
			}
		    }
		    # Deallocate FGSL interpolators.
		    if
			(
			 $declaration->{'intrinsic'} eq "type"
			 &&
			 $declaration->{'type'     } =~ m/^\s*fgsl_interp\s*$/i
			) {
			    $assignments .= "destination%".$_."=fgsl_interp()\n"
				foreach ( @{$declaration->{'variables'}} );
		    }
		    if
			(
			 $declaration->{'intrinsic'} eq "type"
			 &&
			 $declaration->{'type'     } =~ m/^\s*fgsl_interp_accel\s*$/i
			) {
			    $assignments .= "destination%".$_."=fgsl_interp_accel()\n"
				foreach ( @{$declaration->{'variables'}} );
		    }
		}
		# Add any objects declared in the functionClassType class.
		if ( defined($functionClassType) ) {
		    # Search the node for declarations.
		    my $node = $functionClassType->{'node'}->{'firstChild'};
		    while ( $node ) {
			if ( $node->{'type'} eq "declaration" ) {
			    foreach my $declaration ( @{$node->{'declarations'}} ) {
				# Deep copy of functionClass objects.
				(my $type = $declaration->{'type'}) =~ s/(^\s*|\s*$)//g
				    if ( $declaration->{'intrinsic'} eq "class" );
				if
				    (
				     $declaration->{'intrinsic'} eq "class"
				     &&
				     (grep {$_ eq $type    } (@{$stateStorables->{'functionClasses'}},@{$stateStorables->{'functionClassInstances'}}))
				     &&
				     grep {$_ eq "pointer"}  @{$declaration   ->{'attributes'     }}
				    ) 
				{
				    foreach my $object ( @{$declaration->{'variables'}} ) {
					(my $name = $object) =~ s/^([a-zA-Z0-9_]+).*/$1/; # Strip away anything (e.g. assignment operators) after the variable name.
					$assignments .= "nullify(destination%".$name.")\n";
					$assignments .= "allocate(destination%".$name.",mold=self%".$name.")\n";
					$assignments .= "call self%".$name."%deepCopy(destination%".$name.")\n";
					$assignments .= "if (mpiSelf\%isMaster()) call Galacticus_Display_Message(var_str('functionClass[own] (class : ownerName : ownerLoc : objectLoc : sourceLoc): ".$name." : [destination] : ')//loc(destination)//' : '//loc(destination%".$name.")//' : '//".&Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($node,$lineNumber,compact => 1).",verbositySilent)\n"
					    if ( $debugging );
					$assignments .= "call destination%".$name."%autoHook()\n";
				    }
				};
				# Deep copy of HDF5 objects.
				if
				    (
				     $declaration->{'intrinsic'} eq "type"
				     &&
				     $declaration->{'type'     } =~ m/^\s*hdf5object\s*$/i
				    ) {
					$deepCopyModules{'IO_HDF5'} = 1;
					$assignments .= "!\$ call hdf5Access%set  ()\n";
					$assignments .= "call self%".$_."%deepCopy(destination%".$_.")\n"
					    foreach ( @{$declaration->{'variables'}} );
					$assignments .= "!\$ call hdf5Access%unset()\n";
				}
				# Deep copy of non-(class,pointer) functionClass objects.
				if ( exists($class->{'deepCopy'}->{'functionClass'}) ) {
				    foreach my $object ( @{$declaration->{'variables'}} ) {
					(my $name = $object) =~ s/^([a-zA-Z0-9_]+).*/$1/; # Strip away anything (e.g. assignment operators) after the variable name.
					if ( grep {lc($_) eq lc($name)} split(/\s*,\s*/,$class->{'deepCopy'}->{'functionClass'}->{'variables'}) ) {
					    if ( grep {$_ eq "pointer"}  @{$declaration->{'attributes'}} ) {
						$assignments .= "nullify(destination%".$name.")\n";
						$assignments .= "if (associated(self%".$name.")) then\n";
						$assignments .= "allocate(destination%".$name.",mold=self%".$name.")\n";
					    }
					    $assignments .= "call self%".$name."%deepCopy(destination%".$name.")\n";
					    $assignments .= "if (mpiSelf\%isMaster()) call Galacticus_Display_Message(var_str('functionClass[own] (class : ownerName : ownerLoc : objectLoc : sourceLoc): ".$name." : [destination] : ')//loc(destination)//' : '//loc(destination%".$name.")//' : '//".&Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($node,$lineNumber,compact => 1).",verbositySilent)\n"
					    if ( $debugging );
				   	    $assignments .= "call destination%".$name."%autoHook()\n";
					    if ( grep {$_ eq "pointer"}  @{$declaration->{'attributes'}} ) {
						$assignments .= "end if\n";
					    }
					}
				    }
				}
				# Deallocate FGSL interpolators.
				if
				    (
				     $declaration->{'intrinsic'} eq "type"
				     &&
				     $declaration->{'type'     } =~ m/^\s*fgsl_interp\s*$/i
				    ) {
					$assignments .= "destination%".$_."=fgsl_interp()\n"
					    foreach ( @{$declaration->{'variables'}} );
				}
				if
				    (
				     $declaration->{'intrinsic'} eq "type"
				     &&
				     $declaration->{'type'     } =~ m/^\s*fgsl_interp_accel\s*$/i
				    ) {
					$assignments .= "destination%".$_."=fgsl_interp_accel()\n"
					    foreach ( @{$declaration->{'variables'}} );
				}
				# Reinitialize OpenMP locks.
				if
				    (
				     $declaration->{'intrinsic'} eq "integer"
				     &&
				     exists ($declaration->{'type'})
				     &&
				     defined($declaration->{'type'})
				     &&
				     $declaration->{'type'     } =~ m/^\s*omp_lock_kind\s*$/i
				    ) {
					$assignments .= "!\$ call OMP_Init_Lock(destination\%".$_.")\n"
					    foreach ( @{$declaration->{'variables'}} );
				}
				# Reinitialize OpenMP read/write locks.
				if
				    (
				     $declaration->{'intrinsic'} eq "type"
				     &&
				     exists ($declaration->{'type'})
				     &&
				     defined($declaration->{'type'})
				     &&
				     $declaration->{'type'     } =~ m/^\s*ompReadWriteLock\s*$/i
				    ) {
					foreach ( @{$declaration->{'variables'}} ) {
					    my @dimensions =
						exists($declaration->{'attributes'}) 
						?
						map {/^dimension\s*\(([:,]+)\)/} @{$declaration->{'attributes'}} 
					        :
						undef();	    
					    if ( @dimensions ) {
						my @rank = split(",",$dimensions[0]);
						# Add loop index variables.
						$rankMaximum = scalar(@rank)
						    if ( scalar(@rank) > $rankMaximum );
						for(my $i=1;$i<=scalar(@rank);++$i) {
						    $assignments .= "!\$ do i".$i."=lbound(destination\%".$_.",dim=".$i."),ubound(destination\%".$_.",dim=".$i.")\n";
						}
						$assignments .= "!\$    call destination\%".$_."(".join(",",map {"i".$_} 1..scalar(@rank)).")%initialize()\n";
						for(my $i=1;$i<=scalar(@rank);++$i) {
						    $assignments .= "!\$ end do\n";
						}
					    } else {
						# Scalar lock.
						$assignments .= "!\$ call destination\%".$_."%initialize()\n";
					    }
					}
				}		    
			    }
			}
			$node = $node->{'sibling'};
		    }
		}
		# Check that the type of the destination matches, and perform the copy. Reset the reference count to the copy.
		$deepCopyCode .= "type is (".$nonAbstractClass->{'name'}.")\n";
		$deepCopyCode .= "select type (destination)\n";
		$deepCopyCode .= "type is (".$nonAbstractClass->{'name'}.")\n";
		$deepCopyCode .= "destination=self\n";
		$deepCopyCode .= $assignments
		    if ( defined($assignments) );
		$deepCopyCode .= "class default\n";
		$deepCopyCode .= "call Galacticus_Error_Report('destination and source types do not match'//".&Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($nonAbstractClass->{'node'},$nonAbstractClass->{'node'}->{'line'}).")\n";
		$deepCopyCode .= "end select\n";
		# Specify required modules.
		$deepCopyModules{'Galacticus_Error'} = 1;
	    }
            $deepCopyCode .= "end select\n";
            # Reset the reference count to this newly created object.
            $deepCopyCode .= "call destination%referenceCountReset()\n";
            # Reset the state operation ID if necessary.
            $deepCopyCode .= "destination%stateOperationID=0_c_size_t\n";
            # Insert any iterator variables needed.
            $deepCopyCode = "integer :: ".join(",",map {"i".$_} 1..$rankMaximum)."\n".$deepCopyCode
                if ( $rankMaximum > 0 );
	    $methods{'deepCopy'} = 
	    {
		description => "Perform a deep copy of the object.",
		type        => "void",
		recursive   => "yes",
		pass        => "yes",
		modules     => join(" ",keys(%deepCopyModules)),
		argument    => [ "class(".$directive->{'name'}."Class), intent(inout) :: destination" ],
		code        => $deepCopyCode
	    };
	    # Add "stateStore" and "stateRestore" method.
	    my $stateStoreCode;
	    my $stateRestoreCode;
	    my %stateStoreModules   = ( "Galacticus_Display" => 1, "ISO_Varying_String" => 1, "String_Handling" => 1, "FGSL,only:fgsl_file" => 1, "ISO_C_Binding" => 1 );
	    my %stateRestoreModules = ( "Galacticus_Display" => 1, "ISO_Varying_String" => 1, "String_Handling" => 1, "FGSL,only:fgsl_file" => 1, "ISO_C_Binding" => 1 );
	    my @outputUnusedVariables;
	    my @inputUnusedVariables;
	    my $allocatablesFound = 0;
	    my $dimensionalsFound = 0;
	    my $fgslStateFileUsed = 0;
	    my $stateFileUsed     = 0;
	    my $labelUsed         = 0;
	    $rankMaximum          = 0;
	    $stateStoreCode   .= "call Galacticus_Display_Indent(var_str('storing state for \""  .$directive->{'name'}."\" [position: ')//FTell(stateFile)//']',verbosity=verbosityWorking)\n";
	    $stateRestoreCode .= "call Galacticus_Display_Indent(var_str('restoring state for \"".$directive->{'name'}."\" [position: ')//FTell(stateFile)//']',verbosity=verbosityWorking)\n";
	    $stateStoreCode   .= "select type (self)\n";
	    $stateRestoreCode .= "select type (self)\n";
	    foreach my $nonAbstractClass ( @nonAbstractClasses ) {
		# Build the code.
		$stateStoreCode   .= "type is (".$nonAbstractClass->{'name'}.")\n";
		$stateRestoreCode .= "type is (".$nonAbstractClass->{'name'}.")\n";
		$stateStoreCode   .= "if (self%stateOperationID == stateOperationID) then\n"; # If this object was already stored, don't do it again.
		$stateStoreCode   .= " call Galacticus_Display_Unindent('skipping - already stored',verbosity=verbosityWorking)\n";
		$stateStoreCode   .= " return\n";
		$stateStoreCode   .= "end if\n";
		$stateStoreCode   .= "self%stateOperationID=stateOperationID\n";
		$stateRestoreCode .= "if (self%stateOperationID == stateOperationID) then\n"; # If this object was already restored, don't do it again.
		$stateRestoreCode .= " call Galacticus_Display_Unindent('skipping - already restored',verbosity=verbosityWorking)\n";
		$stateRestoreCode .= " return\n";
		$stateRestoreCode .= "end if\n";
		$stateRestoreCode .= "self%stateOperationID=stateOperationID\n";
		$stateStoreCode   .= " call Galacticus_Display_Message('object type \"".$nonAbstractClass->{'name'}."\"',verbosity=verbosityWorking)\n";
		$stateRestoreCode .= " call Galacticus_Display_Message('object type \"".$nonAbstractClass->{'name'}."\"',verbosity=verbosityWorking)\n";
		(my $label = $nonAbstractClass->{'name'}) =~ s/^$directive->{'name'}//;
		$label = lcfirst($label)
		    unless ( $label =~ m/^[A-Z]{2,}/ );
		my $hasCustomStateStore   = 0;
		my $hasCustomStateRestore = 0;
		my $extensionOf;
		# Generate code to output all variables from this class (and any parent class).
		my $outputCode;
		my $inputCode;
		my @staticVariables;
		my $class = $nonAbstractClass;
		while ( $class ) {
		    my $node = $class->{'tree'}->{'firstChild'};
		    $node = $node->{'sibling'}
		    while ( $node && ( $node->{'type'} ne "type" || ( ! exists($node->{'name'}) || $node->{'name'} ne $class->{'name'} ) ) );
		    last
			unless ( $node );
		    # Find the parent class.
		    if ( $class == $nonAbstractClass && $node->{'opener'} =~ m/,\s*extends\s*\(\s*([a-zA-Z0-9_]+)\s*\)/ ) {
			$extensionOf = $1;
		    }
		    # Find any variables to be excluded from state store/restore.
		    my @excludes = exists($class->{'stateStorable'}->{'exclude'}->{'variables'}) ? split(/\s*,\s*/,$class->{'stateStorable'}->{'exclude'}->{'variables'}) : ();
		    # Search the node for declarations.
		    $node = $node->{'firstChild'};
		    while ( $node ) {
			if ( $node->{'type'} eq "declaration" ) {			    
			    foreach my $declaration ( @{$node->{'declarations'}} ) {
				# Identify variable type.
				if ( $declaration->{'intrinsic'} eq "procedure" ) {
				    # Type-bound procedure - nothing to do.
				} elsif ( $declaration->{'intrinsic'} eq "class" || $declaration->{'intrinsic'} eq "type" ) {
				    # Look for pointers to functionClasses.
				    (my $type = $declaration->{'type'}) =~ s/\s//g;
				    if ( 
					$declaration->{'intrinsic'} eq "class"
					&&
					(grep {$_ eq "pointer"} @{$declaration   ->{'attributes'     }})
					&&
					(grep {$_ eq $type    } @{$stateStorables->{'functionClasses'}})
					) {
					# Pointer to a functionClass object.
					foreach ( @{$declaration->{'variables'}} ) {
					    $labelUsed = 1;
					    (my $variableName = $_) =~ s/\s*=.*$//;
					    next
						if ( grep {lc($_) eq lc($variableName)} @excludes );
					    $outputCode .= " if (Galacticus_Verbosity_Level() >= verbosityWorking) then\n";
					    $outputCode .= "  select type (c__ => self%".$variableName.")\n";
					    $outputCode .= "  class is (".$declaration->{'type'}.")\n";
					    $outputCode .= "   write (label,'(i16)') sizeof(c__)\n";
					    $outputCode .= "  end select\n";
					    $outputCode .= "  call Galacticus_Display_Message('storing \"".$variableName."\" with size '//trim(adjustl(label))//' bytes')\n";
					    $outputCode .= " end if\n";
					    $inputCode  .= " call Galacticus_Display_Message('restoring \"".$variableName."\"',verbosity=verbosityWorking)\n";
					    $outputCode .= " call self%".$variableName."%stateStore  (stateFile,fgslStateFile,stateOperationID)\n";
					    $inputCode  .= " call self%".$variableName."%stateRestore(stateFile,fgslStateFile,stateOperationID)\n";
					    $stateFileUsed     = 1;
					    $fgslStateFileUsed = 1;
					}
				    } elsif (
					(
					 (  grep {$_->{'type'} eq $type    } @{$stateStorables->{'stateStorables'        }})
					 ||
					 (  grep {$_           eq $type    } @{$stateStorables->{'functionClassInstances'}})
					)
					&&
					(! grep {$_           eq "pointer"} @{$declaration   ->{'attributes'            }})				       
					){
					# This is a non-pointer object which is explicitly stateStorable.
					# Validate: Currently we do not support store/restore of polymorphic functionClass objects.
					die("Galacticus::Build::SourceTree::Process::FunctionClass::Process_FunctionClass(): storing of polymorphic functionClass objects is not implemented")
					    if
					    (
					     $declaration->{'intrinsic'} eq "class"
					     &&
					     (grep {$_ eq $type} @{$stateStorables->{'functionClassInstances'}})
					    );
					my $isFunctionClass = grep {$_ eq $type} @{$stateStorables->{'functionClassInstances'}};
					# Construct code to output.
					foreach ( @{$declaration->{'variables'}} ) {
					    (my $variableName = $_) =~ s/\s*=.*$//;
					    next
						if ( grep {lc($_) eq lc($variableName)} @excludes );
					    my $rank = 0;
					    if ( grep {$_ =~ m/^dimension\s*\(/} @{$declaration->{'attributes'}} ) {
						my $dimensionDeclarator = join(",",map {/^dimension\s*\(([a-zA-Z0-9_,]+)\)/} @{$declaration->{'attributes'}});
						$rank        = ($dimensionDeclarator =~ tr/,//)+1;
						$rankMaximum = $rank
						    if ( $rank > $rankMaximum );
					    }
					    if ( grep {$_ eq "allocatable"} @{$declaration->{'attributes'}} ) {
						# For allocatable variables we must first store the shape so that they can be reallocated on restore.
						$allocatablesFound  = 1;
						$dimensionalsFound  = 1
						    if ( $rank > 0 );
						$outputCode .= " if (allocated(self%".$variableName.")) then\n";
						$outputCode .= "  write (stateFile) .true.\n";
						$outputCode .= "  write (stateFile) shape(self%".$variableName.",kind=c_size_t)\n"
						    if ( $rank > 0 );
						$inputCode  .= " read (stateFile) wasAllocated\n";
						$inputCode  .= " if (allocated(self%".$variableName.")) deallocate(self%".$variableName.")\n";
						$inputCode  .= " if (wasAllocated) then\n";
						if ( $rank > 0 ) {
						    $inputCode  .= "  allocate(storedShape(".$rank."))\n";
						    $inputCode  .= "  read (stateFile) storedShape\n";
						    $inputCode  .= "  deallocate(storedShape)\n";
						}
						if ( $declaration->{'intrinsic'} eq "class" ) {
						    $inputCode  .= "  call ".$type."ClassRestore(self%".$variableName.",stateFile)\n";
						} else {
						    $inputCode  .= "  allocate(self%".$variableName.($rank > 0 ? "(".join(",",map {"storedShape(".$_.")"} 1..$rank).")" : "").")\n";
						}
					    }
					    for(my $i=1;$i<=$rank;++$i) {
						$outputCode .= (" " x $i)."do i".$i."=1,size(self%".$variableName.",dim=".$i.")\n";
						$inputCode  .= (" " x $i)."do i".$i."=1,size(self%".$variableName.",dim=".$i.")\n";
					    }
					    my $arrayElement = $rank > 0 ? "(".join(",",map {"i".$_} 1..$rank).")" : "";
					    $labelUsed   = 1;
					    $outputCode .= " if (Galacticus_Verbosity_Level() >= verbosityWorking) then\n";
					    if ( $declaration->{'intrinsic'} eq "class" ) {
						$outputCode .= "  select type (c__ => self%".$variableName.")\n";
						$outputCode .= "  class is (".$declaration->{'type'}.")\n";
						$outputCode .= "   write (label,'(i16)') sizeof(c__".$arrayElement.")\n";
						$outputCode .= "  end select\n";
					    } else {
						$outputCode .= "   write (label,'(i16)') sizeof(self%".$variableName.")\n";
					    }
					    $outputCode .= "  call Galacticus_Display_Message('storing \"".$variableName.$arrayElement."\" with size '//trim(adjustl(label))//' bytes')\n";
					    $outputCode .= " end if\n";
					    $inputCode  .= " call Galacticus_Display_Message('restoring \"".$variableName.$arrayElement."\"',verbosity=verbosityWorking)\n";
					    $inputCode  .= (" " x $rank)." call self%".$variableName.$arrayElement."%stateRestore(stateFile,fgslStateFile".($isFunctionClass ? ",stateOperationID" : "").")\n";
					    $outputCode .= (" " x $rank)." call self%".$variableName.$arrayElement."%stateStore  (stateFile,fgslStateFile".($isFunctionClass ? ",stateOperationID" : ",storeIdentifier=".($declaration->{'intrinsic'} eq "class" ? ".true." : ".false.")).")\n";
					    for(my $i=1;$i<=$rank;++$i) {
						$outputCode .= (" " x ($rank+1-$i))."end do\n";
						$inputCode  .= (" " x ($rank+1-$i))."end do\n";
					    }
					    if ( grep {$_ eq "allocatable"} @{$declaration->{'attributes'}} ) {
						$inputCode  .= " end if\n";
						$outputCode .= " else\n";
						$outputCode .= "  write (stateFile) .false.\n";
						$outputCode .= " end if\n";
					    }
					    $stateFileUsed      = 1;
					    $fgslStateFileUsed  = 1;
					}
				    }
				} else {
				    # Intrinsic type.
				    if ( grep {$_ eq "pointer"} @{$declaration->{'attributes'}} ) {
					# Pointers are currently not handled.
				    } elsif ( exists($declaration->{'type'}) && defined($declaration->{'type'}) && $declaration->{'type'} =~ m/^\s*omp_lock_kind\s*/ ) {
					# Do not store OpenMP lock variables.
				    } elsif ( grep {$_ eq "allocatable"} @{$declaration->{'attributes'}} ) {
					# For allocatable variables we must first store the shape so that they can be reallocated on restore.
					my $dimensionDeclarator = join(",",map {/^dimension\s*\(([:,]+)\)/} @{$declaration->{'attributes'}});
					my $rank = ($dimensionDeclarator =~ tr/://);
					foreach my $variableName ( @{$declaration->{'variables'}} ) {
					    next
						if ( grep {lc($_) eq lc($variableName)} @excludes );
					    $allocatablesFound  = 1;
					    $dimensionalsFound  = 1;
					    $stateFileUsed      = 1;
					    $labelUsed          = 1;
					    $outputCode        .= " if (allocated(self%".$variableName.")) then\n";
					    $outputCode        .= "  if (Galacticus_Verbosity_Level() >= verbosityWorking) then\n";
					    $outputCode        .= "   write (label,'(i16)') sizeof(self%".$variableName.")\n";
					    $outputCode        .= "   call Galacticus_Display_Message('storing \"".$variableName."\" with size '//trim(adjustl(label))//' bytes')\n";
					    $outputCode        .= "  end if\n";
					    $outputCode        .= "  write (stateFile) .true.\n";
					    $outputCode        .= "  write (stateFile) shape(self%".$variableName.",kind=c_size_t)\n";
					    $outputCode        .= "  write (stateFile) self%".$variableName."\n";
					    $outputCode        .= " else\n";
					    $outputCode        .= "  write (stateFile) .false.\n";
					    $outputCode        .= " end if\n";
					    $inputCode         .= " read (stateFile) wasAllocated\n";
					    $inputCode         .= " if (allocated(self%".$variableName.")) deallocate(self%".$variableName.")\n";
					    $inputCode         .= " if (wasAllocated) then\n";
					    $inputCode         .= "  call Galacticus_Display_Message('restoring \"".$variableName."\"',verbosity=verbosityWorking)\n";
					    $inputCode         .= "  allocate(storedShape(".$rank."))\n";
					    $inputCode         .= "  read (stateFile) storedShape\n";
					    $inputCode         .= "  allocate(self%".$variableName."(".join(",",map {"storedShape(".$_.")"} 1..$rank)."))\n";
					    $inputCode         .= "  deallocate(storedShape)\n";
					    $inputCode         .= "  read (stateFile) self%".$variableName."\n";
					    $inputCode         .= " end if\n";
					}
				    } else {
					# Statically-sized variable.
					foreach ( @{$declaration->{'variables'}} ) {
					    (my $variableName = $_) =~ s/\s*=.*$//;
					    next
						if ( grep {lc($_) eq lc($variableName)} @excludes );
					    my $store = 1;
					    if ( exists($class->{'stateStorable'}) && exists($class->{'stateStorable'}->{'restoreTo'}) ) {
						foreach ( &List::ExtraUtils::as_array($class->{'stateStorable'}->{'restoreTo'}) ) {
						    my @variables = split(/\s*,\s*/,$_->{'variables'});
						    if ( grep {lc($_) eq lc($variableName)} @variables ) {
							$store = 0;
							$inputCode .= " self%".$variableName."=".$_->{'state'}."\n";
						    }
						}
					    }
					    push(@staticVariables,$variableName)
						if ( $store );
					}
				    }
				}
				$hasCustomStateStore   = 1
				    if
				    (
				     $declaration->{'intrinsic'} eq "procedure"
				     &&
				     $declaration->{'variables'}->[0] =~ m/^stateStore=>/
				    );
				$hasCustomStateRestore = 1
				    if
				    (
				     $declaration->{'intrinsic'} eq "procedure"
				     &&
				     $declaration->{'variables'}->[0] =~ m/^stateRestore=>/
				    );
			    }
			}
			$node = $node->{'type'} eq "contains" ? $node->{'firstChild'} : $node->{'sibling'};
		    }
		    # Move to the parent class.
		    $class = ($class->{'extends'} eq $directive->{'name'}) ? undef() : $classes{$class->{'extends'}};
		}
		# Find any variables to be excluded from state store/restore.
		my @excludes = exists($directive->{'stateStorable'}->{'exclude'}->{'variables'}) ? split(/\s*,\s*/,$directive->{'stateStorable'}->{'exclude'}->{'variables'}) : ();
		# Add any variables declared in the base class.
		foreach my $data ( &List::ExtraUtils::as_array($directive->{'data'}) ) {
		    my $declarationSource;
		    if ( reftype($data) ) {
			$declarationSource = $data->{'content'}
			if ( $data->{'scope'} eq "self" );
		    } else {
			$declarationSource = $data;
		    }
		    next
			unless ( defined($declarationSource) );
		    my $declaration = &Fortran::Utils::Unformat_Variables($declarationSource);
		    die("Galacticus::Build::SourceTree::Process::FunctionClass::Process_FunctionClass(): unable to parse variable declaration")
			unless ( defined($declaration) );
		    # Identify variable type.
		    if ( $declaration->{'intrinsic'} eq "procedure" ) {
			# Type-bound procedure - nothing to do.
		    } elsif ( $declaration->{'intrinsic'} eq "class" || $declaration->{'intrinsic'} eq "type" ) {
			# Look for pointers to functionClasses.
			(my $type = $declaration->{'type'}) =~ s/\s//g;
			if ( 
			    $declaration->{'intrinsic'} eq "class"
			    &&
			    (grep {$_ eq "pointer"} @{$declaration   ->{'attributes'     }})
			    &&
			    (grep {$_ eq $type    } @{$stateStorables->{'functionClasses'}})
			    ) {
			    # Pointer to a functionClass object.
			    foreach ( @{$declaration->{'variables'}} ) {
				(my $variableName = $_) =~ s/\s*=.*$//;
				$labelUsed  = 1;
				$outputCode .= " if (Galacticus_Verbosity_Level() >= verbosityWorking) then\n";
				$outputCode .= "  select type (c__ => self%".$variableName.")\n";
				$outputCode .= "  class is (".$declaration->{'type'}.")\n";
				$outputCode .= "   write (label,'(i16)') sizeof(c__)\n";
				$outputCode .= "  end select\n";
				$outputCode .= "  call Galacticus_Display_Message('storing \"".$variableName."\" with size '//trim(adjustl(label))//' bytes')\n";
				$outputCode .= " end if\n";
				$inputCode  .= " call Galacticus_Display_Message('restoring \"".$variableName."\"',verbosity=verbosityWorking)\n";
				$outputCode .= " call self%".$variableName."%stateStore  (stateFile,fgslStateFile,stateOperationID)\n";
				$inputCode  .= " call self%".$variableName."%stateRestore(stateFile,fgslStateFile,stateOperationID)\n";
				$stateFileUsed     = 1;
				$fgslStateFileUsed = 1;
			    }
			} elsif (
			    (
			     (  grep {$_->{'type'} eq $type    } @{$stateStorables->{'stateStorables'        }})
			     ||
			     (  grep {$_           eq $type    } @{$stateStorables->{'functionClassInstances'}})
			    )
			    &&
			    (! grep {$_           eq "pointer"} @{$declaration   ->{'attributes'            }})
			    ){
			    # This is a non-pointer object which is explicitly stateStorable or implicitly storeable by virtue of being a functionClass.
			    # Validate: Currently we do not support store/restore of polymorphic functionClass objects.
			    die("Galacticus::Build::SourceTree::Process::FunctionClass::Process_FunctionClass(): storing of polymorphic functionClass objects is not implemented")
				if
				(
				 $declaration->{'intrinsic'} eq "class"
				 &&
				 (grep {$_ eq $type} @{$stateStorables->{'functionClassInstances'}})
				);
			    # Construct code to output.
			    foreach ( @{$declaration->{'variables'}} ) {
				(my $variableName = $_) =~ s/\s*=.*$//;
				next
				    if ( grep {lc($_) eq lc($variableName)} @excludes );
				my $rank = 0;
				if ( grep {$_ =~ m/^dimension\s*\(/} @{$declaration->{'attributes'}} ) {
				    my $dimensionDeclarator = join(",",map {/^dimension\s*\(([:,]+)\)/} @{$declaration->{'attributes'}});
				    $rank        = ($dimensionDeclarator =~ tr/,//)+1;
				    $rankMaximum = $rank
					if ( $rank > $rankMaximum );
				}
				if ( grep {$_ eq "allocatable"} @{$declaration->{'attributes'}} ) {
				    # For allocatable variables we must first store the shape so that they can be reallocated on restore.			
				    $allocatablesFound  = 1;
				    $dimensionalsFound  = 1
					if ( $rank > 0 );
				    $outputCode .= " if (allocated(self%".$variableName.")) then\n";
				    $outputCode .= "  write (stateFile) .true.\n";
				    $outputCode .= "  write (stateFile) shape(self%".$variableName.",kind=c_size_t)\n"
					if ( $rank > 0 );
				    $inputCode  .= " read (stateFile) wasAllocated\n";
				    $inputCode  .= " if (allocated(self%".$variableName.")) deallocate(self%".$variableName.")\n";
				    $inputCode  .= " if (wasAllocated) then\n";
				    if ( $rank > 0 ) {
					$inputCode  .= "  allocate(storedShape(".$rank."))\n";
					$inputCode  .= "  read (stateFile) storedShape\n";
					$inputCode  .= "  deallocate(storedShape)\n";
				    }
				    if ( $declaration->{'intrinsic'} eq "class" ) {
					$inputCode  .= "  call ".$type."ClassRestore(self%".$variableName.",stateFile)\n";
				    } else {
					$inputCode  .= "  allocate(self%".$variableName.($rank > 0 ? "(".join(",",map {"storedShape(".$_.")"} 1..$rank).")" : "").")\n";
				    }
				}
				for(my $i=1;$i<=$rank;++$i) {
				    $outputCode .= (" " x $i)."do i".$i."=1,size(self%".$variableName.",dim=".$i.")\n";
				    $inputCode  .= (" " x $i)."do i".$i."=1,size(self%".$variableName.",dim=".$i.")\n";
				}
				my $arrayElement = $rank > 0 ? "(".join(",",map {"i".$_} 1..$rank).")" : "";
				$labelUsed   = 1;
				$outputCode .= " if (Galacticus_Verbosity_Level() >= verbosityWorking) then\n";
				if ( $declaration->{'intrinsic'} eq "class" ) {
				    $outputCode .= "  select type (c__ => self%".$variableName.")\n";
				    $outputCode .= "  class is (".$declaration->{'type'}.")\n";
				    $outputCode .= "   write (label,'(i16)') sizeof(c__".$arrayElement.")\n";
				    $outputCode .= "  end select\n";
				} else {
				    $outputCode .= "   write (label,'(i16)') sizeof(self%".$variableName.$arrayElement.")\n";
				}
				$outputCode .= "  call Galacticus_Display_Message('storing \"".$variableName.$arrayElement."\" with size '//trim(adjustl(label))//' bytes')\n";
				$outputCode .= " end if\n";
				$inputCode  .= " call Galacticus_Display_Message('restoring \"".$variableName.$arrayElement."\"',verbosity=verbosityWorking)\n";
				$inputCode  .= (" " x $rank)." call self%".$variableName.$arrayElement."%stateRestore(stateFile,fgslStateFile)\n";
				$outputCode .= (" " x $rank)." call self%".$variableName.$arrayElement."%stateStore  (stateFile,fgslStateFile,storeIdentifier=".($declaration->{'intrinsic'} eq "class" ? ".true." : ".false.").")\n";
				for(my $i=1;$i<=$rank;++$i) {
				    $outputCode .= (" " x ($rank+1-$i))."end do\n";
				    $inputCode  .= (" " x ($rank+1-$i))."end do\n";
				}
				if ( grep {$_ eq "allocatable"} @{$declaration->{'attributes'}} ) {
				    $inputCode  .= " end if\n";
				    $outputCode .= " else\n";
				    $outputCode .= "  write (stateFile) .false.\n";
				    $outputCode .= " end if\n";
				}
				$stateFileUsed      = 1;
				$fgslStateFileUsed  = 1;
			    }
			}
		    } else {
			# Intrinsic type.
			if ( grep {$_ eq "pointer"} @{$declaration->{'attributes'}} ) {
			    # Pointers are currently not handled.
			} elsif ( exists($declaration->{'type'}) && $declaration->{'type'} =~ m/^\s*omp_lock_kind\s*/ ) {
			    # Do not store OpenMP lock variables.
			} elsif ( grep {$_ eq "allocatable"} @{$declaration->{'attributes'}} ) {
			    # For allocatable variables we must first store the shape so that they can be reallocated on restore.
			    my $dimensionDeclarator = join(",",map {/^dimension\s*\(([:,]+)\)/} @{$declaration->{'attributes'}});
			    my $rank = ($dimensionDeclarator =~ tr/://);
			    foreach my $variableName ( @{$declaration->{'variables'}} ) {
				next
				    if ( grep {lc($_) eq lc($variableName)} @excludes );
				$allocatablesFound  = 1;
				$dimensionalsFound  = 1;
				$stateFileUsed      = 1;
				$labelUsed          = 1;
				$outputCode        .= " if (allocated(self%".$variableName.")) then\n";
				$outputCode        .= "  if (Galacticus_Verbosity_Level() >= verbosityWorking) then\n";
				$outputCode        .= "   write (label,'(i16)') sizeof(self%".$variableName.")\n";
				$outputCode        .= "   call Galacticus_Display_Message('storing \"".$variableName."\" with size '//trim(adjustl(label))//' bytes')\n";
				$outputCode        .= "  end if\n";
				$outputCode        .= "  write (stateFile) .true.\n";
				$outputCode        .= "  write (stateFile) shape(self%".$variableName.",kind=c_size_t)\n";
				$outputCode        .= "  write (stateFile) self%".$variableName."\n";
				$outputCode        .= " else\n";
				$outputCode        .= "  write (stateFile) .false.\n";
				$outputCode        .= " end if\n";
				$inputCode         .= " read (stateFile) wasAllocated\n";
				$inputCode         .= " if (allocated(self%".$variableName.")) deallocate(self%".$variableName.")\n";
				$inputCode         .= " if (wasAllocated) then\n";
				$inputCode         .= "  call Galacticus_Display_Message('restoring \"".$variableName."\"',verbosity=verbosityWorking)\n";
				$inputCode         .= "  allocate(storedShape(".$rank."))\n";
				$inputCode         .= "  read (stateFile) storedShape\n";
				$inputCode         .= "  allocate(self%".$variableName."(".join(",",map {"storedShape(".$_.")"} 1..$rank)."))\n";
				$inputCode         .= "  deallocate(storedShape)\n";
				$inputCode         .= "  read (stateFile) self%".$variableName."\n";
				$inputCode         .= " end if\n";
			    }
			} else {
			    # Statically-sized variable.
			    foreach ( @{$declaration->{'variables'}} ) {
				(my $variableName = $_) =~ s/\s*=.*$//;
				next
				    if ( grep {lc($_) eq lc($variableName)} @excludes );
				push(@staticVariables,$variableName);
			    }
			}
		    }
		}
		# Add any variables declared in the functionClassType class.
		if ( defined($functionClassType) ) {
		    my $node = $functionClassType->{'node'}->{'firstChild'};
		    while ( $node ) {
			if ( $node->{'type'} eq "declaration" ) {			    
			    foreach my $declaration ( @{$node->{'declarations'}} ) {
				# Identify variable type.
				if ( $declaration->{'intrinsic'} eq "procedure" ) {
				    # Type-bound procedure - nothing to do.
				} elsif ( $declaration->{'intrinsic'} eq "class" || $declaration->{'intrinsic'} eq "type" ) {
				    # Look for pointers to functionClasses.
				    (my $type = $declaration->{'type'}) =~ s/\s//g;
				    if ( 
					$declaration->{'intrinsic'} eq "class"
					&&
					(grep {$_ eq "pointer"} @{$declaration   ->{'attributes'     }})
					&&
					(grep {$_ eq $type    } @{$stateStorables->{'functionClasses'}})
					) {
					# Pointer to a functionClass object.
					foreach ( @{$declaration->{'variables'}} ) {
					    $labelUsed = 1;
					    (my $variableName = $_) =~ s/\s*=.*$//;
					    $outputCode .= " if (Galacticus_Verbosity_Level() >= verbosityWorking) then\n";
					    $outputCode .= "  select type (c__ => self%".$variableName.")\n";
					    $outputCode .= "  class is (".$declaration->{'type'}.")\n";
					    $outputCode .= "   write (label,'(i16)') sizeof(c__)\n";
					    $outputCode .= "  end select\n";
					    $outputCode .= "  call Galacticus_Display_Message('storing \"".$variableName."\" with size '//trim(adjustl(label))//' bytes')\n";
					    $outputCode .= " end if\n";
					    $inputCode  .= " call Galacticus_Display_Message('restoring \"".$variableName."\"',verbosity=verbosityWorking)\n";
					    $outputCode .= " call self%".$variableName."%stateStore  (stateFile,fgslStateFile,stateOperationID)\n";
					    $inputCode  .= " call self%".$variableName."%stateRestore(stateFile,fgslStateFile,stateOperationID)\n";
					    $stateFileUsed     = 1;
					    $fgslStateFileUsed = 1;
					}
				    } elsif (
					(
					 (  grep {$_->{'type'} eq $type    } @{$stateStorables->{'stateStorables'        }})
					 ||
					 (  grep {$_           eq $type    } @{$stateStorables->{'functionClassInstances'}})
					)
					&&
					(! grep {$_           eq "pointer"} @{$declaration   ->{'attributes'            }})				       
					){
					# This is a non-pointer object which is explicitly stateStorable.
					# Validate: Currently we do not support store/restore of polymorphic functionClass objects.
					die("Galacticus::Build::SourceTree::Process::FunctionClass::Process_FunctionClass(): storing of polymorphic functionClass objects is not implemented")
					    if
					    (
					     $declaration->{'intrinsic'} eq "class"
					     &&
					     (grep {$_ eq $type} @{$stateStorables->{'functionClassInstances'}})
					    );
					my $isFunctionClass = grep {$_ eq $type} @{$stateStorables->{'functionClassInstances'}};
					# Construct code to output.
					foreach ( @{$declaration->{'variables'}} ) {
					    (my $variableName = $_) =~ s/\s*=.*$//;
					    next
						if ( grep {lc($_) eq lc($variableName)} @excludes );
					    my $rank = 0;
					    if ( grep {$_ =~ m/^dimension\s*\(/} @{$declaration->{'attributes'}} ) {
						my $dimensionDeclarator = join(",",map {/^dimension\s*\(([a-zA-Z0-9_,]+)\)/} @{$declaration->{'attributes'}});
						$rank        = ($dimensionDeclarator =~ tr/,//)+1;
						$rankMaximum = $rank
						    if ( $rank > $rankMaximum );
					    }
					    if ( grep {$_ eq "allocatable"} @{$declaration->{'attributes'}} ) {
						# For allocatable variables we must first store the shape so that they can be reallocated on restore.
						$allocatablesFound  = 1;
						$dimensionalsFound  = 1
						    if ( $rank > 0 );
						$outputCode .= " if (allocated(self%".$variableName.")) then\n";
						$outputCode .= "  write (stateFile) .true.\n";
						$outputCode .= "  write (stateFile) shape(self%".$variableName.",kind=c_size_t)\n"
						    if ( $rank > 0 );
						$inputCode  .= " read (stateFile) wasAllocated\n";
						$inputCode  .= " if (allocated(self%".$variableName.")) deallocate(self%".$variableName.")\n";
						$inputCode  .= " if (wasAllocated) then\n";
						if ( $rank > 0 ) {
						    $inputCode  .= "  allocate(storedShape(".$rank."))\n";
						    $inputCode  .= "  read (stateFile) storedShape\n";
						    $inputCode  .= "  deallocate(storedShape)\n";
						}
						if ( $declaration->{'intrinsic'} eq "class" ) {
						    $inputCode  .= "  call ".$type."ClassRestore(self%".$variableName.",stateFile)\n";
						} else {
						    $inputCode  .= "  allocate(self%".$variableName.($rank > 0 ? "(".join(",",map {"storedShape(".$_.")"} 1..$rank).")" : "").")\n";
						}
					    }
					    for(my $i=1;$i<=$rank;++$i) {
						$outputCode .= (" " x $i)."do i".$i."=1,size(self%".$variableName.",dim=".$i.")\n";
						$inputCode  .= (" " x $i)."do i".$i."=1,size(self%".$variableName.",dim=".$i.")\n";
					    }
					    my $arrayElement = $rank > 0 ? "(".join(",",map {"i".$_} 1..$rank).")" : "";
					    $labelUsed   = 1;
					    $outputCode .= " if (Galacticus_Verbosity_Level() >= verbosityWorking) then\n";
					    if ( $declaration->{'intrinsic'} eq "class" ) {
						$outputCode .= "  select type (c__ => self%".$variableName.")\n";
						$outputCode .= "  class is (".$declaration->{'type'}.")\n";
						$outputCode .= "   write (label,'(i16)') sizeof(c__".$arrayElement.")\n";
						$outputCode .= "  end select\n";
					    } else {
						$outputCode .= "   write (label,'(i16)') sizeof(self%".$variableName.")\n";
					    }
					    $outputCode .= "  call Galacticus_Display_Message('storing \"".$variableName.$arrayElement."\" with size '//trim(adjustl(label))//' bytes')\n";
					    $outputCode .= " end if\n";
					    $inputCode  .= " call Galacticus_Display_Message('restoring \"".$variableName.$arrayElement."\"',verbosity=verbosityWorking)\n";
					    $inputCode  .= (" " x $rank)." call self%".$variableName.$arrayElement."%stateRestore(stateFile,fgslStateFile".($isFunctionClass ? ",stateOperationID" : "").")\n";
					    $outputCode .= (" " x $rank)." call self%".$variableName.$arrayElement."%stateStore  (stateFile,fgslStateFile".($isFunctionClass ? ",stateOperationID" : ",storeIdentifier=".($declaration->{'intrinsic'} eq "class" ? ".true." : ".false.")).")\n";
					    for(my $i=1;$i<=$rank;++$i) {
						$outputCode .= (" " x ($rank+1-$i))."end do\n";
						$inputCode  .= (" " x ($rank+1-$i))."end do\n";
					    }
					    if ( grep {$_ eq "allocatable"} @{$declaration->{'attributes'}} ) {
						$inputCode  .= " end if\n";
						$outputCode .= " else\n";
						$outputCode .= "  write (stateFile) .false.\n";
						$outputCode .= " end if\n";
					    }
					    $stateFileUsed      = 1;
					    $fgslStateFileUsed  = 1;
					}
				    }
				} else {
				    # Intrinsic type.
				    if ( grep {$_ eq "pointer"} @{$declaration->{'attributes'}} ) {
					# Pointers are currently not handled.
				    } elsif ( exists($declaration->{'type'}) && defined($declaration->{'type'}) && $declaration->{'type'} =~ m/^\s*omp_lock_kind\s*/ ) {
					# Do not store OpenMP lock variables.
				    } elsif ( grep {$_ eq "allocatable"} @{$declaration->{'attributes'}} ) {
					# For allocatable variables we must first store the shape so that they can be reallocated on restore.
					my $dimensionDeclarator = join(",",map {/^dimension\s*\(([:,]+)\)/} @{$declaration->{'attributes'}});
					my $rank = ($dimensionDeclarator =~ tr/://);
					foreach my $variableName ( @{$declaration->{'variables'}} ) {
					    next
						if ( grep {lc($_) eq lc($variableName)} @excludes );
					    $allocatablesFound  = 1;
					    $dimensionalsFound  = 1;
					    $stateFileUsed      = 1;
					    $labelUsed          = 1;
					    $outputCode        .= " if (allocated(self%".$variableName.")) then\n";
					    $outputCode        .= "  if (Galacticus_Verbosity_Level() >= verbosityWorking) then\n";
					    $outputCode        .= "   write (label,'(i16)') sizeof(self%".$variableName.")\n";
					    $outputCode        .= "   call Galacticus_Display_Message('storing \"".$variableName."\" with size '//trim(adjustl(label))//' bytes')\n";
					    $outputCode        .= "  end if\n";
					    $outputCode        .= "  write (stateFile) .true.\n";
					    $outputCode        .= "  write (stateFile) shape(self%".$variableName.",kind=c_size_t)\n";
					    $outputCode        .= "  write (stateFile) self%".$variableName."\n";
					    $outputCode        .= " else\n";
					    $outputCode        .= "  write (stateFile) .false.\n";
					    $outputCode        .= " end if\n";
					    $inputCode         .= " read (stateFile) wasAllocated\n";
					    $inputCode         .= " if (allocated(self%".$variableName.")) deallocate(self%".$variableName.")\n";
					    $inputCode         .= " if (wasAllocated) then\n";
					    $inputCode         .= "  call Galacticus_Display_Message('restoring \"".$variableName."\"',verbosity=verbosityWorking)\n";
					    $inputCode         .= "  allocate(storedShape(".$rank."))\n";
					    $inputCode         .= "  read (stateFile) storedShape\n";
					    $inputCode         .= "  allocate(self%".$variableName."(".join(",",map {"storedShape(".$_.")"} 1..$rank)."))\n";
					    $inputCode         .= "  deallocate(storedShape)\n";
					    $inputCode         .= "  read (stateFile) self%".$variableName."\n";
					    $inputCode         .= " end if\n";
					}
				    } else {
					# Statically-sized variable.
					foreach ( @{$declaration->{'variables'}} ) {
					    (my $variableName = $_) =~ s/\s*=.*$//;
					    next
						if ( grep {lc($_) eq lc($variableName)} @excludes );
					    my $store = 1;
					    if ( exists($class->{'stateStorable'}) && exists($class->{'stateStorable'}->{'restoreTo'}) ) {
						foreach ( &List::ExtraUtils::as_array($class->{'stateStorable'}->{'restoreTo'}) ) {
						    my @variables = split(/\s*,\s*/,$_->{'variables'});
						    if ( grep {lc($_) eq lc($variableName)} @variables ) {
							$store = 0;
							$inputCode .= " self%".$variableName."=".$_->{'state'}."\n";
						    }
						}
					    }
					    push(@staticVariables,$variableName)
						if ( $store );
					}
				    }
				}
				$hasCustomStateStore   = 1
				    if
				    (
				     $declaration->{'intrinsic'} eq "procedure"
				     &&
				     $declaration->{'variables'}->[0] =~ m/^stateStore=>/
				    );
				$hasCustomStateRestore = 1
				    if
				    (
				     $declaration->{'intrinsic'} eq "procedure"
				     &&
				     $declaration->{'variables'}->[0] =~ m/^stateRestore=>/
				    );
			    }
			}
			$node = $node->{'type'} eq "contains" ? $node->{'firstChild'} : $node->{'sibling'};
		    }
		}
		# Add code to method.
		$stateFileUsed = 1
		    if ( scalar(@staticVariables) > 0 );
		if ( $hasCustomStateStore   ) {
		    # The class has its own state store function, so we should never arrive at this point in the code.
		    $stateStoreCode .= " call Galacticus_Error_Report('custom state store function exists - this should not happen'//".&Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($nonAbstractClass->{'node'},$nonAbstractClass->{'node'}->{'line'}).")\n";
		    $stateStoreModules{'Galacticus_Error'} = 1;
		} else {
		    foreach ( @staticVariables ) {
			$labelUsed       = 1;
			$stateStoreCode .= " if (Galacticus_Verbosity_Level() >= verbosityWorking) then\n";
			$stateStoreCode .= "  write (label,'(i16)') sizeof(self%".$_.")\n";
			$stateStoreCode .= "  call Galacticus_Display_Message('storing \"".$_."\" with size '//trim(adjustl(label))//' bytes')\n";
			$stateStoreCode .= " end if\n";
		    }
		    $stateStoreCode .= " write (stateFile) ".join(", &\n  & ",map {"self%".$_} @staticVariables)."\n"
			if ( scalar(@staticVariables) > 0 );
		    $stateStoreCode .= $outputCode
			if ( defined($outputCode) );
		}		
		if ( $hasCustomStateRestore ) {
		    # The class has its own state store function, so we should never arrive at this point in the code.
		    $stateRestoreCode .= " call Galacticus_Error_Report('custom state restore function exists - this should not happen'//".&Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($nonAbstractClass->{'node'},$nonAbstractClass->{'node'}->{'line'}).")\n";
		    $stateRestoreModules{'Galacticus_Error'} = 1;
		} else {
		    foreach ( @staticVariables ) {
			$stateRestoreCode .= " call Galacticus_Display_Message('restoring \"".$_."\"',verbosity=verbosityWorking)\n";
		    }
		    $stateRestoreCode .= " read (stateFile) ".join(", &\n  & ",map {"self%".$_} @staticVariables)."\n"
			if ( scalar(@staticVariables) > 0 );
		    $stateRestoreCode .= $inputCode
			if ( defined($inputCode) );
		}		
	    }
	    $stateStoreCode   .= "end select\n";
	    $stateStoreCode   .= "call Galacticus_Display_Unindent('done',verbosity=verbosityWorking)\n";
	    $stateStoreCode   .= "return\n";
	    $stateRestoreCode .= "end select\n";
	    $stateRestoreCode .= "call Galacticus_Display_Unindent('done',verbosity=verbosityWorking)\n";
	    $stateRestoreCode .= "return\n";
	    unless ( $fgslStateFileUsed ) {
		push(@outputUnusedVariables,"fgslStateFile");
		push(@inputUnusedVariables ,"fgslStateFile");
	    }
	    unless ( $stateFileUsed     ) {
		push(@outputUnusedVariables,"stateFile"    );
		push(@inputUnusedVariables ,"stateFile"    );
	    }
	    push(@outputUnusedVariables,"label")
		unless ( $labelUsed );
	    $stateStoreCode   =
		($rankMaximum > 0 ? " integer :: ".join(", ",map {"i".$_} 1..$rankMaximum)."\n" : "").
		(@outputUnusedVariables ? " !GCC\$ attributes unused :: ".join(", ",@outputUnusedVariables)."\n" : "").
		$stateStoreCode  ;
	    $stateRestoreCode = 
		($rankMaximum > 0 ? " integer :: ".join(", ",map {"i".$_} 1..$rankMaximum)."\n" : "").
		(@inputUnusedVariables ? " !GCC\$ attributes unused :: ".join(", ",@inputUnusedVariables)."\n" : "").
		$stateRestoreCode;		
	    if ( $allocatablesFound ) {
		$stateRestoreCode = ($dimensionalsFound ? "integer(c_size_t), allocatable, dimension(:) :: storedShape\n"  : "").
		    ($allocatablesFound ? "logical                                      :: wasAllocated\n" : "").
		    $stateRestoreCode;
	    }
	    $stateStoreCode   = " character(len=16) :: label\n".$stateStoreCode  ;
	    $methods{'stateStore'} = 
	    {
		description => "Store the state of this object to file.",
		type        => "void",
		pass        => "yes",
		modules     => join(" ",keys(%stateStoreModules)),
		argument    => [ "integer, intent(in   ) :: stateFile", "type(fgsl_file), intent(in   ) :: fgslStateFile", "integer(c_size_t), intent(in   ) :: stateOperationID"  ],
		code        => $stateStoreCode
	    };
	    $methods{'stateRestore'} = 
	    {
		description => "Restore the state of this object from file.",
		type        => "void",
		pass        => "yes",
		modules     => join(" ",keys(%stateRestoreModules)),
		argument    => [ "integer, intent(in   ) :: stateFile", "type(fgsl_file), intent(in   ) :: fgslStateFile", "integer(c_size_t), intent(in   ) :: stateOperationID"  ],
		code        => $stateRestoreCode
	    };
	    # Initialize new nodes.
	    my $preContains = 
		[
		 {
		     type       => "code" ,
		     content    => ""     ,
		     firstChild => undef(),
		     sibling    => undef(),
		     parent     => undef(),
		     source     => "Galacticus::Build::SourceTree::Process::FunctionClass::Process_FunctionClass()",
		     line       => 1
		 }
		];
	    my $postContains =
		[
		 {
		     type       => "code" ,
		     content    => ""     ,
		     firstChild => undef(),
		     sibling    => undef(),
		     parent     => undef(),
		     source     => "Galacticus::Build::SourceTree::Process::FunctionClass::Process_FunctionClass()",
		     line       => 1
		 }
		];

	    # Add variable tracking module initialization status.
	    $preContains->[0]->{'content'} .= "   logical, private :: ".$directive->{'name'}."Initialized=.false.\n\n";

	    # Generate the base class.
	    &Galacticus::Build::SourceTree::SetVisibility($node->{'parent'},$directive->{'name'}."Class","public");
	    &Galacticus::Build::SourceTree::SetVisibility($node->{'parent'},$directive->{'name'}        ,"public");
	    my $extends = exists($directive->{'extends'}) ? $directive->{'extends'} : "functionClass";
	    $preContains->[0]->{'content'} .= "   type, extends(".$extends.") :: ".$directive->{'name'}."Class\n";
	    $preContains->[0]->{'content'} .= "    private\n";
	    my $usesNode =
	    {
		type      => "moduleUse",
		moduleUse =>
		{
		    Function_Classes =>
		    {
			intrinsic => 0,
			all       => 1
		    }
		}
	    };
            $usesNode->{'moduleUse'}->{'ISO_C_Binding'} = {intrinsic => 1, all => 1}
                if ( $debugging );
            $preContains->[0]->{'content'} .= "    integer(c_size_t) :: stateOperationID=0\n";
            $usesNode->{'moduleUse'}->{'ISO_C_Binding'} =
                {
		    intrinsic => 1,
		    all       => 1
		};
            &Galacticus::Build::SourceTree::Parse::ModuleUses::AddUses($node->{'parent'},$usesNode);
	    foreach ( &List::ExtraUtils::as_array($directive->{'data'}) ) {
		if ( reftype($_) ) {
		    $_->{'scope'} = "self"
			unless ( exists($_->{'scope'}) );
		    $preContains->[0]->{'content'} .= $_->{'content'}."\n"
			if (  $_->{'scope'} eq "self" );
		} else {
		    $preContains->[0]->{'content'} .= $_."\n";
		}
	    }
	    $preContains->[0]->{'content'} .= "    contains\n";
	    $preContains->[0]->{'content'} .= "    !@ <objectMethods>\n";
	    $preContains->[0]->{'content'} .= "    !@   <object>".$directive->{'name'}."Class</object>\n";
	    my $generics;
	    foreach my $methodName ( keys(%methods) ) {
		my $method = $methods{$methodName};
		my $argumentList = "";
		my @arguments;
		if ( exists($method->{'argument'}) ) {
		    if ( UNIVERSAL::isa($method->{'argument'},"ARRAY") ) {
			push(@arguments,@{$method->{'argument'}});
		    } else {
			push(@arguments,  $method->{'argument'} );
		    }
		}
		my $separator = "";
		foreach my $argument ( @arguments ) {
		    foreach my $intrinsic ( keys(%Fortran::Utils::intrinsicDeclarations) ) {
			my $declarator = $Fortran::Utils::intrinsicDeclarations{$intrinsic};
			if ( my @matches = $argument =~ m/$declarator->{'regEx'}/ ) {
			    my $intrinsicName =                          $declarator->{'intrinsic' }  ;
			    my $type          =                 $matches[$declarator->{'type'      }] ;
			    my $attributeList =                 $matches[$declarator->{'attributes'}] ;
			    $attributeList =~ s/^\s*,?\s*//;
			    $attributeList =~ s/\s*$//;
			    my @attributes = &Fortran::Utils::Extract_Variables($attributeList, keepQualifiers => 1, removeSpaces => 1);
			    my @variables     = split(/\s*,\s*/,$matches[$declarator->{'variables' }]);
			    foreach my $variable ( @variables ) {
				$argumentList .= $separator."\\textcolor{red}{\\textless ".$intrinsicName;
				$argumentList .= "(".$type.")"
				    if ( defined($type) );
				$argumentList .= "\\textgreater} ".$variable;
				foreach my $attribute ( @attributes ) {
				    $argumentList .= "\\argin"
					if ( $attribute eq "intent(in)" );
				    $argumentList .= "\\argout"
					if ( $attribute eq "intent(out)" );
				    $argumentList .= "\\arginout"
					if ( $attribute eq "intent(inout)" );
				}
				$separator     = ",";
			    }
			}
		    }
		}
		$preContains->[0]->{'content'} .= "    !@   <objectMethod>\n";
		$preContains->[0]->{'content'} .= "    !@     <method>".$methodName."</method>\n";
		$preContains->[0]->{'content'} .= "    !@     <type>".latex_encode($method->{'type'})."</type>\n";
		$preContains->[0]->{'content'} .= "    !@     <arguments>".latex_encode($argumentList)."</arguments>\n";
		$preContains->[0]->{'content'} .= "    !@     <description>".$method->{'description'}."</description>\n";
		$preContains->[0]->{'content'} .= "    !@   </objectMethod>\n";
		if ( exists($directive->{'generic'}) ) {
		    foreach my $generic ( &List::ExtraUtils::as_array($directive->{'generic'}) ) {
			if ( grep {$_ eq $methodName} &List::ExtraUtils::as_array($generic->{'method'}) ) {
			    # This method is part of a generic method, store relevant information.
			    $generics->{$generic->{'name'}}->{'type'} = $method->{'type'};
			    push(@{ $generics->{$generic->{'name'}}->{'description'}},$method->{'description'});
			    push(@{ $generics->{$generic->{'name'}}->{'argumentList'}},latex_encode($argumentList));
			}
		    }
		}
	    }
	    foreach my $generic ( &List::ExtraUtils::hashList($generics, keyAs => "name") ) {
		$preContains->[0]->{'content'} .= "    !@   <objectMethod>\n";
		$preContains->[0]->{'content'} .= "    !@     <method>".$generic->{'name'}."</method>\n";
		$preContains->[0]->{'content'} .= "    !@     <type>".latex_encode($generic->{'type'})."</type>\n";
		$preContains->[0]->{'content'} .= "    !@     <arguments>".join(" | ",@{$generic->{'argumentList'}})."</arguments>\n";
		$preContains->[0]->{'content'} .= "    !@     <description>".join(" | ",@{$generic->{'description'}})."</description>\n";
		$preContains->[0]->{'content'} .= "    !@   </objectMethod>\n";
	    }
	    $preContains->[0]->{'content'} .= "    !@ </objectMethods>\n";
	    my $methodTable = Text::Table->new(
		{
		    is_sep => 1,
		    body   => "    procedure"
		},
		{
		    align  => "left"
		},
		{
		    is_sep => 1,
		    body   => " :: "
		},
		{
		    align  => "left"
		},
		{
		    is_sep => 1,
		    body   => " => ",
		},
		{
		    align  => "left"
		}
		);    
	    foreach ( keys(%methods) ) {
		my $method = $methods{$_};
		my $extension = "__";
		$extension = ""
		    if ( exists($method->{'code'}) );
		$methodTable->add("",$_,$directive->{'name'}.ucfirst($_).$extension);
	    }
            $preContains->[0]->{'content'} .= $methodTable->table();
            if ( exists($directive->{'generic'}) ) {
		my $genericTable = Text::Table->new(
		    {
			is_sep => 1,
			body   => "    generic :: "
		    },
		    {
			align  => "left"
		    },
		    {
			is_sep => 1,
			body   => " => ",
		    },
		    {
			align  => "left"
		    }
		    );    
		foreach ( &List::ExtraUtils::as_array($directive->{'generic'}) ) {
		    $genericTable->add($_->{'name'},join(", ",&List::ExtraUtils::as_array($_->{'method'})));
		}
		$preContains->[0]->{'content'} .= $genericTable->table();
            }
	    $preContains->[0]->{'content'} .= "   end type ".$directive->{'name'}."Class\n\n";
	    # Insert any module-scope class content.
	    foreach ( &List::ExtraUtils::as_array($directive->{'data'}) ) {
		if ( reftype($_) ) {
		    if ( exists($_->{'scope'}) && $_->{'scope'} eq "module" ) {
			$preContains->[0]->{'content'} .= $_->{'content'}."\n";
			if ( exists($_->{'threadprivate'}) && $_->{'threadprivate'} eq "yes" && $_->{'content'} =~ m/::\s*(.*)$/ ) {
			    $preContains->[0]->{'content'} .= "   !\$omp threadprivate(".$1.")\n";
			}
		    }
		}
	    }

	    # Generate class constructors
	    $preContains->[0]->{'content'} .= "   interface ".$directive->{'name'}."\n";
	    $preContains->[0]->{'content'} .= "    module procedure ".$directive->{'name'}."CnstrctrDflt\n";
	    $preContains->[0]->{'content'} .= "    module procedure ".$directive->{'name'}."CnstrctrPrmtrs\n";
	    $preContains->[0]->{'content'} .= "   end interface ".$directive->{'name'}."\n";
	    # Add method name parameter.
	    $preContains->[0]->{'content'} .= "   ! Method name parameter.\n";
	    $preContains->[0]->{'content'} .= "   type(varying_string) :: ".$directive->{'name'}."Method\n\n";
	    my $nameUsesNode =
	    {
		type      => "moduleUse",
		moduleUse =>
		{
		    ISO_Varying_String =>
		    {
			intrinsic => 0,
			all       => 1
		    }
		}
	    };
	    &Galacticus::Build::SourceTree::Parse::ModuleUses::AddUses($node->{'parent'},$nameUsesNode);
	    if ( $tree->{'type'} eq "file" ) {
		(my $fileName = $tree->{'name'}) =~ s/\.F90$/.p/;
		open(my $parametersFile,">>".$ENV{'BUILDPATH'}."/".$fileName);
		print $parametersFile $directive->{'name'}."Method\n";
		close($parametersFile);
	    }
	    # Add default implementation.
	    $preContains->[0]->{'content'} .= "   ! Default ".$directive->{'name'}." object.\n";
	    $preContains->[0]->{'content'} .= "   class(".$directive->{'name'}."Class), private , pointer :: ".$directive->{'name'}."Default => null()\n";
	    $preContains->[0]->{'content'} .= "   !\$omp threadprivate(".$directive->{'name'}."Default)\n";
	    $preContains->[0]->{'content'} .= "\n";
	    # Create default constructor.
	    $postContains->[0]->{'content'} .= "   function ".$directive->{'name'}."CnstrctrDflt()\n";
	    $postContains->[0]->{'content'} .= "      !% Return a pointer to the default {\\normalfont \\ttfamily ".$directive->{'name'}."} object.\n";
	    $postContains->[0]->{'content'} .= "      implicit none\n";
	    $postContains->[0]->{'content'} .= "      class(".$directive->{'name'}."Class), pointer :: ".$directive->{'name'}."CnstrctrDflt\n\n";
	    $postContains->[0]->{'content'} .= "      if (.not.associated(".$directive->{'name'}."Default)) call ".$directive->{'name'}."Initialize()\n";
	    $postContains->[0]->{'content'} .= "      ".$directive->{'name'}."CnstrctrDflt => ".$directive->{'name'}."Default\n";
	    $postContains->[0]->{'content'} .= "      return\n";
	    $postContains->[0]->{'content'} .= "   end function ".$directive->{'name'}."CnstrctrDflt\n\n";
	    # Create XML constructor.
	    $postContains->[0]->{'content'} .= "   function ".$directive->{'name'}."CnstrctrPrmtrs(parameters,copyInstance,parameterName)\n";
	    $postContains->[0]->{'content'} .= "      !% Return a pointer to a newly created {\\normalfont \\ttfamily ".$directive->{'name'}."} object as specified by the provided parameters.\n";
	    $postContains->[0]->{'content'} .= "      use Input_Parameters\n";
	    $postContains->[0]->{'content'} .= "      use Galacticus_Error\n";
	    $postContains->[0]->{'content'} .= "      implicit none\n";
	    $postContains->[0]->{'content'} .= "      class    (".$directive->{'name'}."Class), pointer :: ".$directive->{'name'}."CnstrctrPrmtrs\n";
	    $postContains->[0]->{'content'} .= "      type     (inputParameters), intent(inout)           :: parameters\n";
	    $postContains->[0]->{'content'} .= "      integer                   , intent(in   ), optional :: copyInstance\n";
	    $postContains->[0]->{'content'} .= "      character(len=*          ), intent(in   ), optional :: parameterName\n";
	    $postContains->[0]->{'content'} .= "      type     (inputParameters)                          :: subParameters\n";
	    $postContains->[0]->{'content'} .= "      type     (inputParameter ), pointer                 :: parameterNode\n"
                if ( exists($directive->{'default'}) );
	    $postContains->[0]->{'content'} .= "      type     (varying_string )                          :: message      , instanceName, parameterName_\n\n";
	    $postContains->[0]->{'content'} .= "      if (present(parameterName)) then\n";
	    $postContains->[0]->{'content'} .= "        parameterName_=parameterName\n";
	    $postContains->[0]->{'content'} .= "      else\n";
	    $postContains->[0]->{'content'} .= "        parameterName_='".$directive->{'name'}."Method'\n";
	    $postContains->[0]->{'content'} .= "      end if\n";
	    if ( exists($directive->{'default'}) ) {
	        $postContains->[0]->{'content'} .= "      if (parameterName_ == '".$directive->{'name'}."Method' .and. (.not.present(copyInstance) .or. copyInstance == 1) .and. .not.parameters%isPresent(char(parameterName_))) then\n";
	        $postContains->[0]->{'content'} .= "        call parameters%addParameter('".$directive->{'name'}."Method','".$directive->{'default'}."')\n";
	        $postContains->[0]->{'content'} .= "        parameterNode => parameters%node('".$directive->{'name'}."Method',requireValue=.true.)\n";
		$postContains->[0]->{'content'} .= "        subParameters=parameters%subParameters(char(parameterName_))\n";
    		$postContains->[0]->{'content'} .= "        allocate(".$directive->{'name'}.ucfirst($directive->{'default'})." :: ".$directive->{'name'}."CnstrctrPrmtrs)\n";
		$postContains->[0]->{'content'} .= "        select type (".$directive->{'name'}."CnstrctrPrmtrs)\n";
		$postContains->[0]->{'content'} .= "          type is (".$directive->{'name'}.ucfirst($directive->{'default'}).")\n";
		$postContains->[0]->{'content'} .= "            call debugStackPush(loc(".$directive->{'name'}."CnstrctrPrmtrs))\n"
		    if ( $debugging );
		$postContains->[0]->{'content'} .= "            ".$directive->{'name'}."CnstrctrPrmtrs=".$directive->{'name'}.ucfirst($directive->{'default'})."(subParameters)\n";
		$postContains->[0]->{'content'} .= "            call debugStackPop()\n"
		    if ( $debugging );
		$postContains->[0]->{'content'} .= "         end select\n";
                $postContains->[0]->{'content'} .= "         call parameterNode%objectSet(".$directive->{'name'}."CnstrctrPrmtrs)\n";
                $postContains->[0]->{'content'} .= "      else\n";
            }
	    $postContains->[0]->{'content'} .= "      call parameters%value(char(parameterName_),instanceName,copyInstance=copyInstance)\n";
	    $postContains->[0]->{'content'} .= "      subParameters=parameters%subParameters(char(parameterName_),copyInstance=copyInstance)\n";
	    $postContains->[0]->{'content'} .= "      select case (char(instanceName))\n";
	    foreach my $class ( @nonAbstractClasses ) {
		(my $name = $class->{'name'}) =~ s/^$directive->{'name'}//;
		$name = lcfirst($name)
		    unless ( $name =~ m/^[A-Z]{2,}/ );
		$postContains->[0]->{'content'} .= "     case ('".$name."')\n";
		$postContains->[0]->{'content'} .= "        allocate(".$class->{'name'}." :: ".$directive->{'name'}."CnstrctrPrmtrs)\n";
		$postContains->[0]->{'content'} .= "        select type (".$directive->{'name'}."CnstrctrPrmtrs)\n";
		$postContains->[0]->{'content'} .= "          type is (".$class->{'name'}.")\n";
		$postContains->[0]->{'content'} .= "            call debugStackPush(loc(".$directive->{'name'}."CnstrctrPrmtrs))\n"
		    if ( $debugging );
		$postContains->[0]->{'content'} .= "            ".$directive->{'name'}."CnstrctrPrmtrs=".$class->{'name'}."(subParameters)\n";
		$postContains->[0]->{'content'} .= "            call debugStackPop()\n"
		    if ( $debugging );
		$postContains->[0]->{'content'} .= "         end select\n";
	    }
	    $postContains->[0]->{'content'} .= "      case default\n";
	    $postContains->[0]->{'content'} .= "         message='Unrecognized type \"'//trim(instanceName)//'\" Available options are:'\n";
	    my @classNames;
	    push(@classNames,$_->{'name'})
		foreach ( @nonAbstractClasses );
	    foreach ( sort(@classNames) ) {
		(my $name = $_) =~ s/^$directive->{'name'}//;
		$name = lcfirst($name)
		    unless ( $name =~ m/^[A-Z]{2,}/ );
		$postContains->[0]->{'content'} .= "         message=message//char(10)//'   -> ".$name."'\n";
	    }
	    $postContains->[0]->{'content'} .= "         call Galacticus_Error_Report(message//".&Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($node,$node->{'line'}).")\n";
            $postContains->[0]->{'content'} .= "      end select\n";
            $postContains->[0]->{'content'} .= "      end if\n"
                if ( exists($directive->{'default'}) );
 	    $postContains->[0]->{'content'} .= "      return\n";
	    $postContains->[0]->{'content'} .= "   end function ".$directive->{'name'}."CnstrctrPrmtrs\n\n";
	    
	    # Insert class code.
	    foreach my $class ( @classes ) {
		&Galacticus::Build::SourceTree::SetVisibility($node->{'parent'},$class->{'type'},"public");
		my $classTree = $class->{'tree'};
		my $classNode = $classTree->{'firstChild'};
		my $contained = 0;
		while ( $classNode ) {
		    if ( $classNode->{'type'} eq "contains" ) {
			$classNode = $classNode->{'firstChild'};
			$contained = 1;
		    }
		    if ( $contained ) {
			push(@{$postContains},$classNode);
		    } else {
			if ( $classNode->{'type'} eq "moduleUse" ) {
			    &Galacticus::Build::SourceTree::Parse::ModuleUses::AddUses($node->{'parent'},$classNode);
			} else {			    
			    push(@{$preContains },$classNode);
			}
		    }
		    $classNode = $classNode->{'sibling'};
		}
	    }
	    
	    # Create initialization function.
	    $postContains->[0]->{'content'} .= "   subroutine ".$directive->{'name'}."Initialize()\n";
	    $postContains->[0]->{'content'} .= "      !% Initialize the default {\\normalfont \\ttfamily ".$directive->{'name'}."} object.\n";
	    $postContains->[0]->{'content'} .= "      use ISO_Varying_String\n";
	    $postContains->[0]->{'content'} .= "      use Input_Parameters\n";
	    $postContains->[0]->{'content'} .= "      use Galacticus_Error\n";
	    $postContains->[0]->{'content'} .= "      use IO_HDF5\n";
	    $postContains->[0]->{'content'} .= "      implicit none\n";
	    $postContains->[0]->{'content'} .= "      type   (inputParameters) :: subParameters\n";
	    $postContains->[0]->{'content'} .= "      type   (varying_string ) :: message\n";
	    $postContains->[0]->{'content'} .= "      !\$omp critical (".$directive->{'name'}."Initialization)\n";
	    $postContains->[0]->{'content'} .= "      if (.not.".$directive->{'name'}."Initialized) then\n";
	    $postContains->[0]->{'content'} .= "         !@ <inputParameter>\n";
	    $postContains->[0]->{'content'} .= "         !@   <name>".$directive->{'name'}."Method</name>\n";
	    $postContains->[0]->{'content'} .= "         !@   <defaultValue>".$directive->{'default'}."</defaultValue>\n"
              if ( exists($directive->{'default'}) );
	    $postContains->[0]->{'content'} .= "         !@   <attachedTo>module</attachedTo>\n";
	    $postContains->[0]->{'content'} .= "         !@   <description>\n";
	    $postContains->[0]->{'content'} .= "         !@     The method to be used for {\\normalfont \\ttfamily ".$directive->{'name'}."}.\n";
	    $postContains->[0]->{'content'} .= "         !@   </description>\n";
	    $postContains->[0]->{'content'} .= "         !@   <type>string</type>\n";
	    $postContains->[0]->{'content'} .= "         !@   <cardinality>1</cardinality>\n";
	    $postContains->[0]->{'content'} .= "         !@ </inputParameter>\n";
	    $postContains->[0]->{'content'} .= "         call globalParameters%value('".$directive->{'name'}."Method',".$directive->{'name'}."Method";
	    $postContains->[0]->{'content'} .= ",defaultValue=var_str('".$directive->{'default'}."')"
               if ( exists($directive->{'default'}) );
 	    $postContains->[0]->{'content'} .= ")\n";
	    $postContains->[0]->{'content'} .= "         ".$directive->{'name'}."Initialized=.true.\n";
	    $postContains->[0]->{'content'} .= "      end if\n";
	    $postContains->[0]->{'content'} .= "      subParameters=globalParameters%subParameters('".$directive->{'name'}."Method',requirePresent=.false.)\n";
	    $postContains->[0]->{'content'} .= "      select case (char(".$directive->{'name'}."Method))\n";
	    foreach my $class ( @nonAbstractClasses ) {
		(my $name = $class->{'name'}) =~ s/^$directive->{'name'}//;
		$name = lcfirst($name)
		    unless ( $name =~ m/^[A-Z]{2,}/ );
		$postContains->[0]->{'content'} .= "     case ('".$name."')\n";
		$postContains->[0]->{'content'} .= "        allocate(".$class->{'name'}." :: ".$directive->{'name'}."Default)\n";
		$postContains->[0]->{'content'} .= "        select type (".$directive->{'name'}."Default)\n";
		$postContains->[0]->{'content'} .= "        type is (".$class->{'name'}.")\n";

		$postContains->[0]->{'content'} .= "        !# <referenceConstruct ownerLoc=\"module:".$node->{'parent'}->{'name'}."\" object=\"".$directive->{'name'}."Default\" constructor=\"".$class->{'name'}."(subParameters)\" />\n";
		$postContains->[0]->{'content'} .= "        end select\n";
		$postContains->[0]->{'content'} .= "        call ".$directive->{'name'}."Default%autoHook()\n"
		    if ( grep {exists($_->{'autoHook'}) && $_->{'autoHook'} eq "yes"} @classes );
	    }
	    $postContains->[0]->{'content'} .= "      case default\n";
	    $postContains->[0]->{'content'} .= "         message='Unrecognized option for [".$directive->{'name'}."Method](='//".$directive->{'name'}."Method//'). Available options are:'\n";
	    foreach ( sort(@classNames) ) {
		(my $name = $_) =~ s/^$directive->{'name'}//;
		$name = lcfirst($name)
		    unless ( $name =~ m/^[A-Z]{2,}/ );
		$postContains->[0]->{'content'} .= "        message=message//char(10)//'   -> ".$name."'\n";
	    }
	    $postContains->[0]->{'content'} .= "         call Galacticus_Error_Report(message//".&Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($node,$node->{'line'}).")\n";
	    $postContains->[0]->{'content'} .= "      end select\n";
            $postContains->[0]->{'content'} .= "      ".$directive->{'name'}."Default%isDefaultOfClass=.true.\n";
	    $postContains->[0]->{'content'} .= "      !\$omp end critical (".$directive->{'name'}."Initialization)\n";
	    $postContains->[0]->{'content'} .= "      return\n";
	    $postContains->[0]->{'content'} .= "   end subroutine ".$directive->{'name'}."Initialize\n\n";

	    # Create global state store/restore functions.
	    &Galacticus::Build::SourceTree::SetVisibility($node->{'parent'},$directive->{'name'}.$_,"public")
		foreach ( "DoStateStore", "DoStateRetrieve" );
	    $postContains->[0]->{'content'} .= "  !# <galacticusStateStoreTask>\n";
	    $postContains->[0]->{'content'} .= "  !#  <unitName>".$directive->{'name'}."DoStateStore</unitName>\n";
	    $postContains->[0]->{'content'} .= "  !# </galacticusStateStoreTask>\n";
	    $postContains->[0]->{'content'} .= "  subroutine ".$directive->{'name'}."DoStateStore(stateFile,fgslStateFile,stateOperationID)\n";
	    $postContains->[0]->{'content'} .= "    !% Store the state to file.\n";
	    $postContains->[0]->{'content'} .= "    use, intrinsic :: ISO_C_Binding     , only : c_size_t\n";
	    $postContains->[0]->{'content'} .= "    use            :: FGSL              , only : fgsl_file\n";
	    $postContains->[0]->{'content'} .= "    use            :: ISO_Varying_String, only : var_str\n";
	    $postContains->[0]->{'content'} .= "    use            :: String_Handling   , only : operator(//)\n";
	    $postContains->[0]->{'content'} .= "    use            :: Galacticus_Display, only : Galacticus_Display_Message, verbosityWorking\n";
	    $postContains->[0]->{'content'} .= "    implicit none\n";
	    $postContains->[0]->{'content'} .= "    integer           , intent(in   ) :: stateFile\n";
	    $postContains->[0]->{'content'} .= "    integer(c_size_t ), intent(in   ) :: stateOperationID\n";
	    $postContains->[0]->{'content'} .= "    type   (fgsl_file), intent(in   ) :: fgslStateFile\n";
	    $postContains->[0]->{'content'} .= "    if (associated(".$directive->{'name'}."Default)) then\n";
	    $postContains->[0]->{'content'} .= "     write (stateFile) .true.\n";
	    $postContains->[0]->{'content'} .= "     call Galacticus_Display_Message(var_str('storing default object of \""  .$directive->{'name'}."\" class [position: ')//FTell(stateFile)//']',verbosity=verbosityWorking)\n";
	    $postContains->[0]->{'content'} .= "     call ".$directive->{'name'}."Default%stateStore(stateFile,fgslStateFile,stateOperationID)\n";
	    $postContains->[0]->{'content'} .= "    else\n";
	    $postContains->[0]->{'content'} .= "     write (stateFile) .false.\n";
	    $postContains->[0]->{'content'} .= "     call Galacticus_Display_Message(var_str('skipping default object of \""  .$directive->{'name'}."\" class [position: ')//FTell(stateFile)//']',verbosity=verbosityWorking)\n";
	    $postContains->[0]->{'content'} .= "    end if\n";
	    $postContains->[0]->{'content'} .= "    return\n";
	    $postContains->[0]->{'content'} .= "  end subroutine ".$directive->{'name'}."DoStateStore\n\n";
	    $postContains->[0]->{'content'} .= "  !# <galacticusStateRetrieveTask>\n";
	    $postContains->[0]->{'content'} .= "  !#  <unitName>".$directive->{'name'}."DoStateRetrieve</unitName>\n";
	    $postContains->[0]->{'content'} .= "  !# </galacticusStateRetrieveTask>\n";
	    $postContains->[0]->{'content'} .= "  subroutine ".$directive->{'name'}."DoStateRetrieve(stateFile,fgslStateFile,stateOperationID)\n";
	    $postContains->[0]->{'content'} .= "    !% Retrieve the state from file.\n";
	    $postContains->[0]->{'content'} .= "    use, intrinsic :: ISO_C_Binding     , only : c_size_t\n";
	    $postContains->[0]->{'content'} .= "    use            :: FGSL              , only : fgsl_file\n";
	    $postContains->[0]->{'content'} .= "    use            :: ISO_Varying_String, only : var_str\n";
	    $postContains->[0]->{'content'} .= "    use            :: String_Handling   , only : operator(//)\n";
	    $postContains->[0]->{'content'} .= "    use            :: Galacticus_Display, only : Galacticus_Display_Message, verbosityWorking\n";
	    $postContains->[0]->{'content'} .= "    implicit none\n";
	    $postContains->[0]->{'content'} .= "    integer           , intent(in   ) :: stateFile\n";
	    $postContains->[0]->{'content'} .= "    integer(c_size_t ), intent(in   ) :: stateOperationID\n";
	    $postContains->[0]->{'content'} .= "    type   (fgsl_file), intent(in   ) :: fgslStateFile\n";
	    $postContains->[0]->{'content'} .= "    class  (".$directive->{'name'}."Class), pointer :: default\n";
	    $postContains->[0]->{'content'} .= "    logical                                         :: initialized\n\n";
	    $postContains->[0]->{'content'} .= "    read (stateFile) initialized\n";
	    $postContains->[0]->{'content'} .= "    if (initialized) then\n";
	    $postContains->[0]->{'content'} .= "     call Galacticus_Display_Message(var_str('restoring default object of \""  .$directive->{'name'}."\" class [position: ')//FTell(stateFile)//']',verbosity=verbosityWorking)\n";
	    $postContains->[0]->{'content'} .= "     default => ".$directive->{'name'}."()\n";
	    $postContains->[0]->{'content'} .= "     call default%stateRestore(stateFile,fgslStateFile,stateOperationID)\n";
	    $postContains->[0]->{'content'} .= "    else\n";
	    $postContains->[0]->{'content'} .= "     call Galacticus_Display_Message(var_str('skipping default object of \""  .$directive->{'name'}."\" class [position: ')//FTell(stateFile)//']',verbosity=verbosityWorking)\n";
	    $postContains->[0]->{'content'} .= "    end if\n";
	    $postContains->[0]->{'content'} .= "    return\n";
	    $postContains->[0]->{'content'} .= "  end subroutine ".$directive->{'name'}."DoStateRetrieve\n\n";

	    # Create global calculation reset function.
	    if ( exists($directive->{'calculationReset'}) && $directive->{'calculationReset'} eq "yes" ) {
		&Galacticus::Build::SourceTree::SetVisibility($node->{'parent'},$directive->{'name'}."DoCalculationReset","public");
		$postContains->[0]->{'content'} .= "  !# <calculationResetTask>\n";
		$postContains->[0]->{'content'} .= "  !#  <unitName>".$directive->{'name'}."DoCalculationReset</unitName>\n";
		$postContains->[0]->{'content'} .= "  !# </calculationResetTask>\n";
		$postContains->[0]->{'content'} .= "  subroutine ".$directive->{'name'}."DoCalculationReset(node)\n";
		$postContains->[0]->{'content'} .= "    !% Reset calculations.\n";
		$postContains->[0]->{'content'} .= "    implicit none\n";
		$postContains->[0]->{'content'} .= "    type (treeNode), intent(inout) :: node\n";
		$postContains->[0]->{'content'} .= "    class(".$directive->{'name'}."Class), pointer :: default\n\n";
		$postContains->[0]->{'content'} .= "    if (associated(".$directive->{'name'}."Default)) then\n";
		$postContains->[0]->{'content'} .= "      default => ".$directive->{'name'}."()\n";
		$postContains->[0]->{'content'} .= "      call default%calculationReset(node)\n";
		$postContains->[0]->{'content'} .= "    end if\n";
		$postContains->[0]->{'content'} .= "    return\n";
		$postContains->[0]->{'content'} .= "  end subroutine ".$directive->{'name'}."DoCalculationReset\n\n";
            }

	    # Create global destroy function.
	    &Galacticus::Build::SourceTree::SetVisibility($node->{'parent'},$directive->{'name'}."DoDestroy","public");
	    $postContains->[0]->{'content'} .= "  !# <functionClassDestroyTask>\n";
	    $postContains->[0]->{'content'} .= "  !#  <unitName>".$directive->{'name'}."DoDestroy</unitName>\n";
	    $postContains->[0]->{'content'} .= "  !# </functionClassDestroyTask>\n";
	    $postContains->[0]->{'content'} .= "  subroutine ".$directive->{'name'}."DoDestroy()\n";
	    $postContains->[0]->{'content'} .= "    !% Destroy the default object.\n";
	    $postContains->[0]->{'content'} .= "    implicit none\n";
	    $postContains->[0]->{'content'} .= "    !# <objectDestructor owner=\"module:".$node->{'parent'}->{'name'}."\" name=\"".$directive->{'name'}."Default\"/>\n";
	    $postContains->[0]->{'content'} .= "    return\n";
	    $postContains->[0]->{'content'} .= "  end subroutine ".$directive->{'name'}."DoDestroy\n\n";

	    # Create functions.
	    foreach my $methodName ( keys(%methods) ) {
		my $method = $methods{$methodName};
		# Insert arguments.
		my @arguments;
		if ( exists($method->{'argument'}) ) {
		    if ( UNIVERSAL::isa($method->{'argument'},"ARRAY") ) {
			push(@arguments,@{$method->{'argument'}});
		    } else {
			push(@arguments,  $method->{'argument'} );
		    }
		}
		my $pass = "yes";
		$pass = $method->{'pass'}
		    if ( exists($method->{'pass'}) );
		my $argumentList = "";
		my $argumentCode;
		if ( $pass eq "yes" ) {
		    $argumentCode .= "      class(".$directive->{'name'}."Class), intent(inout)";
		    $argumentCode .= ", target"
			if ( exists($method->{'selfTarget'}) && $method->{'selfTarget'} eq "yes" );
		    $argumentCode .= " :: self\n";
		}
		my $separator = "";
		my $unusedCode = "if (sizeof(self)<0.and.sizeof(self)>0) then\nend if\n";
		foreach my $argument ( @arguments ) {
		    (my $variables = $argument) =~ s/^.*::\s*(.*?)\s*$/$1/;
		    $argumentList .= $separator.$variables;
		    $argumentCode .= "      ".$argument."\n";
		    $separator     = ",";
		    my $declaration = &Fortran::Utils::Unformat_Variables($argument);
		    foreach ( @{$declaration->{'variables'}} ) {
			my $function = $declaration->{'intrinsic'} eq "procedure" ? "loc" : "sizeof";
			$unusedCode .= "if (".$function."(".$_.")<0.and.".$function."(".$_.")>0) then\nend if\n";
		    }
		}
		my $type;
		my $category;
		my $self;
		my $extension = "__";	       
		$extension = ""
		    if ( exists($method->{'code'}) );
		my $recursive = exists($method->{'recursive'}) && $method->{'recursive'} eq "yes" ? "recursive " : "";
		if ( $method->{'type'} eq "void" ) {
		    $category = "subroutine";
		    $type     = "";
		    $self     = "";
		} elsif ( $method->{'type'} =~ m/^class/ ) {
		    $category = "function";
		    $type     = "";
		    $self     = "      ".$method->{'type'}.", pointer :: ".$directive->{'name'}.ucfirst($methodName).$extension."\n";
		} elsif ( $method->{'type'} =~ m/^type/ || $method->{'type'} =~ m/,/ ) {
		    $category = "function";
		    $type     = "";
		    $self     = "      ".$method->{'type'}." :: ".$directive->{'name'}.ucfirst($methodName).$extension."\n";
		} else {
		    $category = "function";
		    $type     = $method->{'type'}." ";
		    $self     = "";
		}
		$postContains->[0]->{'content'} .= "   ".$recursive.$type.$category." ".$directive->{'name'}.ucfirst($methodName).$extension."(self";
		$postContains->[0]->{'content'} .= ",".$argumentList
		    unless ( $argumentList eq "" );
		$postContains->[0]->{'content'} .= ")\n";
		$postContains->[0]->{'content'} .= "      !% ".$method->{'description'}."\n";
		if ( exists($method->{'code'}) ) {
		    if ( exists($method->{'modules'}) ) {
			if ( reftype($method->{'modules'}) ) {
			    # Array of modules, with possible "only" clauses.
			    foreach my $module ( @{$method->{'modules'}} ) {
				$postContains->[0]->{'content'} .= "      use ".$module->{'name'}.(exists($module->{'only'}) ? ", only : ".join(",",@{$module->{'only'}}) : "")."\n";
			    }
			} else {
			    # Simple space-separated list of modules.
			    foreach ( split(/\s+/,$method->{'modules'}) ) {			    
				$postContains->[0]->{'content'} .= "      use".($_ eq "ISO_C_Binding" ? ", intrinsic :: " : "")." ".$_."\n";
			    }
			}
		    }
		} else {
		    $postContains->[0]->{'content'} .= "      use Galacticus_Error\n";
		}
		$postContains->[0]->{'content'} .= "      implicit none\n";
		$postContains->[0]->{'content'} .= $argumentCode;
		$postContains->[0]->{'content'} .= $self;
		if ( exists($method->{'code'}) ) {
		    my $code = "      ".$method->{'code'};
		    $code =~ s/\n/\n      /g;
		    $postContains->[0]->{'content'} .= $code."\n";
		} else {
		    $postContains->[0]->{'content'} .= "      call Galacticus_Error_Report('this is a null method - initialize the ".$directive->{'name'}." object before use'//".&Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($node,$node->{'line'}).")\n";
		    if ( $category eq "function" ) {
			# Avoid warnings about unset function values.
			$postContains->[0]->{'content'} .= "      ".$directive->{'name'}.ucfirst($methodName).$extension."=";
			my $setValue;
			if ( $method->{'type'} =~ m/^class/ ) {
			    $setValue = "> null()";
			} elsif ( $method->{'type'} =~ m/^type\s*\(\s*(.*)\s*\)/ ) {
			    if ( $method->{'type'} =~ m/,\s*pointer/ ) {
				$setValue = ">null()";
			    } else {
				$setValue = $1."()";
			    }
			} elsif ( $method->{'type'} =~ m/^integer/ ) {
			    $setValue = "0";
			} elsif ( $method->{'type'} =~ m/^double\s+precision/ ) {
			    $setValue = "0.0d0";
			} elsif ( $method->{'type'} =~ m/^logical/ ) {
			    $setValue = ".false.";
			}
			die("Process_FunctionClass(): do not know how to set '".$method->{'type'}."'")
			    unless ( defined($setValue) );
			$postContains->[0]->{'content'} .= $setValue."\n";
		    }
		    $postContains->[0]->{'content'} .= "      return\n";
		    # <workaround type="gfortran" PR="41209" url="https://gcc.gnu.org/bugzilla/show_bug.cgi?id=41209"/>
		    $postContains->[0]->{'content'} .= $unusedCode;
		}
		$postContains->[0]->{'content'} .= "   end ".$category." ".$directive->{'name'}.ucfirst($methodName).$extension."\n\n";
	    }
 	    
	    # Generate documentation.
	    my $documentation = "\\subsubsection{".$directive->{'descriptiveName'}."}\\label{sec:methods".ucfirst($directive->{'name'})."}\n\n";
	    $documentation   .= "Additional implementations for ".lc($directive->{'descriptiveName'})." are added using the {\\normalfont \\ttfamily ".$directive->{'name'}."} class.\n";
	    $documentation   .= "The implementation should be placed in a file containing the directive:\n";
	    $documentation   .= "\\begin{verbatim}\n";
	    $documentation   .= "!# <".$directive->{'name'}." name=\"".$directive->{'name'}."MyImplementation\">\n";
	    $documentation   .= "!# <description>A short description of the implementation.</description>\n";
	    $documentation   .= "!# </".$directive->{'name'}.">\n";
	    $documentation   .= "\\end{verbatim}\n";
	    $documentation   .= "where {\\normalfont \\ttfamily MyImplementation} is an appropriate name for the implemention. This file should be treated as a regular Fortran module, but without the initial {\\normalfont \\ttfamily module} and final {\\normalfont \\ttfamily end module} lines. That is, it may contain {\\normalfont \\ttfamily use} statements and variable declarations prior to the {\\normalfont \\ttfamily contains} line, and should contain all functions required by the implementation after that line. Function names should begin with {\\normalfont \\ttfamily ".&LaTeX_Breakable($directive->{'name'}."MyImplementation")."}. The file \\emph{must} define a type that extends the {\\normalfont \\ttfamily ".$directive->{'name'}."Class} class (or extends another type which is itself an extension of the {\\normalfont \\ttfamily ".$directive->{'name'}."Class} class), containing any data needed by the implementation along with type-bound functions required by the implementation. The following type-bound functions are required (unless inherited from the parent type):\n";
	    $documentation   .= "\\begin{description}\n";
	    # Create functions.
	    foreach my $methodName ( keys(%methods) ) {
		my $method = $methods{$methodName};
		$documentation   .= "\\item[{\\normalfont \\ttfamily ".$methodName."}] ".$method->{'description'};
		if ( exists($method->{'code'}) ) {
		    $documentation .= " A default implementation exists. If overridden the following interface must be used:\n";
		} else {
		    $documentation .= " Must have the following interface:\n";
		}
		$documentation   .= "\\begin{lstlisting}[language=Fortran,basicstyle=\\small\\ttfamily,escapechar=@,breaklines,prebreak=\\&,postbreak=\\&\\space\\space,columns=flexible,keepspaces=true,breakautoindent=true,breakindent=10pt]\n";
		# Insert arguments.
		my @arguments;
		if ( exists($method->{'argument'}) ) {
		    if ( UNIVERSAL::isa($method->{'argument'},"ARRAY") ) {
			push(@arguments,@{$method->{'argument'}});
		    } else {
			push(@arguments,  $method->{'argument'} );
		    }
		}
		unshift(@arguments,"class(".$directive->{'name'}."Class), intent(inout) :: self");
		my $argumentList = "";
		my $separator    = "";
		my @argumentDefinitions;
		foreach my $argument ( @arguments ) {
		    if ( $argument =~ $Fortran::Utils::variableDeclarationRegEx ) {
			my $intrinsic     = $1;
			my $type          = $2;
			my $attributeList = $3;
			my $variableList  = $4;
			my @variables  = &Fortran::Utils::Extract_Variables($variableList,keepQualifiers => 1,lowerCase => 0);
			my $declaration =
			{
			    intrinsic  => $intrinsic,
			    attributes => $attributeList,
			    variables  => \@variables
			}; 
			if ( defined($type) ) {
			    $type =~ s/\((.*)\)/$1/;
			    $declaration->{'type'} = $type;
			}
			if ( defined($attributeList) ) {
			    $attributeList =~ s/^\s*,\s*//;
			    my @attributes = &Fortran::Utils::Extract_Variables($attributeList,keepQualifiers => 1);
			    $declaration->{'attributes'} = \@attributes;
			}
			push(@argumentDefinitions,$declaration);
		    } else {
			print "Argument does not match expected pattern:\n\t".$argument."\n";
			die("Functions_Generate_Output: argument parse error");
		    }
		    (my $variables = $argument) =~ s/^.*::\s*(.*?)\s*$/$1/;
		    $argumentList .= $separator.$variables;
		    $separator     = ",";
		}
		my $type;
		my $category;
		if ( $method->{'type'} eq "void" ) {
		    $category = "subroutine";
		    $type     = "";
		} else {
		    $category = "function";
		    $type     = $method->{'type'}." ";
		}
		$documentation .= "   ".$type.$category." myImplementation".ucfirst($methodName)."(";
		$documentation .= $argumentList
		    unless ( $argumentList eq "" );
		$documentation .= ")\n";
		$documentation .= &Fortran::Utils::Format_Variable_Definitions(\@argumentDefinitions);
		$documentation .= "   end ".$type.$category." myImplementation".ucfirst($methodName)."\n";
		$documentation .= "\\end{lstlisting}\n\n";
	    }
	    $documentation   .= "\\end{description}\n\n";	    
	    $documentation   .= "Existing implementations are:\n";
	    $documentation   .= "\\begin{description}\n";
	    foreach my $className ( keys(%classes) ) {
		my $class = $classes{$className};
		$documentation   .= "\\item[{\\normalfont \\ttfamily ".$class->{'name'}."}] ".$class->{'description'};
		$documentation   .= " \\iflabelexists{phys:".$directive->{'name'}.":".$class->{'name'}."}{See \\S\\ref{phys:".$directive->{'name'}.":".$class->{'name'}."}.}{}\n";
	    }
	    $documentation   .= "\\end{description}\n\n";
	    system("mkdir -p doc/methods");
	    open(my $docHndl,">doc/methods/".$directive->{'name'}.".tex");
	    print $docHndl $documentation;
	    close($docHndl);
	    # Insert into tree.	
	    # <workaround type="gfortran" PR="41209" url="https://gcc.gnu.org/bugzilla/show_bug.cgi?id=41209"/>		    
	    # To allow processing of "GCC attributes unused" directives by our preprocessor (since GCC does not support them yet),
	    # we parse and process our generated code here, before serializing it back into the original node. Should we need to
	    # retain this behavior permanently it would be cleaner to just generate the code as text (i.e. not in a node), then
	    # parse into a tree and unshift() it to the start of the postcontains array.	    
	    my $treeTmp = &Galacticus::Build::SourceTree::ParseCode ($postContains->[0]->{'content'},'Galacticus::Build::SourceTree::Process::FunctionClass::Process_FunctionClass()');
	    &Galacticus::Build::SourceTree::ProcessTree($treeTmp);
	    $postContains->[0]->{'content'} = &Galacticus::Build::SourceTree::Serialize($treeTmp);
	    # </workaround>
	    &Galacticus::Build::SourceTree::InsertAfterNode   ($node            ,$preContains );
	    &Galacticus::Build::SourceTree::InsertPostContains($node->{'parent'},$postContains);
	}
	$node = &Galacticus::Build::SourceTree::Walk_Tree($node,\$depth);
    }
}

sub LaTeX_Breakable {
    my $text = shift;
    $text =~ s/([a-z])([A-Z])/$1\\-$2/g;
    return $text;
}

sub trimlc {
    (my $result = lc(shift())) =~ s/^\s+|\s+$//g;
    return $result;
}

sub striplc {
    (my $result = lc(shift())) =~ s/\s//g;
    return $result;
}

1;

#  LocalWords:  nonAbstractClass
