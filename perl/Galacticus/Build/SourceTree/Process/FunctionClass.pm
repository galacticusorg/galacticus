# Contains a Perl module which implements processing of functionClass directives.

package Galacticus::Build::SourceTree::Process::FunctionClass;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Data::Dumper;
use XML::Simple;
use Sort::Topo;
use LaTeX::Encode;
use Scalar::Util qw(reftype);
use List::Util;
use List::MoreUtils qw(first_index);
use List::ExtraUtils;
use List::Uniq ':all';
use File::Changes;
use Fortran::Utils;
use Text::Levenshtein;
use Text::Template 'fill_in_string';
use Storable qw(dclone);
use Galacticus::Build::SourceTree::Process::SourceIntrospection;
use Galacticus::Build::SourceTree::Process::FunctionClass::Utils;
use Galacticus::Build::SourceTree::Parse::Declarations;

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
    our $stateStorables;
    # Initialize deep copy actions database.
    our $deepCopyActions;
    # Determine if debugging output is required.
    our $debugging = exists($ENV{'GALACTICUS_OBJECTS_DEBUG'}) && $ENV{'GALACTICUS_OBJECTS_DEBUG'} eq "yes";
    # Get state storables database if we do not have it.
    $stateStorables = $xml->XMLin($ENV{'BUILDPATH'}."/stateStorables.xml")
	unless ( $stateStorables );
    # Walk the tree, looking for code blocks.
    my $node  = $tree;
    my $depth = 0;
    while ( $node ) {	
	if ( grep {$node->{'type'}."Class" eq $_} keys(%{$stateStorables->{'functionClasses'}}) ) {
	    $node->{'directive'}->{'processed'} = 1;
	}
	if ( $node->{'type'} eq "functionClass" ) {
	    # Assert that our parent is a module.
	    die("Process_FunctionClass: parent node must be a module")
		unless ( $node->{'parent'}->{'type'} eq "module" );
	    my $lineNumber = $node->{'line'};
	    # Extract the directive and mark as processed.
	    my $directive = $node->{'directive'};
	    $directive->{'processed'} = 1;
	    # Get code directive locations if we do not have them.
	    $directiveLocations = $xml->XMLin($ENV{'BUILDPATH'}."/directiveLocations.xml")
		unless ( $directiveLocations );
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
	    my %classCounts;
	    foreach my $classLocation ( @classLocations ) {
		my $classTree = &Galacticus::Build::SourceTree::ParseFile($classLocation);
		&Galacticus::Build::SourceTree::ProcessTree($classTree, errorTolerant => 1);
		my $classNode = $classTree;
		(my $class, my @classDependencies) = &Galacticus::Build::SourceTree::Process::FunctionClass::Utils::Class_Dependencies($classNode,$directive->{'name'});
		foreach my $classDependency ( @classDependencies ) {
		    push(@{$dependencies{$classDependency}},$class->{'type'})
			unless ( $classDependency eq $class->{'type'});
		}
		# Store tree and file location.
		$class->{'file'} = $classLocation;
		$class->{'tree'} = $classTree;
		# Set defaults.
		$class->{'abstract'} = "no"
		    unless ( exists($class->{'abstract'}) );
	        # Store to set of all classes.
		die('Galacticus::Build::SourceTree::Process::FunctionClass::Process_FunctionClass: class is undefined in file "'.$classLocation.'"')
		    unless ( defined($class->{'type'}) );
		$classes{$class->{'type'}} = $class;
		push(@{$classCounts{$class->{'type'}}},$class->{'file'});
	    }
	    # Check for duplicated class names.
	    my $duplicatesFound = 0;
	    foreach my $className ( sort(keys(%classCounts)) ) {
		next
		    unless ( scalar(@{$classCounts{$className}}) > 1 );
		$duplicatesFound = 1;
		print "Duplicate class '".$className."' found in:\n".join("\n",map {"\t".$_} @{$classCounts{$className}})."\n";
	    }
	    die("ERROR: Galacticus::Build::SourceTree::Process::FunctionClass::Process_FunctionClass(): duplicate classes found")
		if ( $duplicatesFound );
	    # Sort classes. We first impose an alphanumeric sort to ensure that the ordering is always identical for each build,
	    # and then impose a topological sort to ensure that dependencies are correctly handled.
	    my @unsortedClasses = sort(keys(%classes));
	    my @sortedClasses   = &Sort::Topo::sort(\@unsortedClasses,\%dependencies);
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
	    # Construct short names for each non-abstract class.
	    foreach my $nonAbstractClass ( @nonAbstractClasses ) {
		if ( $nonAbstractClass->{'name'} =~ m/^$directive->{'name'}([a-zA-Z0-9]+)/ ) {
		    $nonAbstractClass->{'shortName'} = $1;
		    $nonAbstractClass->{'shortName'} = lcfirst($nonAbstractClass->{'shortName'})
			unless ( $nonAbstractClass->{'shortName'} =~ m/^[A-Z]{2,}/ );
		} else {
		    die("class name has incorrect format");
		}
	    }
	    # Validate any default class.
	    if ( exists($directive->{'default'}) ) {
		unless ( grep {$directive->{'default'} eq $_->{'shortName'}} @nonAbstractClasses ) {
		    print "ERROR: unrecognized default '".$directive->{'default'}."' for class '".$directive->{'name'}."'\n";
		    print "  allowed defaults are:\n";
		    print join("\n",map {"    ".$_->{'shortName'}} @nonAbstractClasses)."\n";
		    exit 1;
		}
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
			 name      => "ISO_C_Binding",
			 intrinsic => 1              ,
			 only      => [ "c_ptr" ]
		     }
		    ],
		argument    => [ "integer, intent(in   ) :: stateFile", "type(c_ptr), intent(in   ) :: gslStateFile" ],
	    };
	    $methods{'stateRestore'} =
	    {
		description => "Restore the state of the object to file.",
		type        => "void",
		pass        => "yes",
		modules     =>
		    [
		     {
			 name      => "ISO_C_Binding",
			 intrinsic => 1              ,
			 only      => [ "c_ptr" ]
		     }
		    ],
		argument    => [ "integer, intent(in   ) :: stateFile", "type(c_ptr), intent(in   ) :: gslStateFile" ],
	    };
	    # Add auto-hook function.
	    $methods{'autoHook'} =
	    {
		description => "Insert any event hooks required by this object.",
		type        => "void",
		pass        => "yes",
		code        => "!\$GLC attributes unused :: self\n\n! Nothing to do by default.\n"
	    };
	    if ( exists($directive->{'autoHook'}) ) {
		foreach my $module ( &List::ExtraUtils::as_array($directive->{'autoHook'}->{'modules'}) ) {
		    my $moduleName = $module->{'name'};
		    my @only       = split(/\s*,\s*/,$module->{'only'});
		    push(@{$methods{'autoHook'}->{'modules'}},{name => $moduleName, only => \@only});
		}
		$methods{'autoHook'}->{'code'} = $directive->{'autoHook'}->{'code'};
	    }
	    # Add destructor function.
	    if ( exists($directive->{'destructor'}) ) {
		$methods{'destructor'} =
		{
		    description => "Destructor for this class.",
		    type        => "void",
		    pass        => "yes",
		    code        => $directive->{'destructor'}->{'code'}
		};
		foreach my $module ( &List::ExtraUtils::as_array($directive->{'destructor'}->{'modules'}) ) {
		    my $moduleName = $module->{'name'};
		    my @only       = split(/\s*,\s*/,$module->{'only'});
		    push(@{$methods{'destructor'}->{'modules'}},{name => $moduleName, only => \@only});
		}
	    }
	    # Add "descriptor" method.
	    my $descriptorCode;
	    my %descriptorModules = ( "Input_Parameters" => 1 );
	    my %addSubParameters;
	    my $addLabel                  = 0;
	    my $rankMaximum               = 0;
	    my $descriptorUsed            = 0;
	    my $fileModificationCodeAdded = 0;
	    my $descriptorLinkedListVariables;
	    @{$descriptorLinkedListVariables} = ();
	    $descriptorCode .= "logical :: includeFileModificationTimes_\n";
	    $descriptorCode .= "if (present(includeFileModificationTimes)) then\n";
	    $descriptorCode .= " includeFileModificationTimes_=includeFileModificationTimes\n";
	    $descriptorCode .= "else\n";
	    $descriptorCode .= " includeFileModificationTimes_=.false.\n";
	    $descriptorCode .= "end if\n";
	    $descriptorCode .= "select type (self)\n";
	    foreach my $nonAbstractClass ( @nonAbstractClasses ) {
		(my $label = $nonAbstractClass->{'name'}) =~ s/^$directive->{'name'}//;
		$label = lcfirst($label)
		    unless ( $label =~ m/^[A-Z]{2,}/ );
		$nonAbstractClass->{'hasCustomDescriptor'} = 0;
		my $extensionOf;
		# Build lists of all potential parameter and object names for this class, including any from parent classes.
		my $potentialNames = {};
		my $class          = $nonAbstractClass;
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
			&potentialDescriptorParameters($node->{'declarations'},$nonAbstractClass,$potentialNames)
			    if ( $node->{'type'} eq "declaration" );
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
		    &potentialDescriptorParameters($declaration,$nonAbstractClass,$potentialNames);
		}
		# Add any names declared in the functionClassType.
		if ( defined($functionClassType) ) {
		    # Search the node for declarations.
		    my $node = $functionClassType->{'node'}->{'firstChild'};
		    while ( $node ) {
			&potentialDescriptorParameters($node->{'declarations'},$nonAbstractClass,$potentialNames)
			    if ( $node->{'type'} eq "declaration" );
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
		    if ( $node->{'type'} eq "function" && (grep {$_ eq $node->{'name'}} @constructors) && $node->{'opener'} =~ m/^\s*(recursive\s+)??function\s+$node->{'name'}\s*\(\s*parameters\s*(,\s*recursiveConstruct\s*,\s*recursiveSelf\s*)??\)/ ) {
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
					    if ( $constructorNode->{'directive'}->{'variable'} =~ m/(.*)\%(.*)/ ) {
						my $object  = $1;
						my $element = $2;
						if ( lc($object) eq lc($result) ) {
						    # Direct read into an element of the object being constructed.
						    $name = $element;
						} else {
						    # Read of some other derived-type component. Use the name of the derived type
						    # variable in case it is of a type that we can handle.
						    $name = $object;
						}
					    } else {
						$name = $constructorNode->{'directive'}->{'variable'};
					    }
					} else {
					    $name = $constructorNode->{'directive'}->{'name'};
					}
					if ( grep {lc($_) eq lc($name)} (map {@{$_->{'variableNames'}}} @{$potentialNames->{'parameters'}}) ) {
					    push(@{$descriptorParameters->{'parameters'}},{name => $name, inputName => $constructorNode->{'directive'}->{'name'}, source => $constructorNode->{'directive'}->{'source'}});
					    # Find the matched variable.
					    my $descriptor;
					    foreach my $potentialDescriptor ( @{$potentialNames->{'parameters'}} ) {
						$descriptor = $potentialDescriptor
						    if ( grep {lc($_) eq lc($name)} @{$potentialDescriptor->{'variableNames'}} );
					    }
					} elsif ( grep {lc($_) eq lc($name)."_"} (map {@{$_->{'variableNames'}}} @{$potentialNames->{'parameters'}}) ) {
					    push(@{$descriptorParameters->{'parameters'}},{name => $name."_", inputName => $constructorNode->{'directive'}->{'name'}, source => $constructorNode->{'directive'}->{'source'}});
					    # Find the matched variable.
					    my $descriptor;
					    foreach my $potentialDescriptor ( @{$potentialNames->{'parameters'}} ) {
						$descriptor = $potentialDescriptor
						    if ( grep {lc($_) eq lc($name)."_"} @{$potentialDescriptor->{'variableNames'}} );
					    }
					} elsif ( grep {lc($_) eq lc($name)} (map {@{$_->{'variableNames'}}} @{$potentialNames->{'enumerations'}}) ) {
					    push(@{$descriptorParameters->{'enumerations'}},{name => $name, inputName => $constructorNode->{'directive'}->{'name'}, source => $constructorNode->{'directive'}->{'source'}});
					    # Find the matched variable.
					    my $descriptor;
					    foreach my $potentialDescriptor ( @{$potentialNames->{'enumerations'}} ) {
						$descriptor = $potentialDescriptor
						    if ( grep {lc($_) eq lc($name)} @{$potentialDescriptor->{'variableNames'}} );
					    }
					} elsif ( grep {$_ eq lc($name)} (map {@{$_->{'variables'}}} @{$potentialNames->{'statefulTypes'}}) ) {
					    push(@{$descriptorParameters->{'statefulTypes'}},{name => $name, inputName => $constructorNode->{'directive'}->{'name'}, source => $constructorNode->{'directive'}->{'source'}});
					    # Find the matched variable.
					    my $descriptor;
					    foreach my $potentialDescriptor ( @{$potentialNames->{'statefulTypes'}} ) {
						$descriptor = $potentialDescriptor
						    if ( grep {$_ eq lc($name)} @{$potentialDescriptor->{'variables'}} );
					    }
					} else {
					    $supported = -1;
					    my $message = "could not find a matching internal variable for parameter [".$name."]";
					    my @potentialNames = map {@{$_->{'variableNames'}}} @{$potentialNames->{'parameters'}};
					    if ( scalar(@potentialNames) > 0 ) {
						my @distances      = &Text::Levenshtein::distance(lc($name),map {lc($_)} @potentialNames);
						my $indexMinimum   = first_index {$_ == &List::Util::min(@distances)} @distances;
						unless ( $indexMinimum == -1 ) {
						    (my $nameGuess = $potentialNames[$indexMinimum]) =~ s/_//;
						    $message .= " - did you mean [".$nameGuess."]";
						}
					    }
					    push(@failureMessage,$message);
					}
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
				    } elsif ( exists($nonAbstractClass->{'linkedList'}) && grep {$_ eq $name} split(" ",$nonAbstractClass->{'linkedList'}->{'object'}) ) {
					push(@{$descriptorParameters->{'linkedLists'}},$nonAbstractClass->{'linkedList'});
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
		foreach my $subParameterName ( sort(keys(%subParameters)) ) {
		    unless ( exists($subParameters{$subParameters{$subParameterName}->{'parent'}}) || $subParameters{$subParameterName}->{'parent'} eq "parameters" ) {
			$supported = -7;
			push(@failureMessage,"subparameter hierarchy failure");
		    }
		}
		# Build the code.
		$descriptorCode .= "type is (".$nonAbstractClass->{'name'}.")\n";
		if ( $nonAbstractClass->{'hasCustomDescriptor'} ) {
		    # The class has its own descriptor function, so we should never arrive at this point in the code.
		    $descriptorCode .= " call Error_Report('custom descriptor exists - this should not happen'//".&Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($nonAbstractClass->{'node'},$nonAbstractClass->{'node'}->{'line'}).")\n";
		    $descriptorModules{'Error'} = 1;
		} else{
		    # Build an auto-descriptor function.
		    if ( $declarationMatches && ( $supported == 1 || exists($nonAbstractClass->{'descriptorSpecial'}) ) ) {
			$descriptorUsed = 1;
			$descriptorCode .= " if (present(includeClass)) then\n";
			$descriptorCode .= "  includeClass_=includeClass\n";
			$descriptorCode .= " else\n";
			$descriptorCode .= "  includeClass_=.true.\n";
			$descriptorCode .= " end if\n";
			$descriptorCode .= " if (includeClass_) call descriptor%addParameter('".$directive->{'name'}."','".$label."')\n";
			if ( defined($descriptorParameters) ) {
			    # Get subparameters.
			    $addSubParameters{'parameters'} = 1;
			    $descriptorCode   .= "parameters=descriptor%subparameters('".$directive->{'name'}."')\n";
			    foreach my $subParameterName ( sort(keys(%subParameters)) ) {
				$addSubParameters{$subParameterName} = 1;
				$descriptorCode .= $subParameters{$subParameterName}->{'source'};
			    }
			    # Handle parameters set via inputParameter directives.
			    if ( defined($descriptorParameters->{'parameters'}) ) {
				foreach my $parameter ( @{$descriptorParameters->{'parameters'}} ) {
				    foreach my $declaration ( @{$potentialNames->{'parameters'}} ) {
					if ( grep {$_ eq lc($parameter->{'name'})} @{$declaration->{'variables'}} ) {
					    my $format;
					    my $function;
					    my $isLogical = 0;
					    if      ( $declaration->{'intrinsic'} eq "type"             ) {
						# Varying string type.
						$function  = "char";
					    } elsif ( $declaration->{'intrinsic'} eq "logical"          ) {
						# Logical.
						$addLabel  = 1;
						$isLogical = 1;
					    } elsif ( $declaration->{'intrinsic'} eq "double precision" ) {
						$addLabel  = 1;
						$format    = "e17.10";
					    } elsif ( $declaration->{'intrinsic'} eq "integer"          ) {
						$addLabel  = 1;
						$format    = "i17";
					    } elsif ( $declaration->{'intrinsic'} eq "character"        ) {
						$function  = "trim";
					    }
					    if ( grep {$_ =~ m/^dimension\s*\([a-z0-9_:,\s]+\)/} @{$declaration->{'attributes'}} ) {
						# Non-scalar parameter - values must be concatenated.
						my $dimensionDeclarator = join(",",map {/^dimension\s*\(([a-zA-Z0-9_,:\s]+)\)/} @{$declaration->{'attributes'}});
						my $rank                = ($dimensionDeclarator =~ tr/,//)+1;
						$rankMaximum            = $rank
						    if ( $rank > $rankMaximum );
						$descriptorCode .= "parameterValues=''\n";
						for(my $i=1;$i<=$rank;++$i) {
						    $descriptorCode .= " parameterValues=parameterValues//'['\n";
						    $descriptorCode .= "do i".$i."=lbound(self%".$parameter->{'name'}.",dim=".$i."),ubound(self%".$parameter->{'name'}.",dim=".$i.")\n";
						}
						if ( $function ) {
						    $descriptorCode .= " parameterValues=parameterValues//".$function."(self%".$parameter->{'name'}."(".join(",",map {"i".$_} 1..$rank)."))\n";
						} else {
						    if ( $isLogical ) {
							$descriptorCode .= "if (self%".$parameter->{'name'}."(".join(",",map {"i".$_} 1..$rank).") then\n";
							$descriptorCode .= "  parameterLabel='true'\n";
							$descriptorCode .= "else\n";
							$descriptorCode .= "  parameterLabel='false'\n";
							$descriptorCode .= "end if\n";
						    } else {
							$descriptorCode .= "write (parameterLabel,'(".$format.")') self%".$parameter->{'name'}."(".join(",",map {"i".$_} 1..$rank).")\n";
						    }
						    $descriptorCode .= " parameterValues=parameterValues//trim(adjustl(parameterLabel))\n";
						}
						$descriptorCode .= " if (i".$rank." /= size(self%".$parameter->{'name'}.",dim=".$rank.")) parameterValues=parameterValues//','\n";
						for(my $i=1;$i<=$rank;++$i) {
						    $descriptorCode .= "end do\n";
						    $descriptorCode .= " parameterValues=parameterValues//']'\n";
						    $descriptorCode .= " if (i".($i-1)." /= size(self%".$parameter->{'name'}.",dim=".($i-1).")) parameterValues=parameterValues//','\n"
							unless ( $i == 1 );
						}
						$descriptorCode .= "call ".$parameter->{'source'}."%addParameter('".$parameter->{'inputName'}."',char(parameterValues))\n";
					    } else {
						# Scalar parameter.
						if ( $function ) {
						    $descriptorCode .= "call ".$parameter->{'source'}."%addParameter('".$parameter->{'inputName'}."',".$function."(self%".$parameter->{'name'}."))\n";
						} else {
						    if ( $isLogical ) {
							$descriptorCode .= "if (self%".$parameter->{'name'}.") then\n";
							$descriptorCode .= "  parameterLabel='true'\n";
							$descriptorCode .= "else\n";
							$descriptorCode .= "  parameterLabel='false'\n";
							$descriptorCode .= "end if\n";
						    } else {
							$descriptorCode .= "write (parameterLabel,'(".$format.")') self%".$parameter->{'name'}."\n";
						    }
						    $descriptorCode .= "call ".$parameter->{'source'}."%addParameter('".$parameter->{'inputName'}."',trim(adjustl(parameterLabel)))\n";
						}
					    }
					}
				    }
				}
			    }
			    # Enumerations.
			    if ( defined($descriptorParameters->{'enumerations'}) ) {
				foreach my $parameter ( @{$descriptorParameters->{'enumerations'}} ) {
				    foreach my $declaration ( @{$potentialNames->{'enumerations'}} ) {
					if ( grep {$_ eq lc($parameter->{'name'})} @{$declaration->{'variables'}} ) {
					    my $format = "i17";
					    my $isLogical = 0;
					    $addLabel  = 1;
					    if ( grep {$_ =~ m/^dimension\s*\([a-z0-9_:,\s]+\)/} @{$declaration->{'attributes'}} ) {
						# Non-scalar parameter - values must be concatenated.
						my $dimensionDeclarator = join(",",map {/^dimension\s*\(([a-zA-Z0-9_,:\s]+)\)/} @{$declaration->{'attributes'}});
						my $rank                = ($dimensionDeclarator =~ tr/,//)+1;
						$rankMaximum            = $rank
						    if ( $rank > $rankMaximum );
						$descriptorCode .= "parameterValues=''\n";
						for(my $i=1;$i<=$rank;++$i) {
						    $descriptorCode .= " parameterValues=parameterValues//'['\n";
						    $descriptorCode .= "do i".$i."=lbound(self%".$parameter->{'name'}.",dim=".$i."),ubound(self%".$parameter->{'name'}.",dim=".$i.")\n";
						}
						$descriptorCode .= "write (parameterLabel,'(".$format.")') self%".$parameter->{'name'}."(".join(",",map {"i".$_} 1..$rank).")%ID\n";
						$descriptorCode .= " parameterValues=parameterValues//trim(adjustl(parameterLabel))\n";
						$descriptorCode .= " if (i".$rank." /= size(self%".$parameter->{'name'}.",dim=".$rank.")) parameterValues=parameterValues//','\n";
						for(my $i=1;$i<=$rank;++$i) {
						    $descriptorCode .= "end do\n";
						    $descriptorCode .= " parameterValues=parameterValues//']'\n";
						    $descriptorCode .= " if (i".($i-1)." /= size(self%".$parameter->{'name'}.",dim=".($i-1).")) parameterValues=parameterValues//','\n"
							unless ( $i == 1 );
						}
						$descriptorCode .= "call ".$parameter->{'source'}."%addParameter('".$parameter->{'inputName'}."',char(parameterValues))\n";
					    } else {
						# Scalar parameter.
						$descriptorCode .= "write (parameterLabel,'(".$format.")') self%".$parameter->{'name'}."%ID\n";
						$descriptorCode .= "call ".$parameter->{'source'}."%addParameter('".$parameter->{'inputName'}."',trim(adjustl(parameterLabel)))\n";
					    }
					}
				    }
				}
			    }
			    # Stateful types.
			    if ( defined($descriptorParameters->{'statefulTypes'}) ) {
				foreach my $parameter ( @{$descriptorParameters->{'statefulTypes'}} ) {
				    foreach my $declaration ( @{$potentialNames->{'statefulTypes'}} ) {
					if ( grep {$_ eq lc($parameter->{'name'})} @{$declaration->{'variables'}} ) {
					    my $format;
					    my $isLogical = 0;
					    if      ( $declaration->{'type'} eq "statefulInteger" ) {
						$format     = "i17";
					    } elsif ( $declaration->{'type'} eq "statefulDouble"  ) {
						$format     = "e17.10";
					    } elsif ( $declaration->{'type'} eq "statefulLogical" ) {
						$isLogical = 1;
					    } else {
						die("unknown stateful-type");
					    }
					    $addLabel  = 1;
					    if ( grep {$_ =~ m/^dimension\s*\([a-z0-9_:,\s]+\)/} @{$declaration->{'attributes'}} ) {
						# Non-scalar parameter - values must be concatenated.
						my $dimensionDeclarator = join(",",map {/^dimension\s*\(([a-zA-Z0-9_,:\s]+)\)/} @{$declaration->{'attributes'}});
						my $rank                = ($dimensionDeclarator =~ tr/,//)+1;
						$rankMaximum            = $rank
						    if ( $rank > $rankMaximum );
						$descriptorCode .= "parameterValues=''\n";
						for(my $i=1;$i<=$rank;++$i) {
						    $descriptorCode .= " parameterValues=parameterValues//'['\n";
						    $descriptorCode .= "do i".$i."=lbound(self%".$parameter->{'name'}.",dim=".$i."),ubound(self%".$parameter->{'name'}.",dim=".$i.")\n";
						}
						$descriptorCode .= "if (self%".$parameter->{'name'}."(".join(",",map {"i".$_} 1..$rank).")%isSet) then\n";
						if ( $isLogical ) {
						    $descriptorCode .= "if (self%".$parameter->{'name'}."(".join(",",map {"i".$_} 1..$rank)."%value) then\n";
						    $descriptorCode .= "  parameterLabel='true'\n";
						    $descriptorCode .= "else\n";
						    $descriptorCode .= "  parameterLabel='false'\n";
						    $descriptorCode .= "end if\n";
						} else {
						    $descriptorCode .= "write (parameterLabel,'(".$format.")') self%".$parameter->{'name'}."(".join(",",map {"i".$_} 1..$rank).")%value\n";
						}
						$descriptorCode .= " else\n";
						$descriptorCode .= "  parameterLabel='?'\n";
						$descriptorCode .= " end if\n";
						$descriptorCode .= " parameterValues=parameterValues//trim(adjustl(parameterLabel))\n";
						$descriptorCode .= " if (i".$rank." /= size(self%".$parameter->{'name'}.",dim=".$rank.")) parameterValues=parameterValues//','\n";
						for(my $i=1;$i<=$rank;++$i) {
						    $descriptorCode .= "end do\n";
						    $descriptorCode .= " parameterValues=parameterValues//']'\n";
						    $descriptorCode .= " if (i".($i-1)." /= size(self%".$parameter->{'name'}.",dim=".($i-1).")) parameterValues=parameterValues//','\n"
							unless ( $i == 1 );
						}
						$descriptorCode .= "call ".$parameter->{'source'}."%addParameter('".$parameter->{'inputName'}."',char(parameterValues))\n";
					    } else {
						# Scalar parameter.
						$descriptorCode .= "if (self%".$parameter->{'name'}."%isSet) then\n";
						if ( $isLogical ) {
						    $descriptorCode .= "if (self%".$parameter->{'name'}."%value) then\n";
						    $descriptorCode .= "  parameterLabel='true'\n";
						    $descriptorCode .= "else\n";
						    $descriptorCode .= "  parameterLabel='false'\n";
						    $descriptorCode .= "end if\n";
						} else {
						    $descriptorCode .= "write (parameterLabel,'(".$format.")') self%".$parameter->{'name'}."%value\n";
						}
						$descriptorCode .= " else\n";
						$descriptorCode .= "  parameterLabel='?'\n";
						$descriptorCode .= " end if\n";
						$descriptorCode .= "call ".$parameter->{'source'}."%addParameter('".$parameter->{'inputName'}."',trim(adjustl(parameterLabel)))\n";
					    }
					 }
				    }
				}
			    }
			    # Handle objects built via objectBuilder directives.
			    if ( defined($descriptorParameters->{'objects'}) ) {
				foreach ( @{$descriptorParameters->{'objects'}} ) {
                                    # Always include the class for composited objects - this ensures that the object is actually created.
				    $descriptorCode .= "if (associated(self%".$_->{'name'}.")) call self%".$_->{'name'}."%descriptor(parameters,includeClass=.true.,includeFileModificationTimes=includeFileModificationTimes)\n";
				}
			    }
			    # Handle linked lists.
			    if ( defined($descriptorParameters->{'linkedLists'}) ) {
				foreach ( @{$descriptorParameters->{'linkedLists'}} ) {
				    (my $linkedListCode, my $linkedListModule) = &autoDescriptorLinkedList($_,$descriptorLinkedListVariables);
				    $descriptorCode .= $linkedListCode;
				    $descriptorModules{$linkedListModule} = 1
					if ( $linkedListModule );
				}
			    }
			}
			# If the parent constructor was used, call its descriptor method.
			if ( $parentConstructorUsed ) {
			    $descriptorCode .= "call self%".$extensionOf."%descriptor(descriptor,includeClass=.false.,includeFileModificationTimes=includeFileModificationTimes)\n";
			}
		    } elsif ( ! $declarationMatches     && ! exists($nonAbstractClass->{'descriptorSpecial'}) ) {
			die("Automatic descriptor can not be built for class '".$nonAbstractClass->{'name'}."': parameter-based constructor not found");
		    } elsif (   $supported         != 1 && ! exists($nonAbstractClass->{'descriptorSpecial'}) ) {
			die("Automatic descriptor can not be built for class '".$nonAbstractClass->{'name'}."' because:\n   ".join("\n   ",@failureMessage));
		    }
		}
		# Add run-time file dependency modification times if needed.
		{
		    $code::type = $nonAbstractClass->{'name'};
		    my $class = $nonAbstractClass;
		    while ( $class ) {
			if ( exists($class->{'runTimeFileDependencies'}) ) {
			    unless ( $fileModificationCodeAdded ) {
				$descriptorCode                      = "integer :: status\ncharacter(len=30) :: timeModification\ninteger :: countRunTimeFileDependency\ntype(varying_string) :: fileDependencyParameterName\n".$descriptorCode;
				$descriptorModules{'File_Utilities' } = 1;
				$descriptorModules{'String_Handling'} = 1;
				$descriptorModules{'Error'          } = 1;
				$fileModificationCodeAdded            = 1;
			    }
			    $descriptorCode .= "if (includeFileModificationTimes_) then\ncountRunTimeFileDependency=0\n";
			    my @paths = split(" ",$class->{'runTimeFileDependencies'}->{'paths'});
			    foreach $code::path ( @paths ) {
				$code::introspection = &Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($nonAbstractClass->{'node'},$nonAbstractClass->{'node'}->{'line'});
				$descriptorCode .= fill_in_string(<<'CODE', PACKAGE => 'code');
timeModification=File_Modification_Time(self%{$path},status)
if (status == errorStatusSuccess) then
 countRunTimeFileDependency=countRunTimeFileDependency+1
 fileDependencyParameterName=var_str("runTimeFileDependency")//countRunTimeFileDependency
 call descriptor%addParameter(char(fileDependencyParameterName),char(self%{$path}//": "//trim(timeModification)))
else if (status /= errorStatusNotExist) then
 call Error_Report('unable to get file modification time'//{$introspection})
end if
CODE
			    }
			    $descriptorCode .= "end if\n";
			}
			$class = ($class->{'extends'} eq $directive->{'name'}) ? undef() : $classes{$class->{'extends'}};
		    }
		}
		# Call any special descriptor function.
		$descriptorCode .= " call self%".$nonAbstractClass->{'descriptorSpecial'}."(parameters)\n"
		    if ( exists($nonAbstractClass->{'descriptorSpecial'}) );
	    }
	    $descriptorCode .= "end select\n";
	    if ( scalar(@{$descriptorLinkedListVariables}) > 0 ) {
           	$descriptorCode = &Fortran::Utils::Format_Variable_Definitions($descriptorLinkedListVariables).$descriptorCode;
	    }
	    if ( $descriptorUsed ) {
		$descriptorCode = "logical :: includeClass_\n".$descriptorCode;
	    } else {
		$descriptorCode  = " !\$GLC attributes unused :: descriptor, includeClass\n".$descriptorCode;
	    }
 	    $descriptorCode  = "type(inputParameters) :: ".join(",",sort(keys(%addSubParameters)))."\n".$descriptorCode
		if ( %addSubParameters );
 	    $descriptorCode  = "character(len=18) :: parameterLabel\n".$descriptorCode
		if ( $addLabel );
	    if ( $rankMaximum > 0 ) {
		$descriptorCode  = "integer :: ".join(",",map {"i".$_} 1..$rankMaximum)."\ntype(varying_string) :: parameterValues\n".$descriptorCode
	    }
	    $methods{'descriptor'} =
	    {
		description => "Return an input parameter list descriptor which could be used to recreate this object.",
		type        => "void",
		pass        => "yes",
		modules     => join(" ",sort(keys(%descriptorModules))),
		argument    => [ "type(inputParameters), intent(inout) :: descriptor", "logical, intent(in   ), optional :: includeClass, includeFileModificationTimes" ],
		code        => $descriptorCode
	    };
	    # Add a "hashedDescriptor" method.
	    $code::directiveName = $directive->{'name'};
	    # <workaround type="gfortran" PR="102845" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=102845">
	    #  <description>
	    #   Nested parallelism results in memory leaks.
	    #  </description>
	    # </workaround>
	    my $hashedDescriptorCode = fill_in_string(<<'CODE', PACKAGE => 'code');
logical                        :: includeSourceDigest_
type   (inputParameters)       :: descriptor
type   (varying_string )       :: descriptorString
!   Workaround starts here.
! type   (varying_string ), save :: descriptorStringPrevious, hashedDescriptorPrevious
! !$omp threadprivate(descriptorStringPrevious,hashedDescriptorPrevious)
! Workaround ends here.
descriptor=inputParameters()
! Disable live nodeLists in FoX as updating these nodeLists leads to memory leaks.
call setLiveNodeLists(descriptor%document%document,.false.)
call self%descriptor(descriptor,includeClass=.true.,includeFileModificationTimes=includeFileModificationTimes)
descriptorString=descriptor%serializeToString()
call descriptor%destroy()
if (present(includeSourceDigest)) then
 includeSourceDigest_=includeSourceDigest
else
 includeSourceDigest_=.false.
end if
if (includeSourceDigest_) then
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
		$hashedDescriptorCode .= fill_in_string(<<'CODE', PACKAGE => 'code');
type is ({$type})
descriptorString=descriptorString//":sourceDigest\{"//String_C_To_Fortran({$type}5)//"\}"
CODE
	    }
	    # <workaround type="gfortran" PR="102845" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=102845">
	    #  <description>
	    #   Nested parallelism results in memory leaks.
	    #  </description>
	    # </workaround>
	    $hashedDescriptorCode .= fill_in_string(<<'CODE', PACKAGE => 'code');
end select
end if
!   Workaround starts here.
!   if (descriptorString /= descriptorStringPrevious) then
!      descriptorStringPrevious=         descriptorString
!      hashedDescriptorPrevious=Hash_MD5(descriptorString)
!   end if
!   {$directiveName}HashedDescriptor=hashedDescriptorPrevious
   {$directiveName}HashedDescriptor=Hash_MD5(descriptorString)
! Workaround ends here.
CODE
	    $methods{'hashedDescriptor'} =
	    {
		description => "Return a hash of the descriptor for this object, optionally include the source code digest in the hash.",
		type        => "type(varying_string)",
		pass        => "yes",
		modules     => "ISO_Varying_String String_Handling Input_Parameters Hashes_Cryptographic FoX_DOM",
		argument    => [ "logical, intent(in   ), optional :: includeSourceDigest, includeFileModificationTimes" ],
		code        => $hashedDescriptorCode
	    };
	    # Add a "objectType" method.
	    $code::directiveName = $directive->{'name'};
	    my $objectTypeCode = fill_in_string(<<'CODE', PACKAGE => 'code');
logical :: short_
short_=.false.
if (present(short)) short_=short
select type (self)
CODE
	    foreach my $nonAbstractClass ( @nonAbstractClasses ) {
		$code::type       = $nonAbstractClass->{'name'};
		($code::typeShort = $nonAbstractClass->{'name'}) =~ s/^$directive->{'name'}//;
		$code::typeShort  = lcfirst($code::typeShort);
		$objectTypeCode .= fill_in_string(<<'CODE', PACKAGE => 'code');
type is ({$type})
if (short_) then
 {$directiveName}ObjectType='{$typeShort}'
else
 {$directiveName}ObjectType='{$type}'
end if
CODE
	    }
	    $objectTypeCode .= fill_in_string(<<'CODE', PACKAGE => 'code');
end select
CODE
	    $methods{'objectType'} =
	    {
		description => "Return the type of the object.",
		type        => "type(varying_string)",
		pass        => "yes",
		modules     => "ISO_Varying_String",
		argument    => [ "logical, intent(in   ), optional :: short" ],
		code        => $objectTypeCode
	    };
	    # Add "allowedParameters" method.
	    my $allowedParametersCode;
	    my $allowedParameters;
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
		$allowedParameters->{$class->{'name'}}->{'declarationMatches'} = 0;
		while ( $node ) {
		    if ( $node->{'type'} eq "function" && (grep {$_ eq $node->{'name'}} @constructors) && $node->{'opener'} =~ m/^\s*(recursive)??\s+function\s+$node->{'name'}\s*\(\s*parameters\s*(\s*,\s*recursiveConstruct\s*,\s*recursiveSelf\s*)??\)/ ) {
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
				    $allowedParameters->{$class->{'name'}}->{'declarationMatches'} = 1
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
				# Look for calls to a parent class' parameter constructor.
				my $newContent;
				my $modified = 0;
				open(my $code,"<",\$constructorNode->{'content'});
				do {
				    # Get a line.
				    &Fortran::Utils::Get_Fortran_Line($code, my $rawLine, my $processedLine, my $bufferedComments);
				    if ( $processedLine =~ m/^\s*$result%([a-zA-Z0-9_]+)\s*=([a-zA-Z0-9_]+)\(\s*parameters\s*\)/ ) {
					$allowedParameters->{$class->{'name'}}->{'classParent'} = $1;
					$newContent .= $directive->{'name'}."DsblVldtn=.true.\n";
					$newContent .= $rawLine;
					$newContent .= $directive->{'name'}."DsblVldtn=.false.\n";
					$modified    = 1;
				    } else {
					$newContent .= $rawLine;
				    }
				} until ( eof($code) );
				close($code);
				$constructorNode->{'content'} = $newContent
				    if ( $modified );
			    }			    
			    if ( $constructorNode->{'type'} eq "inputParameter" ) {
				my $source = $constructorNode->{'directive'}->{'source'};
				if ( exists($constructorNode->{'directive'}->{'name'}) ) {
				    # A regular parameter, defined by its name.
				    push(@{$allowedParameters->{$class->{'name'}}->{'parameters'}->{$source}->{'all'}},$constructorNode->{'directive'}->{'name' });
				}
			    }
			    if ( $constructorNode->{'type'} eq "objectBuilder"  ) {
				my $source = $constructorNode->{'directive'}->{'source'};
				push(@{$allowedParameters->{$class->{'name'}}->{'parameters'}->{$source}->{'all'    }},exists($constructorNode->{'directive'}->{'parameterName'}) ? $constructorNode->{'directive'}->{'parameterName'} : $constructorNode->{'directive'}->{'class'});
				push(@{$allowedParameters->{$class->{'name'}}->{'parameters'}->{$source}->{'classes'}},exists($constructorNode->{'directive'}->{'parameterName'}) ? $constructorNode->{'directive'}->{'parameterName'} : $constructorNode->{'directive'}->{'class'});
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
							     @{$allowedParameters->{$class->{'name'}}->{'parameters'}->{$source}->{'objects'}},
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
	    }
	    my $allowedParametersLinkedListVariables;
	    @{$allowedParametersLinkedListVariables} = ();
	    $allowedParametersCode .= "select type (self)\n";
	    my %allowedParametersModules;
	    $allowedParametersModules{'ISO_Varying_String'} = 1;
	    foreach my $class ( @classes ) {
		if ( $allowedParameters->{$class->{'name'}}->{'declarationMatches'} ) {
		    my $className = $class->{'name'};
		    $allowedParametersCode .= "type is (".$className.")\n";
		    # Include the class and all parent classes for which the parent class parameter constructor is called.
		    while ( defined($className) ) {
			foreach my $source ( sort(keys(%{$allowedParameters->{$className}->{'parameters'}})) ) {
			    $allowedParametersCode .= "  if (objectsOnly) then\n";
			    {
				my $parameterCount = exists($allowedParameters->{$className}->{'parameters'}->{$source}->{'classes'}) ? scalar(@{$allowedParameters->{$className}->{'parameters'}->{$source}->{'classes'}}) : 0;
				if ( $parameterCount > 0 ) {
				    $parametersPresent          = 1;
				    $allowedParametersCode     .= "   if (sourceName == '".$source."') then\n";
				    $allowedParametersCode     .= "     countNew=0\n";
				    $allowedParametersCode     .= "     if (allocated(allowedParameters)) then\n";
				    for(my $i=0;$i<$parameterCount;++$i) {
					# <workaround type="gfortran" PR="37336" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=37336">
					#   <description>
					#     Array constructors are not correctly finalized. So, avoid using them.
					#   </description>
					# $allowedParametersCode .= "       if (.not.any(trim(allowedParameters) == '".$allowedParameters->{$className}->{'parameters'}->{$source}->{'classes'}->[$i]."')) countNew=countNew+1\n";
					$allowedParametersCode .= "       isNew=.true.\n";
					$allowedParametersCode .= "       do j=1,size(allowedParameters)\n";
					$allowedParametersCode .= "          if (allowedParameters(j) == '".$allowedParameters->{$className}->{'parameters'}->{$source}->{'classes'}->[$i]."') then\n";
					$allowedParametersCode .= "             isNew=.false.\n";
					$allowedParametersCode .= "             exit\n";
					$allowedParametersCode .= "          end if\n";
					$allowedParametersCode .= "       end do\n";
					$allowedParametersCode .= "       if (isNew) countNew=countNew+1\n";
					# </workaround>
				    }
				    $allowedParametersCode     .= "       if (countNew > 0) then\n";
				    $allowedParametersCode     .= "         call move_alloc(allowedParameters,allowedParametersTmp)\n";
				    $allowedParametersCode     .= "         allocate(allowedParameters(size(allowedParametersTmp)+countNew))\n";
				    $allowedParametersCode     .= "         allowedParameters(1:size(allowedParametersTmp))=allowedParametersTmp\n";
				    $allowedParametersCode     .= "         deallocate(allowedParametersTmp)\n";
				    for(my $i=0;$i<$parameterCount;++$i) {
					# <workaround type="gfortran" PR="37336" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=37336">
					#   <description>
					#     Array constructors are not correctly finalized. So, avoid using them.
					#   </description>
					# $allowedParametersCode .= "         if (.not.any(trim(allowedParameters(1:size(allowedParameters)-countNew)) == '".$allowedParameters->{$className}->{'parameters'}->{$source}->{'classes'}->[$i]."')) then\n";
					$allowedParametersCode .= "       isNew=.true.\n";
					$allowedParametersCode .= "       do j=1,size(allowedParameters)-countNew\n";
					$allowedParametersCode .= "          if (allowedParameters(j) == '".$allowedParameters->{$className}->{'parameters'}->{$source}->{'classes'}->[$i]."') then\n";
					$allowedParametersCode .= "             isNew=.false.\n";
					$allowedParametersCode .= "             exit\n";
					$allowedParametersCode .= "          end if\n";
					$allowedParametersCode .= "       end do\n";
					$allowedParametersCode .= "       if (isNew) then\n";
					# </workaround>
					$allowedParametersCode .= "           countNew=countNew-1\n";
					$allowedParametersCode .= "           allowedParameters(size(allowedParameters)-countNew)='".$allowedParameters->{$className}->{'parameters'}->{$source}->{'classes'}->[$i]."'\n";
					$allowedParametersCode .= "         end if\n";
				    }
				    $allowedParametersCode     .= "       end if\n";
				    $allowedParametersCode     .= "     else\n";
				    $allowedParametersCode     .= "       allocate(allowedParameters(".$parameterCount."))\n";
				    for(my $i=0;$i<$parameterCount;++$i) {
					$allowedParametersCode .= "       allowedParameters(".($i+1).")='".$allowedParameters->{$className}->{'parameters'}->{$source}->{'classes'}->[$i]."'\n";
				    }
				    $allowedParametersCode     .= "     end if\n";
				    $allowedParametersCode     .= "   end if\n";
				}
			    }
			    $allowedParametersCode .= "  else\n";
			    {
				my $parameterCount = scalar(@{$allowedParameters->{$className}->{'parameters'}->{$source}->{'all'}});
				if ( $parameterCount > 0 ) {
				    $parametersPresent          = 1;
				    $allowedParametersCode     .= "   if (sourceName == '".$source."') then\n";
				    $allowedParametersCode     .= "     countNew=0\n";
				    $allowedParametersCode     .= "     if (allocated(allowedParameters)) then\n";
				    for(my $i=0;$i<$parameterCount;++$i) {
					# <workaround type="gfortran" PR="37336" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=37336">
					#   <description>
					#     Array constructors are not correctly finalized. So, avoid using them.
					#   </description>
					# $allowedParametersCode .= "       if (.not.any(trim(allowedParameters) == '".$allowedParameters->{$className}->{'parameters'}->{$source}->{'all'}->[$i]."')) countNew=countNew+1\n";
					$allowedParametersCode .= "       isNew=.true.\n";
					$allowedParametersCode .= "       do j=1,size(allowedParameters)\n";
					$allowedParametersCode .= "          if (allowedParameters(j) == '".$allowedParameters->{$className}->{'parameters'}->{$source}->{'all'}->[$i]."') then\n";
					$allowedParametersCode .= "             isNew=.false.\n";
					$allowedParametersCode .= "             exit\n";
					$allowedParametersCode .= "          end if\n";
					$allowedParametersCode .= "       end do\n";
					$allowedParametersCode .= "       if (isNew) countNew=countNew+1\n";
				    }	
				    $allowedParametersCode     .= "       if (countNew > 0) then\n";
				    $allowedParametersCode     .= "         call move_alloc(allowedParameters,allowedParametersTmp)\n";
				    $allowedParametersCode     .= "         allocate(allowedParameters(size(allowedParametersTmp)+countNew))\n";
				    $allowedParametersCode     .= "         allowedParameters(1:size(allowedParametersTmp))=allowedParametersTmp\n";
				    $allowedParametersCode     .= "         deallocate(allowedParametersTmp)\n";
				    for(my $i=0;$i<$parameterCount;++$i) {
					# <workaround type="gfortran" PR="37336" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=37336">
					#   <description>
					#     Array constructors are not correctly finalized. So, avoid using them.
					#   </description>
					# $allowedParametersCode .= "         if (.not.any(trim(allowedParameters(1:size(allowedParameters)-countNew)) == '".$allowedParameters->{$className}->{'parameters'}->{$source}->{'all'}->[$i]."')) then\n";
					$allowedParametersCode .= "       isNew=.true.\n";
					$allowedParametersCode .= "       do j=1,size(allowedParameters)-countNew\n";
					$allowedParametersCode .= "          if (allowedParameters(j) == '".$allowedParameters->{$className}->{'parameters'}->{$source}->{'all'}->[$i]."') then\n";
					$allowedParametersCode .= "             isNew=.false.\n";
					$allowedParametersCode .= "             exit\n";
					$allowedParametersCode .= "          end if\n";
					$allowedParametersCode .= "       end do\n";
					$allowedParametersCode .= "       if (isNew) then\n";
					# </workaround>
					$allowedParametersCode .= "           countNew=countNew-1\n";
					$allowedParametersCode .= "           allowedParameters(size(allowedParameters)-countNew)='".$allowedParameters->{$className}->{'parameters'}->{$source}->{'all'}->[$i]."'\n";
					$allowedParametersCode .= "         end if\n";
				    }
				    $allowedParametersCode     .= "       end if\n";
				    $allowedParametersCode     .= "     else\n";
				    $allowedParametersCode     .= "       allocate(allowedParameters(".$parameterCount."))\n";
				    for(my $i=0;$i<$parameterCount;++$i) {
					$allowedParametersCode .= "       allowedParameters(".($i+1).")='".$allowedParameters->{$className}->{'parameters'}->{$source}->{'all'}->[$i]."'\n";
				    }
				    $allowedParametersCode     .= "     end if\n";
				    $allowedParametersCode     .= "   end if\n";
				}
			    }
			    $allowedParametersCode .= "  end if\n";
			    # Call the allowedParameters() method of any stored obejcts.
			    if ( $className eq $class->{'name'} ) {
				$parametersPresent      = 1
				    if ( exists($allowedParameters->{$className}->{'parameters'}->{$source}->{'objects'}) );
				foreach ( @{$allowedParameters->{$className}->{'parameters'}->{$source}->{'objects'}} ) {
				    $allowedParametersCode .= "  if (associated(self%".$_.")) call self%".$_."%allowedParameters(allowedParameters,'".$source."',.true.)\n";
				}
				# Handle any linked lists.
				(my $linkedListCode, my $linkedListModule) = &allowedParametersLinkedList($class,$allowedParametersLinkedListVariables,$source);
				$allowedParametersCode .= $linkedListCode;
				$allowedParametersModules{$linkedListModule} = 1
				    if ( $linkedListModule );
			    }
			}
			if ( defined($allowedParameters->{$className}->{'classParent'}) ) {
			    $className = $allowedParameters->{$className}->{'classParent'};
			} else {
			    undef($className);
			}
		    }
		}
	    }
	    $allowedParametersCode .= "end select\n";
	    if ( $parametersPresent ) {
		$allowedParametersCode = "type   (varying_string), allocatable, dimension(:) :: allowedParametersTmp\n".$directive->{'name'}."DsblVldtn=".$directive->{'name'}."DsblVldtn\n".$allowedParametersCode;
		$allowedParametersCode = "integer                                            :: countNew, j\n"                                                                              .$allowedParametersCode;
		$allowedParametersCode = "logical                                            :: isNew\n"                                                                                    .$allowedParametersCode;
	    } else {
		$allowedParametersCode = "!\$GLC attributes unused :: self, allowedParameters, sourceName\n".$directive->{'name'}."DsblVldtn=".$directive->{'name'}."DsblVldtn\n";
	    }
            $allowedParametersCode  = &Fortran::Utils::Format_Variable_Definitions($allowedParametersLinkedListVariables).$allowedParametersCode;
	    $methods{'allowedParameters'} =
	    {
		description => "Return a list of parameter names allowed for this object.",
		type        => "void",
		recursive   => "yes",
		pass        => "yes",
		modules     => join(" ",keys(%allowedParametersModules)),
		argument    => [
		    "type     (varying_string), dimension(:), allocatable, intent(inout) :: allowedParameters",
		    "character(len=*         )                           , intent(in   ) :: sourceName"       ,
		    "logical                                             , intent(in   ) :: objectsOnly"
		    ],
		code        => $allowedParametersCode
	    };
	    
	    # Add "assignment(=)" operator.
	    my $assignment;
            $assignment->{'code'        } .= "select type (self)\n";
	    foreach my $nonAbstractClass ( @nonAbstractClasses ) {
		# Add type guards.
		$assignment->{'code'} .= "type is (".$nonAbstractClass->{'name'}.")\n";
		$assignment->{'code'} .= "  select type (from)\n";
		$assignment->{'code'} .= "  type is (".$nonAbstractClass->{'name'}.")\n";
		# Search the tree for this class.
		my $class = $nonAbstractClass;
		while ( $class ) {
		    my $node = $class->{'tree'}->{'firstChild'};
		    $node = $node->{'sibling'}
		        while ( $node && ( $node->{'type'} ne "type" || ( ! exists($node->{'name'}) || $node->{'name'} ne $class->{'name'} ) ) );
		    last
			unless ( $node );
		    # Search the node for declarations.
		    $node = $node->{'firstChild'};
		    while ( $node ) {
			# Process declarations.
			if ( $node->{'type'} eq "declaration" ) {
			    foreach my $declaration ( &List::ExtraUtils::as_array($node->{'declarations'}) ) {
				my $isPointer   = grep {$_ eq "pointer"} @{$declaration->{'attributes'}};
				my $assigner    = $isPointer ? "=>" : "=";
				my $allocatable = grep {$_ eq "allocatable"} @{$declaration->{'attributes'}};
				(my $type = $declaration->{'type'}) =~ s/(^\s*|\s*$)//g
				    if ( $declaration->{'intrinsic'} eq "class" || $declaration->{'intrinsic'} eq "type" );
				my $referenceCount= 
				    ($declaration->{'intrinsic'} eq "class" || $declaration->{'intrinsic'} eq "type")
				    &&
				    (grep {$_ eq $type    } (keys(%{$stateStorables->{'functionClasses'}}),@{$stateStorables->{'functionClassInstances'}}))
				    &&
				    grep {$_ eq "pointer"}  @{$declaration   ->{'attributes'     }};
				foreach my $object ( @{$declaration->{'variables'}} ) {
				    (my $name = $object) =~ s/^([a-zA-Z0-9_]+).*/$1/; # Strip away anything (e.g. assignment operators) after the variable name.
				    if ( $allocatable ) {
					$assignment->{'code'} .= "    if (allocated(self%".$name.")) deallocate(self%".$name."                      )\n";
					$assignment->{'code'} .= "    if (allocated(from%".$name."))   allocate(self%".$name.",source=from%".$name.")\n";
				    } else {
					$assignment->{'code'} .= "    self%".$name.$assigner."from%".$name."\n";
					$assignment->{'code'} .= "    ".($isPointer ? "if (associated(self%".$name.")) " : "")."call self%".$name."%referenceCountIncrement()\n"
					    if ( $referenceCount );
				    }
				}
			    }
			}
			$node = $node->{'sibling'};
		    }
		    # Move to the parent class.
		    $class = ($class->{'extends'} eq $directive->{'name'}) ? undef() : $classes{$class->{'extends'}};
		}
		$assignment->{'code'} .= "  class default\n";
		$assignment->{'code'} .= "    call Error_Report('self and from types do not match'//".&Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($nonAbstractClass->{'node'},$nonAbstractClass->{'node'}->{'line'}).")\n";
		$assignment->{'code'} .= "  end select\n";
	    }
	    $assignment->{'code'} .= "end select\n";
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
		my $isPointer   = grep {$_ eq "pointer"} @{$declaration->{'attributes'}};
		my $assigner    = $isPointer ? "=>" : "=";
		my $allocatable = grep {$_ eq "allocatable"} @{$declaration->{'attributes'}};
		(my $type = $declaration->{'type'}) =~ s/(^\s*|\s*$)//g
							  if ( $declaration->{'intrinsic'} eq "class" || $declaration->{'intrinsic'} eq "type" );
		my $referenceCount= 
		    ($declaration->{'intrinsic'} eq "class" || $declaration->{'intrinsic'} eq "type")
		    &&
		    (grep {$_ eq $type    } (keys(%{$stateStorables->{'functionClasses'}}),@{$stateStorables->{'functionClassInstances'}}))
		    &&
		    grep {$_ eq "pointer"}  @{$declaration   ->{'attributes'     }};
		foreach my $object ( @{$declaration->{'variables'}} ) {
		    (my $name = $object) =~ s/^([a-zA-Z0-9_]+).*/$1/; # Strip away anything (e.g. assignment operators) after the variable name.
		    if ( $allocatable ) {
			$assignment->{'code'} .= "    if (allocated(self%".$name.")) deallocate(self%".$name."                      )\n";
			$assignment->{'code'} .= "    if (allocated(from%".$name."))   allocate(self%".$name.",source=from%".$name.")\n";
		    } else {
			$assignment->{'code'} .= "    self%".$name.$assigner."from%".$name."\n";
		    }
		    $assignment->{'code'} .= "    ".($isPointer ? "if (associated(self%".$name.")) " : "")."call self%".$name."%referenceCountIncrement()\n"
			if ( $referenceCount );
		}
	    }
	    # Add any objects declared in the functionClassType class.
	    if ( defined($functionClassType) ) {
		# Search the node for declarations.
		my @ignore = ();
		my $node   = $functionClassType->{'node'}->{'firstChild'};
		while ( $node ) {
		    if ( $node->{'type'} eq "declaration" ) {
			foreach my $declaration ( &List::ExtraUtils::as_array($node->{'declarations'}) ) {
			    my $isPointer   = grep {$_ eq "pointer"} @{$declaration->{'attributes'}};
			    my $assigner    = $isPointer ? "=>" : "=";
			    my $allocatable = grep {$_ eq "allocatable"} @{$declaration->{'attributes'}};
			    (my $type = $declaration->{'type'}) =~ s/(^\s*|\s*$)//g
								      if ( $declaration->{'intrinsic'} eq "class" || $declaration->{'intrinsic'} eq "type" );
			    my $referenceCount= 
				($declaration->{'intrinsic'} eq "class" || $declaration->{'intrinsic'} eq "type")
				&&
				(grep {$_ eq $type    } (keys(%{$stateStorables->{'functionClasses'}}),@{$stateStorables->{'functionClassInstances'}}))
				&&
				grep {$_ eq "pointer"}  @{$declaration   ->{'attributes'     }};
			    foreach my $object ( @{$declaration->{'variables'}} ) {
				(my $name = $object) =~ s/^([a-zA-Z0-9_]+).*/$1/; # Strip away anything (e.g. assignment operators) after the variable name.
				if ( $allocatable ) {
				    $assignment->{'code'} .= "    if (allocated(self%".$name.")) deallocate(self%".$name."                      )\n";
				    $assignment->{'code'} .= "    if (allocated(from%".$name."))   allocate(self%".$name.",source=from%".$name.")\n";
				} else {
				    $assignment->{'code'} .= "    self%".$name.$assigner."from%".$name."\n";
				}
				$assignment->{'code'} .= "    ".($isPointer ? "if (associated(self%".$name.")) " : "")."call self%".$name."%referenceCountIncrement()\n"
				    if ( $referenceCount );
			    }
			}			}
		    $node = $node->{'sibling'};
		}
	    }
	    # Add objects from the functionClass class.
	    $assignment->{'code'} .= "self%isDefaultOfClass=from%isDefaultOfClass\n";
	    $assignment->{'code'} .= "self%referenceCount=from%referenceCount\n";
	    $assignment->{'code'} .= "return\n";
	    $methods{'assignment(=)'} =
	    {
		description => "Assign the object.",
		type        => "void",
		recursive   => "yes",
		pass        => "yes",
		selfIntent  => "out",
		modules     => "Error",
		argument    => [ "class(".$directive->{'name'}."Class), intent(in   ) :: from" ],
		code        => $assignment->{'code'}
	    };
	    
	    # Add "deepCopy" method.
	    my $deepCopy;
            if ( $debugging ) {
		$deepCopy->{'modules'}->{'MPI_Utilities'     } = 1;
		$deepCopy->{'modules'}->{'ISO_Varying_String'} = 1;
		$deepCopy->{'modules'}->{'String_Handling'   } = 1;
		$deepCopy->{'modules'}->{'Display'           } = 1;
            }
	    $deepCopy->{'rankMaximum'} = 0;
            my $linkedListVariables;
            my $linkedListResetVariables;
            my $linkedListFinalizeVariables;
            @{$linkedListVariables        } = ();
            @{$linkedListResetVariables   } = ();
            @{$linkedListFinalizeVariables} = ();
            $deepCopy->{'resetCode'   } .= "self%copiedSelf => null()\n";
            $deepCopy->{'resetCode'   } .= "select type (self)\n";
	    $deepCopy->{'finalizeCode'} .= "self%copiedSelf => null()\n";
            $deepCopy->{'finalizeCode'} .= "select type (self)\n";
            $deepCopy->{'code'        } .= "select type (self)\n";
	    foreach my $nonAbstractClass ( @nonAbstractClasses ) {
		# Search the tree for this class.
		my $class = $nonAbstractClass;
		undef($deepCopy->{'assignments'});
		# Add a class guard for resets.
		$deepCopy->{'resetCode'   } .= "type is (".$nonAbstractClass->{'name'}.")\n";
		$deepCopy->{'finalizeCode'} .= "type is (".$nonAbstractClass->{'name'}.")\n";
		# Initialize a list of explicity-deep-copied variables that have been found.
		my $foundDeepCopyNames;
		@{$foundDeepCopyNames} = ();
		while ( $class ) {
		    my $node = $class->{'tree'}->{'firstChild'};
		    $node = $node->{'sibling'}
		        while ( $node && ( $node->{'type'} ne "type" || ( ! exists($node->{'name'}) || $node->{'name'} ne $class->{'name'} ) ) );
		    last
			unless ( $node );
		    # Handle linked lists.
		    (my $linkedListCode, my $linkedListResetCode, my $linkedListFinalizeCode, my $linkedListModule) = &deepCopyLinkedList($class,$nonAbstractClass,$linkedListVariables,$linkedListResetVariables,$linkedListFinalizeVariables,$debugging);
		    $deepCopy->{'assignments' } .= $linkedListCode;
		    $deepCopy->{'resetCode'   } .= $linkedListResetCode;
		    $deepCopy->{'finalizeCode'} .= $linkedListFinalizeCode;
		    if ( $linkedListModule ) {
			$deepCopy->{'modules'        }->{$linkedListModule} = 1;
			$deepCopy->{'resetModules'   }->{$linkedListModule} = 1;
			$deepCopy->{'finalizeModules'}->{$linkedListModule} = 1;
		    }
		    # Search the node for declarations.
		    my @ignore = exists($class->{'deepCopy'}->{'ignore'}) ? split(/\s*,\s*/,$class->{'deepCopy'}->{'ignore'}->{'variables'}) : ();
		    $node = $node->{'firstChild'};
		    while ( $node ) {
			&deepCopyDeclarations($class,$nonAbstractClass,$node,$node->{'declarations'},\@ignore,$lineNumber,$deepCopy,$foundDeepCopyNames)
			    if ( $node->{'type'} eq "declaration" );
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
		    my @ignore      = ();
		    &deepCopyDeclarations($class,$nonAbstractClass,$node,$declaration,\@ignore,$lineNumber,$deepCopy,$foundDeepCopyNames);
		}
		# Add any objects declared in the functionClassType class.
		if ( defined($functionClassType) ) {
		    # Search the node for declarations.
		    my @ignore = ();
		    my $node   = $functionClassType->{'node'}->{'firstChild'};
		    while ( $node ) {
			&deepCopyDeclarations($class,$nonAbstractClass,$node,$node->{'declarations'},\@ignore,$lineNumber,$deepCopy,$foundDeepCopyNames)
			    if ( $node->{'type'} eq "declaration" );
			$node = $node->{'sibling'};
		    }
		}
		# Check that the type of the destination matches, and perform the copy. Reset the reference count to the copy.
		$deepCopy->{'code'} .= "type is (".$nonAbstractClass->{'name'}.")\n";
		$deepCopy->{'code'} .= "select type (destination)\n";
		$deepCopy->{'code'} .= "type is (".$nonAbstractClass->{'name'}.")\n";
		$deepCopy->{'code'} .= "destination=self\n";
		$deepCopy->{'code'} .= $deepCopy->{'assignments'}
		    if ( defined($deepCopy->{'assignments'}) );
		$deepCopy->{'code'} .= "class default\n";
		$deepCopy->{'code'} .= "call Error_Report('destination and source types do not match'//".&Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($nonAbstractClass->{'node'},$nonAbstractClass->{'node'}->{'line'}).")\n";
		$deepCopy->{'code'} .= "end select\n";
		# Specify required modules.
		$deepCopy->{'modules'}->{'Error'} = 1;
		# Check that all explicit variables were found.
		{
		    my $class = $nonAbstractClass;
		    while ( $class ) {
			if ( exists($class->{'deepCopy'}->{'functionClass'}) ) {
			    foreach my $variable ( split(/\s*,\s*/,$class->{'deepCopy'}->{'functionClass'}->{'variables'}) ) {				
				die("Error: unable to find variable '".$variable."' marked for deep copy in class '".$class->{'name'}."'")
				    unless ( grep {$_ eq lc($variable)} @{$foundDeepCopyNames} );
			    }
			}
			$class = ($class->{'extends'} eq $directive->{'name'}) ? undef() : $classes{$class->{'extends'}};
		    }
		}
	    }
            $deepCopy->{'code'        } .= "end select\n";
	    $deepCopy->{'resetCode'   } .= "end select\n";
	    $deepCopy->{'finalizeCode'} .= "end select\n";
            # Reset the reference count to this newly created object.
            $deepCopy->{'code'} .= "call destination%referenceCountReset()\n";
            # Reset the state operation ID if necessary.
            $deepCopy->{'code'} .= "destination%stateOperationID=0_c_size_t\n";
            # Insert variables declarations.
            $deepCopy->{'code'        } = &Fortran::Utils::Format_Variable_Definitions($linkedListVariables        ).$deepCopy->{'code'}        ;
            $deepCopy->{'resetCode'   } = &Fortran::Utils::Format_Variable_Definitions($linkedListResetVariables   ).$deepCopy->{'resetCode'}   ;
            $deepCopy->{'finalizeCode'} = &Fortran::Utils::Format_Variable_Definitions($linkedListFinalizeVariables).$deepCopy->{'finalizeCode'};
            # Insert any iterator variables needed.
            $deepCopy->{'code'} = "integer :: ".join(",",map {"i".$_} 1..$deepCopy->{'rankMaximum'})."\n".$deepCopy->{'code'}
                if ( $deepCopy->{'rankMaximum'} > 0 );
	    $methods{'deepCopy'} =
	    {
		description => "Perform a deep copy of the object. This is a wrapper around the actual deep-copy code.",
		type        => "void",
		recursive   => "yes",
		pass        => "yes",
		selfTarget  => "yes",
		argument    => [ "class(".$directive->{'name'}."Class), intent(inout) :: destination" ],
		code        => "call self%deepCopy_(destination)"
	    };
	    $methods{'deepCopy_'} =
	    {
		description => "Perform a deep copy of the object.",
		type        => "void",
		recursive   => "yes",
		pass        => "yes",
		modules     => join(" ",sort(keys(%{$deepCopy->{'modules'}}))),
		argument    => [ "class(".$directive->{'name'}."Class), intent(inout) :: destination" ],
		code        => $deepCopy->{'code'}
	    };
	    $methods{'deepCopyReset'} =
	    {
		description => "Reset deep copy pointers in this object and any objects that it uses.",
		type        => "void",
		recursive   => "yes",
		pass        => "yes",
		code        => $deepCopy->{'resetCode'}
	    };
	    $methods{'deepCopyFinalize'} =
	    {
		description => "Finalize a deep copy in this object and any objects that it uses.",
		type        => "void",
		recursive   => "yes",
		pass        => "yes",
		code        => $deepCopy->{'finalizeCode'}
	    };
	    $methods{'deepCopyReset'   }->{'modules'} = join(" ",sort(keys(%{$deepCopy->{'resetModules'   }})))
		if ( scalar(keys(%{$deepCopy->{'resetModules'   }})) > 0 );
	    $methods{'deepCopyFinalize'}->{'modules'} = join(" ",sort(keys(%{$deepCopy->{'finalizeModules'}})))
		if ( scalar(keys(%{$deepCopy->{'finalizeModules'}})) > 0 );
	    # Add "stateStore" and "stateRestore" method.
	    my $stateStores =
	    {
		stateFileUsed          => 0,
		gslStateFileUsed       => 0,
		rankMaximum            => 0,
		allocatablesFound      => 0,
		explicitFunctionsFound => 0,
		dimensionalsFound      => 0,
		labelUsed              => 0		
	    };
	    my $stateStoreCode;
	    my $stateRestoreCode;
	    my $stateLinkedListVariables;
            @{$stateLinkedListVariables} = ();
	    %{$stateStores->{'stateStoreModules'  }} = ( "Display" => 1, "ISO_Varying_String" => 1, "String_Handling" => 1, "ISO_C_Binding" => 1 );
	    %{$stateStores->{'stateRestoreModules'}} = ( "Display" => 1, "ISO_Varying_String" => 1, "String_Handling" => 1, "ISO_C_Binding" => 1 );
	    my @outputUnusedVariables;
	    my @inputUnusedVariables;
	    $stateStoreCode   .= "position=FTell(stateFile)\n";
	    $stateRestoreCode .= "position=FTell(stateFile)\n";
	    $stateStoreCode   .= "call displayIndent(var_str('storing state for \""  .$directive->{'name'}."\" [position: ')//position//']',verbosity=verbosityLevelWorking)\n";
	    $stateRestoreCode .= "call displayIndent(var_str('restoring state for \"".$directive->{'name'}."\" [position: ')//position//']',verbosity=verbosityLevelWorking)\n";
	    $stateStoreCode   .= "select type (self)\n";
	    $stateRestoreCode .= "select type (self)\n";
	    foreach my $nonAbstractClass ( @nonAbstractClasses ) {
		# Build the code.
		my $stateStore = $stateStores->{$nonAbstractClass->{'name'}};
		$stateStoreCode   .= "type is (".$nonAbstractClass->{'name'}.")\n";
		$stateRestoreCode .= "type is (".$nonAbstractClass->{'name'}.")\n";
		if ( exists($nonAbstractClass->{'recursive'}) && $nonAbstractClass->{'recursive'} eq "yes" ) {
		    # Object allows recursion. If this object is a recursive copy, call the state store/restore functions on the actual copy.
		    $stateStoreCode   .= "if (self%isRecursive) then\n";
		    $stateStoreCode   .= " call displayUnindent('recursive copy - moving to actual',verbosity=verbosityLevelWorking)\n";
		    $stateStoreCode   .= " call self%recursiveSelf%stateStore  (stateFile,gslStateFile,stateOperationID)\n";
		    $stateStoreCode   .= " return\n";
		    $stateStoreCode   .= "end if\n";
		    $stateRestoreCode .= "if (self%isRecursive) then\n";
		    $stateRestoreCode .= " call displayUnindent('recursive copy - moving to actual',verbosity=verbosityLevelWorking)\n";
		    $stateRestoreCode .= " call self%recursiveSelf%stateRestore(stateFile,gslStateFile,stateOperationID)\n";
		    $stateRestoreCode .= " return\n";
		    $stateRestoreCode .= "end if\n";
		}
		$stateStoreCode   .= "if (self%stateOperationID == stateOperationID) then\n"; # If this object was already stored, don't do it again.
		$stateStoreCode   .= " call displayUnindent('skipping - already stored',verbosity=verbosityLevelWorking)\n";
		$stateStoreCode   .= " return\n";
		$stateStoreCode   .= "end if\n";
		$stateStoreCode   .= "self%stateOperationID=stateOperationID\n";
		$stateRestoreCode .= "if (self%stateOperationID == stateOperationID) then\n"; # If this object was already restored, don't do it again.
		$stateRestoreCode .= " call displayUnindent('skipping - already restored',verbosity=verbosityLevelWorking)\n";
		$stateRestoreCode .= " return\n";
		$stateRestoreCode .= "end if\n";
		$stateRestoreCode .= "self%stateOperationID=stateOperationID\n";
		$stateStoreCode   .= " call displayMessage('object type \"".$nonAbstractClass->{'name'}."\"',verbosity=verbosityLevelWorking)\n";
		$stateRestoreCode .= " call displayMessage('object type \"".$nonAbstractClass->{'name'}."\"',verbosity=verbosityLevelWorking)\n";
		(my $label = $nonAbstractClass->{'name'}) =~ s/^$directive->{'name'}//;
		$label = lcfirst($label)
		    unless ( $label =~ m/^[A-Z]{2,}/ );
		$stateStore->{'hasCustomStateStore'  } = 0;
		$stateStore->{'hasCustomStateRestore'} = 0;
		my $extensionOf;
		# Generate code to output all variables from this class (and any parent class).
		@{$stateStore->{'staticVariables'}} = ();
		my $explicitNamesFound;
		@{$explicitNamesFound} = ();
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
		    @{$stateStore->{'excludes'}} = exists($class->{'stateStorable'}->{'exclude'}->{'variables'}) ? split(/\s*,\s*/,$class->{'stateStorable'}->{'exclude'}->{'variables'}) : ();
		    # Search the node for declarations.
		    $node = $node->{'firstChild'};
		    while ( $node ) {
			&stateStoreVariables($stateStores,$stateStore,$class,$node->{'declarations'},$explicitNamesFound)
			    if ( $node->{'type'} eq "declaration" );
			$node = $node->{'type'} eq "contains" ? $node->{'firstChild'} : $node->{'sibling'};
		    }
		    # Handle linked lists.
		    (my $linkedListInputCode, my $linkedListOutputCode, my $linkedListModule) = &stateStoreLinkedList($class,$nonAbstractClass,$stateLinkedListVariables);
		    $stateStore->{'inputCode' } .= $linkedListInputCode;
		    $stateStore->{'outputCode'} .= $linkedListOutputCode;
		    $stateStores->{'stateStoreModules'}->{$linkedListModule} = 1
			if ( $linkedListModule );
		    # Handle explicit state store functions.
		    $stateStores->{'explicitFunctionsFound'} = 1
			if ( exists($nonAbstractClass->{'stateStore'}->{'stateStore'}->{'restore'}) );
		    (my $stateStoreExplicitInputCode, my $stateStoreExplicitOutputCode, my %stateStoreExplicitModules) = &stateStoreExplicitFunction($nonAbstractClass);
		    $stateStore->{'inputCode'}  .= $stateStoreExplicitInputCode;
		    $stateStore->{'outputCode'} .= $stateStoreExplicitOutputCode;
		    foreach my $module ( sort(keys(%stateStoreExplicitModules)) ) {
			$stateStores->{'stateStoreModules'  }->{$module} = 1;
			$stateStores->{'stateRestoreModules'}->{$module} = 1;
		    }
		    # Move to the parent class.
		    $class = ($class->{'extends'} eq $directive->{'name'}) ? undef() : $classes{$class->{'extends'}};
		}
		# Find any variables to be excluded from state store/restore.
		@{$stateStore->{'excludes'}} = exists($directive->{'stateStorable'}->{'exclude'}->{'variables'}) ? split(/\s*,\s*/,$directive->{'stateStorable'}->{'exclude'}->{'variables'}) : ();
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
		    &stateStoreVariables($stateStores,$stateStore,undef(),$declaration,$explicitNamesFound);
		}
		# Add any variables declared in the functionClassType class.
		if ( defined($functionClassType) ) {
		    my $node = $functionClassType->{'node'}->{'firstChild'};
		    while ( $node ) {
			&stateStoreVariables($stateStores,$stateStore,undef(),$node->{'declarations'},$explicitNamesFound)
			    if ( $node->{'type'} eq "declaration" );
			$node = $node->{'type'} eq "contains" ? $node->{'firstChild'} : $node->{'sibling'};
		    }
		}
		# Check that all explicit variables were found.
		{
		    my $class = $nonAbstractClass;
		    while ( $class ) {
			if ( exists($class->{'stateStorable'}->{'functionClass'}) && exists($class->{'stateStorable'}->{'functionClass'}->{'variables'}) ) {
			    foreach my $variable ( split(/\s*,\s*/,$class->{'stateStorable'}->{'functionClass'}->{'variables'}) ) {
				die("Error: unable to find variable '".$variable."' marked as state storable in class '".$class->{'name'}."'")
				    unless ( grep {$_ eq lc($variable)} @{$explicitNamesFound} );
			    }
			}
			# Move to the parent class.
			$class = ($class->{'extends'} eq $directive->{'name'}) ? undef() : $classes{$class->{'extends'}};
		    }
		}
		# Add code to method.
		$stateStores->{'stateFileUsed'} = 1
		    if ( scalar(@{$stateStore->{'staticVariables'}}) > 0 );
		if ( $stateStore->{'hasCustomStateStore'  } ) {
		    # The class has its own state store function, so we should never arrive at this point in the code.
		    $stateStoreCode .= " call Error_Report('custom state store function exists - this should not happen'//".&Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($nonAbstractClass->{'node'},$nonAbstractClass->{'node'}->{'line'}).")\n";
		    $stateStores->{'stateStoreModules'}->{'Error'} = 1;
		} else {
		    foreach ( @{$stateStore->{'staticVariables'}} ) {
			$stateStores->{'labelUsed'}       = 1;
			$stateStoreCode .= " if (displayVerbosity() >= verbosityLevelWorking) then\n";
			# <workaround type="gfortran" PR="94446" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=94446">
			#  <description>
			#   Using the sizeof() intrinsic on a treeNode object causes a bogus "type mismatch" error when this module is used.
			#  </description>
			# </workaround>
			$stateStoreCode .= "   write (label,'(i16)') 0\n";
			#$stateStoreCode .= "  write (label,'(i16)') sizeof(self%".$_.")\n";
			$stateStoreCode .= "  call displayMessage('storing \"".$_."\" with size '//trim(adjustl(label))//' bytes')\n";
			$stateStoreCode .= " end if\n";
		    }
		    $stateStoreCode .= " write (stateFile) ".join(", &\n  & ",map {"self%".$_} @{$stateStore->{'staticVariables'}})."\n"
			if ( scalar(@{$stateStore->{'staticVariables'}}) > 0 );
		    $stateStoreCode .= $stateStore->{'outputCode'}
			if ( defined($stateStore->{'outputCode'}) );
		}
		if ( $stateStore->{'hasCustomStateRestore'} ) {
		    # The class has its own state store function, so we should never arrive at this point in the code.
		    $stateRestoreCode .= " call Error_Report('custom state restore function exists - this should not happen'//".&Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($nonAbstractClass->{'node'},$nonAbstractClass->{'node'}->{'line'}).")\n";
		    $stateStores->{'stateRestoreModules'}->{'Error'} = 1;
		} else {
		    foreach ( @{$stateStore->{'staticVariables'}} ) {
			$stateRestoreCode .= " call displayMessage('restoring \"".$_."\"',verbosity=verbosityLevelWorking)\n";
		    }
		    $stateRestoreCode .= " read (stateFile) ".join(", &\n  & ",map {"self%".$_} @{$stateStore->{'staticVariables'}})."\n"
			if ( scalar(@{$stateStore->{'staticVariables'}}) > 0 );
		    $stateRestoreCode .= $stateStore->{'inputCode'}
			if ( defined($stateStore->{'inputCode'}) );
		}
	    }
	    $stateStoreCode   .= "end select\n";
	    $stateStoreCode   .= "call displayUnindent('done',verbosity=verbosityLevelWorking)\n";
	    $stateStoreCode   .= "return\n";
	    $stateRestoreCode .= "end select\n";
	    $stateRestoreCode .= "call displayUnindent('done',verbosity=verbosityLevelWorking)\n";
	    $stateRestoreCode .= "return\n";
	    unless ( $stateStores->{'gslStateFileUsed'} ) {
		push(@outputUnusedVariables,"gslStateFile");
		push(@inputUnusedVariables ,"gslStateFile");
	    }
	    unless ( $stateStores->{'stateFileUsed'}     ) {
		push(@outputUnusedVariables,"stateFile"    );
		push(@inputUnusedVariables ,"stateFile"    );
	    }
	    $stateStoreCode   =
		($stateStores->{'rankMaximum'} > 0 ? " integer :: ".join(", ",map {"i".$_} 1..$stateStores->{'rankMaximum'})."\n" : "").
		(@outputUnusedVariables ? " !\$GLC attributes unused :: ".join(", ",@outputUnusedVariables)."\n" : "").
		$stateStoreCode  ;
	    $stateRestoreCode =
		($stateStores->{'rankMaximum'} > 0 ? " integer :: ".join(", ",map {"i".$_} 1..$stateStores->{'rankMaximum'})."\n" : "").
		(@inputUnusedVariables ? " !\$GLC attributes unused :: ".join(", ",@inputUnusedVariables)."\n" : "").
		$stateRestoreCode;
	    if ( $stateStores->{'allocatablesFound'} ) {
		$stateRestoreCode = ($stateStores->{'dimensionalsFound'} ? "integer(c_size_t), allocatable, dimension(:) :: storedShape\n"  : "").
		    ($stateStores->{'allocatablesFound'} ? "logical                                      :: wasAllocated\n" : "").
		    $stateRestoreCode;
	    }
	    if ( $stateStores->{'explicitFunctionsFound'} ) {
		$stateRestoreCode = "logical :: wasAssociated\n".$stateRestoreCode;
	    }
	    $stateStoreCode   = " character(len=16) :: label\n".$stateStoreCode
                 if ( $stateStores->{'labelUsed'} );
            $stateStoreCode   = " integer(c_size_t) :: position\n".$stateStoreCode;
            $stateRestoreCode = " integer(c_size_t) :: position\n".$stateRestoreCode;
	    $stateStoreCode   = &Fortran::Utils::Format_Variable_Definitions($stateLinkedListVariables).$stateStoreCode;
	    $stateRestoreCode = &Fortran::Utils::Format_Variable_Definitions($stateLinkedListVariables).$stateRestoreCode;
	    $methods{'stateStore'} =
	    {
		description => "Store the state of this object to file.",
		type        => "void",
		pass        => "yes",
		argument    => [ "integer, intent(in   ) :: stateFile", "type(c_ptr), intent(in   ) :: gslStateFile", "integer(c_size_t), intent(in   ) :: stateOperationID"  ],
		code        => "call self%stateStore_(stateFile,gslStateFile,stateOperationID)"
	    };
	    $methods{'stateStore_'} =
	    {
		description => "Store the state of this object to file.",
		type        => "void",
		pass        => "yes",
		modules     => join(" ",sort(keys(%{$stateStores->{'stateStoreModules'}}))),
		argument    => [ "integer, intent(in   ) :: stateFile", "type(c_ptr), intent(in   ) :: gslStateFile", "integer(c_size_t), intent(in   ) :: stateOperationID"  ],
		code        => $stateStoreCode
	    };
	    $methods{'stateRestore'} =
	    {
		description => "Restore the state of this object from file.",
		type        => "void",
		pass        => "yes",
		argument    => [ "integer, intent(in   ) :: stateFile", "type(c_ptr), intent(in   ) :: gslStateFile", "integer(c_size_t), intent(in   ) :: stateOperationID"  ],
		code        => "call self%stateRestore_(stateFile,gslStateFile,stateOperationID)"
	    };
	    $methods{'stateRestore_'} =
	    {
		description => "Restore the state of this object from file.",
		type        => "void",
		pass        => "yes",
		modules     => join(" ",sort(keys(%{$stateStores->{'stateRestoreModules'}}))),
		argument    => [ "integer, intent(in   ) :: stateFile", "type(c_ptr), intent(in   ) :: gslStateFile", "integer(c_size_t), intent(in   ) :: stateOperationID"  ],
		code        => $stateRestoreCode
	    };
	    # Initialize structure that will hold all generated code.
	    my $codeContent;
	    $codeContent->{'module'}->{'preContains'} =
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
	    $codeContent->{'module'}->{'postContains'} =
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
	    $codeContent->{'module'}->{'interfaces'} =
 	       [
	      ];
	    my $modulePreContains  = $codeContent->{'module'}->{'preContains' }->[0];
	    my $modulePostContains = $codeContent->{'module'}->{'postContains'}->[0];

	    # Generate the base class.
	    &Galacticus::Build::SourceTree::SetVisibility($node->{'parent'},$directive->{'name'}."Class","public");
	    &Galacticus::Build::SourceTree::SetVisibility($node->{'parent'},$directive->{'name'}        ,"public");
	    my $extends = exists($directive->{'extends'}) ? $directive->{'extends'} : "functionClass";
	    $modulePreContains->{'content'} .= "   type, extends(".$extends.") :: ".$directive->{'name'}."Class\n";
	    $modulePreContains->{'content'} .= "    private\n";
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
		},
		source     => "Galacticus::Build::SourceTree::Process::FunctionClass::Process_FunctionClass()",
		line       => 1
	    };
            $usesNode->{'moduleUse'}->{'ISO_C_Binding'} = {intrinsic => 1, all => 1}
                if ( $debugging );
            $modulePreContains->{'content'} .= "    integer(c_size_t) :: stateOperationID=0\n";
            $modulePreContains->{'content'} .= "    class(".$directive->{'name'}."Class), public, pointer :: copiedSelf => null()\n";
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
		    $modulePreContains->{'content'} .= $_->{'content'}."\n"
			if (  $_->{'scope'} eq "self" );
		} else {
		    $modulePreContains->{'content'} .= $_."\n";
		}
	    }
	    $modulePreContains->{'content'} .= "    contains\n";
	    $modulePreContains->{'content'} .= "    !![\n";
	    $modulePreContains->{'content'} .= "    <methods>\n";
	    my $generics;
            foreach my $methodName ( sort(keys(%methods)) ) {
                next
                    if ( $methodName eq "destructor" );
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
		    foreach my $intrinsic ( sort(keys(%Fortran::Utils::intrinsicDeclarations)) ) {
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
				$argumentList .= $separator."\\textcolor{red}{\\textless ".latex_encode($intrinsicName);
				$argumentList .= latex_encode($type)
				    if ( defined($type) );
				$argumentList .= "\\textgreater} ".latex_encode($variable);
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
			    last;
			}
		    }
		}
		$modulePreContains->{'content'} .= "     <method method=\"".$methodName."\">\n";
		$modulePreContains->{'content'} .= "      <description>\n";
                $modulePreContains->{'content'} .= join("\n",map {"       ".$_} split("\n",$method->{'description'}))."\n";
                $modulePreContains->{'content'} .= "      </description>\n";
		$modulePreContains->{'content'} .= "     </method>\n";
		if ( exists($directive->{'generic'}) ) {
		    foreach my $generic ( &List::ExtraUtils::as_array($directive->{'generic'}) ) {
			if ( grep {$_ eq $methodName} &List::ExtraUtils::as_array($generic->{'method'}) ) {
			    # This method is part of a generic method, store relevant information.
			    $generics->{$generic->{'name'}}->{'type'} = $method->{'type'};
			    push(@{ $generics->{$generic->{'name'}}->{'description'}},$method->{'description'});
			    push(@{ $generics->{$generic->{'name'}}->{'argumentList'}},$argumentList);
			}
		    }
		}
	    }
	    foreach my $generic ( &List::ExtraUtils::hashList($generics, keyAs => "name") ) {
		$modulePreContains->{'content'} .= "    <method method=\"".$generic->{'name'}."\">\n";
		$modulePreContains->{'content'} .= "       <description>\n";
		$modulePreContains->{'content'} .= "        ".join(" | ",@{$generic->{'description'}})."\n";
		$modulePreContains->{'content'} .= "       </description>\n";
		$modulePreContains->{'content'} .= "    </method>\n";
	    }
	    $modulePreContains->{'content'} .= "    </methods>\n";
	    $modulePreContains->{'content'} .= "    !!]\n";
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
            foreach ( sort(keys(%methods)) ) {
                next
                    if ( $_ eq "destructor" || $_ eq "assignment(=)" );
                my $method = $methods{$_};
                my $functionName;
                if ( exists($method->{'function'}) ) {
		    $functionName = $method->{'function'};
                } else {
		    $functionName = $directive->{'name'}.ucfirst($_).(exists($method->{'code'}) ? "" : "__");
                }
		$methodTable->add("",$_,$functionName);
	    }
	    $modulePreContains->{'content'} .= $methodTable->table();
	    $modulePreContains->{'content'} .= "procedure :: ".$directive->{'name'}."Assignment\n";
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
	    if ( exists($directive->{'generic'}) ) {
		foreach ( &List::ExtraUtils::as_array($directive->{'generic'}) ) {
		    $genericTable->add($_->{'name'},join(", ",&List::ExtraUtils::as_array($_->{'method'})));
		}
            }
	    $genericTable->add("assignment(=)",$directive->{'name'}."Assignment");
	    $modulePreContains->{'content'} .= $genericTable->table();
	    $modulePreContains->{'content'} .= "   final :: ".$directive->{'name'}."Destructor\n"
                if ( exists($methods{'destructor'}) );
	    $modulePreContains->{'content'} .= "   end type ".$directive->{'name'}."Class\n\n";
	    # Insert any module-scope class content.
	    foreach ( &List::ExtraUtils::as_array($directive->{'data'}) ) {
		if ( reftype($_) ) {
		    if ( exists($_->{'scope'}) && $_->{'scope'} eq "module" ) {
			$modulePreContains->{'content'} .= $_->{'content'}."\n";
			if ( exists($_->{'threadprivate'}) && $_->{'threadprivate'} eq "yes" && $_->{'content'} =~ m/::\s*(.*)$/ ) {
			    my @declarations = split(/\s*,\s*/,$1);
			    foreach my $declaration ( @declarations ) {
				$declaration =~ s/\s*=.*//;
			    }
			    $modulePreContains->{'content'} .= "   !\$omp threadprivate(".join(",",@declarations).")\n";
			}
		    }
		}
	    }

	    # Insert state variable for input parameter validation.
	    $modulePreContains->{'content'} .= "   logical :: ".$directive->{'name'}."DsblVldtn=.false.\n";
	    $modulePreContains->{'content'} .= "   !\$omp threadprivate(".$directive->{'name'}."DsblVldtn)\n";

	    # Generate class constructors
	    $modulePreContains->{'content'} .= "   interface ".$directive->{'name'}."\n";
	    $modulePreContains->{'content'} .= "    module procedure ".$directive->{'name'}."CnstrctrPrmtrs\n";
	    $modulePreContains->{'content'} .= "   end interface ".$directive->{'name'}."\n";
	    # Add a variable which records whether construction of the default object is underway, and for detecting and
	    # performing recursive build of the default object from a parameter list.
	    my $allowRecursion = grep {exists($_->{'recursive'}) && $_->{'recursive'} eq "yes"} @classes;
	    if ( $allowRecursion ) {
		$modulePreContains->{'content'} .= "   type(inputParameter), pointer :: ".$directive->{'name'}."RecursiveBuildNode => null()\n";
		$modulePreContains->{'content'} .= "   class(".$directive->{'name'}."Class), pointer :: ".$directive->{'name'}."RecursiveBuildObject => null()\n";
		$modulePreContains->{'content'} .= "   !\$omp threadprivate(".$directive->{'name'}."RecursiveBuildNode,".$directive->{'name'}."RecursiveBuildObject)\n\n";
		my $usesNode =
		{
		    type      => "moduleUse",
		    moduleUse =>
		    {
			Input_Parameters =>
			{
			    intrinsic => 0,
			    all       => 1
			}
		    },
		    source     => "Galacticus::Build::SourceTree::Process::FunctionClass::Process_FunctionClass()",
		    line       => 1
		};
		&Galacticus::Build::SourceTree::Parse::ModuleUses::AddUses($node->{'parent'},$usesNode);
 	    }
	    # Add method name parameter.
	    if ( $tree->{'type'} eq "file" ) {
		(my $fileName = $tree->{'name'}) =~ s/\.F90$/.p/;
		open(my $parametersFile,">>".$ENV{'BUILDPATH'}."/".$fileName);
		print $parametersFile $directive->{'name'}."\n";
		close($parametersFile);
	    }
	    # Create XML constructor.
	    $modulePostContains->{'content'} .= "   ".($allowRecursion ? "recursive " : "")."function ".$directive->{'name'}."CnstrctrPrmtrs(parameters,copyInstance,parameterName) result(self)\n";
	    $modulePostContains->{'content'} .= "      !!{\n";
	    $modulePostContains->{'content'} .= "      Return a pointer to a newly created {\\normalfont \\ttfamily ".$directive->{'name'}."} object as specified by the provided parameters.\n";
	    $modulePostContains->{'content'} .= "      !!}\n";
	    $modulePostContains->{'content'} .= "      use :: Input_Parameters  , only : inputParameter         , inputParameters\n";
	    $modulePostContains->{'content'} .= "      use :: Locks  , only : ompLock\n";
	    $modulePostContains->{'content'} .= "      use :: Error  , only : Error_Report\n";
	    $modulePostContains->{'content'} .= "      use :: ISO_Varying_String, only : varying_string         , char           , trim, operator(//), operator(==), assignment(=)\n";
	    $modulePostContains->{'content'} .= "      implicit none\n";
	    $modulePostContains->{'content'} .= "      class    (".$directive->{'name'}."Class), pointer :: self\n";
	    $modulePostContains->{'content'} .= "      type     (inputParameters), intent(inout)           :: parameters\n";
	    $modulePostContains->{'content'} .= "      integer                   , intent(in   ), optional :: copyInstance\n";
	    $modulePostContains->{'content'} .= "      character(len=*          ), intent(in   ), optional :: parameterName\n";
	    $modulePostContains->{'content'} .= "      type     (inputParameters)                          :: subParameters\n";
	    if ( exists($directive->{'default'}) ) {
		$modulePostContains->{'content'} .= "      type     (inputParameter ), pointer                 :: parameterNode\n";
		$modulePostContains->{'content'} .= "      type     (ompLock        ), save                    :: addLock\n";
		$modulePostContains->{'content'} .= "      logical                   , save                    :: addLockInitialized=.false.\n";
		$modulePostContains->{'content'} .= "      logical                                             :: needLock\n";
	    }
	    $modulePostContains->{'content'} .= "      type     (varying_string )                          :: message      , instanceName, parameterName_\n";
	    $modulePostContains->{'content'} .= "      integer                                             :: copyInstance_\n\n";
	    $modulePostContains->{'content'} .= "      if (present(parameterName)) then\n";
	    $modulePostContains->{'content'} .= "        parameterName_=parameterName\n";
	    $modulePostContains->{'content'} .= "      else\n";
	    $modulePostContains->{'content'} .= "        parameterName_='".$directive->{'name'}."'\n";
	    $modulePostContains->{'content'} .= "      end if\n";
	    $modulePostContains->{'content'} .= "      if (present(copyInstance)) then\n";
	    $modulePostContains->{'content'} .= "        copyInstance_=copyInstance\n";
	    $modulePostContains->{'content'} .= "      else\n";
	    $modulePostContains->{'content'} .= "        copyInstance_=1\n";
	    $modulePostContains->{'content'} .= "      end if\n";
            if ( exists($directive->{'default'}) ) {
                (my $class) = grep {$_->{'name'} eq $directive->{'name'}.ucfirst($directive->{'default'})} @nonAbstractClasses;
		if ( exists($directive->{'default'}) ) {
		    $modulePostContains->{'content'} .= "      if (.not.addLockInitialized) then\n";
		    $modulePostContains->{'content'} .= "      !\$omp critical (addLockInitialize".ucfirst($directive->{'default'}).")\n";
		    $modulePostContains->{'content'} .= "          if (.not.addLockInitialized) then\n";
		    $modulePostContains->{'content'} .= "          addLockInitialized=.true.\n";
		    $modulePostContains->{'content'} .= "          addLock=ompLock()\n";
		    $modulePostContains->{'content'} .= "      end if\n";
		    $modulePostContains->{'content'} .= "      !\$omp end critical (addLockInitialize".ucfirst($directive->{'default'}).")\n";
		    $modulePostContains->{'content'} .= "      end if\n";
		    $modulePostContains->{'content'} .= "      needLock=.not.addLock%ownedByThread()\n";
		    $modulePostContains->{'content'} .= "      if (needLock) call addLock%set()\n";
		}
	        $modulePostContains->{'content'} .= "      if (parameterName_ == '".$directive->{'name'}."' .and. copyInstance_ == 1 .and. .not.parameters%isPresent(char(parameterName_))) then\n";
	        $modulePostContains->{'content'} .= "        call parameters%addParameter('".$directive->{'name'}."','".$directive->{'default'}."')\n";
	        $modulePostContains->{'content'} .= "        parameterNode => parameters%node('".$directive->{'name'}."',requireValue=.true.)\n";
		$modulePostContains->{'content'} .= "        subParameters=parameters%subParameters(char(parameterName_))\n";
		$modulePostContains->{'content'} .= "        allocate(".$directive->{'name'}.ucfirst($directive->{'default'})." :: self)\n";
		if ( exists($class->{'recursive'}) && $class->{'recursive'} eq "yes" ) {
		    $modulePostContains->{'content'} .= "        ".$directive->{'name'}."RecursiveBuildNode   => parameterNode\n";
		    $modulePostContains->{'content'} .= "        ".$directive->{'name'}."RecursiveBuildObject => self\n";
		}
		$modulePostContains->{'content'} .= "        select type (self)\n";
		$modulePostContains->{'content'} .= "          type is (".$directive->{'name'}.ucfirst($directive->{'default'}).")\n";
		$modulePostContains->{'content'} .= "            call debugStackPush(loc(self))\n"
		    if ( $debugging );
		$modulePostContains->{'content'} .= "            self=".$directive->{'name'}.ucfirst($directive->{'default'})."(subParameters)\n";
		$modulePostContains->{'content'} .= "            call debugStackPop()\n"
		    if ( $debugging );
		$modulePostContains->{'content'} .= "         end select\n";
		if ( exists($class->{'recursive'}) && $class->{'recursive'} eq "yes" ) {
		    $modulePostContains->{'content'} .= "        ".$directive->{'name'}."RecursiveBuildNode   => null()\n";
		    $modulePostContains->{'content'} .= "        ".$directive->{'name'}."RecursiveBuildObject => null()\n";
		}
                $modulePostContains->{'content'} .= "         call parameterNode%objectSet(self)\n";
                $modulePostContains->{'content'} .= "      else\n";
            }
	    # Detect recursive builds if any class member allows it.
	    if ( $allowRecursion ) {
		$modulePostContains->{'content'} .= "        parameterNode => parameters%node(char(parameterName_),requireValue=.true.)\n";
		$modulePostContains->{'content'} .= "        if (associated(parameterNode,".$directive->{'name'}."RecursiveBuildNode)) then\n";
		foreach my $class ( @nonAbstractClasses ) {
		    next
			unless ( exists($class->{'recursive'}) && $class->{'recursive'} eq "yes" );
		    $modulePostContains->{'content'} .= "           select type (".$directive->{'name'}."RecursiveBuildObject)\n";
		    $modulePostContains->{'content'} .= "              type is (".$class->{'name'}.")\n";
		    $modulePostContains->{'content'} .= "              allocate(".$class->{'name'}." :: self)\n";
		    $modulePostContains->{'content'} .= "              select type (self)\n";
		    $modulePostContains->{'content'} .= "              type is (".$class->{'name'}.")\n";
		    $modulePostContains->{'content'} .= "                 self%isRecursive=.true.\n";
		    $modulePostContains->{'content'} .= "                 self%recursiveSelf => ".$directive->{'name'}."RecursiveBuildObject\n";
		    $modulePostContains->{'content'} .= "              end select\n";
		    $modulePostContains->{'content'} .= "           end select\n";
		}
		$modulePostContains->{'content'} .= "           if (needLock) call addLock%unset()\n";
		$modulePostContains->{'content'} .= "           return\n";
		$modulePostContains->{'content'} .= "        end if\n";
	    }
	    # Build the object.
	    $modulePostContains->{'content'} .= "      call parameters%value(char(parameterName_),instanceName,copyInstance=copyInstance_)\n";
	    $modulePostContains->{'content'} .= "      subParameters=parameters%subParameters(char(parameterName_),copyInstance=copyInstance_)\n";
	    $modulePostContains->{'content'} .= "      select case (char(instanceName))\n";
	    foreach my $class ( @nonAbstractClasses ) {
		(my $name = $class->{'name'}) =~ s/^$directive->{'name'}//;
		$name = lcfirst($name)
		    unless ( $name =~ m/^[A-Z]{2,}/ );
		$modulePostContains->{'content'} .= "     case ('".$name."')\n";
		$modulePostContains->{'content'} .= "        allocate(".$class->{'name'}." :: self)\n";
		if ( exists($class->{'recursive'}) && $class->{'recursive'} eq "yes" ) {
		    $modulePostContains->{'content'} .= "        ".$directive->{'name'}."RecursiveBuildNode   => parameterNode\n";
		    $modulePostContains->{'content'} .= "        ".$directive->{'name'}."RecursiveBuildObject => self\n";
		}
		$modulePostContains->{'content'} .= "        select type (self)\n";
		$modulePostContains->{'content'} .= "          type is (".$class->{'name'}.")\n";
		$modulePostContains->{'content'} .= "            call debugStackPush(loc(self))\n"
		    if ( $debugging );
		$modulePostContains->{'content'} .= "            self=".$class->{'name'}."(subParameters)\n";
		$modulePostContains->{'content'} .= "            call debugStackPop()\n"
		    if ( $debugging );
		$modulePostContains->{'content'} .= "         end select\n";
		if ( exists($class->{'recursive'}) && $class->{'recursive'} eq "yes" ) {
		    $modulePostContains->{'content'} .= "        ".$directive->{'name'}."RecursiveBuildNode   => null()\n";
		    $modulePostContains->{'content'} .= "        ".$directive->{'name'}."RecursiveBuildObject => null()\n";
		}
	    }
	    $modulePostContains->{'content'} .= "      case default\n";
	    $modulePostContains->{'content'} .= "         message='Unrecognized type \"'//trim(instanceName)//'\" Available options are:'\n";
	    my @classNames;
	    push(@classNames,$_->{'name'})
		foreach ( @nonAbstractClasses );
	    foreach ( sort(@classNames) ) {
		(my $name = $_) =~ s/^$directive->{'name'}//;
		$name = lcfirst($name)
		    unless ( $name =~ m/^[A-Z]{2,}/ );
		$modulePostContains->{'content'} .= "         message=message//char(10)//'   -> ".$name."'\n";
	    }
	    $modulePostContains->{'content'} .= "         call Error_Report(message//".&Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($node,$node->{'line'}).")\n";
            $modulePostContains->{'content'} .= "      end select\n";
	    if ( exists($directive->{'default'}) ) {
		$modulePostContains->{'content'} .= "      end if\n";
		$modulePostContains->{'content'} .= "      if (needLock) call addLock%unset()\n";
	    }
	    $modulePostContains->{'content'} .= "      return\n";
	    $modulePostContains->{'content'} .= "   end function ".$directive->{'name'}."CnstrctrPrmtrs\n\n";

	    # Insert class code. Each class goes into its own submodule, so we build this code inside a separate part of the code
	    # structure for each class. This is initialized below.
            foreach my $class ( @classes ) {		  
                &Galacticus::Build::SourceTree::SetVisibility($node->{'parent'},$class->{'type'},"public");
                ($codeContent->{'submodule'}->{$class->{'type'}}->{'fileName'   } = $class->{'file'}) =~ s/^.*\/(.*)\.F90$/$1.p.F90/;
                $codeContent ->{'submodule'}->{$class->{'type'}}->{'fileName'   } = $ENV{'BUILDPATH'}."/".$codeContent->{'submodule'}->{$class->{'type'}}->{'fileName'};
                $codeContent ->{'submodule'}->{$class->{'type'}}->{'preContains'} =
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
                $codeContent->{'submodule'}->{$class->{'type'}}->{'postContains'} =
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
                @{$codeContent->{'submodule'}->{$class->{'type'}}->{'interfaces'}} = ();
                my $submodulePreContains  = $codeContent->{'submodule'}->{$class->{'type'}}->{'preContains' };
                my $submodulePostContains = $codeContent->{'submodule'}->{$class->{'type'}}->{'postContains'};
                my @moduleScoped;
		my @moduleSymbols;
		my @moduleUseNodes;
		my $classTree = $class    ->{'tree'      };
		my $classNode = $classTree->{'firstChild'};
		my $contained = 0;
                # Walk the tree containing the class node.
		while ( $classNode ) {
		    if ( $classNode->{'type'} eq "contains" ) {
			# Record that we've found the "contains" divider, and move to the contained content.
			$classNode = $classNode->{'firstChild'};
			$contained = 1;
		    }
		    if ( $contained ) {
			# Contained content always goes into the submodule.
			push(@{$submodulePostContains},$classNode);
			if ( $classNode->{'type'} eq "function" || $classNode->{'type'} eq "subroutine" ) {
			    # For functions and subroutines we must determine if an interface must be inserted into the parent
			    # module. We have built a list of all named functions/subroutines that require such interfaces, so
			    # simply look for the current name in that list.
			    if ( grep {$_ eq lc($classNode->{'name'})} @{$codeContent->{'submodule'}->{$class->{'type'}}->{'interfaces'}} ) {
				# Must add the interface for this function to the module.
				my @interfaceNodes;
				my $interfaceOpener =
				{
				    type       => "code" ,
				    content    => ""     ,
				    firstChild => undef(),
				    sibling    => undef(),
				    parent     => undef(),
				    source     => "Galacticus::Build::SourceTree::Process::FunctionClass::Process_FunctionClass()",
				    line       => 1
				};
				$interfaceOpener->{'content'} .= "interface\n";
				$interfaceOpener->{'content'} .= "module ".$classNode->{'opener'};
				push(@interfaceNodes,$interfaceOpener);
				my $returnName   = $classNode->{'opener'    } =~ m/result\s*\(\s*([a-zA-Z0-9_]+)\s*\)/ ? $1 : $classNode->{'name'};
				my $functionNode = $classNode->{'firstChild'};
				# Walk through the function/subroutine and process declarations.
				while ( $functionNode ) {
				    if ( $functionNode->{'type'} eq "declaration" ) {
					my @declarationsLocal;
					foreach my $declaration ( @{$functionNode->{'declarations'}} ) {
					    if ( ( grep {$_ eq lc($returnName)} @{$declaration->{'variables'}}) || ( grep {$_ =~ m/intent\s*\(\s*(in|out|inout)\s*\)/i } @{$declaration->{'attributes'}} ) ) {
						# Function result variables and any variables with "intent()" attributes must
						# always be specified in the interface.
						die('Galacticus::Build::SourceTree::Process::FunctionClass::Process_FunctionClass: expected a "code" node as first child')
						    unless ( $functionNode->{'firstChild'}->{'type'} eq "code" );
						my $interfaceDeclaration =
						{
						    type       => "code" ,
						    content    => ""     ,
						    firstChild => undef(),
						    sibling    => undef(),
						    parent     => undef(),
						    source     => $functionNode->{'firstChild'}->{'source'},
						    line       => $declaration                 ->{'line'  }
						};
						$interfaceDeclaration->{'content'} .= &Fortran::Utils::Format_Variable_Definitions([$declaration]);
						push(@interfaceNodes,$interfaceDeclaration);
						# Capture any type names - these will need to be imported into the parent module.
						if  ( defined($declaration->{'type'}) ) {
						    # A type is defined - strip initial part before any "=" and any whitespace.
						    (my $type = $declaration->{'type'}) =~ s/^.+=//;
						    $type =~ s/\s//g;
						    # Ignore purely numerical types (e.g. the "32" in a "len=32' character type).
						    push(@moduleSymbols,$type)
							unless ( $type =~ m/^\d+$/ );
						}
					    } elsif ( (grep {$_ =~ m/external/i } @{$declaration->{'attributes'}}) || $declaration->{'intrinsic'} eq "procedure" ) {
						# For declarations of "external" symbols we must determine if they appear as a
						# dummy argument - if they do they must go into the interface.
						if ( my @matches = $functionNode->{'parent'}->{'opener'} =~ $Fortran::Utils::unitOpeners{$functionNode->{'parent'}->{'type'}}->{'regEx'} ) {
						    my @arguments = split(/\s*,\s*/,lc($matches[$Fortran::Utils::unitOpeners{$functionNode->{'parent'}->{'type'}}->{'arguments'}]))
							if (
							    exists (         $Fortran::Utils::unitOpeners{$functionNode->{'parent'}->{'type'}}->{'arguments'} )
							    &&
							    defined($matches[$Fortran::Utils::unitOpeners{$functionNode->{'parent'}->{'type'}}->{'arguments'}])
							);
						    my $moduleScope    = dclone($declaration);
						    my $submoduleScope = dclone($declaration);
						    @{$moduleScope   ->{'variables'    }} = ();
						    @{$moduleScope   ->{'variableNames'}} = ();
						    @{$submoduleScope->{'variables'    }} = ();
						    @{$submoduleScope->{'variableNames'}} = ();
						    for(my $i=0;$i<scalar(@{$declaration->{'variableNames'}});++$i) {
							my $variableName = $declaration->{'variableNames'}->[$i];
							if ( grep {$_ eq lc($variableName)} @arguments ) {
							    push(@{$moduleScope   ->{'variables'    }},$declaration->{'variables'    }->[$i]);
							    push(@{$moduleScope   ->{'variableNames'}},$declaration->{'variableNames'}->[$i]);
							} else {
							    push(@{$submoduleScope->{'variables'    }},$declaration->{'variables'    }->[$i]);
							    push(@{$submoduleScope->{'variableNames'}},$declaration->{'variableNames'}->[$i]);
							}
						    }
						    push(@declarationsLocal,$submoduleScope)
							if ( scalar(@{$submoduleScope->{'variables'}}) > 0 );
						    if ( scalar(@{$moduleScope   ->{'variables'}}) > 0 ) {
							die('Galacticus::Build::SourceTree::Process::FunctionClass::Process_FunctionClass: expected a "code" node as first child')
							    unless ( $functionNode->{'firstChild'}->{'type'} eq "code" );
							my $interfaceDeclaration =
							{
							    type       => "code" ,
							    content    => ""     ,
							    firstChild => undef(),
							    sibling    => undef(),
							    parent     => undef(),
							    source     => $functionNode->{'firstChild'}->{'source'},
							    line       => $declaration                 ->{'line'  }
							};
							$interfaceDeclaration->{'content'} .= &Fortran::Utils::Format_Variable_Definitions([$moduleScope]);
							push(@interfaceNodes,$interfaceDeclaration);
							# Capture any type names - these will need to be imported into the parent module.
						   	if  ( defined($declaration->{'type'}) ) {
							    # A type is defined - strip initial part before any "=" and any whitespace.
							    (my $type = $declaration->{'type'}) =~ s/^.+=//;
							    $type =~ s/\s//g;
							    # Ignore purely numerical types (e.g. the "32" in a "len=32' character type).
							    push(@moduleSymbols,$type)
								unless ( $type =~ m/^\d+$/ );
							}
						    }
						}
					    } else {
						# Any other declarations are local to the function/subroutine in the module.
						push(@declarationsLocal,$declaration);
					    }
					}
					# Insert submodule and module declarations.
					@{$functionNode->{'declarations'}} = @declarationsLocal;
					&Galacticus::Build::SourceTree::Parse::Declarations::BuildDeclarations($functionNode);
				    } elsif ( $functionNode->{'type'} eq "moduleUse" ) {
					# Record this moduleUse node for later inclusion into the module with a reduced set of
					# symbols, corresponding to types that are used in interfaces.
					my $content =
					{
					    type       => "code"        ,
					    content    => $functionNode->{'firstChild'}->{'content'},
					    sibling    => undef()       ,
					    parent     => undef()       ,
					    firstChild => undef()       ,
					    source     => "Galacticus::Build::SourceTree::Process::FunctionClass::Process_FunctionClass()",
					    line       => 1
					};
					my $moduleUseNode =
					{
					    type       => "moduleUse",
					    sibling    => undef()    ,
					    parent     => undef()    ,
					    firstChild => $content   ,
					    source     => "Galacticus::Build::SourceTree::Process::FunctionClass::Process_FunctionClass()",
					    line       => 1          ,
					    moduleUse  => dclone($functionNode->{'moduleUse' })
					};
					$content->{'parent'} = $moduleUseNode;
					push(@moduleUseNodes,$moduleUseNode);
				    }
				    $functionNode = $functionNode->{'sibling'};
				}
				my $interfaceCloser =
				{
				    type       => "code" ,
				    content    => ""     ,
				    firstChild => undef(),
				    sibling    => undef(),
				    parent     => undef(),
				    source     => "Galacticus::Build::SourceTree::Process::FunctionClass::Process_FunctionClass()",
				    line       => 1
				};
				$interfaceCloser->{'content'} .= $classNode->{'closer'};
				$interfaceCloser->{'content'} .= "end interface\n";
				$classNode      ->{'opener' } = "module procedure ".$classNode->{'name'}."\n";
				$classNode      ->{'closer' } = "end procedure "   .$classNode->{'name'}."\n";
				push(@interfaceNodes,$interfaceCloser);
				# Insert this new interfaces into the interface list.
				push(@{$codeContent->{'module'}->{'interfaces'}},@interfaceNodes);
			    }
			}
		    } else {
			if ( $classNode->{'type'} eq "scoping" ) {
			    # Capture any variables that are to have module-scope.
			    push(@moduleScoped,split(/\s*,\s*/,$classNode->{'directive'}->{'module'}->{'variables'}))
				if ( exists($classNode->{'directive'}->{'module'}->{'variables'}) );
			} elsif ( $classNode->{'type'} eq "moduleUse" ) {
			    # Any module use statements must be placed in the parent module.
			    &Galacticus::Build::SourceTree::Parse::ModuleUses::AddUses($node->{'parent'},$classNode);
			} elsif (
			    $classNode->{'type'} eq "type"
			    ||
			    $classNode->{'type'} eq "interface"
			    ||
			    $classNode->{'type'} eq $directive->{'name'}
			    ) {
			    # Type declarations, interface, and the functionClass directive itself are placed into the module.
			    push(@{$codeContent->{'module'}->{'interfaces'}},$classNode);
			    if ( $classNode->{'type'} eq "interface" ) {
				# Look for module procedures defined in this interface - an explicit interface for these
				# functions/subroutines must be built.
				my $interfaceNode = $classNode->{'firstChild'};
				while ( $interfaceNode ) {
				    if ( $interfaceNode->{'type'} eq "moduleProcedure" ) {
					push(@{$codeContent->{'submodule'}->{$class->{'type'}}->{'interfaces'}},map {lc($_)} @{$interfaceNode->{'names'}});
				    }
				    $interfaceNode = $interfaceNode->{'sibling'};
				}
			    } elsif ( $classNode->{'type'} eq "type" ) {
				# Look for any type-bound functions as these require an interface.
				my $postContains = 0;
				my $typeNode = $classNode->{'firstChild'};
				while ( $typeNode ) {
				    if ( $typeNode->{'type'} eq "contains" ) {
					$postContains = 1;
					$typeNode = $typeNode->{'firstChild'};
					last
					    unless ( $typeNode );
				    }
				    if ( $postContains ) {
					if ( $typeNode->{'type'} eq "declaration" ) {
					    # Handle regular type-bindings.
					    foreach my $declaration ( @{$typeNode->{'declarations'}} ) {
						push(@{$codeContent->{'submodule'}->{$class->{'type'}}->{'interfaces'}},map {$_ =~ m/=>([a-z0-9_]+)/ ? $1 : $_} @{$declaration->{'variables'}});
					    }
					} elsif ( $typeNode->{'type'} eq "code" && $typeNode->{'content'} =~ m/final\s*::\s*([a-zA-Z0-9_]+)/ ) {
					    # Handle "final" bindings separately.
					    push(@{$codeContent->{'submodule'}->{$class->{'type'}}->{'interfaces'}},lc($1));
					}
				    } else {
					# Pre-contains.
					if ( $typeNode->{'type'} eq "declaration" ) {
					    # Any types defined in declarations must be included at module-scope.
					    foreach my $declaration ( @{$typeNode->{'declarations'}} ) {
						if  ( defined($declaration->{'type'}) ) {
						    # A type is defined - strip initial part before any "=" and any whitespace.
						    (my $type = $declaration->{'type'}) =~ s/^.+=//;
						    $type =~ s/\s//g;
						    # Ignore purely numerical types (e.g. the "32" in a "len=32' character type).
						    push(@moduleSymbols,$type)
							unless ( $type =~ m/^\d+$/ );
						}
					    }
					}
				    }
				    $typeNode = $typeNode->{'sibling'};
				}
			    }
			} elsif ( $classNode->{'type'} eq "visibility" ) {
			    # Visibility statements must go in the module.
			    push(@{$codeContent->{'submodule'}->{$class->{'type'}}->{'interfaces'}},map {lc($_)} sort(keys(%{$classNode->{'visibility'}->{'public'}})));
			    delete($classNode->{'visibility'}->{'private'})
				if ( exists($classNode->{'visibility'}->{'private'}) );
			    &Galacticus::Build::SourceTree::Parse::Visibilities::UpdateVisibilities($classNode);
			    push(@{$codeContent->{'module'}->{'interfaces'}},$classNode);
			} elsif ( $classNode->{'type'} eq "declaration" ) {
			    # Variables declared pre-contains are usually placed in the submodule, unless they are publicly
			    # visible, or if explicitly stated to have module-scope.
			    my @declarationsSubmodule;
			    foreach my $declaration ( @{$classNode->{'declarations'}} ) {
				my $moduleScope = 0;
				if ( grep {lc($_) eq "public"} @{$declaration->{'attributes'}} ) {
				    # Public variables must go in the module.
				    my $declarationNew = {
					type       => "declaration",
					sibling    => undef()      ,
					parent     => undef()      ,
					source     => "Galacticus::Build::SourceTree::Process::FunctionClass::Process_FunctionClass()",
					line       => 1
				    };
				    $declarationNew->{'firstChild'} =
				    {
					type       => "code"         ,
					content    => ""             ,
					sibling    => undef()        ,
					parent     => $declarationNew,
					firstChild => undef()        ,
					source     => "Galacticus::Build::SourceTree::Process::FunctionClass::Process_FunctionClass()",
					line       => 1
				    };
				    push(@{$declarationNew->{'declarations'}},$declaration);
				    push(@{$codeContent->{'module'}->{'interfaces'}},$declarationNew);
				    &Galacticus::Build::SourceTree::Parse::Declarations::BuildDeclarations($declarationNew);
				    $moduleScope = 1;
				} else {
				    # Private variables go in the module only if explicitly scoped.
				    my @moduleVariables;
				    my @submoduleVariables;
				    for(my $i=0;$i<scalar(@{$declaration->{'variables'}});++$i) {
					my $variable = $declaration->{'variableNames'}->[$i];
					if ( grep {lc($_) eq lc($variable)} @moduleScoped ) {
					    push(@moduleVariables   ,{variable => $declaration->{'variables'}->[$i], variableName => $declaration->{'variableNames'}->[$i]});
					} else {
					    push(@submoduleVariables,{variable => $declaration->{'variables'}->[$i], variableName => $declaration->{'variableNames'}->[$i]});
					}
				    }
				    if ( scalar(@moduleVariables) > 0 ) {
					my $moduleDeclaration = dclone($declaration);
					@{$moduleDeclaration->{'variables'    }} = map {$_->{'variable'    }} @moduleVariables;
					@{$moduleDeclaration->{'variableNames'}} = map {$_->{'variableName'}} @moduleVariables;
					my $declarationNew = {
					    type       => "declaration",
					    sibling    => undef()      ,
					    parent     => undef()      ,
					    source     => "Galacticus::Build::SourceTree::Process::FunctionClass::Process_FunctionClass()",
					    line       => 1
					};
					$declarationNew->{'firstChild'} =
					{
					    type       => "code"         ,
					    content    => ""             ,
					    sibling    => undef()        ,
					    parent     => $declarationNew,
					    firstChild => undef()        ,
					    source     => "Galacticus::Build::SourceTree::Process::FunctionClass::Process_FunctionClass()",
					    line       => 1
					};
					push(@{$declarationNew->{'declarations'}},$moduleDeclaration);
					push(@{$codeContent->{'module'}->{'interfaces'}},$declarationNew);
					&Galacticus::Build::SourceTree::Parse::Declarations::BuildDeclarations($declarationNew);
					$moduleScope = 1;
				    }
				    if ( scalar(@submoduleVariables) > 0 ) {
					# For submodule variables we strip any visibility attribute since submodule variables are
					# private by definition and so can not have such attributes.
					my $submoduleDeclaration = dclone($declaration);
					@{$submoduleDeclaration->{'variables'    }} = map {$_->{'variable'    }} @submoduleVariables;
					@{$submoduleDeclaration->{'variableNames'}} = map {$_->{'variableName'}} @submoduleVariables;
					@{$submoduleDeclaration->{'attributes'   }} = grep {$_ ne "public" && $_ ne "private"} @{$declaration->{'attributes'}};
					push(@declarationsSubmodule,$submoduleDeclaration);
				    }
				}
				if ( $moduleScope ) {
				    # Some variables from this declaration have module-scope - record their type.
				    if  ( defined($declaration->{'type'}) ) {
					# A type is defined - strip initial part before any "=" and any whitespace.
					(my $type = $declaration->{'type'}) =~ s/^.+=//;
					$type =~ s/\s//g;
					# Ignore purely numerical types (e.g. the "32" in a "len=32' character type).
					push(@moduleSymbols,$type)
					    unless ( $type =~ m/^\d+$/ );
				    }
				}
			    }
			    # Add all accumulated declarations to the submodule.
			    @{$classNode->{'declarations'}} = @declarationsSubmodule;
			    &Galacticus::Build::SourceTree::Parse::Declarations::BuildDeclarations($classNode);
			    push(@{$submodulePreContains},$classNode);
			} else {
			    # Anything else pre-contains goes in the submodule.
			    push(@{$submodulePreContains},$classNode);
			}
		    }
		    $classNode = $classNode->{'sibling'};
		}
		# Insert module use statements to import any type symbols needed by interfaces into the module.
		@moduleSymbols = uniq(sort(map {lc($_)} @moduleSymbols));
		foreach my $moduleUseNode ( @moduleUseNodes ) {
		    # Remove symbols not used.
		    foreach my $moduleName ( sort(keys(%{$moduleUseNode->{'moduleUse'}})) ) {
			if ( exists($moduleUseNode->{'moduleUse'}->{$moduleName}->{'all'}) ) {
			    # All symbols imported - keep this module for now - in principle we would like to have no cases were all symbols are imported though!
			} else {
			    foreach my $symbolName ( sort(keys(%{$moduleUseNode->{'moduleUse'}->{$moduleName}->{'only'}})) ) {
				delete($moduleUseNode->{'moduleUse'}->{$moduleName}->{'only'}->{$symbolName})
				    unless ( grep {$_ eq lc($symbolName)} @moduleSymbols );
			    }
			    delete($moduleUseNode->{'moduleUse'}->{$moduleName})
				if ( scalar(keys(%{$moduleUseNode->{'moduleUse'}->{$moduleName}->{'only'}})) == 0 );
			}
		    }
		    &Galacticus::Build::SourceTree::Parse::ModuleUses::AddUses($node->{'parent'},$moduleUseNode);
		}
	    }
	    # Insert bindings for any hashed source code descriptors.
	    foreach my $nonAbstractClass ( @nonAbstractClasses ) {
                $modulePreContains->{'content'} .= $code::digest = &Galacticus::Build::SourceTree::Process::SourceDigest::Binding($nonAbstractClass->{'name'});
	    }
	    # C-types are required for the bindings.
	    my $bindingNode =
	    {
		type       => "moduleUse",
	        sibling    => undef()    ,
	        parent     => undef()    ,
		firstChild => undef()    ,
                moduleUse  => {ISO_C_Binding => {intrinsic => 1, only => {C_Char => 1}}},
		source     => "Galacticus::Build::SourceTree::Process::FunctionClass::Process_FunctionClass()",
		line       => 1
	    };
            &Galacticus::Build::SourceTree::Parse::ModuleUses::AddUses($node->{'parent'},$bindingNode);
	    # Create functions.
	    foreach my $methodName ( sort(keys(%methods)) ) {
                my $method = $methods{$methodName};
                next
                    if ( exists($method->{'function'}) );
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
                    my $intrinsic = $methodName eq "destructor" ? "type" : "class";
		    $argumentCode .= "      ".$intrinsic."(".$directive->{'name'}."Class), intent(".(exists($method->{'selfIntent'}) ? $method->{'selfIntent'} : "inout").")";
		    $argumentCode .= ", target"
			if ( exists($method->{'selfTarget'}) && $method->{'selfTarget'} eq "yes" );
		    $argumentCode .= " :: self\n";
		}
		my $separator = "";
                my @unusedVariables = ( "self" );
		foreach my $argument ( @arguments ) {
		    (my $variables = $argument) =~ s/^.*::\s*(.*?)\s*$/$1/;
		    my $isOpenMP = $argument =~ m/^\s*!\$/;
		    $argumentList .= ($isOpenMP ? " &\n!\$ & " : "").$separator.$variables.($isOpenMP ? " &\n& " : "");
		    $argumentCode .= "      ".$argument."\n";
		    $separator     = ",";
		    my $declaration = &Fortran::Utils::Unformat_Variables($argument);
		    push(@unusedVariables,@{$declaration->{'variables'}});
		}
		my $type;
		my $category;
		my $self;
		my $extension = "__";
		$extension = ""
		    if ( exists($method->{'code'}) );
		my $nonOperatorName;
		if ( $methodName eq "assignment(=)" ) {
		    $nonOperatorName = "assignment";
		} else {
		    $nonOperatorName = $methodName;
		}
		my $functionName = $directive->{'name'}.ucfirst($nonOperatorName).$extension;
		my $recursive = exists($method->{'recursive'}) && $method->{'recursive'} eq "yes" ? "recursive " : "";
		my $elemental = exists($method->{'elemental'}) && $method->{'elemental'} eq "yes" ? "elemental " : "";
		if ( $method->{'type'} eq "void" ) {
		    $category = "subroutine";
		    $type     = "";
		    $self     = "";
		} elsif ( $method->{'type'} =~ m/^class/ ) {
		    $category = "function";
		    $type     = "";
		    $self     = "      ".$method->{'type'}.", pointer :: ".$functionName."\n";
		} elsif ( $method->{'type'} =~ m/^type/ || $method->{'type'} =~ m/,/ ) {
		    $category = "function";
		    $type     = "";
		    $self     = "      ".$method->{'type'}." :: ".$functionName."\n";
		} else {
		    $category = "function";
		    $type     = $method->{'type'}." ";
		    $self     = "";
		}
		$modulePostContains->{'content'} .= "   ".$recursive.$elemental.$type.$category." ".$functionName."(self";
		$modulePostContains->{'content'} .= ",".$argumentList
		    unless ( $argumentList eq "" );
		$modulePostContains->{'content'} .= ")\n";
                $modulePostContains->{'content'} .= "      !!{\n";
                $modulePostContains->{'content'} .= "      Default implementation of the {\\normalfont \\ttfamily ".$methodName."} method for the {\\normalfont \\ttfamily ".$directive->{'name'}."} class.\n";
                $modulePostContains->{'content'} .= "      !!}\n";
		if ( exists($method->{'code'}) ) {
		    if ( exists($method->{'modules'}) ) {
			if ( reftype($method->{'modules'}) ) {
			    # Array of modules, with possible "only" clauses.
			    my @modules = reftype($method->{'modules'}) eq "ARRAY" ? @{$method->{'modules'}} : &List::ExtraUtils::hashList($method->{'modules'},keyAs => "name");
			    foreach my $module ( @modules ) {
				my $moduleName = $module->{'name'};
				my $prefix = $moduleName =~ m/OMP_Lib/ ? "!\$ " : "";
				$modulePostContains->{'content'} .= "      ".$prefix."use ".$moduleName.(exists($module->{'only'}) ? ", only : ".join(",",&List::ExtraUtils::as_array($module->{'only'})) : "")."\n";
			    }
			} else {
			    # Simple space-separated list of modules.
			    foreach ( split(/\s+/,$method->{'modules'}) ) {
				my $prefix = $_ =~ m/OMP_Lib/ ? "!\$ " : "";
				$modulePostContains->{'content'} .= "      ".$prefix."use".($_ eq "ISO_C_Binding" ? ", intrinsic :: " : "")." ".$_."\n";
			    }
			}
		    }
		} else {
		    $modulePostContains->{'content'} .= "      use Error             , only : Error_Report\n";
		    $modulePostContains->{'content'} .= "      use ISO_Varying_String, only : char\n";
		}
		$modulePostContains->{'content'} .= "      implicit none\n";
		$modulePostContains->{'content'} .= $argumentCode;
		$modulePostContains->{'content'} .= $self;
		$modulePostContains->{'content'} .= "!\$GLC attributes unused :: ".join(", ",@unusedVariables)."\n";
		if ( exists($method->{'code'}) ) {
		    my $code = "      ".$method->{'code'};
		    $code =~ s/\n/\n      /g;
		    $modulePostContains->{'content'} .= $code."\n";
		} else {
		    $modulePostContains->{'content'} .= "      call Error_Report('this is a null method - initialize the ".$directive->{'name'}." object before use and/or check that the \"'//char(self%objectType())//'\" class implements this method'//".&Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($node,$node->{'line'}).")\n";
		    if ( $category eq "function" ) {
			# Avoid warnings about unset function values.
			$modulePostContains->{'content'} .= "      ".$functionName."=";
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
			$modulePostContains->{'content'} .= $setValue."\n";
		    }
		    $modulePostContains->{'content'} .= "      return\n";
		}
		$modulePostContains->{'content'} .= "   end ".$category." ".$functionName."\n\n";
	    }

	    # Generate documentation. We construct two sets of documentation, one describing the physics models, and one describing the code implementation.
            my $documentationPhysics = "\\section{"      .$directive->{'descriptiveName'}."}\\label{phys:".$directive->{'name'}."}\\hyperdef{physics}{".$directive->{'name'}."}{}\n\n";
	    if ( exists($directive->{'default'}) ) {
		$documentationPhysics .= "Default implementation: \\refPhysics{".$directive->{'name'}.ucfirst($directive->{'default'})."}\n\n";
	    } else {
		$documentationPhysics .= "No default implementation\n\n";
	    }
	    foreach my $className ( sort {lc($a) cmp lc($b)} keys(%classes) ) {
		my $class = $classes{$className};
                (my $suffix = $class->{'name'}) =~ s/^$directive->{'name'}//;
                $suffix = lcfirst($suffix)
                    unless ( $suffix =~ m/^[A-Z]{2}/ );
                $documentationPhysics .= "\\subsection{\\normalfont \\ttfamily ".$suffix."}\\label{phys:".$class->{'name'}."}\\hyperdef{physics}{".$class->{'name'}."}{}\n\n";
                $documentationPhysics .= $class->{'description'}."\n\n";
		$documentationPhysics .= "\\noindent \\textbf{(Default)}\n\n"
		    if ( exists($directive->{'default'}) && $directive->{'name'}.ucfirst($directive->{'default'}) eq $class->{'name'} );
                $documentationPhysics .= "\\noindent \\emph{Implemented by} \\refClass{".$class->{'name'}."}\n";
		# Search the tree for this class to find the interface to the parameters constructor.
		my $node = $classes{$className}->{'tree'}->{'firstChild'};
		$node = $node->{'sibling'}
		    while ( $node && ( $node->{'type'} ne "interface" || ( ! exists($node->{'name'}) || $node->{'name'} ne $className) ) );
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
		$node = $classes{$className}->{'tree'}->{'firstChild'};
		my @objects;
		my @parameters;
		while ( $node ) {
		    # Identify constructor functions.
		    if ( $node->{'type'} eq "function" && (grep {$_ eq $node->{'name'}} @constructors) ) {
			my $constructorNode = $node->{'firstChild'};
			my $depth           = 0;
			while ( $constructorNode ) {
			    # Process node.			  
			    if ( $constructorNode->{'type'} eq "inputParameter" ) {
				# Get the associated variable declaration.
				my $declaration;
				my $variableName = exists($constructorNode->{'directive'}->{'variable'}) ? $constructorNode->{'directive'}->{'variable'} : $constructorNode->{'directive'}->{'name'};
				if ( $variableName =~ m/([a-zA-Z0-9_]+)(\s*\(\s*[a-zA-Z0-9_:,]\s*\)\s*)??\%([a-zA-Z0-9_]+)/ ) {
				    my $objectName         = $1;
				    my $objectVariableName = $3;
				    if ( $objectName eq "self" || $objectName eq $node->{'name'} ) {
					my $parentClass = $class;
					while ( $parentClass ) {
					    my $node = $parentClass->{'tree'}->{'firstChild'};
					    $node = $node->{'sibling'}
					    while ( $node && ( $node->{'type'} ne "type" || ( ! exists($node->{'name'}) || $node->{'name'} ne $parentClass->{'name'} ) ) );
					    last
						unless ( $node );
					    $declaration = &Galacticus::Build::SourceTree::Parse::Declarations::GetDeclaration($node,$objectVariableName)
						if ( &Galacticus::Build::SourceTree::Parse::Declarations::DeclarationExists($node,$objectVariableName) );
					    last
						if ( $declaration );
					    # Move to the parent class.
					    $parentClass = ($parentClass->{'extends'} eq $directive->{'name'}) ? undef() : $classes{$parentClass->{'extends'}};
					}
				    } else {
					# The object being read into is not the object being constructed. We can handle a few special types here.
					if ( &Galacticus::Build::SourceTree::Parse::Declarations::DeclarationExists($node,$objectName) ) {
					    my $declarationTmp = &Galacticus::Build::SourceTree::Parse::Declarations::GetDeclaration($node,$objectName);
					    if ( $declarationTmp->{'intrinsic'} eq "type" ) {
						if      ( $declarationTmp->{'type'} eq "statefulDouble"  ) {
						    $declaration = dclone($declarationTmp);
						    $declaration->{'intrinsic'} = "double precision";
						    $declaration->{'type'     } = undef();
						} elsif ( $declarationTmp->{'type'} eq "statefulInteger" ) {
						    $declaration = dclone($declarationTmp);
						    $declaration->{'intrinsic'} = "integer";
						    $declaration->{'type'     } = undef();
						} elsif ( $declarationTmp->{'type'} eq "statefulLogical" ) {
						    $declaration = dclone($declarationTmp);
						    $declaration->{'intrinsic'} = "logical";
						    $declaration->{'type'     } = undef();
						}
					    }
					}
				    }
				} else {
				    $declaration = &Galacticus::Build::SourceTree::Parse::Declarations::GetDeclaration($node,$variableName)
					if ( &Galacticus::Build::SourceTree::Parse::Declarations::DeclarationExists($node,$variableName) );
				}
				unless ( $declaration ) {
				    if ( exists($constructorNode->{'directive'}->{'type'}) && exists($constructorNode->{'directive'}->{'cardinality'}) ) {
					$declaration->{'parameterType'       } = $constructorNode->{'directive'}->{'type'       };
					$declaration->{'parameterCardinality'} = $constructorNode->{'directive'}->{'cardinality'};
				    } else {
					print "Galacticus::Build::SourceTree::Process::FunctionClass::Process_FunctionClass(): unable to find parameter variable declaration for \"".$constructorNode->{'directive'}->{'name'}."\" in class \"".$className."\"\n";
					die('abort');
				    }
				}
				## Determine the type of the variable.
				my $type;
				if ( exists($declaration->{'parameterType'}) ) {
				    $type = $declaration->{'parameterType'};
				} else {
				    if      ( $declaration->{'intrinsic'} eq "double precision"                                               ) {
					$type = "real";
				    } elsif ( $declaration->{'intrinsic'} eq "integer"                                                        ) {
					$type = "integer";
				    } elsif ( $declaration->{'intrinsic'} eq "logical"                                                        ) {
					$type = "boolean";
				    } elsif ( $declaration->{'intrinsic'} eq "character"                                                      ) {
					$type = "string";
				    } elsif ( $declaration->{'intrinsic'} eq "type"             && $declaration->{'type'} eq "varying_string" ) {
					$type = "string";
				    } else {
					print "Galacticus::Build::SourceTree::Process::FunctionClass::Process_FunctionClass(): unable to determine parameter type - declaration follows:\n";
					print Dumper($declaration);
					die('aborting');
				    }
				}
				## Determine the allowed cardinality.
				my $cardinality;
				if ( exists($declaration->{'parameterCardinality'}) ) {
				    $cardinality = $declaration->{'parameterCardinality'};
				} else {
				    my $cardinalityMinimum;
				    my $cardinalityMaximum;
				    my @dimension = grep {$_ =~ /dimension/} @{$declaration->{'attributes'}};
				    if ( @dimension ) {
					# Non-scalar.
					if ( $dimension[0] =~ m/^dimension\s*\(\s*(.+?)\s*\)/ ) {
					    my @shape = split(/\s*,\s*/,$1);
					    $cardinalityMaximum = 1;
					    for(my $i=0;$i<scalar(@shape);++$i) {
						if ( $shape[$i] =~ m/^\d+$/ ) {
						    $cardinalityMaximum *= $shape[$i];
						} else {
						    $cardinalityMaximum = "*";
						    last;
						}
					    }
					} else {
					    print "Galacticus::Build::SourceTree::Process::FunctionClass::Process_FunctionClass(): unable to parse dimension attribute - declaration follows:\n";
					    print Dumper($declaration);
					    die('aborting');
					}
				    } else {
					# Scalar variable.
					$cardinalityMaximum = 1;
				    }
				    if ( $cardinalityMaximum eq "*" ) {
					# No upper limit.
					if ( exists($constructorNode->{'directive'}->{'defaultValue'}) ) {
					    $cardinality = "0..*"; # Default is present, no upper limit - can have zero or more.
					} else {
					    $cardinality = "1..*"; # Default not present, no upper limit - can have one or more.
					}
				    } else {
					# Upper limit exists.
					if ( exists($constructorNode->{'directive'}->{'defaultValue'}) ) {
					    $cardinality = "0,".$cardinalityMaximum; # Default is present, has upper limit - can have zero or upper limit.
					} else {
					    $cardinality = $cardinalityMaximum; # Default not present, has upper limit - must have precisely the upper limit.
					}
				    }
				}
				my $description =  "\\item[{\\normalfont \\ttfamily [".latex_encode($constructorNode->{'directive'}->{'name'})."]}] ";
				$description .= "(".$type."; ".$cardinality.") ";
				if ( exists($constructorNode->{'directive'}->{'defaultValue'}) ) {
				    my $value = latex_encode($constructorNode->{'directive'}->{'defaultValue'});
				    $value = "true"
					if ( $value eq ".true." );
				    $value = "false"
					if ( $value eq ".false." );
				    if ( $type eq "real" ) {
					$value =~ s/(\d)d([\+\-0-9])/$1e$2/g;
				    }
				    if ( $type eq "integer" ) {
					$value =~ s/\\_c\\_size\\_t//;
				    }
				    if ( $type eq "string" ) {
					$value =~ s/^var\\_str\(['"](.*)['"]\)/$1/;
				    }
				    $description .= " \\{{\\normalfont \\ttfamily ".$value."}".(exists($constructorNode->{'directive'}->{'defaultSource'}) ? "; ".$constructorNode->{'directive'}->{'defaultSource'} : "")."\\} ";
				}
				$description .= $constructorNode->{'directive'}->{'description'};
				push(
				    @parameters,
				    $description	   
				    );
			    }
			    if ( $constructorNode->{'type'} eq "objectBuilder" ) {				
				push(
				    @objects,
				    $constructorNode->{'directive'}->{'class'}
				    );
			    }
			    $constructorNode = &Galacticus::Build::SourceTree::Walk_Tree($constructorNode,\$depth);
			    last
				if ( $depth < 0 );
			}
		    }
		    $node = $node->{'type'} eq "contains" ? $node->{'firstChild'} : $node->{'sibling'};
		}
		$documentationPhysics .= "\n\n\\noindent\\emph{Parameters}\n\\begin{description}\n".join("\n",@parameters)."\n\\end{description}\n"
		    if ( @parameters );
                if ( @objects ) {
		    my @sortedObjects = sort(@objects);
		    $documentationPhysics .= "\n\\noindent\\emph{Classes used}\n\n\\begin{tabular}{ll}\n";
		    for(my $i=0;$i<scalar(@sortedObjects);$i+=2) {
			$documentationPhysics .=    "\\refPhysics{".$sortedObjects[$i  ]."}";
			$documentationPhysics .= " & \\refPhysics{".$sortedObjects[$i+1]."}"
			    if ( $i+1 < scalar(@sortedObjects) );
			$documentationPhysics .= "\\\\\n";
		    }
		    $documentationPhysics .= "\\end{tabular}\n\n";
                }
            }
            (my $descriptiveName = lc($directive->{'descriptiveName'})) =~ s/\s/_/g;
	    system("mkdir -p doc/physics");
	    open(my $docPhysics,">doc/physics/".$descriptiveName.".tex");
	    print $docPhysics $documentationPhysics;
	    close($docPhysics);
	    # Insert into tree.	  
            ## To allow processing of directives by our preprocessor, we parse and process our generated code here.
	    my $treePrecontainsTmp  = &Galacticus::Build::SourceTree::ParseCode($modulePreContains ->{'content'},'Galacticus::Build::SourceTree::Process::FunctionClass::Process_FunctionClass()');
	    my $treePostcontainsTmp = &Galacticus::Build::SourceTree::ParseCode($modulePostContains->{'content'},'Galacticus::Build::SourceTree::Process::FunctionClass::Process_FunctionClass()');
	    &Galacticus::Build::SourceTree::ProcessTree       (                   $treePrecontainsTmp                              );
	    &Galacticus::Build::SourceTree::ProcessTree       (                   $treePostcontainsTmp                             );
	    &Galacticus::Build::SourceTree::InsertAfterNode   ($node            ,[$treePrecontainsTmp                             ]);
	    &Galacticus::Build::SourceTree::InsertPreContains ($node->{'parent'}, $codeContent        ->{'module'}->{'interfaces'} );
	    &Galacticus::Build::SourceTree::InsertPostContains($node->{'parent'},[$treePostcontainsTmp                            ]);	   
	    # Generate submodule files.
	    foreach my $className ( sort(keys(%{$codeContent->{'submodule'}})) ) {
                # Submodule names are just the class name with an underscore appended.
                my $submoduleName = $className."_";
		# Build a file node.
                my $file =
                 {
		     type       => "file" ,
		     parent     => undef(),
		     firstChild => undef(),
		     sibling    => undef(),
		     source     => "Galacticus::Build::SourceTree::Process::FunctionClass::Process_FunctionClass()",
		     line       => 1
		 };
                # The parent (sub)module is either the main module, or the submodule associated with whatever class this class extends.
                my $parentName = $node->{'parent'}->{'name'}.($classes{$className}->{'extends'} eq $directive->{'name'}."Class" ? "" : ":".$classes{$className}->{'extends'}."_");
                # Build a submodule node.
                my $submodule =
                 {
		     type       => "submodule"                                       ,
		     opener     => "submodule (".$parentName.") ".$submoduleName."\n",
		     closer     => "end submodule "              .$submoduleName."\n",
		     parent     => undef(),
		     firstChild => undef(),
		     sibling    => undef(),
		     source     => "Galacticus::Build::SourceTree::Process::FunctionClass::Process_FunctionClass()",
		     line       => 1
		 };
		&Galacticus::Build::SourceTree::PrependChildToNode ($file,[$submodule]);
		# Remove existing links from all submodule nodes.
		foreach ( @{$codeContent->{'submodule'}->{$className}->{'preContains' }}, @{$codeContent->{'submodule'}->{$className}->{'postContains' }} ) {
		    $_->{'parent' } = undef();
		    $_->{'sibling'} = undef();
		}
		# Insert the pre- and post-contains content into the submodule.
		&Galacticus::Build::SourceTree::PrependChildToNode($submodule,$codeContent->{'submodule'}->{$className}->{'preContains' });
		&Galacticus::Build::SourceTree::InsertPostContains($submodule,$codeContent->{'submodule'}->{$className}->{'postContains'});
		# Write the submodule to a temporary file, and update the actual file only if it has changed (to avoid recompilation cascades).
		(my $submoduleContent, my $submoduleMappings) = &Galacticus::Build::SourceTree::Serialize($file, stripMappings => 1);
		open(my $submoduleFile,">",$codeContent->{'submodule'}->{$className}->{'fileName'}.".tmp");
		print $submoduleFile $submoduleContent;
		close($submoduleFile);
		&File::Changes::Update($codeContent->{'submodule'}->{$className}->{'fileName'},$codeContent->{'submodule'}->{$className}->{'fileName'}.".tmp", proveUpdate => "yes");
		open(my $mappingFile,">",$codeContent->{'submodule'}->{$className}->{'fileName'}.".lmap");
		print $mappingFile $submoduleMappings;
		close($mappingFile);
	       }
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

sub deepCopyLinkedList {
    # Create deep-copy instructions for linked list objects.
    my $class                       = shift();
    my $nonAbstractClass            = shift();
    my $linkedListVariables         = shift();
    my $linkedListResetVariables    = shift();
    my $linkedListFinalizeVariables = shift();
    my $debugging                   = shift();
    return ("","","",undef())
	unless ( exists($class->{'linkedList'}) );
    my $linkedList = $class->{'linkedList'};
    # Get object names and types.
    my @objects     = split(" ",$linkedList->{'object'    });
    my @objectTypes = split(" ",$linkedList->{'objectType'});
    # Add variables needed for linked list processing.
    push(
	@{$linkedListVariables},
	{
	    intrinsic  => 'type',
	    type       => $linkedList->{'type'},
	    attributes => [ 'pointer' ],
	    variables  => [ $linkedList->{'type'}.'item', $linkedList->{'type'}.'destination', $linkedList->{'type'}.'itemNew' ]
	}
	)
	unless ( grep {$_->{'type'} eq $linkedList->{'type'}} @{$linkedListVariables} );
    push(
	@{$linkedListResetVariables},
	{
	    intrinsic  => 'type',
	    type       => $linkedList->{'type'},
	    attributes => [ 'pointer' ],
	    variables  => [ $linkedList->{'type'}.'item' ]
	}
	)
	unless ( grep {$_->{'type'} eq $linkedList->{'type'}} @{$linkedListResetVariables} );
    push(
	@{$linkedListFinalizeVariables},
	{
	    intrinsic  => 'type',
	    type       => $linkedList->{'type'},
	    attributes => [ 'pointer' ],
	    variables  => [ $linkedList->{'type'}.'item' ]
	}
	)
	unless ( grep {$_->{'type'} eq $linkedList->{'type'}} @{$linkedListFinalizeVariables} );
    # Generate code for the walk through the linked list.
    my $deepCopyCode;
    my $deepCopyResetCode;
    my $deepCopyFinalizeCode;
    for(my $i=0;$i<scalar(@objects);++$i) {
	$code::type            =                                            $linkedList->{'type'    };
	$code::variable        =                                            $linkedList->{'variable'};
	$code::next            =                                            $linkedList->{'next'    };
	$code::object          =                                            $objects    [$i]         ;
	$code::objectType      =                                            $objectTypes[$i]         ;
	$code::location        = &Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($class->{'node'},$class->{'node'}->{'line'});
	$code::debugCode       = $debugging ? "if (debugReporting.and.mpiSelf\%isMaster()) call displayMessage(var_str('functionClass[own] (class : ownerName : ownerLoc : objectLoc : sourceLoc): [".$code::objectType."] : ".$code::object." : ')//loc(".$code::type."itemNew)//' : '//loc(".$code::type."itemNew%".$code::object.")//' : '//".&Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($class->{'node'},$class->{'node'}->{'line'},compact => 1).",verbosityLevelSilent)\n" : "";
	if ( $i == 0 ) {
	    $deepCopyCode .= fill_in_string(<<'CODE', PACKAGE => 'code');
destination%{$variable} => null            ()
CODE
	}
	$deepCopyCode .= fill_in_string(<<'CODE', PACKAGE => 'code');
{$type}destination      => null            ()
{$type}item             => self%{$variable}
do while (associated({$type}item))
CODE
	if ( $i == 0 ) {
	    $deepCopyCode .= fill_in_string(<<'CODE', PACKAGE => 'code');
   allocate({$type}itemNew)
   if (associated({$type}destination)) then
      {$type}destination%{$next}     => {$type}itemNew
      {$type}destination             => {$type}itemNew
   else
      destination       %{$variable} => {$type}itemNew
      {$type}destination             => {$type}itemNew
   end if
CODE
	} else {
	    $deepCopyCode .= fill_in_string(<<'CODE', PACKAGE => 'code');
   if (associated({$type}destination)) then
      {$type}itemNew     => {$type}destination%{$next}
      {$type}destination => {$type}destination%{$next}
   else
      {$type}itemNew     => destination       %{$variable}
      {$type}destination => destination       %{$variable}
   end if
CODE
	}
	$deepCopyCode .= fill_in_string(<<'CODE', PACKAGE => 'code');
      nullify({$type}itemNew%{$object})
      if (associated({$type}item%{$object})) then
       if (associated({$type}item%{$object}%copiedSelf)) then
        select type(s => {$type}item%{$object}%copiedSelf)
        class is ({$objectType})
         {$type}itemNew%{$object} => s
        class default
         call Error_Report('copiedSelf has incorrect type'//{$location})
        end select
        call {$type}item%{$object}%copiedSelf%referenceCountIncrement()
       else
        allocate({$type}itemNew%{$object},mold={$type}item%{$object})
        call {$type}item%{$object}%deepCopy({$type}itemNew%{$object})
        {$type}item%{$object}%copiedSelf => {$type}itemNew%{$object}
        call {$type}itemNew%{$object}%autoHook()
       end if
       {$debugCode}
      end if
   {$type}item => {$type}item%{$next}
end do
CODE
	$deepCopyResetCode .= fill_in_string(<<'CODE', PACKAGE => 'code');
{$type}item => self%{$variable}
do while (associated({$type}item))
   call {$type}item%{$object}%deepCopyReset()
   {$type}item => {$type}item%{$next}
end do
CODE
	$deepCopyFinalizeCode .= fill_in_string(<<'CODE', PACKAGE => 'code');
{$type}item => self%{$variable}
do while (associated({$type}item))
   call {$type}item%{$object}%deepCopyFinalize()
   {$type}item => {$type}item%{$next}
end do
CODE
    }
    my $deepCopyModule = exists($linkedList->{'module'}) ? $linkedList->{'module'} : undef();
    return ($deepCopyCode,$deepCopyResetCode,$deepCopyFinalizeCode,$deepCopyModule);
}

sub stateStoreLinkedList {
    # Create state store/restore instructions for linked list objects.
    my $class               = shift();
    my $nonAbstractClass    = shift();
    my $linkedListVariables = shift();
    return ("","","",undef())
	unless ( exists($class->{'linkedList'}) );
    my $linkedList = $class->{'linkedList'};
    # Get object names.
    my @objects = split(" ",$linkedList->{'object'});
    # Add variables needed for linked list processing.
    push(
	@{$linkedListVariables},
	{
	    intrinsic  => 'type',
	    type       => $linkedList->{'type'},
	    attributes => [ 'pointer' ],
	    variables  => [ $linkedList->{'type'}.'item' ]
	}
	)
	unless ( grep {$_->{'type'} eq $linkedList->{'type'}} @{$linkedListVariables} );
    # Generate code for the walk through the linked list.
    my $inputCode;
    my $outputCode;
    for(my $i=0;$i<scalar(@objects);++$i) {
	$code::type     = $linkedList->{'type'    };
	$code::variable = $linkedList->{'variable'};
	$code::next     = $linkedList->{'next'    };
	$code::object   = $objects[$i];
	$inputCode .= fill_in_string(<<'CODE', PACKAGE => 'code');
{$type}item => self%{$variable}
do while (associated({$type}item))
   call {$type}item%{$object}%stateRestore(stateFile,gslStateFile,stateOperationID)
   {$type}item => {$type}item%{$next}
end do
CODE
	$outputCode .= fill_in_string(<<'CODE', PACKAGE => 'code');
{$type}item => self%{$variable}
do while (associated({$type}item))
   call {$type}item%{$object}%stateStore(stateFile,gslStateFile,stateOperationID)
   {$type}item => {$type}item%{$next}
end do
CODE
    }
    my $deepCopyModule = exists($linkedList->{'module'}) ? $linkedList->{'module'} : undef();
    return ($inputCode,$outputCode,$deepCopyModule);
}

sub allowedParametersLinkedList {
    # Create allowed parameter instructions for linked list objects.
    my $nonAbstractClass    = shift();
    my $linkedListVariables = shift();
    my $source              = shift();
    return ("",undef())
	unless ( exists($nonAbstractClass->{'linkedList'}) );
    my $linkedList = $nonAbstractClass->{'linkedList'};
    # Get object names.
    my @objects = split(" ",$linkedList->{'object'});
    # Add variables needed for linked list processing.
    push(
	@{$linkedListVariables},
	{
	    intrinsic  => 'type',
	    type       => $linkedList->{'type'},
	    attributes => [ 'pointer' ],
	    variables  => [ $linkedList->{'type'}.'item' ]
	}
	)
	unless ( grep {$_->{'type'} eq $linkedList->{'type'}} @{$linkedListVariables} );
    # Generate code for the walk through the linked list.
    my $iterator;
    for(my $i=0;$i<scalar(@objects);++$i) {
	$code::type     = $linkedList->{'type'    };
	$code::variable = $linkedList->{'variable'};
	$code::next     = $linkedList->{'next'    };
	$code::object   = $objects[$i];
	$code::source   = $source;
	$iterator .= fill_in_string(<<'CODE', PACKAGE => 'code');
{$type}item => self%{$variable}
do while (associated({$type}item))
   call {$type}item%{$object}%allowedParameters(allowedParameters,'{$source}',.true.)
   {$type}item => {$type}item%{$next}
end do
CODE
    }
    my $deepCopyModule = exists($linkedList->{'module'}) ? $linkedList->{'module'} : undef();
    return ($iterator,$deepCopyModule);
}

sub autoDescriptorLinkedList {
    # Create auto-descriptor instructions for linked list objects.
    my $linkedList          = shift();
    my $linkedListVariables = shift();
    # Get object names.
    my @objects = split(" ",$linkedList->{'object'});
    # Add variables needed for linked list processing.
    push(
	@{$linkedListVariables},
	{
	    intrinsic  => 'type',
	    type       => $linkedList->{'type'},
	    attributes => [ 'pointer' ],
	    variables  => [ $linkedList->{'type'}.'item' ]
	}
	)
	unless ( grep {$_->{'type'} eq $linkedList->{'type'}} @{$linkedListVariables} );
    # Generate code for the walk through the linked list.
    my $iterator;
    for(my $i=0;$i<scalar(@objects);++$i) {
	$code::type     = $linkedList->{'type'    };
	$code::variable = $linkedList->{'variable'};
	$code::next     = $linkedList->{'next'    };
	$code::object   = $objects[$i];
	$iterator .= fill_in_string(<<'CODE', PACKAGE => 'code');
{$type}item => self%{$variable}
do while (associated({$type}item))
   call {$type}item%{$object}%descriptor(parameters)
   {$type}item => {$type}item%{$next}
end do
CODE
    }
    my $deepCopyModule = exists($linkedList->{'module'}) ? $linkedList->{'module'} : undef();
    return ($iterator,$deepCopyModule);
}

sub stateStoreExplicitFunction {
    # Create state store/restore instructions for objects with explicit functions.
    my $nonAbstractClass  = shift();
    my $inputCode  = "";
    my $outputCode = "";
    my %modules;
    # Handle store.
    if ( exists($nonAbstractClass->{'stateStore'}->{'stateStore'}->{'variables'}) ) {
	foreach my $explicitVariable ( split(" ",$nonAbstractClass->{'stateStore'}->{'stateStore'}->{'variables'}) ) {
	    if ( exists($nonAbstractClass->{'stateStore'}->{'stateStore'}->{'store'  }) ) {
		$outputCode .= "if (associated(self%".$explicitVariable.")) then\n";
		$outputCode .= " write (stateFile) .true.\n";
		$outputCode .= " call ".$nonAbstractClass->{'stateStore'}->{'stateStore'}->{'store'  }."(self%".$explicitVariable.",stateFile,gslStateFile,stateOperationID)\n";
		$outputCode .= "else\n";
		$outputCode .= " write (stateFile) .false.\n";
		$outputCode .= "end if\n";
	    }
	    if ( exists($nonAbstractClass->{'stateStore'}->{'stateStore'}->{'restore'}) ) {
		$inputCode  .= "read (stateFile) wasAssociated\n";
		$inputCode  .= "if (wasAssociated) then\n";
		$inputCode  .= " call ".$nonAbstractClass->{'stateStore'}->{'stateStore'}->{'restore'}."(self%".$explicitVariable.",stateFile,gslStateFile,stateOperationID)\n";
		$inputCode  .= "else\n";
		$inputCode  .= " nullify(self%".$explicitVariable.")\n";
		$inputCode  .= "end if\n";
	    }
	    $modules{$nonAbstractClass->{'stateStore'}->{'stateStore'}->{'module'}} = 1
		if ( exists($nonAbstractClass->{'stateStore'}->{'stateStore'}->{'module' }) );   
	}
    }
    return ($inputCode,$outputCode,%modules);
}

sub potentialDescriptorParameters {
    # Process variable declarations for potential parameters to include in descriptors.
    my $declarations   = shift();
    my $class          = shift();
    my $potentialNames = shift();
    our $stateStorables;
    foreach my $declaration ( &List::ExtraUtils::as_array($declarations) ) {
	# Identify object pointers.
	push(@{$potentialNames->{'objects'}},map {$_ =~ s/\s*([a-zA-Z0-9_]+).*/$1/; $_} @{$declaration->{'variables'}})
	    if
	    (
	     $declaration->{'intrinsic'} eq "class"
	     &&
	     (
	      (grep {&lctrim($declaration->{'type'}) eq lc($_)} keys                       (%{$stateStorables->{'functionClasses'       }}))
	      ||
	      (grep {&lctrim($declaration->{'type'}) eq lc($_)} &List::ExtraUtils::as_array(  $stateStorables->{'functionClassInstances'} ))
	     )
	     &&
	     grep {$_ eq "pointer"} @{$declaration->{'attributes'}}
	    );
	# Identify stateful types.
	push(@{$potentialNames->{'statefulTypes'}},$declaration)
	    if
	    (
	     $declaration->{'intrinsic'} eq "type"
	     &&
	     $declaration->{'type'     } =~ m/^stateful(Integer|Double|Logical)\s*$/i
	    );
	# Identify enumerations.
	push(@{$potentialNames->{'enumerations'}},$declaration)
	    if
	    (
	     $declaration->{'intrinsic'} eq "type"
	     &&
	     $declaration->{'type'     } =~ m/^enumeration[a-z0-9_]+type\s*$/i
	    );
	# Identify regular parameters.
	push(@{$potentialNames->{'parameters'}},$declaration)
	    if
	    (
	     (grep {$_ eq $declaration->{'intrinsic'}} ( "integer", "logical", "double precision", "character" ))
	     ||
	     (
	      $declaration->{'intrinsic'}  eq "type"
	      &&
	      trimlc($declaration->{'type'     }) eq "varying_string"
	     )
	    );
	$class->{'hasCustomDescriptor'} = 1
	    if
	    (
	     $declaration->{'intrinsic'} eq "procedure"
	     &&
	     $declaration->{'variables'}->[0] =~ m/^descriptor=>/
	    );
    }
}
    
sub deepCopyDeclarations {
    # Process variable declarations from a node for deep copy.
    my $class              =   shift() ;
    my $nonAbstractClass   =   shift() ;
    my $node               =   shift() ;
    my $declarations       =   shift() ;
    my @ignore             = @{shift()};
    my $lineNumber         =   shift() ;
    my $deepCopy           =   shift() ;
    my $foundDeepCopyNames =   shift() ;
    our $stateStorables;
    our $debugging;
    our $deepCopyActions;
    foreach my $declaration ( &List::ExtraUtils::as_array($declarations) ) {
	# Deep copy of functionClass objects.
	(my $type = $declaration->{'type'}) =~ s/(^\s*|\s*$)//g
	    if ( $declaration->{'intrinsic'} eq "class" || $declaration->{'intrinsic'} eq "type" );
	if
	    (
	     $declaration->{'intrinsic'} eq "class"
	     &&
	     (grep {$_ eq $type    } (keys(%{$stateStorables->{'functionClasses'}}),@{$stateStorables->{'functionClassInstances'}}))
	     &&
	     grep {$_ eq "pointer"}  @{$declaration   ->{'attributes'     }}
	    )
	{
	    foreach my $object ( @{$declaration->{'variables'}} ) {
		(my $name = $object) =~ s/^([a-zA-Z0-9_]+).*/$1/; # Strip away anything (e.g. assignment operators) after the variable name.
		next
		    if ( grep {lc($_) eq lc($name)} @ignore );
		$deepCopy->{'resetCode'   } .= "if (associated(self%".$name.")) call self%".$name."%deepCopyReset   ()\n";
		$deepCopy->{'finalizeCode'} .= "if (associated(self%".$name.")) call self%".$name."%deepCopyFinalize()\n";
		$deepCopy->{'assignments' } .= "nullify(destination%".$name.")\n";
		$deepCopy->{'assignments' } .= "if (associated(self%".$name.")) then\n";
		$deepCopy->{'assignments' } .= " if (associated(self%".$name."\%copiedSelf)) then\n";
		$deepCopy->{'assignments' } .= "  select type(s => self%".$name."\%copiedSelf)\n";
		$deepCopy->{'assignments' } .= "  ".$declaration->{'intrinsic'}." is (".$declaration->{'type'}.")\n";
		$deepCopy->{'assignments' } .= "   destination%".$name." => s\n";
		$deepCopy->{'assignments' } .= "  class default\n";
		$deepCopy->{'assignments' } .= "   call Error_Report('copiedSelf has incorrect type'//".&Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($nonAbstractClass->{'node'},$nonAbstractClass->{'node'}->{'line'}).")\n";
		$deepCopy->{'assignments' } .= "  end select\n";
		$deepCopy->{'assignments' } .= "  call self%".$name."\%copiedSelf\%referenceCountIncrement()\n";
		$deepCopy->{'assignments' } .= " else\n";
		$deepCopy->{'assignments' } .= "  allocate(destination%".$name.",mold=self%".$name.")\n";
		$deepCopy->{'assignments' } .= "  call self%".$name."%deepCopy(destination%".$name.")\n";
		$deepCopy->{'assignments' } .= "  self%".$name."%copiedSelf => destination%".$name."\n";
		$deepCopy->{'assignments' } .= "  call destination%".$name."%autoHook()\n";
		$deepCopy->{'assignments' } .= " end if\n";
		$deepCopy->{'assignments' } .= " if (debugReporting.and.mpiSelf\%isMaster()) call displayMessage(var_str('functionClass[own] (class : ownerName : ownerLoc : objectLoc : sourceLoc): ".$name." : [destination] : ')//loc(destination)//' : '//loc(destination%".$name.")//' : '//".&Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($node,$lineNumber,compact => 1).",verbosityLevelSilent)\n"
		    if ( $debugging );
		$deepCopy->{'assignments' } .= "end if\n";
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
	     (grep {$_->{'type'} eq $type} &List::ExtraUtils::as_array($deepCopyActions->{'deepCopyActions'}))
	    ) {
		my $isAllocatable = grep {$_ eq "allocatable"} @{$declaration->{'attributes'}};
		my $isPointer     = grep {$_ eq "pointer"    } @{$declaration->{'attributes'}};
		my $rank = 0;
		if ( grep {$_ =~ m/^dimension\s*\(/} @{$declaration->{'attributes'}} ) {
		    my $dimensionDeclarator = join(",",map {/^dimension\s*\(([a-zA-Z0-9_,:\s]+)\)/} @{$declaration->{'attributes'}});
		    $rank        = ($dimensionDeclarator =~ tr/,//)+1;
		    $deepCopy->{'rankMaximum'} = $rank
			if ( $rank > $deepCopy->{'rankMaximum'} );					    
		}
		foreach my $variableName ( @{$declaration->{'variableNames'}} ) {
		    $deepCopy->{'assignments'} .= "if (allocated(self%".$variableName.")) then\n"
			if ( $isAllocatable );
		    $deepCopy->{'assignments'} .= "if (associated(self%".$variableName.")) then\n"
			if ( $isPointer     );
		    for(my $i=1;$i<=$rank;++$i) {
			$deepCopy->{'assignments'} .= (" " x $i)."do i".$i."=lbound(self%".$variableName.",dim=".$i."),ubound(self%".$variableName.",dim=".$i.")\n";
		    }
		    my $arrayElement = $rank > 0 ? "(".join(",",map {"i".$_} 1..$rank).")" : "";
		    $deepCopy->{'assignments'} .= (" " x $rank)."call destination%".$variableName.$arrayElement."%deepCopyActions()\n";
		    for(my $i=1;$i<=$rank;++$i) {
			$deepCopy->{'assignments'} .= (" " x ($rank+1-$i))."end do\n";
		    }
		    $deepCopy->{'assignments'} .= "end if\n"
			if ( $isAllocatable || $isPointer );
		}
	}
	# Deep copy of HDF5 objects.
	if
	    (
	     $declaration->{'intrinsic'} eq "type"
	     &&
	     $declaration->{'type'     } =~ m/^\s*hdf5object\s*$/i
	    ) {
		$deepCopy->{'modules'}->{'HDF5_Access'} = 1;
		$deepCopy->{'assignments'} .= "!\$ call hdf5Access%set  ()\n";
		$deepCopy->{'assignments'} .= "call self%".$_."%deepCopy(destination%".$_.")\n"
		    foreach ( @{$declaration->{'variables'}} );
		$deepCopy->{'assignments'} .= "!\$ call hdf5Access%unset()\n";
	}
	# Deep copy of non-(class,pointer) functionClass objects.
	if ( exists($class->{'deepCopy'}->{'functionClass'}) ) {
	    foreach my $object ( @{$declaration->{'variables'}} ) {
		(my $name = $object) =~ s/^([a-zA-Z0-9_]+).*/$1/; # Strip away anything (e.g. assignment operators) after the variable name.
		if ( grep {lc($_) eq lc($name)} split(/\s*,\s*/,$class->{'deepCopy'}->{'functionClass'}->{'variables'}) ) {
		    push(@{$foundDeepCopyNames},$name);
		    if ( grep {$_ eq "pointer"}  @{$declaration->{'attributes'}} ) {
			$deepCopy->{'assignments' } .= "nullify(destination%".$name.")\n";
			$deepCopy->{'assignments' } .= "if (associated(self%".$name.")) then\n";
			$deepCopy->{'resetCode'   } .= "if (associated(self%".$name.")) then\n";
			$deepCopy->{'finalizeCode'} .= "if (associated(self%".$name.")) then\n";
			$deepCopy->{'assignments' } .= "if (associated(self%".$name."\%copiedSelf)) then\n";
			$deepCopy->{'assignments' } .= "  select type(s => self%".$name."\%copiedSelf)\n";
			$deepCopy->{'assignments' } .= "  ".$declaration->{'intrinsic'}." is (".$declaration->{'type'}.")\n";
			$deepCopy->{'assignments' } .= "   destination%".$name." => s\n";
			$deepCopy->{'assignments' } .= "  class default\n";
			$deepCopy->{'assignments' } .= "   call Error_Report('copiedSelf has incorrect type'//".&Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($nonAbstractClass->{'node'},$nonAbstractClass->{'node'}->{'line'}).")\n";
			$deepCopy->{'assignments' } .= "  end select\n";
			$deepCopy->{'assignments' } .= "  call self%".$name."\%copiedSelf\%referenceCountIncrement()\n";
			$deepCopy->{'assignments' } .= "else\n";
			$deepCopy->{'assignments' } .= " allocate(destination%".$name.",mold=self%".$name.")\n";
		    }
		    $deepCopy->{'resetCode'   } .= "call self%".$name."%deepCopyReset   ()\n";
		    $deepCopy->{'finalizeCode'} .= "call self%".$name."%deepCopyFinalize()\n";
		    $deepCopy->{'assignments' } .= "call self%".$name."%deepCopy(destination%".$name.")\n";
		    $deepCopy->{'assignments' } .= "self%".$name."%copiedSelf => destination%".$name."\n";
		    $deepCopy->{'assignments' } .= "call destination%".$name."%autoHook()\n";
		    $deepCopy->{'assignments' } .= "end if\n"
			if ( grep {$_ eq "pointer"}  @{$declaration->{'attributes'}} );
		    $deepCopy->{'assignments' } .= "if (debugReporting.and.mpiSelf\%isMaster()) call displayMessage(var_str('functionClass[own] (class : ownerName : ownerLoc : objectLoc : sourceLoc): ".$name." : [destination] : ')//loc(destination)//' : '//loc(destination%".$name.")//' : '//".&Galacticus::Build::SourceTree::Process::SourceIntrospection::Location($node,$lineNumber,compact => 1).",verbosityLevelSilent)\n"
			if ( $debugging );
		    if ( grep {$_ eq "pointer"}  @{$declaration->{'attributes'}} ) {
			$deepCopy->{'assignments' } .= "end if\n";
			$deepCopy->{'resetCode'   } .= "end if\n";
			$deepCopy->{'finalizeCode'} .= "end if\n";
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
			$deepCopy->{'assignments'} .= "!\$omp atomic\n"
			    if ( exists($class->{'deepCopy'}->{'increment'}->{'atomic'}) && $class->{'deepCopy'}->{'increment'}->{'atomic'} eq "yes" );
			$deepCopy->{'assignments'} .= "destination\%".$increment->{'variable'}."=destination\%".$increment->{'variable'}."+1\n";
		    }
		}
	    }
	}
	# Perform any sets.
	if ( exists($class->{'deepCopy'}->{'setTo'}) ) {
	    my @setTos = map {{variable => $_}} split(/\s*,\s*/,$class->{'deepCopy'}->{'setTo'}->{'variables'});
	    foreach ( @setTos ) {
		($_->{'host'} = $_->{'variable'}) =~ s/^([^%]+)%.+/$1/;
	    }
	    foreach my $object ( @{$declaration->{'variables'}} ) {
		(my $name = $object) =~ s/^([a-zA-Z0-9_]+).*/$1/; # Strip away anything (e.g. assignment operators) after the variable name.
		foreach my $setTo ( @setTos ) {
		    if ( lc($setTo->{'host'}) eq lc($name) ) {
			$deepCopy->{'assignments'} .= "destination\%".$setTo->{'variable'}."=".$class->{'deepCopy'}->{'setTo'}->{'value'}."\n";
		    }
		}
	    }
	}
	# Perform any explicit deep copies.
	if ( exists($class->{'deepCopy'}->{'deepCopy'}) ) {
	    my @deepCopies = split(/\s*,\s*/,$class->{'deepCopy'}->{'deepCopy'}->{'variables'});
	    foreach my $object ( @{$declaration->{'variableNames'}} ) {
		foreach my $deepCopy ( @deepCopies ) {
		    if ( lc($object) eq lc($deepCopy) ) {
			$deepCopy->{'assignments'} .= "nullify(destination\%".$object.")\n";
			$deepCopy->{'assignments'} .= "allocate(destination\%".$object.",mold=self\%".$object.")\n";
			if ( exists($class->{'deepCopy'}->{'deepCopy'}->{'copy'}) ) {
			    $deepCopy->{'modules'}->{$class->{'deepCopy'}->{'deepCopy'}->{'module'}} = 1
				if ( exists($class->{'deepCopy'}->{'deepCopy'}->{'module'}) );
			    $deepCopy->{'assignments'} .= "if (associated(self\%".$object.")) call ".$class->{'deepCopy'}->{'deepCopy'}->{'copy'}."(self\%".$object.",destination\%".$object.")\n";
			} else {
			    $deepCopy->{'assignments'} .= "if (associated(self\%".$object.")) call self\%".$object."\%deepCopy(destination\%".$object.")\n";
			}
			if ( exists($class->{'deepCopy'}->{'deepCopy'}->{'reset'}) ) {
			    $deepCopy->{'resetModules'}->{$class->{'deepCopy'}->{'deepCopy'}->{'module'}} = 1
				if ( exists($class->{'deepCopy'}->{'deepCopy'}->{'module'}) );
			    $deepCopy->{'resetCode'} .= "if (associated(self\%".$object.")) call ".$class->{'deepCopy'}->{'deepCopy'}->{'reset'}."(self\%".$object.")\n";
			}
			if ( exists($class->{'deepCopy'}->{'deepCopy'}->{'finalize'}) ) {
			    $deepCopy->{'finalizeModules'}->{$class->{'deepCopy'}->{'deepCopy'}->{'module'}} = 1
				if ( exists($class->{'deepCopy'}->{'deepCopy'}->{'module'}) );
			    $deepCopy->{'finalizeCode'} .= "if (associated(self\%".$object.")) call ".$class->{'deepCopy'}->{'deepCopy'}->{'finalize'}."(self\%".$object.")\n";
			}
		    }
		}
	    }
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
		$deepCopy->{'assignments'} .= "!\$ call OMP_Init_Lock(destination\%".$_.")\n"
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
			$deepCopy->{'rankMaximum'} = scalar(@rank)
			    if ( scalar(@rank) > $deepCopy->{'rankMaximum'} );
			for(my $i=1;$i<=scalar(@rank);++$i) {
			    $deepCopy->{'assignments'} .= "!\$ do i".$i."=lbound(destination\%".$_.",dim=".$i."),ubound(destination\%".$_.",dim=".$i.")\n";
			}
			$deepCopy->{'assignments'} .= "!\$    call destination\%".$_."(".join(",",map {"i".$_} 1..scalar(@rank)).")%initialize()\n";
			for(my $i=1;$i<=scalar(@rank);++$i) {
			    $deepCopy->{'assignments'} .= "!\$ end do\n";
			}
		    } else {
			# Scalar lock.
			$deepCopy->{'assignments'} .= "!\$ call destination\%".$_."%initialize()\n";
		    }
		}
	}
    }
}

sub stateStoreVariables {
    # Generate code to store/restore variables in functionClass objects.
    my  $stateStores        = shift();
    my  $stateStore         = shift();
    my  $class              = shift();
    my  $declarations       = shift();
    my  $explicitNamesFound = shift();
    our $stateStorables              ;
    foreach my $declaration ( &List::ExtraUtils::as_array($declarations) ) {
	# Identify variable type.
	if ( $declaration->{'intrinsic'} eq "procedure" || $declaration->{'intrinsic'} eq "final" ) {
	    # Type-bound procedure - nothing to do.
	} elsif ( $declaration->{'intrinsic'} eq "class" || $declaration->{'intrinsic'} eq "type" ) {
	    # Look for pointers to functionClasses.
	    (my $type = $declaration->{'type'}) =~ s/\s//g;
	    if (
		$declaration->{'intrinsic'} eq "class"
		&&
		(grep {$_ eq "pointer"}      @{$declaration   ->{'attributes'     }} )
		&&
		(grep {$_ eq $type    } keys(%{$stateStorables->{'functionClasses'}}))
		) {
		# Pointer to a functionClass object.
		foreach ( @{$declaration->{'variables'}} ) {
		    $stateStores->{'labelUsed'} = 1;
		    (my $variableName = $_) =~ s/\s*=.*$//;
		    next
			if ( grep {lc($_) eq lc($variableName)} @{$stateStore->{'excludes'}} );
		    $stateStore->{'outputCode'} .= " if (displayVerbosity() >= verbosityLevelWorking) then\n";
		    $stateStore->{'outputCode'} .= "  select type (c__ => self%".$variableName.")\n";
		    $stateStore->{'outputCode'} .= "  class is (".$declaration->{'type'}.")\n";
		    # <workaround type="gfortran" PR="94446" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=94446">
		    #  <description>
		    #   Using the sizeof() intrinsic on a treeNode object causes a bogus "type mismatch" error when this module is used.
		    #  </description>
		    # </workaround>
		    $stateStore->{'outputCode'} .= "   write (label,'(i16)') 0\n";
		    #$stateStore->{'outputCode'} .= "   write (label,'(i16)') sizeof(c__)\n";
		    $stateStore->{'outputCode'} .= "  end select\n";
		    $stateStore->{'outputCode'} .= "  call displayMessage('storing \"".$variableName."\" with size '//trim(adjustl(label))//' bytes')\n";
		    $stateStore->{'outputCode'} .= " end if\n";
		    $stateStore->{'inputCode'}  .= " call displayMessage('restoring \"".$variableName."\"',verbosity=verbosityLevelWorking)\n";
		    $stateStore->{'outputCode'} .= " call self%".$variableName."%stateStore  (stateFile,gslStateFile,stateOperationID)\n";
		    $stateStore->{'inputCode'}  .= " call self%".$variableName."%stateRestore(stateFile,gslStateFile,stateOperationID)\n";
		    $stateStores->{'stateFileUsed'}    = 1;
		    $stateStores->{'gslStateFileUsed'} = 1;
		}
	    } elsif (
		$declaration->{'intrinsic'} eq "type"
		&&
		$declaration->{'type'     } =~ m/^enumeration[a-z0-9_]+type/i
		) {
		# Enumeration.
		if ( grep {$_ eq "allocatable"} @{$declaration->{'attributes'}} ) {
		    # For allocatable variables we must first store the shape so that they can be reallocated on restore.
		    my $dimensionDeclarator = join(",",map {/^dimension\s*\(([:,]+)\)/} @{$declaration->{'attributes'}});
		    my $rank = ($dimensionDeclarator =~ tr/://);
		    foreach my $variableName ( @{$declaration->{'variables'}} ) {
			next
			    if ( grep {lc($_) eq lc($variableName)} @{$stateStore->{'excludes'}} );
			$stateStores->{'allocatablesFound'}  = 1;
			$stateStores->{'dimensionalsFound'}  = 1;
			$stateStores->{'stateFileUsed'}      = 1;
			$stateStores->{'labelUsed'}          = 1;
			$stateStore->{'outputCode'}        .= " if (allocated(self%".$variableName.")) then\n";
			$stateStore->{'outputCode'}        .= "  if (displayVerbosity() >= verbosityLevelWorking) then\n";
			# <workaround type="gfortran" PR="94446" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=94446">
			#  <description>
			#   Using the sizeof() intrinsic on a treeNode object causes a bogus "type mismatch" error when this module is used.
			#  </description>
			# </workaround>
			$stateStore->{'outputCode'} .= "   write (label,'(i16)') 0\n";
			#$stateStore->{'outputCode'}        .= "   write (label,'(i16)') sizeof(self%".$variableName.")\n";
			$stateStore->{'outputCode'}        .= "   call displayMessage('storing \"".$variableName."\" with size '//trim(adjustl(label))//' bytes')\n";
			$stateStore->{'outputCode'}        .= "  end if\n";
			$stateStore->{'outputCode'}        .= "  write (stateFile) .true.\n"
			    . "  write (stateFile) shape(self%".$variableName.",kind=c_size_t)\n"
			    . "  write (stateFile) self%".$variableName."%ID\n";
			$stateStore->{'outputCode'}        .= " else\n";
			$stateStore->{'outputCode'}        .= "  write (stateFile) .false.\n";
			$stateStore->{'outputCode'}        .= " end if\n";
			$stateStore->{'inputCode'}         .= " read (stateFile) wasAllocated\n";
			$stateStore->{'inputCode'}         .= " if (allocated(self%".$variableName.")) deallocate(self%".$variableName.")\n";
			$stateStore->{'inputCode'}         .= " if (wasAllocated) then\n";
			$stateStore->{'inputCode'}         .= "  call displayMessage('restoring \"".$variableName."\"',verbosity=verbosityLevelWorking)\n";
			$stateStore->{'inputCode'}         .= "  allocate(storedShape(".$rank."))\n";
			$stateStore->{'inputCode'}         .= "  read (stateFile) storedShape\n";
			$stateStore->{'inputCode'}         .= "  allocate(self%".$variableName."(".join(",",map {"storedShape(".$_.")"} 1..$rank)."))\n";
			$stateStore->{'inputCode'}         .= "  deallocate(storedShape)\n";
			$stateStore->{'inputCode'}         .= "  read (stateFile) self%".$variableName."%ID\n";
			$stateStore->{'inputCode'}         .= " end if\n";
		    }
		} else {
		    $stateStore->{'outputCode'} .= " if (displayVerbosity() >= verbosityLevelWorking) then\n";
		    foreach ( @{$declaration->{'variableNames'}} ) {
			# <workaround type="gfortran" PR="94446" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=94446">
			#  <description>
			#   Using the sizeof() intrinsic on a treeNode object causes a bogus "type mismatch" error when this module is used.
			#  </description>
			# </workaround>
			$stateStore->{'outputCode'} .= "   write (label,'(i16)') 0\n";
			#$stateStore->{'outputCode'} .= "   write (label,'(i16)') sizeof(".$_.")\n";
			$stateStore->{'outputCode'} .= "  call displayMessage('storing \"".$_."\" with size '//trim(adjustl(label))//' bytes')\n";
		    }
		    $stateStore->{'outputCode'} .= " end if\n";
		    $stateStore->{'outputCode'} .= "  write (stateFile) ".join(",",map {"self%".$_."%ID"} @{$declaration->{'variableNames'}})."\n";
		    $stateStore->{'inputCode'}  .= "  read  (stateFile) ".join(",",map {"self%".$_."%ID"} @{$declaration->{'variableNames'}})."\n";
		}
	    } elsif (
		(  grep {$_->{'type'} eq $type    } &List::ExtraUtils::as_array($stateStorables->{'stateStorables'        }))
		||
		(  grep {$_           eq $type    } &List::ExtraUtils::as_array($stateStorables->{'functionClassInstances'}))
		){
		# This is a non-pointer object which is explicitly stateStorable.
		# Get presence of pointer attribute.
		my $isPointer = grep {$_ eq "pointer"} @{$declaration->{'attributes'}};
		# Get list of named pointer variables that are allowed.
		my @explicits = (defined($class) && exists($class->{'stateStorable'}->{'functionClass'}->{'variables'})) ? split(/\s*,\s*/,$class->{'stateStorable'}->{'functionClass'}->{'variables'}) : ();
		my $isFunctionClass = grep {$_ eq $type} @{$stateStorables->{'functionClassInstances'}};
		# Construct code to output.
		foreach ( @{$declaration->{'variables'}} ) {
		    (my $variableName = $_) =~ s/\s*=.*$//;
		    next
			if ( grep {lc($_) eq lc($variableName)} @{$stateStore->{'excludes'}} );
		    my $isExplicit = grep {lc($_) eq lc($variableName)} @explicits;
		    next
			unless ( (! $isPointer) || $isExplicit );
		    push(@{$explicitNamesFound},lc($variableName))
			if ( $isExplicit );
		    my $rank = 0;
		    if ( grep {$_ =~ m/^dimension\s*\(/} @{$declaration->{'attributes'}} ) {
			my $dimensionDeclarator = join(",",map {/^dimension\s*\(([a-zA-Z0-9_,:\s]+)\)/} @{$declaration->{'attributes'}});
			$rank        = ($dimensionDeclarator =~ tr/,//)+1;
			$stateStores->{'rankMaximum'} = $rank
			    if ( $rank > $stateStores->{'rankMaximum'} );
		    }
		    if ( grep {$_ eq "allocatable"} @{$declaration->{'attributes'}} ) {
			# For allocatable variables we must first store the shape so that they can be reallocated on restore.
			$stateStores->{'allocatablesFound'}  = 1;
			$stateStores->{'dimensionalsFound'}  = 1
			    if ( $rank > 0 );
			$stateStore->{'outputCode'} .= " if (allocated(self%".$variableName.")) then\n";
			$stateStore->{'outputCode'} .= "  write (stateFile) .true.\n";
			$stateStore->{'outputCode'} .= "  write (stateFile) shape(self%".$variableName.",kind=c_size_t)\n"
			    if ( $rank > 0 );
			$stateStore->{'inputCode'}  .= " read (stateFile) wasAllocated\n";
			$stateStore->{'inputCode'}  .= " if (allocated(self%".$variableName.")) deallocate(self%".$variableName.")\n";
			$stateStore->{'inputCode'}  .= " if (wasAllocated) then\n";
			if ( $rank > 0 ) {
			    $stateStore->{'inputCode'}  .= "  allocate(storedShape(".$rank."))\n";
			    $stateStore->{'inputCode'}  .= "  read (stateFile) storedShape\n";
			}
			if ( $declaration->{'intrinsic'} eq "class" ) {
			    (my $storable) = grep {$_->{'type'} eq $type} @{$stateStorables->{'stateStorables'}};
			    my $functionName = $type."ClassRestore".($rank > 0 ? $rank."D" : "");
			    $stateStores->{'stateRestoreModules'}->{$storable->{'module'}.",only:".$functionName} = 1;
			    $stateStore->{'inputCode'}  .= "  call ".$functionName."(self%".$variableName.",stateFile".($rank > 0 ? ",storedShape" : "").")\n";
			} else {
			    $stateStore->{'inputCode'}  .= "  allocate(self%".$variableName.($rank > 0 ? "(".join(",",map {"storedShape(".$_.")"} 1..$rank).")" : "").")\n";
			}
			if ( $rank > 0 ) {
			    $stateStore->{'inputCode'}  .= "  deallocate(storedShape)\n";
			}
		    }
		    for(my $i=1;$i<=$rank;++$i) {
			$stateStore->{'outputCode'} .= (" " x $i)."do i".$i."=lbound(self%".$variableName.",dim=".$i."),ubound(self%".$variableName.",dim=".$i.")\n";
			$stateStore->{'inputCode'}  .= (" " x $i)."do i".$i."=lbound(self%".$variableName.",dim=".$i."),ubound(self%".$variableName.",dim=".$i.")\n";
		    }
		    my $arrayElement = $rank > 0 ? "(".join(",",map {"i".$_} 1..$rank).")" : "";
		    $stateStores->{'labelUsed'}   = 1;
		    $stateStore->{'outputCode'} .= " if (displayVerbosity() >= verbosityLevelWorking) then\n";
		    if ( $declaration->{'intrinsic'} eq "class" ) {
			$stateStore->{'outputCode'} .= "  select type (c__ => self%".$variableName.$arrayElement.")\n";
			$stateStore->{'outputCode'} .= "  class is (".$declaration->{'type'}.")\n";
			# <workaround type="gfortran" PR="94446" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=94446">
			#  <description>
			#   Using the sizeof() intrinsic on a treeNode object causes a bogus "type mismatch" error when this module is used.
			#  </description>
			# </workaround>
			$stateStore->{'outputCode'} .= "   write (label,'(i16)') 0\n";
			#$stateStore->{'outputCode'} .= "   write (label,'(i16)') sizeof(c__)\n";
			$stateStore->{'outputCode'} .= "  end select\n";
		    } else {
			# <workaround type="gfortran" PR="94446" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=94446">
			#  <description>
			#   Using the sizeof() intrinsic on a treeNode object causes a bogus "type mismatch" error when this module is used.
			#  </description>
			# </workaround>
			$stateStore->{'outputCode'} .= "   write (label,'(i16)') 0\n";
			#$stateStore->{'outputCode'} .= "   write (label,'(i16)') sizeof(self%".$variableName.")\n";
		    }
		    $stateStore->{'outputCode'} .= "  call displayMessage('storing \"".$variableName.$arrayElement."\" with size '//trim(adjustl(label))//' bytes')\n";
		    $stateStore->{'outputCode'} .= " end if\n";
		    $stateStore->{'inputCode'}  .= " call displayMessage('restoring \"".$variableName.$arrayElement."\"',verbosity=verbosityLevelWorking)\n";
		    $stateStore->{'inputCode'}  .= (" " x $rank)." call self%".$variableName.$arrayElement."%stateRestore(stateFile,gslStateFile".($isFunctionClass ? ",stateOperationID" : "").")\n";
		    $stateStore->{'outputCode'} .= (" " x $rank)." call self%".$variableName.$arrayElement."%stateStore  (stateFile,gslStateFile".($isFunctionClass ? ",stateOperationID" : ",storeIdentifier=".($declaration->{'intrinsic'} eq "class" ? ".true." : ".false.")).")\n";
		    for(my $i=1;$i<=$rank;++$i) {
			$stateStore->{'outputCode'} .= (" " x ($rank+1-$i))."end do\n";
			$stateStore->{'inputCode'}  .= (" " x ($rank+1-$i))."end do\n";
		    }
		    if ( grep {$_ eq "allocatable"} @{$declaration->{'attributes'}} ) {
			$stateStore->{'inputCode'}  .= " end if\n";
			$stateStore->{'outputCode'} .= " else\n";
			$stateStore->{'outputCode'} .= "  write (stateFile) .false.\n";
			$stateStore->{'outputCode'} .= " end if\n";
		    }
		    $stateStores->{'stateFileUsed'}    = 1;
		    $stateStores->{'gslStateFileUsed'} = 1;
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
			if ( grep {lc($_) eq lc($variableName)} @{$stateStore->{'excludes'}} );
		    $stateStores->{'allocatablesFound'}  = 1;
		    $stateStores->{'dimensionalsFound'}  = 1;
		    $stateStores->{'stateFileUsed'}      = 1;
		    $stateStores->{'labelUsed'}          = 1;
		    $stateStore->{'outputCode'}        .= " if (allocated(self%".$variableName.")) then\n";
		    $stateStore->{'outputCode'}        .= "  if (displayVerbosity() >= verbosityLevelWorking) then\n";
		    # <workaround type="gfortran" PR="94446" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=94446">
		    #  <description>
		    #   Using the sizeof() intrinsic on a treeNode object causes a bogus "type mismatch" error when this module is used.
		    #  </description>
		    # </workaround>
		    $stateStore->{'outputCode'} .= "   write (label,'(i16)') 0\n";
		    #$stateStore->{'outputCode'}        .= "   write (label,'(i16)') sizeof(self%".$variableName.")\n";
		    $stateStore->{'outputCode'}        .= "   call displayMessage('storing \"".$variableName."\" with size '//trim(adjustl(label))//' bytes')\n";
		    $stateStore->{'outputCode'}        .= "  end if\n";
		    $stateStore->{'outputCode'}        .= "  write (stateFile) .true.\n"
			. "  write (stateFile) shape(self%".$variableName.",kind=c_size_t)\n"
			. "  write (stateFile) self%".$variableName."\n";
		    $stateStore->{'outputCode'}        .= " else\n";
		    $stateStore->{'outputCode'}        .= "  write (stateFile) .false.\n";
		    $stateStore->{'outputCode'}        .= " end if\n";
		    $stateStore->{'inputCode'}         .= " read (stateFile) wasAllocated\n";
		    $stateStore->{'inputCode'}         .= " if (allocated(self%".$variableName.")) deallocate(self%".$variableName.")\n";
		    $stateStore->{'inputCode'}         .= " if (wasAllocated) then\n";
		    $stateStore->{'inputCode'}         .= "  call displayMessage('restoring \"".$variableName."\"',verbosity=verbosityLevelWorking)\n";
		    $stateStore->{'inputCode'}         .= "  allocate(storedShape(".$rank."))\n";
		    $stateStore->{'inputCode'}         .= "  read (stateFile) storedShape\n";
		    $stateStore->{'inputCode'}         .= "  allocate(self%".$variableName."(".join(",",map {"storedShape(".$_.")"} 1..$rank)."))\n";
		    $stateStore->{'inputCode'}         .= "  deallocate(storedShape)\n";
		    $stateStore->{'inputCode'}         .= "  read (stateFile) self%".$variableName."\n";
		    $stateStore->{'inputCode'}         .= " end if\n";
		}
	    } else {
		# Statically-sized variable.
		foreach ( @{$declaration->{'variables'}} ) {
		    (my $variableName = $_) =~ s/\s*=.*$//;
		    next
			if ( grep {lc($_) eq lc($variableName)} @{$stateStore->{'excludes'}} );
		    my $store = 1;
		    if ( defined($class) && exists($class->{'stateStorable'}) && exists($class->{'stateStorable'}->{'restoreTo'}) ) {
			foreach ( &List::ExtraUtils::as_array($class->{'stateStorable'}->{'restoreTo'}) ) {
			    my @variables = split(/\s*,\s*/,$_->{'variables'});
			    if ( grep {lc($_) eq lc($variableName)} @variables ) {
				$store = 0;
				$stateStore->{'inputCode'} .= " self%".$variableName."=".$_->{'state'}."\n";
			    }
			}
		    }
		    push(@{$stateStore->{'staticVariables'}},$variableName)
			if ( $store );
		}
	    }
	}
	# Check for a custom state store/restore.
	$stateStore->{'hasCustomStateStore'  } = 1
	    if
	    (
	     $declaration->{'intrinsic'} eq "procedure"
	     &&
	     $declaration->{'variables'}->[0] =~ m/^stateStore=>/
	    );
	$stateStore->{'hasCustomStateRestore'} = 1
	    if
	    (
	     $declaration->{'intrinsic'} eq "procedure"
	     &&
	     $declaration->{'variables'}->[0] =~ m/^stateRestore=>/
	    );
    }
}

sub lctrim {
    # Trim trailing whitespace and return lowercased.
    my $string = shift();
    $string =~ s/\s*$//;
    return lc($string);
}
    
1;

#  LocalWords:  nonAbstractClass
