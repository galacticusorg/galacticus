# Contains a Perl module which implements processing of "component" directives in the Galacticus build system.

package Galacticus::Build::Components;
use strict;
use warnings;
use utf8;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use DateTime;
use Data::Dumper;
use Text::Table;
use Text::Template 'fill_in_string';
use Sort::Topological qw(toposort);
use Scalar::Util 'reftype';
use XML::SAX::ParserFactory;
use XML::Validator::Schema;
use LaTeX::Encode;
use Carp 'verbose';
use Sub::Identify ':all';
use File::Changes;
use Fortran::Utils;
use Galacticus::Path;
use Galacticus::Build::Hooks;
use Galacticus::Build::Components::Utils qw(@booleanLabel $implementationPropertyNameLengthMax $fullyQualifiedNameLengthMax padClass padLinkedData padImplementationPropertyName padFullyQualified isIntrinsic isOutputIntrinsic offsetName);
use Galacticus::Build::Components::CodeGeneration;
use Galacticus::Build::Components::Hierarchy;
use Galacticus::Build::Components::Hierarchy::ODESolver;
use Galacticus::Build::Components::Hierarchy::State;
use Galacticus::Build::Components::Hierarchy::Utils;
use Galacticus::Build::Components::TreeNodes;
use Galacticus::Build::Components::TreeNodes::CreateDestroy;
use Galacticus::Build::Components::TreeNodes::Classes;
use Galacticus::Build::Components::TreeNodes::ODESolver;
use Galacticus::Build::Components::TreeNodes::Map;
use Galacticus::Build::Components::TreeNodes::Serialization;
use Galacticus::Build::Components::TreeNodes::Utils;
use Galacticus::Build::Components::NodeEvents;
use Galacticus::Build::Components::BaseTypes;
use Galacticus::Build::Components::Classes;
use Galacticus::Build::Components::Classes::Names;
use Galacticus::Build::Components::Classes::CreateDestroy;
use Galacticus::Build::Components::Classes::Defaults;
use Galacticus::Build::Components::Classes::State;
use Galacticus::Build::Components::Classes::Serialization;
use Galacticus::Build::Components::Classes::Utils;
use Galacticus::Build::Components::Implementations;
use Galacticus::Build::Components::Implementations::Utils;
use Galacticus::Build::Components::Implementations::Names;
use Galacticus::Build::Components::Implementations::CreateDestroy;
use Galacticus::Build::Components::Implementations::State;
use Galacticus::Build::Components::Implementations::Serialization;
use Galacticus::Build::Components::Implementations::ODESolver;
use Galacticus::Build::Components::Properties;
use Galacticus::Build::Components::Properties::Set;
use Galacticus::Build::Components::Properties::Utils;
use Galacticus::Build::Components::Attributes;
use Galacticus::Build::Components::DataTypes;
$SIG{ __DIE__ } = sub { Carp::confess( @_ ) };

# Insert hooks for our functions.
%Galacticus::Build::Hooks::moduleHooks = 
    (
     %Galacticus::Build::Hooks::moduleHooks,
     component => 
     {
	 parse    => \&Components_Parse_Directive,
	 generate => \&Components_Generate_Output,
	 validate => \&Components_Validate
     }
    );

# Include debugging code.
my $debugging = 0;

# Adjectives for attributes.
my %attributeAdjective =
    (
     get  => "isGettable" ,
     set  => "isSettable" ,
     rate => "isEvolvable"
    );

sub Components_Validate {
    # Validate a component document.
    my $document  = shift;
    my $file      = shift;
    my $validator = XML::Validator::Schema->new(file => &galacticusPath()."schema/componentSchema.xsd");
    my $parser    = XML::SAX::ParserFactory->parser(Handler => $validator); 
    eval { $parser->parse_string($document) };
    die "Galacticus::Build::Components::Components_Validate(): validation failed in file ".$file.":\n".$@
	if ( $@ );
}

sub Components_Parse_Directive {
    # Parse content for a "component" directive.
    my $build = shift;

    # Assert that we have a name and class.
    die("Galacticus::Build::Components::Components_Parse_Directive: no currentDocument present")
	unless ( exists($build->{'currentDocument'}           ) );
    die("Galacticus::Build::Components::Components_Parse_Directive: no name present"           )
	unless ( exists($build->{'currentDocument'}->{'name' }) );
    die("Galacticus::Build::Components::Components_Parse_Directive: no class present"          )
	unless ( exists($build->{'currentDocument'}->{'class'}) );
    # Construct an ID for this component.
    my $componentID = ucfirst($build->{'currentDocument'}->{'class'}).ucfirst($build->{'currentDocument'}->{'name'});    
    # Store a copy of the component's defining document.
    $build->{'components'}->{$componentID} = $build->{'currentDocument'};
}

sub Components_Generate_Output {
    # Generate output for a "component" directive.
    my $build = shift();

    # Sort hooks.
    my @hooks = map
    {{name => $_, hook => $Galacticus::Build::Component::Utils::componentUtils{$_}}}
    &List::ExtraUtils::sortedKeys(\%Galacticus::Build::Component::Utils::componentUtils);

    # Iterate over phases.
    print "--> Phase:\n";
    foreach my $phase ( "preValidate", "default", "gather", "scatter", "postValidate", "content", "types", "interfaces", "functions" ) {
	print "   --> ".ucfirst($phase)."...\n";
	foreach my $hook ( @hooks ) {	
	    if ( exists($hook->{'hook'}->{$phase}) ) {
		my @functions = &List::ExtraUtils::as_array($hook->{'hook'}->{$phase});
		foreach my $function ( @functions ) {
		    print "      --> ".$hook->{'name'}.(scalar(@functions) > 1 ? " {".sub_name($function)."}" : "")."\n";
		    &{$function}($build);
		}
	    }
	}
    }
    # Iterate over all functions, calling them with the build data object.
    &{$_}($build)
	foreach (
	    # Generate functions to output nodes.
	    \&Generate_Node_Output_Functions                         ,
	    # Generate output functions for each component class.
	    \&Generate_Component_Class_Output_Functions              ,
	    # Generate output functions for each implementation.
	    \&Generate_Implementation_Output_Functions               ,
	    # Generate functions for getting/setting/rating value via a deferred function.
	    \&Generate_Deferred_GSR_Function                         ,
	    # Generate functions for getting/setting/rating value directly.
	    \&Generate_GSR_Functions                                 ,
	    # Generate functions for returning which components support getting/setting/rating.
	    \&Generate_GSR_Availability_Functions                    ,
	    # Generate deferred procedure pointers.
	    \&Generate_Deferred_Procedure_Pointers                   ,
	    # Generate deferred binding procedure pointers.
	    \&Generate_Deferred_Binding_Procedure_Pointers           ,
	    # Generate deferred binding procedure pointers.
	    \&Generate_Deferred_Binding_Functions                    ,
	    # Generate required null binding functions.
	    \&Generate_Null_Binding_Functions
	);

    # Insert all derived-type variable definitions.
    &derivedTypesSerialize($build);
    # Insert all interfaces.
    &interfacesSerialize  ($build);
    # Insert all module scope variables.
    $build->{'content'} .= &Fortran::Utils::Format_Variable_Defintions($build->{'variables'})."\n";
    # Insert the "contains" line.
    $build->{'content'} .= "contains\n\n";
    # Serialize all functions.
    &functionsSerialize   ($build);
    
    # Insert all legacy-defined functions into content.
    if ( scalar(@{$build->{'code'}->{'functions'}}) > 0 ) {
	print "--> Legacy functions:\n";
	foreach my $function ( @{$build->{'code'}->{'functions'}} ) {
	    foreach my $line ( split("\n",$function) ) {
		foreach my $type ( "subroutine", "function" ) {
		    if 
			( 
			  my @matches = $line =~ $Fortran::Utils::unitOpeners{$type}->{'regEx'}
			)
		    {
			print "   --> ".$type.": ".$matches[$Fortran::Utils::unitOpeners{$type}->{'unitName'}]."\n";
		    }
		}
	    }
	}
	print "--> ".scalar(@{$build->{'code'}->{'functions'}})." legacy functions remain.....\n";
	sleep(10);
    } else {
	print "--> No legacy functions remain.....functionality should be removed\n";
	sleep(10);
    }
    $build->{'content'} .= join("\n",@{$build->{'code'}->{'functions'}})."\n";
    
    # Insert include statements to bring in all functions associated with components.
    my @includeDependencies;
    foreach my $component ( @{$build->{'componentIdList'}} ) {
     	if ( exists($build->{'components'}->{$component}->{'functions'}) ) {
     	    $build->{'content'} .= "  include \"".$build->{'components'}->{$component}->{'functions'}."\"\n";
     	    push(@includeDependencies,$build->{'components'}->{$component}->{'functions'});
     	}
    }

    # Create a Makefile to specify dependencies on these include files.
    open(makeFile,">".$ENV{'BUILDPATH'}."/Makefile_Component_Includes.tmp");
    print makeFile $ENV{'BUILDPATH'}."/objects.nodes.o:".join("",map {" ".$ENV{'BUILDPATH'}."/".$_} @includeDependencies)
	if ( scalar(@includeDependencies) > 0 );
    close(makeFile);
    &File::Changes::Update($ENV{'BUILDPATH'}."/Makefile_Component_Includes" ,$ENV{'BUILDPATH'}."/Makefile_Component_Includes.tmp" );

}

sub Generate_Deferred_Binding_Procedure_Pointers {
    # Generate deferred binding procedure pointers.
    my $build = shift;
    # Initialize data content.
    my @dataContent;
    # Initialize class pointers.
    my %classPointers;
    # Iterate over component implementations.
    foreach my $componentID ( @{$build->{'componentIdList'}} ) {
	# Get the component.
	my $component = $build->{'components'}->{$componentID};
	# Iterate over bindings.
	foreach my $binding ( @{$component->{'bindings'}->{'binding'}} ) {
	    if ( $binding->{'isDeferred'} ) {
		# Create a pointer for the component class level if needed.
		my $classFunctionName = $component->{'class'}.ucfirst($binding->{'method'});
		if ( $binding->{'bindsTo'} eq 'componentClass' && ! exists($classPointers{$classFunctionName}) ) {
		    push(
			@dataContent,
			{
			    intrinsic  => "procedure",
			    type       => $component->{'class'}.ucfirst($binding->{'method'})."Interface",
			    attributes => [ "pointer" ],
			    variables  => [ $classFunctionName."Deferred" ]
			},
			{
			    intrinsic  => "logical",
			    variables  => [ $classFunctionName."IsSetValue=.false." ]
			}
			);
		    $classPointers{$classFunctionName} = 1;
		}
		# Create a pointer for the component level.
		my $componentFunctionName = $componentID.ucfirst($binding->{'method'});
		push(
		    @dataContent,
		    {
			intrinsic  => "procedure",
			type       => $component->{'class'}.ucfirst($binding->{'method'})."Interface",
			attributes => [ "pointer" ],
			variables  => [ $componentFunctionName."Deferred" ]
		    },
		    {
			intrinsic  => "logical",
			variables  => [ $componentFunctionName."IsSetValue=.false." ]
		    }
		    );
	    }
	}
    }
    $build->{'content'} .= &Fortran::Utils::Format_Variable_Defintions(\@dataContent, indent => 2)."\n";
}

sub Generate_Deferred_Binding_Functions {
    # Generate deferred binding functions.
    my $build = shift;
    # Initialize class functions.
    my %classFunctions;
    # Initialize interfaces.
    my %interfaces;
    # Iterate over component implementations.
    foreach my $componentID ( @{$build->{'componentIdList'}} ) {
	# Get the component.
	my $component = $build->{'components'}->{$componentID};
	# Iterate over bindings.
	foreach my $binding ( @{$component->{'bindings'}->{'binding'}} ) {
	    if ( $binding->{'isDeferred'} ) {
		# Determine type and arguments of the function.
		my $type = $binding->{'interface'}->{'type'};
		($type, my $name, my $attributeList) = &Galacticus::Build::Components::DataTypes::dataObjectPrimitiveName($binding->{'interface'})
		    unless ( $type eq "void" );
		my $endType;
		my $openType;
		if ( $type eq "void" ) {
		    $openType = "subroutine";
		    $endType  = "subroutine";
		} else {
		    $openType = $type." function";
		    $endType  = "function";
		}
		my @arguments;
		my @selflessArguments;
		my @selfishArguments;
		if ( exists($binding->{'interface'}->{'argument'}) ) {
		    foreach ( @{$binding->{'interface'}->{'argument'}} ) {
			if ( $_ =~ m/::\s*([a-zA-Z0-9_,\s\(\):]+)\s*$/ ) {
			    push(
				@selflessArguments,
				&Fortran::Utils::Extract_Variables($1)
				);
			} else {
			    die "Generate_Deferred_Binding_Functions: unrecognized argument format"
			}
		    }
		}
		push(@arguments       ,"self",@selflessArguments);
		push(@selfishArguments,"self"                   )
		    if ( $binding->{'interface'}->{'self'}->{'pass'} eq "true" );
		push(@selfishArguments       ,@selflessArguments);
		# Find the highest level at which this method binds.
		my $highLevel;
		if ( $binding->{'bindsTo'} eq "componentClass" ) {
		    $highLevel = $component->{'class'};
		} else {
		    $highLevel = $componentID;
		    my $parentComponent = $component;
		    while ( exists($parentComponent->{'extends'}) ) {
			$highLevel = ucfirst($parentComponent->{'extends'}->{'class'}).ucfirst($parentComponent->{'extends'}->{'name'});
			$parentComponent = $build->{'components'}->{$highLevel};
		    }
		}
		# Create an abstract interface for the deferred function.
		my $interfaceName = $component->{'class'}.ucfirst($binding->{'method'});
		unless ( exists($build->{'interfaces'}->{$interfaceName}) ) {
		    my @data;
		    push
			(
			 @data,
			 {
			     intrinsic  => "class"                   ,
			     type       => "nodeComponent".$highLevel,
			     attributes => [ "intent(".$binding->{'interface'}->{'self'}->{'intent'}.")" ]       ,
			     variables  => [ "self" ]
			 }
			)
			if ( $binding->{'interface'}->{'self'}->{'pass'} eq "true" );
		    push(@data,&Fortran::Utils::Unformat_Variables($_))
			foreach ( @{$binding->{'interface'}->{'argument'}} );
		    $build->{'interfaces'}->{$interfaceName} =
		    {
			name      => $interfaceName."Interface",
			intrinsic => $type,
			comment   => "Interface for deferred function for {\\normalfont \\ttfamily ".$binding->{'method'}."} method of the {\\normalfont \\ttfamily ".$component->{'class'}."} class."
		    };
		    $build->{'interfaces'}->{$interfaceName}->{'data'} = \@data
			if ( scalar(@data) > 0 );
		}
		# Create functions for the component class level if needed.
		my $classFunctionName = $component->{'class'}.ucfirst($binding->{'method'});
		if ( $binding->{'bindsTo'} eq 'componentClass' && ! exists($classFunctions{$classFunctionName}) ) {
		    my @dataContent =
			(
			 {
			     intrinsic  => "procedure",
			     type       => $component->{'class'}.ucfirst($binding->{'method'})."Interface",
			     variables  => [ "deferredFunction" ]
			 },
			);
		    my $functionCode;
		    $functionCode  = "  subroutine ".$classFunctionName."DeferredFunctionSet(deferredFunction)\n";
		    $functionCode .= "    !% Set the function to be used for the {\\normalfont \\ttfamily ".$binding->{'method'}."} method of the {\\normalfont \\ttfamily ".$component->{'class'}."} component class.\n";
		    $functionCode .= "    implicit none\n";
		    $functionCode .= &Fortran::Utils::Format_Variable_Defintions(\@dataContent)."\n";
		    $functionCode .= "    ".$classFunctionName."Deferred   => deferredFunction\n";
		    $functionCode .= "    ".$classFunctionName."IsSetValue =  .true.\n";
		    $functionCode .= "    return\n";
		    $functionCode .= "  end subroutine ".$classFunctionName."DeferredFunctionSet\n";
		    # Insert into the function list.
		    push(
			@{$build->{'code'}->{'functions'}},
			$functionCode
			);
		    $functionCode  = "  logical function ".$classFunctionName."DfrrdFnctnIsSet()\n";
		    $functionCode .= "    !% Return true if the deferred function for the {\\normalfont \\ttfamily ".$binding->{'method'}."} method of the {\\normalfont \\ttfamily ".$component->{'class'}."} component class has been set.\n";
		    $functionCode .= "    implicit none\n";
		    $functionCode .= "    ".$classFunctionName."DfrrdFnctnIsSet=".$classFunctionName."IsSetValue\n";
		    $functionCode .= "    return\n";
		    $functionCode .= "  end function ".$classFunctionName."DfrrdFnctnIsSet\n";
		    # Insert into the function list.
		    push(
			@{$build->{'code'}->{'functions'}},
			$functionCode
			);
		    # Create a function to call the deferred function.
		    $functionCode  = "  ".$openType." ".$classFunctionName."(".join(",",@arguments).")\n";
		    $functionCode .= "    !% Call the deferred function for the {\\normalfont \\ttfamily ".$binding->{'method'}."} method of the {\\normalfont \\ttfamily ".$component->{'class'}."} component class.\n";
		    $functionCode .= "    use Galacticus_Error\n";
		    $functionCode .= "    implicit none\n";
		    $functionCode .= "    class(nodeComponent".ucfirst($component->{'class'})."), intent(".$binding->{'interface'}->{'self'}->{'intent'}.") :: self\n";
		    $functionCode .= "    ".$_."\n"
			foreach ( @{$binding->{'interface'}->{'argument'}} );
		    $functionCode .= "    if (self%".$binding->{'method'}."FunctionIsSet()) then\n";
		    if ( $openType eq "subroutine" ) {
			$functionCode .= "       call ".$classFunctionName."Deferred(".join(",",@selfishArguments).")\n";
		    } else {
			$functionCode .= "       ".$classFunctionName."=".$classFunctionName."Deferred(".join(",",@selfishArguments).")\n";
		    }
		    $functionCode .= "    else\n";
		    if ( $openType eq "double precision" ) {
			$functionCode .= "       ".$classFunctionName."=0.0d0\n";
		    }
		    $functionCode .= "       call Galacticus_Error_Report('".$classFunctionName."','deferred function has not been assigned')\n";
		    $functionCode .= "    end if\n";
		    $functionCode .= "    return\n";
		    $functionCode .= "  end ".$endType." ".$classFunctionName."\n";
		    # Insert into the function list.
		    push(
			@{$build->{'code'}->{'functions'}},
			$functionCode
			);
		    # Record that we have created functions for this class.
		    $classFunctions{$classFunctionName} = 1;
		}
		# Create functions for the component level.
		my $componentFunctionName = $componentID.ucfirst($binding->{'method'});
		my @dataContent =
		    (
		     {
			 intrinsic  => "procedure",
			 type       => $component->{'class'}.ucfirst($binding->{'method'})."Interface",
			 variables  => [ "deferredFunction" ]
		     },
		    );
		my $functionCode;
		$functionCode  = "  subroutine ".$componentFunctionName."DeferredFunctionSet(deferredFunction)\n";
		$functionCode .= "    !% Set the function to be used for the {\\normalfont \\ttfamily ".$binding->{'method'}."} method of the {\\normalfont \\ttfamily ".$componentID."} component.\n";
		$functionCode .= "    implicit none\n";
		$functionCode .= &Fortran::Utils::Format_Variable_Defintions(\@dataContent)."\n";
		$functionCode .= "    ".$componentFunctionName."Deferred   => deferredFunction\n";
		$functionCode .= "    ".$componentFunctionName."IsSetValue =  .true.\n";
		$functionCode .= "    return\n";
		$functionCode .= "  end subroutine ".$componentFunctionName."DeferredFunctionSet\n";
		# Insert into the function list.
		push(
		    @{$build->{'code'}->{'functions'}},
		    $functionCode
		    );
		$functionCode  = "  logical function ".$componentFunctionName."DfrrdFnctnIsSet()\n";
		$functionCode .= "    !% Return true if the deferred function for the {\\normalfont \\ttfamily ".$binding->{'method'}."} method of the {\\normalfont \\ttfamily ".$componentID."} component has been set.\n";
		$functionCode .= "    implicit none\n";
		$functionCode .= "    ".$componentFunctionName."DfrrdFnctnIsSet=".$componentFunctionName."IsSetValue\n";
		$functionCode .= "    return\n";
		$functionCode .= "  end function ".$componentFunctionName."DfrrdFnctnIsSet\n";
		# Insert into the function list.
		push(
		    @{$build->{'code'}->{'functions'}},
		    $functionCode
		    );
		# Create a function that calls the deferred function.
		$functionCode  = "  ".$openType." ".$componentFunctionName."(".join(",",@arguments).")\n";
		$functionCode .= "    !% Call the deferred function for the {\\normalfont \\ttfamily ".$binding->{'method'}."} method of the {\\normalfont \\ttfamily ".$componentID."} component.\n";
		$functionCode .= "    use Galacticus_Error\n";
		$functionCode .= "    implicit none\n";
		$functionCode .= "    class(nodeComponent".ucfirst($componentID)."), intent(".$binding->{'interface'}->{'self'}->{'intent'}.") :: self\n";
		$functionCode .= "    ".$_."\n"
		    foreach ( @{$binding->{'interface'}->{'argument'}} );
		$functionCode .= "    select type (self)\n";
		$functionCode .= "    class is (nodeComponent".ucfirst($componentID).")\n";
		$functionCode .= "       if (self%".$binding->{'method'}."FunctionIsSet()) then\n";
		if ( $openType eq "subroutine" ) {
		    $functionCode .= "          call ".$componentFunctionName."Deferred(".join(",",@selfishArguments).")\n";
		} else {
		    $functionCode .= "          ".$componentFunctionName."=".$componentFunctionName."Deferred(".join(",",@selfishArguments).")\n";
		}
		$functionCode .= "       else\n";
		my $parentType;
		if ( exists($component->{'extends'}) ) {
		    $parentType = "nodeComponent".ucfirst($component->{'extends'}->{'class'}).ucfirst($component->{'extends'}->{'name'});
		} elsif ( $binding->{'bindsTo'} eq "componentClass" ) {
		    $parentType = "nodeComponent".ucfirst($component->{'class'});
		}
		if ( defined($parentType) ) {
		    if ( $openType eq "subroutine" ) {
			$functionCode .= "          call self%".$parentType."%".$binding->{'method'}."(".join(",",@selflessArguments).")\n";
		    } else {
			$functionCode .= "          ".$componentFunctionName."=self%".$parentType."%".$binding->{'method'}."(".join(",",@selflessArguments).")\n";
		    }
		} else {
		    if ( $openType eq "double precision function" ) {
			$functionCode .= "       ".$componentFunctionName."=0.0d0\n";
		    }
		    $functionCode .= "          call Galacticus_Error_Report('".$componentFunctionName."','deferred function has not been assigned')\n";
		}
		$functionCode .= "       end if\n";
		$functionCode .= "    class default\n";
		if ( $openType eq "double precision function" ) {
		    $functionCode .= "       ".$componentFunctionName."=0.0d0\n";
		}
		$functionCode .= "       call Galacticus_Error_Report('".$componentFunctionName."','incorrect class - this should not happen')\n";
		$functionCode .= "    end select\n";
		$functionCode .= "    return\n";
		$functionCode .= "  end ".$endType." ".$componentFunctionName."\n";
		# Insert into the function list.
		push(
		    @{$build->{'code'}->{'functions'}},
		    $functionCode
		    );
	    }
	}
    }
}

sub Generate_Deferred_Procedure_Pointers {
    # Generate deferred procedure pointers.
    my $build = shift;
    # Initialize record of pointers which have been created.
    my %createdPointers;
    # Insert comment.
    $build->{'content'} .= "  ! Procedure pointers for deferred custom functions.\n";
    # Initialize data content.
    my @dataContent;
    # Iterate over component implementations.
    foreach my $componentID ( @{$build->{'componentIdList'}} ) {
	# Get the component.
	my $component = $build->{'components'}->{$componentID};
	# Get the component class name.
	my $componentClassName = $component->{'class'};
	# Iterate over properties.
	foreach my $propertyName ( &List::ExtraUtils::sortedKeys($component->{'properties'}->{'property'}) ) {
	    my $property = $component->{'properties'}->{'property'}->{$propertyName};
	    unless ( $property->{'attributes' }->{'isDeferred'} eq "" ) {
		my $selfType = "generic";
		$selfType = $component->{'class'}
		   unless ( $property->{'attributes'}->{'bindsTo'} eq "top" );
		(my $dataObject, my $label) = &Galacticus::Build::Components::DataTypes::dataObjectDefinition($property);
		my $dataType = $label.$property->{'rank'};
		# Determine where to attach.
		my $attachTo = $componentID;
		$attachTo = $componentClassName
		    if ( $property->{'attributes'}->{'bindsTo'} eq "top" );
		# Iterate over attributes.
		foreach ( "get", "set", "rate" ) {		    
		    # Determine function name.
		    my $functionLabel = lcfirst($attachTo).ucfirst($propertyName).ucfirst($_);
		    # Determine if this attribute is deferred and has not yet had a procedure pointer created.
		    if (
			$property->{'attributes' }->{'isDeferred'} =~ m/$_/ 
			&& $property->{'attributes' }->{$attributeAdjective{$_}}
			&& ! exists($createdPointers{$functionLabel})
			) {
			# Construct the template function.
			my $template = $selfType."NullBinding".ucfirst($_).$dataType."InOut";
			if ( $_ eq "get" ) {
			    $template = lcfirst($componentID).ucfirst($propertyName).ucfirst($_);
			} else {
			    # Record that a null binding function was used.
			    $build->{'nullBindingsUsed'}->{lc($template)} = 1;
			}
			# Generate the procedure pointer and a boolean to indicate if is has been attached.
			push(
			    @dataContent,
			    {
				intrinsic  => "procedure",
				type       => $template,
				attributes => [ "pointer" ],
				variables  => [ $functionLabel."Deferred" ]
			    },
			    {
				intrinsic  => "logical",
				variables  => [ $functionLabel."IsAttachedValue=.false." ]
			    },
			    );
			# Add the required null property to the list.
			$build->{'nullProperties'}->{$selfType}->{$dataType."InOut"} =
			{
			    type   => $property->{'type'},
			    rank   => $property->{'rank'},
			    intent => "inout"
			};
			# Record that this procedure pointer has been created.
			$createdPointers{$functionLabel} = 1;
		    }		    
		}
	    }
	}
    }
    # Insert data content.
    $build->{'content'} .= &Fortran::Utils::Format_Variable_Defintions(\@dataContent, indent => 2)."\n";
}

sub Generate_Node_Output_Functions {
    # Generate functions to output node properties.
    my $build = shift;

    # Create an output count function.
    my @dataContent =
	(
	 {
	     intrinsic  => "class",
	     type       => "treeNode",
	     attributes => [ "intent(inout)" ],
	     variables  => [ "self" ]
	 },
	 {
	     intrinsic  => "integer",
	     attributes => [ "intent(inout)" ],
	     variables  => [ "integerPropertyCount", "doublePropertyCount" ]
	 },
	 {
	     intrinsic  => "double precision",
	     attributes => [ "intent(in   )" ],
	     variables  => [ "time" ]
	 },
	 {
	     intrinsic  => "integer",
	     variables  => [ "i" ]
	 }
	);
    my $functionCode;
    $functionCode  = "  subroutine Node_Output_Count(self,integerPropertyCount,doublePropertyCount,time)\n";
    $functionCode .= "    !% Increment the count of properties to output for this node.\n";
    $functionCode .= "    implicit none\n";
    $functionCode .= &Fortran::Utils::Format_Variable_Defintions(\@dataContent)."\n";
    # Iterate over all component classes
    foreach ( @{$build->{'componentClassList'}} ) {	    
	$functionCode .= "    if (allocated(self%component".&padClass(ucfirst($_),[0,0]).")) then\n";
	$functionCode .= "      do i=1,size(self%component".&padClass(ucfirst($_),[0,0]).")\n";
	$functionCode .= "        call self%component".&padClass(ucfirst($_),[0,0])."(i)%outputCount(integerPropertyCount,doublePropertyCount,time,instance=i)\n";
	$functionCode .= "      end do\n";
	$functionCode .= "    end if\n";
    }
    $functionCode .= "    return\n";
    $functionCode .= "  end subroutine Node_Output_Count\n";
    # Insert into the function list.
    push(
	@{$build->{'code'}->{'functions'}},
	$functionCode
	);
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{type => "procedure", name => "outputCount", function => "Node_Output_Count", description => "Increment the count of properties to output for a node.", returnType => "\\void", arguments => "\\intzero\\ integerPropertyCount\\arginout, \\intzero\\ doublePropertyCount\\arginout, \\doublezero\\ time\\argin"}
	);
    # Create an output property names function.
    @dataContent =
	(
	 {
	     intrinsic  => "class",
	     type       => "treeNode",
	     attributes => [ "intent(inout)" ],
	     variables  => [ "self" ]
	 },
	 {
	     intrinsic  => "double precision",
	     attributes => [ "intent(in   )" ], 
	     variables  => [ "time" ]
	 },
	 {
	     intrinsic  => "integer", 
	     attributes => [ "intent(inout)" ], 
	     variables  => [ "integerProperty", "doubleProperty" ]
	 },
	 {
	     intrinsic  => "character",
	     type       => "len=*",
	     attributes => [ "intent(inout)", "dimension(:)" ], 
	     variables  => [ "integerPropertyNames", "integerPropertyComments", "doublePropertyNames", "doublePropertyComments" ]
	 },
	 {
	     intrinsic  => "double precision",
	     attributes => [ "intent(inout)", "dimension(:)" ],
	     variables  => [ "integerPropertyUnitsSI", "doublePropertyUnitsSI" ]
	 },
	 {
	     intrinsic  => "integer",
	     variables  => [ "i" ]
	 }
	);
    undef($functionCode);
    $functionCode  = "  subroutine Node_Output_Names(self,integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI,doubleProperty,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)\n";
    $functionCode .= "    !% Establish the names of properties to output for this node.\n";
    $functionCode .= "    implicit none\n";
    $functionCode .= &Fortran::Utils::Format_Variable_Defintions(\@dataContent)."\n";
    # Iterate over all component classes
    foreach ( @{$build->{'componentClassList'}} ) {	    
	$functionCode .= "    if (allocated(self%component".&padClass(ucfirst($_),[0,0]).")) then\n";
	$functionCode .= "      do i=1,size(self%component".&padClass(ucfirst($_),[0,0]).")\n";
	$functionCode .= "        call self%component".&padClass(ucfirst($_),[0,0])."(i)%outputNames(integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI,doubleProperty,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time,instance=i)\n";
	$functionCode .= "      end do\n";
	$functionCode .= "    end if\n";
    }
    $functionCode .= "    return\n";
    $functionCode .= "  end subroutine Node_Output_Names\n";
    # Insert into the function list.
    push(
	@{$build->{'code'}->{'functions'}},
	$functionCode
	);
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{type => "procedure", name => "outputNames", function => "Node_Output_Names", description => "Establish the names of properties to output for a node.", returnType  => "\\void", arguments   => "\\intzero\\ integerProperty\\arginout, \\textcolor{red}{\\textless char[*](:)\\textgreater} integerPropertyNames\\arginout, \\textcolor{red}{\\textless char[*](:)\\textgreater} integerPropertyComments\\arginout, \\doubleone\\ integerPropertyUnitsSI\\arginout, \\intzero\\ doubleProperty\\arginout, \\textcolor{red}{\\textless char[*](:)\\textgreater} doublePropertyNames\\arginout, \\textcolor{red}{\\textless char[*](:)\\textgreater} doublePropertyComments\\arginout, \\doubleone\\ doublePropertyUnitsSI\\arginout, \\doublezero\\ time\\argin"}
	);
    # Create an output function.
    @dataContent =
	(
	 {
	     intrinsic  => "class",
	     type       => "treeNode",
	     attributes => [ "intent(inout)" ],
	     variables  => [ "self" ]
	 },
	 {
	     intrinsic  => "double precision",
	     attributes => [ "intent(in   )" ], 
	     variables  => [ "time" ]
	 },
	 {
	     intrinsic  => "integer", 
	     attributes => [ "intent(inout)" ], 
	     variables  => [ "integerProperty", "integerBufferCount", "doubleProperty", "doubleBufferCount" ]
	 },
	 {
	     intrinsic  => "integer",
	     type       => "kind=kind_int8",
	     attributes => [ "intent(inout)", "dimension(:,:)" ],
	     variables  => [ "integerBuffer" ]
	 },
	 {
	     intrinsic  => "double precision",
	     attributes => [ "intent(inout)", "dimension(:,:)" ],
	     variables  => [ "doubleBuffer" ]
	 },
	 {
	     intrinsic  => "integer",
	     variables  => [ "i" ]
	 }
	);
    undef($functionCode);
    $functionCode  = "  subroutine Node_Output(self,integerProperty,integerBufferCount,integerBuffer,doubleProperty,doubleBufferCount,doubleBuffer,time)\n";
    $functionCode .= "    ! Output properties for this node.\n";
    $functionCode .= "    implicit none\n";
    $functionCode .= &Fortran::Utils::Format_Variable_Defintions(\@dataContent)."\n";
    # Iterate over all component classes
    foreach ( @{$build->{'componentClassList'}} ) {	    
	$functionCode .= "    if (allocated(self%component".&padClass(ucfirst($_),[0,0]).")) then\n";
	$functionCode .= "      do i=1,size(self%component".&padClass(ucfirst($_),[0,0]).")\n";
	$functionCode .= "        call self%component".&padClass(ucfirst($_),[0,0])."(i)%output(integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time,instance=i)\n";
	$functionCode .= "      end do\n";
	$functionCode .= "    end if\n";
    }
    $functionCode .= "    return\n";
    $functionCode .= "  end subroutine Node_Output\n";
    # Insert into the function list.
    push(
	@{$build->{'code'}->{'functions'}},
	$functionCode
	);
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{type => "procedure", name => "output", function => "Node_Output", description => "Populate output buffers with properties for a node.", returnType  => "\\void", arguments   => "\\intzero\\ integerProperty\\arginout, \\intzero\\ integerBufferCount\\arginout, \\inttwo\\ integerBuffer\\arginout, \\intzero doubleProperty\\arginout, \\intzero\\ doubleBufferCount\\arginout, \\doubletwo\\ doubleBuffer\\arginout, \\doublezero\\ time\\argin"}
	);
}

sub Generate_Implementation_Output_Functions {
    # Generate output functions for each component implementation.
    my $build = shift;
    # Iterate over component implementations.
    foreach my $componentID ( @{$build->{'componentIdList'}} ) {
	# Get the component.
	my $component = $build->{'components'}->{$componentID};
	# Find modules required.
	my %modulesRequired;
	foreach my $propertyName ( &List::ExtraUtils::sortedKeys($component->{'properties'}->{'property'}) ) {
	    my $property = $component->{'properties'}->{'property'}->{$propertyName};
	    # Check if this property is to be output.
	    if ( exists($property->{'output'}) ) {
		if ( exists($property->{'output'}->{'modules'}) ) {
		    my $moduleList = $property->{'output'}->{'modules'};
		    $moduleList =~ s/^\s*//;
		    $moduleList =~ s/\s*$//;
		    my @modules = split(/\s*,\s*/,$moduleList);
		    foreach ( @modules ) {
			$modulesRequired{$_} = 1;
		    }
		}
	    }
	}
	# Initialize function code.
	my $functionCode;
	# Create property count function.
	my @dataContent =
	    (
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($componentID),
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => ["intent(inout)" ],
		 variables  => [ "integerPropertyCount", "doublePropertyCount" ]
	     },
	     {
		 intrinsic  => "double precision",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "time" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => ["intent(in   )" ],
		 variables  => [ "instance" ]
	     }
	    );
	# Check for rank-1 outputs.
	my $counterAdded = 0;
	foreach my $propertyName ( &List::ExtraUtils::sortedKeys($component->{'properties'}->{'property'}) ) {
	    my $property = $component->{'properties'}->{'property'}->{$propertyName};
	    # Check if this property is to be output.
	    if ( exists($property->{'output'}) ) {
		# Define rank, type and value.
		my $rank;
		my $type;
		# Check if this property has any linked data in this component.
		if ( exists($property->{'linkedData'}) ) {
		    my $linkedDataName = $property->{'linkedData'};
		    my $linkedData     = $component->{'content'}->{'data'}->{$linkedDataName};
		    $rank   = $linkedData->{'rank'};
		    $type   = $linkedData->{'type'};
		} elsif ( $property->{'attributes'}->{'isVirtual'} && $property->{'attributes'}->{'isGettable'} ) {
		    $rank = $property->{'rank'};
		    $type = $property->{'type'};
		} else {
		    die("Generate_Implementation_Output_Functions(): can not output [".$propertyName."]");
		}
		# Increment the counters.
		if (&isOutputIntrinsic($type)) {
		    if ( $rank == 1 && exists($property->{'output'}->{'condition'}) && $counterAdded == 0 ) {
			push(
			    @dataContent,
			    {
				intrinsic  => "integer",
				variables  => [ "i" ]
			    }
			    );
			$counterAdded = 1;
		    }
		}
	    }
	}
	# Find all derived types to be output.
	my %outputTypes;
	foreach my $propertyName ( &List::ExtraUtils::sortedKeys($component->{'properties'}->{'property'}) ) {
	    my $property = $component->{'properties'}->{'property'}->{$propertyName};
	    # Check if property is to be output.
	    if ( exists($property->{'output'}) ) {
		# Get the type of this component.
		my $type;
		if ( exists($property->{'linkedData'}) ) {
		    my $linkedDataName = $property->{'linkedData'};
		    my $linkedData     = $component->{'content'}->{'data'}->{$linkedDataName};
		    $type   = $linkedData->{'type'};
		} else {
		    $type = $property->{'type'};		
		}
		$outputTypes{$type} = 1
		    unless (&isOutputIntrinsic($type));
	    }
	}
	my @outputTypes;
	foreach ( &List::ExtraUtils::sortedKeys(\%outputTypes) ){
	    push(
		@outputTypes,
		{
		    intrinsic => "type",
		    type      => $_,
		    variables => [ "output".ucfirst($_) ]
		}
		);
	}
	push(@dataContent,@outputTypes);
	# Generate output count function.
	undef($functionCode);
	$functionCode  = "  subroutine Node_Component_".ucfirst($componentID)."_Output_Count(self,integerPropertyCount,doublePropertyCount,time,instance)\n";
	$functionCode .= "    !% Increment output property count for a ".$component->{'name'}." implementation of the ".$component->{'class'}." component.\n";
	$functionCode .= "    use ".$_."\n"
	    foreach ( &List::ExtraUtils::sortedKeys(\%modulesRequired) );
	$functionCode .= "    implicit none\n";
	$functionCode .= &Fortran::Utils::Format_Variable_Defintions(\@dataContent)."\n";
	my %typeUsed =
	    (
	     integer => 0,
	     double  => 0
	    );
	my $extraUsed    = 0;
	my $instanceUsed = 0;
	my $functionBody = "";
	unless ( $component->{'name'} eq "null" ) {
	    # Act on the parent type if necessary.
	    if ( exists($component->{'extends'}) ) {
		$functionBody .= "    call self%nodeComponent".ucfirst($component->{'extends'}->{'class'}).ucfirst($component->{'extends'}->{'name'})."%outputCount(integerPropertyCount,doublePropertyCount,time,instance)\n";
		$extraUsed               = 1;
		$instanceUsed            = 1;
		$typeUsed    {'integer'} = 1;
		$typeUsed    {'double' } = 1;
	    }
	    # Check that this instance is to be output.
	    my $checkAdded = 0;
	    if (
		exists($component->{'output'}               )           &&
		exists($component->{'output'}->{'instances'})           &&
		$component->{'output'}->{'instances'} eq "first"
		) {
		$functionBody .= "    if (instance == 1) then\n";
		$checkAdded   = 1;
		$instanceUsed = 1;
	    }
	    # Initialize counts.
	    my %typeCount =
		(
		 double => 0,
		 integer         => 0
		);
	    my %typeMap =
		(
		 double      => "double" ,
		 integer     => "integer",
		 longInteger => "integer"
		);
	    foreach my $propertyName ( &List::ExtraUtils::sortedKeys($component->{'properties'}->{'property'}) ) {
		my $property = $component->{'properties'}->{'property'}->{$propertyName};
		# Check if this property is to be output.
		if ( exists($property->{'output'}) ) {
		    # Define rank, type and value.
		    my $rank;
		    my $type;
		    # Check if this property has any linked data in this component.
		    if ( exists($property->{'linkedData'}) ) {
			my $linkedDataName = $property->{'linkedData'};
			my $linkedData     = $component->{'content'}->{'data'}->{$linkedDataName};
			$rank   = $linkedData->{'rank'};
			$type   = $linkedData->{'type'};
		    } elsif ( $property->{'attributes'}->{'isVirtual'} && $property->{'attributes'}->{'isGettable'} ) {
			$rank = $property->{'rank'};
			$type = $property->{'type'};
		    } else {
			die("Generate_Implementation_Output_Functions(): can not output [".$propertyName."]");
		    }
		    # Determine count.
		    my $count;
		    if ( $rank == 0 ) {
			$count = 1;
		    } elsif ( $rank ==1 ) {
			die("Generate_Implementation_Output_Functions(): output of rank>0 objects requires a labels attribute")
			    unless ( exists($property->{'output'}->{'labels'}) );	
			if ( $property->{'output'}->{'labels'} =~ m/^\[(.*)\]$/ ) {
			    my $labelText = $1;
			    $labelText    =~ s/^\s*//;
			    $labelText    =~ s/\s*$//;
			    my @labels    = split(/\s*,\s*/,$labelText);
			    $count = scalar(@labels);
			} elsif ( exists($property->{'output'}->{'count'}) ) {
			    $count = $property->{'output'}->{'count'};
			} else {
			    die('Generate_Implementation_Output_Functions(): no method to determine output count of property');
			}
		    } else {
			die("Generate_Implementation_Output_Functions(): output of rank>1 arrays not supported");
		    }
		    # Increment the counters.
		    if (&isOutputIntrinsic($type)) {
			if ( $rank == 0 ) {
			    if ( exists($property->{'output'}->{'condition'}) ) {
				my $condition = $property->{'output'}->{'condition'};
				$condition =~ s/\[\[([^\]]+)\]\]/$1/g;
				$functionBody .= "    if (".$condition.") ".$typeMap{$type}."PropertyCount=".$typeMap{$type}."PropertyCount+".$count."\n";
				$typeUsed{$typeMap{$type}} = 1;
			    } elsif ( $count =~ m/^\d/ ) {
				$typeCount{$typeMap{$type}} += $count;
				$typeUsed{$typeMap{$type}} = 1;
			    } else {
				$functionBody .= "    ".$typeMap{$type}."PropertyCount=".$typeMap{$type}."PropertyCount+".$count."\n";
				$typeUsed{$typeMap{$type}} = 1;
			    }
			} elsif ( $rank == 1 ) {
			    if ( exists($property->{'output'}->{'condition'}) ) {
				my $condition = $property->{'output'}->{'condition'};
				$condition =~ s/\[\[([^\]]+)\]\]/$1/g;
				$condition =~ s/\{i\}/i/g;
				$functionBody .= "    do i=1,".$count."\n";
				$functionBody .= "    if (".$condition.") ".$typeMap{$type}."PropertyCount=".$typeMap{$type}."PropertyCount+1\n";
				$functionBody .= "    end do\n";
				$typeUsed{$typeMap{$type}} = 1;
			    } else {
				if ( $count =~ m/^\d/ ) {
				    $typeCount{$typeMap{$type}} += $count;
				} else {
				    $functionBody .= "    ".$typeMap{$type}."PropertyCount=".$typeMap{$type}."PropertyCount+".$count."\n";
				    $typeUsed{$typeMap{$type}} = 1;
				}  
			    }
			}
		    }
		    else {
			if ( exists($property->{'output'}->{'condition'}) ) {
			    my $condition = $property->{'output'}->{'condition'};
			    $condition =~ s/\[\[([^\]]+)\]\]/$1/g;
			    $functionBody .= "    if (".$condition.") then\n";
			}
			$functionBody .= "    output".ucfirst($type)."=self%".$propertyName."()\n";
			$functionBody .= "    call output".ucfirst($type)."%outputCount(integerPropertyCount,doublePropertyCount,time)\n";
			$functionBody .= "    end if\n"
			    if ( exists($property->{'output'}->{'condition'}) );
			$typeUsed    {'integer'} = 1;
			$typeUsed    {'double' } = 1;
			$extraUsed               = 1;
		    }
		}
	    }
	    $functionBody .= "    doublePropertyCount =doublePropertyCount +".$typeCount{'double' }."\n"
		unless ( $typeCount{'double' } == 0 );
	    $functionBody .= "    integerPropertyCount=integerPropertyCount+".$typeCount{'integer'}."\n"
		unless ( $typeCount{'integer'} == 0 );
	    $functionBody .= "    end if\n"
		if ( $checkAdded == 1 );
	}
	$functionBody .= "    return\n";
	$functionBody .= "  end subroutine Node_Component_".ucfirst($componentID)."_Output_Count\n";
	foreach my $type ( sort(keys(%typeUsed)) ) {
	    $functionCode .= "  !GCC\$ attributes unused :: ".$type."PropertyCount\n"
		unless ( $typeUsed{$type} );
	}	
	$functionCode .= "  !GCC\$ attributes unused :: self, time\n"
		unless ( $extraUsed );
	$functionCode .= "  !GCC\$ attributes unused :: instance\n"
		unless ( $instanceUsed );
	$functionCode .= $functionBody;
	# Insert into the function list.
	push(
	    @{$build->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Insert a type-binding for this function into the implementation type.
	push(
	    @{$build->{'types'}->{'nodeComponent'.ucfirst($componentID)}->{'boundFunctions'}},
	    {type => "procedure", name => "outputCount", function => "Node_Component_".ucfirst($componentID)."_Output_Count"},
	    );
	# Create property names function.
	@dataContent =
	    (
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($componentID),
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "double precision",
		 attributes => [ "intent(in   )" ], 
		 variables  => [ "time" ]
	     },
	     {
		 intrinsic  => "integer", 
		 attributes => [ "intent(inout)" ], 
		 variables  => [ "integerProperty", "doubleProperty" ]
	     },
	     {
		 intrinsic  => "character",
		 type       => "len=*",
		 attributes => [ "intent(inout)", "dimension(:)" ], 
		 variables  => [ "integerPropertyNames", "integerPropertyComments", "doublePropertyNames", "doublePropertyComments" ]
	     },
	     {
		 intrinsic  => "double precision",
		 attributes => [ "intent(inout)", "dimension(:)" ],
		 variables  => [ "integerPropertyUnitsSI", "doublePropertyUnitsSI" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => ["intent(in   )" ],
		 variables  => [ "instance" ]
	     }
	    );
	push(@dataContent,@outputTypes);
	undef($functionCode);
	$functionCode  = "  subroutine Node_Component_".ucfirst($componentID)."_Output_Names(self,integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI,doubleProperty,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time,instance)\n";
	$functionCode .= "    !% Establish property names for a ".$component->{'name'}." implementation of the ".$component->{'class'}." component.\n";
	$functionCode .= "    use ".$_."\n"
	    foreach ( &List::ExtraUtils::sortedKeys(\%modulesRequired) );
	$functionCode .= "    implicit none\n";
	$functionBody = "";
	my $nameCounterAdded = 0;
	$typeUsed {'integer'} = 0;
	$typeUsed {'double' } = 0;
	$extraUsed            = 0;
	$instanceUsed         = 0;
	unless ( $component->{'name'} eq "null" ) {
	    # Act on the parent type if necessary.
	    if ( exists($component->{'extends'}) ) {
		$functionBody .= "    call self%nodeComponent".ucfirst($component->{'extends'}->{'class'}).ucfirst($component->{'extends'}->{'name'})."%outputNames(integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI,doubleProperty,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time,instance)\n";
		$extraUsed               = 1;
		$instanceUsed            = 1;
		$typeUsed    {'integer'} = 1;
		$typeUsed    {'double' } = 1;
	    }
	    # Check that this instance is to be output.
	    my $checkAdded = 0;
	    if (
		exists($component->{'output'}               )           &&
		exists($component->{'output'}->{'instances'})           &&
		$component->{'output'}->{'instances'} eq "first"
		) {
		$functionBody .= "    if (instance == 1) then\n";
		$checkAdded   = 1;
		$instanceUsed = 1;
	    }
	    my %typeMap =
		(
		 double      => "double" ,
		 integer     => "integer",
		 longInteger => "integer"
		);
	    foreach my $propertyName ( &List::ExtraUtils::sortedKeys($component->{'properties'}->{'property'}) ) {
		my $property = $component->{'properties'}->{'property'}->{$propertyName};
		# Check if this property is to be output.
		if ( exists($property->{'output'}) ) {
		    # Define rank, type and value.
		    my $rank;
		    my $type;
		    # Check if this property has any linked data in this component.
		    if ( exists($property->{'linkedData'}) ) {
			my $linkedDataName = $property->{'linkedData'};
			my $linkedData     = $component->{'content'}->{'data'}->{$linkedDataName};
			$rank   = $linkedData->{'rank'};
			$type   = $linkedData->{'type'};
		    } elsif ( $property->{'attributes'}->{'isVirtual'} && $property->{'attributes'}->{'isGettable'} ) {
			$rank = $property->{'rank'};
			$type = $property->{'type'};
		    } else {
			die("Generate_Implementation_Output_Functions(): can not output [".$propertyName."]");
		    }		   
		    # Increment the counters.
		    if (&isOutputIntrinsic($type)) {
			$typeUsed{$typeMap{$type}} = 1;
			if ( $rank == 0 ) {
			    if ( exists($property->{'output'}->{'condition'}) ) {
				my $condition = $property->{'output'}->{'condition'};
				$condition =~ s/\[\[([^\]]+)\]\]/$1/g;
				$functionBody .= "    if (".$condition.") then\n";
			    }
			    $functionBody .= "       ".$typeMap{$type}."Property=".$typeMap{$type}."Property+1\n";
			    $functionBody .= "       ".$typeMap{$type}."PropertyNames   (".$typeMap{$type}."Property)='".$component->{'class'}.ucfirst($propertyName)."'\n";
			    $functionBody .= "       ".$typeMap{$type}."PropertyComments(".$typeMap{$type}."Property)='".$property->{'output'}->{'comment'  }."'\n";
			    $functionBody .= "       ".$typeMap{$type}."PropertyUnitsSI (".$typeMap{$type}."Property)=".$property->{'output'}->{'unitsInSI'}."\n";
			    $functionBody .= "    end if\n"
				if ( exists($property->{'output'}->{'condition'}) );
			} elsif ( $rank == 1 ) {
			    die("Generate_Implementation_Output_Functions(): output of rank>0 objects requires a labels attribute")
				unless ( exists($property->{'output'}->{'labels'}) );
			    
			    if ( $property->{'output'}->{'labels'} =~ m/^\[(.*)\]$/ ) {
				if ( exists($property->{'output'}->{'condition'}) ) {
				    my $condition = $property->{'output'}->{'condition'};
				    $condition =~ s/\[\[([^\]]+)\]\]/$1/g;
				    $functionBody .= "    if (".$condition.") then\n";
				}
				my $labelText = $1;
				$labelText =~ s/^\s*//;
				$labelText =~ s/\s*$//;
				my @labels    = split(/\s*,\s*/,$labelText);
				foreach my $label ( @labels ) {
				    $functionBody .= "       ".$typeMap{$type}."Property=".$typeMap{$type}."Property+1\n";
				    $functionBody .= "       ".$typeMap{$type}."PropertyNames   (".$typeMap{$type}."Property)='".$component->{'class'}.ucfirst($propertyName).$label."'\n";
				    $functionBody .= "       ".$typeMap{$type}."PropertyComments(".$typeMap{$type}."Property)='".$property->{'output'}->{'comment'  }." [".$label."]'\n";
				    $functionBody .= "       ".$typeMap{$type}."PropertyUnitsSI (".$typeMap{$type}."Property)=".$property->{'output'}->{'unitsInSI'}."\n";	
				}
				$functionBody .= "    end if\n"
				    if ( exists($property->{'output'}->{'condition'}) );
			    } elsif ( exists($property->{'output'}->{'count'}) ) {
				my $count = $property->{'output'}->{'count' };
				my $label = $property->{'output'}->{'labels'};
				$label =~ s/\{i\}/i/g;
				if ( $nameCounterAdded == 0 ) {
				    push(
					@dataContent,
					{
					    intrinsic  => "integer",
					    variables  => [ "i" ]
					}
					);
				    $nameCounterAdded = 1;
				}
				$functionBody .= "       do i=1,".$count."\n";
				if ( exists($property->{'output'}->{'condition'}) ) {
				    my $condition = $property->{'output'}->{'condition'};
				    $condition =~ s/\[\[([^\]]+)\]\]/$1/g;
				    $condition =~ s/\{i\}/i/g;
				    $functionBody .= "    if (".$condition.") then\n";
				}
				$functionBody .= "         ".$typeMap{$type}."Property=".$typeMap{$type}."Property+1\n";
				$functionBody .= "         ".$typeMap{$type}."PropertyNames   (".$typeMap{$type}."Property)='".$component->{'class'}.ucfirst($propertyName)."'//".$label."\n";
				$functionBody .= "         ".$typeMap{$type}."PropertyComments(".$typeMap{$type}."Property)='".$property->{'output'}->{'comment'  }." ['//".$label."//']'\n";
				$functionBody .= "         ".$typeMap{$type}."PropertyUnitsSI (".$typeMap{$type}."Property)=".$property->{'output'}->{'unitsInSI'}."\n";	
				$functionBody .= "    end if\n"
				    if ( exists($property->{'output'}->{'condition'}) );
				$functionBody .= "end do\n";
			    } else {
				die('Generate_Implementation_Output_Functions(): no method to determine output count of property');
			    }
			} else {
			    die("Generate_Implementation_Output_Functions(): can not output rank>1 properties");
			}
		    }
		    else {
			my $unitsInSI = "0.0d0";
			$unitsInSI = $property->{'output'}->{'unitsInSI'}
			if ( exists($property->{'output'}->{'unitsInSI'}) );
			if ( exists($property->{'output'}->{'condition'}) ) {
			    my $condition = $property->{'output'}->{'condition'};
			    $condition =~ s/\[\[([^\]]+)\]\]/$1/g;
			    $functionBody .= "    if (".$condition.") then\n";
			}
			$functionBody .= "    output".ucfirst($type)."=self%".$propertyName."()\n";			   
			$functionBody .= "    call output".ucfirst($type)."%outputNames(integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI,doubleProperty,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time,'".$component->{'class'}.ucfirst($propertyName)."','".$property->{'output'}->{'comment'}."',".$unitsInSI.")\n";
			$functionBody .= "    end if\n"
			    if ( exists($property->{'output'}->{'condition'}) );
			$typeUsed{'integer'} = 1;
			$typeUsed{'double' } = 1;
			$extraUsed           = 1;
		    }
		}
	    }
	    $functionBody .= "    end if\n"
		if ( $checkAdded == 1 );
	}
	$functionBody .= "    return\n";
	$functionBody .= "  end subroutine Node_Component_".ucfirst($componentID)."_Output_Names\n";
	$functionCode .= &Fortran::Utils::Format_Variable_Defintions(\@dataContent)."\n";
	foreach my $type ( sort(keys(%typeUsed)) ) {
	    $functionCode .= "  !GCC\$ attributes unused :: ".join(", ",map {$type.$_} ('Property','PropertyNames','PropertyComments','PropertyUnitsSI'))."\n"
		unless ( $typeUsed{$type} );
	}	
	$functionCode .= "  !GCC\$ attributes unused :: self, time\n"
		unless ( $extraUsed );
	$functionCode .= "  !GCC\$ attributes unused :: instance\n"
		unless ( $instanceUsed );
	$functionCode .= $functionBody;
	# Insert into the function list.
	push(
	    @{$build->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Insert a type-binding for this function into the implementation type.
	push(
	    @{$build->{'types'}->{'nodeComponent'.ucfirst($componentID)}->{'boundFunctions'}},
	    {type => "procedure", name => "outputNames", function => "Node_Component_".ucfirst($componentID)."_Output_Names"},
	    );
    }
}

sub Generate_Deferred_Function_Attacher {
    # Generate functions to attach a function to a deferred method and to query the attachment state.
    my $component = shift;
    my $property  = shift;
    my $build = shift;
    my $gsr       = shift;
    my $gsrSuffix = "";
    $gsrSuffix = ucfirst($gsr)
	unless ( $gsr eq "get" );
    # Get the component fully-qualified, class and metho names.
    my $componentClassName = $component->{'class'             };
    my $componentName      = $component->{'fullyQualifiedName'};
    my $propertyName       = $property->{'name'};
    # Determine where to attach.
    my $attachTo = $componentName;
    $attachTo = $componentClassName
	if ( $property->{'attributes'}->{'bindsTo'} eq "top" );
    # Define the function name.
    my $functionLabel = lcfirst($attachTo).ucfirst($propertyName).ucfirst($gsr);
    my $functionName = $functionLabel."Function";    
    # Skip if this function was already created.
    unless ( exists($build->{'deferredFunctionComponentClassMethodsMade'}->{$functionLabel}) ) {
	# Define the data content.
	my $selfType = "generic";
	$selfType = $component->{'class'}
	   unless ( $property->{'attributes'}->{'bindsTo'} eq "top" );
	(my $dataObject, my $label) = &Galacticus::Build::Components::DataTypes::dataObjectDefinition($property);
	my $dataType = $label.$property->{'rank'};
	my $type = $selfType."NullBinding".ucfirst($gsr).$dataType."InOut";
	$type = $componentName.ucfirst($propertyName).ucfirst($gsr)
	    if ( $gsr eq "get" );
	my @dataContent =
	    (
	     {
		 intrinsic  => "procedure",
		 type       => $type,
		 variables  => [ "deferredFunction" ]
	     },
	    );
	# Construct the function code.
	my $functionCode;
	$functionCode  = "  subroutine ".$functionName."(deferredFunction)\n";
	$functionCode .= "    !% Set the function to be used for ".$gsr." of the {\\normalfont \\ttfamily ".$propertyName."} property of the {\\normalfont \\ttfamily ".$attachTo."} component class.\n";
	$functionCode .= "    implicit none\n";
	$functionCode .= &Fortran::Utils::Format_Variable_Defintions(\@dataContent)."\n";
	$functionCode .= "    ".$functionLabel."Deferred       => deferredFunction\n";
	$functionCode .= "    ".$functionLabel."IsAttachedValue=  .true.\n";
	$functionCode .= "    return\n";
	$functionCode .= "  end subroutine ".$functionName."\n";
	# Insert into the function list.
	push(
	    @{$build->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Bind this function to the relevant type.
	if ( 
	    ( $property->{'attributes'}->{'bindsTo'} ne "top" && ( $gsr eq "get" || $gsr eq "set" ) ) ||
	    (                                                      $gsr eq "rate"                   )
	    ) {
	    my $functionType = "\\void";
	    $functionType = &Galacticus::Build::Components::DataTypes::dataObjectDocName($property)
		if ( $gsr eq "get" );
	    push(
		@{$build->{'types'}->{"nodeComponent".ucfirst($attachTo)}->{'boundFunctions'}},
		{type => "procedure", pass => "nopass", name => $propertyName.$gsrSuffix."Function", function => $functionName, description => "Set the function to be used for the {\\normalfont \\ttfamily ".$gsr."} method of the {\\normalfont \\ttfamily ".$propertyName."} property of the {\\normalfont \\ttfamily ".$attachTo."} component.", returnType => "\\void", arguments => "\\textcolor{red}{\\textless function()\\textgreater} deferredFunction"}
		);
	}
	# Also create a function to return whether or not the deferred function has been attached.
	$functionCode  = "  logical function ".$functionLabel."IsAttached()\n";
	$functionCode .= "    !% Return true if the deferred function used to ".$gsr." the {\\normalfont \\ttfamily ".$propertyName."} property of the {\\normalfont \\ttfamily ".$attachTo."} component class has been attached.\n";
	$functionCode .= "    implicit none\n";
	$functionCode .= "    ".$functionLabel."IsAttached=".$functionLabel."IsAttachedValue\n";
	$functionCode .= "    return\n";
	$functionCode .= "  end function ".$functionLabel."IsAttached\n";
	# Insert into the function list.
	push(
	    @{$build->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Bind this function to the relevant type.
	if ( 
	    ( $property->{'attributes'}->{'bindsTo'} ne "top" && ( $gsr eq "get" || $gsr eq "set" ) ) ||
	    (                                                      $gsr eq "rate"                   )
	    ) {
	    push(
		@{$build->{'types'}->{"nodeComponent".ucfirst($attachTo)}->{'boundFunctions'}},
		{type => "procedure", pass => "nopass", name => $propertyName.$gsrSuffix."IsAttached", function => $functionLabel."IsAttached", description => "Return whether the ".$gsr." method of the ".$propertyName." property of the {\\normalfont \\ttfamily ".$attachTo."} component has been attached to a function.", returnType => "\\logicalzero", arguments => ""}
		);
	}
	# Record that these functions have now been created.
	$build->{'deferredFunctionComponentClassMethodsMade'}->{$functionLabel} = 1;
    }
}

sub Generate_Deferred_GSR_Function {
    # Generate function to get/set/rate the value of a property via a deferred function.
    my $build = shift;
    # Record bindings already made.
    my %bindings;
    # Iterate over component implementations
    foreach my $componentID ( @{$build->{'componentIdList'}} ) {
	# Get the component.
	my $component = $build->{'components'}->{$componentID};
	# Iterate over properties.
	foreach my $propertyName ( &List::ExtraUtils::sortedKeys($component->{'properties'}->{'property'}) ) {
	    # Get the property.
	    my $property = $component->{'properties'}->{'property'}->{$propertyName};
	    # Get the component fully-qualified and class names.
	    my $componentClassName = $component->{'class'             };
	    my $componentName      = $component->{'fullyQualifiedName'};
	    # Ignore non-deferred functions.
	    unless ( $property->{'attributes'}->{'isDeferred'} eq "" ) {
		# Function code data.
		my $functionCode;
		# Get the name of the property.
		my $propertyName = $property->{'name'};
		# Get properties of the data type needed.
		(my $dataDefinition, my $label) = &Galacticus::Build::Components::DataTypes::dataObjectDefinition($property,matchOnly => 1);
		# Identify properties with a deferred get function to be built.
		if (
		    $property->{'attributes' }->{'isDeferred'} =~ m/get/ &&
		    $property->{'attributes' }->{'isGettable'}           &&
		    $property->{'getFunction'}->{'build'     } eq "true"
		    )
		{
		    # Define data content of this function.
		    @{$dataDefinition->{'variables'}} = ( $componentName.ucfirst($propertyName)."Get" );
		    my @dataContent =
			(
			 $dataDefinition,
			 {
			     intrinsic  => "class",
			     type       => "nodeComponent".ucfirst($componentName),
			     attributes => [ "intent(inout)" ],
			     variables  => [ "self" ]
			 },
			);
		    # Construct the function code.
		    $functionCode  = "  function ".$componentName.ucfirst($propertyName)."Get(self)\n";
		    $functionCode .= "    !% Get the value of the {\\normalfont \\ttfamily ".$propertyName."} property of the {\\normalfont \\ttfamily ".$componentName."} component using a deferred function.\n";
		    $functionCode .= "    implicit none\n";
		    $functionCode .= &Fortran::Utils::Format_Variable_Defintions(\@dataContent)."\n";
		    $functionCode .= "    ".$componentName.ucfirst($propertyName)."Get=".$componentName.ucfirst($propertyName)."GetDeferred(self)\n";
		    $functionCode .= "    return\n";
		    $functionCode .= "  end function ".$componentName.ucfirst($propertyName)."Get\n";
		    # Insert into the function list.
		    push(
			@{$build->{'code'}->{'functions'}},
			$functionCode
			);
		    # Bind this function to the relevant type.
		    push(
			@{$build->{'types'}->{"nodeComponent".ucfirst($componentName)}->{'boundFunctions'}},
			{type => "procedure", name => $propertyName, function => $componentName.ucfirst($propertyName)."Get"}
			);
		    # Generate an attacher function.
		    &Generate_Deferred_Function_Attacher($component,$property,$build,"get");
		}
		# Add an "intent(in)" attribute to the data definition for set and rate functions.
		push(@{$dataDefinition->{'attributes'}},"intent(in   )");
		# Identify properties with a deferred set function to be built.
		if (
		    $property->{'attributes' }->{'isDeferred'} =~ m/set/
		    && $property->{'attributes' }->{'isSettable'} 
		    && $property->{'setFunction'}->{'build'     } eq "true"
		    )
		{
		    @{$dataDefinition->{'variables'}} = ( "setValue" );
		    my @dataContent =
			(
			 $dataDefinition,
			 {
			     intrinsic  => "class",
			     type       => "nodeComponent".ucfirst($componentName),
			     attributes => [ "intent(inout)" ],
			     variables  => [ "self" ]
			 },
			);
		    $functionCode  = "  subroutine ".$componentName.ucfirst($propertyName)."Set(self,setValue)\n";
		    $functionCode .= "    !% Set the value of the {\\normalfont \\ttfamily ".$propertyName."} property of the {\\normalfont \\ttfamily ".$componentName."} component using a deferred function.\n";
		    $functionCode .= "    implicit none\n";
		    $functionCode .= &Fortran::Utils::Format_Variable_Defintions(\@dataContent)."\n";
		    $functionCode .= "    call ".$componentName.ucfirst($propertyName)."SetDeferred(self,setValue)\n";
		    $functionCode .= "    return\n";
		    $functionCode .= "  end subroutine ".$componentName.ucfirst($propertyName)."Set\n\n";
		    # Insert into the function list.
		    push(
			@{$build->{'code'}->{'functions'}},
			$functionCode
			);
		    # Bind this function to the relevant type.
		    push(
			@{$build->{'types'}->{"nodeComponent".ucfirst($componentName)}->{'boundFunctions'}},
			{type => "procedure", name => $propertyName."Set", function => $componentName.ucfirst($propertyName)."Set"}
			);
		    # Generate an attacher function.
		    &Generate_Deferred_Function_Attacher($component,$property,$build,"set");		  
		}
		# Identify properties with a deferred rate function to be built.
		if (
		    $property->{'attributes' }->{'isDeferred' } =~ m/rate/
		    && $property->{'attributes' }->{'isEvolvable'} 
		    )
		{
		    # Define data content of this function.
		    @{$dataDefinition->{'variables'}} = ( "setValue" );
		    my $type = "nodeComponent";
		    $type .= ucfirst($componentName)
			unless ( $property->{'attributes' }->{'bindsTo'} eq "top" );
		    my $attachTo = $componentName;
		    $attachTo = $componentClassName
			if ( $property->{'attributes' }->{'bindsTo'} eq "top" );
		    my $functionLabel = $componentName.ucfirst($propertyName)."Rate";
		    $functionLabel = "node".ucfirst($propertyName)."Rate"
			if ( $property->{'attributes' }->{'bindsTo'} eq "top" );
		    unless ( exists($build->{'topLevelDeferredFunctionsCreated'}->{$functionLabel}) ) {
			$build->{'topLevelDeferredFunctionsCreated'}->{$functionLabel} = 1;
			my @dataContent =
			    (
			     $dataDefinition,
			     {
				 intrinsic  => "class",
				 type       => $type,
				 attributes => [ "intent(inout)" ],
				 variables  => [ "self" ]
			     },
			     {
				 intrinsic  => "logical",
				 attributes => [ "intent(inout)", "optional" ],
				 variables  => [ "interrupt" ]
			     },
			     {
				 intrinsic  => "procedure",
			     type       => "interruptTask",
				 attributes => [ "pointer", "optional", "intent(inout)" ],
				 variables  => [ "interruptProcedure" ]
			     }
			    );
			$functionCode  = "  subroutine ".$functionLabel."(self,setValue,interrupt,interruptProcedure)\n";
			$functionCode .= "    !% Set the rate of the {\\normalfont \\ttfamily ".$propertyName."} property of the {\\normalfont \\ttfamily ".$componentName."} component using a deferred function.\n";
			$functionCode .= "    implicit none\n";
			$functionCode .= &Fortran::Utils::Format_Variable_Defintions(\@dataContent)."\n";
			$functionCode .= "    call ".$attachTo.ucfirst($propertyName)."RateDeferred(self,setValue,interrupt,interruptProcedure)\n";
			$functionCode .= "    return\n";
			$functionCode .= "  end subroutine ".$functionLabel."\n\n";
			# Insert into the function list.
			push(
			@{$build->{'code'}->{'functions'}},
			    $functionCode
			    );
			# Bind this function to the relevant type.
			my $bindingName = $type.$propertyName."Rate";
			unless ( exists($bindings{$bindingName}) && $property->{'attributes' }->{'bindsTo'} eq "top" ) {
			    push(
			    @{$build->{'types'}->{$type}->{'boundFunctions'}},
			    {type => "procedure", name => $propertyName."Rate", function => $functionLabel, description => "Cumulate to the rate of the {\\normalfont \\ttfamily ".$propertyName."} property of the {\\normalfont \\ttfamily ".$componentClassName."} component.", returnType => "\\void", arguments => &Galacticus::Build::Components::DataTypes::dataObjectDocName($property)."\\ value"}
				);
			    $bindings{$bindingName} = 1;
			}
			# Generate an attacher function.
		    &Generate_Deferred_Function_Attacher($component,$property,$build,"rate");
		    }
		}
	    }
	}
    }
}

sub Generate_GSR_Functions {
    # Generate functions to get/set/rate the value of a property directly.
    my $build = shift;
    # Initialize function code.
    my $functionCode;
    # Initialize data content.
    my @dataContent;
    # Initialize records of functions created.
    my %classRatesCreated;
    my %deferredFunctionComponentClassMethodsMade;
    # Iterate over component implementations.
    foreach my $componentID ( @{$build->{'componentIdList'}} ) {
	my $component = $build->{'components'}->{$componentID};
	# Get the parent class.
	my $componentClassName = $component->{'class'};
	# Iterate over properties.
	foreach my $propertyName ( &List::ExtraUtils::sortedKeys($component->{'properties'}->{'property'}) ) {
	    my $property = $component->{'properties'}->{'property'}->{$propertyName};
	    # Handle cases where a get function is explicitly specified for a non-deferred virtual property.
	    if (
		$property->{'attributes' }->{'isGettable'}                &&
		$property->{'getFunction'}->{'build'     } eq "false"     &&
		$property->{'getFunction'}->{'bindsTo'   } eq "component" &&
		$property->{'attributes' }->{'isDeferred'} !~ m/get/
		) {
		# No need to build the function - just insert a type-binding into the implementation type.
		push(
		    @{$build->{'types'}->{'nodeComponent'.ucfirst($componentID)}->{'boundFunctions'}},
		    {type => "procedure", name => $propertyName, function => $property->{'getFunction'}->{'content'}}
		    );
	    }
	    # Handle cases where a set function is explicitly specified for a non-deferred virtual property.
	    if (
		$property->{'attributes' }->{'isSettable'}                &&
		$property->{'setFunction'}->{'build'     } eq "false"     &&
		$property->{'setFunction'}->{'bindsTo'   } eq "component" &&
		$property->{'attributes' }->{'isDeferred'} !~ m/set/
		) {
		# No need to build the function - just insert a type-binding into the implementation type.
		push(
		    @{$build->{'types'}->{'nodeComponent'.ucfirst($componentID)}->{'boundFunctions'}},
		    {type => "procedure", name => $propertyName."Set", function => $property->{'setFunction'}->{'content'}}
		    );	
	    }
   	    # Check if this property has any linked data in this component.
	    if ( exists($property->{'linkedData'}) ) {
		# Get the linked data.
		my $linkedDataName = $property->{'linkedData'};
		my $linkedData     = $component->{'content'}->{'data'}->{$linkedDataName};
		# Create a "get" function if the property is gettable.
		if ( $property->{'attributes'}->{'isGettable'} ) {
		    # Skip get function creation if a custom function which binds at the component level has been specified.
		    unless (
			$property->{'getFunction'}->{'build'  } eq "false"     &&
			$property->{'getFunction'}->{'bindsTo'} eq "component"
			)
		    {
			# Determine the suffix for this function.
			my $suffix = "";
			$suffix = "Value"
			    if ( $property->{'attributes' }->{'isDeferred'} =~ m/get/ );
			# Specify the data content.
			(my $dataDefinition,my $label) = &Galacticus::Build::Components::DataTypes::dataObjectDefinition($linkedData);
			push(@{$dataDefinition->{'variables'}},$componentID.ucfirst($propertyName)."Get".$suffix);
			@dataContent = (
			    $dataDefinition,
			    {
				intrinsic  => "class",
				type       => "nodeComponent".ucfirst($componentID),
				attributes => [ "intent(inout)" ],
				variables  => [ "self" ]
			    }
			    );
			# Generate the code.
			$functionCode  = "  function ".$componentID.ucfirst($propertyName)."Get".$suffix."(self)\n";
			$functionCode .= "    !% Return the {\\normalfont \\ttfamily ".$propertyName."} property of the {\\normalfont \\ttfamily ".$componentID."} component implementation.\n";
			$functionCode .= "    implicit none\n";
			$functionCode .= &Fortran::Utils::Format_Variable_Defintions(\@dataContent)."\n";
			$functionCode .= "    ".$componentID.$propertyName."Get".$suffix."=self%".$linkedDataName."\n";
			$functionCode .= "    return\n";
			$functionCode .= "  end function ".$componentID.ucfirst($propertyName)."Get".$suffix."\n\n";
			# Insert into the function list.
			push(
			    @{$build->{'code'}->{'functions'}},
			    $functionCode
			    );
			# Insert a type-binding for this function into the implementation type.
			push(
			    @{$build->{'types'}->{'nodeComponent'.ucfirst($componentID)}->{'boundFunctions'}},
			    {type => "procedure", name => $propertyName.$suffix, function => $componentID.ucfirst($propertyName)."Get".$suffix, description => "Get the {\\normalfont \\ttfamily ".$propertyName."} property of the {\\normalfont \\ttfamily ".$componentClassName."} component.", returnType => &Galacticus::Build::Components::DataTypes::dataObjectDocName($property), arguments => ""}
			    );
		    }
		}
		# Create a "set" method unless the property is not settable or a custom set function has been specified.
		if ( $property->{'attributes' }->{'isSettable'} ) {
		    if ( $property->{'setFunction'}->{'build'} eq "true" ) {
			# Determine the suffix for this function.
			my $suffix = "";
			$suffix = "Value"
			    if ( $property->{'attributes' }->{'isDeferred'} =~ m/set/ );
			# Specify the data content.
			(my $dataDefinition,my $label) = &Galacticus::Build::Components::DataTypes::dataObjectDefinition($linkedData,matchOnly => 1);
			push(@{$dataDefinition->{'attributes'}},"intent(in   )");
			push(@{$dataDefinition->{'variables' }},"setValue"     );
			@dataContent = (
			    $dataDefinition,
			    {
				intrinsic  => "class",
				type       => "nodeComponent".ucfirst($componentID),
				attributes => [ "intent(inout)" ],
				variables  => [ "self" ]
			    }
			    );
			# Generate the function code.
			$functionCode  = "  subroutine ".$componentID.ucfirst($propertyName)."Set".$suffix."(self,setValue)\n";
			$functionCode .= "    !% Set the {\\normalfont \\ttfamily ".$propertyName."} property of the {\\normalfont \\ttfamily ".$componentID."} component implementation.\n";
			$functionCode .= "    use Memory_Management\n";
			$functionCode .= "    implicit none\n";
			$functionCode .= &Fortran::Utils::Format_Variable_Defintions(\@dataContent)."\n";
			# For non-real properties we also set the rate and scale content. This ensures that they get reallocated to
			# the correct size.
			if ( $linkedData->{'rank'} == 0 ) {
			    $functionCode .= "    self%".$linkedDataName."=setValue\n";
			}
			elsif ( $linkedData->{'rank'} == 1 ) {
			    $functionCode .= "    if (.not.allocated(self%".$linkedDataName.")) then\n";
			    $functionCode .= "       call    Alloc_Array  (self%".$linkedDataName.",shape(setValue))\n";
			    $functionCode .= "    else\n";
			    $functionCode .= "       if (size(self%".$linkedDataName.") /= size(setValue)) then\n";
			    $functionCode .= "          call Dealloc_Array(self%".$linkedDataName."                )\n";
			    $functionCode .= "          call Alloc_Array  (self%".$linkedDataName.",shape(setValue))\n";
			    $functionCode .= "       end if\n";
			    $functionCode .= "    end if\n";
			    $functionCode .= "    self%".$linkedDataName."=setValue\n";
			}
			$functionCode .= "    return\n";
			$functionCode .= "  end subroutine ".$componentID.ucfirst($propertyName)."Set".$suffix."\n\n";
			# Insert into the function list.
			push(
			    @{$build->{'code'}->{'functions'}},
			    $functionCode
			    );
			# Insert a type-binding for this function into the implementation type.
			push(
			    @{$build->{'types'}->{'nodeComponent'.ucfirst($componentID)}->{'boundFunctions'}},
			    {type => "procedure", name => $propertyName."Set".$suffix, function => $componentID.ucfirst($propertyName)."Set".$suffix, description => "Set the {\\normalfont \\ttfamily ".$propertyName."} property of the {\\normalfont \\ttfamily ".$componentClassName."} component.", returnType => "\\void", arguments => &Galacticus::Build::Components::DataTypes::dataObjectDocName($property)."\\ value"}
			    );
		    }
		}
		# Create "count", "rate" and "scale" functions if the property is evolvable.
		if ( $property->{'attributes'}->{'isEvolvable'} ) {
		    # Specify the "count" function data content.
		    @dataContent = (
			{
			    intrinsic  => "class",
			    type       => "nodeComponent".ucfirst($componentID),
			    attributes => [ "intent(in   )" ],
			    variables  => [ "self" ]
			}
			);
		    # Generate the "count" function code.
		    my $selfUsed   = 0;
		    $functionCode  = "  integer function ".$componentID.ucfirst($propertyName)."Count(self)\n";
		    $functionCode .= "    !% Return a count of the number of scalar properties in the {\\normalfont \\ttfamily ".$propertyName."} property of the {\\normalfont \\ttfamily ".lcfirst($componentID)."} component implementation.\n";
		    $functionCode .= "    implicit none\n";
		    $functionCode .= &Fortran::Utils::Format_Variable_Defintions(\@dataContent)."\n";
		    my $functionBody = "";
		    if ( $linkedData->{'rank'} ==  0 ) {
			$functionBody .= "    ".$componentID.$propertyName."Count=1\n";
		    }
		    elsif ( $linkedData->{'rank'} == 1 ) {
			$selfUsed      = 1;
			$functionBody .= "    if (allocated(self%".$linkedDataName.")) then\n";
			$functionBody .= "    ".$componentID.$propertyName."Count=size(self%".$linkedDataName.")\n";
			$functionBody .= "    else\n";
			$functionBody .= "    ".$componentID.$propertyName."Count=0\n";
			$functionBody .= "    end if\n";
		    }
		    $functionBody .= "    return\n";
		    $functionBody .= "  end function ".$componentID.ucfirst($propertyName)."Count\n\n";
		    $functionCode .= "   !GCC\$ attributes unused :: self\n"
			unless ( $selfUsed );
		    $functionCode .= $functionBody;
		    # Insert into the function list.
		    push(
			@{$build->{'code'}->{'functions'}},
			$functionCode
			);
		    # Insert a type-binding for this function into the implementation type.
		    push(
			@{$build->{'types'}->{'nodeComponent'.ucfirst($componentID)}->{'boundFunctions'}},
			{type => "procedure", name => $propertyName."Count", function => $componentID.ucfirst($propertyName)."Count"}
			);
		    # Get the data content for remaining functions.
		    (my $dataDefinition,my $label) = &Galacticus::Build::Components::DataTypes::dataObjectDefinition($linkedData,matchOnly => 1);
		    push(@{$dataDefinition->{'variables' }},"setValue"     );
		    push(@{$dataDefinition->{'attributes'}},"intent(in   )");
		    (my $currentDefinition,my $currentLabel) = &Galacticus::Build::Components::DataTypes::dataObjectDefinition($linkedData,matchOnly => 1);
		    push(@{$currentDefinition->{'variables' }},"current"     );
		    # If rate function is deferred, then create an intrinsic version.
		    my $rateSuffix = "";
		    $rateSuffix = "intrinsic"
			if ( $property->{'attributes'}->{'isDeferred'} =~ m/rate/ );
		    # Specify the "rate" function data content.
		    @dataContent = (
			$dataDefinition,
			{
			    intrinsic  => "class",
			    type       => "nodeComponent".ucfirst($componentID),
			    attributes => [ "intent(inout)" ],
			    variables  => [ "self" ]
			},
			{
			    intrinsic  => "logical",
			    attributes => [ "optional", "intent(inout)" ],
			    variables  => [ "interrupt" ]
			},
			{
			    intrinsic  => "procedure",
			    type       => "interruptTask",
			    attributes => [ "optional", "intent(inout)", "pointer" ],
			    variables  => [ "interruptProcedure" ]
			}
			);
		    push(@dataContent,$currentDefinition)
			unless ( $linkedData->{'type'} eq "double" );
		    # Generate the rate function code.
		    $functionCode  = "  subroutine ".$componentID.ucfirst($propertyName)."Rate".ucfirst($rateSuffix)."(self,setValue,interrupt,interruptProcedure)\n";
		    $functionCode .= "    !% Accumulate to the {\\normalfont \\ttfamily ".$propertyName."} property rate of change of the {\\normalfont \\ttfamily ".$componentID."} component implementation.\n";
		    $functionCode .= "    implicit none\n";
		    my $rateSetCode;
		    if ( $linkedData->{'type'} eq "double" ) {
			if ( $linkedData->{'rank'} == 0 ) {
			    $rateSetCode .= "    nodeRates(".&offsetName($componentID,$propertyName).")=nodeRates(".&offsetName($componentID,$propertyName).")+setValue\n";
			} else {
			    push(
				@dataContent,
				{
				    intrinsic  => "integer",
				    variables  => [ "count" ]
				}				
				);
			    $rateSetCode .= "    count=size(setValue)\n";
			    $rateSetCode .= "    nodeRates(".&offsetName($componentID,$propertyName).":".&offsetName($componentID,$propertyName)."+count-1)=nodeRates(".&offsetName($componentID,$propertyName).":".&offsetName($componentID,$propertyName)."+count-1)+setValue\n";
			}
		    } else {
			push(
			    @dataContent,
			    {
				intrinsic  => "integer",
				variables  => [ "count" ]
			    }				
			    );
			$rateSetCode .= "    count=self%".$propertyName."Data%serializeCount()\n";
			$rateSetCode .= "    if (count > 0) then\n";
			$rateSetCode .= "       current=self%".$propertyName."Data\n";
			$rateSetCode .= "       call current%deserialize(nodeRates(".&offsetName($componentID,$propertyName).":".&offsetName($componentID,$propertyName)."+count-1))\n";
			$rateSetCode .= "       call current%increment(setValue)\n";
			$rateSetCode .= "       call current%serialize(nodeRates(".&offsetName($componentID,$propertyName).":".&offsetName($componentID,$propertyName)."+count-1))\n";
			$rateSetCode .= "    end if\n";
		    }
		    $functionCode .= &Fortran::Utils::Format_Variable_Defintions(\@dataContent)."\n";
		    $functionCode .= "   !GCC\$ attributes unused :: self, interrupt, interruptProcedure\n";
		    $functionCode .= $rateSetCode;
		    $functionCode .= "    return\n";
		    $functionCode .= "  end subroutine ".$componentID.ucfirst($propertyName)."Rate".ucfirst($rateSuffix)."\n\n";
		    # Insert into the function list.
		    push(
			@{$build->{'code'}->{'functions'}},
			$functionCode
			);
		    # Insert a type-binding for this function into the implementation type.
		    my %typeDefinition = (
			type     => "procedure",
			name     => $propertyName."Rate".ucfirst($rateSuffix),
			function => $componentID.ucfirst($propertyName)."Rate".ucfirst($rateSuffix)
			);
		    if ( $rateSuffix eq "intrinsic" ) {
			$typeDefinition{'description'} = "Cumulate directly (i.e. circumventing any deferred function binding) to the rate of the {\\normalfont \\ttfamily ".$propertyName."} property of the {\\normalfont \\ttfamily ".$componentID."} component.";
			$typeDefinition{'returnType' } = "\\void";
			$typeDefinition{'arguments'  } = &Galacticus::Build::Components::DataTypes::dataObjectDocName($property)."\\ value";
		    }
		    push(
			@{$build->{'types'}->{'nodeComponent'.ucfirst($componentID)}->{'boundFunctions'}},
			\%typeDefinition
			);
		    if ( $property->{'attributes' }->{'makeGeneric'} ) {
			# Create a version of this rate function which binds to the top-level class, and so is suitable for
			# attaching to inter-component pipes.
			# Specify the data content.		
			@dataContent = (
			    $dataDefinition,
			    {
				intrinsic  => "class",
				type       => "nodeComponent",
				attributes => [ "intent(inout)" ],
				variables  => [ "self" ]
			    },
			    {
				intrinsic  => "logical",
				attributes => [ "optional", "intent(inout)" ],
				variables  => [ "interrupt" ]
			    },
			    {
				intrinsic  => "procedure",
				type       => "interruptTask",
				attributes => [ "optional", "intent(inout)", "pointer" ],
				variables  => [ "interruptProcedure" ]
			    },
			    {
				intrinsic  => "class",
				type       => "nodeComponent".ucfirst($componentClassName),
				attributes => [ "pointer" ],
				variables  => [ "this".ucfirst($componentClassName) ]
			    },
			    {
				intrinsic  => "type",
				type       => "treeNode",
				attributes => [ "pointer" ],
				variables  => [ "thisNode" ]
			    }
			    );
			# Generate the function code.
			$functionCode  = "  subroutine ".$componentID.ucfirst($propertyName)."RateGeneric(self,setValue,interrupt,interruptProcedure)\n";
			$functionCode .= "    !% Set the rate of the {\\normalfont \\ttfamily ".$propertyName."} property of the {\\normalfont \\ttfamily ".$componentID."} component via a generic {\\normalfont \\ttfamily nodeComponent}.\n";
			$functionCode .= "    use Galacticus_Error\n"
			    if ( $property->{'attributes'}->{'createIfNeeded'} );
			$functionCode .= "    implicit none\n";
			$functionCode .= &Fortran::Utils::Format_Variable_Defintions(\@dataContent)."\n";
			$functionCode .= "    thisNode => self%host()\n";
			$functionCode .= "    this".ucfirst($componentClassName)." => thisNode%".$componentClassName."()\n";
			if ( $property->{'attributes'}->{'createIfNeeded'} ) {
			    $functionCode .= "    select type (this".ucfirst($componentClassName).")\n";
			    $functionCode .= "    type is (nodeComponent".ucfirst($componentClassName).")\n";
			    $functionCode .= "      ! No specific component exists, we must interrupt and create one.\n";
			    if ( $linkedData->{'rank'} == 0 ) {
				if    ( $linkedData->{'type'} eq "double" ) {
				    $functionCode .= "   if (setValue == 0.0d0) return\n";
				}
				else {
				    $functionCode .= "   if (setValue%isZero()) return\n";
				}
			    } else {
				if    ( $linkedData->{'type'} eq "double" ) {
				    $functionCode .= "   if (all(setValue == 0.0d0)) return\n";
				}
				else {
				    die('auto-create of rank>0 objects not supported');
				}
			    }
			    $functionCode .= "      if (.not.(present(interrupt).and.present(interruptProcedure))) call Galacticus_Error_Report('".$componentID.ucfirst($propertyName)."RateGeneric','interrupt required, but optional arguments missing')\n";
			    $functionCode .= "      interrupt=.true.\n";
			    $functionCode .= "      interruptProcedure => ".$componentClassName."CreateByInterrupt\n";
			    $functionCode .= "      return\n";
			    $functionCode .= "    end select\n";
			}
			$functionCode .= "    call this".ucfirst($componentClassName)."%".$propertyName."Rate(setValue,interrupt,interruptProcedure)\n";
			$functionCode .= "    return\n";
			$functionCode .= "  end subroutine ".$componentID.ucfirst($propertyName)."RateGeneric\n\n";
			# Insert into the function list.
			push(
			    @{$build->{'code'}->{'functions'}},
			    $functionCode
			    );
		    }
		    if ( $property->{'attributes' }->{'createIfNeeded'} ) {
			# Create a version of this rate function which binds to the component class, and so can auto-create the component as needed.
			my $label = $componentClassName.ucfirst($propertyName);
			unless ( exists($classRatesCreated{$label}) ) {
			    $classRatesCreated{$label} = 1;
			    # Specify the data content.		
			    @dataContent = (
				$dataDefinition,
				{
				    intrinsic  => "class",
				    type       => "nodeComponent".ucfirst($componentClassName),
				    attributes => [ "intent(inout)" ],
				    variables  => [ "self" ]
				},
				{
				    intrinsic  => "logical",
				    attributes => [ "optional", "intent(inout)" ],
				    variables  => [ "interrupt" ]
				},
				{
				    intrinsic  => "procedure",
				    type       => "interruptTask",
				    attributes => [ "optional", "intent(inout)", "pointer" ],
				    variables  => [ "interruptProcedure" ]
				}
				);
			    # Generate the function code.
			    $functionCode  = "  subroutine ".$componentClassName.ucfirst($propertyName)."Rate(self,setValue,interrupt,interruptProcedure)\n";
			    $functionCode .= "    !% Accept a rate set for the {\\normalfont \\ttfamily ".$propertyName."} property of the {\\normalfont \\ttfamily ".$componentClassName."} component class. Trigger an interrupt to create the component.\n";
			    $functionCode .= "    use Galacticus_Error\n";
			    $functionCode .= "    implicit none\n";
			    $functionCode .= &Fortran::Utils::Format_Variable_Defintions(\@dataContent)."\n";
			    $functionCode .= "    !GCC\$ attributes unused :: self\n";
			    $functionCode .= "    ! No specific component exists, so we must interrupt and create one unless the rate is zero.\n";
			    if ( $linkedData->{'rank'} == 0 ) {
				if    ( $linkedData->{'type'} eq"double" ) {
				    $functionCode .= "   if (setValue == 0.0d0) return\n";
				}
				else {
				    $functionCode .= "   if (setValue%isZero()) return\n";
				}
			    } else {
				if    ( $linkedData->{'type'} eq"double" ) {
				    $functionCode .= "   if (all(setValue == 0.0d0)) return\n";
				}
				else {
				    die('auto-create of rank>0 objects not supported');
				}
			    }
			    $functionCode .= "    if (.not.(present(interrupt).and.present(interruptProcedure))) call Galacticus_Error_Report('".$componentClassName.ucfirst($propertyName)."Rate','interrupt required, but optional arguments missing')\n";
			    $functionCode .= "    interrupt=.true.\n";
			    $functionCode .= "    interruptProcedure => ".$componentClassName."CreateByInterrupt\n";
			    $functionCode .= "    return\n";
			    $functionCode .= "  end subroutine ".$componentClassName.ucfirst($propertyName)."Rate\n\n";
			    # Insert into the function list.
			    push(
				@{$build->{'code'}->{'functions'}},
				$functionCode
				);
			}
		    }
		    # Specify the data content for the "scale" function.		
		    @dataContent = (
			$dataDefinition,
			{
			    intrinsic  => "class",
			    type       => "nodeComponent".ucfirst($componentID),
			    attributes => [ "intent(inout)" ],
			    variables  => [ "self" ]
			}
			);
		    # Generate a function to set the "scale".
		    $functionCode  = "  subroutine ".$componentID.ucfirst($propertyName)."Scale(self,setValue)\n";
		    $functionCode .= "    !% Set the {\\normalfont \\ttfamily ".$propertyName."} property scale of the {\\normalfont \\ttfamily ".$componentID."} component implementation.\n";
		    $functionCode .= "    implicit none\n";
		    my $scaleSetCode;
		    if ( $linkedData->{'type'} eq "double" ) {
			if ( $linkedData->{'rank'} == 0 ) {
			    $scaleSetCode .= "    nodeScales(".&offsetName($componentID,$propertyName).")=setValue\n";
			} else {
			    $scaleSetCode .= "    nodeScales(".&offsetName($componentID,$propertyName).":".&offsetName($componentID,$propertyName)."+size(setValue))=setValue\n";
			}
		    } else {
			push(
			    @dataContent,
			    {
				intrinsic  => "integer",
				variables  => [ "count" ]
			    }
			    );
			$scaleSetCode .= "    count=setValue%serializeCount()\n";
			$scaleSetCode .= "    if (count > 0) call setValue%serialize(nodeScales(".&offsetName($componentID,$propertyName).":".&offsetName($componentID,$propertyName)."+count-1))\n";
		    }
		    $functionCode .= &Fortran::Utils::Format_Variable_Defintions(\@dataContent)."\n";
		    $functionCode .= "    !GCC\$ attributes unused :: self\n";
		    $functionCode .= $scaleSetCode;
		    $functionCode .= "    return\n";
		    $functionCode .= "  end subroutine ".$componentID.ucfirst($propertyName)."Scale\n\n";
		    # Insert into the function list.
		    push(
			@{$build->{'code'}->{'functions'}},
			$functionCode
			);
		    # Insert a type-binding for this function into the implementation type.
		    push(
			@{$build->{'types'}->{'nodeComponent'.ucfirst($componentID)}->{'boundFunctions'}},
			{type => "procedure", name => $propertyName."Scale", function => $componentID.ucfirst($propertyName)."Scale"}
			);
		}
	    }
	}
    }
}

sub Generate_GSR_Availability_Functions {
    # Generate functions to return text described which components support setting/getting/rating of a particular property.
    my $build = shift;
    # Iterate over classes.
    foreach my $componentClassName ( @{$build->{'componentClassList'}} ) {
	# Initialize a structure of properties.
	my $properties;
	# Iterate over class members.
	foreach my $componentName ( @{$build->{'componentClasses'}->{$componentClassName}->{'memberNames'}} ) {
	    # Get the component.
	    my $componentID = ucfirst($componentClassName).ucfirst($componentName);
	    my $component   = $build->{'components'}->{$componentID};
	    # Iterate over component and parents.
	    while ( defined($component) ) {
		# Iterate over the properties of this implementation.
	        foreach my $propertyName ( &List::ExtraUtils::sortedKeys($component->{'properties'}->{'property'}) ) {
		    # Get the property.
		    my $property = $component->{'properties'}->{'property'}->{$propertyName};
		    # Record attributes.
		    $properties->{$propertyName}->{$componentName}->{'set' } = $property->{'attributes'}->{'isSettable' }; 
		    $properties->{$propertyName}->{$componentName}->{'get' } = $property->{'attributes'}->{'isGettable' }; 
		    $properties->{$propertyName}->{$componentName}->{'rate'} = $property->{'attributes'}->{'isEvolvable'}; 
		}
		if ( exists($component->{'extends'}) ) {
		    my $parentID = ucfirst($component->{'extends'}->{'class'}).ucfirst($component->{'extends'}->{'name'});
		    $component = $build->{'components'}->{$parentID};
		} else {
		    undef($component);
		}
	    }
	}
	# Iterate over properties, creating a function for each.
	foreach my $propertyName ( &List::ExtraUtils::sortedKeys($properties) ) {
	    my $property = $properties->{$propertyName};
	    my $functionName = $componentClassName.ucfirst($propertyName)."AttributeMatch";
	    my $functionCode;
	    $functionCode  = "  function ".$functionName."(requireSettable,requireGettable,requireEvolvable)\n";
	    $functionCode .= "   !% Return a text list of component implementations in the {\\normalfont \\ttfamily ".$componentClassName."} class that have the desired attributes for the {\\normalfont \\ttfamily ".$propertyName."} property\n";
	    $functionCode .= "   use ISO_Varying_String\n";
	    $functionCode .= "   implicit none\n";
	    $functionCode .= "   type   (varying_string), allocatable  , dimension(:) :: ".$functionName."\n";
	    $functionCode .= "   logical                , intent(in   ), optional     :: requireSettable      , requireGettable      , requireEvolvable\n";
	    $functionCode .= "   logical                                              :: requireSettableActual, requireGettableActual, requireEvolvableActual\n";
	    $functionCode .= "   type   (varying_string), allocatable  , dimension(:) :: temporaryList\n\n";
	    $functionCode .= "   requireSettableActual =.false.\n";
	    $functionCode .= "   requireGettableActual =.false.\n";
	    $functionCode .= "   requireEvolvableActual=.false.\n";
	    $functionCode .= "   if (present(requireSettable )) requireSettableActual =requireSettable\n";
	    $functionCode .= "   if (present(requireGettable )) requireGettableActual =requireGettable\n";
	    $functionCode .= "   if (present(requireEvolvable)) requireEvolvableActual=requireEvolvable\n";
	    # Iterate over component implementations.
	    foreach my $componentName ( &List::ExtraUtils::sortedKeys($property) ) {
		my $component = $property->{$componentName};
		my @logic;
		push(@logic,".not.requireSettableActual" )
		    unless ( $component->{'set' } );
		push(@logic,".not.requireGettableActual" )
		    unless ( $component->{'get' } );
		push(@logic,".not.requireEvolvableActual")
		    unless ( $component->{'rate'} );
		my $logicCode;
		if ( @logic ) {
		    $logicCode .= "   if (".join(".and.",@logic).") then\n";
		}
		$functionCode .= $logicCode
		    if ( defined($logicCode) );
		$functionCode .= "    if (allocated(".$functionName.")) then\n";
		$functionCode .= "     call Move_Alloc(".$functionName.",temporaryList)\n";
		$functionCode .= "     allocate(".$functionName."(size(temporaryList)+1))\n";
		$functionCode .= "     ".$functionName."(1:size(temporaryList))=temporaryList\n";
		$functionCode .= "     deallocate(temporaryList)\n";
		$functionCode .= "    else\n";
		$functionCode .= "     allocate(".$functionName."(1))\n";
		$functionCode .= "    end if\n";
		$functionCode .= "    ".$functionName."(size(".$functionName."))='".$componentName."'\n";
		$functionCode .= "   end if\n"
		    if ( defined($logicCode) );
	    }
	    $functionCode .= "   return\n";
	    $functionCode .= "  end function ".$functionName."\n\n";
	    # Insert into the function list.
	    push(
		@{$build->{'code'}->{'functions'}},
		$functionCode
		);
	    # Bind this function to the relevant type.
	    push(
		@{$build->{'types'}->{'nodeComponent'.ucfirst($componentClassName)}->{'boundFunctions'}},
		{type => "procedure", pass => "nopass", name => $propertyName."AttributeMatch", function => $functionName, description => "Return a list of implementations that provide the given list off attributes for the {\\normalfont \\ttfamily ".$propertyName."} property of the {\\normalfont \\ttfamily ".$componentClassName."} component", returnType => "\\textcolor{red}{\\textless type(varying\\_string)(:)\\textgreater}", arguments => "\\logicalzero [requireGettable]\\argin, \\logicalzero [requireSettable]\\argin, \\logicalzero [requireEvolvable]\\argin"}
		);
	}
    }
}

sub Generate_Component_Class_Output_Functions {
    # Generate output for each component class.
    my $build = shift;
    # Iterate over component classes.
    foreach my $componentClassName ( @{$build->{'componentClassList'}} ) {
	# Initialize function code.
	my $functionCode;
	# Create property count function.
	my @dataContent =
	    (
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($componentClassName),
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => [ "intent(inout)" ],
		 variables  => [ "integerPropertyCount", "doublePropertyCount" ]
	     },
	     {
		 intrinsic  => "double precision",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "time" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => ["intent(in   )" ],
		 variables  => [ "instance" ]
	     },
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($componentClassName),
		 attributes => [ "allocatable" ],
		 variables  => [ "selfDefault" ]
	     }
	    );
	undef($functionCode);
	$functionCode  = "  subroutine Node_Component_".ucfirst($componentClassName)."_Output_Count(self,integerPropertyCount,doublePropertyCount,time,instance)\n";
	$functionCode .= "    !% Increment the count of properties to output for a generic ".$componentClassName." component.\n";
	$functionCode .= "    implicit none\n";
	$functionCode .= &Fortran::Utils::Format_Variable_Defintions(\@dataContent)."\n";
	$functionCode .= "    allocate(selfDefault,source=default".ucfirst($componentClassName)."Component)\n";
	$functionCode .= "    selfDefault%hostNode => self%hostNode\n";
	$functionCode .= "    call selfDefault%outputCount(integerPropertyCount,doublePropertyCount,time,instance)\n";
	$functionCode .= "    return\n";
	$functionCode .= "  end subroutine Node_Component_".ucfirst($componentClassName)."_Output_Count\n";
	# Insert into the function list.
	push(
	    @{$build->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Insert a type-binding for this function into the implementation type.
	push(
	    @{$build->{'types'}->{'nodeComponent'.ucfirst($componentClassName)}->{'boundFunctions'}},
	    {type => "procedure", name => "outputCount", function => "Node_Component_".ucfirst($componentClassName)."_Output_Count"},
	    );
	# Create property names function.
	@dataContent =
	    (
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($componentClassName),
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "double precision",
		 attributes => [ "intent(in   )" ], 
		 variables  => [ "time" ]
	     },
	     {
		 intrinsic  => "integer", 
		 attributes => [ "intent(inout)" ], 
		 variables  => [ "integerProperty", "doubleProperty" ]
	     },
	     {
		 intrinsic  => "character",
		 type       => "len=*",
		 attributes => [ "intent(inout)", "dimension(:)" ], 
		 variables  => [ "integerPropertyNames", "integerPropertyComments", "doublePropertyNames", "doublePropertyComments" ]
	     },
	     {
		 intrinsic  => "double precision",
		 attributes => [ "intent(inout)", "dimension(:)" ],
		 variables  => [ "integerPropertyUnitsSI", "doublePropertyUnitsSI" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => ["intent(in   )" ],
		 variables  => [ "instance" ]
	     },
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($componentClassName),
		 attributes => [ "allocatable" ],
		 variables  => [ "selfDefault" ]
	     }
	    );
	undef($functionCode);
	$functionCode  = "  subroutine Node_Component_".ucfirst($componentClassName)."_Output_Names(self,integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI,doubleProperty,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time,instance)\n";
	$functionCode .= "    !% Establish property names for a generic ".$componentClassName." component.\n";
	$functionCode .= "    implicit none\n";
	$functionCode .= &Fortran::Utils::Format_Variable_Defintions(\@dataContent)."\n";
	$functionCode .= "    allocate(selfDefault,source=default".ucfirst($componentClassName)."Component)\n";
	$functionCode .= "    selfDefault%hostNode => self%hostNode\n";
	$functionCode .= "    call selfDefault%outputNames(integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI,doubleProperty,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time,instance)\n";
	$functionCode .= "    return\n";
	$functionCode .= "  end subroutine Node_Component_".ucfirst($componentClassName)."_Output_Names\n";
	# Insert into the function list.
	push(
	    @{$build->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Insert a type-binding for this function into the implementation type.
	push(
	    @{$build->{'types'}->{'nodeComponent'.ucfirst($componentClassName)}->{'boundFunctions'}},
	    {type => "procedure", name => "outputNames", function => "Node_Component_".ucfirst($componentClassName)."_Output_Names"},
	    );
	# Create output function.
	my %typeMap =
	    (
	     double => "double" ,
	     integer         => "integer",
	     longInteger     => "integer"
	    );
	@dataContent =
	    (
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($componentClassName),
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "double precision",
		 attributes => [ "intent(in   )" ], 
		 variables  => [ "time" ]
	     },
	     {
		 intrinsic  => "integer", 
		 attributes => [ "intent(inout)" ], 
		 variables  => [ "integerProperty", "integerBufferCount", "doubleProperty", "doubleBufferCount" ]
	     },
	     {
		 intrinsic  => "integer",
		 type       => "kind=kind_int8",
		 attributes => [ "intent(inout)", "dimension(:,:)" ],
		 variables  => [ "integerBuffer" ]
	     },
	     {
		 intrinsic  => "double precision",
		 attributes => [ "intent(inout)", "dimension(:,:)" ],
		 variables  => [ "doubleBuffer" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => ["intent(in   )" ],
		 variables  => [ "instance" ]
	     }
	    );
	# Find all derived types to be output.
	my %outputTypes;
	my %rank1OutputTypes;
	foreach my $componentName ( @{$build->{'componentClasses'}->{$componentClassName}->{'memberNames'}} ) {
	    # Get the component.
	    my $componentID  = ucfirst($componentClassName).ucfirst($componentName);
	    my $component    = $build->{'components'}->{$componentID};
	    foreach my $propertyName ( &List::ExtraUtils::sortedKeys($component->{'properties'}->{'property'}) ) {
		my $property = $component->{'properties'}->{'property'}->{$propertyName};
		# Check if property is to be output.
		if ( exists($property->{'output'}) ) {
		    # Get the type of this component.
		    my $type;
		    if ( exists($property->{'linkedData'}) ) {
			my $linkedDataName = $property->{'linkedData'};
			my $linkedData     = $component->{'content'}->{'data'}->{$linkedDataName};
			$type   = $linkedData->{'type'};
		    } else {
			$type = $property->{'type'};		
		    }
		    $outputTypes{$type} = 1
			unless (&isOutputIntrinsic($type));
		    $rank1OutputTypes{$type} = 1
			if ( &isOutputIntrinsic($type) && $property->{'rank'} == 1 && exists($property->{'output'}->{'condition'}) );
		}
	    }
	}
	my %intrinsicMap = 
	    (
	     integer         => "integer(kind=kind_int8)",
	     longInteger     => "integer(kind=kind_int8)",
	     double => "double precision"
	    );
	foreach ( &List::ExtraUtils::sortedKeys(\%rank1OutputTypes) ) {
	    push(
		@dataContent,
		{
		    intrinsic  => $intrinsicMap{$_},
		    attributes => [ "allocatable", "dimension(:)" ],
		    variables  => [ "outputRank1".ucfirst($typeMap{$_}) ]
		}
		);
	}
	push(
	    @dataContent,
	    {
		intrinsic  => "integer",
		variables  => [ "i" ]
	    }
	    )
	    if ( scalar(keys(%rank1OutputTypes)) > 0 );
	my @outputTypes;
	foreach ( &List::ExtraUtils::sortedKeys(\%outputTypes) ){
	    push(
		@outputTypes,
		{
		    intrinsic => "type",
		    type      => $_,
		    variables => [ "output".ucfirst($_) ]
		}
		);
	}
	push(@dataContent,@outputTypes);
	undef($functionCode);
	# Find modules required.
	my %modulesRequired;
	foreach my $componentName ( @{$build->{'componentClasses'}->{$componentClassName}->{'memberNames'}} ) {
	    # Get the component.
	    my $componentID  = ucfirst($componentClassName).ucfirst($componentName);
	    my $component    = $build->{'components'}->{$componentID};
	    foreach my $propertyName ( &List::ExtraUtils::sortedKeys($component->{'properties'}->{'property'}) ) {
		my $property = $component->{'properties'}->{'property'}->{$propertyName};
		# Check if this property is to be output.
		if ( exists($property->{'output'}) ) {
		    if ( exists($property->{'output'}->{'modules'}) ) {
			my $moduleList = $property->{'output'}->{'modules'};
			$moduleList =~ s/^\s*//;
			$moduleList =~ s/\s*$//;
			my @modules = split(/\s*,\s*/,$moduleList);
			foreach ( @modules ) {
			    $modulesRequired{$_} = 1;
			}
		    }
		}
	    }
	}
	$functionCode  = "  subroutine Node_Component_".ucfirst($componentClassName)."_Output(self,integerProperty,integerBufferCount,integerBuffer,doubleProperty,doubleBufferCount,doubleBuffer,time,instance)\n";
	$functionCode .= "    !% Output properties for a ".$componentClassName." component.\n";
	$functionCode .= "    use ".$_."\n" 
	    foreach ( &List::ExtraUtils::sortedKeys(\%modulesRequired) );
	$functionCode .= "    implicit none\n";
	$functionCode .= &Fortran::Utils::Format_Variable_Defintions(\@dataContent)."\n";
	my $functionBody = "";
	my $selfUsed     = 0;
	my $instanceUsed = 0;
	my $timeUsed     = 0;
	my %typeUsed     =
	    (
	     integer => 0,
	     double  => 0
	    );
        foreach my $componentName ( @{$build->{'componentClasses'}->{$componentClassName}->{'memberNames'}} ) {
	    # Get the component.
	    my $componentID  = ucfirst($componentClassName).ucfirst($componentName);
	    my $component    = $build->{'components'}->{$componentID};
	    my $activeCheck  = "    if (default".ucfirst($componentClassName)."Component%".$componentName."IsActive()";
	    if (
		exists($component->{'output'}               )           &&
		exists($component->{'output'}->{'instances'})           &&
		$component->{'output'}->{'instances'} eq "first"
		) {
		$activeCheck .= ".and.instance == 1";
		$instanceUsed = 1;
	    }
	    $activeCheck .= ") then\n";
	    my $outputsFound = 0;
	    foreach my $propertyName ( &List::ExtraUtils::sortedKeys($component->{'properties'}->{'property'}) ) {
		my $property = $component->{'properties'}->{'property'}->{$propertyName};
		# Check if this property is to be output.
		if ( exists($property->{'output'}) ) {
		    $selfUsed = 1;
		    # Add conditional statement if necessary.
		    if ( $outputsFound == 0 ) {
			$functionBody .= $activeCheck;
			$outputsFound  = 1;
		    }
		    # Define rank, type and value.
		    my $rank;
		    my $type;
		    # Check if this property has any linked data in this component.
		    if ( exists($property->{'linkedData'}) ) {
			my $linkedDataName = $property->{'linkedData'};
			my $linkedData     = $component->{'content'}->{'data'}->{$linkedDataName};
			$rank   = $linkedData->{'rank'};
			$type   = $linkedData->{'type'};
		    } elsif ( $property->{'attributes'}->{'isVirtual'} && $property->{'attributes'}->{'isGettable'} ) {
			$rank = $property->{'rank'};
			$type = $property->{'type'};
		    } else {
			die("Generate_Component_Class_Output_Functions(): can not output [".$propertyName."]");
		    }		   
		    # Determine count.
		    my $count;
		    if ( $rank == 0 ) {
			$count = 1;
		    } elsif ( $rank ==1 ) {
			die("Generate_Component_Class_Output_Functions(): output of rank>0 objects requires a labels attribute")
			    unless ( exists($property->{'output'}->{'labels'}) );	
			if ( $property->{'output'}->{'labels'} =~ m/^\[(.*)\]$/ ) {
			    my $labelText = $1;
			    $labelText    =~ s/^\s*//;
			    $labelText    =~ s/\s*$//;
			    my @labels    = split(/\s*,\s*/,$labelText);
			    $count = scalar(@labels);
			} elsif ( exists($property->{'output'}->{'count'}) ) {
			    $count = $property->{'output'}->{'count'};
			} else {
			    die('Generate_Component_Class_Output_Functions(): no method to determine output size for rank-1 property');
			}
		    } else {
			die("Generate_Component_Class_Output_Functions(): output of rank>1 arrays not supported");
		    }
		    # Increment the counters.
		    if (&isOutputIntrinsic($type)) {
			$typeUsed{$typeMap{$type}} = 1;
			if ( $rank == 0 ) {
			    if ( exists($property->{'output'}->{'condition'}) ) {
				my $condition = $property->{'output'}->{'condition'};
				$condition =~ s/\[\[([^\]]+)\]\]/$1/g;
				$functionBody .= "    if (".$condition.") then\n";
			    }
			    $functionBody .= "       ".$typeMap{$type}."Property=".$typeMap{$type}."Property+1\n";
			    $functionBody .= "       ".$typeMap{$type}."Buffer(".$typeMap{$type}."BufferCount,".$typeMap{$type}."Property)=self%".$propertyName."()\n";
			    $functionBody .= "    end if\n"
				if ( exists($property->{'output'}->{'condition'}) );
			} else {
			    if ( exists($property->{'output'}->{'condition'}) ) {
				die("Generate_Component_Class_Output_Functions(): conditions for rank>1 properties not supported")
				    unless ( $rank == 1 );
				my $condition = $property->{'output'}->{'condition'};
				$condition =~ s/\[\[([^\]]+)\]\]/$1/g;
				$condition =~ s/\{i\}/i/g;
				$functionBody .= "    outputRank1".ucfirst($typeMap{$type})."=self%".$propertyName."()\n";
				$functionBody .= "    do i=1,".$property->{'output'}->{'count'}."\n";
				$functionBody .= "      if (".$condition.") then\n";
				$functionBody .= "        ".$typeMap{$type}."Property=".$typeMap{$type}."Property+1\n";
				$functionBody .= "        ".$typeMap{$type}."Buffer(".$typeMap{$type}."BufferCount,".$typeMap{$type}."Property)=outputRank1".ucfirst($typeMap{$type})."(i)\n";
				$functionBody .= "      end if\n";
				$functionBody .= "    end do\n";
				$functionBody .= "    deallocate(outputRank1".ucfirst($typeMap{$type}).")\n";
			    } else {
				$functionBody .= "       ".$typeMap{$type}."Buffer(".$typeMap{$type}."BufferCount,".$typeMap{$type}."Property+1:".$typeMap{$type}."Property+".$count.")=reshape(self%".$propertyName."(),[".$count."])\n";
				$functionBody .= "       ".$typeMap{$type}."Property=".$typeMap{$type}."Property+".$count."\n";
			    }
			}
		    }
		    else {
			if ( exists($property->{'output'}->{'condition'}) ) {
			    my $condition = $property->{'output'}->{'condition'};
			    $condition =~ s/\[\[([^\]]+)\]\]/$1/g;
			    $functionBody .= "    if (".$condition.") then\n";
			}
			$functionBody .= "      output".ucfirst($type)."=self%".$propertyName."()\n";
			$functionBody .= "      call output".ucfirst($type)."%output(integerProperty,integerBufferCount,integerBuffer,doubleProperty,doubleBufferCount,doubleBuffer,time)\n";
			$functionBody .= "      call self%".$propertyName."Set(output".ucfirst($type).")\n";
			$functionBody .= "    end if\n"
			    if ( exists($property->{'output'}->{'condition'}) );
			$timeUsed            = 1;
			$typeUsed{'integer'} = 1;
			$typeUsed{'double' } = 1;
		    }
		}
	    }
	    $functionBody .= "    end if\n"
		if ( $outputsFound == 1 );
	}
	$functionBody .= "    return\n";
	$functionBody .= "  end subroutine Node_Component_".ucfirst($componentClassName)."_Output\n";
	$functionCode .= "    !GCC\$ attributes unused :: self\n"
	    unless ( $selfUsed );
	$functionCode .= "    !GCC\$ attributes unused :: instance\n"
	    unless ( $instanceUsed );
	$functionCode .= "    !GCC\$ attributes unused :: time\n"
	    unless ( $timeUsed );
	foreach my $type ( sort(keys(%typeUsed)) ) {
	    $functionCode .= "    !GCC\$ attributes unused :: ".join(", ",map {$type.$_} ("Property", "BufferCount", "Buffer"))."\n"
		unless ( $typeUsed{$type} );
	}
	$functionCode .= $functionBody;
	# Insert into the function list.
	push(
	    @{$build->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Insert a type-binding for this function into the implementation type.
	push(
	    @{$build->{'types'}->{'nodeComponent'.ucfirst($componentClassName)}->{'boundFunctions'}},
	    {type => "procedure", name => "output", function => "Node_Component_".ucfirst($componentClassName)."_Output"},
	    );
    }
}

sub Generate_Null_Binding_Functions {
    # Generate null binding functions.
    my $build = shift;
    # Iterate over component classes.
    foreach my $componentClassName ( &List::ExtraUtils::sortedKeys($build->{'nullProperties'}) ) {
	# Get the null functions required for this component class.
	my $componentClass = $build->{'nullProperties'}->{$componentClassName};
	# Iterate over required null functions for this component class.
	foreach my $nullFunctionName ( &List::ExtraUtils::sortedKeys($componentClass) ) {
	    # Get the null function definition.
	    my $nullFunction = $componentClass->{$nullFunctionName};
	    # Construct a datatype for this null function.
	    (my $dataDefinition, my $label) = &Galacticus::Build::Components::DataTypes::dataObjectDefinition($nullFunction,matchOnly => 1);
	    my $labelRaw = $label;
	    # Build a label describing the intrinsic type of the data.
	    my $intrinsicType = $dataDefinition->{'intrinsic'};
	    $intrinsicType .= $dataDefinition->{'type'}
	    if ( exists($dataDefinition->{'type'}) );
	    # Append rank to the label for this data.
	    $label .= $nullFunction->{'rank'};
	    # Extract the intent of the function.
	    my $intent = $nullFunction->{'intent'};
	    # Construct the type of "self" for this function.
	    my $selfType = "nodeComponent";
	    $selfType .= ucfirst($componentClassName)
		unless ( $componentClassName eq "generic" );
	    # Add an intent to the attributes of the datatype and specify the name of the input datatype.
	    push(@{$dataDefinition->{'attributes'}},"intent(in   )");
	    @{$dataDefinition->{'variables'}} = "setValue";
	    # Build code for the null set function.
	    my @dataContent =
		(
		 $dataDefinition,
		 {
		     intrinsic  => "class",
		     type       => $selfType,
		     attributes => [ "intent(".$intent.")" ],
		     variables  => [ "self" ]
		 }
		);
	    my $functionCode;
	    my $nullBindingName = $componentClassName."NullBindingSet".$label.$intent;
	    $functionCode  = "  subroutine ".$nullBindingName."(self,setValue)\n";
	    $functionCode .= "    !% A null set function for rank ".$nullFunction->{'rank'}." ".latex_encode(lc($intrinsicType))."s.\n";
	    $functionCode .= "    implicit none\n";
	    $functionCode .= &Fortran::Utils::Format_Variable_Defintions(\@dataContent)."\n";
	    $functionCode .= "   !GCC\$ attributes unused :: ".join(", ",@{$_->{'variables'}})."\n"
		foreach ( @dataContent );
	    $functionCode .= "    return\n";
	    $functionCode .= "  end subroutine ".$nullBindingName."\n";
	    # Determine if this null rate function is actually used.
	    my $nullBindingUsed = 0;
	    foreach my $typeName ( keys(%{$build->{'types'}}) ) {
		foreach my $typeBoundFunction ( @{$build->{'types'}->{$typeName}->{'boundFunctions'}} ) {
		    if (
			(
			 exists($typeBoundFunction->{'function'})
			 &&
			 lc($typeBoundFunction->{'function'}) eq lc($nullBindingName) 
			)
			||
			(
			 exists($typeBoundFunction->{'descriptor'})
			 &&
			 lc($typeBoundFunction->{'descriptor'}->{'name'}) eq lc($nullBindingName) 
			 )
			) {
			$nullBindingUsed = 1;
			last;
		    }
		}
		last
		    if ( $nullBindingUsed );
	    }
	    $nullBindingUsed = grep {$_ eq lc($nullBindingName)} keys(%{$build->{'nullBindingsUsed'}})
		unless ( $nullBindingUsed );
	    # Insert into the function list.
	    push(
		@{$build->{'code'}->{'functions'}},
		$functionCode
		)
		if ( $nullBindingUsed && $intent ne "in" );
	    # Build code for the null rate function.
	    @dataContent =
		(
		 $dataDefinition,
		 {
		     intrinsic  => "class",
		     type       => $selfType,
		     attributes => [ "intent(".$intent.")" ],
		     variables  => [ "self" ]
		 },
		 {
		     intrinsic  => "logical",
		     attributes => [ "intent(inout)", "optional" ],
		     variables  => [ "interrupt" ]
		 },
		 {
		     intrinsic  => "procedure",
		     type       => "interruptTask", 
		     attributes => [ "intent(inout)", "optional", "pointer" ],
		     variables  => [ "interruptProcedure" ]
		 }
		);
	    $nullBindingName = $componentClassName."NullBindingRate".$label.$intent;
	    $functionCode  = "  subroutine ".$nullBindingName."(self,setValue,interrupt,interruptProcedure)\n";
	    $functionCode .= "    !% A null rate function for rank ".$nullFunction->{'rank'}." ".latex_encode(lc($intrinsicType))."s.\n";
	    $functionCode .= "    implicit none\n";
	    $functionCode .= &Fortran::Utils::Format_Variable_Defintions(\@dataContent)."\n";
	    $functionCode .= "   !GCC\$ attributes unused :: ".join(", ",@{$_->{'variables'}})."\n"
		foreach ( @dataContent );
	    $functionCode .= "    return\n";
	    $functionCode .= "  end subroutine ".$nullBindingName."\n";
	    # Determine if this null rate function is actually used.
	    $nullBindingUsed = 0;
	    foreach my $typeName ( keys(%{$build->{'types'}}) ) {
		foreach my $typeBoundFunction ( @{$build->{'types'}->{$typeName}->{'boundFunctions'}} ) {
		    if (
			(
			 exists($typeBoundFunction->{'function'})
			 &&
			 lc($typeBoundFunction->{'function'}) eq lc($nullBindingName) 
			)
			||
			(
			 exists($typeBoundFunction->{'descriptor'})
			 &&
			 lc($typeBoundFunction->{'descriptor'}->{'name'}) eq lc($nullBindingName) 
			 )
			) {
			$nullBindingUsed = 1;
			last;
		    }
		}
		last
		    if ( $nullBindingUsed );
	    }
	    $nullBindingUsed = grep {$_ eq lc($nullBindingName)} keys(%{$build->{'nullBindingsUsed'}})
		unless ( $nullBindingUsed );
	    # Insert into the function list (if it is used).
	    push(
		@{$build->{'code'}->{'functions'}},
		$functionCode
		)
		if ( $nullBindingUsed && $intent ne "in" );
	    # Build code for the null get function.
	    pop(@{$dataDefinition->{'attributes'}});
	    push(@{$dataDefinition->{'attributes'}},"allocatable")
		if ( $nullFunction->{'rank'} > 0 );
	    @{$dataDefinition->{'variables'}} = $componentClassName."NullBinding".$label.$intent;
	    @dataContent =
		(
		 $dataDefinition,
		 {
		     intrinsic  => "class",
		     type       => $selfType,
		     attributes => [ "intent(".$intent.")" ],
		     variables  => [ "self" ]
		 }		 
		);
	    my $functionName = $componentClassName."NullBinding".$label.$intent;
	    $functionCode  = "  function ".$functionName."(self)\n";
	    $functionCode .= "    !% A null get function for rank ".$nullFunction->{'rank'}." ".latex_encode(lc($intrinsicType))."s.\n";
	    $functionCode .= "    implicit none\n";
	    $functionCode .= &Fortran::Utils::Format_Variable_Defintions(\@dataContent)."\n";
	    $functionCode .= "   !GCC\$ attributes unused :: self\n";
	    if ( $nullFunction->{'rank'} == 0 ) {
		if    ( $labelRaw eq "Double" ) {
		    $functionCode .= "    ".$functionName."=0.0d0\n";
		}
		elsif ( $labelRaw eq "Integer"     ) {
		    $functionCode .= "    ".$functionName."=0\n";
		}
		elsif ( $labelRaw eq "LongInteger" ) {
		    $functionCode .= "    ".$functionName."=0_kind_int8\n";
		}
		elsif ( $labelRaw eq "Logical"     ) {
		    $functionCode .= "    ".$functionName."=.false.\n";
		}
		else {
		    $functionCode .= "     call ".$functionName."%reset()\n";
		}
	    } else {
		$functionCode .= "    ".$functionName."=null".$labelRaw.$nullFunction->{'rank'}."d\n";
	    }
	    $functionCode .= "    return\n";
	    $functionCode .= "  end function ".$functionName."\n\n";
	    # Insert into the function list.
	    push(
		@{$build->{'code'}->{'functions'}},
		$functionCode
		);
	}
    }
}

sub boundFunctionTable {
    # Get the list of type-bound functions.
    my $objectName         = shift();
    my @typeBoundFunctions = 
	sort
         {$a->{'functionName'} cmp $b->{'functionName'}}
         map
          {$_->{'functionName'} = exists($_->{'descriptor'}) ? $_->{'descriptor'}->{'name'} : $_->{'function'}; $_}
          @{shift()};
    # Create a text table object suitable for type-bound function definitions.
    my $table =  Text::Table->new(
	{
	    is_sep => 1,
	    body   => "     "
	},
	{
	    align  => "left"
	},
	{
	    align  => "left"
	},
	{
	    align  => "left"
	},
	{
	    align  => "left"
	},
 	{
	    align  => "left"
	},
        {
	    align  => "left"
	},
	{
	    align  => "left"
	}
	);
    # Iterate over type-bound functions and insert them into the table.
    foreach ( @typeBoundFunctions ) {
	# Determine pass status.
	my $pass = "";
	$pass = ", ".$_->{'pass'}
	   if ( exists($_->{'pass'}) );
	# Determine the connector to use.
	my $connector = "";
	$connector = " => "
	    if ( exists($_->{'name'}) );
	# Determine the name to use.
	my $name = "";
	$name = $_->{'name'}
	    if ( exists($_->{'name'}) );
	# Add a row to the table.
	if ( defined(reftype($_->{'function'})) && reftype($_->{'function'}) eq "ARRAY" ) {
	    # Multiple functions specified. List them, one per row.
	    my $i = 0;
	    foreach my $function ( sort(@{$_->{'function'}}) ) {
		++$i;
		# Determine a suitable suffix for this line.
		my $suffix = ", &";
		$suffix = ""
		    if ( $i == scalar(@{$_->{'function'}}) );
		# Add the line to the table.
		if ( $i == 1 ) {
		    $table->add($_->{'type'},$pass," :: ",$name,$connector,$function,$suffix);
		} else {
		    $table->add("     &"    ,""   ,""    ,""   ,""        ,$function,$suffix);
		}		
	    }
	} elsif ( exists($_->{'descriptor'}) ) {
	    # Single function specified, and a function descriptor was provided. Take the function name directly from the descriptor.
	    $table->add($_->{'type'},$pass," :: ",$name,$connector,$_->{'descriptor'}->{'name'},"");
	} else {
	    # Single function specified. Simply add to the table.
	    $table->add($_->{'type'},$pass," :: ",$name,$connector,$_->{'function'  }          ,"");
	}
    }
    # Add any descriptions.
    my $description;
    my $methodCount = 0;
    foreach ( @typeBoundFunctions ) {
	my $descriptionText;
	if ( exists($_->{'descriptor'}) ) {
	    $descriptionText = $_->{'descriptor'}->{'description'};
	} elsif ( exists($_->{'description'}) ) {
	    $descriptionText = $_                ->{'description'};
	}
	if ( defined($descriptionText) ) {
	    ++$methodCount;
	    $description .= "     !@  <objectMethod>\n     !@   <method>".$_->{'name'}."</method>\n     !@   <description>".$descriptionText."</description>\n";
	    if ( exists($_->{'descriptor'}) ) {
		# Build type and argument descriptions directly from the function descriptor.
		my $returnType;	
		if ( $_->{'descriptor'}->{'type'} =~ m/^([a-zA-Z0-9_\(\)\s]+)\s+=>\s+([a-zA-Z0-9_]+)/ ) {
		    $returnType = latex_encode($1);
		} else {
		    $returnType = $_->{'descriptor'}->{'type'};
		}
		my @arguments;
		if ( exists($_->{'descriptor'}->{'variables'}) ) {
		    foreach my $declaration ( @{$_->{'descriptor'}->{'variables'}} ) {
			my $intent = join(" ",map {$_ =~ m/intent\(\s*((in|out)+)\s*\)/ ? $1 : ()} @{$declaration->{'attributes'}});
			if ( $intent ne "" ) {
			    my $type = $declaration->{'intrinsic'};
			    $type .= "(".$declaration->{'type'}.")"
				if ( exists($declaration->{'type'}) );
			    $type = "*".$type
				if ( grep {$_ eq "pointer"} @{$declaration->{'attributes'}} );
			    my $dimension = join(",",map {$_ =~ m/dimension\(\s*(.*)\s*\)/ ? $1 : ()} @{$declaration->{'attributes'}});
			    $type .= "[".$dimension."]"
				unless ( $dimension eq "" );
			    my $optional =  grep {$_ eq "optional"} @{$declaration->{'attributes'}};
			    foreach ( @{$declaration->{'variables'}} ) {
				push(
				    @arguments,
				    "\\textcolor{red}{\\textless ".$type."\\textgreater} ".($optional ? "[" : "").$_.($optional ? "]" : "")."\\arg".$intent
				    )
				    unless ( $_ eq "self" );
			    }
			}
		    }
		}
		$description .= "     !@    <type>\\textcolor{red}{\\textless ".$returnType."\\textgreater}</type>\n";
		$description .= "     !@    <arguments>".join(", ",@arguments)."</arguments>\n";
	    } else {
		# Build type and arguments from those directly supplied.
		$description .= "     !@    <type>".$_->{'returnType'}."</type>\n"
		    if ( exists($_->{'returnType'}) );
		$description .= "     !@    <arguments>".$_->{'arguments'}."</arguments>\n"
		    if ( exists($_->{'arguments'}) );
	    }
	    $description .= "     !@  </objectMethod>\n";
	}
    }
    if ( $methodCount == 1 ) {
	$description =~ s/(\!\@\s+\<method\>)/!@   <object>$objectName<\/object>\n     $1/;
    } elsif ( $methodCount > 1 ) {
	$description = "     !@ <objectMethods>\n     !@  <object>".$objectName."</object>\n".$description."     !@ </objectMethods>\n";
    }
    # Construct final product.
    my $product = "";
    $product .= $description
	if ( defined($description) );
    $product .= $table;
    # Return the table.
    return $product;
}

sub derivedTypesSerialize {
    # Generate and insert code for all type definitions.
    my $build = shift();
    # Sort types into dependency order.
    my %typeDependencies;
    foreach ( &List::ExtraUtils::hashList($build->{'types'}) ) {
	# Types are dependent on their parent type.
	push(@{$typeDependencies{$_->{'extends'}}},$_->{'name'})
	    if ( exists($_->{'extends'}) );
	# Types are also dependent on any types used as components, unless that component is a pointer.
	foreach my $dataContent ( @{$_->{'dataContent'}} ) {
	    push(@{$typeDependencies{$dataContent->{'type'}}},$_->{'name'})
		 if
		 (
		  (
		   $dataContent->{'intrinsic'} eq "type"
		   ||
		   $dataContent->{'intrinsic'} eq "class"
		  )
		  && 
		  exists($build->{'types'}->{$dataContent->{'type'}})
		  &&
		  ! grep {$_ eq "pointer"} @{$dataContent->{'attributes'}}
		 );
	}
    }
    my @typeSort  = &List::ExtraUtils::sortedKeys($build->{'types'});
    my @typeOrder =
	toposort
	(
	 sub { @{$typeDependencies{$_[0]} || []}; },
	 \@typeSort
	);
    # Iterate over types.
    foreach ( @typeOrder ) {
	# Get the type.
	my $type = $build->{'types'}->{$_};
	# Insert the type opening.
	$build->{'content'} .= "  type";
	$build->{'content'} .= ", public"
	    if ( exists($type->{'isPublic'}) && $type->{'isPublic'} );
	$build->{'content'} .= ", extends(".$type->{'extends'}.")"
	    if ( exists($type->{'extends'}) );
	$build->{'content'} .= " :: ".$type->{'name'}."\n";
	# Insert any comment.
	$build->{'content'} .= "  !% ".$type->{'comment'}."\n"
	    if ( exists($type->{'comment'}) );
	# Declare contents private.
	$build->{'content'} .= "    private\n";
	# Process any data content.
	$build->{'content'} .= &Fortran::Utils::Format_Variable_Defintions($type->{'dataContent'})
	    if ( exists($type->{'dataContent'}) );
	# Generate and insert a type-bound function table.
	if ( exists($type->{'boundFunctions'}) ) {
	    $build->{'content'} .= "   contains\n";
	    $build->{'content'} .= &boundFunctionTable($type->{'name'},$type->{'boundFunctions'});
	}
	# Insert the type closing.
	$build->{'content'} .= "  end type ".$type->{'name'}."\n\n";
    }
}

sub interfacesSerialize {
    # Generate and insert code for all interfaces.
    my $build = shift();
    # Iterate over interfaces.
    print "   --> Serialize interfaces...\n";
    foreach ( &List::ExtraUtils::hashList($build->{'interfaces'}) ) {
	print "      ---> ".$_->{'name'}."\n";
	$Galacticus::Build::Components::CodeGeneration::interface = $_;
	$build->{'content'} .= 
	    fill_in_string(<<'CODE', PACKAGE => 'Galacticus::Build::Components::CodeGeneration')."\n";
! {$interface->{'comment'}}
abstract interface
  {$interface->{'intrinsic'} eq "void" ? "subroutine" : $interface->{'intrinsic'}." function"} {$interface->{'name'}}({join(",",&Function_Arguments($interface->{'data'}))})
    {&Importables($interface->{'data'}) ? "import ".join(", ",&Importables($interface->{'data'})) : ""}
{&Fortran::Utils::Format_Variable_Defintions($interface->{'data'}, indent => 4)}
  end {$interface->{'intrinsic'} eq "void" ? "subroutine" : "function"} {$interface->{'name'}}
end interface
CODE
    }
    
}

sub functionsSerialize {
    # Serialize function code.
    my $build = shift();
    print "   --> Serialize functions...\n";
    # Iterate over functions.
    foreach my $function 
	(
	 # Include all functions explicitly added to the functions list.
	 @{$build->{'functions'}},
	 # Also include functions which were bound to a derived type. We use two map operations to do this. In the first we
	 # iterate over all derived types, selecting out those with bound functions. For each such derived type, we use a second
	 # map to select out each bound function with an included function descriptor, and return a list of those descriptors.
	 map
	  {exists($_->{'boundFunctions'}) ? map {exists($_->{'descriptor'}) ? $_->{'descriptor'} : ()} @{$_->{'boundFunctions'}} : ()}
	  &List::ExtraUtils::hashList($build->{'types'})
	) {
	# Report.
	print "      --> ".$function->{'name'}."\n";
	# Build function type definition.
	my $form;
	my $type;
	my $result;
	if ( $function->{'type'} =~ m/^([a-zA-Z0-9_\(\),:=\s]+)\s+=>\s+([a-zA-Z0-9_]+)/ ) {
	    my $returnDescriptor = $1;
	    my $returnName       = $2;
	    $result              = "result(".$2.")";
	    $type                = "";
	    $form                = "function";
	    my @attributes;
	    my $returnType;
	    if ( $returnDescriptor =~ m/([a-zA-Z0-9_\(\)\s]+)\s*,\s*([a-zA-Z0-9_,\(\):=\s]+)/ ) {
		$returnType = $1;
		@attributes = split(/\s*,\s*/,$2);
	    } else {
		$returnType = $returnDescriptor;
	    }
	    my $declaration      =
	    {
		attributes => \@attributes,
		variables  => [ $returnName ]
	    };
	    if ( $returnType =~ m/^(type|class)\s*\(\s*([a-zA-Z0-9_]+)\s*\)/ ) {
		$declaration->{'intrinsic'} = $1;
		$declaration->{'type'     } = $2;
	    } else {
		$declaration->{'intrinsic'} = $returnType;
	    }
	    push(@{$function->{'variables'}},$declaration);
	} else {
	    $form   = $function->{'type'} eq "void" ? "subroutine" : "function";
	    $type   = $function->{'type'} eq "void" ? ""           : $function->{'type'};
	    $result = "";
	}
	my @functionAttributes;
	push(@functionAttributes,"recursive")
	    if ( exists($function->{'recursive'}) && $function->{'recursive'} );
	# Build a list of arguments.
	my @arguments;
	foreach my $variables ( @{$function->{'variables'}} ) {
	    push(
		@arguments,
		@{$variables->{'variables'}}
		)
		if (
		    # Regular variables with intent.
		    (
		     exists($variables->{'attributes'})
		     &&
		     grep {$_ =~ m/^intent\s*\(\s*(in|inout|out)\s*\)/} @{$variables->{'attributes'}}		      
		    )
		    ||
		    # Procedure pointers
		    (
		     $variables->{'intrinsic'} eq "procedure"
		     &&
		     exists($variables->{'attributes'})
		     &&
		     grep {$_ eq "pointer"} @{$variables->{'attributes'}}		      
		    )
		    ||
		    # External variables.
		    $variables->{'intrinsic'} eq "external"
		);
	}

	# Serialize function opener.
	$build->{'content'} .= join(" ",@functionAttributes)." ".$type." ".$form." ".$function->{'name'}."(".join(",",@arguments).") ".$result."\n";
	# Serialize description.
	$build->{'content'} .= "   !% ".$function->{'description'}."\n";
	# Serialize module uses.
	$build->{'content'} .= "   use ".$_."\n"
	    foreach ( @{$function->{'modules'}} );
	# Serialize variable definitions.
	$build->{'content'} .= "   implicit none\n";
	$build->{'content'} .= &Fortran::Utils::Format_Variable_Defintions($function->{'variables'})
	    if ( exists($function->{'variables'}) );
	# Serialize content.
	$build->{'content'} .= $function->{'content'};
	# Serialize function closer.
	$build->{'content'} .= "   return\n";
	$build->{'content'} .= "end ".$form." ".$function->{'name'}."\n\n";
    }
}

1;
