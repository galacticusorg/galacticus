# Contains a Perl module which implements processing of "component" directives in the Galacticus build system.

package Component;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V094"}) ) {
    $galacticusPath = $ENV{"GALACTICUS_ROOT_V094"};
    $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
    $galacticusPath = "./";
}
unshift(@INC, $galacticusPath."perl"); 
use strict;
use warnings;
use utf8;
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
require File::Changes;
require Fortran::Utils;
require Galacticus::Build::Hooks;
require Galacticus::Build::Components::Utils;
require Galacticus::Build::Components::CodeGeneration;
require Galacticus::Build::Components::Hierarchy;
require Galacticus::Build::Components::TreeNodes;
require Galacticus::Build::Components::NodeEvents;
require Galacticus::Build::Components::BaseTypes;
require Galacticus::Build::Components::Classes;
require Galacticus::Build::Components::Implementations;
require Galacticus::Build::Components::Properties;
require Galacticus::Build::Components::Properties::Set;
require Galacticus::Build::Components::Attributes;
require Galacticus::Build::Components::DataTypes;
$SIG{ __DIE__ } = sub { Carp::confess( @_ ) };

# Insert hooks for our functions.
%Hooks::moduleHooks = 
    (
     %Hooks::moduleHooks,
     component => 
     {
	 parse    => \&Components_Parse_Directive,
	 generate => \&Components_Generate_Output,
	 validate => \&Components_Validate
     }
    );

# Include debugging code.
my $debugging                           = 0;

# Switch to control gfortran workarounds.
my $workaround                         = 1;

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
    my $validator = XML::Validator::Schema->new(file => $galacticusPath."schema/componentSchema.xsd");
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
    my $build = shift;

    # Sort hooks.
    my @hooks = map
    {{name => $_, hook => $Galacticus::Build::Component::Utils::componentUtils{$_}}}
    &ExtraUtils::sortedKeys(\%Galacticus::Build::Component::Utils::componentUtils);

    # Iterate over phases.
    print "--> Phase:\n";
    foreach my $phase ( "preValidate", "default", "gather", "scatter", "postValidate", "content", "types", "functions" ) {
	print "   --> ".ucfirst($phase)."...\n";
	foreach my $hook ( @hooks ) {	
	    if ( exists($hook->{'hook'}->{$phase}) ) {
		foreach my $function ( &ExtraUtils::as_array($hook->{'hook'}->{$phase}) ) {
		    print "      --> ".$hook->{'name'}."\n";
		    &{$function}($build);
		}
	    }
	}
    }
    # Iterate over all functions, calling them with the build data object.
    &{$_}($build)
	foreach (
	    # Generate a finalization method.
	    \&Generate_Finalization_Function                         ,
	    # Generate functions to map other functions over components.
	    \&Generate_Map_Functions                                 ,
	    # Generate functions to dump nodes.
	    \&Generate_Node_Dump_Function                            ,
	    # Generate functions to output nodes.
	    \&Generate_Node_Output_Functions                         ,
	    # Generate functions to serialize/deserialize nodes to/from arrays.
	    \&Generate_Node_Serialization_Functions                  ,
	    # Generate functions compute offsets into serialization arrays.
	    \&Generate_Node_Offset_Functions                         ,
	    # Generate functions to get property names from a supplied index.
	    \&Generate_Node_Property_Name_From_Index_Function        ,
	    # Generate function to copy one node to another.
	    \&Generate_Node_Copy_Function                            ,
	    # Generate function to move one node to another.
	    \&Generate_Node_Move_Function                            ,
	    # Generate a tree node creation function.
	    \&Generate_Tree_Node_Creation_Function                   ,
	    # Generate a tree node destruction function.
	    \&Generate_Tree_Node_Destruction_Function                ,
	    # Generate a tree node builder function.
	    \&Generate_Tree_Node_Builder_Function                    ,
	    # Generate type name functions.
	    \&Generate_Type_Name_Functions                           ,
	    # Generate component assignment function.
	    \&Generate_Component_Assignment_Function                 ,
	    # Generate component class destruction functions.
	    \&Generate_Component_Class_Destruction_Functions         ,
	    # Generate component class removal functions.
	    \&Generate_Component_Class_Removal_Functions             ,
	    # Generate component class move functions.
	    \&Generate_Component_Class_Move_Functions                ,
	    # Generate dump functions for each component class.
	    \&Generate_Component_Class_Dump_Functions                ,
	    # Generate builder functions for each component class.
	    \&Generate_Component_Class_Builder_Functions             ,
	    # Generate initializor functions for each component class.
	    \&Generate_Component_Class_Initializor_Functions         ,
	    # Generate output functions for each component class.
	    \&Generate_Component_Class_Output_Functions              ,
	    # Generate functions to determine if an implementation is active.
	    \&Generate_Is_Active_Functions                           ,
	    # Generate component implementation destruction functions.
	    \&Generate_Component_Implementation_Destruction_Functions,
	    # Generate ODE solver initialization functions.
	    \&Generate_Node_ODE_Initialization_Functions             ,
	    # Generate dump functions for each implementation.
	    \&Generate_Implementation_Dump_Functions                 ,
	    # Generate initializor functions for each implementation.
	    \&Generate_Implementation_Initializor_Functions          ,
	    # Generate builder functions for each implementation.
	    \&Generate_Implementation_Builder_Functions              ,
	    # Generate output functions for each implementation.
	    \&Generate_Implementation_Output_Functions               ,
	    # Generate functions to get the name of properties from an index.
	    \&Generate_Implementation_Name_From_Index_Functions      ,
	    # Generate serialization/deserialization functions for each implementation.
	    \&Generate_Implementation_Serialization_Functions        ,
	    # Generate serialization offset functions for each implementation.
	    \&Generate_Implementation_Offset_Functions               ,
	    # Generate component count methods.
	    \&Generate_Component_Count_Functions                     ,
	    # Generate component get methods.
	    \&Generate_Component_Get_Functions                       ,
	    # Generate component destruction functions.
	    \&Generate_Component_Destruction_Functions               ,
	    # Generate component creation functions.
	    \&Generate_Component_Creation_Functions                  ,
	    # Generate component class default value functions.
	    \&Generate_Component_Class_Default_Value_Functions       ,
	    # Generate functions for getting/setting/rating value via a deferred function.
	    \&Generate_Deferred_GSR_Function                         ,
	    # Generate functions for getting/setting/rating value directly.
	    \&Generate_GSR_Functions                                 ,
	    # Generate functions for returning which components support getting/setting/rating.
	    \&Generate_GSR_Availability_Functions                    ,
	    # Insert code for type-definitions.
	    \&Insert_Type_Definitions                                ,
	    # Generate module status.
	    \&Generate_Initialization_Status                         ,
	    # Generate objects that record which type of component will be used by default.
	    \&Generate_Default_Component_Sources                     ,
	    # Generate records of which component implementations are selected.
	    \&Generate_Active_Implementation_Records                 ,
	    # Generate variables that record offsets for serialization.
	    \&Generate_Serialization_Offset_Variables                ,
	    # Generate deferred procedure pointers.
	    \&Generate_Deferred_Procedure_Pointers                   ,
	    # Generate deferred binding procedure pointers.
	    \&Generate_Deferred_Binding_Procedure_Pointers           ,
	    # Generate deferred binding procedure pointers.
	    \&Generate_Deferred_Binding_Functions                    ,
	    # Generate required null binding functions.
	    \&Generate_Null_Binding_Functions                        ,
	    # Generate all interfaces.
	    \&Generate_Interfaces                                    ,
	);

    # Insert all module scope variables.
    $build->{'content'} .= &Fortran_Utils::Format_Variable_Defintions($build->{'variables'})."\n";

    # Insert the "contains" line.
    $build->{'content'} .= "contains\n\n";

    # Serialize all functions.
    &functionsSerialize($build);
    
    # Insert all functions into content.
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
    &File_Changes::Update($ENV{'BUILDPATH'}."/Makefile_Component_Includes" ,$ENV{'BUILDPATH'}."/Makefile_Component_Includes.tmp" );

}

sub Get_Type {
    # Returns the type of a method of pipe.
    my $build->{'currentDocument'} = shift;
    # Assume scalar type by default
    my $type = "scalar";
    # If a type is specified, then return it instead.
    if ( exists($build->{'currentDocument'}->{'type'}) ) {$type = $build->{'currentDocument'}->{'type'}};
    return $type;
}

sub pad {
    # Pad a string to give nicely aligned formatting in the output code.
    die("pad() requires two arguments")
	unless (scalar(@_) == 2);
    my $text       = shift;
    my $padLength  = $_[0];

    die("pad(): text is too long to pad: '".$text."'")
	if ( length($text) > $padLength );
    
    my $paddedText = $text." " x ($padLength-length($text));
    return $paddedText;
}

sub Generate_Implementations {
    # Generate a type for each component implementation.
    my $build = shift;
    # Create classes for each specific implementation.
    foreach my $componentID ( @{$build->{'componentIdList'}} ) {
	# Get the implementation.
	my $component          = $build->{'components'}->{$componentID};
	# Get the parent class.
	my $componentClassName = $component->{'class'};
    	# Determine the name of the class which this component extends (use the "nodeComponent" class by default).
    	my $extensionOf = 
	    exists($component->{'extends'})
	    ?
	    "nodeComponent".ucfirst($component->{'extends'}->{'class'}).ucfirst($component->{'extends'}->{'name'})
	    : 
	    "nodeComponent".ucfirst($componentClassName);
     	# Create data objects to store all of the linked data for this component.
	my @dataContent;
    	foreach ( &ExtraUtils::sortedKeys($component->{'content'}->{'data'}) ) {
    	    my $type = &DataTypes::dataObjectName($component->{'content'}->{'data'}->{$_});
	    (my $typeDefinition, my $typeLabel) = &DataTypes::dataObjectDefinition($component->{'content'}->{'data'}->{$_});
	    $typeDefinition->{'variables'} = [ $_ ];
	    push(
		@dataContent,
		$typeDefinition
		);
    	}
	# Create a list for type-bound functions.
	my @typeBoundFunctions;
	# Add binding for deferred create function set function.
	push(
	    @typeBoundFunctions,
	    {type => "procedure", pass => "nopass", name => "createFunctionSet", function => $componentID."CreateFunctionSet", description => "Set the function used to create {\\normalfont \\ttfamily ".$componentID."} components.", returnType => "\\void", arguments => "\\textcolor{red}{\\textless function()\\textgreater}"}
	    )
	    if ( 
		exists($component->{'createFunction'})
		&&
		$component->{'createFunction'}->{'isDeferred'}
	    );
     	# If this component has bindings defined, scan through them and create an appropriate method.
    	if ( exists($component->{'bindings'}) ) {
    	    foreach ( @{$component->{'bindings'}->{'binding'}} ) {
		my %function = (
		    type        => "procedure",
		    name        => $_->{'method'}
		    );
		if ( ! $_->{'isDeferred'} ) {
		    # Binding is not deferred, simply map to the given function.
		    $function{'function'} = $_->{'function'};
		} else {
		    # Binding is deferred, map to a suitable wrapper function.
		    $function{'function'   } = $componentID.$_->{'method'};
		    $function{'returnType' } = &DataTypes::dataObjectDocName($_->{'interface'});
		    $function{'arguments'  } = "";
		    $function{'description'} = "Get the {\\normalfont \\ttfamily ".$_->{'method'}."} property of the {\\normalfont \\ttfamily ". $componentID."} component.";
		    # Also add bindings to functions to set and test the deferred function.
		    my %setFunction = (
			type        => "procedure",
			pass        => "nopass",
			name        => $_->{'method'}."Function",
			function    => $componentID.$_->{'method'}."DeferredFunctionSet",
			returnType  => "\\void",
			arguments   => "procedure(".$componentID.$_->{'method'}."Interface) deferredFunction",
			description => "Set the function for the deferred {\\normalfont \\ttfamily ".$_->{'method'}."} propert of the {\\normalfont \\ttfamily ". $componentID."} component."
			);
		    my %testFunction = (
			type        => "procedure",
			pass        => "nopass",
			name        => $_->{'method'}."FunctionIsSet",
			function    => $componentID.$_->{'method'}."DfrrdFnctnIsSet",
			returnType  => "\\logicalzero",
			arguments   => "",
			description => "Specify whether the deferred function for the {\\normalfont \\ttfamily ".$_->{'method'}."} property of the {\\normalfont \\ttfamily ". $componentID."} component has been set."
			);
		    push(@typeBoundFunctions,\%setFunction,\%testFunction);
		}
		foreach my $attribute ( "description", "returnType", "arguments" ) {
		    $function{$attribute} = $_->{$attribute}
		    if ( exists($_->{$attribute}) );
		}
		push(@typeBoundFunctions,\%function);
	    }
	}
	# Iterate over properties.
	foreach my $propertyName ( &ExtraUtils::sortedKeys($component->{'properties'}->{'property'}) ) {
	    # Get the property.
	    my $property = $component->{'properties'}->{'property'}->{$propertyName};
	    push(
		@typeBoundFunctions,
		{type => "procedure", pass => "nopass", name => $propertyName."IsGettable", function => "Boolean_".ucfirst($Utils::booleanLabel[$property->{'attributes'}->{'isGettable'}])},
		{type => "procedure", pass => "nopass", name => $propertyName."IsSettable", function => "Boolean_".ucfirst($Utils::booleanLabel[$property->{'attributes'}->{'isSettable'}])}
		);
	}
	# Create the type.
	$build->{'types'}->{'nodeComponent'.ucfirst($componentID)} = {
	    name           => "nodeComponent".ucfirst($componentID),
	    comment        => "Class for the ".$component->{'name'}." implementation of the ".$componentClassName." component.",
	    isPublic       => 1,
	    extends        => $extensionOf,
	    boundFunctions => \@typeBoundFunctions,
	    dataContent    => \@dataContent
	};
    }
}

sub Generate_Active_Implementation_Records{
    # Generate records of which component implementations are selected.
    my $build = shift;
    # Create a table.
    my $recordTable = Text::Table->new(
	{
	    is_sep => 1,
	    body   => "  logical :: "
	},
	{
	    align  => "left"
	},
	{
	    is_sep => 1,
	    body   => "=.false."
	}
	);
    # Iterate over all component implementations.
    foreach ( @{$build->{'componentIdList'}} ) {
	$recordTable->add("nodeComponent".$_."IsActive");
    }
    # Insert into the document.
    $build->{'content'} .= "  ! Records of which component implementations are active.\n";
    $build->{'content'} .= $recordTable->table()."\n";
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
    $build->{'content'} .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent, indent => 2)."\n";
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
		($type, my $name, my $attributeList) = &DataTypes::dataObjectPrimitiveName($binding->{'interface'})
		    unless ( $type eq "void" );
		my $endType;
		if ( $type eq "void" ) {
		    $type    = "subroutine";
		    $endType = "subroutine";
		} else {
		    $type    .= " function";
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
				&Fortran_Utils::Extract_Variables($1)
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
		unless ( exists($interfaces{$interfaceName}) ) {
		    $build->{'content'} .= "abstract interface\n";
		    $build->{'content'} .= "  ".$type." ".$interfaceName."Interface(".join(",",@arguments).")\n";
		    $build->{'content'} .= "    import nodeComponent".$highLevel."\n"
			if ( $binding->{'interface'}->{'self'}->{'pass'} eq "true" );
		    $build->{'content'} .= "    class(nodeComponent".$highLevel."), intent(".$binding->{'interface'}->{'self'}->{'intent'}.") :: self\n"
			if ( $binding->{'interface'}->{'self'}->{'pass'} eq "true" );
		    $build->{'content'} .= "    ".$_."\n"
			foreach ( @{$binding->{'interface'}->{'argument'}} );
		    $build->{'content'} .= "  end ".$endType." ".$interfaceName."Interface\n";
		    $build->{'content'} .= "end interface\n\n";
		    $interfaces{$interfaceName} = 1;
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
		    $functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
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
		    $functionCode  = "  ".$type." ".$classFunctionName."(".join(",",@arguments).")\n";
		    $functionCode .= "    !% Call the deferred function for the {\\normalfont \\ttfamily ".$binding->{'method'}."} method of the {\\normalfont \\ttfamily ".$component->{'class'}."} component class.\n";
		    $functionCode .= "    use Galacticus_Error\n";
		    $functionCode .= "    implicit none\n";
		    $functionCode .= "    class(nodeComponent".ucfirst($component->{'class'})."), intent(".$binding->{'interface'}->{'self'}->{'intent'}.") :: self\n";
		    $functionCode .= "    ".$_."\n"
			foreach ( @{$binding->{'interface'}->{'argument'}} );
		    $functionCode .= "    if (self%".$binding->{'method'}."FunctionIsSet()) then\n";
		    if ( $type eq "subroutine" ) {
			$functionCode .= "       call ".$classFunctionName."Deferred(".join(",",@selfishArguments).")\n";
		    } else {
			$functionCode .= "       ".$classFunctionName."=".$classFunctionName."Deferred(".join(",",@selfishArguments).")\n";
		    }
		    $functionCode .= "    else\n";
		    if ( $type eq "double precision" ) {
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
		$functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
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
		$functionCode  = "  ".$type." ".$componentFunctionName."(".join(",",@arguments).")\n";
		$functionCode .= "    !% Call the deferred function for the {\\normalfont \\ttfamily ".$binding->{'method'}."} method of the {\\normalfont \\ttfamily ".$componentID."} component.\n";
		$functionCode .= "    use Galacticus_Error\n";
		$functionCode .= "    implicit none\n";
		$functionCode .= "    class(nodeComponent".ucfirst($componentID)."), intent(".$binding->{'interface'}->{'self'}->{'intent'}.") :: self\n";
		$functionCode .= "    ".$_."\n"
		    foreach ( @{$binding->{'interface'}->{'argument'}} );
		$functionCode .= "    select type (self)\n";
		$functionCode .= "    class is (nodeComponent".ucfirst($componentID).")\n";
		$functionCode .= "       if (self%".$binding->{'method'}."FunctionIsSet()) then\n";
		if ( $type eq "subroutine" ) {
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
		    if ( $type eq "subroutine" ) {
			$functionCode .= "          call self%".$parentType."%".$binding->{'method'}."(".join(",",@selflessArguments).")\n";
		    } else {
			$functionCode .= "          ".$componentFunctionName."=self%".$parentType."%".$binding->{'method'}."(".join(",",@selflessArguments).")\n";
		    }
		} else {
		    if ( $type eq "double precision function" ) {
			$functionCode .= "       ".$componentFunctionName."=0.0d0\n";
		    }
		    $functionCode .= "          call Galacticus_Error_Report('".$componentFunctionName."','deferred function has not been assigned')\n";
		}
		$functionCode .= "       end if\n";
		$functionCode .= "    class default\n";
		if ( $type eq "double precision function" ) {
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
	# Create pointer for deferred create functions.
	push(
	    @dataContent,
	    {
		intrinsic  => "procedure",
		type       => "",
		attributes => [ "pointer" ],
		variables  => [ $componentID."CreateFunction" ]
	    }
	    )
	    if (
		exists($component->{'createFunction'})
		&&
		$component->{'createFunction'}->{'isDeferred'}
	    );
	# Iterate over properties.
	foreach my $propertyName ( &ExtraUtils::sortedKeys($component->{'properties'}->{'property'}) ) {
	    my $property = $component->{'properties'}->{'property'}->{$propertyName};
	    unless ( $property->{'attributes' }->{'isDeferred'} eq "" ) {
		my $selfType = "generic";
		$selfType = $component->{'class'}
		   unless ( $property->{'attributes'}->{'bindsTo'} eq "top" );
		(my $dataObject, my $label) = &DataTypes::dataObjectDefinition($property);
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
    $build->{'content'} .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent, indent => 2)."\n";
}

sub Generate_Node_Event_Interface {
    # Generate interace for node event tasks.
    my $build = shift;
    # Initialize data content.
    @code::dataContent = 
	(
	 {
	     intrinsic  => "class",
	     type       => "nodeEvent",
	     attributes => [ "intent(in   )" ],
	     variables  => [ "thisEvent" ]
	 },
	 {
	     intrinsic  => "type",
	     type       => "treeNode",
	     attributes => [ "pointer", "intent(inout)" ],
	     variables  => [ "thisNode" ]
	 },
	 {
	     intrinsic  => "integer",
	     attributes => [ "intent(inout)" ],
	     variables  => [ "deadlockStatus" ]
	 }
	);
    # Insert interface.
    $build->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
! Interface for node event tasks.
abstract interface
  logical function nodeEventTask(thisEvent,thisNode,deadlockStatus)
    import nodeEvent,treeNode
{&Fortran_Utils::Format_Variable_Defintions(\@dataContent, indent => 4)}
  end function nodeEventTask
end interface
CODE
}

sub Generate_Default_Component_Sources{
    # Generate records of which component implementations are selected.
    my $build = shift;
    # Create default objects for each class.
    $build->{'content'} .= "  ! Objects that will record which type of each component is to be used by default.\n";
    $build->{'content'} .= &Fortran_Utils::Format_Variable_Defintions
	(
	 [
	  map
	  {
	      {
		  intrinsic  => "class"                              ,
		  type       => "nodeComponent".ucfirst($_)          ,
		  attributes => [ "public", "allocatable" ]          ,
		  variables  => [ "default".ucfirst($_)."Component" ]
	      }
	  } 
	  @{$build->{'componentClassList'}}
	 ],
	 indent => 2
	);
    # Create source objects for each class.
    $build->{'content'} .= "  ! Objects used to allocate components of given class..\n";
    $build->{'content'} .= &Fortran_Utils::Format_Variable_Defintions
	(
	 [
	  map
	  {
	      {
		  intrinsic  => "type"                               ,
		  type       => "nodeComponent".ucfirst($_)          ,
		  attributes => [ ]                                  ,
		  variables  => [ $_."Class" ]
	      }
	  } 
	  @{$build->{'componentClassList'}}
	 ],
	 indent => 2
	);
}

sub Generate_Initialization_Status {
    # Generate a variable that stores initialization status..
    my $build = shift;
    # Insert into the document.
    $build->{'content'} .= "  ! Record of module initialization status.\n";
    $build->{'content'} .= "  logical :: moduleIsInitialized=.false.\n";
}

sub Generate_Finalization_Function {
    # Generate a finalization function.
    my $build = shift;
    # Create a table for the deallocation code.
    my $table = Text::Table->new(
	{
	    is_sep => 1,
	    body   => "    deallocate(default"
	},
	{
	    align  => "left"
	},
	{
	    is_sep => 1,
	    body   => ")"
	},
	);
    # Populate the table.
    foreach ( @{$build->{'componentClassList'}} ) {
	$table->add(ucfirst($_)."Component");
    }
    # Generate the function code.
    my $functionCode;
    $functionCode .= "  subroutine Galacticus_Nodes_Finalize()\n";
    $functionCode .= "    !% Finialize the \\glc\\ object system.\n";
    $functionCode .= "    implicit none\n\n";
    $functionCode .= "    if (.not.moduleIsInitialized) return\n";
    $functionCode .= $table->table();
    $functionCode .= "    return\n";
    $functionCode .= "  end subroutine Galacticus_Nodes_Finalize\n";
    # Insert into the function list.
    push(
	@{$build->{'code'}->{'functions'}},
	$functionCode
	);
}

sub Generate_Map_Functions {
    # Generate functions to map other functions over components.
    my $build = shift;

    # Function for mapping a void function.
    # Generate variables.
    my @dataContent =
	(
	 {
	     intrinsic  => "class",
	     type       => "treeNode",
	     attributes => [ "intent(inout)" ],
	     variables  => [ "self" ]
	 },
	 {
	     intrinsic  => "procedure",
	     type       => "Node_Component_Null_Void0_InOut",
	     attributes => [ "pointer" ],
	     variables  => [ "mapFunction" ]
	 },
	 {
	     intrinsic  => "integer",
	     variables  => [ "i" ]
	 }
	);
    my $functionCode;
    $functionCode .= "  subroutine mapComponentsVoid(self,mapFunction)\n";
    $functionCode .= "    !% Map a void function over components.\n";
    $functionCode .= "    implicit none\n";
    $functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
    foreach ( @{$build->{'componentClassList'}} ) {	    
     	$functionCode .= "    if (allocated(self%component".&Utils::padClass(ucfirst($_),[19,0]).")) then\n";
	$functionCode .= "      do i=1,size(self%component".&Utils::padClass(ucfirst($_),[19,0]).")\n";
	$functionCode .= "        call mapFunction(self%component".&Utils::padClass(ucfirst($_),[19,0])."(i))\n";
	$functionCode .= "      end do\n";
	$functionCode .= "    end if\n";
    }
    $functionCode .= "    return\n";
    $functionCode .= "  end subroutine mapComponentsVoid\n\n";
    # Insert into the function list.
    push(
	@{$build->{'code'}->{'functions'}},
	$functionCode
	);
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{type => "procedure", name => "mapVoid", function => "mapComponentsVoid", description => "Map a void function over components.", returnType => "\\void", arguments => "\\textcolor{red}{\\textless *function()\\textgreater} mapFunction"}
	);

    # Function for mapping a scalar double function.
    @dataContent =
	(
	 {
	     intrinsic  => "class",
	     type       => "treeNode",
	     attributes => [ "intent(inout)" ],
	     variables  => [ "self" ]
	 },
	 {
	     intrinsic  => "procedure",
	     type       => "Node_Component_Null_Double0_InOut",
	     attributes => [ "pointer" ],
	     variables  => [ "mapFunction" ]
	 },
	 {
	     intrinsic  => "integer",
	     attributes => [ "intent(in   )" ],
	     variables  => [ "reduction" ]
	 },
	 {
	     intrinsic  => "integer",
	     attributes => [ "intent(in   )", "optional" ],
	     variables  => [ "optimizeFor" ]
	 },
	 {
	     intrinsic  => "double precision",
	     variables  => [ "componentValue" ]
	 },
	 {
	     intrinsic  => "integer",
	     variables  => [ "i" ]
	 }
	);
    $functionCode  = "  double precision function mapComponentsDouble0(self,mapFunction,reduction,optimizeFor)\n";
    $functionCode .= "    !% Map a scalar double function over components with a specified {\\normalfont \\ttfamily reduction}.\n";
    $functionCode .= "    use Galacticus_Error\n";
    $functionCode .= "    implicit none\n";
    $functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
    # Scan through available node component methods and find ones which are mappable. Create optimized versions of this function
    # for them.
    my $optimizationLabel      = -1;
    my $optimizationsGenerated = 0;
    foreach my $boundFunction ( @{$build->{'types'}->{'nodeComponent'}->{'boundFunctions'}} ) {
	if ( exists($boundFunction->{'mappable'}) ) {
	    my @reductions = split(/:/,$boundFunction->{'mappable'});
	    foreach my $reduction ( @reductions ) {
		# Record that optimized versions were generated.
		$optimizationsGenerated = 1;
		# Insert test for optimized case.
		++$optimizationLabel;
		$build->{'content'} .= "   integer, public, parameter :: optimizeFor".ucfirst($boundFunction->{'name'}).ucfirst($reduction)."=".$optimizationLabel."\n";
		$functionCode .= "   ";
		$functionCode .= "else"
		    unless ( $optimizationLabel == 0 );
		$functionCode .= "if (present(optimizeFor).and.optimizeFor == optimizeFor".ucfirst($boundFunction->{'name'}).ucfirst($reduction).") then\n";
		# Initialize reduction.
		if ( $reduction eq "summation" ) {
		    $functionCode .= "      if (reduction /= reductionSummation) call Galacticus_Error_Report('mapComponentsDouble0','reduction mismatch')\n";
		    $functionCode .= "      mapComponentsDouble0=0.0d0\n";
		} elsif ( $reduction eq "product" ) {
		    $functionCode .= "      if (reduction /= reductionProduct  ) call Galacticus_Error_Report('mapComponentsDouble0','reduction mismatch')\n";
		    $functionCode .= "      mapComponentsDouble0=1.0d0\n";
		} else {
		    die("Generate_Map_Functions(): unrecognized reduction");
		}
		# Iterate over available types.
		foreach my $type ( &ExtraUtils::sortedKeys($build->{'types'}) ) {
		    if ( $type =~ m/^nodeComponent.+/  && grep {$_->{'name'} eq $boundFunction->{'name'}} @{$build->{'types'}->{$type}->{'boundFunctions'}} ) {
			# Determine the class of this component.
			my $baseClass = $type;
			while ( exists($build->{'types'}->{$baseClass}->{'extends'}) && $build->{'types'}->{$baseClass}->{'extends'} ne "nodeComponent" ) {
			    $baseClass = $build->{'types'}->{$baseClass}->{'extends'};
			}
			$baseClass =~ s/^nodeComponent//;
			$baseClass = lc($baseClass);
			# Construct code for this component.
			$functionCode .= "    if (allocated(self%component".&Utils::padClass(ucfirst($baseClass),[0,0]).")) then\n";
			$functionCode .= "      select type (c => self%component".&Utils::padClass(ucfirst($baseClass),[0,0]).")\n";
			$functionCode .= "      type is (".$type.")\n";
			$functionCode .= "         do i=1,size(self%component".&Utils::padClass(ucfirst($baseClass),[0,0]).")\n";
			$functionCode .= "            mapComponentsDouble0=mapComponentsDouble0";
			if ( $reduction eq "summation" ) {
			    $functionCode .= "+";
			} elsif ( $reduction eq "product" ) {
			    $functionCode .= "*";
			}
			$functionCode .= "mapFunction(self%component".&Utils::padClass(ucfirst($baseClass),[0,0])."(i))\n";
			$functionCode .= "         end do\n";
			$functionCode .= "      end select\n";
			$functionCode .= "    end if\n";
		    }
		}
	    }
	}
    }
    $build->{'content'} .= "\n";
    # Generate the generic, unoptimized function.
    $functionCode .= "    else\n"
	if ( $optimizationsGenerated == 1 );
    $functionCode .= "    select case (reduction)\n";
    $functionCode .= "    case (reductionSummation)\n";
    $functionCode .= "      mapComponentsDouble0=0.0d0\n";
    $functionCode .= "    case (reductionProduct  )\n";
    $functionCode .= "      mapComponentsDouble0=1.0d0\n";
    $functionCode .= "    case default\n";
    $functionCode .= "      mapComponentsDouble0=1.0d0\n";
    $functionCode .= "      call Galacticus_Error_Report('mapComponentsDouble0','unknown reduction')\n";
    $functionCode .= "    end select\n";
    foreach ( @{$build->{'componentClassList'}} ) {	    
     	$functionCode .= "    if (allocated(self%component".&Utils::padClass(ucfirst($_),[0,0]).")) then\n";
	$functionCode .= "      do i=1,size(self%component".&Utils::padClass(ucfirst($_),[0,0]).")\n";
     	$functionCode .= "        componentValue=mapFunction(self%component".&Utils::padClass(ucfirst($_),[0,0])."(i))\n";
     	$functionCode .= "        select case (reduction)\n";
     	$functionCode .= "        case (reductionSummation)\n";
     	$functionCode .= "          mapComponentsDouble0=mapComponentsDouble0+componentValue\n";
     	$functionCode .= "        case (reductionProduct  )\n";
     	$functionCode .= "          mapComponentsDouble0=mapComponentsDouble0*componentValue\n";
     	$functionCode .= "        end select\n";
	$functionCode .= "      end do\n";
     	$functionCode .= "    end if\n";
    }
    $functionCode .= "    end if\n"
	if ( $optimizationsGenerated == 1 );
    $functionCode .= "    return\n";
    $functionCode .= "  end function mapComponentsDouble0\n\n";
    # Insert into the function list.
    push(
	@{$build->{'code'}->{'functions'}},
	$functionCode
	);
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{type => "procedure", name => "mapDouble0", function => "mapComponentsDouble0", description => "Map a scalar double function over components.", returnType => "\\doublezero", arguments => "\\textcolor{red}{\\textless *function()\\textgreater} mapFunction"}
	);
}

sub Generate_Node_Dump_Function {
    # Generate function to dump node properties.
    my $build = shift;
    # Create the function.
    my @dataContent =
	(
	 {
	     intrinsic  => "class",
	     type       => "treeNode",
	     attributes => [ "intent(in   )" ],
	     variables  => [ "self" ]
	 },
	 {
	     intrinsic  => "integer",
	     variables  => [ "i" ]
	 },
	 {
	     intrinsic  => "type",
	     type       => "varying_string",
	     variables  => [ "message" ]
	 }
	);
    my $functionCode;
    $functionCode  = "  subroutine Node_Dump(self)\n";
    $functionCode .= "    !% Dump node content.\n";
    $functionCode .= "    use ISO_Varying_String\n";
    $functionCode .= "    use Galacticus_Display\n";
    $functionCode .= "    use String_Handling\n";
    $functionCode .= "    implicit none\n";
    $functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
    # Iterate over pointers.
    $functionCode .= "    message='Dumping node '\n";
    $functionCode .= "    message=message//self%index()\n";
    $functionCode .= "    call Galacticus_Display_Indent(message)\n";
    $functionCode .= "    message='host tree: '\n";
    $functionCode .= "    message=message//self%hostTree%index\n";
    $functionCode .= "    call Galacticus_Display_Message(message)\n";
    $functionCode .= "    call Galacticus_Display_Indent('pointers')\n";
    foreach my $pointer ( "parent", "firstChild", "sibling", "firstSatellite", "mergeTarget", "firstMergee", "siblingMergee", "formationNode" ) {
	$functionCode .= "   message='".(" " x (14-length($pointer))).$pointer.": '\n";
	$functionCode .= "   message=message//self%".pad($pointer,14)."%index()\n";
	$functionCode .= "   call Galacticus_Display_Message(message)\n";
    }
    $functionCode .= "   call Galacticus_Display_Unindent('done')\n";
    # Iterate over all component classes
    foreach ( @{$build->{'componentClassList'}} ) {	    
	$functionCode .= "    if (allocated(self%component".&Utils::padClass(ucfirst($_),[0,0]).")) then\n";
	$functionCode .= "      do i=1,size(self%component".&Utils::padClass(ucfirst($_),[0,0]).")\n";
	$functionCode .= "        call self%component".&Utils::padClass(ucfirst($_),[0,0])."(i)%dump()\n";
	$functionCode .= "      end do\n";
	$functionCode .= "    end if\n";
    }
    $functionCode .= "    call Galacticus_Display_Unindent('done')\n";
    $functionCode .= "    return\n";
    $functionCode .= "  end subroutine Node_Dump\n";
    # Insert into the function list.
    push(
	@{$build->{'code'}->{'functions'}},
	$functionCode
	);
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{type => "procedure", name => "dump", function => "Node_Dump", description => "Generate an ASCII dump of all content of a node.", returnType => "\\void", arguments => ""}
	);
    # Create the function.
    @dataContent =
	(
	 {
	     intrinsic  => "class",
	     type       => "treeNode",
	     attributes => [ "intent(inout)" ],
	     variables  => [ "self" ]
	 },
	 {
	     intrinsic  => "integer",
	     attributes => [ "intent(in   )" ],
	     variables  => [ "fileHandle" ]
	 },
	 {
	     intrinsic  => "integer",
	     variables  => [ "i" ]
	 },
	 {
	     intrinsic  => "character",
	     type       => "len=20",
	     variables  => [ "idLabel", "treeLabel" ]
	 }
	);
    $functionCode  = "  subroutine Node_Dump_XML(self,fileHandle)\n";
    $functionCode .= "    !% Dump node content.\n";
    $functionCode .= "    use ISO_Varying_String\n";
    $functionCode .= "    use Galacticus_Display\n";
    $functionCode .= "    use String_Handling\n";
    $functionCode .= "    implicit none\n";
    $functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
    # Iterate over pointers.
    $functionCode .= "    !\$omp critical(Node_XML_Dump)\n";
    $functionCode .= "    write (  idLabel,'(i20)') self         %index()\n";
    $functionCode .= "    write (treeLabel,'(i20)') self%hostTree%index\n";
    $functionCode .= "    write (fileHandle,'(a,a,a,a,a)') ' <node tree=\"',trim(adjustl(treeLabel)),'\" id=\"',trim(adjustl(idLabel)),'\" >'\n";
    $functionCode .= "    write (fileHandle,'(a)') '  <pointer>'\n";
    foreach my $pointer ( "parent", "firstChild", "sibling", "firstSatellite", "mergeTarget", "firstMergee", "siblingMergee", "formationNode" ) {
	$functionCode .= "   write (idLabel,'(i20)') self%".pad($pointer,14)."%index()\n";
	$functionCode .= "   write (fileHandle,'(a,a,a)') '   <".$pointer.">',trim(adjustl(idLabel)),'</".$pointer.">'\n"
    }
    $functionCode .= "    write (fileHandle,'(a)') '  </pointer>'\n";
    # Iterate over all component classes
    foreach ( @{$build->{'componentClassList'}} ) {	    
	$functionCode .= "    if (allocated(self%component".&Utils::padClass(ucfirst($_),[0,0]).")) then\n";
	$functionCode .= "      do i=1,size(self%component".&Utils::padClass(ucfirst($_),[0,0]).")\n";
	$functionCode .= "        call self%component".&Utils::padClass(ucfirst($_),[0,0])."(i)%dumpXML(fileHandle)\n";
	$functionCode .= "      end do\n";
	$functionCode .= "    end if\n";
    }
    $functionCode .= "    write (fileHandle,*) ' </node>'\n";
    $functionCode .= "    !\$omp end critical(Node_XML_Dump)\n";
    $functionCode .= "    return\n";
    $functionCode .= "  end subroutine Node_Dump_XML\n";
    # Insert into the function list.
    push(
	@{$build->{'code'}->{'functions'}},
	$functionCode
	);
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{type => "procedure", name => "dumpXML", function => "Node_Dump_XML", description => "Generate an XML dump of all content of a node.", returnType => "\\void", arguments => ""}
	);
    # Create a function for doing a raw (binary) dump.
    @dataContent =
	(
	 {
	     intrinsic  => "class",
	     type       => "treeNode",
	     attributes => [ "intent(in   )" ],
	     variables  => [ "self" ]
	 },
	 {
	     intrinsic  => "integer",
	     attributes => [ "intent(in   )" ],
	     variables  => [ "fileHandle" ]
	 },
	 {
	     intrinsic  => "integer",
	     variables  => [ "i" ]
	 }
	);
    $functionCode  = "  subroutine Node_Dump_Raw(self,fileHandle)\n";
    $functionCode .= "    !% Dump node content in binary.\n";
    $functionCode .= "    implicit none\n";
    $functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
    $functionCode .= "    write (fileHandle) self%isPhysicallyPlausible\n";
    # Iterate over all component classes
    foreach ( @{$build->{'componentClassList'}} ) {	    
	$functionCode .= "    write (fileHandle) allocated(self%component".&Utils::padClass(ucfirst($_),[0,0]).")\n";
	$functionCode .= "    if (allocated(self%component".&Utils::padClass(ucfirst($_),[0,0]).")) then\n";
	$functionCode .= "      select type (component => self%component".ucfirst($_)."(1))\n";
	$functionCode .= "      type is (nodeComponent".ucfirst($_).")\n";
	$functionCode .= "        write (fileHandle) .false.\n";
	$functionCode .= "      class is (nodeComponent".ucfirst($_).")\n";
	$functionCode .= "        write (fileHandle) .true.\n";
	$functionCode .= "      end select\n";
	$functionCode .= "      write (fileHandle) size(self%component".&Utils::padClass(ucfirst($_),[0,0]).")\n";
	$functionCode .= "      do i=1,size(self%component".&Utils::padClass(ucfirst($_),[0,0]).")\n";
	$functionCode .= "        call self%component".&Utils::padClass(ucfirst($_),[0,0])."(i)%dumpRaw(fileHandle)\n";
	$functionCode .= "      end do\n";
	$functionCode .= "    end if\n";
    }
    $functionCode .= "    return\n";
    $functionCode .= "  end subroutine Node_Dump_Raw\n";
    # Insert into the function list.
    push(
	@{$build->{'code'}->{'functions'}},
	$functionCode
	);
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{type => "procedure", name => "dumpRaw", function => "Node_Dump_Raw", description => "Generate a binary dump of all content of a node.", returnType => "\\void", arguments => "\\intzero\\ fileHandle\\argin"}
	);
    # Create a function for doing a raw (binary) read.
    @dataContent =
	(
	 {
	     intrinsic  => "class",
	     type       => "treeNode",
	     attributes => [ "intent(inout)", "target" ],
	     variables  => [ "self" ]
	 },
	 {
	     intrinsic  => "integer",
	     attributes => [ "intent(in   )" ],
	     variables  => [ "fileHandle" ]
	 },
	 {
	     intrinsic  => "integer",
	     variables  => [ "i", "componentCount" ]
	 },
	 {
	     intrinsic  => "logical",
	     variables  => [ "isAllocated" ]
	 }
	);
    $functionCode  = "  subroutine Node_Read_Raw(self,fileHandle)\n";
    $functionCode .= "    !% Dump node content in binary.\n";
    $functionCode .= "    implicit none\n";
    $functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
    $functionCode .= "    read (fileHandle) self%isPhysicallyPlausible\n";
    # Iterate over all component classes
    foreach ( @{$build->{'componentClassList'}} ) {	    
	$functionCode .= "    read (fileHandle) isAllocated\n";
	$functionCode .= "    if (isAllocated) then\n";
	$functionCode .= "      read (fileHandle) isAllocated\n";
	$functionCode .= "      read (fileHandle) componentCount\n";
	$functionCode .= "      if (allocated(self%component".ucfirst($_).")) deallocate(self%component".ucfirst($_).")\n";
	$functionCode .= "      if (isAllocated) then\n";
	$functionCode .= "        allocate(self%component".ucfirst($_)."(componentCount),source=default".ucfirst($_)."Component)\n";
	$functionCode .= "      else\n";
	$functionCode .= "        allocate(self%component".ucfirst($_)."(componentCount),source=".ucfirst($_)."Class)\n";
	$functionCode .= "      end if\n";
	$functionCode .= "      select type (self)\n";
	$functionCode .= "      type is (treeNode)\n";
	$functionCode .= "        do i=1,componentCount\n";
	$functionCode .= "          self%component".ucfirst($_)."(i)%hostNode => self\n";
	$functionCode .= "        end do\n";
	$functionCode .= "      end select\n";
	$functionCode .= "      do i=1,componentCount\n";
	$functionCode .= "        call self%component".&Utils::padClass(ucfirst($_),[0,0])."(i)%readRaw(fileHandle)\n";
	$functionCode .= "      end do\n";
	$functionCode .= "    else\n";
	$functionCode .= "       if (allocated(self%component".&Utils::padClass(ucfirst($_),[0,0]).")) deallocate(self%component".&Utils::padClass(ucfirst($_),[0,0]).")\n";
	$functionCode .= "       allocate(self%component".&Utils::padClass(ucfirst($_),[0,0])."(1))\n";
	$functionCode .= "    end if\n";
    }
    $functionCode .= "    return\n";
    $functionCode .= "  end subroutine Node_Read_Raw\n";
    # Insert into the function list.
    push(
	@{$build->{'code'}->{'functions'}},
	$functionCode
	);
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{type => "procedure", name => "readRaw", function => "Node_Read_Raw", description => "Read a binary dump of all content of a node.", returnType => "\\void", arguments => "\\intzero\\ fileHandle\\argin"}
	);
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
    $functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
    # Iterate over all component classes
    foreach ( @{$build->{'componentClassList'}} ) {	    
	$functionCode .= "    if (allocated(self%component".&Utils::padClass(ucfirst($_),[0,0]).")) then\n";
	$functionCode .= "      do i=1,size(self%component".&Utils::padClass(ucfirst($_),[0,0]).")\n";
	$functionCode .= "        call self%component".&Utils::padClass(ucfirst($_),[0,0])."(i)%outputCount(integerPropertyCount,doublePropertyCount,time,instance=i)\n";
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
    $functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
    # Iterate over all component classes
    foreach ( @{$build->{'componentClassList'}} ) {	    
	$functionCode .= "    if (allocated(self%component".&Utils::padClass(ucfirst($_),[0,0]).")) then\n";
	$functionCode .= "      do i=1,size(self%component".&Utils::padClass(ucfirst($_),[0,0]).")\n";
	$functionCode .= "        call self%component".&Utils::padClass(ucfirst($_),[0,0])."(i)%outputNames(integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI,doubleProperty,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time,instance=i)\n";
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
    $functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
    # Iterate over all component classes
    foreach ( @{$build->{'componentClassList'}} ) {	    
	$functionCode .= "    if (allocated(self%component".&Utils::padClass(ucfirst($_),[0,0]).")) then\n";
	$functionCode .= "      do i=1,size(self%component".&Utils::padClass(ucfirst($_),[0,0]).")\n";
	$functionCode .= "        call self%component".&Utils::padClass(ucfirst($_),[0,0])."(i)%output(integerProperty,integerBufferCount,integerBuffer,doubleProperty&
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

sub Generate_Node_Property_Name_From_Index_Function {
    # Generate function to get the name of a property given an index.
    my $build = shift;

    # Define variables.
    my @dataContent =
	(
	 {
	     intrinsic  => "class",
	     type       => "treeNode",
	     attributes => [ "intent(in   )" ],
	     variables  => [ "self" ]
	 },
	 {
	     intrinsic  => "integer",
	     attributes => [ "intent(in   )" ],
	     variables  => [ "index" ]
	 },
	 {
	     intrinsic  => "type",
	     type       => "varying_string",
	     variables  => [ "name" ]
	 },
	 {
	     intrinsic  => "integer",
	     variables  => [ "count", "i" ]
	 }
	);
    my $functionCode;
    $functionCode .= "  function Node_Property_Name_From_Index(self,index) result (name)\n";
    $functionCode .= "    !% Return the name of a property given its index.\n";
    $functionCode .= "    use ISO_Varying_String\n";
    $functionCode .= "    implicit none\n";
    $functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";

    # Loop over all component classes
    $functionCode .= "  name='unknown'\n";
    $functionCode .= "  count=index\n";
    foreach ( @{$build->{'componentClassList'}} ) {	    
     	$functionCode .= "    if (allocated(self%component".&Utils::padClass(ucfirst($_),[0,0]).")) then\n";
	$functionCode .= "      do i=1,size(self%component".&Utils::padClass(ucfirst($_),[0,0]).")\n";
	$functionCode .= "        call self%component".&Utils::padClass(ucfirst($_),[0,0])."(i)%nameFromIndex(count,name)\n";
	$functionCode .= "        if (count <= 0) return\n";
	$functionCode .= "      end do\n";
	$functionCode .= "    end if\n";
    }
    $functionCode .= "    if (name == 'unknown') call Galacticus_Error_Report('Node_Property_Name_From_Index','property index out of range')\n";
    $functionCode .= "    return\n";
    $functionCode .= "  end function Node_Property_Name_From_Index\n";
    # Insert into the function list.
    push(
	@{$build->{'code'}->{'functions'}},
	$functionCode
	);
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{type => "procedure", name => "nameFromIndex", function => "Node_Property_Name_From_Index", description => "Return the name of a property given its index in a node.", returnType => "\\textcolor{red}{\\textless varying\\_string\\textgreater}", arguments => "\\intzero\\ index\\argin"}
	);
}

sub Generate_Node_Serialization_Functions {
    # Generate functions to serialize/deserialize nodes to/from arrays.
    my $build = shift;

    # Function computing a count of the serialization length.
    my @dataContent =
	(
	 {
	     intrinsic  => "class",
	     type       => "treeNode",
	     attributes => [ "intent(in   )" ],
	     variables  => [ "self" ]
	 },
	 {
	     intrinsic  => "integer",
	     variables  => [ "count", "i" ]
	 }
	);
    my $functionCode;
    $functionCode .= "  function SerializeToArrayCount(self) result (count)\n";
    $functionCode .= "    !% Return a count of the size of the serialized {\\normalfont \\ttfamily treeNode} object.\n";
    $functionCode .= "    implicit none\n";
    $functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
    $functionCode .= "    count=0\n";
    # Loop over all component classes
    foreach ( @{$build->{'componentClassList'}} ) {	    
     	$functionCode .= "    if (allocated(self%component".&Utils::padClass(ucfirst($_),[0,0]).")) then\n";
	$functionCode .= "      do i=1,size(self%component".&Utils::padClass(ucfirst($_),[0,0]).")\n";
	$functionCode .= "        count=count+self%component".&Utils::padClass(ucfirst($_),[0,0])."(i)%serializeCount()\n";
	$functionCode .= "      end do\n";
	$functionCode .= "    end if\n";
    }
    $functionCode .= "    return\n";
    $functionCode .= "  end function SerializeToArrayCount\n\n";
    # Insert into the function list.
    push(
	@{$build->{'code'}->{'functions'}},
	$functionCode
	);
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{type => "procedure", name => "serializeCount", function => "serializeToArrayCount", description => "Return a count of the number of evolvable properties of the serialized object.", returnType => "\\intzero", arguments => ""}
	);
    # Create the serialization function.
    @dataContent =
	(
	 {
	     intrinsic  => "class",
	     type       => "treeNode",
	     attributes => [ "intent(in   )" ],
	     variables  => [ "self" ]
	 },
	 {
	     intrinsic  => "integer",
	     variables  => [ "count", "offset", "i" ],
	 },
	 {
	     intrinsic  => "double precision",
	     attributes => [ "dimension(:)", "intent(  out)" ],
	     variables  => [ "array" ]
	 }
	);
    $functionCode  = "  subroutine SerializeToArrayValues(self,array)\n";
    $functionCode .= "    !% Serialize values to array.\n";
    $functionCode .= "    use Memory_Management\n";
    $functionCode .= "    implicit none\n";
    $functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
    $functionCode .= "    offset=1\n";
    # Loop over all component classes
    foreach ( @{$build->{'componentClassList'}} ) {	    
	$functionCode .= "    if (allocated(self%component".&Utils::padClass(ucfirst($_),[0,0]).")) then\n";
	$functionCode .= "      do i=1,size(self%component".&Utils::padClass(ucfirst($_),[0,0]).")\n";
	$functionCode .= "        count=self%component".&Utils::padClass(ucfirst($_),[0,0])."(i)%serializeCount()\n";
	$functionCode .= "        if (count > 0) call self%component".&Utils::padClass(ucfirst($_),[0,0])."(i)%serializeValues(array(offset:))\n";
	$functionCode .= "        offset=offset+count\n";
	$functionCode .= "      end do\n";
	$functionCode .= "    end if\n";
    }
    $functionCode .= "    return\n";
    $functionCode .= "  end subroutine SerializeToArrayValues\n\n";
    # Insert into the function list.
    push(
	@{$build->{'code'}->{'functions'}},
	$functionCode
	);
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{type => "procedure", name => "serializeValues", function => "serializeToArrayValues", description => "Serialize values to {\\normalfont \\ttfamily array}.", returnType => "\\void", arguments => "\\doubleone\\ array\\argout"}
	);
    # Create the deserialization function.
    @dataContent =
	(
	 {
	     intrinsic  => "class",
	     type       => "treeNode",
	     attributes => [ "intent(inout)" ],
	     variables  => [ "self" ]
	 },
	 {
	     intrinsic  => "integer",
	     variables  => [ "count", "offset", "i" ],
	 },
	 {
	     intrinsic  => "double precision",
	     attributes => [ "dimension(:)", "intent(in   )" ],
	     variables  => [ "array" ]
	 }
	);
    $functionCode  = "  subroutine DeserializeFromArrayValues(self,array)\n";
    $functionCode .= "    !% Deserialize values from {\\normalfont \\ttfamily array}.\n";
    $functionCode .= "    use Memory_Management\n";
    $functionCode .= "    implicit none\n";
    $functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
    $functionCode .= "    offset=1\n";
    # Loop over all component classes
    foreach ( @{$build->{'componentClassList'}} ) {	    
	$functionCode .= "    if (allocated(self%component".&Utils::padClass(ucfirst($_),[0,0]).")) then\n";
	$functionCode .= "      do i=1,size(self%component".&Utils::padClass(ucfirst($_),[0,0]).")\n";
	$functionCode .= "        count=self%component".&Utils::padClass(ucfirst($_),[0,0])."(i)%serializeCount()\n";
	$functionCode .= "        if (count > 0) call self%component".&Utils::padClass(ucfirst($_),[0,0])."(i)%deserializeValues(array(offset:))\n";
	$functionCode .= "        offset=offset+count\n";
	$functionCode .= "      end do\n";
	$functionCode .= "    end if\n";
    }
    $functionCode .= "    return\n";
    $functionCode .= "  end subroutine DeserializeFromArrayValues\n\n";
    # Insert into the function list.
    push(
	@{$build->{'code'}->{'functions'}},
	$functionCode
	);
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{type => "procedure", name => "deserializeValues", function => "deserializeFromArrayValues", description => "Deserialize values from {\\normalfont \\ttfamily array}.", returnType => "\\void", arguments => "\\doubleone\\ array\\argin"}
	);
    # Generate serialization functions for scales and rates.
    foreach my $content ( "scale", "rate" ) {
	# Create the serialization function.
	@dataContent =
	    (
	     {
		 intrinsic  => "class",
		 type       => "treeNode",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "double precision",
		 attributes => [ "dimension(:)", "intent(  out)" ],
		 variables  => [ "array" ]
	     }
	    );
	$functionCode  = "  subroutine SerializeToArray".pad(ucfirst($content)."s",6)."(self,array)\n";
	$functionCode .= "    !% Serialize ".$content."s to array.\n";
	$functionCode .= "    implicit none\n";
	$functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
	$functionCode .= "    !GCC\$ attributes unused :: self\n";
	$functionCode .= "    array(1:nodeSerializationCount)=node".ucfirst($content)."s(1:nodeSerializationCount)\n";
	$functionCode .= "    return\n";
	$functionCode .= "  end subroutine SerializeToArray".ucfirst($content)."s\n\n";
	# Insert into the function list.
	push(
	    @{$build->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Insert a type-binding for this function into the treeNode type.
	push(
	    @{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	    {type => "procedure", name => "serialize".ucfirst($content)."s", function => "serializeToArray".ucfirst($content)."s", description => "Serialize ".$content."s to {\\normalfont \\ttfamily array}.", returnType => "\\void", arguments => "\\doubleone\\ array\\argout"}
	    );
    }
}

sub Generate_Node_ODE_Initialization_Functions {
    # Generate functions initialize a node for an ODE step.
    my $build = shift;
    # Create functions to initialize property rates for an ODE step.
    my @dataContent =
	(
	 {
	     intrinsic  => "class",
	     type       => "treeNode",
	     attributes => [ "intent(inout)" ],
	     variables  => [ "self" ]
	 }
	);
    my $functionCode;
    $functionCode  = "  subroutine Tree_Node_ODE_Step_Rates_Initialize(self)\n";
    $functionCode .= "    !% Initialize the rates in components of tree node {\\normalfont \\ttfamily self} in preparation for an ODE solver step.\n";
    $functionCode .= "    implicit none\n";
    $functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
    $functionCode .= "    !GCC\$ attributes unused :: self\n";
    $functionCode .= "    nodeRates=0.0d0\n";
    $functionCode .= "    return\n";
    $functionCode .= "  end subroutine Tree_Node_ODE_Step_Rates_Initialize\n\n";
    # Insert into the function list.
    push(
	@{$build->{'code'}->{'functions'}},
	$functionCode
	);
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{type => "procedure", name => "odeStepRatesInitialize", function => "Tree_Node_ODE_Step_Rates_Initialize", description => "Initialize rates of evolvable properties.", returnType => "\\void", arguments => ""},
	);    
    # Create functions to initialize property scales for an ODE step.
    $functionCode  = "  subroutine Tree_Node_ODE_Step_Scales_Initialize(self)\n";
    $functionCode .= "    !% Initialize the scales in components of tree node {\\normalfont \\ttfamily self} in preparation for an ODE solver step.\n";
    $functionCode .= "    implicit none\n";
    $functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
    $functionCode .= "    !GCC\$ attributes unused :: self\n";
    $functionCode .= "    nodeScales=1.0d0\n";
    $functionCode .= "    return\n";
    $functionCode .= "  end subroutine Tree_Node_ODE_Step_Scales_Initialize\n\n";
    # Insert into the function list.
    push(
	@{$build->{'code'}->{'functions'}},
	$functionCode
	);
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{type => "procedure", name => "odeStepScalesInitialize" , function => "Tree_Node_ODE_Step_Scales_Initialize", description => "Initialize tolerance scales of evolvable properties.", returnType => "\\void", arguments => ""}
	);
}

sub Generate_Implementation_Dump_Functions {
    # Generate dump for each component implementation.
    my $build = shift;
    # Iterate over component implementations.
    foreach my $componentID ( @{$build->{'componentIdList'}} ) {
	# Get the component.
	my $component = $build->{'components'}->{$componentID};
	# Initialize function code.
	my $functionCode;
	# Initialize data content.
	my @dataContent =
	    (
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($componentID),
		 attributes => [ "intent(in   )" ],
		 variables  => [ "self" ]
	     }
	    );
	my $counterAdded = 0;
	unless ( $component->{'name'} eq "null" ) {
	    push(
		@dataContent,
		{
		    intrinsic  => "type",
		    type       => "varying_string",
		    variables  => [ "message" ]
		},
		{
		    intrinsic  => "character",
		    type       => "len=18",
		    variables  => [ "label" ]
		}
		);
	    foreach my $propertyName ( &ExtraUtils::sortedKeys($component->{'properties'}->{'property'}) ) {
		my $property = $component->{'properties'}->{'property'}->{$propertyName};
		# Check if this property has any linked data in this component.
		if ( exists($property->{'linkedData'}) ) {
		    my $linkedDataName = $property->{'linkedData'};
		    my $linkedData     = $component->{'content'}->{'data'}->{$linkedDataName};
		    if ( $linkedData->{'rank'} == 1 && $counterAdded == 0 ) {
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
	# Format labels for different data types.
	my %formatLabel = 
	    (
	     "double" => "'(e12.6)'",
	     "integer"         => "'(i8)'"   ,
	     "longInteger"     => "'(i16)'"  ,
	     "logical"         => "'(l1)'"
	    );
	# Generate dump function.
	$functionCode  = "  subroutine Node_Component_".ucfirst($componentID)."_Dump(self)\n";
	$functionCode .= "    !% Dump the contents of a ".$component->{'name'}." implementation of the ".$component->{'class'}." component.\n";
	$functionCode .= "    use Galacticus_Display\n";
	$functionCode .= "    use ISO_Varying_String\n";
	$functionCode .= "    use String_Handling\n";
	$functionCode .= "    implicit none\n";
	$functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
	$functionCode .= "    !GCC\$ attributes unused :: self\n"
	    if ( $component->{'name'} eq "null" );
	unless ( $component->{'name'} eq "null" ) {
	    # Dump the parent type if necessary.
	    $functionCode .= "    call self%nodeComponent".ucfirst($component->{'extends'}->{'class'}).ucfirst($component->{'extends'}->{'name'})."%dump()\n"
		if ( exists($component->{'extends'}) );
	    $functionCode .= "    call Galacticus_Display_Indent('".$component->{'class'}.": ".(" " x ($Utils::fullyQualifiedNameLengthMax-length($component->{'class'}))).$component->{'name'}."')\n";
	    foreach my $propertyName ( &ExtraUtils::sortedKeys($component->{'properties'}->{'property'}) ) {
		my $property = $component->{'properties'}->{'property'}->{$propertyName};
		# Check if this property has any linked data in this component.
		if ( exists($property->{'linkedData'}) ) {
		    my $linkedDataName = $property->{'linkedData'};
		    my $linkedData     = $component->{'content'}->{'data'}->{$linkedDataName};
		    if ( $linkedData->{'rank'} == 0 ) {
			if (&Utils::isIntrinsic($linkedData->{'type'})) {
			    $functionCode .= "    write (label,".$formatLabel{$linkedData->{'type'}}.") self%".&Utils::padLinkedData($linkedDataName,[0,0])."\n";
			    $functionCode .= "    message='".$propertyName.": ".(" " x ($Utils::implementationPropertyNameLengthMax-length($propertyName)))."'//label\n";
			    $functionCode .= "    call Galacticus_Display_Message(message)\n";
			}
			else {
			    $functionCode .= "    message='".$propertyName.":'\n";
			    $functionCode .= "    call Galacticus_Display_Indent(message)\n";
			    $functionCode .= "    call self%".&Utils::padLinkedData($linkedDataName,[0,0])."%dump()\n";
			    $functionCode .= "    call Galacticus_Display_Unindent('end')\n";
			}
		    } elsif ( $linkedData->{'rank'} == 1 ) {
			if (&Utils::isIntrinsic($linkedData->{'type'})) {
			    $functionCode .= "    do i=1,size(self%".$linkedDataName.")\n";
			    $functionCode .= "       write (label,'(i3)') i\n";
			    $functionCode .= "       message='".$propertyName.": ".(" " x ($Utils::implementationPropertyNameLengthMax-length($propertyName)))." '//trim(label)\n";
			    $functionCode .= "       write (label,".$formatLabel{$linkedData->{'type'}}.") self%".$linkedDataName."(i)\n";
			    $functionCode .= "       message=message//': '//label\n";
			    $functionCode .= "       call Galacticus_Display_Message(message)\n";
			    $functionCode .= "    end do\n";
			}
			else {
			    $functionCode .= "    do i=1,size(self%".$linkedDataName.")\n";
			    $functionCode .= "       write (label,'(i3)') i\n";
			    $functionCode .= "       message='".$propertyName.": ".(" " x ($Utils::implementationPropertyNameLengthMax-length($propertyName)))." '//trim(label)\n";
			    $functionCode .= "       call Galacticus_Display_Indent(message)\n";
			    $functionCode .= "       call self%".$linkedDataName."(i)%dump()\n";
			    $functionCode .= "       call Galacticus_Display_Unindent('end')\n";
			    $functionCode .= "    end do\n";
			}
		    }
		}
	    }
	    $functionCode .= "    call Galacticus_Display_Unindent('done')\n";
	}
	$functionCode .= "    return\n";
	$functionCode .= "  end subroutine Node_Component_".ucfirst($componentID)."_Dump\n";
	# Insert into the function list.
	push(
	    @{$build->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Insert a type-binding for this function into the implementation type.
	push(
	    @{$build->{'types'}->{'nodeComponent'.ucfirst($componentID)}->{'boundFunctions'}},
	    {type => "procedure", name => "dump", function => "Node_Component_".ucfirst($componentID)."_Dump"},
	    );
	# Initialize data content.
	@dataContent =
	    (
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($componentID),
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "fileHandle" ]
	     }
	    );
	push(
	    @dataContent,
	    {
		intrinsic  => "integer",
		variables  => [ "i" ]
	    }
	    )
	    if ( $counterAdded == 1 );
	# Generate XML dump function.
	$functionCode  = "  subroutine Node_Component_".ucfirst($componentID)."_Dump_XML(self,fileHandle)\n";
	$functionCode .= "    !% Dump the contents of a ".$component->{'name'}." implementation of the ".$component->{'class'}." component to XML.\n";
	$functionCode .= "    implicit none\n";
	$functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
	my $selfUsed = 0;
	my $fileUsed = 0;
	my $functionBody = "";
	unless ( $component->{'name'} eq "null" ) {
	    $fileUsed = 1;
	    # Dump the parent type if necessary.
	    if ( exists($component->{'extends'}) ) {
		$functionBody .= "    call self%nodeComponent".ucfirst($component->{'extends'}->{'class'}).ucfirst($component->{'extends'}->{'name'})."%dumpXML(fileHandle)\n";
		$selfUsed = 1;
	    }
	    $functionBody .= "    write (fileHandle,'(a)') '  <".$component->{'class'}." type=\"".$component->{'name'}."\">'\n";
            foreach my $propertyName ( &ExtraUtils::sortedKeys($component->{'properties'}->{'property'}) ) {
		my $property = $component->{'properties'}->{'property'}->{$propertyName};
		# Check if this property has any linked data in this component.
		if ( exists($property->{'linkedData'}) ) {
		    $selfUsed = 1;
		    my $linkedDataName = $property->{'linkedData'};
		    my $linkedData     = $component->{'content'}->{'data'}->{$linkedDataName};
		    if ( $linkedData->{'rank'} == 0 ) {
			if (&Utils::isIntrinsic($linkedData->{'type'})) {
			    (my $typeFormat = $formatLabel{$linkedData->{'type'}}) =~ s/^\'\((.*)\)\'$/$1/g;
				$functionBody .= "    write (fileHandle,'(a,".$typeFormat.",a)') '   <".$propertyName.">',self%".&Utils::padLinkedData($linkedDataName,[0,0]).",'</".$propertyName.">'\n";
			}
			else {
			    $functionBody .= "    write (fileHandle,'(a)') '   <".$propertyName.">'\n";
			    $functionBody .= "    write (fileHandle,'(a)') '   </".$propertyName.">'\n";
			}
		    } elsif ( $linkedData->{'rank'} == 1 ) {
			if (&Utils::isIntrinsic($linkedData->{'type'})) {
			    (my $typeFormat = $formatLabel{$linkedData->{'type'}}) =~ s/^\'\((.*)\)\'$/$1/g;
			    $functionBody .= "    do i=1,size(self%".$linkedDataName.")\n";
			    $functionBody .= "       write (fileHandle,'(a,".$typeFormat.",a)') '   <".$propertyName.">',self%".$linkedDataName."(i),'</".$propertyName.">'\n";
			    $functionBody .= "    end do\n";
			}
			else {
			    $functionBody .= "    do i=1,size(self%".$linkedDataName.")\n";
			    $functionBody .= "       write (fileHandle,'(a)') '   <".$propertyName.">'\n";
			    $functionBody .= "       write (fileHandle,'(a)') '   </".$propertyName.">'\n";
			    $functionBody .= "    end do\n";
			}			
		    }
		} elsif ( $property->{'isVirtual'} && $property->{'rank'} == 0 ) {
		    if (&Utils::isIntrinsic($property->{'type'})) {
			(my $typeFormat = $formatLabel{$property->{'type'}}) =~ s/^\'\((.*)\)\'$/$1/g;
			$functionBody .= "    write (fileHandle,'(a,".$typeFormat.",a)') '   <".$propertyName.">',self%".&Utils::padImplementationPropertyName($propertyName,[0,0])."(),'</".$propertyName.">'\n";
		    }
		}
	    }
	    $functionBody .= "    write (fileHandle,'(a)') '  </".$component->{'class'}.">'\n";
	}
	$functionBody .= "    return\n";
	$functionBody .= "  end subroutine Node_Component_".ucfirst($componentID)."_Dump_XML\n";
	$functionCode .= "  !GCC\$ attributes unused :: self\n"
	    unless ( $selfUsed );
	$functionCode .= "  !GCC\$ attributes unused :: fileHandle\n"
	    unless ( $fileUsed );
	$functionCode .= $functionBody;
	# Insert into the function list.
	push(
	    @{$build->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Insert a type-binding for this function into the implementation type.
	push(
	    @{$build->{'types'}->{'nodeComponent'.ucfirst($componentID)}->{'boundFunctions'}},
	    {type => "procedure", name => "dumpXML", function => "Node_Component_".ucfirst($componentID)."_Dump_XML"},
	    );
	# Create function to do a raw (binary) dump.
	# Initialize data content.
	@dataContent =
	    (
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($componentID),
		 attributes => [ "intent(in   )" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "fileHandle" ]
	     }
	    );
	$counterAdded = 0;
	foreach my $propertyName ( &ExtraUtils::sortedKeys($component->{'properties'}->{'property'}) ) {
	    my $property = $component->{'properties'}->{'property'}->{$propertyName};
	    # Check if this property has any linked data in this component.
	    if ( exists($property->{'linkedData'}) ) {
		my $linkedDataName = $property->{'linkedData'};
		my $linkedData     = $component->{'content'}->{'data'}->{$linkedDataName};
		if ( $linkedData->{'rank'} == 1 && $counterAdded == 0) {
		    unless (&Utils::isIntrinsic($linkedData->{'type'})) {
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
	# Generate raw dump function.
	$functionCode  = "  subroutine Node_Component_".ucfirst($componentID)."_Dump_Raw(self,fileHandle)\n";
	$functionCode .= "    !% Dump the contents of a ".$component->{'name'}." implementation of the ".$component->{'class'}." component in binary.\n";
	$functionCode .= "    implicit none\n";
	$functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
	$functionBody = "";
	$selfUsed     = 0;
	$fileUsed     = 0;	
	unless ( $component->{'name'} eq "null" ) {
	    # Dump the parent type if necessary.
	    if ( exists($component->{'extends'}) ) {
		$functionBody .= "    call self%nodeComponent".ucfirst($component->{'extends'}->{'class'}).ucfirst($component->{'extends'}->{'name'})."%dumpRaw(fileHandle)\n";
		$selfUsed = 1;
		$fileUsed = 1;
	    }
	    foreach my $propertyName ( &ExtraUtils::sortedKeys($component->{'properties'}->{'property'}) ) {
		my $property = $component->{'properties'}->{'property'}->{$propertyName};
		# Check if this property has any linked data in this component.
		if ( exists($property->{'linkedData'}) ) {
		    $selfUsed = 1;
		    $fileUsed = 1;
		    my $linkedDataName = $property->{'linkedData'};
		    my $linkedData     = $component->{'content'}->{'data'}->{$linkedDataName};
		    if ( $linkedData->{'rank'} == 0 ) {
			if (&Utils::isIntrinsic($linkedData->{'type'})) {
			    $functionCode .= "    write (fileHandle) self%".&Utils::padLinkedData($linkedDataName,[0,0])."\n";
			    $functionBody .= "    write (fileHandle) self%".&Utils::padLinkedData($linkedDataName,[0,0])."\n";
			}
			else {
			    $functionBody .= "    call self%".&Utils::padLinkedData($linkedDataName,[0,0])."%dumpRaw(fileHandle)\n";
			}
		    } elsif ( $linkedData->{'rank'} == 1 ) {
			$functionBody .= "    write (fileHandle) allocated(self%".$linkedDataName.")\n";
			$functionBody .= "    if (allocated(self%".$linkedDataName.")) then\n";
			$functionBody .= "       write (fileHandle) size(self%".$linkedDataName.")\n";
			if (&Utils::isIntrinsic($linkedData->{'type'})) {
			    $functionBody .= "      write (fileHandle) self%".$linkedDataName."\n";
			}
			else {
			    $functionBody .= "       do i=1,size(self%".$linkedDataName.")\n";
			    $functionBody .= "          call self%".$linkedDataName."(i)%dumpRaw(fileHandle)\n";
			    $functionBody .= "       end do\n";
			}
			$functionBody .= "    end if\n";
		    }
		}
	    }
	}
	$functionBody .= "    return\n";
	$functionBody .= "  end subroutine Node_Component_".ucfirst($componentID)."_Dump_Raw\n";
	$functionCode .= "  !GCC\$ attributes unused :: self\n"
	    unless ( $selfUsed );
	$functionCode .= "  !GCC\$ attributes unused :: fileHandle\n"
	    unless ( $fileUsed );
	$functionCode .= $functionBody;
	# Insert into the function list.
	push(
	    @{$build->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Insert a type-binding for this function into the implementation type.
	push(
	    @{$build->{'types'}->{'nodeComponent'.ucfirst($componentID)}->{'boundFunctions'}},
	    {type => "procedure", name => "dumpRaw", function => "Node_Component_".ucfirst($componentID)."_Dump_Raw"},
	    );
	# Create function to do a raw (binary) read.
	# Initialize data content.
	@dataContent =
	    (
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($componentID),
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "fileHandle" ]
	     }
	    );
	my $readCounterAdded = 0;
	my $readArraysAdded  = 0;
	foreach my $propertyName ( &ExtraUtils::sortedKeys($component->{'properties'}->{'property'}) ) {
	    my $property = $component->{'properties'}->{'property'}->{$propertyName};
	    # Check if this property has any linked data in this component.
	    if ( exists($property->{'linkedData'}) ) {
		my $linkedDataName = $property->{'linkedData'};
		my $linkedData     = $component->{'content'}->{'data'}->{$linkedDataName};
		if ( $linkedData->{'rank'} == 1 ) {
		    if ( $readArraysAdded == 0 ) {
			push(
			    @dataContent,
			    {
				intrinsic  => "integer",
				variables  => [ "arraySize" ]
			    },
			    {
				intrinsic  => "logical",
				variables  => [ "isAllocated" ]
			    }
			    );
			$readArraysAdded = 1;
		    }
		    if ( $readCounterAdded == 0 ) {
			unless (&Utils::isIntrinsic($linkedData->{'type'})) {
			    push(
				@dataContent,
				{
				    intrinsic  => "integer",
				    variables  => [ "i" ]
				}
				);
			    $readCounterAdded = 1;
			}
		    }
		}
	    }
	}
	# Generate raw read function.
	$functionCode  = "  subroutine Node_Component_".ucfirst($componentID)."_Read_Raw(self,fileHandle)\n";
	$functionCode .= "    !% Read the contents of a ".$component->{'name'}." implementation of the ".$component->{'class'}." component in binary.\n";
	$functionCode .= "    use Memory_Management\n";
	$functionCode .= "    implicit none\n";
	$functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
	$selfUsed     = 0;
	$fileUsed     = 0;
	$functionBody = "";
	unless ( $component->{'name'} eq "null" ) {
	    # Dump the parent type if necessary.
	    if ( exists($component->{'extends'}) ) {
		$functionBody .= "    call self%nodeComponent".ucfirst($component->{'extends'}->{'class'}).ucfirst($component->{'extends'}->{'name'})."%readRaw(fileHandle)\n";
		$selfUsed = 1;
		$fileUsed = 1;
	    }
            foreach my $propertyName ( &ExtraUtils::sortedKeys($component->{'properties'}->{'property'}) ) {
		my $property = $component->{'properties'}->{'property'}->{$propertyName};
		# Check if this property has any linked data in this component.
		if ( exists($property->{'linkedData'}) ) {
		    $selfUsed = 1;
		    $fileUsed = 1;
		    my $linkedDataName = $property->{'linkedData'};
		    my $linkedData     = $component->{'content'}->{'data'}->{$linkedDataName};
		    if ( $linkedData->{'rank'} == 0 ) {
			if (&Utils::isIntrinsic($linkedData->{'type'})) {
			    $functionBody .= "    read (fileHandle) self%".&Utils::padLinkedData($linkedDataName,[0,0])."\n";
			} else {
			    $functionBody .= "    call self%".&Utils::padLinkedData($linkedDataName,[0,0])."%readRaw(fileHandle)\n";
			}
		    } elsif ( $linkedData->{'rank'} == 1 ) {
			$functionBody .= "    read (fileHandle) isAllocated\n";
			$functionBody .= "    if (isAllocated) then\n";
			$functionBody .= "       read (fileHandle) arraySize\n";
			if (&Utils::isIntrinsic($linkedData->{'type'})) {
			    $functionBody .= "      call Alloc_Array(self%".$linkedDataName.",[arraySize])\n";
			    $functionBody .= "      read (fileHandle) self%".$linkedDataName."\n";
			}
			else {
			    $functionBody .= "       allocate(self%".$linkedDataName."(arraySize))\n";
			    $functionBody .= "       do i=1,arraySize)\n";
			    $functionBody .= "          call self%".$linkedDataName."(i)%readRaw(fileHandle)\n";
			    $functionBody .= "       end do\n";
			}
			$functionBody .= "    end if\n";
		    }
		}
	    }
	}
	$functionBody .= "    return\n";
	$functionBody .= "  end subroutine Node_Component_".ucfirst($componentID)."_Read_Raw\n";
	$functionCode .= "  !GCC\$ attributes unused :: self\n"
	    unless ( $selfUsed );
	$functionCode .= "  !GCC\$ attributes unused :: fileHandle\n"
	    unless ( $fileUsed );
	$functionCode .= $functionBody;
	# Insert into the function list.
	push(
	    @{$build->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Insert a type-binding for this function into the implementation type.
	push(
	    @{$build->{'types'}->{'nodeComponent'.ucfirst($componentID)}->{'boundFunctions'}},
	    {type => "procedure", name => "readRaw", function => "Node_Component_".ucfirst($componentID)."_Read_Raw", description => "Read a binary dump of the {\\normalfont \\ttfamily nodeComponent} from the given {\\normalfont \\ttfamily fileHandle}.", returnType => "\\void", arguments => "\\intzero\\ fileHandle\\argin"},
	    );
    }
}

sub Generate_Implementation_Initializor_Functions {
    # Generate initializor for each component implementation.
    my $build = shift;
    # Iterate over component implementations.
    foreach my $componentID ( @{$build->{'componentIdList'}} ) {
	# Get the component.
	my $component = $build->{'components'}->{$componentID};
	# Initialize function code.
	my $functionCode;
	# Initialize data content.
	my @dataContent =
	    (
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($componentID),
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     }
	    );
	# Generate the initialization code.
	my %requiredComponents;
	my $initializeCode = "";
	foreach my $propertyName ( &ExtraUtils::sortedKeys($component->{'properties'}->{'property'}) ) {
	    my $property = $component->{'properties'}->{'property'}->{$propertyName};
	    if ( exists($property->{'linkedData'}) ) {
		my $linkedDataName = $property->{'linkedData'};
		my $linkedData     = $component->{'content'}->{'data'}->{$linkedDataName};
		# Set to a class default value if available.
		if ( exists($property->{'classDefault'}) ) {
		    my $default = $property->{'classDefault'}->{'code'};
		    while ( $default =~ m/self([a-zA-Z]+)Component\s*%/ ) {
			$requiredComponents{$1} = 1;
			$default =~ s/self([a-zA-Z]+)Component\s*%//;
		    }
		    $default = $property->{'classDefault'}->{'code'};
		    if ( exists($property->{'classDefault'}->{'count'}) ) {
			$initializeCode .= "           call Alloc_Array(self%".&Utils::padLinkedData($linkedDataName,[0,0]).",[".$property->{'classDefault'}->{'count'}."])\n";
		    }
		    $initializeCode .= "            self%".&Utils::padLinkedData($linkedDataName,[0,0])."=".$default."\n";
		} else {
		    # Set to null.
		    if    ( $linkedData->{'type'} eq"double" ) {
			if ( $linkedData->{'rank'} == 0 ) {
			    $initializeCode .= "            self%".&Utils::padLinkedData($linkedDataName,[0,0])."=0.0d0\n";
			} else {
			    $initializeCode .= "            call Alloc_Array(self%".&Utils::padLinkedData($linkedDataName,[0,0]).",[".join(",","0" x $linkedData->{'rank'})."])\n";
			}
		    }
		    elsif ( $linkedData->{'type'} eq"integer"     ) {
			$initializeCode .= "            self%".&Utils::padLinkedData($linkedDataName,[0,0])."=0\n";
		    }
		    elsif ( $linkedData->{'type'} eq"longInteger" ) {
			$initializeCode .= "            self%".&Utils::padLinkedData($linkedDataName,[0,0])."=0_kind_int8\n";
		    }
		    elsif ( $linkedData->{'type'} eq"logical"     ) {
			$initializeCode .= "            self%".&Utils::padLinkedData($linkedDataName,[0,0])."=.false.\n";
		    }
		    else {
			$initializeCode .= "       call self%".&Utils::padLinkedData($linkedDataName,[0,0])."%reset()\n";			    
		    }
		}
	    }
	}
	# Add pointers for each required component.
	push(
	    @dataContent,
	    {
		intrinsic  => "class",
		type       => "nodeComponent".ucfirst($_),
		attributes => [ "pointer" ],
		variables  => [ "self".ucfirst($_)."Component" ]
	    }
	    )
	    foreach ( &ExtraUtils::sortedKeys(\%requiredComponents) );
	# Generate initializor function.
	my $selfUsed     = 0;
	my $functionBody = "";
	$functionCode  = "  subroutine Node_Component_".ucfirst($componentID)."_Initializor(self)\n";
	$functionCode .= "    !% Initialize a ".$component->{'name'}." implementation of the ".$component->{'class'}." component.\n";
	$functionCode .= "    use Memory_Management\n";
	# Insert any required modules.
	my %requiredModules;
	foreach my $propertyName ( &ExtraUtils::sortedKeys($component->{'properties'}->{'property'}) ) {
	    my $property = $component->{'properties'}->{'property'}->{$propertyName};
	    if ( exists($property->{'classDefault'}) && exists($property->{'classDefault'}->{'modules'}) ) {
		foreach ( @{$property->{'classDefault'}->{'modules'}} ) {
		    $requiredModules{$_} = 1;
		}
	    }
	}
	foreach ( &ExtraUtils::sortedKeys(\%requiredModules) ) {
	    $functionCode .= "    use ".$_."\n";
	}
	$functionCode .= "    implicit none\n";
	$functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";	
	unless ( $component->{'name'} eq "null" ) {
	    # Initialize the parent type if necessary.
	    if ( exists($component->{'extends'}) ) {
		$functionBody .= "    call self%nodeComponent".ucfirst($component->{'extends'}->{'class'}).ucfirst($component->{'extends'}->{'name'})."%initialize()\n";
		$selfUsed = 1;
	    }
	}
	foreach my $requiredComponent ( &ExtraUtils::sortedKeys(\%requiredComponents) ) {
	    $functionBody .= "     self".$requiredComponent."Component => self%hostNode%".lc($requiredComponent)."()\n";
	    $selfUsed = 1;
	}
	$functionBody .= $initializeCode;
	$functionBody .= "    return\n";
	$functionBody .= "  end subroutine Node_Component_".ucfirst($componentID)."_Initializor\n";
	$functionCode .= "  !GCC\$ attributes unused :: self\n"
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
	    {type => "procedure", name => "initialize", function => "Node_Component_".ucfirst($componentID)."_Initializor"},
	    );
    }
}

sub Generate_Implementation_Builder_Functions {
    # Generate builder for each component implementation.
    my $build = shift;
    # Iterate over component implementations.
    foreach my $componentID ( @{$build->{'componentIdList'}} ) {
	# Get the component.
	my $component = $build->{'components'}->{$componentID};
	# Initialize function code.
	my $functionCode;
	# Initialize data content.
	my @dataContent =
	    (
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($componentID),
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "type",
		 type       => "node",
		 attributes => [ "intent(in   )", "pointer" ],
		 variables  => [ "componentDefinition" ]
	     }
	    );
	unless ( $component->{'name'} eq "null" ) {
	    push(
		@dataContent,
		{
		    intrinsic  => "type",
		    type       => "node",
		    attributes => [ "pointer" ],
		    variables  => [ "property" ]
		},
		{
		    intrinsic  => "type",
		    type       => "nodeList",
		    attributes => [ "pointer" ],
		    variables  => [ "propertyList" ]
		}
		);
	    my $counterAdded = 0;
	    foreach my $propertyName ( &ExtraUtils::sortedKeys($component->{'properties'}->{'property'}) ) {
		my $property = $component->{'properties'}->{'property'}->{$propertyName};
		# Check if this property has any linked data in this component.
		if ( exists($property->{'linkedData'}) ) {
		    my $linkedDataName = $property->{'linkedData'};
		    my $linkedData     = $component->{'content'}->{'data'}->{$linkedDataName};
		    if ( $linkedData->{'rank'} == 1 && $counterAdded == 0 ) {
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
	# Generate builder function.
	$functionCode  = "  subroutine Node_Component_".ucfirst($componentID)."_Builder(self,componentDefinition)\n";
	$functionCode .= "    !% Build a ".$component->{'name'}." implementation of the ".$component->{'class'}." component.\n";
	$functionCode .= "    use FoX_DOM\n";
	$functionCode .= "    use Galacticus_Error\n";
	$functionCode .= "    use Memory_Management\n";
	$functionCode .= "    implicit none\n";
	$functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
	$functionCode .= "    !GCC\$ attributes unused :: self, componentDefinition\n"
	    if ( $component->{'name'} eq "null" );
	unless ( $component->{'name'} eq "null" ) {
	    # Initialize the component.
	    $functionCode .= "    call self%initialize()\n";
	    # Build the parent type if necessary.
	    $functionCode .= "    call self%nodeComponent".ucfirst($component->{'extends'}->{'class'}).ucfirst($component->{'extends'}->{'name'})."%builder(componentDefinition)\n"
		if ( exists($component->{'extends'}) );
	    # Enter critical section.
	    $functionCode .= "    !\$omp critical (FoX_DOM_Access)\n";
	    foreach my $propertyName ( sort(keys(%{$component->{'properties'}->{'property'}})) ) {
		my $property = $component->{'properties'}->{'property'}->{$propertyName};
		# Check if this property has any linked data in this component.
		if ( exists($property->{'linkedData'}) ) {
		    my $linkedDataName = $property->{'linkedData'};
		    my $linkedData     = $component->{'content'}->{'data'}->{$linkedDataName};
		    $functionCode .= "    propertyList => getElementsByTagName(componentDefinition,'".$propertyName."')\n";
		    if ( $linkedData->{'rank'} == 0 ) {
			$functionCode .= "    if (getLength(propertyList) > 1) call Galacticus_Error_Report('Node_Component_".ucfirst($componentID)."_Builder','scalar property must have precisely one value')\n";
			$functionCode .= "    if (getLength(propertyList) == 1) then\n";
			$functionCode .= "      property => item(propertyList,0)\n";
			if (
			    $linkedData->{'type'} eq "double"
			    ||
			    $linkedData->{'type'} eq "integer"
			    ||
			    $linkedData->{'type'} eq "logical"
			    ) {
			    $functionCode .= "      call extractDataContent(property,self%".&Utils::padLinkedData($linkedDataName,[0,0]).")\n";
			} elsif ( $linkedData->{'type'} eq "longInteger" ) {
			    $functionCode .= "      call Galacticus_Error_Report('Node_Component_".ucfirst($componentID)."_Builder','building of long integer properties currently not supported')\n";
			}
			else {
			    $functionCode .= "      call self%".&Utils::padLinkedData($linkedDataName,[0,0])."%builder(property)\n";
			}
			$functionCode .= "    end if\n";
		    } elsif ( $linkedData->{'rank'} == 1 ) {
			$functionCode .= "    if (getLength(propertyList) >= 1) then\n";
			if (&Utils::isIntrinsic($linkedData->{'type'})) {
			    $functionCode .= "      call Alloc_Array(self%".$linkedDataName.",[getLength(propertyList)])\n";
			    $functionCode .= "      do i=1,getLength(propertyList)\n";
			    $functionCode .= "        property => item(propertyList,i-1)\n";
			    $functionCode .= "        call extractDataContent(property,self%".$linkedDataName."(i))\n";
			    $functionCode .= "      end do\n";
			}
			else {
			    $functionCode .= "      allocate(self%".$linkedDataName."(getLength(propertyList)))\n";
			    $functionCode .= "      do i=1,getLength(propertyList)\n";
			    $functionCode .= "        property => item(propertyList,i-1)\n";
			    $functionCode .= "        call self%".$linkedDataName."(i)%builder(property)\n";
			    $functionCode .= "      end do\n";
			}
			$functionCode .= "    end if\n";
		    }
		}
	    }
	    $functionCode .= "    !\$omp end critical (FoX_DOM_Access)\n";
	}
	$functionCode .= "    return\n";
	$functionCode .= "  end subroutine Node_Component_".ucfirst($componentID)."_Builder\n";
	# Insert into the function list.
	push(
	    @{$build->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Insert a type-binding for this function into the implementation type.
	push(
	    @{$build->{'types'}->{'nodeComponent'.ucfirst($componentID)}->{'boundFunctions'}},
	    {type => "procedure", name => "builder", function => "Node_Component_".ucfirst($componentID)."_Builder"},
	    );
    }
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
	foreach my $propertyName ( &ExtraUtils::sortedKeys($component->{'properties'}->{'property'}) ) {
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
	foreach my $propertyName ( &ExtraUtils::sortedKeys($component->{'properties'}->{'property'}) ) {
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
		} elsif ( $property->{'isVirtual'} && $property->{'attributes'}->{'isGettable'} ) {
		    $rank = $property->{'rank'};
		    $type = $property->{'type'};
		} else {
		    die("Generate_Implementation_Output_Functions(): can not output [".$propertyName."]");
		}
		# Increment the counters.
		if (&Utils::isOutputIntrinsic($type)) {
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
	foreach my $propertyName ( &ExtraUtils::sortedKeys($component->{'properties'}->{'property'}) ) {
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
		    unless (&Utils::isOutputIntrinsic($type));
	    }
	}
	my @outputTypes;
	foreach ( &ExtraUtils::sortedKeys(\%outputTypes) ){
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
	    foreach ( &ExtraUtils::sortedKeys(\%modulesRequired) );
	$functionCode .= "    implicit none\n";
	$functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
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
	    foreach my $propertyName ( &ExtraUtils::sortedKeys($component->{'properties'}->{'property'}) ) {
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
		    } elsif ( $property->{'isVirtual'} && $property->{'attributes'}->{'isGettable'} ) {
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
		    if (&Utils::isOutputIntrinsic($type)) {
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
	    foreach ( &ExtraUtils::sortedKeys(\%modulesRequired) );
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
	    foreach my $propertyName ( &ExtraUtils::sortedKeys($component->{'properties'}->{'property'}) ) {
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
		    } elsif ( $property->{'isVirtual'} && $property->{'attributes'}->{'isGettable'} ) {
			$rank = $property->{'rank'};
			$type = $property->{'type'};
		    } else {
			die("Generate_Implementation_Output_Functions(): can not output [".$propertyName."]");
		    }		   
		    # Increment the counters.
		    if (&Utils::isOutputIntrinsic($type)) {
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
	$functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
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

sub Generate_Implementation_Name_From_Index_Functions {
    # Generate serialization/deserialization functions for each component implementation.
    my $build = shift;
    # Initialize function code.
    my $functionCode;
    # Initialize data content.
    my @dataContent;
    # Iterate over component implementations.
    foreach my $componentID ( @{$build->{'componentIdList'}} ) {
	# Get the component.
	my $component = $build->{'components'}->{$componentID};
	# Generate data content.
	@dataContent =
	    (
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($componentID),
		 attributes => [ "intent(in   )" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => [ "intent(inout)" ],
		 variables  => [ "count" ]
	     },
	     {
		 intrinsic  => "type",
		 type       => "varying_string",	
		 attributes => [ "intent(inout)" ],
		 variables  => [ "name" ]
	     }
	    );
	# Generate the function.
	my $selfUsed  = 0;
	my $countUsed = 0;
  	$functionCode  = "  subroutine Node_Component_".ucfirst($componentID)."_Name_From_Index(self,count,name)\n";
	$functionCode .= "    !% Return the name of the property of given index for a ".$component->{'name'}." implementation of the ".$component->{'class'}." component.\n";
	$functionCode .= "    use ISO_Varying_String\n";
	$functionCode .= "    implicit none\n";
	$functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
	my $functionBody = "";
	# If this component is an extension, first call on the extended type.
	if ( exists($build->{'components'}->{$componentID}->{'extends'}) ) {
	    my $extends = $build->{'components'}->{$componentID}->{'extends'};
	    $functionBody .= "    call self%nodeComponent".ucfirst($extends->{'class'}).ucfirst($extends->{'name'})."%nameFromIndex(count,name)\n";
	    $functionBody .= "    if (count <= 0) return\n";	    
	    $selfUsed  = 1;
	    $countUsed = 1;
	}
	# Iterate over properties.
	foreach my $propertyName ( &ExtraUtils::sortedKeys($component->{'properties'}->{'property'}) ) {
	    my $property = $component->{'properties'}->{'property'}->{$propertyName};
   	    # Check if this property has any linked data in this component.
	    if ( exists($property->{'linkedData'}) ) {
		# For each linked datum count if necessary.
		my $linkedDataName = $property->{'linkedData'};
		my $linkedData     = $component->{'content'}->{'data'}->{$linkedDataName};
		if ( $linkedData->{'isEvolvable'} ) {
		    $countUsed = 1;
		    if ( $linkedData->{'rank'} == 0 ) {
			if ( $linkedData->{'type'} eq "double" ) {
			    $functionBody .= "    count=count-1\n";
			}
			else {
			    $functionBody .= "    count=count-self%".&Utils::padLinkedData($linkedDataName,[0,0])."%serializeCount()\n";
			    $selfUsed  = 1;
			}
		    } else {
			$functionBody .= "    if (allocated(self%".&Utils::padLinkedData($linkedDataName,[0,0]).")) count=count-size(self%".&Utils::padLinkedData($linkedDataName,[0,0]).")\n";
			$selfUsed  = 1;
		    }
		    $functionBody .= "    if (count <= 0) then\n";
		    $functionBody .= "      name='".$component->{'class'}.":".$component->{'name'}.":".$propertyName."'\n";
		    $functionBody .= "      return\n";
		    $functionBody .= "    end if\n";
		}
	    }
	}
	$functionBody .= "    name='?'\n";
	$functionBody .= "    return\n";
	$functionBody .= "  end subroutine Node_Component_".ucfirst($componentID)."_Name_From_Index\n\n";
	$functionCode .= "   !GCC\$ attributes unused :: self\n"
	    unless ( $selfUsed  );
	$functionCode .= "   !GCC\$ attributes unused :: count\n"
	    unless ( $countUsed );
	$functionCode .= $functionBody;
	# Insert into the function list.
	push(
	    @{$build->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Insert a type-binding for this function into the implementation type.
	push(
	    @{$build->{'types'}->{'nodeComponent'.ucfirst($componentID)}->{'boundFunctions'}},
	    {type => "procedure", name => "nameFromIndex", function => "Node_Component_".ucfirst($componentID)."_Name_From_Index"}
	    );
    }
    # Generate data content.
    @dataContent =
	(
	 {
	     intrinsic  => "class",
	     type       => "nodeComponent",
	     attributes => [ "intent(in   )" ],
	     variables  => [ "self" ]
	 },
	 {
	     intrinsic  => "integer",
	     attributes => [ "intent(inout)" ],
	     variables  => [ "count" ]
	 },
	 {
	     intrinsic  => "type",
	     type       => "varying_string",	
	     attributes => [ "intent(inout)" ],
	     variables  => [ "name" ]
	 }
	);
    # Generate the function.
    $functionCode  = "  subroutine Node_Component_Name_From_Index(self,count,name)\n";
    $functionCode .= "    !% Return the name of the property of given index.\n";
    $functionCode .= "    use ISO_Varying_String\n";
    $functionCode .= "    implicit none\n";
    $functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
    $functionCode .= "    !GCC\$ attributes unused :: self, count\n";
    $functionCode .= "    name='?'\n";
    $functionCode .= "    return\n";
    $functionCode .= "  end subroutine Node_Component_Name_From_Index\n\n";
    # Insert into the function list.
    push(
	@{$build->{'code'}->{'functions'}},
	$functionCode
	);
    # Insert a type-binding for this function into the implementation type.
    push(
	@{$build->{'types'}->{'nodeComponent'}->{'boundFunctions'}},
	{type => "procedure", name => "nameFromIndex", function => "Node_Component_Name_From_Index", description => "Return the name of a property given is index.", returnType => "\\void", arguments => "\\intzero\\ count\\argin, \\textcolor{red}{\\textless varying\\_string\\textgreater}name\\argout"}
	);
}

sub Generate_Implementation_Serialization_Functions {
    # Generate serialization/deserialization functions for each component implementation.
    my $build = shift;
    # Initialize function code.
    my $functionCode;
    # Initialize data content.
    my @dataContent;
    # Iterate over component implementations.
    foreach my $componentID ( @{$build->{'componentIdList'}} ) {
	# Get the component.
	my $component = $build->{'components'}->{$componentID};
	# Generate data content.
	@dataContent =
	    (
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($componentID),
		 attributes => [ "intent(in   )" ],
		 variables  => [ "self" ]
	     }
	    );
	# Generate a count function.
  	$functionCode  = "  integer function Node_Component_".ucfirst($componentID)."_Count(self)\n";
	$functionCode .= "    !% Return a count of the serialization of a ".$component->{'name'}." implementation of the ".$component->{'class'}." component.\n";
	$functionCode .= "    implicit none\n";
	$functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
	my $functionBody = "";
	my $selfUsed = 0;
	# If this component is an extension, get the count of the extended type.
	$functionBody .= "    Node_Component_".ucfirst($componentID)."_Count=";
	if ( exists($build->{'components'}->{$componentID}->{'extends'}) ) {
	    my $extends = $build->{'components'}->{$componentID}->{'extends'};
	    $functionBody .= "self%nodeComponent".ucfirst($extends->{'class'}).ucfirst($extends->{'name'})."%serializeCount()\n";
	    $selfUsed = 1;
	} else {
	    $functionBody .= "0\n";
	}
	# Initialize a count of scalar properties.
	my $scalarPropertyCount = 0;
	# Iterate over properties.
	foreach my $propertyName ( &ExtraUtils::sortedKeys($component->{'properties'}->{'property'}) ) {
	    my $property = $component->{'properties'}->{'property'}->{$propertyName};
   	    # Check if this property has any linked data in this component.
	    if ( exists($property->{'linkedData'}) ) {
		# For each linked datum count if necessary.
		my $linkedDataName = $property->{'linkedData'};
		my $linkedData     = $component->{'content'}->{'data'}->{$linkedDataName};
		if ( $linkedData->{'isEvolvable'} ) {
		    if ( $linkedData->{'rank'} == 0 ) {
			if ( $linkedData->{'type'} eq "double" ) {
			    ++$scalarPropertyCount;
			}
			else {
			    $functionBody .= "    Node_Component_".ucfirst($componentID)."_Count=Node_Component_".ucfirst($componentID)."_Count+self%".&Utils::padLinkedData($linkedDataName,[0,0])."%serializeCount()\n";
			    $selfUsed = 1;
			}
		    } else {
			$functionBody .= "    if (allocated(self%".&Utils::padLinkedData($linkedDataName,[0,0]).")) Node_Component_".ucfirst($componentID)."_Count=Node_Component_".ucfirst($componentID)."_Count+size(self%".&Utils::padLinkedData($linkedDataName,[0,0]).")\n";
			$selfUsed = 1;
		    }
		}
	    }
	}
	# Insert the final count of scalar properties.
	$functionBody .= "    Node_Component_".ucfirst($componentID)."_Count=Node_Component_".ucfirst($componentID)."_Count+".$scalarPropertyCount."\n"
	    if ($scalarPropertyCount > 0);
	$functionBody .= "    return\n";
	$functionBody .= "  end function Node_Component_".ucfirst($componentID)."_Count\n\n";
	$functionCode .= "   !GCC\$ attributes unused :: self\n"
	    unless ( $selfUsed  );		    
	$functionCode .= $functionBody;
	# Insert into the function list.
	push(
	    @{$build->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Insert a type-binding for this function into the treeNode type.
	push(
	    @{$build->{'types'}->{'nodeComponent'.ucfirst($componentID)}->{'boundFunctions'}},
	    {type => "procedure", name => "serializeCount", function => "Node_Component_".ucfirst($componentID)."_Count"}
	    );
	# Specify data content for serialization functions.
	@dataContent =
	    (
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($componentID),
		 attributes => [ "intent(in   )" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "double precision",
		 attributes => [ "intent(  out)", "dimension(:)" ],
		 variables  => [ "array" ]
	     }
	    );
	# Generate serialization function.
	$functionCode  = "  subroutine Node_Component_".ucfirst($componentID)."_Serialize_Values(self,array)\n";
	$functionCode .= "    !% Serialize values of a ".$component->{'name'}." implementation of the ".$component->{'class'}." component.\n";
	$functionCode .= "    implicit none\n";
	my $serializationCode;
	my $needCount = 0;
	my $arrayUsed = 0;
	$selfUsed     = 0;
	# If this component is an extension, call serialization on the extended type.
	if ( exists($build->{'components'}->{$componentID}->{'extends'}) ) {
	    my $extends = $build->{'components'}->{$componentID}->{'extends'};
	    $serializationCode .= " count=self%nodeComponent".ucfirst($extends->{'class'}).ucfirst($extends->{'name'})."%serializeCount()\n";
	    $serializationCode .= " if (count > 0) then\n";
	    $serializationCode .= "  call self%nodeComponent".ucfirst($extends->{'class'}).ucfirst($extends->{'name'})."%serializeValues(array)\n";
	    $serializationCode .= "  offset=offset+count\n";
	    $serializationCode .= " end if\n";
	    $needCount = 1;
	    $arrayUsed = 1;
	    $selfUsed  = 1;
	}
	foreach my $propertyName ( &ExtraUtils::sortedKeys($component->{'properties'}->{'property'}) ) {
	    my $property = $component->{'properties'}->{'property'}->{$propertyName};
	    # Check if this property has any linked data in this component.
	    if ( exists($property->{'linkedData'}) ) {
		my $linkedDataName = $property->{'linkedData'};
		my $linkedData     = $component->{'content'}->{'data'}->{$linkedDataName};
		if ( $linkedData->{'isEvolvable'} ) {
		    $arrayUsed = 1;
		    $selfUsed  = 1;
		    if ( $linkedData->{'rank'} == 0 ) {
			if ( $linkedData->{'type'} eq "double" ) {
			    $serializationCode .= "    write (0,*) 'DEBUG -> Node_Component_".ucfirst($componentID)."_Serialize_Values -> ".$linkedDataName."',offset,size(array)\n"
				if ( $debugging == 1 );
			    $serializationCode .= "    array(offset)=self%".&Utils::padLinkedData($linkedDataName,[0,0])."\n";
			    $serializationCode .= "    offset=offset+1\n";
			}
			else {
			    $serializationCode .= "    count=self%".&Utils::padLinkedData($linkedDataName,[0,0])."%serializeCount(                            )\n";
			    $serializationCode .= "    write (0,*) 'DEBUG -> Node_Component_".ucfirst($componentID)."_Serialize_Values -> ".$linkedDataName."',offset,count,size(array)\n"
				if ( $debugging == 1 );
			    $serializationCode .= "    if (count > 0) call  self%".&Utils::padLinkedData($linkedDataName,[0,0])."%serialize     (array(offset:offset+count-1))\n";
			    $serializationCode .= "    offset=offset+count\n";
			    $needCount = 1;
			}
		    } else {
			$serializationCode .= "    if (allocated(self%".&Utils::padLinkedData($linkedDataName,[0,0]).")) then\n";
			$serializationCode .= "       count=size(self%".&Utils::padLinkedData($linkedDataName,[0,0]).")\n";
			$serializationCode .= "    write (0,*) 'DEBUG -> Node_Component_".ucfirst($componentID)."_Serialize_Values -> ".$linkedDataName."',offset,count,size(array)\n"
			    if ( $debugging == 1 );
			$serializationCode .= "       array(offset:offset+count-1)=reshape(self%".&Utils::padLinkedData($linkedDataName,[0,0]).",[count])\n";
			$serializationCode .= "       offset=offset+count\n";
			$serializationCode .= "    end if\n";
			$needCount = 1;
		    }
		}
	    }
	}
	if ( defined($serializationCode) ) {
	    my @variables = ( "offset" );
	    push(@variables,"count") 
		if ( $needCount == 1 );
	    push(
		@dataContent,
		{
		    intrinsic  => "integer",
		    variables  => \@variables
		}    
		);
	    $serializationCode = "    offset=1\n".$serializationCode;
	}
	$functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
	$functionCode .= "   !GCC\$ attributes unused :: array\n"
	    unless ( $arrayUsed );		    
	$functionCode .= "   !GCC\$ attributes unused :: self\n"
	    unless ( $selfUsed  );		    
	$functionCode .= $serializationCode
	    if ( defined($serializationCode) );
	$functionCode .= "    return\n";
	$functionCode .= "  end subroutine Node_Component_".ucfirst($componentID)."_Serialize_Values\n\n";
	# Insert into the function list.
	push(
	    	@{$build->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Insert a type-binding for this function into the implementation type.
	push(
	    @{$build->{'types'}->{'nodeComponent'.ucfirst($componentID)}->{'boundFunctions'}},
	    {type => "procedure", name => "serializeValues", function => "Node_Component_".ucfirst($componentID)."_Serialize_Values"},
	    );
	# Specify data content for deserialization functions.
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
		 attributes => [ "intent(in   )", "dimension(:)" ],
		 variables  => [ "array" ]
	     }
	    );
	# Generate deserialization function.
	$functionCode  = "  subroutine Node_Component_".ucfirst($componentID)."_Deserialize_Values(self,array)\n";
	$functionCode .= "    !% Serialize values of a ".$component->{'name'}." implementation of the ".$component->{'class'}." component.\n";
	$functionCode .= "    implicit none\n";
	my $deserializationCode;
	$selfUsed  = 0;
	$arrayUsed = 0;
	$needCount = 0;
	# If this component is an extension, call deserialization on the extended type.
	if ( exists($build->{'components'}->{$componentID}->{'extends'}) ) {
	    my $extends = $build->{'components'}->{$componentID}->{'extends'};
	    $deserializationCode .= " count=self%nodeComponent".ucfirst($extends->{'class'}).ucfirst($extends->{'name'})."%serializeCount()\n";
	    $deserializationCode .= " if (count > 0) then\n";
	    $deserializationCode .= "  call self%nodeComponent".ucfirst($extends->{'class'}).ucfirst($extends->{'name'})."%deserializeValues(array)\n";
	    $deserializationCode .= "  offset=offset+count\n";
	    $deserializationCode .= " end if\n";
	    $needCount = 1;
	    $arrayUsed = 1;
	    $selfUsed  = 1;
	}
	foreach my $propertyName ( &ExtraUtils::sortedKeys($component->{'properties'}->{'property'}) ) {
	    my $property = $component->{'properties'}->{'property'}->{$propertyName};
	    # Check if this property has any linked data in this component.
	    if ( exists($property->{'linkedData'}) ) {
		# For each linked datum count if necessary.
		my $linkedDataName = $property->{'linkedData'};
		my $linkedData     = $component->{'content'}->{'data'}->{$linkedDataName};
		if ( $linkedData->{'isEvolvable'} ) {
		    $arrayUsed = 1;
		    $selfUsed  = 1;
		    if ( $linkedData->{'rank'} == 0 ) {
			if ( $linkedData->{'type'} eq  "double" ) {
			    $deserializationCode .= "    self%".&Utils::padLinkedData($linkedDataName,[0,0])."=array(offset)\n";
			    $deserializationCode .= "    offset=offset+1\n";
			}
			else {
			    $deserializationCode .= "    count=self%".&Utils::padLinkedData($linkedDataName,[0,0])."%serializeCount(                            )\n";
			    $deserializationCode .= "    call  self%".&Utils::padLinkedData($linkedDataName,[0,0])."%deserialize   (array(offset:offset+count-1))\n";
			    $deserializationCode .= "    offset=offset+count\n";
			    $needCount = 1;
			}
		    } else {
			$deserializationCode .= "    if (allocated(self%".&Utils::padLinkedData($linkedDataName,[0,0]).")) then\n";
			$deserializationCode .= "       count=size(self%".&Utils::padLinkedData($linkedDataName,[0,0]).")\n";
			$deserializationCode .= "       self%".&Utils::padLinkedData($linkedDataName,[0,0])."=reshape(array(offset:offset+count-1),shape(self%".&Utils::padLinkedData($linkedDataName,[0,0])."))\n";
			$deserializationCode .= "       offset=offset+count\n";
			$deserializationCode .= "    end if\n";
			$needCount = 1;
		    }
		}
	    }
	}
	if ( defined($deserializationCode) ) {
	    my @variables = ( "offset" );
	    push(@variables,"count") 
		if ( $needCount == 1 );
	    push(
		@dataContent,
		{
		    intrinsic  => "integer",
		    variables  => \@variables
		}    
		);
	    $deserializationCode = "    offset=1\n".$deserializationCode;
	}
	$functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
	$functionCode .= "   !GCC\$ attributes unused :: array\n"
	    unless ( $arrayUsed );		    
	$functionCode .= "   !GCC\$ attributes unused :: self\n"
	    unless ( $selfUsed  );		    
	$functionCode .= $deserializationCode
	    if ( defined($deserializationCode) );
	$functionCode .= "    return\n";
	$functionCode .= "  end subroutine Node_Component_".ucfirst($componentID)."_Deserialize_Values\n\n";
	# Insert into the function list.
	push(
	    @{$build->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Insert a type-binding for this function into the implementation type.
	push(
	    @{$build->{'types'}->{'nodeComponent'.ucfirst($componentID)}->{'boundFunctions'}},
	    {type => "procedure", name => "deserializeValues", function => "Node_Component_".ucfirst($componentID)."_Deserialize_Values"},
	    );
    }
}

sub Generate_Serialization_Offset_Variables {
    # Generate variables which store offsets into arrays for serialization.
    my $build = shift;
    # Create a table.
    my $offsetTable = Text::Table->new(
	{
	    is_sep => 1,
	    body   => "  integer :: "
	},
	{
	    align  => "left"
	}
	);
    my $privateTable = Text::Table->new(
	{
	    is_sep => 1,
	    body   => "  !\$omp threadprivate("
	},
	{
	    align  => "left"
	},
	{
	    is_sep  => 1,
	    body    => ")"
	}
	);
    # Iterate over component implementations.
    foreach my $componentID ( @{$build->{'componentIdList'}} ) {
	# Get the component.
	my $component = $build->{'components'}->{$componentID};
	# Iterate over properties.
	foreach my $propertyName ( &ExtraUtils::sortedKeys($component->{'properties'}->{'property'}) ) {
	    my $property = $component->{'properties'}->{'property'}->{$propertyName};
   	    # Check if this property has any linked data in this component.
	    if ( exists($property->{'linkedData'}) ) {
		# For each linked datum count if necessary.
		my $linkedDataName = $property->{'linkedData'};
		my $linkedData     = $component->{'content'}->{'data'}->{$linkedDataName};
		if ( $linkedData->{'isEvolvable'} ) {
		    my $offsetName = &offsetName($componentID,$propertyName);
		    $offsetTable ->add($offsetName);
		    $privateTable->add($offsetName);
		}
	    }
	}
    }
    # Insert into the document.
    $build->{'content'} .= "  ! Offsets into serialization arrays.\n";
    $build->{'content'} .= $offsetTable ->table()."\n";
    $build->{'content'} .= $privateTable->table()."\n";
    $build->{'content'} .= " integer                                     :: nodeSerializationCount\n";
    $build->{'content'} .= " double precision, allocatable, dimension(:) :: nodeScales, nodeRates, nodeRatesIncrement\n";
    $build->{'content'} .= " !\$omp threadprivate(nodeScales,nodeRates,nodeRatesIncrement,nodeSerializationCount)\n";
}

sub offsetName {
    my $componentName = shift();
    my $propertyName  = shift();
    return "offset".ucfirst($componentName).ucfirst($propertyName);
}

sub Generate_Node_Offset_Functions {
    # Generate functions to compute offsets into serialization arrays.
    my $build = shift;

    # Function computing a count of the serialization length.
    my @dataContent =
	(
	 {
	     intrinsic  => "class",
	     type       => "treeNode",
	     attributes => [ "intent(in   )" ],
	     variables  => [ "self" ]
	 },
	 {
	     intrinsic  => "integer",
	     variables  => [ "i", "count" ]
	 }
	);
    my $functionCode;
    $functionCode .= "  subroutine SerializationOffsets(self)\n";
    $functionCode .= "    !% Compute offsets into serialization arrays for {\\normalfont \\ttfamily treeNode} object.\n";
    $functionCode .= "    implicit none\n";
    $functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
    $functionCode .= "    count=0\n";
    # Loop over all component classes
    foreach ( @{$build->{'componentClassList'}} ) {	    
     	$functionCode .= "    if (allocated(self%component".&Utils::padClass(ucfirst($_),[0,0]).")) then\n";
	$functionCode .= "      do i=1,size(self%component".&Utils::padClass(ucfirst($_),[0,0]).")\n";
	$functionCode .= "        call self%component".&Utils::padClass(ucfirst($_),[0,0])."(i)%serializationOffsets(count)\n";
	$functionCode .= "      end do\n";
	$functionCode .= "    end if\n";
    }
    $functionCode .= "    if (.not.allocated(nodeScales)) then\n";
    $functionCode .= "       allocate  (nodeScales        (count))\n";
    $functionCode .= "       allocate  (nodeRates         (count))\n";
    $functionCode .= "       allocate  (nodeRatesIncrement(count))\n";
    $functionCode .= "    else if (size(nodeScales) < count) then\n";
    $functionCode .= "       deallocate(nodeScales               )\n";
    $functionCode .= "       deallocate(nodeRates                )\n";
    $functionCode .= "       deallocate(nodeRatesIncrement       )\n";
    $functionCode .= "       allocate  (nodeScales        (count))\n";
    $functionCode .= "       allocate  (nodeRates         (count))\n";
    $functionCode .= "       allocate  (nodeRatesIncrement(count))\n";
    $functionCode .= "    end if\n";
    $functionCode .= "    nodeSerializationCount=count\n";
    $functionCode .= "    return\n";
    $functionCode .= "  end subroutine SerializationOffsets\n\n";
    # Insert into the function list.
    push(
	@{$build->{'code'}->{'functions'}},
	$functionCode
	);
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{type => "procedure", name => "serializationOffsets", function => "SerializationOffsets", description => "Compute offsets into serialization arrays for all properties", returnType => "\\void", arguments => ""}
	);
}

sub Generate_Implementation_Offset_Functions {
    # Generate serialization offset functions for each component implementation.
    my $build = shift;
    # Initialize function code.
    my $functionCode;
    # Initialize data content.
    my @dataContent;
    # Iterate over component implementations.
    foreach my $componentID ( @{$build->{'componentIdList'}} ) {
	# Get the component.
	my $component = $build->{'components'}->{$componentID};
	# Generate data content.
	@dataContent =
	    (
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($componentID),
		 attributes => [ "intent(in   )" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => [ "intent(inout)" ],
		 variables  => [ "count" ]
	     }
	    );
	# Generate a count function.
  	$functionCode  = "  subroutine Node_Component_".ucfirst($componentID)."_Offsets(self,count)\n";
	$functionCode .= "    !% Return a count of the serialization of a ".$component->{'name'}." implementation of the ".$component->{'class'}." component.\n";
	$functionCode .= "    implicit none\n";
	$functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
	my $functionBody = "";
	# If this component is an extension, compute offsets of the extended type.
	if ( exists($build->{'components'}->{$componentID}->{'extends'}) ) {
	    my $extends = $build->{'components'}->{$componentID}->{'extends'};
	    $functionBody .= "call self%nodeComponent".ucfirst($extends->{'class'}).ucfirst($extends->{'name'})."%serializationOffsets(count)\n";
	}
	# Iterate over properties.
	my $countUsed = 0;
	my $selfUsed  =0;
	foreach my $propertyName ( &ExtraUtils::sortedKeys($component->{'properties'}->{'property'}) ) {
	    my $property = $component->{'properties'}->{'property'}->{$propertyName};
   	    # Check if this property has any linked data in this component.
	    if ( exists($property->{'linkedData'}) ) {
		# For each linked datum count if necessary.
		my $linkedDataName = $property->{'linkedData'};
		my $linkedData     = $component->{'content'}->{'data'}->{$linkedDataName};
		if ( $linkedData->{'isEvolvable'} ) {
		    $countUsed = 1;
		    my $offsetName = &offsetName($componentID,$propertyName);
		    $functionBody .= "    ".$offsetName."=count+1\n";
		    if ( $linkedData->{'rank'} == 0 ) {
			if ( $linkedData->{'type'} eq "double" ) {
			    $functionBody .= "    count=count+1\n";
			}
			else {
			    $selfUsed = 1;
			    $functionBody .= "    count=count+self%".&Utils::padLinkedData($linkedDataName,[0,0])."%serializeCount()\n";
			}
		    } else {
			$selfUsed = 1;
			$functionBody .= "    if (allocated(self%".&Utils::padLinkedData($linkedDataName,[0,0]).")) count=count+size(self%".&Utils::padLinkedData($linkedDataName,[0,0]).")\n";
		    }
		}
	    }
	}
	$functionBody .= "    return\n";
	$functionBody .= "  end subroutine Node_Component_".ucfirst($componentID)."_Offsets\n\n";
	$functionCode .= "   !GCC\$ attributes unused :: count\n"
	    unless ( $countUsed );		    
	$functionCode .= "   !GCC\$ attributes unused :: self\n"
	    unless ( $selfUsed  );		    
	$functionCode .= $functionBody;	    
	# Insert into the function list.
	push(
	    @{$build->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Insert a type-binding for this function into the treeNode type.
	push(
	    @{$build->{'types'}->{'nodeComponent'.ucfirst($componentID)}->{'boundFunctions'}},
	    {type => "procedure", name => "serializationOffsets", function => "Node_Component_".ucfirst($componentID)."_Offsets"}
	    );
    }
}

sub Generate_Component_Count_Functions {
    # Generate component count functions.
    my $build = shift;
    # Initialize function code.
    my $functionCode;
    # Initialize data content.
    my @dataContent;
    # Create methods to get components.
    foreach my $componentClassName ( @{$build->{'componentClassList'}} ) {
	# Specify data content for component get function.
	@dataContent =
	    (
	     {
		 intrinsic  => "class",
		 type       => "treeNode",
		 attributes => [ "intent(inout)", "target" ],
		 variables  => [ "self" ]
	     }
	    );
	# Generate code for the count function.
    	$functionCode  = "  integer function ".$componentClassName."CountLinked(self)\n";
	$functionCode .= "    !% Returns the number of {\\normalfont \\ttfamily ".$componentClassName."} components in {\\normalfont \\ttfamily self}.\n";
    	$functionCode .= "    implicit none\n";
	$functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
    	$functionCode .= "    select type (self)\n";
    	$functionCode .= "    class is (treeNode)\n";
	$functionCode .= "     if (allocated(self%component".ucfirst($componentClassName).")) then\n";
	$functionCode .= "       select type (component => self%component".ucfirst($componentClassName)."(1))\n";
	$functionCode .= "       type is (nodeComponent".ucfirst($componentClassName).")\n";
	$functionCode .= "         ".$componentClassName."CountLinked=0\n";
	$functionCode .= "       class default\n";
	$functionCode .= "         ".$componentClassName."CountLinked=size(self%component".ucfirst($componentClassName).")\n";
	$functionCode .= "       end select\n";
	$functionCode .= "     else\n";
	$functionCode .= "        ".$componentClassName."CountLinked=0\n";
	$functionCode .= "     end if\n";
	$functionCode .= "    class default\n";
	$functionCode .= "     ".$componentClassName."CountLinked=0\n";
	$functionCode .= "     call Galacticus_Error_Report('".$componentClassName."CountLinked\','treeNode of unknown class')\n";
    	$functionCode .= "    end select\n";
    	$functionCode .= "    return\n";
    	$functionCode .= "  end function ".$componentClassName."CountLinked\n";
	# Insert into the function list.
	push(
	    @{$build->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Bind this function to the treeNode type.
	push(
	    @{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	    {type => "procedure", name => $componentClassName."Count", function => $componentClassName."CountLinked", description => "Returns the number of {\\normalfont \\ttfamily ".$componentClassName."} components in the node.", returnType => "\\intzero", arguments => ""}
	    );
    }
}

sub Generate_Component_Get_Functions {
    # Generate component get methods.
    my $build = shift;
    # Initialize function code.
    my $functionCode;
    # Initialize data content.
    my @dataContent;
    # Create methods to get components.
    foreach my $componentClassName ( @{$build->{'componentClassList'}} ) {
	# Specify data content for component get function.
	@dataContent =
	    (
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($componentClassName),
		 attributes => [ "pointer" ],
		 variables  => [ $componentClassName."Get" ]
	     },
	     {
		 intrinsic  => "class",
		 type       => "treeNode",
		 attributes => [ "intent(inout)", "target" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => [ "intent(in   )", "optional" ],
		 variables  => [ "instance" ]
	     },
	     {
		 intrinsic  => "logical",
		 attributes => [ "intent(in   )", "optional" ],
		 variables  => [ "autoCreate" ]
	     },
	     {
		 intrinsic  => "integer",
		 variables  => [ "instanceActual" ]
	     },
	     {
		 intrinsic  => "logical",
		 variables  => [ "autoCreateActual" ]
	     }
	    );
	# Generate code for the get function.
    	$functionCode  = "  recursive function ".$componentClassName."Get(self,instance,autoCreate)\n";
	$functionCode .= "    !% Returns the {\\normalfont \\ttfamily ".$componentClassName."} component of {\\normalfont \\ttfamily self}.\n";
	$functionCode .= "    use Galacticus_Error\n";   
 	$functionCode .= "    implicit none\n";
	$functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
    	$functionCode .= "    instanceActual=1\n";
    	$functionCode .= "    if (present(instance)) instanceActual=instance\n";
    	$functionCode .= "    autoCreateActual=.false.\n";
    	$functionCode .= "    if (present(autoCreate)) autoCreateActual=autoCreate\n";
	$functionCode .= "    if (autoCreateActual.and.allocated(self%component".ucfirst($componentClassName).")) then\n";
	# If we are allowed to autocreate the component and it still has generic type then deallocate it to
	# force it to be created later.
	$functionCode .= "      if (same_type_as(self%component".ucfirst($componentClassName)."(1),".ucfirst($componentClassName)."Class)) deallocate(self%component".ucfirst($componentClassName).")\n";
	$functionCode .= "    end if\n";
    	$functionCode .= "    if (.not.allocated(self%component".ucfirst($componentClassName).")) then\n";
    	$functionCode .= "      if (autoCreateActual) then\n";
	$functionCode .= "         call self%".lc($componentClassName)."Create"."()\n";
	$functionCode .= "      else\n";
	$functionCode .= "         call Galacticus_Error_Report('".$componentClassName."Get','component is not allocated')\n";
	$functionCode .= "      end if\n";
    	$functionCode .= "    end if\n";
    	$functionCode .= "    ".$componentClassName."Get => self%component".ucfirst($componentClassName)."(instanceActual)\n";
    	$functionCode .= "    return\n";
    	$functionCode .= "  end function ".$componentClassName."Get\n";
	# Insert into the function list.
	push(
	    @{$build->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Bind this function to the treeNode type.
	push(
	    @{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	    {type => "procedure", name => $componentClassName, function => $componentClassName."Get", description => "Return a ".$componentClassName." component member of the node. If no {\\normalfont \\ttfamily instance} is specified, return the first instance. If {\\normalfont \\ttfamily autoCreate} is {\\normalfont \\ttfamily true} then create a single instance of the component if none exists in the node.", returnType => "\\textcolor{red}{\\textless *class(nodeComponent".ucfirst($componentClassName).")\\textgreater}", arguments => "\\intzero\\ [instance]\\argin, \\logicalzero\\ [autoCreate]\\argin"}
	    );
	# Specify data content for create-by-interrupt function.
	@dataContent =
	    (
	     {
		 intrinsic  => "type",
		 type       => "treeNode",
		 attributes => [ "pointer", "intent(inout)" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($componentClassName),
		 attributes => [ "pointer" ],
		 variables  => [ $componentClassName ]
	     }
	    );
	# Generate function to create component via an interrupt.
	my $createIfNeeded = 0;
	# Iterate over component implementations.
	foreach my $implementationName ( @{$build->{'componentClasses'}->{$componentClassName}->{'memberNames'}} ) {
	    my $implementationID = ucfirst($componentClassName).ucfirst($implementationName);
	    my $component = $build->{'components'}->{$implementationID};
	    # Iterate over properties.
	    foreach my $propertyName ( keys(%{$component->{'properties'}->{'property'}}) ) {
		my $property = $component->{'properties'}->{'property'}->{$propertyName};
		$createIfNeeded = 1
		    if ( $property->{'attributes'}->{'createIfNeeded'} );
	    }
	}
	if ( $createIfNeeded ) {
	    $functionCode  = "  subroutine ".$componentClassName."CreateByInterrupt(self)\n";
	    $functionCode .= "    !% Create the {\\normalfont \\ttfamily ".$componentClassName."} component of {\\normalfont \\ttfamily self} via an interrupt.\n";
	    $functionCode .= "    implicit none\n";
	    $functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
	    $functionCode .= "    ".$componentClassName." => self%".$componentClassName."(autoCreate=.true.)\n";
	    # Loop over instances of this class, and call custom create routines if necessary.
	    my $foundCreateFunctions = 0;
    	foreach my $componentName ( @{$build->{'componentClasses'}->{$componentClassName}->{'memberNames'}} ) {
		my $componentID = ucfirst($componentClassName).ucfirst($componentName);
	    my $component = $build->{'components'}->{$componentID};
		if ( exists($component->{'createFunction'}) ) {
		    if ( $foundCreateFunctions == 0 ) {
			$functionCode .= "    select type (".$componentClassName.")\n";
			$foundCreateFunctions = 1;
		    }
		$functionCode .= "    type is (nodeComponent".&Utils::padFullyQualified(ucfirst($componentID),[0,0]).")\n";
		    my $createFunction = $component->{'createFunction'};
		    $createFunction = $component->{'createFunction'}->{'content'}
		    if ( exists($component->{'createFunction'}->{'content'}) );
		    $createFunction = $componentID."CreateFunction"
			if (
			exists($component->{'createFunction'}->{'isDeferred'})
			&&
			$component->{'createFunction'}->{'isDeferred'}
			);
		    $functionCode .= "       call ".$createFunction."(".$componentClassName.")\n";
		}
	    }
	    $functionCode .= "    end select\n"
		unless ( $foundCreateFunctions == 0 );
	    $functionCode .= "    return\n";
	    $functionCode .= "  end subroutine ".$componentClassName."CreateByInterrupt\n";
	    # Insert into the function list.
	    push(
	    @{$build->{'code'}->{'functions'}},
		$functionCode
		);
	}
	# If any create function is deferred, create a function to set it at runt time.
    	foreach my $componentName ( @{$build->{'componentClasses'}->{$componentClassName}->{'memberNames'}} ) {
	    my $componentID = ucfirst($componentClassName).ucfirst($componentName);
	    my $component = $build->{'components'}->{$componentID};
	    if (
		exists($component->{'createFunction'}                           ) && 
		exists($component->{'createFunction'}->{'isDeferred'}           ) &&
		$component->{'createFunction'}->{'isDeferred'}
		) {
		$functionCode  = "   subroutine ".$componentID."CreateFunctionSet(createFunction)\n";
		$functionCode .= "     !% Set the create function for the {\\normalfont \\ttfamily ".$componentID."} component.\n";
		$functionCode .= "     implicit none\n";
		$functionCode .= "     external createFunction\n";
		my $createFunction = $componentID."CreateFunction"; 
		$functionCode .= "     ".$createFunction." => createFunction\n";
		$functionCode .= "     return\n";
		$functionCode .= "   end subroutine ".$componentID."CreateFunctionSet\n";
		# Insert into the function list.
		push(
		    @{$build->{'code'}->{'functions'}},
		    $functionCode
		    );
	    }
	}
    }
}

sub Generate_Component_Destruction_Functions {
    # Generate component destruction functions.
    my $build = shift;
    # Iterate over component classes.
    foreach my $componentClassName ( @{$build->{'componentClassList'}} ) {
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
		 variables  => [ "i" ]
	     }
	    );
    	my $functionCode = "  subroutine ".$componentClassName."DestroyLinked(self)\n";
	$functionCode   .= "    !% Destroy the {\\normalfont \\ttfamily ".$componentClassName."} component of {\\normalfont \\ttfamily self}.\n";
    	$functionCode   .= "    implicit none\n";
	$functionCode   .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
    	$functionCode   .= "    if (allocated(self%component".ucfirst($componentClassName).")) then\n";
	$functionCode   .= "      do i=1,size(self%component".ucfirst($componentClassName).")\n";
	$functionCode   .= "        call        self%component".ucfirst($componentClassName)."(i)%destroy()\n";
	$functionCode   .= "      end do\n";
	$functionCode   .= "      deallocate (self%component".ucfirst($componentClassName).")\n";
	$functionCode   .= "    end if\n";
    	$functionCode   .= "    return\n";
    	$functionCode   .= "  end subroutine ".$componentClassName."DestroyLinked\n\n";
	# Insert into the function list.
	push(
	    @{$build->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Bind this function to the treeNode type.
	push(
	    @{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	    {type => "procedure", name => $componentClassName."Destroy" , function => $componentClassName."DestroyLinked", description => "Destroy the {\\normalfont \\ttfamily ".$componentClassName."} component(s) of the node.", returnType => "\\void", arguments => ""}
	    );
    }
}

sub Generate_Component_Creation_Functions {
    # Generate component creation functions.
    my $build = shift;
    # Iterate over component classes.
    foreach my $componentClassName ( @{$build->{'componentClassList'}} ) {
	# Specify data content.
	my @dataContent =
	    (
	     {
		 intrinsic  => "class",
		 type       => "treeNode",
		 attributes => [ "target", "intent(inout)" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($componentClassName),
		 attributes => [ "intent(in   )", "optional" ],
		 variables  => [ "template" ]
	     },
	     {
		 intrinsic  => "type",
		 type       => "varying_string",
		 variables  => [ "message" ]
	     },
	     {
		 intrinsic  => "integer",
		 variables  => [ "i" ]
	     }
	    );
	# Generate function code.
	my $functionCode;
    	$functionCode .= "  subroutine ".$componentClassName."CreateLinked(self,template)\n";
	$functionCode .= "    !% Create the {\\normalfont \\ttfamily ".$componentClassName."} component of {\\normalfont \\ttfamily self}.\n";
	$functionCode .= "    use ISO_Varying_String\n";
	$functionCode .= "    use Galacticus_Display\n";
	$functionCode .= "    use Galacticus_Error\n";
	$functionCode .= "    use String_Handling\n";
    	$functionCode .= "    implicit none\n";
 	$functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
	$functionCode .= "    if (Galacticus_Verbosity_Level() >= verbosityInfo) then\n";
	$functionCode .= "      message='Creating ".$componentClassName." in node '\n";
	$functionCode .= "      message=message//self%index()\n";
	$functionCode .= "      call Galacticus_Display_Message(message,verbosityInfo)\n";
	$functionCode .= "    end if\n";
    	$functionCode .= "    if (present(template)) then\n";
	$functionCode .= "       allocate(self%component".ucfirst($componentClassName)."(1),source=template)\n";
    	$functionCode .= "    else\n";
	$functionCode .= "       select type (default".ucfirst($componentClassName)."Component)\n";
	$functionCode .= "       type is (nodeComponent".ucfirst($componentClassName)."Null)\n";
	$functionCode .= "          call Galacticus_Error_Report('".$componentClassName."CreateLinked','refusing to create null instance')\n";
	$functionCode .= "       class default\n";
	$functionCode .= "          allocate(self%component".ucfirst($componentClassName)."(1),source="."default".ucfirst($componentClassName)."Component)\n";
   	$functionCode .= "       end select\n";
   	$functionCode .= "    end if\n";
     	$functionCode .= "    select type (self)\n";
	$functionCode .= "    type is (treeNode)\n";
	$functionCode .= "      do i=1,size(self%component".ucfirst($componentClassName).")\n";
	$functionCode .= "        self%component".ucfirst($componentClassName)."(i)%hostNode => self\n";
	$functionCode .= "        call self%component".ucfirst($componentClassName)."(i)%initialize()\n";
	$functionCode .= "      end do\n";
    	$functionCode .= "    end select\n";
    	$functionCode .= "    return\n";
    	$functionCode .= "  end subroutine ".$componentClassName."CreateLinked\n\n";
	# Insert into the function list.
	push(
	    @{$build->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Bind this function to the treeNode type.
	push(
	    @{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	    {type => "procedure", name => $componentClassName."Create" , function => $componentClassName."CreateLinked", description => "Create a {\\normalfont \\ttfamily ".$componentClassName."} component in the node. If no {\\normalfont \\ttfamily template} is specified use the active implementation of this class.", returnType => "\\void", arguments => "\\textcolor{red}{\\textless class(nodeComponent".ucfirst($componentClassName).")\\textgreater}\\ [template]\\argin"}
	    );
    }
}

sub Generate_Node_Copy_Function {
    # Generate function to copy one node to another.
    my $build = shift;
    # Specify variables.
    my @dataContent =
	(
	 {
	     intrinsic  => "class",
	     type       => "treeNode",
	     attributes => [ "intent(in   )" ],
	     variables  => [ "self" ]
	 },
	 {
	     intrinsic  => "class",
	     type       => "treeNode",
	     attributes => [ "intent(inout)" ],
	     variables  => [ "targetNode" ]
	 },
	 {
	     intrinsic  => "logical",
	     attributes => [ "intent(in   )", "optional" ],
	     variables  => [ "skipFormationNode" ]
	 },
	 {
	     intrinsic  => "logical",
	     variables  => [ "skipFormationNodeActual" ]
	 },
	 {
	     intrinsic  => "integer",
	     variables  => [ "i" ]
	 }
	);
    # Generate the code.
    my $functionCode;
    $functionCode .= "  subroutine Tree_Node_Copy_Node_To(self,targetNode,skipFormationNode)\n";
    $functionCode .= "    !% Make a copy of {\\normalfont \\ttfamily self} in {\\normalfont \\ttfamily targetNode}.\n";
    $functionCode .= "    implicit none\n";
    $functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
    $functionCode .= "    skipFormationNodeActual=.false.\n";
    $functionCode .= "    if (present(skipFormationNode)) skipFormationNodeActual=skipFormationNode\n";
    $functionCode .= "    targetNode%".&Utils::padClass($_,[8,14])." =  self%".$_."\n"
	foreach ( "uniqueIdValue", "indexValue", "timeStepValue" );
    $functionCode .= "    targetNode%".&Utils::padClass($_,[8,14])." => self%".$_."\n"
	foreach ( "parent", "firstChild", "sibling", "firstSatellite", "mergeTarget", "firstMergee", "siblingMergee", "event", "hostTree" );
    $functionCode .= "    if (.not.skipFormationNodeActual) targetNode%formationNode => self%formationNode\n";
    # Loop over all component classes
    if ( $workaround == 1 ) { # Workaround "Assignment to an allocatable polymorphic variable is not yet supported"
	foreach my $componentClassName ( @{$build->{'componentClassList'}} ) {
	    $functionCode .= "    if (allocated(targetNode%component".&Utils::padClass(ucfirst($componentClassName),[0,0]).")) deallocate(targetNode%component".&Utils::padClass(ucfirst($componentClassName),[0,0]).")\n";
	    $functionCode .= "    allocate(targetNode%component".&Utils::padClass(ucfirst($componentClassName),[0,0])."(size(self%component".&Utils::padClass(ucfirst($componentClassName),[0,0]).")),source=self%component".&Utils::padClass(ucfirst($componentClassName),[0,0])."(1))\n";
	    $functionCode .= "    do i=1,size(self%component".&Utils::padClass(ucfirst($componentClassName),[0,0]).")\n";
	    foreach my $implementationName ( @{$build->{'componentClasses'}->{$componentClassName}->{'memberNames'}} ) {
		$functionCode .= "      select type (from => self%component".&Utils::padClass(ucfirst($componentClassName),[0,0]).")\n";
		$functionCode .= "      type is (nodeComponent".&Utils::padFullyQualified(ucfirst($componentClassName).ucfirst($implementationName),[0,0]).")\n";
		$functionCode .= "        select type (to => targetNode%component".&Utils::padClass(ucfirst($componentClassName),[0,0]).")\n";
		$functionCode .= "        type is (nodeComponent".&Utils::padFullyQualified(ucfirst($componentClassName).ucfirst($implementationName),[0,0]).")\n";
		$functionCode .= "          to=from\n";
		$functionCode .= "        end select\n";
		$functionCode .= "      end select\n";
	    }
	    $functionCode .= "    end do\n";
	}
    } else {
	foreach ( @{$build->{'componentClassList'}} ) {
	    $functionCode .= "    targetNode%component".&Utils::padClass(ucfirst($_),[0,14])."=  self%component".ucfirst($_)."\n";
	}
    }
    # Update target node pointers.
    $functionCode .= "    select type (targetNode)\n";
    $functionCode .= "    type is (treeNode)\n";
    foreach ( @{$build->{'componentClassList'}} ) {
	$functionCode .= "      do i=1,size(self%component".&Utils::padClass(ucfirst($_),[0,0]).")\n";
	
	$functionCode .= "        targetNode%component".&Utils::padClass(ucfirst($_),[0,14])."(i)%hostNode =>  targetNode\n";
	$functionCode .= "      end do\n";
    }
    $functionCode .= "    end select\n";
    $functionCode .= "    return\n";
    $functionCode .= "  end subroutine Tree_Node_Copy_Node_To\n\n";
    # Insert into the function list.
    push(
	@{$build->{'code'}->{'functions'}},
	$functionCode
	);
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{type => "procedure", name => "copyNodeTo", function => "Tree_Node_Copy_Node_To", description => "Make a copy of the node in {\\normalfont \\ttfamily targetNode}. If {\\normalfont \\ttfamily skipFormationNode} is {\\normalfont \\ttfamily true} then do not copy any pointer to the formation node.", returnType => "\\void", arguments => "\\textcolor{red}{\\textless class(treeNode)\\textgreater} targetNode\\arginout, \\logicalzero\\ [skipFormationNode]\\argin"}
	);
}

sub Generate_Node_Move_Function {
    # Generate function to move one node to another.
    my $build = shift;
    # Specify variables.
    my @dataContent =
	(
	 {
	     intrinsic  => "class",
	     type       => "treeNode",
	     attributes => [ "intent(inout)" ],
	     variables  => [ "self" ]
	 },
	 {
	     intrinsic  => "type",
	     type       => "treeNode",
	     attributes => [ "intent(inout)", "target" ],
	     variables  => [ "targetNode" ]
	 },
	 {
	     intrinsic  => "integer",
	     variables  => [ "i" ]
	 }
	);
    # Generate the code.
    my $functionCode;
    # Create functions for moving node components.
    $functionCode .= "  subroutine Tree_Node_Move_Components(self,targetNode)\n";
    $functionCode .= "    !% Move components from {\\normalfont \\ttfamily self} to {\\normalfont \\ttfamily targetNode}.\n";
    $functionCode .= "    implicit none\n";
    $functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
    # Loop over all component classes
    foreach ( @{$build->{'componentClassList'}} ) {	    
	$functionCode .= "    if (allocated(targetNode%component".&Utils::padClass(ucfirst($_),[0,0]).")) then\n";
	$functionCode .= "      do i=1,size(targetNode%component".&Utils::padClass(ucfirst($_),[0,0]).")\n";
	$functionCode .= "        call targetNode%component".&Utils::padClass(ucfirst($_),[0,0])."(i)%destroy()\n";
	$functionCode .= "      end do\n";
	$functionCode .= "      deallocate(targetNode%component".&Utils::padClass(ucfirst($_),[0,0]).")\n";
	$functionCode .= "    end if\n";
	$functionCode .= "    if (allocated(self      %component".&Utils::padClass(ucfirst($_),[0,0]).")) then\n";
	$functionCode .= "       call Move_Alloc(self%component".&Utils::padClass(ucfirst($_),[0,0]).",targetNode%component".&Utils::padClass(ucfirst($_),[0,0]).")\n";
	$functionCode .= "       do i=1,size(targetNode%component".&Utils::padClass(ucfirst($_),[0,0]).")\n";
	$functionCode .= "         targetNode%component".&Utils::padClass(ucfirst($_),[0,0])."(i)%hostNode => targetNode\n";
	$functionCode .= "       end do\n";
	$functionCode .= "    end if\n";
    }
    $functionCode .= "    return\n";
    $functionCode .= "  end subroutine Tree_Node_Move_Components\n\n";
    # Insert into the function list.
    push(
	@{$build->{'code'}->{'functions'}},
	$functionCode
	);
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{type => "procedure", name => "moveComponentsTo", function => "Tree_Node_Move_Components", description => "Move components from a node to {\\normalfont \\ttfamily targetNode}.", returnType => "\\void", arguments => "\\textcolor{red}{\\textless class(treeNode)\\textgreater} targetNode\\arginout"}
	);
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
	(my $dataObject, my $label) = &DataTypes::dataObjectDefinition($property);
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
	$functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
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
	    $functionType = &DataTypes::dataObjectDocName($property)
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
	foreach my $propertyName ( &ExtraUtils::sortedKeys($component->{'properties'}->{'property'}) ) {
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
		(my $dataDefinition, my $label) = &DataTypes::dataObjectDefinition($property,matchOnly => 1);
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
		    $functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
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
		    $functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
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
			$functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
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
			    {type => "procedure", name => $propertyName."Rate", function => $functionLabel, description => "Cumulate to the rate of the {\\normalfont \\ttfamily ".$propertyName."} property of the {\\normalfont \\ttfamily ".$componentClassName."} component.", returnType => "\\void", arguments => &DataTypes::dataObjectDocName($property)."\\ value"}
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
	foreach my $propertyName ( &ExtraUtils::sortedKeys($component->{'properties'}->{'property'}) ) {
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
			(my $dataDefinition,my $label) = &DataTypes::dataObjectDefinition($linkedData);
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
			$functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
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
			    {type => "procedure", name => $propertyName.$suffix, function => $componentID.ucfirst($propertyName)."Get".$suffix, description => "Get the {\\normalfont \\ttfamily ".$propertyName."} property of the {\\normalfont \\ttfamily ".$componentClassName."} component.", returnType => &DataTypes::dataObjectDocName($property), arguments => ""}
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
			(my $dataDefinition,my $label) = &DataTypes::dataObjectDefinition($linkedData,matchOnly => 1);
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
			$functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
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
			    {type => "procedure", name => $propertyName."Set".$suffix, function => $componentID.ucfirst($propertyName)."Set".$suffix, description => "Set the {\\normalfont \\ttfamily ".$propertyName."} property of the {\\normalfont \\ttfamily ".$componentClassName."} component.", returnType => "\\void", arguments => &DataTypes::dataObjectDocName($property)."\\ value"}
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
		    $functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
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
		    (my $dataDefinition,my $label) = &DataTypes::dataObjectDefinition($linkedData,matchOnly => 1);
		    push(@{$dataDefinition->{'variables' }},"setValue"     );
		    push(@{$dataDefinition->{'attributes'}},"intent(in   )");
		    (my $currentDefinition,my $currentLabel) = &DataTypes::dataObjectDefinition($linkedData,matchOnly => 1);
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
		    $functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
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
			$typeDefinition{'arguments'  } = &DataTypes::dataObjectDocName($property)."\\ value";
		    }
		    push(
			@{$build->{'types'}->{'nodeComponent'.ucfirst($componentID)}->{'boundFunctions'}},
			\%typeDefinition
			);
		    if ( $property->{'attributes' }->{'makeGeneric'} eq "true" ) {
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
			$functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
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
			    $functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
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
		    $functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
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

sub Generate_Tree_Node_Creation_Function {
    # Generate a tree node creation function.
    my $build = shift;
    # Specify data content.
    my @dataContent =
	(
	 {
	     intrinsic  => "class",
	     type       => "treeNode",
	     attributes => [ "intent(inout)", "target" ],
	     variables  => [ "self" ]
	 },
	 {
	     intrinsic  => "integer",
	     type       => "kind=kind_int8",
	     attributes => [ "intent(in   )", "optional" ],
	     variables  => [ "index" ]
	 },
	 {
	     intrinsic  => "type",
	     type       => "mergerTree",
	     attributes => [ "intent(in   )", "optional", "target" ],
	     variables  => [ "hostTree" ]
	 }
	);
    # Create the function code.
    my $functionCode;
    $functionCode .= "  subroutine treeNodeInitialize(self,index,hostTree)\n";
    $functionCode .= "    !% Initialize a {\\normalfont \\ttfamily treeNode} object.\n";
    $functionCode .= "    use Galacticus_Error\n";
    $functionCode .= "    implicit none\n";
    $functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
    $functionCode .= "    ! Ensure pointers are nullified.\n";
    $functionCode .= "    nullify (self%".&Utils::padClass($_,[9,14]).")\n"
	foreach ( "parent", "firstChild", "sibling", "firstSatellite", "mergeTarget", "firstMergee", "siblingMergee", "formationNode", "event" );
    foreach ( @{$build->{'componentClassList'}} ) {
    	$functionCode .= "    allocate(self%".&Utils::padClass("component".ucfirst($_),[9,14])."(1))\n";
    }
    $functionCode .= "    select type (self)\n";
    $functionCode .= "    type is (treeNode)\n";
    foreach ( @{$build->{'componentClassList'}} ) {
	$functionCode .= "       self%component".&Utils::padClass(ucfirst($_),[0,0])."(1)%hostNode => self\n";
    }
    $functionCode .= "    end select\n";
    $functionCode .= "    ! Assign a host tree if supplied.\n";
    $functionCode .= "    if (present(hostTree)) self%hostTree => hostTree\n";
    $functionCode .= "    ! Assign index if supplied.\n";
    $functionCode .= "    if (present(index)) call self%indexSet(index)\n";
    $functionCode .= "    ! Assign a unique ID.\n";
    $functionCode .= "    !\$omp critical(UniqueID_Assign)\n";
    $functionCode .= "    uniqueIDCount=uniqueIDCount+1\n";
    $functionCode .= "    if (uniqueIDCount <= 0) call Galacticus_Error_Report('treeNodeInitialize','ran out of unique ID numbers')\n";
    $functionCode .= "    self%uniqueIdValue=uniqueIDCount\n";
    $functionCode .= "    !\$omp end critical(UniqueID_Assign)\n";
    $functionCode .= "    ! Assign a timestep.\n";
    $functionCode .= "    self%timeStepValue=-1.0d0\n";
    $functionCode .= "    return\n";
    $functionCode .= "  end subroutine treeNodeInitialize\n";	
    # Insert into the function list.
    push(
	@{$build->{'code'}->{'functions'}},
	$functionCode
	);
}

sub Generate_Tree_Node_Destruction_Function {
    # Generate a tree node destruction function.
    my $build = shift;
    # Specify data content.
    my @dataContent =
	(
	 {
	     intrinsic  => "class",
	     type       => "treeNode",
	     attributes => [ "intent(inout)" ],
	     variables  => [ "self" ]
	 },
	 {
	     intrinsic  => "class",
	     type       => "nodeEvent",
	     attributes => [ "pointer" ],
	     variables  => [ "thisEvent", "pairEvent", "lastEvent", "nextEvent" ]
	 },
	 {
	     intrinsic  => "logical",
	     variables  => [ "pairMatched" ]
	 }
	);
    # Create the function code.
    my $functionCode;
    $functionCode .= "  subroutine treeNodeDestroy(self)\n";
    $functionCode .= "    !% Destroy a {\\normalfont \\ttfamily treeNode} object.\n";
    $functionCode .= "    implicit none\n";
    $functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
    foreach my $componentClass ( @{$build->{'componentClassList'}} ) {
     	$functionCode .= "    call self%".&Utils::padClass(lc($componentClass)."Destroy",[7,0])."()\n";
    }
    # Remove any events attached to the node, along with their paired event in other nodes.
    $functionCode .= "    ! Iterate over all attached events.\n";
    $functionCode .= "    thisEvent => self%event\n";
    $functionCode .= "    do while (associated(thisEvent))\n";
    $functionCode .= "        ! Locate the paired event and remove it.\n";
    $functionCode .= "        pairEvent => thisEvent%node%event\n";
    $functionCode .= "        lastEvent => thisEvent%node%event\n";
    $functionCode .= "        ! Iterate over all events.\n";
    $functionCode .= "        pairMatched=.false.\n";
    $functionCode .= "        do while (associated(pairEvent).and..not.pairMatched)\n";
    $functionCode .= "           ! Match the paired event ID with the current event ID.\n";
    $functionCode .= "           if (pairEvent%ID == thisEvent%ID) then\n";
    $functionCode .= "              pairMatched=.true.\n";
    $functionCode .= "              if (associated(pairEvent,thisEvent%node%event)) then\n";
    $functionCode .= "                 thisEvent%node  %event => pairEvent%next\n";
    $functionCode .= "                 lastEvent       => thisEvent%node %event\n";
    $functionCode .= "              else\n";
    $functionCode .= "                 lastEvent%next  => pairEvent%next\n";
    $functionCode .= "              end if\n";
    $functionCode .= "              nextEvent => pairEvent%next\n";
    $functionCode .= "              deallocate(pairEvent)\n";
    $functionCode .= "              pairEvent => nextEvent\n";
    $functionCode .= "           else\n";
    $functionCode .= "              lastEvent => pairEvent\n";
    $functionCode .= "              pairEvent => pairEvent%next\n";
    $functionCode .= "           end if\n";
    $functionCode .= "        end do\n";
    $functionCode .= "        if (.not.pairMatched) call Galacticus_Error_Report('treeNodeDestroy','unable to find paired event')\n";
    $functionCode .= "        nextEvent => thisEvent%next\n";
    $functionCode .= "        deallocate(thisEvent)\n";
    $functionCode .= "        thisEvent => nextEvent\n";
    $functionCode .= "    end do\n";
    $functionCode .= "    return\n";
    $functionCode .= "  end subroutine treeNodeDestroy\n";	
    # Insert into the function list.
    push(
	@{$build->{'code'}->{'functions'}},
	$functionCode
	);
}

sub Generate_Tree_Node_Builder_Function {
    # Generate a tree node builder function.
    my $build = shift;
    # Specify data content.
    my @dataContent =
	(
	 {
	     intrinsic  => "class",
	     type       => "treeNode",
	     attributes => [ "intent(inout)" ],
	     variables  => [ "self" ]
	 },
	 {
	     intrinsic  => "type",
	     type       => "node",
	     attributes => [ "intent(in   )", "pointer" ],
	     variables  => [ "nodeDefinition" ]
	 },
	 {
	     intrinsic  => "type",
	     type       => "node",
	     attributes => [ "pointer" ],
	     variables  => [ "componentDefinition" ]
	 },
	 {
	     intrinsic  => "type",
	     type       => "nodeList",
	     attributes => [ "pointer" ],
	     variables  => [ "componentList" ]
	 },
	 {
	     intrinsic  => "integer",
	     variables  => [ "i", "j", "componentCount" ]
	 },
	 {
	     intrinsic  => "type",
	     type       => "integerScalarHash",
	     variables  => [ "componentIndex" ]
	 },
	 {
	     intrinsic  => "character",
	     type       => "len=128",
	     variables  => [ "nodeName" ]
	 }
	);
    # Create the function code.
    my $functionCode;
    $functionCode .= "  subroutine Tree_Node_Component_Builder(self,nodeDefinition)\n";
    $functionCode .= "    !% Build components in a {\\normalfont \\ttfamily treeNode} object given an XML definition.\n";
    $functionCode .= "    use FoX_Dom\n";
    $functionCode .= "    use Hashes\n";
    $functionCode .= "    implicit none\n";
    $functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
    $functionCode .= "    select type (self)\n";
    $functionCode .= "    type is (treeNode)\n";
    $functionCode .= "       call componentIndex%initialize()\n";
    $functionCode .= "       !\$omp critical (FoX_DOM_Access)\n";
   foreach my $componentClass ( @{$build->{'componentClassList'}} ) {
	$functionCode .= "    componentList => getChildNodes(nodeDefinition)\n";
	$functionCode .= "    componentCount=0\n";
	$functionCode .= "    do i=0,getLength(componentList)-1\n";
	$functionCode .= "      componentDefinition => item(componentList,i)\n";
	$functionCode .= "      if (getNodeName(componentDefinition) == '".$componentClass."') componentCount=componentCount+1\n";
	$functionCode .= "    end do\n";
	$functionCode .= "    if (componentCount > 0) then\n";
	$functionCode .= "      if (allocated(self%component".ucfirst($componentClass).")) deallocate(self%component".ucfirst($componentClass).")\n";
	$functionCode .= "      allocate(self%component".ucfirst($componentClass)."(componentCount),source=default".ucfirst($componentClass)."Component)\n";
	$functionCode .= "      call componentIndex%set('".$componentClass."',0)\n";
	$functionCode .= "    end if\n";
    }
    $functionCode .= "    componentCount=getLength(componentList)\n";
    $functionCode .= "    !\$omp end critical (FoX_DOM_Access)\n";
    $functionCode .= "    do i=0,componentCount-1\n";
    foreach my $componentClass ( @{$build->{'componentClassList'}} ) {
	$functionCode .= "     !\$omp critical (FoX_DOM_Access)\n";
	$functionCode .= "     componentDefinition => item(componentList,i)\n";
	$functionCode .= "     nodeName=getNodeName(componentDefinition)\n";
	$functionCode .= "     !\$omp end critical (FoX_DOM_Access)\n";
	$functionCode .= "     if (trim(nodeName) == '".$componentClass."') then\n";
	$functionCode .= "       j=componentIndex%value('".$componentClass."')\n";
	$functionCode .= "       j=j+1\n";
	$functionCode .= "       self%component".ucfirst($componentClass)."(j)%hostNode => self\n";
	$functionCode .= "       call self%component".ucfirst($componentClass)."(j)%builder(componentDefinition)\n";
	$functionCode .= "       call componentIndex%set('".$componentClass."',j)\n";
	$functionCode .= "     end if\n";
    }
    $functionCode .= "      end do\n";
    $functionCode .= "      call componentIndex%destroy()\n";
    $functionCode .= "    end select\n";
    $functionCode .= "    return\n";
    $functionCode .= "  end subroutine Tree_Node_Component_Builder\n";	
    # Insert into the function list.
    push(
	@{$build->{'code'}->{'functions'}},
	$functionCode
	);
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
	        foreach my $propertyName ( &ExtraUtils::sortedKeys($component->{'properties'}->{'property'}) ) {
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
	foreach my $propertyName ( &ExtraUtils::sortedKeys($properties) ) {
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
	    foreach my $componentName ( &ExtraUtils::sortedKeys($property) ) {
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

sub Generate_Type_Name_Functions {
    # Generate a type name functions.
    my $build = shift;
    # Initialize data content.
    my @dataContent;
    # Initialize the function code.
    my $functionCode;
    # Iterate over component classes.
    foreach ( @{$build->{'componentClassList'}} ) {
	# Specify data content.
	@dataContent =
	    (
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($_),
		 attributes => [ "intent(in   )" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "type",
		 type       => "varying_string",
		 variables  => [ "Node_Component_".ucfirst($_)."_Type" ]
	     }
	    );
	# Create the function code.
	$functionCode  = "  function Node_Component_".ucfirst($_)."_Type(self)\n";
	$functionCode .= "     !% Returns the type for the ".$_." component.\n";
	$functionCode .= "     implicit none\n";
	$functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
	$functionCode .= "     !GCC\$ attributes unused :: self\n";
	$functionCode .= "     ".&Utils::padClass("Node_Component_".ucfirst($_)."_Type",[20,0])."='nodeComponent:".$_."'\n";
	$functionCode .= "     return\n";
	$functionCode .= "  end function Node_Component_".ucfirst($_)."_Type\n\n";
	# Insert into the function list.
	push(
	    @{$build->{'code'}->{'functions'}},
	    $functionCode
	    );

	# Bind this function to the relevant type.
	push(
	    @{$build->{'types'}->{'nodeComponent'.ucfirst($_)}->{'boundFunctions'}},
	    {type => "procedure", name => "type", function => "Node_Component_".ucfirst($_)."_Type"}
	    );
    }
    # Iterate over implementations.
    foreach my $componentName ( @{$build->{'componentIdList'}} ) {
	my $component = $build->{'components'}->{$componentName};
	# Specify data content.
	@dataContent =
	    (
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($componentName),
		 attributes => [ "intent(in   )" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "type",
		 type       => "varying_string",
		 variables  => [ "Node_Component_".ucfirst($componentName)."_Type" ]
	     }
	    );
	# Create the function code.
  	$functionCode  = "  function Node_Component_".ucfirst($componentName)."_Type(self)\n";
	$functionCode .= "    !% Returns the type for the ".$component->{'name'}." implementation of the ".$component->{'class'}." component.\n";
	$functionCode .= "    implicit none\n";
	$functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
	$functionCode .= "    !GCC\$ attributes unused :: self\n";
	$functionCode .= "    ".&Utils::padImplementationPropertyName("Node_Component_".ucfirst($componentName)."_Type",[20,0])."='nodeComponent:".$component->{'class'}.":".$component->{'name'}."'\n";
	$functionCode .= "    return\n";
	$functionCode .= "  end function Node_Component_".ucfirst($componentName)."_Type\n\n";
	# Insert into the function list.
	push(
	    @{$build->{'code'}->{'functions'}},
	    $functionCode
	    );	
	# Bind this function to the relevant type.
	push(
	    @{$build->{'types'}->{'nodeComponent'.ucfirst($componentName)}->{'boundFunctions'}},
	    {type => "procedure", name => "type", function => "Node_Component_".ucfirst($componentName)."_Type"}
	    );
    }
}

sub Generate_Component_Assignment_Function {
    # Generate a type name functions.
    my $build = shift;
    # Specify data content.
    my @dataContent =
	(
	 {
	     intrinsic  => "class",
	     type       => "nodeComponent",
	     attributes => [ "intent(  out)" ],
	     variables  => [ "to" ]
	 },
	 {
	     intrinsic  => "class",
	     type       => "nodeComponent",
	     attributes => [ "intent(in   )" ],
	     variables  => [ "from" ]
	 }
	);
    # Generate the function code.
    my $functionCode;
    # Component assignment functions.
    $functionCode  = "  subroutine Node_Component_Assign(to,from)\n";
    $functionCode .= "    !% Assign a node component to another node component.\n";
    $functionCode .= "    implicit none\n";
    $functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
    $functionCode .= "    to%hostNode => from%hostNode\n";
    $functionCode .= "    select type (to)\n";
    foreach my $componentName ( @{$build->{'componentIdList'}} ) {
	my $component = $build->{'components'}->{$componentName};
	$functionCode .= "    type is (nodeComponent".&Utils::padFullyQualified($componentName,[0,0]).")\n";
	$functionCode .= "       select type (from)\n";
	$functionCode .= "       type is (nodeComponent".&Utils::padFullyQualified($componentName,[0,0]).")\n";
	foreach my $propertyName ( &ExtraUtils::sortedKeys($component->{'properties'}->{'property'}) ) {
	    my $property = $component->{'properties'}->{'property'}->{$propertyName};
	    if ( exists($property->{'linkedData'}) ) {
		my $linkedDataName = $property->{'linkedData'};		
		my $linkedData     = $component->{'content'}->{'data'}->{$linkedDataName};
		if ( $linkedData->{'type'} eq"double" ) {
		    # Deallocate if necessary.
		    $functionCode .= "   if (allocated(to%".&Utils::padLinkedData($linkedDataName,[0,0]).")) call Dealloc_Array(to%".&Utils::padLinkedData($linkedDataName,[0,0]).") \n"
			if ( $linkedData->{'rank'} > 0 );
		}
		elsif ( $linkedData->{'type'} eq"integer"     ) {
		    # Nothing to do in this case.
		}
		elsif ( $linkedData->{'type'} eq"longInteger" ) {
		    # Nothing to do in this case.
		}
		elsif ( $linkedData->{'type'} eq"logical"     ) {
		    # Nothing to do in this case.
		}
		else {
		    $functionCode .= "    call to%".&Utils::padLinkedData($linkedDataName,[0,0])."%destroy()\n";
		}
		$functionCode .= "          to%".&Utils::padLinkedData($linkedDataName,[0,0])."=from%".&Utils::padLinkedData($linkedDataName,[0,0])."\n";
	    }
	}
	$functionCode .= "       end select\n";
    }
    $functionCode .= "    end select\n";
    $functionCode .= "    return\n";
    $functionCode .= "  end subroutine Node_Component_Assign\n\n";
    # Insert into the function list.
    push(
	@{$build->{'code'}->{'functions'}},
	$functionCode
	);	
    # Bind this function to the nodeComponent type.
    push(
	@{$build->{'types'}->{'nodeComponent'}->{'boundFunctions'}},
	{type => "procedure", name => "assign"       , function => "Node_Component_Assign", description => "Assign a {\\normalfont \\ttfamily nodeComponent} to another {\\normalfont \\ttfamily nodeComponent}.", returnType => "\\textcolor{red}{\\textless class(nodeComponent)\\textgreater}", arguments => "\\textcolor{red}{\\textless class(nodeComponent)\\textgreater} from\\argin"},
	{type => "generic"  , name => "assignment(=)", function => "assign" }
	);
}

sub Generate_Component_Class_Destruction_Functions {
    # Generate class destruction functions.
    my $build = shift;
    # Initialize data content.
    my @dataContent;
    # Generate the function code.
    my $functionCode;
    foreach my $componentClassName ( @{$build->{'componentClassList'}} ) {
	# Specify data content.
	@dataContent =
	    (
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($componentClassName),
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     }
	    );
	# Generate the function code.
	$functionCode  = "  subroutine Node_Component_".ucfirst($componentClassName)."_Destroy(self)\n";
	$functionCode .= "    !% Destroys a ".$componentClassName." component.\n";
	$functionCode .= "    implicit none\n";
	$functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent);
	$functionCode .= "    !GCC\$ attributes unused :: self\n\n";
	$functionCode .= "    ! Do nothing.\n";
	$functionCode .= "    return\n";
	$functionCode .= "  end subroutine Node_Component_".ucfirst($componentClassName)."_Destroy\n\n";
	# Insert into the function list.
	push(
	    @{$build->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Bind this function to the treeNode type.
	push(
	    @{$build->{'types'}->{'nodeComponent'.ucfirst($componentClassName)}->{'boundFunctions'}},
	    {type => "procedure", name => "destroy", function => "Node_Component_".ucfirst($componentClassName)."_Destroy"}
	    );
    }
}

sub Generate_Component_Class_Removal_Functions {
    # Generate class removal functions.
    my $build = shift;

    # Initialize data content.
    my @dataContent;
    # Generate the function code.
    my $functionCode;
    foreach my $componentClassName ( @{$build->{'componentClassList'}} ) {
	# Specify data content.
	@dataContent =
	    (
	     {
		 intrinsic  => "class",
		 type       => "treeNode",
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "instance" ]
	     },
	     {
		 intrinsic  => "integer",
		 variables  => [ "instanceCount" ]
	     },
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($componentClassName),
		 attributes => [ "allocatable, dimension(:)" ],
		 variables  => [ "instancesTemporary" ]
	     }
	    );
	# Generate the function code.
	$functionCode  = "  subroutine Node_Component_".ucfirst($componentClassName)."_Remove(self,instance)\n";
	$functionCode .= "    !% Removes an instance of the ".$componentClassName." component, shifting other instances to keep the array contiguous.\n";
	$functionCode .= "    use Galacticus_Error\n";
	$functionCode .= "    implicit none\n";
	$functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
	$functionCode .= "    instanceCount=self%".$componentClassName."count()\n";
	$functionCode .= "    if (instance < 1 .or. instance > instanceCount) call Galacticus_Error_Report('Node_Component_".ucfirst($componentClassName)."_Remove','instance out of range')\n";
	$functionCode .= "    call self%component".ucfirst($componentClassName)."(instance)%destroy()\n";
	$functionCode .= "    if (instanceCount == 1) then\n";
	$functionCode .= "      ! Only one instance of this component. Deallocate it and reallocate with generic type.\n";
	$functionCode .= "      deallocate(self%component".ucfirst($componentClassName).")\n";
	$functionCode .= "      allocate(self%component".ucfirst($componentClassName)."(1))\n";
	$functionCode .= "    else\n";
	$functionCode .= "      ! Multiple instances, so remove the specified instance.\n";
	$functionCode .= "      allocate(instancesTemporary(instanceCount-1),source=self%component".ucfirst($componentClassName)."(1))\n";
	foreach my $implementationName ( @{$build->{'componentClasses'}->{$componentClassName}->{'memberNames'}} ) {
	    $functionCode .= "      select type (from => self%component".ucfirst($componentClassName).")\n";
	    $functionCode .= "      type is (nodeComponent".&Utils::padFullyQualified(ucfirst($componentClassName).ucfirst($implementationName),[0,0]).")\n";
	    $functionCode .= "        select type (to => instancesTemporary)\n";
	    $functionCode .= "        type is (nodeComponent".&Utils::padFullyQualified(ucfirst($componentClassName).ucfirst($implementationName),[0,0]).")\n";
	    $functionCode .= "          if (instance >             1) to(       1:instance     -1)=from(         1:instance     -1)\n";
	    $functionCode .= "          if (instance < instanceCount) to(instance:instanceCount-1)=from(instance+1:instanceCount  )\n";
	    $functionCode .= "        end select\n";
	    $functionCode .= "      end select\n";
	}
	$functionCode .= "      deallocate(self%component".ucfirst($componentClassName).")\n";
	$functionCode .= "      call Move_Alloc(instancesTemporary,self%component".ucfirst($componentClassName).")\n";
	
	$functionCode .= "    end if\n";
	$functionCode .= "    return\n";
	$functionCode .= "  end subroutine Node_Component_".ucfirst($componentClassName)."_Remove\n";
	# Insert into the function list.
	push(
	    @{$build->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Bind this function to the treeNode type.
	push(
	    @{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	    {type => "procedure", name => $componentClassName."Remove", function => "Node_Component_".ucfirst($componentClassName)."_Remove", description => "Remove an instance of the ".$componentClassName." component, shifting other instances to keep the array contiguous. If no {\\normalfont \\ttfamily instance} is specified, the first instance is assumed.", returnType => "\\void", arguments => "\\intzero\\ [instance]\\argin"}
	    );
    }
}

sub Generate_Component_Class_Move_Functions {
    # Generate class move functions.
    my $build = shift;
    # Initialize data content.
    my @dataContent;
    # Generate the function code.
    my $functionCode;
    foreach my $componentClassName ( @{$build->{'componentClassList'}} ) {
	# Specify data content.
	@dataContent =
	    (
	     {
		 intrinsic  => "class",
		 type       => "treeNode",
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "type",
		 type       => "treeNode",
		 attributes => [ "intent(inout)", "target" ],
		 variables  => [ "targetNode" ]
	     },
	     {
		 intrinsic  => "logical",
		 attributes => [ "intent(in   )", "optional" ],
		 variables  => [ "overwrite" ]
	     },
	     {
		 intrinsic  => "integer",
		 variables  => [ "instanceCount", "targetCount", "i" ]
	     },
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($componentClassName),
		 attributes => [ "allocatable, dimension(:)" ],
		 variables  => [ "instancesTemporary" ]
	     },
	     {
		 intrinsic  => "logical",
		 variables  => [ "overwrite_" ]
	     }
	    );
	# Generate the function code.
	$functionCode  = "  subroutine Node_Component_".ucfirst($componentClassName)."_Move(self,targetNode,overwrite)\n";
	$functionCode .= "    !% Move instances of the ".$componentClassName." component, from one node to another.\n";
	$functionCode .= "    use Galacticus_Error\n";
	$functionCode .= "    implicit none\n";
	$functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
	$functionCode .= "    overwrite_=.false.\n";
	$functionCode .= "    if (present(overwrite)) overwrite_=overwrite\n";
	$functionCode .= "    instanceCount=self      %".$componentClassName."count()\n";
	$functionCode .= "    targetCount  =targetNode%".$componentClassName."count()\n";
	$functionCode .= "    if (overwrite_ .and. targetCount > 0) then\n";
	$functionCode .= "      do i=1,targetCount\n";
	$functionCode .= "        call targetNode%component".ucfirst($componentClassName)."(i)%destroy()\n";
	$functionCode .= "      end do \n";
        $functionCode .= "      targetCount=0\n";
	$functionCode .= "      deallocate(targetNode%component".ucfirst($componentClassName).")\n";
	$functionCode .= "      allocate(targetNode%component".ucfirst($componentClassName)."(1))\n";
        $functionCode .= "    end if\n";	
	$functionCode .= "    if (instanceCount == 0) return\n";
	$functionCode .= "    if (targetCount == 0) then\n";
	$functionCode .= "      deallocate(targetNode%component".ucfirst($componentClassName).")\n";
	$functionCode .= "      call Move_Alloc(self%component".ucfirst($componentClassName).",targetNode%component".ucfirst($componentClassName).")\n";
	$functionCode .= "    else\n";
	$functionCode .= "      ! Multiple instances, so remove the specified instance.\n";
	$functionCode .= "      allocate(instancesTemporary(instanceCount+targetCount),source=self%component".ucfirst($componentClassName)."(1))\n";
	foreach my $implementationName ( @{$build->{'componentClasses'}->{$componentClassName}->{'memberNames'}} ) {
	    $functionCode .= "      select type (from => targetNode%component".ucfirst($componentClassName).")\n";
	    $functionCode .= "      type is (nodeComponent".&Utils::padFullyQualified(ucfirst($componentClassName).ucfirst($implementationName),[0,0]).")\n";
	    $functionCode .= "        select type (to => instancesTemporary)\n";
	    $functionCode .= "        type is (nodeComponent".&Utils::padFullyQualified(ucfirst($componentClassName).ucfirst($implementationName),[0,0]).")\n";
	    $functionCode .= "          to(1:targetCount)=from\n";
	    $functionCode .= "        end select\n";
	    $functionCode .= "      end select\n";
	}
	foreach my $implementationName ( @{$build->{'componentClasses'}->{$componentClassName}->{'memberNames'}} ) {
	    $functionCode .= "      select type (from => self%component".ucfirst($componentClassName).")\n";
	    $functionCode .= "      type is (nodeComponent".&Utils::padFullyQualified(ucfirst($componentClassName).ucfirst($implementationName),[0,0]).")\n";
	    $functionCode .= "        select type (to => instancesTemporary)\n";
	    $functionCode .= "        type is (nodeComponent".&Utils::padFullyQualified(ucfirst($componentClassName).ucfirst($implementationName),[0,0]).")\n";
	    $functionCode .= "          to(targetCount+1:targetCount+instanceCount)=from\n";
	    $functionCode .= "        end select\n";
	    $functionCode .= "      end select\n";
	}
	$functionCode .= "      call targetNode%".$componentClassName."Destroy()\n";
	$functionCode .= "      call self      %".$componentClassName."Destroy()\n";
	$functionCode .= "      call Move_Alloc(instancesTemporary,targetNode%component".ucfirst($componentClassName).")\n";
	$functionCode .= "      allocate(self%component".ucfirst($componentClassName)."(1))\n";
	$functionCode .= "    end if\n";
	$functionCode .= "    do i=1,size(targetNode%component".ucfirst($componentClassName).")\n";
	$functionCode .= "       targetNode%component".ucfirst($componentClassName)."(i)%hostNode => targetNode\n";
	$functionCode .= "    end do\n";
	$functionCode .= "    return\n";
	$functionCode .= "  end subroutine Node_Component_".ucfirst($componentClassName)."_Move\n";
	# Insert into the function list.
	push(
	    @{$build->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Bind this function to the treeNode type.
	push(
	    @{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	    {type => "procedure", name => $componentClassName."Move", function => "Node_Component_".ucfirst($componentClassName)."_Move", description => "", returnType => "\\void", arguments => "\\textcolor{red}{\\textless type(treeNode)\\textgreater} targetNode\\arginout"}
	    );
    }
}

sub Generate_Component_Class_Dump_Functions {
    # Generate dump for each component class.
    my $build = shift;
    # Iterate over component classes.
    foreach my $componentClassName ( @{$build->{'componentClassList'}} ) {
	# Initialize function code.
	my $functionCode;
	# Initialize data content.
	my @dataContent =
	    (
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($componentClassName),
		 attributes => [ "intent(in   )" ],
		 variables  => [ "self" ]
	     }
	    );
	# Generate dump function.
	$functionCode  = "  subroutine Node_Component_".ucfirst($componentClassName)."_Dump(self)\n";
	$functionCode .= "    !% Dump the contents of a generic ".$componentClassName." component.\n";
	$functionCode .= "    use Galacticus_Display\n";
	$functionCode .= "    use ISO_Varying_String\n";
	$functionCode .= "    implicit none\n";
	$functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
	$functionCode .= "    !GCC\$ attributes unused :: self\n";
	$functionCode .= "    call Galacticus_Display_Indent('".$componentClassName.": ".(" " x ($Utils::fullyQualifiedNameLengthMax-length($componentClassName)))."generic')\n";
	$functionCode .= "    call Galacticus_Display_Unindent('done')\n";
	$functionCode .= "    return\n";
	$functionCode .= "  end subroutine Node_Component_".ucfirst($componentClassName)."_Dump\n";
	# Insert into the function list.
	push(
	    @{$build->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Insert a type-binding for this function into the implementation type.
	push(
	    @{$build->{'types'}->{'nodeComponent'.ucfirst($componentClassName)}->{'boundFunctions'}},
	    {type => "procedure", name => "dump", function => "Node_Component_".ucfirst($componentClassName)."_Dump"},
	    );
    }
}

sub Generate_Component_Class_Initializor_Functions {
    # Generate initializor for each component class.
    my $build = shift;
    # Iterate over component classes.
    foreach my $componentClassName ( @{$build->{'componentClassList'}} ) {
	# Initialize function code.
	my $functionCode;
	# Initialize data content.
	my @dataContent =
	    (
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($componentClassName),
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     }
	    );
	# Generate initializor function.
	$functionCode  = "  subroutine Node_Component_".ucfirst($componentClassName)."_Initializor(self)\n";
	$functionCode .= "    !% Initialize a generic ".$componentClassName." component.\n";
	$functionCode .= "    use Galacticus_Error\n";
	$functionCode .= "    implicit none\n";
	$functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
	$functionCode .= "    !GCC\$ attributes unused :: self\n";
	$functionCode .= "    call Galacticus_Error_Report('Node_Component_".ucfirst($componentClassName)."_Initializor','can not initialize a generic component')\n";
	$functionCode .= "    return\n";
	$functionCode .= "  end subroutine Node_Component_".ucfirst($componentClassName)."_Initializor\n";
	# Insert into the function list.
	push(
	    @{$build->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Insert a type-binding for this function into the implementation type.
	push(
	    @{$build->{'types'}->{'nodeComponent'.ucfirst($componentClassName)}->{'boundFunctions'}},
	    {type => "procedure", name => "initialize", function => "Node_Component_".ucfirst($componentClassName)."_Initializor", description => "Initialize the object.", returnType => "\\void", arguments => ""},
	    );
    }
}

sub Generate_Component_Class_Builder_Functions {
    # Generate builder for each component class.
    my $build = shift;
    # Iterate over component classes.
    foreach my $componentClassName ( @{$build->{'componentClassList'}} ) {
	# Initialize function code.
	my $functionCode;
	# Initialize data content.
	my @dataContent =
	    (
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($componentClassName),
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "type",
		 type       => "node",
		 attributes => [ "intent(in   )", "pointer" ],
		 variables  => [ "componentDefinition" ]
	     }
	    );
	# Generate dump function.
	$functionCode  = "  subroutine Node_Component_".ucfirst($componentClassName)."_Builder(self,componentDefinition)\n";
	$functionCode .= "    !% Build a generic ".$componentClassName." component.\n";
	$functionCode .= "    use FoX_DOM\n";
	$functionCode .= "    use Galacticus_Error\n";
	$functionCode .= "    implicit none\n";
	$functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
	$functionCode .= "    !GCC\$ attributes unused :: self, componentDefinition\n";
	$functionCode .= "    call Galacticus_Error_Report('Node_Component_".ucfirst($componentClassName)."_Builder','can not build a generic component')\n";
	$functionCode .= "    return\n";
	$functionCode .= "  end subroutine Node_Component_".ucfirst($componentClassName)."_Builder\n";
	# Insert into the function list.
	push(
	    @{$build->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Insert a type-binding for this function into the implementation type.
	push(
	    @{$build->{'types'}->{'nodeComponent'.ucfirst($componentClassName)}->{'boundFunctions'}},
	    {type => "procedure", name => "builder", function => "Node_Component_".ucfirst($componentClassName)."_Builder", description => "Build a {\\normalfont \\ttfamily nodeComponent} from a supplied XML definition.", returnType => "\\void", arguments => "\\textcolor{red}{\\textless *type(node)\\textgreater}componentDefinition\\argin"},
	    );
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
	$functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
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
	$functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
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
	    foreach my $propertyName ( &ExtraUtils::sortedKeys($component->{'properties'}->{'property'}) ) {
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
			unless (&Utils::isOutputIntrinsic($type));
		    $rank1OutputTypes{$type} = 1
			if ( &Utils::isOutputIntrinsic($type) && $property->{'rank'} == 1 && exists($property->{'output'}->{'condition'}) );
		}
	    }
	}
	my %intrinsicMap = 
	    (
	     integer         => "integer(kind=kind_int8)",
	     longInteger     => "integer(kind=kind_int8)",
	     double => "double precision"
	    );
	foreach ( &ExtraUtils::sortedKeys(\%rank1OutputTypes) ) {
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
	foreach ( &ExtraUtils::sortedKeys(\%outputTypes) ){
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
	    foreach my $propertyName ( &ExtraUtils::sortedKeys($component->{'properties'}->{'property'}) ) {
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
	    foreach ( &ExtraUtils::sortedKeys(\%modulesRequired) );
	$functionCode .= "    implicit none\n";
	$functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
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
	    foreach my $propertyName ( &ExtraUtils::sortedKeys($component->{'properties'}->{'property'}) ) {
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
		    } elsif ( $property->{'isVirtual'} && $property->{'attributes'}->{'isGettable'} ) {
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
		    if (&Utils::isOutputIntrinsic($type)) {
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

sub Generate_Is_Active_Functions {
    # Generate "isActive" functions.
    my $build = shift;
    # Initialize function code.
    my $functionCode;
    # Iterate over component implementations.
    foreach my $componentID ( @{$build->{'componentIdList'}} ) {
	my $component = $build->{'components'}->{$componentID};
	$functionCode  = "  logical function Node_Component_".ucfirst($componentID)."_Is_Active()\n";
	$functionCode .= "    !% Return true if the ".$component->{'name'}." implementation of the ".$component->{'class'}." component is the active choice.\n";
	$functionCode .= "    implicit none\n\n";
	$functionCode .= "    Node_Component_".ucfirst($componentID)."_Is_Active=nodeComponent".ucfirst($componentID)."IsActive\n";
	$functionCode .= "    return\n";
	$functionCode .= "  end function Node_Component_".ucfirst($componentID)."_Is_Active\n";
	# Insert into the function list.
	push(
	    @{$build->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Bind this function to the implementation type.
	push(
	    @{$build->{'types'}->{'nodeComponent'.ucfirst($component->{'class'})}->{'boundFunctions'}},
	    {type => "procedure", pass => "nopass", name => lcfirst($component->{'name'})."IsActive", function => "Node_Component_".ucfirst($componentID)."_Is_Active", description => "Return whether the ".$component->{'name'}." implementation of the ".$component->{'class'}." component class is active.", returnType => "\\logicalzero", arguments => ""}
	    );
    }
}

sub Generate_Component_Implementation_Destruction_Functions {
    # Generate component implementation destruction functions.
    my $build = shift;
    # Initialize function code.
    my $functionCode;
    foreach my $componentID ( @{$build->{'componentIdList'}} ) {
	my $component = $build->{'components'}->{$componentID};
	# Specify data content.
	my @dataContent =
	    (
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($componentID),
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     }
	    );
	# Generate function code.
  	$functionCode  = "  subroutine Node_Component_".ucfirst($componentID)."_Destroy(self)\n";
	$functionCode .= "    !% Destroy a ".$component->{'name'}." implementation of the ".$component->{'class'}." component.\n";
	$functionCode .= "    use Memory_Management\n";
	$functionCode .= "    implicit none\n";
	$functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
	# Iterate over properties.
	my $selfUsed     = 0;
	my $functionBody = "";
	foreach my $propertyName ( &ExtraUtils::sortedKeys($component->{'properties'}->{'property'}) ) {
	    my $property = $component->{'properties'}->{'property'}->{$propertyName};
   	    # Check if this property has any linked data in this component.
	    if ( exists($property->{'linkedData'}) ) {
		# For each linked datum deallocate if necessary.
		my $linkedDataName = $property->{'linkedData'};
		my $linkedData     = $component->{'content'}->{'data'}->{$linkedDataName};
		if    ( $linkedData->{'type'} eq"double" ) {
		    # Nothing to do in this case.
		}
		elsif ( $linkedData->{'type'} eq"integer"     ) {
		    # Nothing to do in this case.
		}
		elsif ( $linkedData->{'type'} eq"longInteger" ) {
		    # Nothing to do in this case.
		}
		elsif ( $linkedData->{'type'} eq"logical"     ) {
		    # Nothing to do in this case.
		}
		else {
		    $functionBody .= "    call self%".&Utils::padLinkedData($linkedDataName,[0,0])."%destroy()\n";
		    $selfUsed = 1;
		}
		if ( $linkedData->{'rank'} > 0 ) {
		    $functionCode .= "    if (allocated(self%".&Utils::padLinkedData($linkedDataName,[0,0]).")) call Dealloc_Array(self%".&Utils::padLinkedData($linkedDataName,[0,0]).")\n";
		    $selfUsed = 1;
		}
	    }
	}
	$functionBody .= "    return\n";
	$functionBody .= "  end subroutine Node_Component_".ucfirst($componentID)."_Destroy\n\n";
	$functionCode .= "    !GCC\$ attributes unused :: self\n"
	    unless ( $selfUsed );
	$functionCode .= $functionBody;
	# Insert into the function list.
	push(
	    @{$build->{'code'}->{'functions'}},
	    $functionCode
	    );
	# Bind this function to the implementation type.
	push(
	    @{$build->{'types'}->{'nodeComponent'.ucfirst($componentID)}->{'boundFunctions'}},
	    {type => "procedure", name => "destroy", function => "Node_Component_".ucfirst($componentID)."_Destroy"}
	    );    
    }
}

sub Generate_Null_Binding_Functions {
    # Generate null binding functions.
    my $build = shift;
    # Iterate over component classes.
    foreach my $componentClassName ( &ExtraUtils::sortedKeys($build->{'nullProperties'}) ) {
	# Get the null functions required for this component class.
	my $componentClass = $build->{'nullProperties'}->{$componentClassName};
	# Iterate over required null functions for this component class.
	foreach my $nullFunctionName ( &ExtraUtils::sortedKeys($componentClass) ) {
	    # Get the null function definition.
	    my $nullFunction = $componentClass->{$nullFunctionName};
	    # Construct a datatype for this null function.
	    (my $dataDefinition, my $label) = &DataTypes::dataObjectDefinition($nullFunction,matchOnly => 1);
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
	    $functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
	    $functionCode .= "   !GCC\$ attributes unused :: ".join(", ",@{$_->{'variables'}})."\n"
		foreach ( @dataContent );
	    $functionCode .= "    return\n";
	    $functionCode .= "  end subroutine ".$nullBindingName."\n";
	    # Determine if this null rate function is actually used.
	    my $nullBindingUsed = 0;
	    foreach my $typeName ( keys(%{$build->{'types'}}) ) {
		foreach my $typeBoundFunction ( @{$build->{'types'}->{$typeName}->{'boundFunctions'}} ) {
		    if ( lc($typeBoundFunction->{'function'}) eq lc($nullBindingName) ) {
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
	    $functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
	    $functionCode .= "   !GCC\$ attributes unused :: ".join(", ",@{$_->{'variables'}})."\n"
		foreach ( @dataContent );
	    $functionCode .= "    return\n";
	    $functionCode .= "  end subroutine ".$nullBindingName."\n";
	    # Determine if this null rate function is actually used.
	    $nullBindingUsed = 0;
	    foreach my $typeName ( keys(%{$build->{'types'}}) ) {
		foreach my $typeBoundFunction ( @{$build->{'types'}->{$typeName}->{'boundFunctions'}} ) {
		    if ( lc($typeBoundFunction->{'function'}) eq lc($nullBindingName) ) {
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
	    $functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
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

sub Generate_Component_Class_Default_Value_Functions {
    # Generate component class default value functions.
    my $build = shift;
    # Iterate over component classes.
    foreach my $componentClassName ( @{$build->{'componentClassList'}} ) {
	# Initialize hash to track which property have been created already.
	my %propertiesCreated;
	# Iterate over implementations in this class.
    	foreach my $componentName ( @{$build->{'componentClasses'}->{$componentClassName}->{'memberNames'}} ) {
	    # Get the component.
	    my $componentID = ucfirst($componentClassName).ucfirst($componentName);
	    my $component   = $build->{'components'}->{$componentID};
	    # Iterate over the properties of this implementation.
	    foreach my $propertyName ( &ExtraUtils::sortedKeys($component->{'properties'}->{'property'}) ) {
		# Get the property.
		my $property = $component->{'properties'}->{'property'}->{$propertyName};
		# Get the linked data.
		my $linkedData;
		if ( exists($property->{'linkedData'}) ) {
		    my $linkedDataName = $property->{'linkedData'};
		    $linkedData = $component->{'content'}->{'data'}->{$linkedDataName};
		} else {
		    $linkedData = $property;
		}
		# Specify required data content.
		(my $dataDefinition, my $label ) = &DataTypes::dataObjectDefinition($linkedData);
		push(@{$dataDefinition->{'variables' }},ucfirst($componentClassName).ucfirst($propertyName));
		my @dataContent = (
		    $dataDefinition,
		    {
			intrinsic  => "class",
			type       => "nodeComponent".ucfirst($componentClassName),
			attributes => [ "intent(inout)" ],
			variables   => [ "self" ]
		    }
		    );
		# Skip if this property has already been created.
		unless ( exists($propertiesCreated{$propertyName}) ) {
		    # Generate code for "isGettable" function.
		    my $functionCode;
		    $functionCode  = "   logical function ".ucfirst($componentClassName).ucfirst($propertyName)."IsGettable()\n";
		    $functionCode .= "     !% Returns true if the {\\normalfont \\ttfamily ".$propertyName."} property is gettable for the {\\normalfont \\ttfamily ".$componentClassName."} component class.\n\n"; 
		    $functionCode .= "     implicit none\n";
		    $functionCode .= "     ".ucfirst($componentClassName).ucfirst($propertyName)."IsGettable=.false.\n";
		    foreach my $componentName2 ( @{$build->{'componentClasses'}->{$componentClassName}->{'memberNames'}} ) {
			my $component2ID = ucfirst($componentClassName).ucfirst($componentName2);
			my $component2   = $build->{'components'}->{$component2ID};
			$functionCode .= "     if (nodeComponent".ucfirst($component2ID)."IsActive) ".ucfirst($componentClassName).ucfirst($propertyName)."IsGettable=.true.\n"
			    if (
				exists($component2->{'properties'}->{'property'}->{$propertyName}                  ) && 
				exists($component2->{'properties'}->{'property'}->{$propertyName}->{'classDefault'})
			    );
		    }
		    $functionCode .= "     return\n";
		    $functionCode .= "   end function ".ucfirst($componentClassName).ucfirst($propertyName)."IsGettable\n";
		    # Insert into the function list.
		    push(
			@{$build->{'code'}->{'functions'}},
			$functionCode
			);
		    # Bind this function to the implementation type.
		    push(
			@{$build->{'types'}->{'nodeComponent'.ucfirst($componentClassName)}->{'boundFunctions'}},
			{type => "procedure", pass => "nopass", name => $propertyName."IsGettable", function => ucfirst($componentClassName).ucfirst($propertyName)."IsGettable", description => "Get the {\\normalfont \\ttfamily ".$propertyName."} property of the {\\normalfont \\ttfamily ".$componentClassName."} component.", returnType => &DataTypes::dataObjectDocName($property), arguments => ""}
			);
		    # Generate code for default value function.
		    $functionCode  = "  function ".ucfirst($componentClassName).ucfirst($propertyName)."(self)\n";
		    $functionCode .= "    !% Returns the default value for the {\\normalfont \\ttfamily ".$propertyName."} property for the {\\normalfont \\ttfamily ".$componentClassName."} component class.\n";
		    # Insert any required modules.
		    if ( exists($property->{'classDefault'}) && exists($property->{'classDefault'}->{'modules'}) ) {
			foreach ( @{$property->{'classDefault'}->{'modules'}} ) {
			    $functionCode .= "    use ".$_."\n";
			}
		    }
		    $functionCode .= "    implicit none\n";
		    # Build default value code, and accumulate which additional components are needed.
		    my $defaultLines = "";
		    my %requiredComponents;
		    foreach my $componentName2 ( @{$build->{'componentClasses'}->{$componentClassName}->{'memberNames'}} ) {
			my $component2ID = ucfirst($componentClassName).ucfirst($componentName2);
			my $component2   = $build->{'components'}->{$component2ID};
			if ( exists($component2->{'properties'}->{'property'}->{$propertyName}) ) {
			    my $property2 = $component2->{'properties'}->{'property'}->{$propertyName};
			    if ( exists($property2->{'classDefault'}) ) {
				$defaultLines .= "     if (nodeComponent".ucfirst($component2ID)."IsActive) then\n";
				my %selfComponents;
				my $default = $property2->{'classDefault'}->{'code'};
				while ( $default =~ m/self([a-zA-Z]+)Component\s*%/ ) {
				    $selfComponents{$1} = 1;
				    $requiredComponents{$1} = 1;
				    $default =~ s/self([a-zA-Z]+)Component\s*%//;
				}
				$defaultLines .= "    selfNode => self%host()\n" if ( scalar(keys(%selfComponents)) > 0 );
				foreach my $selfComponent ( &ExtraUtils::sortedKeys(\%selfComponents) ) {
				    $defaultLines .= "     self".$selfComponent."Component => selfNode%".lc($selfComponent)."()\n";
				}
				$defaultLines .= "       call Alloc_Array(".ucfirst($componentClassName).ucfirst($propertyName).",[".$property2->{'classDefault'}->{'count'}."])\n"
				    if ( exists($property2->{'classDefault'}->{'count'}) );
				$defaultLines .= "       ".ucfirst($componentClassName).ucfirst($propertyName)."=".$property2->{'classDefault'}->{'code'}."\n";
				$defaultLines .= "       return\n";
				$defaultLines .= "     end if\n";
			    }
			}
		    }
		    # Add a self node pointer if other components are required.
		    push(
			@dataContent,
			{
			    intrinsic  => "type",
			    type       => "treeNode",
			    attributes => [ "pointer" ],
			    variables  => [ "selfNode" ]
			}
			) if ( scalar(keys(%requiredComponents)) > 0 );
		    # Add pointers for each required component.
		    push(
			@dataContent,
			{
			    intrinsic  => "class",
			    type       => "nodeComponent".ucfirst($_),
			    attributes => [ "pointer" ],
			    variables  => [ "self".ucfirst($_)."Component" ]
			}
			)
			foreach ( &ExtraUtils::sortedKeys(\%requiredComponents) );
		    # Insert data content.
		    $functionCode .= &Fortran_Utils::Format_Variable_Defintions(\@dataContent)."\n";
		    $functionCode .= "   !GCC\$ attributes unused :: self\n";
		    # Insert code to set required default.
		    $functionCode .= $defaultLines;
		    # Insert code to return zero values by default.
		    if ( $linkedData->{'rank'} == 0 ) {
			if    ( $linkedData->{'type'} eq"double" ) {
			    $functionCode .= "    ".ucfirst($componentClassName).ucfirst($propertyName)."=0.0d0\n";
			}
			elsif ( $linkedData->{'type'} eq"integer"     ) {
			    $functionCode .= "    ".ucfirst($componentClassName).ucfirst($propertyName)."=0\n";
			}
			elsif ( $linkedData->{'type'} eq"longInteger" ) {
			    $functionCode .= "    ".ucfirst($componentClassName).ucfirst($propertyName)."=0_kind_int8\n";
			}
			elsif ( $linkedData->{'type'} eq"logical"     ) {
			    $functionCode .= "    ".ucfirst($componentClassName).ucfirst($propertyName)."=.false.\n";
			}
			else {
			    $functionCode .= "     call ".ucfirst($componentClassName).ucfirst($propertyName)."%reset()\n";
			}
		    } else {
			$functionCode .= "    ".ucfirst($componentClassName).ucfirst($propertyName)."=null".$label.$linkedData->{'rank'}."d\n";
		    }
		    # Close the function.
		    $functionCode .= "    return\n";
		    $functionCode .= "  end function ".ucfirst($componentClassName).ucfirst($propertyName)."\n";
		    # Insert into the function list.
		    push(
			@{$build->{'code'}->{'functions'}},
			$functionCode
			);
		    # Bind this function to the implementation type.
		    push(
			@{$build->{'types'}->{'nodeComponent'.ucfirst($componentClassName)}->{'boundFunctions'}},
			{type => "procedure", name => $propertyName, function => ucfirst($componentClassName).ucfirst($propertyName), description => "Get the {\\normalfont \\ttfamily ".$propertyName."} property of the {\\normalfont \\ttfamily ".$componentClassName."} component.", returnType => &DataTypes::dataObjectDocName($property), arguments => ""}
			);
		    # Record that this property has been created.
		    $propertiesCreated{$propertyName} = 1;
		}
	    }
	}
    }
}

sub Bound_Function_Table {
    # Get the list of type-bound functions.
    my $objectName         = shift;
    my @typeBoundFunctions = sort {$a->{'function'} cmp $b->{'function'}} @{$_[0]};
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
	} else {
	    # Single function specified. Simply add to the table.
	    $table->add($_->{'type'},$pass," :: ",$name,$connector,$_->{'function'},"");
	}
    }
    # Add any descriptions.
    my $description;
    my $methodCount = 0;
    foreach ( @typeBoundFunctions ) {
	if ( exists($_->{'description'}) ) {
	    ++$methodCount;
	    $description .= "     !@  <objectMethod>\n     !@   <method>".$_->{'name'}."</method>\n     !@   <description>".$_->{'description'}."</description>\n";
	    $description .= "     !@    <type>".$_->{'returnType'}."</type>\n"
		if ( exists($_->{'returnType'}) );
	    $description .= "     !@    <arguments>".$_->{'arguments'}."</arguments>\n"
		if ( exists($_->{'arguments'}) );
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

sub Insert_Type_Definitions {
    # Generate and insert code for all type definitions.
    my $build = shift;
    # Sort types into dependency order.
    my %typeDependencies;
    foreach ( &ExtraUtils::hashList($build->{'types'}) ) {
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
    my @typeSort  = &ExtraUtils::sortedKeys($build->{'types'});
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
	$build->{'content'} .= &Fortran_Utils::Format_Variable_Defintions($type->{'dataContent'})
	    if ( exists($type->{'dataContent'}) );
	# Generate and insert a type-bound function table.
	if ( exists($type->{'boundFunctions'}) ) {
	    $build->{'content'} .= "   contains\n";
	    my $boundFunctionTable = &Bound_Function_Table($type->{'name'},$type->{'boundFunctions'});   
	    $build->{'content'} .= $boundFunctionTable;
	}
	# Insert the type closing.
	$build->{'content'} .= "  end type ".$type->{'name'}."\n\n";
    }
}

sub Generate_Interfaces {
    # Generate and insert code for all interfaces.
    my $build = shift();
    # Iterate over interfaces.
    print "   --> Interfaces...\n";
    foreach ( &ExtraUtils::hashList($build->{'interfaces'}) ) {
	print "      ---> ".$_->{'name'}."\n";
	$CodeGeneration::interface = $_;
	$build->{'content'} .= 
	    fill_in_string(<<'CODE', PACKAGE => 'CodeGeneration')."\n";
! {$interface->{'comment'}}
abstract interface
  {$interface->{'intrinsic'} eq "void" ? "subroutine" : $interface->{'intrinsic'}." function"} {$interface->{'name'}}({join(",",&Function_Arguments($interface->{'data'}))})
    {&Importables($interface->{'data'}) ? "import ".join(", ",&Importables($interface->{'data'})) : ""}
{&Fortran_Utils::Format_Variable_Defintions($interface->{'data'}, indent => 4)}
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
    foreach my $function ( @{$build->{'functions'}} ) {
	# Report.
	print "      --> ".$function->{'name'}."\n";
	# Build function type definition.
	$function->{'type'} = $function->{'type'} eq "void" ? "subroutine" : $function->{'type'}." function";
	# Build a list of arguments.
	my @arguments;
	foreach my $variables ( @{$function->{'variables'}} ) {
	    push(
		@arguments,
		@{$variables->{'variables'}}
		)
		if ( exists($variables->{'attributes'}) && grep {$_ =~ m/^intent\s*\(\s*(in|inout|out)\s*\)/} @{$variables->{'attributes'}} );
	}
	# Serialize function opener.
	$build->{'content'} .= $function->{'type'}." ".$function->{'name'}."(".join(",",@arguments).")\n";
	# Serialize description.
	$build->{'content'} .= "   !% ".$function->{'description'}."\n";
	# Serialize module uses.
	$build->{'content'} .= "   use ".$_."\n"
	    foreach ( @{$function->{'modules'}} );
	# Serialize variable definitions.
	$build->{'content'} .= "   implicit none\n";
	$build->{'content'} .= &Fortran_Utils::Format_Variable_Defintions($function->{'variables'})
	    if ( exists($function->{'variables'}) );
	# Serialize content.
	$build->{'content'} .= $function->{'content'};
	# Serialize function closer.
	$build->{'content'} .= "   return\n";
	$build->{'content'} .= "end ".$function->{'type'}." ".$function->{'name'}."\n\n";
    }
}

1;
