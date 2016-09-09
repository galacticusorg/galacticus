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
use Galacticus::Build::Components::NullFunctions qw(createNullFunction);
use Galacticus::Build::Components::CodeGeneration;
use Galacticus::Build::Components::Hierarchy;
use Galacticus::Build::Components::Hierarchy::ODESolver;
use Galacticus::Build::Components::Hierarchy::State;
use Galacticus::Build::Components::Hierarchy::Utils;
use Galacticus::Build::Components::TreeNodes;
use Galacticus::Build::Components::TreeNodes::CreateDestroy;
use Galacticus::Build::Components::TreeNodes::Classes;
use Galacticus::Build::Components::TreeNodes::ODESolver;
use Galacticus::Build::Components::TreeNodes::Output;
use Galacticus::Build::Components::TreeNodes::Map;
use Galacticus::Build::Components::TreeNodes::Serialization;
use Galacticus::Build::Components::TreeNodes::Utils;
use Galacticus::Build::Components::NodeEvents;
use Galacticus::Build::Components::BaseTypes;
use Galacticus::Build::Components::Classes;
use Galacticus::Build::Components::Classes::Names;
use Galacticus::Build::Components::Classes::CreateDestroy;
use Galacticus::Build::Components::Classes::Deferred;
use Galacticus::Build::Components::Classes::Defaults;
use Galacticus::Build::Components::Classes::Output;
use Galacticus::Build::Components::Classes::State;
use Galacticus::Build::Components::Classes::Serialization;
use Galacticus::Build::Components::Classes::Utils;
use Galacticus::Build::Components::Implementations;
use Galacticus::Build::Components::Implementations::Utils;
use Galacticus::Build::Components::Implementations::Names;
use Galacticus::Build::Components::Implementations::Deferred;
use Galacticus::Build::Components::Implementations::CreateDestroy;
use Galacticus::Build::Components::Implementations::State;
use Galacticus::Build::Components::Implementations::Serialization;
use Galacticus::Build::Components::Implementations::Output;
use Galacticus::Build::Components::Implementations::ODESolver;
use Galacticus::Build::Components::Properties;
use Galacticus::Build::Components::Properties::Attributes;
use Galacticus::Build::Components::Properties::Deferred;
use Galacticus::Build::Components::Properties::Get;
use Galacticus::Build::Components::Properties::Set;
use Galacticus::Build::Components::Properties::Evolve;
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
    my @hooks =
	map
        {{name => $_, hook => $Galacticus::Build::Component::Utils::componentUtils{$_}}}
        &List::ExtraUtils::sortedKeys(\%Galacticus::Build::Component::Utils::componentUtils);
    # Iterate over phases.
    print "--> Phase:\n";
    foreach my $phase ( "preValidate", "default", "gather", "scatter", "postValidate", "content", "types", "interfaces", "functions" ) {
	print "   --> ".ucfirst($phase)."...\n";
	foreach my $hook ( grep {exists($_->{'hook'}->{$phase})} @hooks ) {	
	    my @functions = &List::ExtraUtils::as_array($hook->{'hook'}->{$phase});
	    foreach my $function ( @functions ) {
		print "      --> ".$hook->{'name'}.(scalar(@functions) > 1 ? " {".sub_name($function)."}" : "")."\n";
		&{$function}($build);
	    }
	}
    }
    # Iterate over all functions, calling them with the build data object.
    &{$_}($build)
	foreach (
	    # Generate functions for getting/setting/rating value via a deferred function.
	    \&Generate_Deferred_GSR_Function
	);

    # Insert all derived-type variable definitions.
    &derivedTypesSerialize($build);
    # Insert all interfaces.
    &interfacesSerialize  ($build);
    # Insert all module scope variables.
    $build->{'content'} .= &Fortran::Utils::Format_Variable_Definitions($build->{'variables'})."\n";
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
	sleep(2);
    } else {
	print "--> No legacy functions remain.....functionality should be removed\n";
	sleep(2);
    }
    $build->{'content'} .= join("\n",@{$build->{'code'}->{'functions'}})."\n";
    
    # Insert include statements to bring in all functions associated with components.
    my @includeDependencies  = map {exists($_->{'functions'}) ? $_->{'functions'} : ()} &List::ExtraUtils::hashList($build->{'components'});
    $build->{'content'}     .= join("\n",map {"  include \"".$_."\"\n"} @includeDependencies)."\n";
    # Create a Makefile to specify dependencies on these include files.
    open(makeFile,">".$ENV{'BUILDPATH'}."/Makefile_Component_Includes.tmp");
    print makeFile $ENV{'BUILDPATH'}."/objects.nodes.o:".join("",map {" ".$ENV{'BUILDPATH'}."/".$_} @includeDependencies)
	if ( scalar(@includeDependencies) > 0 );
    close(makeFile);
    &File::Changes::Update($ENV{'BUILDPATH'}."/Makefile_Component_Includes" ,$ENV{'BUILDPATH'}."/Makefile_Component_Includes.tmp" );
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
    # Get the component fully-qualified, class and method names.
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
	my $type = $gsr eq "get" ? $componentName.ucfirst($propertyName).ucfirst($gsr) : &createNullFunction($build,{selfType => $selfType, attribute => $gsr, property => $property, intent => "inout"});
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
	$functionCode .= &Fortran::Utils::Format_Variable_Definitions(\@dataContent)."\n";
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
		    $property->{'getFunction'}->{'build'     }
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
		    $functionCode .= &Fortran::Utils::Format_Variable_Definitions(\@dataContent)."\n";
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
		    && $property->{'setFunction'}->{'build'     }
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
		    $functionCode .= &Fortran::Utils::Format_Variable_Definitions(\@dataContent)."\n";
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
			$functionCode .= &Fortran::Utils::Format_Variable_Definitions(\@dataContent)."\n";
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
	(
	 {
	     align  => "left"
	 }
	) x 7
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
	$build->{'content'} .= &Fortran::Utils::Format_Variable_Definitions($type->{'dataContent'})
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
{&Fortran::Utils::Format_Variable_Definitions($interface->{'data'}, indent => 4)}
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
	my @arguments = &Galacticus::Build::Components::Utils::argumentList(@{$function->{'variables'}});
	# Serialize function opener.
	$build->{'content'} .= join(" ",@functionAttributes)." ".$type." ".$form." ".$function->{'name'}."(".join(",",@arguments).") ".$result."\n";
	# Serialize description.
	$build->{'content'} .= "   !% ".$function->{'description'}."\n";
	# Serialize module uses.
	$build->{'content'} .= "   use ".$_."\n"
	    foreach ( @{$function->{'modules'}} );
	# Serialize variable definitions.
	$build->{'content'} .= "   implicit none\n";
	$build->{'content'} .= &Fortran::Utils::Format_Variable_Definitions($function->{'variables'})
	    if ( exists($function->{'variables'}) );
	# Serialize content.
	$build->{'content'} .= $function->{'content'};
	# Serialize function closer.
	$build->{'content'} .= "   return\n";
	$build->{'content'} .= "end ".$form." ".$function->{'name'}."\n\n";
    }
}

1;
