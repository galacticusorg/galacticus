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
use Galacticus::Build::Components::Implementations::ODESolver;
use Galacticus::Build::Components::Properties;
use Galacticus::Build::Components::Properties::Attributes;
use Galacticus::Build::Components::Properties::Deferred;
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
	    # Generate output functions for each implementation.
	    \&Generate_Implementation_Output_Functions               ,
	    # Generate functions for getting/setting/rating value via a deferred function.
	    \&Generate_Deferred_GSR_Function                         ,
	    # Generate functions for getting/setting/rating value directly.
	    \&Generate_GSR_Functions
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
    my @includeDependencies  = map {exists($_->{'functions'}) ? $_->{'functions'} : ()} &List::ExtraUtils::hashList($build->{'components'});
    $build->{'content'}     .= join("\n",map {"  include \"".$_."\"\n"} @includeDependencies)."\n";
    # Create a Makefile to specify dependencies on these include files.
    open(makeFile,">".$ENV{'BUILDPATH'}."/Makefile_Component_Includes.tmp");
    print makeFile $ENV{'BUILDPATH'}."/objects.nodes.o:".join("",map {" ".$ENV{'BUILDPATH'}."/".$_} @includeDependencies)
	if ( scalar(@includeDependencies) > 0 );
    close(makeFile);
    &File::Changes::Update($ENV{'BUILDPATH'}."/Makefile_Component_Includes" ,$ENV{'BUILDPATH'}."/Makefile_Component_Includes.tmp" );
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
