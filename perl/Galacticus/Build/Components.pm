# Contains a Perl module which implements processing of "component" directives in the Galacticus build system.

package Galacticus::Build::Components;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Data::Dumper;
use Text::Table;
use Text::Template 'fill_in_string';
use Sort::Topo;
use Scalar::Util 'reftype';
use XML::SAX::ParserFactory;
use XML::Validator::Schema;
use LaTeX::Encode;
use Carp 'verbose';
use Sub::Identify ':all';
use File::Changes;
use Fortran::Utils;
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
use Galacticus::Build::Components::TreeNodes::State;
use Galacticus::Build::Components::TreeNodes::Utils;
use Galacticus::Build::Components::NodeEvents;
use Galacticus::Build::Components::Components;
use Galacticus::Build::Components::Classes;
use Galacticus::Build::Components::Classes::Names;
use Galacticus::Build::Components::Classes::CreateDestroy;
use Galacticus::Build::Components::Classes::Evolve;
use Galacticus::Build::Components::Classes::Deferred;
use Galacticus::Build::Components::Classes::Defaults;
use Galacticus::Build::Components::Classes::MetaProperties;
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
    my $validator = XML::Validator::Schema->new(file => $ENV{'GALACTICUS_EXEC_PATH'}."/schema/componentSchema.xsd");
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
    # Check for pre-existing component with identical name.
    die("Galacticus::Build::Components::Components_Parse_Directive: multiple components with ID '".$componentID."'")
	if ( exists($build->{'components'}->{$componentID}) );
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
    # Insert include statements to bring in all functions associated with components.
    my @includeDependencies  = map {exists($_->{'functions'}) ? $_->{'functions'} : ()} &List::ExtraUtils::hashList($build->{'components'});
    $build->{'content'}     .= join("\n",map {"  include \"".$_."\"\n"} @includeDependencies)."\n";
    # Create a Makefile to specify dependencies on these include files.
    open(makeFile,">".$ENV{'BUILDPATH'}."/Makefile_Component_Includes.tmp");
    print makeFile $ENV{'BUILDPATH'}."/objects.nodes.o ".$ENV{'BUILDPATH'}."/objects.nodes.p.F90:".join("",map {" ".$ENV{'BUILDPATH'}."/".$_} @includeDependencies)
	if ( scalar(@includeDependencies) > 0 );
    close(makeFile);
    &File::Changes::Update($ENV{'BUILDPATH'}."/Makefile_Component_Includes" ,$ENV{'BUILDPATH'}."/Makefile_Component_Includes.tmp" ); 
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
	    my $methodName = (exists($_->{'descriptor'}) && exists($_->{'descriptor'}->{'methodName'})) ? $_->{'descriptor'}->{'methodName'} : $_->{'name'};
	    $description .= "      <method method=\"".$methodName."\" description=\"".$descriptionText."\"/>\n";
	}
    }
    if ( $methodCount >= 1 ) {
	$description = "     !![\n     <methods>\n".$description."     </methods>\n     !!]\n\n";
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
    my @typeOrder = &Sort::Topo::sort(\@typeSort,\%typeDependencies);
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
	if ( exists($type->{'comment'}) ) {
	    $build->{'content'} .= "  !!{\n";
	    $build->{'content'} .= "  ".$type->{'comment'}."\n";
	    $build->{'content'} .= "  !!}\n";
	}
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
  {$interface->{'intrinsic'} eq "void" ? "subroutine" : (exists($interface->{'result'}) && $interface->{'result'} eq "yes" ? "" : $interface->{'intrinsic'})." function"} {$interface->{'name'}}({join(",",&Function_Arguments($interface->{'data'}))})
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
	    if ( $returnDescriptor =~ m/([a-zA-Z0-9_=\(\)\s]+)\s*,\s*([a-zA-Z0-9_,\(\):=\s]+)/ ) {
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
	$build->{'content'} .= "   !!{\n";
	$build->{'content'} .= "   ".$function->{'description'}."\n";
	$build->{'content'} .= "   !!}\n";
	# Serialize module uses.
	my @intrinsicModules = ( "iso_c_binding" );
	foreach my $module ( @{$function->{'modules'}} ) {
	    my $intrinsic = (grep {lc($module) eq $_} @intrinsicModules) ? ", intrinsic" : "";
	    $build->{'content'} .= "   use".$intrinsic." :: ".$module."\n";
	}
	# Serialize variable definitions.
	$build->{'content'} .= "   implicit none\n";
	$build->{'content'} .= &Fortran::Utils::Format_Variable_Definitions($function->{'variables'})
	    if ( exists($function->{'variables'}) );
	# Serialize content.
	$build->{'content'} .= $function->{'content'}
	    if ( exists($function->{'content'  }) );
	# Serialize function closer.
	$build->{'content'} .= "   return\n";
	$build->{'content'} .= "end ".$form." ".$function->{'name'}."\n\n";
    }
}

1;
