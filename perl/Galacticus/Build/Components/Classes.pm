# Contains a Perl module which handles component classes during build.

package Galacticus::Build::Components::Classes;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Data::Dumper;
use List::ExtraUtils;
use Galacticus::Build::Components::Utils qw($verbosityLevel);
use Galacticus::Build::Components::NullFunctions qw(createNullFunction);
use Galacticus::Build::Components::DataTypes;

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     classes => {
	 gather =>
	     [
	      \&Gather_Classes
	     ],
	 types =>
	     [
	      \&Build_Component_Classes
	     ]
     }
    );

sub Gather_Classes {
    # Determine class membership of each component implementation, and build a list of component class names.
    my $build = shift();
    # Iterate over all component instances.
    foreach my $implementation ( &List::ExtraUtils::hashList($build->{'components'}) ) {
	# Get the name of the component class.
	my $className          = $implementation->{'class'};
	# Get the name of the implementation.
	my $implementationName = $implementation->{'name' };
	# Append the component name to the list of members for its class.
	push(
	    @{$build->{'componentClasses'}->{$className}->{'memberNames'}},
	    $implementationName
	    );
	push(
	    @{$build->{'componentClasses'}->{$className}->{'members'    }},
	    $implementation
	    );
    }
    # Construct a list of component class names.
    @{$build->{'componentClassList'}} = &List::ExtraUtils::sortedKeys($build->{'componentClasses'});
    # Build a list of active classes. By default this is all available classes, but if the GALACTICUS_ACTIVE_COMPONENTS
    # environment variable is set, then use it as a list of active components.
    if ( exists($ENV{'GALACTICUS_ACTIVE_COMPONENTS'}) ) {
	my @activeClasses = split(" ",$ENV{'GALACTICUS_ACTIVE_COMPONENTS'});
	foreach my $className ( @{$build->{'componentClassList'}} ) {
	    push(@{$build->{'componentClassListActive'}},$className)
		if ( grep {$className eq $_} @activeClasses );
	}
    } else {
	@{$build->{'componentClassListActive'}} = @{$build->{'componentClassList'}};
    }    
    # Output report if sufficiently verbose.
    if ( $verbosityLevel > 0 ) {
	print "         --> Found the following component classes and implementations:\n";
	foreach my $className ( @{$build->{'componentClassList'}} ) {
	    print "            --> ".$className."\n";
	    foreach my $implementationName ( @{$build->{'componentClasses'}->{$className}->{'memberNames'}} ) {
		print "               --> ".$implementationName."\n";
	    }
	}
    }
}

sub Build_Component_Classes {
    # Generate type definitions for each component class.
    my $build = shift();    
    # Iterate over all component classes.
    my %classGetDefaults;
    foreach my $class ( &List::ExtraUtils::hashList($build->{'componentClasses'}, keyAs => 'name') ) {
	my $className  = $class->{'name'};
	# Define a hash to record which properties have already been created.
	my %propertiesCreated;
	# Create a list for type-bound functions.
	my @typeBoundFunctions;
  	# Insert definitions for each method associated with a component implementation of this component class.
    	foreach my $implementationName ( @{$class->{'memberNames'}} ) {
	    # Construct a fully-qualified name for this implementation.
	    my $implementationIdentifier = ucfirst($className).ucfirst($implementationName);
	    # Get the implementation.
	    my $implementation           = $build->{'components'}->{$implementationIdentifier};
	    # Iterate over properties beloning to this implementation.
	    foreach my $property ( &List::ExtraUtils::hashList($implementation->{'properties'}->{'property'}, keyAs => 'name') ) {
		# Create functions to set/get/evolve each property as necessary.
		if ( 
		    $property->{'attributes'}->{'isGettable' }
		    || 
		    $property->{'attributes'}->{'isSettable' }
		    || 
		    $property->{'attributes'}->{'isEvolvable'}
		    )
		{
		    # Name of function being processed.
		    my $functionName;
		    # Get a fully-qualified type identfier for this property.
		    (my $intrinsic,my $type,my $attributes) = &Galacticus::Build::Components::DataTypes::dataObjectPrimitiveName($property);
		    $type .= $property->{'rank'}."InOut";
		    # Handle set functions and related functions.
		    if ( $property->{'attributes'}->{'isSettable'} ) {
			# Create a "set" function if one does not already exist.
			$functionName = $property->{'name'}."Set";
			unless ( exists($propertiesCreated{$functionName}) ) {
				# Create a binding to a null function here. 
			    my $boundTo = &createNullFunction($build,{selfType => $className, attribute => "set", property => $property, intent => "inout"});
			    push(
				@typeBoundFunctions,
				{
				    type        => "procedure"  , 
				    name        => $functionName, 
				    function    => $boundTo     ,
				    returnType  => "\\void"     ,
				    arguments   => &Galacticus::Build::Components::DataTypes::dataObjectDocName($property)."\\ value",
				    description => "Set the {\\normalfont \\ttfamily ".$property->{'name'}."} property of the {\\normalfont \\ttfamily ".$className."} component."
				}
				);
			    $propertiesCreated{$functionName} = 1;
			}
		    }
		    # Handle evolve functions.
		    if ( $property->{'attributes'}->{'isEvolvable'} ) {
			# Create the "count" function.
			$functionName = $property->{'name'}."Count";
			unless ( exists($propertiesCreated{$functionName}) ) {
			    push(
				@typeBoundFunctions,
				{
				    type        => "procedure"                       ,
				    name        => $functionName                     ,
				    function    => &createNullFunction($build,{selfType => $className, attribute => "get", property => {type => "integer", rank => 0}, intent => "in"}),
				    returnType  => "\\intzero"                       ,
				    arguments   => ""                                ,
				    description => "Compute the count of evolvable quantities in the {\\normalfont \\ttfamily ".$property->{'name'}."} property of the {\\normalfont \\ttfamily ".$implementationIdentifier."} component."
				}
				);
			    $propertiesCreated{$functionName} = 1;
			}
			# Create the "rate" function.
			$functionName = $property->{'name'}."Rate";
			unless (  exists($propertiesCreated{$functionName}) ) {
			    # Do not create a rate function if it is to be deferred, or if it binds at the top level.
			    push(
				@typeBoundFunctions,
				{
				    type        => "procedure"  ,
				    name        => $functionName,
				    function    => &createNullFunction($build,{selfType => $className, attribute => "rate", property => $property, intent => "inout"}),
				    returnType  => "\\void"     ,
				    arguments   => &Galacticus::Build::Components::DataTypes::dataObjectDocName($property)."\\ value",
				    description => "Cumulate to the rate of the {\\normalfont \\ttfamily ".$property->{'name'}."} property of the {\\normalfont \\ttfamily ".$implementationIdentifier."} component."
				}
				)
				unless ( $property->{'attributes'}->{'createIfNeeded'} );
			    # Create an analytic function unless this is a virtual property.
			    push(
				@typeBoundFunctions,
				{
				    type        => "procedure"  ,
				    name        => $property->{'name'}."Analytic",
				    function    => &createNullFunction($build,{selfType => $className, attribute => "analytic", property => $property, intent => "inout"}),
				    returnType  => "\\void"     ,
				    arguments   => "",
				    description => "Mark the {\\normalfont \\ttfamily ".$property->{'name'}."} property of the {\\normalfont \\ttfamily ".$implementationIdentifier."} component as analtyically-solvable."
				}
				)
				unless ( $property->{'attributes'}->{'isVirtual'} );
			    # Create an inactive function unless this is a virtual property.
			    push(
				@typeBoundFunctions,
				{
				    type        => "procedure"  ,
				    name        => $property->{'name'}."Inactive",
				    function    => &createNullFunction($build,{selfType => $className, attribute => "inactive", property => $property, intent => "inout"}),
				    returnType  => "\\void"     ,
				    arguments   => "",
				    description => "Mark the {\\normalfont \\ttfamily ".$property->{'name'}."} property of the {\\normalfont \\ttfamily ".$implementationIdentifier."} component as inactive."
				}
				)
				unless ( $property->{'attributes'}->{'isVirtual'} );
			    # Create a "scale" function unless this is a virtual property.
			    push(
				@typeBoundFunctions,
				{
				    type        => "procedure"                      ,
				    name        => $property->{'name'}."Scale"      ,
				    function    => &createNullFunction($build,{selfType => $className, attribute => "scale", property => $property, intent => "inout"}),
				    returnType  => "\\void"                         ,
				    arguments   => &Galacticus::Build::Components::DataTypes::dataObjectDocName($property)."\\ value",
				    description => "Set the scale of the {\\normalfont \\ttfamily ".$property->{'name'}."} property of the {\\normalfont \\ttfamily ".$implementationIdentifier."} component."
				}
				)
				unless ( $property->{'attributes'}->{'isVirtual'} );
			    $propertiesCreated{$functionName} = 1;
			}
		    }
		}
	    }
	}
	# Create the type.
	$build->{'types'}->{'nodeComponent'.ucfirst($className)} = 
	{
	    name           => "nodeComponent".ucfirst($className),
	    comment        => "Type for the {\\normalfont \\ttfamily ".$className."} component class.",
	    isPublic       => 1,
	    extends        => "nodeComponent",
	    boundFunctions => \@typeBoundFunctions,
	};
    }
}

1;
