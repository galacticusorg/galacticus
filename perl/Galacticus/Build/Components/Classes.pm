# Contains a Perl module which handles component classes during build.

package Classes;
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
use Data::Dumper;
require List::ExtraUtils;
require Galacticus::Build::Components::Utils;
require Galacticus::Build::Components::DataTypes;

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
    foreach my $implementation ( &ExtraUtils::hashList($build->{'components'}) ) {
	# Get the name of the component class.
	my $className          = $implementation->{'class'};
	# Get the name of the implementation.
	my $implementationName = $implementation->{'name' };
	# Append the component name to the list of members for its class.
	push(
	    @{$build->{'componentClasses'}->{$className}->{'members'}},
	    $implementationName
	    );
    }
    # Construct a list of component class names.
    @{$build->{'componentClassList'}} = &ExtraUtils::sortedKeys($build->{'componentClasses'});
    # Output report if sufficiently verbose.
    if ( $Utils::verbosityLevel > 0 ) {
	print "         --> Found the following component classes and implementations:\n";
	foreach my $className ( @{$build->{'componentClassList'}} ) {
	    print "            --> ".$className."\n";
	    foreach my $implementationName ( @{$build->{'componentClasses'}->{$className}->{'members'}} ) {
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
    foreach my $class ( &ExtraUtils::hashList($build->{'componentClasses'}, keyAs => 'name') ) {
	my $className  = $class->{'name'};
	# Define a hash to record which properties have already been created.
	my %propertiesCreated;
	# Create a list for type-bound functions.
	my @typeBoundFunctions;
  	# Insert definitions for each method associated with a component implementation of this component class.
    	foreach my $implementationName ( @{$class->{'members'}} ) {
	    # Construct a fully-qualified name for this implementation.
	    my $implementationIdentifier = ucfirst($className).ucfirst($implementationName);
	    # Get the implementation.
	    my $implementation           = $build->{'components'}->{$implementationIdentifier};
	    # Iterate over properties beloning to this implementation.
	    foreach my $property ( &ExtraUtils::hashList($implementation->{'properties'}->{'property'}, keyAs => 'name') ) {
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
		    (my $intrinsic,my $type,my $attributes) = &DataTypes::dataObjectPrimitiveName($property);
		    $type .= $property->{'rank'}."InOut";
		    # Record the null bindings needed.
		    $build->{'nullProperties'}->{$className}->{$type} = 
		    {
			type   => $property->{'type'},
			rank   => $property->{'rank'},
			intent => "inout"
		    };
		    # # Create the "isSettable" function.
		    # $functionName = $property->{'name'}."IsSettable";
		    # unless ( exists($propertiesCreated{$functionName}) ) {
		    # 	push(
		    # 	    @typeBoundFunctions,
		    # 	    {
		    # 		type        => "procedure"    ,
		    # 		pass        => "nopass"       ,
		    # 		name        => $functionName  , 
		    # 		function    => "Boolean_False",
		    # 		returnType  => "\\logicalzero",
		    # 		arguments   => ""             ,
		    # 		description => "Specify whether the {\\normalfont \\ttfamily ".$property->{'name'}."} property of the {\\normalfont \\ttfamily ".$className."} component is settable."
		    # 	    }
		    # 	    );
		    # 	$propertiesCreated{$functionName} = 1;
		    # }
		    # Handle set functions and related functions.
		    if ( $property->{'attributes'}->{'isSettable'} ) {
			# Create a "set" function if one does not already exist.
			$functionName = $property->{'name'}."Set";
			unless ( exists($propertiesCreated{$functionName}) ) {
			    my $boundTo =
				$property->{'setFunction'}->{'bindsTo'} eq "componentClass" 
				?
				# A setFunction was specified that binds to the component class, so bind to it here.
				$property->{'setFunction'}->{'content'}
			        :
				# Create a binding to a null function here. 
				$className."NullBindingSet".$type;
			    push(
				@typeBoundFunctions,
				{
				    type        => "procedure"  , 
				    name        => $functionName, 
				    function    => $boundTo     ,
				    returnType  => "\\void"     ,
				    arguments   => &DataTypes::dataObjectDocName($property)."\\ value",
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
			    $build->{'nullProperties'}->{$className}->{"Integer0In"} =
			    {
				type   => "integer",
				rank   => 0        ,
				intent => "in"
			    };
			    push(
				@typeBoundFunctions,
				{
				    type        => "procedure"                       ,
				    name        => $functionName                     ,
				    function    => $className."NullBindingInteger0In",
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
			    unless ( 
				$property->{'attributes'}->{'isDeferred' } =~ m/rate/
				&&
				$property->{'attributes'}->{'bindsTo'    } eq "top"
				) {
				my $boundTo =
				    $property->{'attributes'}->{'createIfNeeded'} 
				    ?
				    $className.ucfirst($property->{'name'})."Rate"
				    :
				    $className."NullBindingRate".$type;
				push(
				    @typeBoundFunctions,
				    {
					type        => "procedure"  ,
					name        => $functionName,
					function    => $boundTo     , 
					returnType  => "\\void"     ,
					arguments   => &DataTypes::dataObjectDocName($property)."\\ value",
					description => "Cumulate to the rate of the {\\normalfont \\ttfamily ".$property->{'name'}."} property of the {\\normalfont \\ttfamily ".$implementationIdentifier."} component."
				    }
				    );
			    }
			    # Create a "scale" function unless this is a virtual property.
			    push(
				@typeBoundFunctions,
				{
				    type        => "procedure"                      ,
				    name        => $property->{'name'}."Scale"      ,
				    function    => $className."NullBindingSet".$type,
				    returnType  => "\\void"                         ,
				    arguments   => &DataTypes::dataObjectDocName($property)."\\ value",
				    description => "Set the scale of the {\\normalfont \\ttfamily ".$property->{'name'}."} property of the {\\normalfont \\ttfamily ".$implementationIdentifier."} component."
}
				)
				unless ( $property->{'isVirtual'} );
			    $propertiesCreated{$functionName} = 1;
			}
		    }
		    # Add any bindings which bind at the component class level.
		    if ( exists($implementation->{'bindings'}) ) {
			foreach ( @{$implementation->{'bindings'}->{'binding'}} ) {
			    if ( $_->{'bindsTo'} eq "componentClass" ) {
				my %function = 
				    (
				     type => "procedure",
				     name => $_->{'method'},
				    );
				if ( ! $_->{'isDeferred'} ) {
				    # Binding is not deferred, simply map to the given function.
				    $function{'function'} = $_->{'function'};
				} else {
				    # Binding is deferred, map to a suitable wrapper function.
				    $function{'function'} = $className.ucfirst($_->{'method'});
				    # Also add bindings to functions to set and test the deferred function.
				    my %setFunction =
					(
					 type        => "procedure",
					 pass        => "nopass",
					 name        => $_->{'method'}."Function",
					 function    => $className.$_->{'method'}."DeferredFunctionSet",
					 returnType  => "\\void",
					 arguments   => "procedure(".$className.$_->{'method'}."Interface) deferredFunction",
					 description => "Set the function for the deferred {\\normalfont \\ttfamily ".$_->{'method'}."} propert of the {\\normalfont \\ttfamily ". $className."} component."
					);
				    my %testFunction =
					(
					 type        => "procedure",
					 pass        => "nopass",
					 name        => $_->{'method'}."FunctionIsSet",
					 function    => $className.$_->{'method'}."DfrrdFnctnIsSet",
					 returnType  => "\\logicalzero",
					 arguments   => "",
					 description => "Specify whether the deferred function for the {\\normalfont \\ttfamily ".$_->{'method'}."} propert of the {\\normalfont \\ttfamily ". $className."} component has been set."
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
