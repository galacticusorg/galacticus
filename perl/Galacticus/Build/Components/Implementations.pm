# Contains a Perl module which handles component implementations during build.

package Galacticus::Build::Components::Implementations;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Data::Dumper;
use Sort::Topo;
use Scalar::Util;
use NestedMap;
use List::ExtraUtils;
use Galacticus::Build::Components::Utils qw(applyDefaults $verbosityLevel @booleanLabel);
use Galacticus::Build::Components::DataTypes;

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     implementations => 
     {
	 preValidate =>
	     [
	      \&Implementation_ID_List
	      ],
	 default     =>
	     [
	      \&Implementation_Defaults,
	      \&Null_Implementations   ,
	      \&Default_Full_Name
	     ],
	 gather      =>
	     [
	      \&Implementation_Dependencies    ,
	      \&Implementation_Parents         ,
	      \&Implementation_Bindings_Inherit
	     ],
	 types       =>
	     [
	      \&Build_Component_Implementations
	     ]
     }
    );

sub Default_Full_Name {
    # Set a fully-qualified name for each implementation.
    my $build = shift();
    # Iterate over components.
    $_->{'fullyQualifiedName'} = ucfirst($_->{'class'}).ucfirst($_->{'name'})
	foreach ( &List::ExtraUtils::hashList($build->{'components'}) );
}

sub Implementation_ID_List {
    # Create a list of component IDs.
    my $build = shift();
    # Construct a list of all component names.
    @{$build->{'componentIdList'}} = &List::ExtraUtils::sortedKeys($build->{'components'});
}

sub Implementation_Defaults {
    # Set defaults for implementations.
    my $build = shift();
    # Default settings for implementations.
    my %defaults =
	(
	 isDefault     => "booleanFalse",
	 bindings      =>
	 {
	     binding =>
	     {
		 isDeferred => "booleanFalse",
		 interface =>
		 {
		     self =>
		     {
			 pass => "booleanFalse"
		     }
		 }
	     }
	 },
	 createFunction =>
	 {
	     isDeferred => "booleanFalse"
	 },
	 output         =>
	 {
	     instances => "all"
	 }
	);
    # Iterate over implementations and apply all defaults.
    foreach my $implementation ( &List::ExtraUtils::hashList($build->{'components'}) ) {
	&applyDefaults($implementation,$_,$defaults{$_})
	    foreach ( keys(%defaults) );
    }
    # Record the default implementation for each class.
    foreach my $implementation ( &List::ExtraUtils::hashList($build->{'components'}) ) {
	$build->{'componentClasses'}->{$implementation->{'class'}}->{'defaultImplementation'} = $implementation->{'name'}
	    if ( $implementation->{'isDefault'} );
    }
}

sub Null_Implementations {
    # Create a null implementation for each class.
    my $build = shift();
    # Iterate over components to determine which classes need a null case building.
    my %classes;
    foreach my $implementation ( &List::ExtraUtils::hashList($build->{'components'}) ) {
	# Initialize this class if it hasn't been seen before.
	unless ( exists($classes{$implementation->{'class'}}) ) {
	    $classes{$implementation->{'class'}}->{'hasNull'   } = 0;
	    $classes{$implementation->{'class'}}->{'hasDefault'} = 0;
	}
	# Record if a null component already exists.
	$classes{$implementation->{'class'}}->{'hasNull'   } = 1
	    if ( $implementation->{'name'     } eq "null" );
	# Record if a default is already specified.
	$classes{$implementation->{'class'}}->{'hasDefault'} = 1
	    if ( $implementation->{'isDefault'} );
    }
    # Iterate over classes, creating null components as necessary.
    foreach my $class ( &List::ExtraUtils::sortedKeys(\%classes) ) {       
	# Test for pre-existing null component.
	if ( $classes{$class}->{'hasNull'} == 0 ) {
	    # No pre-existing null component is present, so simply insert one into the build data.
	    my $implementationName = ucfirst($class)."Null";
	    my $isDefault          = $classes{$class}->{'hasDefault'} ? 0 : 1;
	    $build->{'components'      }->{$implementationName}->{'class'                } = $class;
	    $build->{'components'      }->{$implementationName}->{'name'                 } = "null";
	    $build->{'components'      }->{$implementationName}->{'isDefault'            } = $isDefault; 
	    $build->{'componentClasses'}->{$class             }->{'defaultImplementation'} = "null"
		if ( $isDefault );
	    # Append this new component ID to the component ID list.
	    push(@{$build->{'componentIdList'}},$implementationName);
	    # Display a message.
	    if ( $verbosityLevel >= 1 ) {
		print "         --> Adding null implementation ";
		print "as default "
		    unless ( $classes{$class}->{'hasDefault'} );
		print "for ".$class." class\n";
	    }
	} elsif ( $verbosityLevel >= 1 ) {
	    # Advise that null components don't need to be explicitly specified.
	    print "         --> INFO: a pre-existing null component was found for the '".$class."' class,\n";
	    print "                   but would be built automatically.\n";
	}
    }
}

sub Implementation_Dependencies {
    # Order class members such that parent classes come before child classes.
    my $build = shift();
    # Iterate over classes
    print "         --> Sorting implentations into parent->child order:\n"
	if ( $verbosityLevel >= 1 );
    foreach my $className ( @{$build->{'componentClassList'}} ) {
	my %dependencies;
	foreach my $implementationName ( @{$build->{'componentClasses'}->{$className}->{'memberNames'}} ) {
	    my $implementationID = ucfirst($className).ucfirst($implementationName);
 	    my $implementation   = $build->{'components'}->{$implementationID};
	    push(@{$dependencies{$implementation->{'extends'}->{'name'}}},$implementationName)
		if ( exists($implementation->{'extends'}) );
	}
	# Sort member names into parent->child order.
	@{$build->{'componentClasses'}->{$className}->{'memberNames'}} =
	    &Sort::Topo::sort($build->{'componentClasses'}->{$className}->{'memberNames'},\%dependencies);
	# Create the same ordering in the list of members.
	my %members;
	$members{$_->{'name'}} = $_
	    foreach ( @{$build->{'componentClasses'}->{$className}->{'members'}} );
	@{$build->{'componentClasses'}->{$className}->{'members'}} = map {$members{$_}} @{$build->{'componentClasses'}->{$className}->{'memberNames'}};
	if ( $verbosityLevel >= 1 ) {
	    print "            --> ".$className.":\n";
	    print "               --> ".$_."\n"
		foreach ( @{$build->{'componentClasses'}->{$className}->{'memberNames'}} );
	}
    }
    # Update the component ID list such that it is in dependency order.
    @{$build->{'componentIdList'}} =
	nestedmap {
	    nestedmap {
		ucfirst($NestedMap::stack[1]).ucfirst($NestedMap::stack[0])		
	    } @{$build->{'componentClasses'}->{$NestedMap::stack[0]}->{'memberNames'}}	    
    } @{$build->{'componentClassList'}};
}

sub Implementation_Parents {
    # Create links to parent implementations.
    my $build = shift();
    # Iterate over implementations.
    foreach my $implementation ( &List::ExtraUtils::hashList($build->{'components'}) ) {
	# For extensions, locate and link to the parent.
	if ( exists($implementation->{'extends'}) ) {
	    my $parentIdentifier = ucfirst($implementation->{'extends'}->{'class'}).ucfirst($implementation->{'extends'}->{'name'});
	    $implementation->{'extends'}->{'implementation'} = $build->{'components'}->{$parentIdentifier};
	}
    }
}

sub Implementation_Bindings_Inherit {
    # Inherit bindings from any parent implementations.
    my $build = shift();
    # Iterate over implementations.
    foreach my $implementation ( &List::ExtraUtils::hashList($build->{'components'}) ) {
	# For extensions, copy any binding from the parent class.
	if ( exists($implementation->{'extends'}) ) {
	    foreach my $parentBinding ( @{$implementation->{'extends'}->{'implementation'}->{'bindings'}->{'binding'}} ) {
		push
		    (
		     @{$implementation->{'bindings'}->{'binding'}},
		     $parentBinding
		    )
		    unless ( 
			! $parentBinding->{'isDeferred'}
			||
			grep {$_->{'method'} eq $parentBinding->{'method'}} @{$implementation->{'bindings'}->{'binding'}}
		    );		    
	    }
	}
    }
}

sub Build_Component_Implementations {
    # Generate a class for each component implementation.
    my $build = shift();
    # Iterate over implementations.
    foreach my $implementation ( &List::ExtraUtils::hashList($build->{'components'}) ) {
    	# Determine the name of the class which this component extends (use the "nodeComponent" class by default).
    	my $extensionOf = 
	    exists($implementation->{'extends'})
	    ?
	    "nodeComponent".$implementation->{'extends'}->{'implementation'}->{'fullyQualifiedName'}
	    : 
	    "nodeComponent".ucfirst($implementation->{'class'});
     	# Create data objects to store all of the linked data for this component.
	my @dataContent;
    	foreach ( &List::ExtraUtils::sortedKeys($implementation->{'content'}->{'data'}) ) {
	    (my $typeDefinition, my $typeLabel) = &Galacticus::Build::Components::DataTypes::dataObjectDefinition($implementation->{'content'}->{'data'}->{$_});
	    $typeDefinition->{'variables'} = [ $_ ];
	    push(
		@dataContent,
		$typeDefinition
		);
    	}
	# Create a list for type-bound functions.
	my @typeBoundFunctions;
     	# If this component has bindings defined, scan through them and create an appropriate method.
    	if ( exists($implementation->{'bindings'}) ) {
    	    foreach ( @{$implementation->{'bindings'}->{'binding'}} ) {
		my %function = (
		    type        => "procedure",
		    name        => $_->{'method'}
		    );
		if ( ! $_->{'isDeferred'} ) {
		    # Binding is not deferred, simply map to the given function. For deferred functions a suitable wrapper will be attached later.
		    $function{'function'} = $_->{'function'};
		    foreach my $attribute ( "description", "returnType", "arguments" ) {
			$function{$attribute} = $_->{$attribute}
		        if ( exists($_->{$attribute}) );
		    }
		    push(@typeBoundFunctions,\%function);
		}
	    }
	}
	# Iterate over properties.
	foreach my $propertyName ( &List::ExtraUtils::sortedKeys($implementation->{'properties'}->{'property'}) ) {
	    # Get the property.
	    my $property = $implementation->{'properties'}->{'property'}->{$propertyName};
	    push
		(
		 @typeBoundFunctions,
		 {
		     type     => "procedure", 
		     pass     => "nopass",
		     name     => $propertyName."IsGettable", 
		     function => "Boolean_".ucfirst($booleanLabel[$property->{'attributes'}->{'isGettable'}])
		 },
		 {
		     type     => "procedure",
		     pass     => "nopass", 
		     name     => $propertyName."IsSettable",
		     function => "Boolean_".ucfirst($booleanLabel[$property->{'attributes'}->{'isSettable'}])
		 }
		);
	}
	# Create the type.
	$build->{'types'}->{'nodeComponent'.ucfirst($implementation->{'fullyQualifiedName'})} = 
	{
	    name           => "nodeComponent".ucfirst($implementation->{'fullyQualifiedName'}),
	    comment        => "Class for the ".$implementation->{'name'}." implementation of the ".$implementation->{'class'}." component.",
	    isPublic       => 1,
	    extends        => $extensionOf,
	    boundFunctions => \@typeBoundFunctions,
	    dataContent    => \@dataContent
	};
    }
}

1;
