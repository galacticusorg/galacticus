# Contains a Perl module which creates get/set/evolve code for component property attributes.

package Galacticus::Build::Components::Attributes;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Data::Dumper;
use Scalar::Util 'reftype';
use List::ExtraUtils;
use Galacticus::Build::Components::Utils qw(isIntrinsic);

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     attributes => {
	 preValidate  =>
	     [
	      \&Validate_Deferreds_Functionless
	     ],
	 default      =>
	     [
	      \&Default_Functions              ,
	     ],
	 postValidate =>
	     [
	      \&Validate_Boolean               ,
	      \&Validate_Evolvable_Intrinsics 
	     ]
     }
    );

sub Validate_Boolean {
    # Validate that property attributes are "true" or "false".
    my $build = shift();
    # Iterate over components.
    foreach my $component ( &List::ExtraUtils::hashList($build->{'components'}) ) {
	# Iterate over all properties belonging to this component.	
	next 
	    unless ( exists($component->{'properties'}) );
	foreach my $property ( &List::ExtraUtils::hashList($component->{'properties'}->{'property'}, keyAs => 'name' ) ) {
	    # Assert that property attributes must be "false" or "true", and convert to boolean form.
	    foreach ( "isSettable", "isGettable", "isEvolvable" ) {
		if ( exists($property->{'attributes'}->{$_}) ) {
		    if      ( $property->{'attributes'}->{$_} eq "true"  ) {
			$property->{'attributes'}->{$_} = 1;
		    } elsif ( $property->{'attributes'}->{$_} eq "false" ) {
			$property->{'attributes'}->{$_} = 0;
		    } else {
			die(
			    "Validate_Boolean: value of '"                     .
			    $_                                                 . 
			    "' attribute of '"                                 .
			    $property->{'name'}                                .
			    "' of '"                                           .
			    $component->{'class'}.ucfirst($component->{'name'}).
			    "' must be either 'true' or 'false'"
			    );
		    }
		}
	    }
	}
    }
}

sub Validate_Evolvable_Intrinsics {
    # Validate that evolvable intrinsic properties are of "double" type only.
    my $build = shift();
    # Iterate over components.
    foreach my $component ( &List::ExtraUtils::hashList($build->{'components'}) ) {
	# Iterate over all properties belonging to this component.	
	next 
	    unless ( exists($component->{'properties'}) );
	foreach my $property ( &List::ExtraUtils::hashList($component->{'properties'}->{'property'}, keyAs => 'name' ) ) {
	    # Assert that non-real intrinsic properties be not evolvable.
	    die(
		"Validate_Evolvable_Intrinsics: non-real intrinsic property '".
		$property->{'name'}                                           .
		"' of '"                                                      .
		$component->{'class'}.ucfirst($component->{'name'})           .
		"' component can not be evolvable"
		)
		if (
		    &isIntrinsic($property->{'type'})
		    &&
		    $property->{'type'      }                  ne "double"
		    &&
		    $property->{'attributes'}->{'isEvolvable'}
		);
	}
    }
}

sub Validate_Deferreds_Functionless {
    # Validate that deferred attributes of properties do not have functions specified at build time.
    my $build = shift();
    # Iterate over components.
    foreach my $component ( &List::ExtraUtils::hashList($build->{'components'}) ) {
	# Iterate over all properties belonging to this component.	
	next 
	    unless ( exists($component->{'properties'}) );
	foreach my $property ( &List::ExtraUtils::hashList($component->{'properties'}->{'property'}, keyAs => 'name' ) ) {
	    my @deferredMethods = split(/:/,$property->{'attributes'}->{'isDeferred'})
		if ( exists($property->{'attributes'}->{'isDeferred'}) );
	    foreach ( @deferredMethods ) {
		die(
		    "Validate_Deferreds_Functionless(): cannot specify '".
		    $_                                                   .
		    "Function' when '"                                   .
		    $_                                                   .
		    "' method is deferred for property '"                .
		    $property->{'name'}                                  .
		    "' of component '"                                   .
		    $component->{'class'}.ucfirst($component->{'name'})  .
		    "'"
		    )
			if ( exists($property->{$_."Function"}) );
		}
	}
    }
}

sub Default_Functions {
    # Set default functions for property attributes.
    my $build = shift();
    # Iterate over implementations.
    foreach my $implementation ( &List::ExtraUtils::hashList($build->{'components'}) ) {
	my $componentIdentifier = 
	    ucfirst($implementation->{'class'}).
	    ucfirst($implementation->{'name' });
	# Iterate over all properties belonging to this component.	
	next 
	    unless ( exists($implementation->{'properties'}) );
	foreach my $property ( &List::ExtraUtils::hashList($implementation->{'properties'}->{'property'}, keyAs => 'name' ) ) {	    
	    # Rate function.
	    $property->{'rateFunction'} = 
		$componentIdentifier         .
		ucfirst($property->{'name' }).
		"Rate"
		unless ( exists($property->{'rateFunction'}) );
	    # Get/Set functions.
	    foreach ( "get", "set" ) {
		if ( exists($property->{$_.'Function'}) ) {
		    # A function element was specified.
		    unless ( reftype($property->{$_.'Function'}) ) {
			# The function element is simply the function name. Replace with a structure.
			$property->{$_.'Function'} = 
			{
			    content => $property->{$_.'Function'}
			};
		    }
		    # Since function was specified, we will not need to build a function.
		    $property->{$_.'Function'}->{'build'} = 0;
		} else {
		    # No function element was specified, assign a default function and record that a function must be built.
		    $property->{$_.'Function'} = 
		    {
			content => 
			    lcfirst($componentIdentifier          ).
			    ucfirst($property           ->{'name'}).
			    ucfirst($_                            ),
			build   => 1
		    };
		}
	    }
	}
    }
}

1;
