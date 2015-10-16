# Contains a Perl module which creates get/set/evolve code for component property attributes.

package Attributes;
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

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     attributes => {
	 postValidate =>
	     [
	      \&Validate_Boolean               ,
	      \&Validate_Evolvable_Intrinsics  ,
	      \&Validate_Deferreds_Functionless
	     ]
     }
    );

sub Validate_Boolean {
    # Validate that property attributes are "true" or "false".
    my $build = shift();
    # Iterate over components.
    foreach my $component ( &ExtraUtils::hashList($build->{'components'}) ) {
	# Iterate over all properties belonging to this component.	
	next 
	    unless ( exists($component->{'properties'}) );
	foreach my $property ( &ExtraUtils::hashList($component->{'properties'}->{'property'}, keyAs => 'name' ) ) {
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
    foreach my $component ( &ExtraUtils::hashList($build->{'components'}) ) {
	# Iterate over all properties belonging to this component.	
	next 
	    unless ( exists($component->{'properties'}) );
	foreach my $property ( &ExtraUtils::hashList($component->{'properties'}->{'property'}, keyAs => 'name' ) ) {
	    # Assert that non-real intrinsic properties be not evolvable.
	    die(
		"Validate_Evolvable_Intrinsics: non-real intrinsic property '".
		$property->{'name'}                                           .
		"' of '"                                                      .
		$component->{'class'}.ucfirst($component->{'name'})           .
		"' component can not be evolvable"
		)
		if (
		    &Utils::isIntrinsic($property->{'type'})
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
    foreach my $component ( &ExtraUtils::hashList($build->{'components'}) ) {
	# Iterate over all properties belonging to this component.	
	next 
	    unless ( exists($component->{'properties'}) );
	foreach my $property ( &ExtraUtils::hashList($component->{'properties'}->{'property'}, keyAs => 'name' ) ) {
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

1;
