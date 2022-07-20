# Contains a Perl module which handles setting class default values for properties.

package Galacticus::Build::Components::Classes::Defaults;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Data::Dumper;
use List::ExtraUtils;
use Galacticus::Build::Components::Utils qw(&isIntrinsic %intrinsicNulls);
use Text::Template 'fill_in_string';

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     classesDefaults => 
     {
	 classIteratedFunctions =>
	     [
	      \&Class_Property_Is_Gettable ,
	      \&Class_Property_Default     ,
	      \&Class_Property_Rate_Default
	     ]
     }
    );

sub Class_Property_Is_Gettable {
    # Build functions that specify if a property value is gettable from the uninstantiated component class.
    my $build       = shift();
    $code::class    = shift();
    # Iterate over properties with class defaults.
    foreach $code::property ( &List::ExtraUtils::hashList(&Class_Members_By_Property($code::class)) ) {
	# Build the function.
	my $function =
	{
	    type        => "logical",
	    name        => $code::class->{'name'}.ucfirst($code::property->{'name'})."IsGettable",
	    description => "Returns true if the {\\normalfont \\ttfamily ".$code::property->{'name'}."} property is gettable for the {\\normalfont \\ttfamily ".$code::class->{'name'}."} component class."
	};
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
{$class->{'name'}.ucfirst($property->{'name'})}IsGettable=.false.
CODE
        # Iterate over class members which provide this property.
        foreach $code::member ( @{$code::property->{'members'}} ) {
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if (nodeComponent{ucfirst($class->{'name'}).ucfirst($member->{'name'})}IsActiveValue) {$class->{'name'}.ucfirst($property->{'name'})}IsGettable=.true.
CODE
        }
	# Insert a type-binding for this function into the treeNode type.
	push(
	    @{$build->{'types'}->{'nodeComponent'.ucfirst($code::class->{'name'})}->{'boundFunctions'}},
	    {
		type        => "procedure", 
		descriptor  => $function,
		name        => $code::property->{'name'}."IsGettable",
		pass        => "nopass"
	    }
	);
    }
}

sub Class_Property_Default {
    # Build functions that return the class default value for properties.
    my $build       = shift();
    $code::class    = shift();
    # Iterate over properties with class defaults.
    foreach $code::property ( &List::ExtraUtils::hashList(&Class_Members_By_Property($code::class)) ) {
	# Build the function.
	(my $functionTypeDescriptor, $code::functionTypeLabel) = &Galacticus::Build::Components::DataTypes::dataObjectDefinition($code::property->{'property'}->{'data'});
	my $function =
	{
	    type        => 
		$functionTypeDescriptor->{'intrinsic'}.
		(
		      exists(       $functionTypeDescriptor->{'type'      } )
		 ?
		             "(" .  $functionTypeDescriptor->{'type'      }  .")"
		 :
		 ""
		).
		(
		      exists(       $functionTypeDescriptor->{'attributes'} )
		 ? 
		 ", ".join  (", ",@{$functionTypeDescriptor->{'attributes'}}) 
		 :
		 ""
		).
		" => classDefault",
	    name        => $code::class->{'name'}.ucfirst($code::property->{'name'}),
	    description => "Returns the default value for the {\\normalfont \\ttfamily ".$code::property->{'name'}."} property for the {\\normalfont \\ttfamily ".$code::class->{'name'}."} component class.",
	    variables =>
		[
		 {
		     intrinsic  => "class",
		     type       => "nodeComponent".ucfirst($code::class->{'name'}),
		     attributes => [ "intent(inout)" ],
		     variables   => [ "self" ]
		 }
		]
	};
	# Discover all components required to set defaults and all modules needed.
	my $requiredComponents;
	my $requiredModules;
	foreach my $member ( @{$code::property->{'members'}} ) {
	    my $implementationProperty = $member->{'properties'}->{'property'}->{$code::property->{'name'}};
	    next
		unless ( exists($implementationProperty->{'classDefault'}->{'code'}) );
	    $requiredModules->{$_} = 1
		foreach ( @{$implementationProperty->{'classDefault'}->{'modules'}} );
	    $requiredModules->{'Memory_Management,only:allocateArray'} = 1
		if ( exists($implementationProperty->{'classDefault'}->{'count'}) );
	    my $classDefault = $implementationProperty->{'classDefault'}->{'code'};
	    while ( $classDefault =~ m/self([a-zA-Z]+)\s*%/ ) {
		$requiredComponents->{'all'            }->{$1} = 1;
		$requiredComponents->{$member->{'name'}}->{$1} = 1;
		$classDefault                                  =~ s/self([a-zA-Z]+)\s*%//;
	    }
	}
	# Insert any required modules.
	@{$function->{'modules'}} = keys(%{$requiredModules})
	    if ( defined($requiredModules) );
	# Insert a pointer to the host node if required.
	push(
	    @{$function->{'variables'}},
	    {
		intrinsic  => "type",
		type       => "treeNode",
		attributes => [ "pointer" ],
		variables  => [ "node" ]
	    }
	    )
	    if ( scalar(keys(%{$requiredComponents->{'all'}})) > 0 );
	# Insert pointers to all node components required.
	push(
	    @{$function->{'variables'}},
	    {
		intrinsic  => "class",
		type       => "nodeComponent".ucfirst($_),
		attributes => [ "pointer" ],
		variables  => [ "self".ucfirst($_) ]
	    }
	    )
	    foreach ( &List::ExtraUtils::sortedKeys($requiredComponents->{'all'}) );
	# Generate the function code.
	$function->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
!$GLC attributes unused :: self
CODE
	# Iterate over class members which provide a class default for this property.
	foreach $code::member ( @{$code::property->{'members'}} ) {
	    $code::implementationProperty = $code::member->{'properties'}->{'property'}->{$code::property->{'name'}};
	    next
		unless ( exists($code::implementationProperty->{'classDefault'}->{'code'}) );
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if (nodeComponent{ucfirst($class->{'name'}).ucfirst($member->{'name'})}IsActiveValue) then
CODE
	    if ( scalar(keys(%{$requiredComponents->{$code::member->{'name'}}})) > 0 ) {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
   node => self%host()
CODE
	    }
	    foreach $code::component ( &List::ExtraUtils::sortedKeys($requiredComponents->{$code::member->{'name'}}) ) {
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
   self{$component} => node%{$component}()
CODE
	    }
	    if ( exists($code::implementationProperty->{'classDefault'}->{'count'}) ) {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
   call allocateArray(classDefault,[{$implementationProperty->{'classDefault'}->{'count'}}])
CODE
	    }
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
   classDefault={$implementationProperty->{'classDefault'}->{'code'}}
   return
end if
CODE
	}
	# Insert code to return zero values by default.	
	if ( $code::property->{'property'}->{'data'}->{'rank'} == 0 ) {
	    if ( &isIntrinsic($code::property->{'property'}->{'data'}->{'type'}) ) {
		%code::intrinsicNulls = %intrinsicNulls;
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
classDefault={$intrinsicNulls{$property->{'property'}->{'data'}->{'type'}}}
CODE
	    } else {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
call classDefault%reset()
CODE
	    }
	} else {
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
    classDefault=null{$functionTypeLabel.$property->{'property'}->{'data'}->{'rank'}}d
CODE
	}
	# Insert a type-binding for this function into the treeNode type.
	push(
	    @{$build->{'types'}->{'nodeComponent'.ucfirst($code::class->{'name'})}->{'boundFunctions'}},
	    {
		type        => "procedure", 
		descriptor  => $function,
		name        => $code::property->{'name'}
	    }
	);
    }
}

sub Class_Property_Rate_Default {
    # Build functions that return the class default (0) rate for properties.
    my $build       = shift();
    $code::class    = shift();
    # Iterate over properties.
    foreach $code::property ( &List::ExtraUtils::hashList(&Class_Members_By_Property($code::class)) ) {
	# Build the function.
	(my $functionTypeDescriptor, $code::functionTypeLabel) = &Galacticus::Build::Components::DataTypes::dataObjectDefinition($code::property->{'property'}->{'data'});
	my $function =
	{
	    type        => 
		$functionTypeDescriptor->{'intrinsic'}.
		(
		      exists(       $functionTypeDescriptor->{'type'      } )
		 ?
		             "(" .  $functionTypeDescriptor->{'type'      }  .")"
		 :
		 ""
		).
		(
		      exists(       $functionTypeDescriptor->{'attributes'} )
		 ? 
		 ", ".join  (", ",@{$functionTypeDescriptor->{'attributes'}}) 
		 :
		 ""
		).
		" => classDefault",
	    name        => $code::class->{'name'}.ucfirst($code::property->{'name'})."RateGet",
	    description => "Returns a zero rate for the {\\normalfont \\ttfamily ".$code::property->{'name'}."} property for the {\\normalfont \\ttfamily ".$code::class->{'name'}."} component class.",
	    variables =>
		[
		 {
		     intrinsic  => "class",
		     type       => "nodeComponent".ucfirst($code::class->{'name'}),
		     attributes => [ "intent(inout)" ],
		     variables   => [ "self" ]
		 }
		]
	};
	# Generate the function code.
	$function->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
!$GLC attributes unused :: self
CODE
	# Insert code to return zero values.
	if ( $code::property->{'property'}->{'data'}->{'rank'} == 0 ) {
	    if ( &isIntrinsic($code::property->{'property'}->{'data'}->{'type'}) ) {
		%code::intrinsicNulls = %intrinsicNulls;
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
classDefault={$intrinsicNulls{$property->{'property'}->{'data'}->{'type'}}}
CODE
	    } else {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
call classDefault%reset()
CODE
	    }
	} else {
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
    classDefault=null{$functionTypeLabel.$property->{'property'}->{'data'}->{'rank'}}d
CODE
	}
	# Insert a type-binding for this function into the treeNode type.
	push(
	    @{$build->{'types'}->{'nodeComponent'.ucfirst($code::class->{'name'})}->{'boundFunctions'}},
	    {
		type        => "procedure", 
		descriptor  => $function,
		name        => $code::property->{'name'}."RateGet"
	    }
	);
    }
}

sub Class_Members_By_Property {
    # Return a data structure containing class members organized by the properties which they provide.
    my $class = shift();
    # Construct lists of all members which provide a given property with a class default.
    my $properties;
    # Iterate over class member implementations, and populate the structure with all class members that provide each property.
    foreach my $member ( @{$class->{'members'}} ) {
	# Iterate over all properties belonging to this member.	
	foreach my $property ( &List::ExtraUtils::hashList($member->{'properties'}->{'property'}, keyAs => 'name' ) ) {
	    $properties->{$property->{'name'}}->{'name'    } = $property->{'name'};
	    $properties->{$property->{'name'}}->{'property'} = $property;
	    push(@{$properties->{$property->{'name'}}->{'members'}},$member)
		if ( exists($property->{'classDefault'}) );
	}
    }
    return $properties;
}

1;
