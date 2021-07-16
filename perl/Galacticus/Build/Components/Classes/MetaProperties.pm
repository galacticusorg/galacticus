# Contains a Perl module which provides handling of meta-properties for component classes.

package Galacticus::Build::Components::Classes::MetaProperties;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Text::Template 'fill_in_string';
use List::ExtraUtils;
use Galacticus::Build::Components::Utils;
use Galacticus::Build::Components::DataTypes;

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     classMetaProperties =>
     {
	 classIteratedFunctions =>
	     [
	      \&Class_Meta_Property_Get,
	      \&Class_Meta_Property_Set
	     ]
     }
    );

sub Class_Meta_Property_Get {
    # Generate a function to return the value of an indexed meta-property.
    my $build            = shift();
    $code::class         = shift();
    $code::classTypeName = "nodeComponent".ucfirst($code::class->{'name'});
    my $function =
    {
	type        => "double precision",
	name        => $code::classTypeName."MetaPropertyGet",
	description => "Return the value of a meta-property of a ".$code::classTypeName." component given its ID.",
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => $code::classTypeName,
		 attributes => [ "intent(in   )" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "metaPropertyID" ]
	     }
	    ],
	content     => $code::classTypeName."MetaPropertyGet=self%metaProperties(metaPropertyID)\n"
    };
    # Insert a type-binding for this function.
    push(
	@{$build->{'types'}->{$code::classTypeName}->{'boundFunctions'}},
	{
	    type        => "procedure",
	    descriptor  => $function,
	    name        => "metaPropertyGet", 
	}
	);
}

sub Class_Meta_Property_Set {
    # Generate a function to set the value of an indexed meta-property.
    my $build            = shift();
    $code::class         = shift();
    $code::classTypeName = "nodeComponent".ucfirst($code::class->{'name'});
    my $function =
    {
	type        => "void",
	name        => $code::classTypeName."MetaPropertySet",
	description => "Set the value of a meta-property of a ".$code::classTypeName." component given its ID.",
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => $code::classTypeName,
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "metaPropertyID" ]
	     },
	     {
		 intrinsic  => "double precision",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "metaPropertyValue" ]
	     }
	    ],
	content     => "self%metaProperties(metaPropertyID)=metaPropertyValue\n"
    };
    # Insert a type-binding for this function.
    push(
	@{$build->{'types'}->{$code::classTypeName}->{'boundFunctions'}},
	{
	    type        => "procedure",
	    descriptor  => $function,
	    name        => "metaPropertySet", 
	}
	);
}

1;
