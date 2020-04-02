# Contains a Perl module which provides naming functions for component implementations.

package Galacticus::Build::Components::Implementations::Names;
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
     implementationNames =>
     {
	 implementationIteratedFunctions =>
	     [
	      \&Implementation_Type
	     ]
     }
    );

sub Implementation_Type {
    # Generate functions to provide type names for component implementations.
    my $build     = shift();
    $code::class  = shift();
    $code::member = shift();
    my $implementationTypeName = "nodeComponent".ucfirst($code::class->{'name'}).ucfirst($code::member->{'name'});
    my $function =
    {
	type        => "type(varying_string) => name",
	name        => $implementationTypeName."Type",
	description => "Returns the type name for the ".$code::member->{'name'}." implementation of the ".$code::class->{'name'}." component class.",
	modules     =>
	    [
	     "ISO_Varying_String"
	    ],
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => $implementationTypeName,
		 attributes => [ "intent(in   )" ],
		 variables  => [ "self" ]
	     }
	    ]
    };    
    $function->{'content'}  = fill_in_string(<<'CODE', PACKAGE => 'code');
!$GLC attributes unused :: self
name='nodeComponent:{$class->{'name'}}:{$member->{'name'}}'
CODE
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{$implementationTypeName}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    descriptor  => $function,
	    name        => "type"
	}
	);
}

1;
