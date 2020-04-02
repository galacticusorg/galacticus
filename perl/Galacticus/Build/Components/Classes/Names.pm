# Contains a Perl module which provides naming functions for component classes.

package Galacticus::Build::Components::Classes::Names;
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
     classNames =>
     {
	 classIteratedFunctions =>
	     [
	      \&Class_Type
	     ]
     }
    );

sub Class_Type {
    # Generate functions to provide type names for component classes.
    my $build    = shift();
    $code::class = shift();
    my $classTypeName = "nodeComponent".ucfirst($code::class->{'name'});
    my $function =
    {
	type        => "type(varying_string) => name",
	name        => $classTypeName."Type",
	description => "Returns the type name for the ".$code::class->{'name'}." component class.",
	modules     =>
	    [
	     "ISO_Varying_String"
	    ],
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => $classTypeName,
		 attributes => [ "intent(in   )" ],
		 variables  => [ "self" ]
	     }
	    ]
    };    
    $function->{'content'}  = fill_in_string(<<'CODE', PACKAGE => 'code');
!$GLC attributes unused :: self
name='nodeComponent:{$class->{'name'}}'
CODE
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{$classTypeName}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    descriptor  => $function,
	    name        => "type"
	}
	);
}

1;
