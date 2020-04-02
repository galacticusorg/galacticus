# Contains a Perl module which provides various ODE solver-related functions for component hierarchy parent classes.

package Galacticus::Build::Components::Hierarchy::ODESolver;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Text::Template 'fill_in_string';
use List::ExtraUtils;
use Fortran::Utils;
use Galacticus::Build::Components::Utils;
use Galacticus::Build::Components::DataTypes;

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     hierarchyODESolver => 
     {
	 functions =>
	     [
	      \&Component_ODE_Name_From_Index
	     ]
     }
    );

sub Component_ODE_Name_From_Index {
    # Generate a function to return the name of a property given the index of that property in a generic component.
    my $build = shift();
    # Generate the function.
    my $function =
    {
	type        => "type(varying_string) => name",
	name        => "nodeComponentNameFromIndex",
	description => "Return the name of the property of given index for a {\\normalfont \\ttfamily nodeComponent} object.",
	modules     => 
	    [ 
	      "ISO_Varying_String"
	    ],
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => [ "intent(inout)" ],
		 variables  => [ "count" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "propertyType" ]
	     }
	    ]
    };
    # This generic (parent) node component class has no properties, so return an unknown name.
    $function->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
!$GLC attributes unused :: self, count, propertyType
name='?'
CODE
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{'nodeComponent'}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    descriptor  => $function,
	    name        => "nameFromIndex"
	}
	);	    
}

1;
