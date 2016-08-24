# Contains a Perl module which provides various utility functions for component hierarchy parent classes.

package Galacticus::Build::Components::Hierarchy::Utils;
use strict;
use warnings;
use utf8;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use Text::Template 'fill_in_string';
use List::ExtraUtils;
use Fortran::Utils;
use Galacticus::Build::Components::Utils qw(isIntrinsic);
use Galacticus::Build::Components::DataTypes;

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     hierarchyUtils => 
     {
	 functions =>
	     [
	      \&Component_Assign
	     ]
     }
    );

sub Component_Assign {
    # Generate a function to assign one component to another.
    my $build = shift();
    # Generate the function.
    my $function =
    {
	type        => "void",
	name        => "nodeComponentAssign",
	description => "Assign a node component to another node component.",
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent",
		 attributes => [ "intent(  out)" ],
		 variables  => [ "to" ]
	     },
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "from" ]
	     }
	    ]
    };
    # Create the code.
    $function->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
to%hostNode => from%hostNode
select type (to)
CODE
    # Iterate over component classes.
    foreach $code::class ( &List::ExtraUtils::hashList($build->{'componentClasses'}) ) {
	# Iterate over class member implementations.
	foreach $code::member ( @{$code::class->{'members'}} ) {
	    $code::implementationTypeName = "nodeComponent".ucfirst($code::class->{'name'}).ucfirst($code::member->{'name'});
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
type is ({$implementationTypeName})
   select type (from)
   type is ({$implementationTypeName})
CODE
	    # Iterate over properties.
	    foreach $code::property ( &List::ExtraUtils::hashList($code::member->{'properties'}->{'property'}) ) {
		# Skip virtual properties.
		next
		    if ( $code::property->{'attributes'}->{'isVirtual'} );		
		if ( &isIntrinsic($code::property->{'data'}->{'type'}) ) {
		    if ( $code::property->{'data'}->{'type'} eq "double" && $code::property->{'data'}->{'rank'} > 0 ) {
			$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
      if (allocated(to%{$property->{'name'}}Data)) call Dealloc_Array(to%{$property->{'name'}}Data)
CODE
		    }
		} else {
			$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
      call to%{$property->{'name'}}Data%destroy()
CODE
		}
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
      to%{$property->{'name'}}Data=from%{$property->{'name'}}Data
CODE
	    }
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
   end select
CODE
	}
    }
    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
end select
CODE
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{'nodeComponent'}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    descriptor  => $function,
	    name        => "assign"
	},
	{
	    type        => "generic",
	    name        => "assignment(=)", 
	    function    => "assign" 
	}
	);	    
}

1;
