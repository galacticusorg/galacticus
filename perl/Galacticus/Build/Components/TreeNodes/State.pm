# Contains a Perl module which provides various state functions for tree nodes.

package Galacticus::Build::Components::TreeNodes::State;
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
     treeNodeState =>
     {
	 functions =>
	     [
	      \&Tree_Node_Size_Of
	     ]
     }
    );

sub Tree_Node_Size_Of {
    # Generate a function to compute the size of a tree node.
    my $build    = shift();
    my $function =
    {
	type        => "integer(c_size_t)",
	name        => "treeNodeSizeOf",
	description => "Compute the size (in bytes) of the tree node.",
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => "treeNode",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "integer",
		 variables  => [ "i" ]
	     },
	     {
		 intrinsic  => "class",
		 type       => "nodeEvent",
		 attributes => [ "pointer" ],
		 variables  => [ "event" ]
	     }
	    ]
    };    
    $function->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
treeNodeSizeOf=sizeof(self)
CODE
    # Iterate over all component classes
    foreach $code::class ( &List::ExtraUtils::hashList($build->{'componentClasses'}) ) {
	next
	    unless ( grep {$code::class->{'name'} eq $_} @{$build->{'componentClassListActive'}} );
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if (allocated(self%component{ucfirst($class->{'name'})})) then
  do i=1,size(self%component{ucfirst($class->{'name'})})
    treeNodeSizeOf=treeNodeSizeOf+self%component{ucfirst($class->{'name'})}(i)%sizeOf()
  end do
end if
CODE
    }
    # Iterate over events.
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
event => self%event
do while (associated(event))
  treeNodeSizeOf=treeNodeSizeOf+sizeof(event)+event%nonStaticSizeOf()
  event => event%next
end do
CODE
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    descriptor  => $function,
	    name        => "sizeOf"
	}
	);
}

1;
