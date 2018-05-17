# Contains a Perl module which handles creation and destruction of the treeNode class.

package Galacticus::Build::Components::TreeNodes::Classes;
use strict;
use warnings;
use utf8;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use Text::Template 'fill_in_string';
use List::ExtraUtils;
use Galacticus::Build::Components::Utils;
use Galacticus::Build::Components::DataTypes;

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     treeNodesClasses =>
     {
	 classIteratedFunctions =>
	     [
	      \&Tree_Node_Class_Count,
	      \&Tree_Node_Class_Get
	     ]
     }
    );

sub Tree_Node_Class_Count {
    # Generate a function to return a count of the number of a given component class attached to a tree node.
    my $build    = shift();
    $code::class = shift();
    my $function =
    {
	type        => "integer",
	name        => "treeNode".ucfirst($code::class->{'name'})."Count",
	description => "Returns the number of {\\normalfont \\ttfamily ".$code::class->{'name'}."} components in the node.",
	modules     =>
	    [
	     "Galacticus_Error"
	    ],
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => "treeNode",
		 attributes => [ "intent(inout)", "target" ],
		 variables  => [ "self" ]
	     }
	    ]
    };
    $function->{'content'}  = fill_in_string(<<'CODE', PACKAGE => 'code');
select type (self)
class is (treeNode)
 if (allocated(self%component{ucfirst($class->{'name'})})) then
   select type (component => self%component{ucfirst($class->{'name'})}(1))
   type is (nodeComponent{ucfirst($class->{'name'})})
     treeNode{ucfirst($class->{'name'})}Count=0
   class default
     treeNode{ucfirst($class->{'name'})}Count=size(self%component{ucfirst($class->{'name'})})
   end select
 else
    treeNode{ucfirst($class->{'name'})}Count=0
 end if
class default
 treeNode{ucfirst($class->{'name'})}Count=0
 call Galacticus_Error_Report('treeNode of unknown class'//\{introspection:location\})
end select
CODE
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    descriptor  => $function,
	    name        => $code::class->{'name'}."Count"
	}
	);
}

sub Tree_Node_Class_Get {
    # Generate a function to return a given component class attached to a tree node.
    my $build    = shift();
    $code::class = shift();
    my $function =
    {
	type        => "class(nodeComponent".ucfirst($code::class->{'name'})."), pointer => component",
	name        => "treeNode".ucfirst($code::class->{'name'})."Get",
	description => "Return a {\\normalfont \\ttfamily ".$code::class->{'name'}."} component member of the node. If no {\\normalfont \\ttfamily instance} is specified, return the first instance. If {\\normalfont \\ttfamily autoCreate} is {\\normalfont \\ttfamily true} then create a single instance of the component if none exists in the node.",
	recursive   => 1,
	modules     =>
	    [
	     "Galacticus_Error",
	     "ISO_Varying_String",
	     "String_Handling"
	    ],
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => "treeNode",
		 attributes => [ "intent(inout)", "target" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => [ "intent(in   )", "optional" ],
		 variables  => [ "instance" ]
	     },
	     {
		 intrinsic  => "logical",
		 attributes => [ "intent(in   )", "optional" ],
		 variables  => [ "autoCreate" ]
	     },
	     {
		 intrinsic  => "integer",
		 variables  => [ "instanceActual" ]
	     },
	     {
		 intrinsic  => "logical",
		 variables  => [ "autoCreateActual" ]
	     },
	     {
		 intrinsic  => "type",
		 type       => "varying_string",
		 variables  => [ "message" ]
	     }
	    ]
    };
    $function->{'content'}  = fill_in_string(<<'CODE', PACKAGE => 'code');
instanceActual=1
if (present(instance)) instanceActual=instance
autoCreateActual=.false.
if (present(autoCreate)) autoCreateActual=autoCreate
if (autoCreateActual.and.allocated(self%component{ucfirst($class->{'name'})})) then
   ! If we are allowed to autocreate the component and it still has generic type then deallocate it to force it to be created later.
   if (same_type_as(self%component{ucfirst($class->{'name'})}(1),{ucfirst($class->{'name'})}Class)) deallocate(self%component{ucfirst($class->{'name'})})
end if
if (.not.allocated(self%component{ucfirst($class->{'name'})})) then
  if (autoCreateActual) then
     call self%{$class->{'name'}}Create()
  else
     message='component is not allocated in node '
     message=message//self%index()
     call Galacticus_Error_Report(message//\{introspection:location\})
  end if
end if
component => self%component{ucfirst($class->{'name'})}(instanceActual)
CODE
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    descriptor  => $function,
	    name        => $code::class->{'name'}
	}
	);
}

1;
