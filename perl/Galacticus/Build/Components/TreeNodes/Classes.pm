# Contains a Perl module which handles creation and destruction of the treeNode class.

package Galacticus::Build::Components::TreeNodes::Classes;
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
	     "Error"
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
    if ( grep {$code::class->{'name'} eq $_} @{$build->{'componentClassListActive'}} ) {
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
 call Error_Report('treeNode of unknown class'//\{introspection:location\})
end select
CODE
    } else {
	$function->{'content'}  = fill_in_string(<<'CODE', PACKAGE => 'code');
!$GLC attributes unused :: self
treeNode{ucfirst($class->{'name'})}Count=0
call Error_Report('Galacticus was not compiled with support for this class'//\{introspection:location\})
CODE
    }
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
	     "Error"
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
	     }
	    ]
    };
    if ( grep {$code::class->{'name'} eq $_} @{$build->{'componentClassListActive'}} ) {
	$function->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
if (.not.present(autoCreate).and..not.present(instance)) then
   component => self%component{ucfirst($class->{'name'})}(1)
else
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
        call nodeComponentGetError('{$class->{'name'}}',self%index())
     end if
   end if
   component => self%component{ucfirst($class->{'name'})}(instanceActual)
end if
CODE
    } else {
	$function->{'content'}  = fill_in_string(<<'CODE', PACKAGE => 'code');
!$GLC attributes unused :: self, instance, instanceActual
autoCreateActual=.false.
if (present(autoCreate)) autoCreateActual=autoCreate
if (autoCreateActual) then
 ! Support for this component was not compiled, so we can not create it.
 component => null()
 call Error_Report('Galacticus was compiled without support for this class'//\{introspection:location\})
else
 ! Support for this component was not compiled, return the default of the class - and trust that the user knows what they are doing.
 component => default{ucfirst($class->{'name'})}Component
end if
CODE
    }
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
