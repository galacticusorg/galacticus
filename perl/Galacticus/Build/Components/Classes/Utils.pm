# Contains a Perl module which provides utility functions for component classes.

package Galacticus::Build::Components::Classes::Utils;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Sub::Identify ':all';
use Text::Template 'fill_in_string';
use List::ExtraUtils;
use Data::Dumper;
use Galacticus::Build::Components::Utils;
use Galacticus::Build::Components::DataTypes;

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     classUtils =>
     {
	 functions =>
	     [
	      \&Class_Move             ,
	      \&Class_Remove           ,
	      \&Class_Function_Iterator
	     ]
     }
    );

sub Class_Move {
    # Generate functions to move component classes.
    my $build = shift();
    # Iterate over component classes.
    foreach $code::class ( &List::ExtraUtils::hashList($build->{'componentClasses'}) ) {
	my $classTypeName = "nodeComponent".ucfirst($code::class->{'name'});
	my $function =
	{
	    type        => "void",
	    name        => $classTypeName."Move",
	    description => "Move instances of the {\\normalfont \\ttfamily ".$code::class->{'name'}."} component, from one node to another.",
	    modules     =>
		[
		 "Galacticus_Error"
		],
	    variables   =>
		[
		 {
		     intrinsic  => "class",
		     type       => "treeNode",
		     attributes => [ "intent(inout)" ],
		     variables  => [ "self" ]
		 },
		 {
		     intrinsic  => "type",
		     type       => "treeNode",
		     attributes => [ "intent(inout)", "target" ],
		     variables  => [ "targetNode" ]
		 },
		 {
		     intrinsic  => "logical",
		     attributes => [ "intent(in   )", "optional" ],
		     variables  => [ "overwrite" ]
		 },
		 {
		     intrinsic  => "integer",
		     variables  => [ "instanceCount", "targetCount", "i" ]
		 },
		 {
		     intrinsic  => "class",
		     type       => $classTypeName,
		     attributes => [ "allocatable, dimension(:)" ],
		     variables  => [ "instancesTemporary" ]
		 },
		 {
		     intrinsic  => "logical",
		     variables  => [ "overwrite_" ]
		 }
		]
	};
	if ( grep {$code::class->{'name'} eq $_} @{$build->{'componentClassListActive'}} ) {
	    $function->{'content'}  = fill_in_string(<<'CODE', PACKAGE => 'code');
overwrite_=.false.
if (present(overwrite)) overwrite_=overwrite
instanceCount=self      %{$class->{'name'}}count()
targetCount  =targetNode%{$class->{'name'}}count()
if (overwrite_ .and. targetCount > 0) then
  do i=1,targetCount
    call targetNode%component{ucfirst($class->{'name'})}(i)%destroy()
  end do 
  targetCount=0
  deallocate(targetNode%component{ucfirst($class->{'name'})})
  allocate(targetNode%component{ucfirst($class->{'name'})}(1))
end if	
if (instanceCount == 0) return
if (targetCount == 0) then
  deallocate(targetNode%component{ucfirst($class->{'name'})})
  call Move_Alloc(self%component{ucfirst($class->{'name'})},targetNode%component{ucfirst($class->{'name'})})
else
  ! Multiple instances, so remove the specified instance.
  allocate(instancesTemporary(instanceCount+targetCount),source=self%component{ucfirst($class->{'name'})}(1))
CODE
	    foreach $code::member ( @{$code::class->{'members'}} ) {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
  select type (from => targetNode%component{ucfirst($class->{'name'})})
  type is (nodeComponent{ucfirst($class->{'name'}).ucfirst($member->{'name'})})
    select type (to => instancesTemporary)
    type is (nodeComponent{ucfirst($class->{'name'}).ucfirst($member->{'name'})})
      to(1:targetCount)=from
    end select
  end select
CODE
	    }
	    foreach $code::member ( @{$code::class->{'members'}} ) {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
   select type (from => self%component{ucfirst($class->{'name'})})
   type is (nodeComponent{ucfirst($class->{'name'}).ucfirst($member->{'name'})})
     select type (to => instancesTemporary)
     type is (nodeComponent{ucfirst($class->{'name'}).ucfirst($member->{'name'})})
       to(targetCount+1:targetCount+instanceCount)=from
     end select
   end select
CODE
	    }
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
  call targetNode%{$class->{'name'}}Destroy()
  call self      %{$class->{'name'}}Destroy()
  call Move_Alloc(instancesTemporary,targetNode%component{ucfirst($class->{'name'})})
  allocate(self%component{ucfirst($class->{'name'})}(1))
end if
do i=1,size(targetNode%component{ucfirst($class->{'name'})})
   targetNode%component{ucfirst($class->{'name'})}(i)%hostNode => targetNode
end do
CODE
	} else {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
!$GLC attributes unused :: self, targetNode, overwrite, instanceCount, targetCount, i, instancesTemporary, overwrite_
call Galacticus_Error_Report('Galacticus was not compiled with support for this class'//\{introspection:location\})
CODE
	}
	# Insert a type-binding for this function into the treeNode type.
	push(
	    @{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	    {
		type        => "procedure", 
		descriptor  => $function,
		name        => $code::class->{'name'}."Move", 
		returnType  => "\\void", 
		arguments   => "\\textcolor{red}{\\textless type(treeNode)\\textgreater} targetNode\\arginout"
	    }
	    );
    }
}

sub Class_Remove {
    # Generate functions to remove component classes.
    my $build = shift();
    # Iterate over component classes.
    foreach $code::class ( &List::ExtraUtils::hashList($build->{'componentClasses'}) ) {
	my $classTypeName = "nodeComponent".ucfirst($code::class->{'name'});
	my $function =
	{
	    type        => "void",
	    name        => $classTypeName."Remove",
	    description => "Remove an instance of the {\\normalfont \\ttfamily ".$code::class->{'name'}."} component from a node.",
	    modules     =>
		[
		 "Galacticus_Error"
		],
	    variables   =>
		[
		 {
		     intrinsic  => "class",
		     type       => "treeNode",
		     attributes => [ "intent(inout)" ],
		     variables  => [ "self" ]
		 },
		 {
		     intrinsic  => "integer",
		     attributes => [ "intent(in   )" ],
		     variables  => [ "instance" ]
		 },
		 {
		     intrinsic  => "integer",
		     variables  => [ "instanceCount" ]
		 },
		 {
		     intrinsic  => "class",
		     type       => $classTypeName,
		     attributes => [ "allocatable, dimension(:)" ],
		     variables  => [ "instancesTemporary" ]
		 }
		]
	};    
	if ( grep {$code::class->{'name'} eq $_} @{$build->{'componentClassListActive'}} ) {
	    $function->{'content'}  = fill_in_string(<<'CODE', PACKAGE => 'code');
instanceCount=self%{$class->{'name'}}count()
if (instance < 1 .or. instance > instanceCount) call Galacticus_Error_Report('instance out of range'//\{introspection:location\})
call self%component{ucfirst($class->{'name'})}(instance)%destroy()
if (instanceCount == 1) then
  ! Only one instance of this component. Deallocate it and reallocate with generic type.
  deallocate(self%component{ucfirst($class->{'name'})})
  allocate(self%component{ucfirst($class->{'name'})}(1))
else
  ! Multiple instances, so remove the specified instance.
  allocate(instancesTemporary(instanceCount-1),source=self%component{ucfirst($class->{'name'})}(1))
CODE
	    foreach $code::member ( @{$code::class->{'members'}} ) {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
  select type (from => self%component{ucfirst($class->{'name'})})
  type is (nodeComponent{ucfirst($class->{'name'}).ucfirst($member->{'name'})})
    select type (to => instancesTemporary)
    type is (nodeComponent{ucfirst($class->{'name'}).ucfirst($member->{'name'})})
      if (instance >             1) to(       1:instance     -1)=from(         1:instance     -1)
      if (instance < instanceCount) to(instance:instanceCount-1)=from(instance+1:instanceCount  )
    end select
  end select
CODE
	    }
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
  deallocate(self%component{ucfirst($class->{'name'})})
  call Move_Alloc(instancesTemporary,self%component{ucfirst($class->{'name'})})
end if
CODE
	} else {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
!$GLC attributes unused :: self, instance, instanceCount, instancesTemporary
call Galacticus_Error_Report('Galacticus was not compiled with support for this class'//\{introspection:location\})
CODE
	}
	# Insert a type-binding for this function into the treeNode type.
	push(
	    @{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	    {
		type        => "procedure", 
		descriptor  => $function,
		name        => $code::class->{'name'}."Remove"
	    }
	    );
    }
}

sub Class_Function_Iterator {
    # Iterates over component classes and calls registered functions for each class.
    my $build = shift();
    my @hooks = &List::ExtraUtils::hashList(\%Galacticus::Build::Component::Utils::componentUtils, keyAs => 'name');
    foreach my $hook ( @hooks ) {
	if ( exists($hook->{'classIteratedFunctions'}) ) {
	    my @functions = &List::ExtraUtils::as_array($hook->{'classIteratedFunctions'});
	    foreach my $function ( @functions ) {
		print "         --> ".$hook->{'name'}.(scalar(@functions) > 1 ? " {".sub_name($function)."}" : "")."\n";
		# Iterate over classes.
		foreach my $class ( &List::ExtraUtils::hashList($build->{'componentClasses'}) ) {
		    &{$function}($build,$class);
		}
	    }
	}
    }
}

1;
