# Contains a Perl module which handles creation and destruction of the component classes.

package Galacticus::Build::Components::Classes::CreateDestroy;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Text::Template 'fill_in_string';
use List::ExtraUtils;
use Galacticus::Build::Components::Utils qw(offsetName);
use Galacticus::Build::Components::DataTypes;

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     classesCreateDestroy =>
     {
	 classIteratedFunctions =>
	     [
	      \&Class_Initialization     ,
	      \&Class_Builder            ,
	      \&Class_Finalization       ,
	      \&Class_Create_By_Interrupt,
	      \&Class_Add_Meta_Property
	     ]
     }
    );

sub Class_Initialization {
    # Generate a function to create/initialize component classes.
    my $build    = shift();
    $code::class = shift();
    my $function =
    {
	type        => "void",
	name        => "nodeComponent".ucfirst($code::class->{'name'})."Initialize",
	description => "Initialize a generic {\\normalfont \\ttfamily ".$code::class->{'name'}."} component.",
	modules     =>
	    [
	     "Galacticus_Error"
	    ],
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($code::class->{'name'}),
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     }
	    ]
    };
    $function->{'content'}  = fill_in_string(<<'CODE', PACKAGE => 'code');
!$GLC attributes unused :: self

call Galacticus_Error_Report('can not initialize a generic component'//\{introspection:location\})
CODE
    # Insert a type-binding for this function.
    push(
	@{$build->{'types'}->{"nodeComponent".ucfirst($code::class->{'name'})}->{'boundFunctions'}},
	{
	    type        => "procedure",
	    descriptor  => $function,
	    name        => "initialize", 
	}
	);
}

sub Class_Finalization {
    # Generate a function to finalize component classes.
    my $build    = shift();
    $code::class = shift();
    my $function =
    {
	type        => "void",
	name        => "nodeComponent".ucfirst($code::class->{'name'})."Finalize",
	description => "Finalize a generic {\\normalfont \\ttfamily ".$code::class->{'name'}."} component.",
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($code::class->{'name'}),
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     }
	    ]
    };
    $function->{'content'}  = fill_in_string(<<'CODE', PACKAGE => 'code');
!$GLC attributes unused :: self

! Nothing to do.
CODE
    # Insert a type-binding for this function.
    push(
	@{$build->{'types'}->{"nodeComponent".ucfirst($code::class->{'name'})}->{'boundFunctions'}},
	{
	    type        => "procedure",
	    descriptor  => $function,
	    name        => "destroy", 
	}
	);
}

sub Class_Builder {
    # Generate a function to build component classes from XML definitions.
    my $build    = shift();
    $code::class = shift();
    my $function =
    {
	type        => "void",
	name        => "nodeComponent".ucfirst($code::class->{'name'})."Builder",
	description => "Build a generic {\\normalfont \\ttfamily ".$code::class->{'name'}."} component from a supplied XML definition.",
	modules     =>
	    [
	     "Galacticus_Error",
	     "FoX_DOM, only : node"
	    ],
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($code::class->{'name'}),
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "type",
		 type       => "node",
		 attributes => [ "intent(in   )", "pointer" ],
		 variables  => [ "componentDefinition" ]
	     }
	    ]
	};
    $function->{'content'}  = fill_in_string(<<'CODE', PACKAGE => 'code');
!$GLC attributes unused :: self, componentDefinition

call Galacticus_Error_Report('can not build a generic component'//\{introspection:location\})
CODE
    # Insert a type-binding for this function.
    push(
	@{$build->{'types'}->{"nodeComponent".ucfirst($code::class->{'name'})}->{'boundFunctions'}},
	{
	    type        => "procedure",
	    descriptor  => $function,
	    name        => "builder", 
	}
	);
}

sub Class_Create_By_Interrupt {
    # Generate a function to create a component of given class in a tree node via an ODE solver interrupt.
    my $build    = shift();
    $code::class = shift();
    # A function is required for this class only if at least one member has at least one property with the "createIfNeeded" attribute.
    my $functionRequired = 0;
    foreach my $member ( @{$code::class->{'members'}} ) {
	$functionRequired = 1
	    if ( grep {$_->{'attributes'}->{'createIfNeeded'}} &List::ExtraUtils::hashList($member->{'properties'}->{'property'}) );	
    }
    return
	unless ( $functionRequired );
    # Create the function.
    my $function =
    {
	type        => "void",
	name        => $code::class->{'name'}."CreateByInterrupt",
	description => "Create the {\\normalfont \\ttfamily ".$code::class->{'name'}."} component of {\\normalfont \\ttfamily self} via an interrupt.",
	variables   =>
	    [
	     {
		 intrinsic  => "type",
		 type       => "treeNode",
		 attributes => [ "target", "intent(inout)" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($code::class->{'name'}),
		 attributes => [ "pointer" ],
		 variables  => [ $code::class->{'name'} ]
	     }
	    ]
    };
    $function->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
{$class->{'name'}} => self%{$class->{'name'}}(autoCreate=.true.)
CODE
    # Iterate over class, and call custom create routines if necessary.
    if ( grep {exists($_->{'createFunction'})}  @{$code::class->{'members'}} ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
select type ({$class->{'name'}})
CODE
	foreach $code::member ( grep {exists($_->{'createFunction'})}  @{$code::class->{'members'}} ) {
	    $code::createFunction = 
		        $code::member->{'createFunction'}->{'isDeferred'} 
	        ? 
		$code::class->{'name'}.ucfirst($code::member->{'name'})."CreateFunction" 
		: 
		(
		 exists($code::member->{'createFunction'}->{'content'   }) 
		 ?
		        $code::member->{'createFunction'}->{'content'   } 
		 :
		        $code::member->{'createFunction'}
		);
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
type is (nodeComponent{ucfirst($code::class->{'name'}).ucfirst($code::member->{'name'})})
   call {$createFunction}({$class->{'name'}})
CODE
	}
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
end select
CODE
    }
    # Insert into the functions list.
    push(@{$build->{'functions'}},$function);    
}

sub Class_Add_Meta_Property {
    # Generate a function to add meta-properties to component classes.
    my $build    = shift();
    $code::class = shift();
    my $function =
    {
	type        => "integer",
	name        => "nodeComponent".ucfirst($code::class->{'name'})."AddMetaProperty",
	description => "Add a meta-property to the generic {\\normalfont \\ttfamily ".$code::class->{'name'}."} component.",
	modules     =>
	    [
	     "ISO_Varying_String"
	    ],
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($code::class->{'name'}),
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "type",
		 type       => "varying_string",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "label" ]
	     },
	     {
		 intrinsic  => "character",
		 type       => "len=*",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "name" ]
	     },
	     {
		 intrinsic  => "logical",
		 attributes => [ "intent(in   )", "optional" ],
		 variables  => [ "isEvolvable" ]
	     },
	     {
		 intrinsic  => "type",
		 type       => "varying_string",
		 attributes => [ "allocatable", "dimension(:)" ],
		 variables  => [ "labelsTmp", "namesTmp" ]
	     },
	     {
		 intrinsic  => "logical",
		 attributes => [ "allocatable", "dimension(:)" ],
		 variables  => [ "evolvableTmp" ]
	     },
	     {
		 intrinsic  => "logical",
		 variables  => [ "found" ]
	     }
	    ]
    };
    if ( grep {$code::class->{'name'} eq $_} @{$build->{'componentClassListActive'}} ) {
	$function->{'content'}  = fill_in_string(<<'CODE', PACKAGE => 'code');
!$GLC attributes unused :: self

!$omp critical ({class->{'name'}}MetaPropertyUpdate)
found=.false.
if (allocated({$class->{'name'}}MetaPropertyLabels)) then
 do nodeComponent{ucfirst($class->{'name'})}AddMetaProperty=1,size({$class->{'name'}}MetaPropertyLabels)
  if ({$class->{'name'}}MetaPropertyLabels(nodeComponent{ucfirst($class->{'name'})}AddMetaProperty) == label) then
   found=.true.
   exit
  end if
 end do
 if (.not.found) then
  call move_alloc({$class->{'name'}}MetaPropertyLabels   ,   labelsTmp)
  call move_alloc({$class->{'name'}}MetaPropertyNames    ,    namesTmp)
  call move_alloc({$class->{'name'}}MetaPropertyEvolvable,evolvableTmp)
  allocate({$class->{'name'}}MetaPropertyLabels   (size(   labelsTmp)+1))
  allocate({$class->{'name'}}MetaPropertyNames    (size(    namesTmp)+1))
  allocate({$class->{'name'}}MetaPropertyEvolvable(size(evolvableTmp)+1))
  {$class->{'name'}}MetaPropertyLabels   (1:size(   labelsTmp))=   labelsTmp
  {$class->{'name'}}MetaPropertyNames    (1:size(    namesTmp))=    namesTmp
  {$class->{'name'}}MetaPropertyEvolvable(1:size(evolvableTmp))=evolvableTmp
  deallocate(   labelsTmp)
  deallocate(    namesTmp)
  deallocate(evolvableTmp)
 end if
else
 allocate({$class->{'name'}}MetaPropertyLabels   (                1))
 allocate({$class->{'name'}}MetaPropertyNames    (                1))
 allocate({$class->{'name'}}MetaPropertyEvolvable(                1))
end if
if (.not.found) then
 nodeComponent{ucfirst($class->{'name'})}AddMetaProperty=size({$class->{'name'}}MetaPropertyLabels)
 {$class->{'name'}}MetaPropertyLabels   (nodeComponent{ucfirst($class->{'name'})}AddMetaProperty)=label
 {$class->{'name'}}MetaPropertyNames    (nodeComponent{ucfirst($class->{'name'})}AddMetaProperty)=name
 if (present(isEvolvable)) then
  {$class->{'name'}}MetaPropertyEvolvable(nodeComponent{ucfirst($class->{'name'})}AddMetaProperty)=isEvolvable
 else
  {$class->{'name'}}MetaPropertyEvolvable(nodeComponent{ucfirst($class->{'name'})}AddMetaProperty)=.false.
 end if
 {$class->{'name'}}MetaPropertyCount={$class->{'name'}}MetaPropertyCount+1
 if ({$class->{'name'}}MetaPropertyEvolvable(nodeComponent{ucfirst($class->{'name'})}AddMetaProperty)) {$class->{'name'}}MetaPropertyEvolvableCount={$class->{'name'}}MetaPropertyEvolvableCount+1
 propertyNameLengthMax=max(len(name),propertyNameLengthMax) 
end if
!$omp end critical ({class->{'name'}}MetaPropertyUpdate)
CODE
    }
    # Insert a type-binding for this function.
    push(
	@{$build->{'types'}->{"nodeComponent".ucfirst($code::class->{'name'})}->{'boundFunctions'}},
	{
	    type        => "procedure",
	    descriptor  => $function,
	    name        => "addMetaProperty", 
	}
	);
}

1;
