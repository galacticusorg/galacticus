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
use Galacticus::Build::Components::Classes::MetaProperties;

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
	      \&Class_Add_Meta_Property  ,
	      \&Class_Count_Meta_Property,
	      \&Class_Name_Meta_Property  
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
	     "Error"
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

call Error_Report('can not initialize a generic component'//\{introspection:location\})
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
	     "Error",
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

call Error_Report('can not build a generic component'//\{introspection:location\})
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
    foreach my $metaPropertyType ( @Galacticus::Build::Components::Classes::MetaProperties::metaPropertyTypes ) {
	$code::label    = $metaPropertyType->{'label'};
	$code::rank     = $metaPropertyType->{'rank' };
	$code::prefix   = ucfirst($metaPropertyType->{'label'})."Rank".$metaPropertyType->{'rank' };
	my $isEvolvable = $metaPropertyType->{'label'} eq "float" && $metaPropertyType->{'rank'} == 0;
	my @options     = ( "isCreator"  );
	my @tmps        = ( "creatorTmp" );
	if ( $isEvolvable ) {
	    unshift(@options,"isEvolvable" );
	    unshift(@tmps   ,"evolvableTmp");
	}
	my $function =
	{
	    type        => "integer",
	    name        => "nodeComponent".ucfirst($code::class->{'name'})."Add".ucfirst($metaPropertyType->{'label'})."Rank".$metaPropertyType->{'rank'}."MetaProperty",
	    description => "Add a rank-".$metaPropertyType->{'rank'}." ".$metaPropertyType->{'label'}."meta-property to the generic {\\normalfont \\ttfamily ".$code::class->{'name'}."} component.",
	    modules     =>
		[
		 "ISO_Varying_String",
		 "Error"
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
		     variables  => \@options
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
		     variables  => \@tmps
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

!$omp critical ({$class->{'name'}.$prefix}MetaPropertyUpdate)
found=.false.
if (allocated({$class->{'name'}.$prefix}MetaPropertyLabels)) then
 do nodeComponent{ucfirst($class->{'name'})}Add{$prefix}MetaProperty=1,size({$class->{'name'}.$prefix}MetaPropertyLabels)
  if ({$class->{'name'}.$prefix}MetaPropertyLabels(nodeComponent{ucfirst($class->{'name'})}Add{$prefix}MetaProperty) == label) then
   found=.true.
   exit
  end if
 end do
 if (.not.found) then
  call move_alloc({$class->{'name'}.$prefix}MetaPropertyLabels , labelsTmp)
  call move_alloc({$class->{'name'}.$prefix}MetaPropertyNames  ,  namesTmp)
  call move_alloc({$class->{'name'}.$prefix}MetaPropertyCreator,creatorTmp)
  allocate({$class->{'name'}.$prefix}MetaPropertyLabels (size( labelsTmp)+1))
  allocate({$class->{'name'}.$prefix}MetaPropertyNames  (size(  namesTmp)+1))
  allocate({$class->{'name'}.$prefix}MetaPropertyCreator(size(creatorTmp)+1))
  {$class->{'name'}.$prefix}MetaPropertyLabels (1:size( labelsTmp))= labelsTmp
  {$class->{'name'}.$prefix}MetaPropertyNames  (1:size(  namesTmp))=  namesTmp
  {$class->{'name'}.$prefix}MetaPropertyCreator(1:size(creatorTmp))=creatorTmp
  deallocate( labelsTmp)
  deallocate(  namesTmp)
  deallocate(creatorTmp)
CODE
	    if ( $isEvolvable ) {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
  call move_alloc({$class->{'name'}.$prefix}MetaPropertyEvolvable,evolvableTmp)
  allocate({$class->{'name'}.$prefix}MetaPropertyEvolvable(size(evolvableTmp)+1))
  {$class->{'name'}.$prefix}MetaPropertyEvolvable(1:size(evolvableTmp))=evolvableTmp
  deallocate(evolvableTmp)
CODE
	    }
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
 end if
else
 allocate({$class->{'name'}.$prefix}MetaPropertyLabels (1))
 allocate({$class->{'name'}.$prefix}MetaPropertyNames  (1))
 allocate({$class->{'name'}.$prefix}MetaPropertyCreator(1))
CODE
	    if ( $isEvolvable ) {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
 allocate({$class->{'name'}.$prefix}MetaPropertyEvolvable(1))
CODE
	    }
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
end if
CODE
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if (.not.found) then
 nodeComponent{ucfirst($class->{'name'})}Add{$prefix}MetaProperty=size({$class->{'name'}.$prefix}MetaPropertyLabels)
 {$class->{'name'}.$prefix}MetaPropertyLabels   (nodeComponent{ucfirst($class->{'name'})}Add{$prefix}MetaProperty)=label
 {$class->{'name'}.$prefix}MetaPropertyNames    (nodeComponent{ucfirst($class->{'name'})}Add{$prefix}MetaProperty)=name
CODE
	    if ( $isEvolvable ) {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
 if (present(isEvolvable)) then
  {$class->{'name'}.$prefix}MetaPropertyEvolvable(nodeComponent{ucfirst($class->{'name'})}Add{$prefix}MetaProperty)=isEvolvable
 else
  {$class->{'name'}.$prefix}MetaPropertyEvolvable(nodeComponent{ucfirst($class->{'name'})}Add{$prefix}MetaProperty)=.false.
 end if
CODE
	    }
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
 if (present(isCreator  )) then
  {$class->{'name'}.$prefix}MetaPropertyCreator  (nodeComponent{ucfirst($class->{'name'})}Add{$prefix}MetaProperty)=isCreator
 else
  {$class->{'name'}.$prefix}MetaPropertyCreator  (nodeComponent{ucfirst($class->{'name'})}Add{$prefix}MetaProperty)=.false.
 end if
 {$class->{'name'}.$prefix}MetaPropertyCount={$class->{'name'}.$prefix}MetaPropertyCount+1
CODE
	    if ( $isEvolvable ) {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
 if ({$class->{'name'}.$prefix}MetaPropertyEvolvable(nodeComponent{ucfirst($class->{'name'})}Add{$prefix}MetaProperty)) {$class->{'name'}.$prefix}MetaPropertyEvolvableCount={$class->{'name'}.$prefix}MetaPropertyEvolvableCount+1
CODE
	    }
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
 propertyNameLengthMax=max(len(name),propertyNameLengthMax) 
else
 if (present(isCreator)) then
  if (isCreator) {$class->{'name'}.$prefix}MetaPropertyCreator(nodeComponent{ucfirst($class->{'name'})}Add{$prefix}MetaProperty)=.true.
 end if
CODE
	    if ( $isEvolvable ) {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
 if (present(isEvolvable)) then
  if ({$class->{'name'}.$prefix}MetaPropertyEvolvable(nodeComponent{ucfirst($class->{'name'})}Add{$prefix}MetaProperty) .neqv. isEvolvable) call Error_Report('inconsistent evolvability for meta-property'//\{introspection:location\})
 else
  if ({$class->{'name'}.$prefix}MetaPropertyEvolvable(nodeComponent{ucfirst($class->{'name'})}Add{$prefix}MetaProperty)                   ) call Error_Report('inconsistent evolvability for meta-property'//\{introspection:location\})
 end if
CODE
	    }
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
end if
!$omp end critical ({$class->{'name'}.$prefix}MetaPropertyUpdate)
CODE
	}
	# Insert a type-binding for this function.
	push(
	    @{$build->{'types'}->{"nodeComponent".ucfirst($code::class->{'name'})}->{'boundFunctions'}},
	    {
		type        => "procedure",
		descriptor  => $function,
		name        => "add".ucfirst($metaPropertyType->{'label'})."Rank".$metaPropertyType->{'rank'}."MetaProperty", 
	    }
	    );
    }
}

sub Class_Count_Meta_Property {
    # Generate a function to return a count of meta-properties for a component classes.
    my $build    = shift();
    $code::class = shift();
    foreach my $metaPropertyType ( @Galacticus::Build::Components::Classes::MetaProperties::metaPropertyTypes ) {
	$code::label    = $metaPropertyType->{'label'};
	$code::rank     = $metaPropertyType->{'rank' };
	$code::prefix   = ucfirst($metaPropertyType->{'label'})."Rank".$metaPropertyType->{'rank' };
	my $function =
	{
	    type        => "integer => countMetaProperties",
	    name        => "component".ucfirst($code::class->{'name'})."Count".ucfirst($metaPropertyType->{'label'})."Rank".$metaPropertyType->{'rank'}."MetaProperties",
	    description => "Return the number of rank-".$metaPropertyType->{'rank'}." ".$metaPropertyType->{'label'}."meta-properties associated with the generic {\\normalfont \\ttfamily ".$code::class->{'name'}."} component.",
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
	if ( grep {$code::class->{'name'} eq $_} @{$build->{'componentClassListActive'}} ) {
	    $function->{'content'}  = fill_in_string(<<'CODE', PACKAGE => 'code');
!$GLC attributes unused :: self

!$omp critical ({$class->{'name'}.$prefix}MetaPropertyUpdate)
if (allocated({$class->{'name'}.$prefix}MetaPropertyNames)) then
 countMetaProperties=size({$class->{'name'}.$prefix}MetaPropertyNames)
else
 countMetaProperties=0
end if
!$omp end critical ({$class->{'name'}.$prefix}MetaPropertyUpdate)
CODE
	}
	# Insert a type-binding for this function.
	push(
	    @{$build->{'types'}->{"nodeComponent".ucfirst($code::class->{'name'})}->{'boundFunctions'}},
	    {
		type        => "procedure",
		descriptor  => $function,
		name        => "count".ucfirst($metaPropertyType->{'label'})."Rank".$metaPropertyType->{'rank'}."MetaProperties", 
	    }
	    );
    }
}

sub Class_Name_Meta_Property {
    # Generate a function to return the name of the indexed meta-property for a component classes.
    my $build    = shift();
    $code::class = shift();
    foreach my $metaPropertyType ( @Galacticus::Build::Components::Classes::MetaProperties::metaPropertyTypes ) {
	$code::label    = $metaPropertyType->{'label'};
	$code::rank     = $metaPropertyType->{'rank' };
	$code::prefix   = ucfirst($metaPropertyType->{'label'})."Rank".$metaPropertyType->{'rank' };
	my $function =
	{
	    type        => "type(varying_string) => nameMetaProperty",
	    name        => "component".ucfirst($code::class->{'name'})."Name".ucfirst($metaPropertyType->{'label'})."Rank".$metaPropertyType->{'rank'}."MetaProperty",
	    description => "Return the name of the indexed of rank-".$metaPropertyType->{'rank'}." ".$metaPropertyType->{'label'}." meta-property associated with the generic {\\normalfont \\ttfamily ".$code::class->{'name'}."} component.",
	    modules     => [ "ISO_Varying_String" ],
	    variables   =>
		[
		 {
		     intrinsic  => "class",
		     type       => "nodeComponent".ucfirst($code::class->{'name'}),
		     attributes => [ "intent(inout)" ],
		     variables  => [ "self" ]
		 },
		 {
		     intrinsic  => "integer",
		     attributes => [ "intent(in   )" ],
		     variables  => [ "index" ]
		 }
		]
	};
	if ( grep {$code::class->{'name'} eq $_} @{$build->{'componentClassListActive'}} ) {
	    push(@{$function->{'modules'}},"Error");
	    $function->{'content'}  = fill_in_string(<<'CODE', PACKAGE => 'code');
!$GLC attributes unused :: self

!$omp critical ({$class->{'name'}.$prefix}MetaPropertyUpdate)
if (index > 0 .and. index <= size({$class->{'name'}.$prefix}MetaPropertyNames)) then
 nameMetaProperty={$class->{'name'}.$prefix}MetaPropertyNames(index)
else
 nameMetaProperty=var_str('')
 call Error_Report('meta-property index is out of range'//\{introspection:location\})
end if
!$omp end critical ({$class->{'name'}.$prefix}MetaPropertyUpdate)
CODE
	}
	# Insert a type-binding for this function.
	push(
	    @{$build->{'types'}->{"nodeComponent".ucfirst($code::class->{'name'})}->{'boundFunctions'}},
	    {
		type        => "procedure",
		descriptor  => $function,
		name        => "name".ucfirst($metaPropertyType->{'label'})."Rank".$metaPropertyType->{'rank'}."MetaProperty", 
	    }
	    );
    }
}

1;
