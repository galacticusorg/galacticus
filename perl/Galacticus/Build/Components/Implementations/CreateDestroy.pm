# Contains a Perl module which handles creation and destruction of the component implementations.

package Galacticus::Build::Components::Implementations::CreateDestroy;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Text::Template 'fill_in_string';
use List::ExtraUtils;
use Galacticus::Build::Components::Utils qw(&isIntrinsic %intrinsicNulls);
use Galacticus::Build::Components::DataTypes;
use Galacticus::Build::Components::Classes::MetaProperties;
use Data::Dumper;

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     implementationsCreateDestroy =>
     {
	 implementationIteratedFunctions =>
	     [
	      \&Implementation_Creation           ,
	      \&Implementation_Builder            ,
	      \&Implementation_Finalization       ,
	      \&Implementation_Deferred_Create_Set
	     ]
     }
    );

sub Implementation_Creation {
    # Generate a function to create/initialize component classes.
    my $build     = shift();
    $code::class  = shift();
    $code::member = shift();
    $code::implementationTypeName = "nodeComponent".ucfirst($code::class->{'name'}).ucfirst($code::member->{'name'});
    my $function  =
    {
	type        => "void",
	name        => $code::implementationTypeName."Initialize",
	description => "Initialize a {\\normalfont \\ttfamily ".$code::member->{'name'}."} member of the {\\normalfont \\ttfamily ".$code::class->{'name'}."} component.",
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => $code::implementationTypeName,
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "integer",
		 variables  => [ "i" ]
	     }
	    ]
    };
    my %modules;
    # Add variables for any other components required for initialization.
    my $componentInitialization;
    my %requiredComponents;
    foreach my $property ( &List::ExtraUtils::hashList($code::member->{'properties'}->{'property'}) ) {
	next
	    if ( $property->{'attributes'}->{'isVirtual'} || ! exists($property->{'classDefault'}->{'code'}) );
	my $classDefault = $property->{'classDefault'}->{'code'};
	while ( $classDefault =~ s/self([a-zA-Z]+)\s*%// ) {
	    my $componentName = $1;
	    unless ( exists($requiredComponents{lc($componentName)}) ) {
		$componentInitialization .= "self".$componentName." => self%hostNode%".lc($componentName)."()\n";
		push(
		    @{$function->{'variables'}},
		    {
			intrinsic  => "class",
			type       => "nodeComponent".ucfirst($componentName),
			attributes => [ "pointer" ],
			variables  => [ "self".ucfirst($componentName) ]
		    }
		    );
	    }
	    $requiredComponents{lc($componentName)} = 1;
	}
    }
    # Insert any required modules.
    push(
	@{$function->{'modules'}},
	map
	{
	    exists($_->{'classDefault'}->{'modules'}) 
		?
		@{$_->{'classDefault'}->{'modules'}} 
	    :
		()
	} 
	&List::ExtraUtils::hashList($code::member->{'properties'}->{'property'})
	);
    # Mark self as unused if necessary.
    $function->{'content'} = "";
    unless (
	defined($componentInitialization)
	||
	exists($code::member->{'extends'})
	|| 
	grep {! $_->{'attributes'}->{'isVirtual'}} &List::ExtraUtils::hashList($code::member->{'properties'}->{'property'}) 
	) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
!$GLC attributes unused :: self
CODE
    }
    # Initialize the parent class.
    if ( exists($code::member->{'extends'}) ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
call self%nodeComponent{ucfirst($member->{'extends'}->{'class'}).ucfirst($member->{'extends'}->{'name'})}%initialize()
CODE
    }
    # Insert component initialization.
    $function->{'content'} .= $componentInitialization
	if ( defined($componentInitialization) );
    # Iterate over properties.
    foreach $code::property ( &List::ExtraUtils::hashList($code::member->{'properties'}->{'property'}) ) {
	# Skip this property if it is virtual.
	next
	    if ( $code::property->{'attributes'}->{'isVirtual'} );
	# Set to a class default value if available.
	if ( exists($code::property->{'classDefault'}->{'code'}) ) {
	    if ( exists($code::property->{'classDefault'}->{'count'}) ) {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
allocate(self%{$property->{'name'}}Data({$property->{'classDefault'}->{'count'}}))
CODE
	    }
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
self%{$property->{'name'}}Data={$property->{'classDefault'}->{'code'}}
CODE
	} else {
	    # Set to appropriate null value.
	    if (
		(
		 $code::property->{'data'}->{'type'} eq "double"
		 ||
		 $code::property->{'data'}->{'type'} eq "longInteger"
		)
		&&
		$code::property->{'data'}->{'rank'} == 1 
		) {
	    	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
allocate(self%{$property->{'name'}}Data({join(",","0" x $property->{'data'}->{'rank'})}))
CODE
	    }
	    if (&isIntrinsic($code::property->{'data'}->{'type'})) {
		%code::nullValues = %intrinsicNulls;
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
self%{$property->{'name'}}Data={$nullValues{$property->{'data'}->{'type'}}}
CODE
	    } else {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
call self%{$property->{'name'}}Data%reset()
CODE
	    }
	}
    }
    # Allocate and initialize meta-properties.
    if ( grep {$code::class->{'name'} eq $_} @{$build->{'componentClassListActive'}} ) {
	foreach my $metaPropertyType ( @Galacticus::Build::Components::Classes::MetaProperties::metaPropertyTypes ) {
	    $code::label    = $metaPropertyType->{'label'};
	    $code::rank     = $metaPropertyType->{'rank' };
	    $code::prefix   = ucfirst($metaPropertyType->{'label'})."Rank".$metaPropertyType->{'rank'};
	    if ( $metaPropertyType->{'rank'} == 0 ) {
		$code::initializer = " self%".$code::prefix."MetaProperties=";
		if      ( $metaPropertyType->{'intrinsic'} eq "double precision" ) {
		    $code::initializer .= "0.0";
		} elsif ( $metaPropertyType->{'intrinsic'} eq "integer"          ) {
		    $code::initializer .= "0";
		} else {
		    die("Galacticus::Build::Components::Implementations::CreateDestroy::Implementation_Creation: unknown meta-property type");
		}
		if      ( exists($metaPropertyType->{'type'})                   ) {
		    $code::initializer .= "_".$metaPropertyType->{'type'};
		} elsif ( $metaPropertyType->{'intrinsic'} eq "double precision" ) {
		    $code::initializer .= "d0";
		}
	    } else {
		$code::initializer = fill_in_string(<<'CODE', PACKAGE => 'code');
 do i=1,size({$class->{'name'}.$prefix}MetaPropertyNames)
  allocate(self%{$prefix}MetaProperties(i)%values(0))
 end do
CODE
	    }	    
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if (allocated({$class->{'name'}.$prefix}MetaPropertyNames).and..not.allocated(self%{$prefix}MetaProperties)) then
 allocate(self%{$prefix}MetaProperties(size({$class->{'name'}.$prefix}MetaPropertyNames)))
{$initializer}
end if
CODE
	}
    }
    # Add required modules to function.
    push(@{$function->{'modules'}},keys(%modules))
	if ( %modules );
    # Insert a type-binding for this function.
    push(
	@{$build->{'types'}->{$code::implementationTypeName}->{'boundFunctions'}},
	{
	    type        => "procedure",
	    descriptor  => $function,
	    name        => "initialize", 
	}
	);
}

sub Implementation_Finalization {
    # Generate a function to build component implementations from XML definitions.
    my $build     = shift();
    $code::class  = shift();
    $code::member = shift();
    # Define the function.
    $code::implementationTypeName = "nodeComponent".ucfirst($code::class->{'name'}).ucfirst($code::member->{'name'});
    my $function =
    {
	type        => "void",
	name        => $code::implementationTypeName."Finalize",
	description => "Finalize a {\\normalfont \\ttfamily ".$code::member->{'name'}."} implementation of the {\\normalfont \\ttfamily ".$code::class->{'name'}."} component.",
	content     => "",
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => $code::implementationTypeName,
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     }
	    ]
    };
    unless
	( 
	  grep 
	  {
	           !              $_->{'attributes'}->{'isVirtual'}
	      &&
		  (
		                  $_->{'data'      }->{'rank'     } > 0 
		   ||
		   ! &isIntrinsic($_->{'data'      }->{'type'     }) 
		  ) 
	  }
	  &List::ExtraUtils::hashList($code::member->{'properties'}->{'property'}) 	
	) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
!$GLC attributes unused :: self
CODE
    }
    # Destroy the parent type.
    if ( exists($code::member->{'extends'}) ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
call self%nodeComponent{ucfirst($member->{'class'}).ucfirst($member->{'extends'}->{'name'})}%destroy()
CODE
    }
    # Iterate over properties.
    foreach $code::property ( &List::ExtraUtils::hashList($code::member->{'properties'}->{'property'}) ) {
	# Skip this property if it is virtual.
	next
	    if ( $code::property->{'attributes'}->{'isVirtual'} );
	unless ( &isIntrinsic($code::property->{'data'}->{'type'}) ) {
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
call self%{$property->{'name'}}Data%destroy()
CODE
	}
	if ( $code::property->{'data'}->{'rank'} > 0 ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if (allocated(self%{$property->{'name'}}Data)) deallocate(self%{$property->{'name'}}Data)
CODE
	}
    }
    # Insert a type-binding for this function.
    push(
	@{$build->{'types'}->{$code::implementationTypeName}->{'boundFunctions'}},
	{
	    type        => "procedure",
	    descriptor  => $function,
	    name        => "destroy", 
	}
	);
}

sub Implementation_Builder {
    # Generate a function to build component implementations from XML definitions.
    my $build     = shift();
    $code::class  = shift();
    $code::member = shift();
    # Define the function.
    $code::implementationTypeName = "nodeComponent".ucfirst($code::class->{'name'}).ucfirst($code::member->{'name'});
    my $function  =
    {
	type        => "void",
	name        => $code::implementationTypeName."Builder",
	description => "Build a {\\normalfont \\ttfamily ".$code::member->{'name'}."} implementation of the {\\normalfont \\ttfamily ".$code::class->{'name'}."} component from a supplied XML definition.",
	modules     =>
	    [
	     "Error, only : Error_Report",
	     "FoX_DOM, only : node, nodeList, getLength, extractDataContent, getElementsByTagName, item",
	     "ISO_Varying_String"
	    ],
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => $code::implementationTypeName,
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "type",
		 type       => "node",
		 attributes => [ "intent(in   )", "pointer" ],
		 variables  => [ "componentDefinition" ]
	     },
	     {
		 intrinsic  => "type",
		 type       => "node",
		 attributes => [ "pointer" ],
		 variables  => [ "property" ]
	     },
	     {
		 intrinsic  => "type",
		 type       => "nodeList",
		 attributes => [ "pointer" ],
		 variables  => [ "propertyList" ]
	     },
	     {
		 intrinsic  => "integer",
		 variables  => [ "i", "j" ]
	     }
	    ]
	};
    # Build the function code.
    if ( $code::member->{'name'} eq "null" ) {
	$function->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
!$GLC attributes unused :: componentDefinition
CODE
    } else {
	# Initialize the component.
 	$function->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
call self%initialize()
CODE
	if ( exists($code::member->{'extends'}) ) {
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
call self%nodeComponent{ucfirst($code::member->{'extends'}->{'class'}).ucfirst($code::member->{'extends'}->{'name'})}%builder(componentDefinition)
CODE
	}
 	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
!$omp critical (FoX_DOM_Access)
CODE
	# Iterate over properties.
	foreach $code::property ( &List::ExtraUtils::hashList($code::member->{'properties'}->{'property'}) ) {
	    # Skip this property if it is virtual.
	    next
		if ( $code::property->{'attributes'}->{'isVirtual'} );
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
propertyList => getElementsByTagName(componentDefinition,'{$property->{'name'}}')
CODE
	    if ( $code::property->{'data'}->{'rank'} == 0 ) {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if (getLength(propertyList) > 1) call Error_Report('scalar property must have precisely one value'//\{introspection:location\})
if (getLength(propertyList) == 1) then
  property => item(propertyList,0)
CODE
		if (
		    $code::property->{'data'}->{'type'} eq "double"
		    ||
		    $code::property->{'data'}->{'type'} eq "integer"
		    ||
		    $code::property->{'data'}->{'type'} eq "logical"
		    ) {
		    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
  call extractDataContent(property,self%{$property->{'name'}}Data)
CODE
		} elsif ( $code::property->{'data'}->{'type'} eq "longInteger" ) {
		    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
  call Error_Report('building of long integer properties currently not supported'//\{introspection:location\})
CODE
		} else {
		    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
  call self%{$property->{'name'}}Data%builder(property)
CODE
		}
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
end if
CODE
	    } elsif ( $code::property->{'data'}->{'rank'} == 1 ) {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if (getLength(propertyList) >= 1) then
CODE
		if (
		    $code::property->{'data'}->{'type'} eq "double"
		    ||
		    $code::property->{'data'}->{'type'} eq "integer"
		    ||
		    $code::property->{'data'}->{'type'} eq "logical"
		    ) {
		    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
  allocate(self%{$property->{'name'}}Data(getLength(propertyList)))
  do i=1,getLength(propertyList)
    property => item(propertyList,i-1)
    call extractDataContent(property,self%{$property->{'name'}}Data(i))
  end do
CODE
		} elsif ( $code::property->{'data'}->{'type'} eq "longInteger" ) {
		    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
  call Error_Report('building of long integer properties currently not supported'//\{introspection:location\})
CODE
		} else {
		    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
  allocate(self%{$property->{'name'}}(getLength(propertyList)))
  do i=1,getLength(propertyList)
    property => item(propertyList,i-1)
    call self%{$property->{'name'}}Data(i)%builder(property)
  end do
CODE
		}
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
end if
CODE
	    }
	}
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
!$omp end critical (FoX_DOM_Access)
CODE
    }
    # Build any meta-properties.
    if ( grep {$code::class->{'name'} eq $_} @{$build->{'componentClassListActive'}} ) {
	foreach my $metaPropertyType ( @Galacticus::Build::Components::Classes::MetaProperties::metaPropertyTypes ) {
	    $code::label    = $metaPropertyType->{'label'};
	    $code::rank     = $metaPropertyType->{'rank' };
	    $code::prefix   = ucfirst($metaPropertyType->{'label'})."Rank".$metaPropertyType->{'rank' };
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if (allocated({$class->{'name'}.$prefix}MetaPropertyNames)) then
 !$omp critical (FoX_DOM_Access)
 do i=1,size(({$class->{'name'}.$prefix}MetaPropertyNames))
  propertyList => getElementsByTagName(componentDefinition,char({$class->{'name'}.$prefix}MetaPropertyNames(i)))
CODE
	    if      ( $metaPropertyType->{'rank'} == 0 ) {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
  if (getLength(propertyList) > 1) call Error_Report('meta-property must have precisely one value'//\{introspection:location\})
CODE
		if ( $metaPropertyType->{'intrinsic'} eq "integer" && defined($metaPropertyType->{'type'}) && $metaPropertyType->{'type'} eq "kind_int8" ) {
		    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
  call Error_Report('building of long integer properties currently not supported'//\{introspection:location\})
CODE
		} else {
		    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
  if (getLength(propertyList) == 1) then
    property => item(propertyList,0)
    call extractDataContent(property,self%{$prefix}MetaProperties(i))
  end if
CODE
		}
	    } elsif ( $metaPropertyType->{'rank'} == 1 ) {
		if ( $metaPropertyType->{'intrinsic'} eq "integer" && defined($metaPropertyType->{'type'}) && $metaPropertyType->{'type'} eq "kind_int8" ) {
		    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
  call Error_Report('building of long integer properties currently not supported'//\{introspection:location\})
CODE
		} else {
		    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
  if (getLength(propertyList) >= 1) then
    allocate(self%{$prefix}MetaProperties(i)%values(getLength(propertyList)))
    do j=1,getLength(propertyList)
     property => item(propertyList,j-1)
     call extractDataContent(property,self%{$prefix}MetaProperties(i)%values(j))
    end do
  end if
CODE
		}
	    } else {
		die("Galacticus::Build::Components::Implementations::CreateDestroy::Implementation_Builder: unsupported meta-property rank");
	    }
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
 end do
 !$omp end critical (FoX_DOM_Access)
end if
CODE
	}
    }
    # Insert a type-binding for this function.
    push(
	@{$build->{'types'}->{$code::implementationTypeName}->{'boundFunctions'}},
	{
	    type        => "procedure",
	    descriptor  => $function,
	    name        => "builder", 
	}
	);
}

sub Implementation_Deferred_Create_Set {
    # Generate a function to set the deferred implementation create function.
    my $build     = shift();
    $code::class  = shift();
    $code::member = shift();
    # If creation function is not deferred, we do not need to build a setter function.
    return
	unless ( $code::member->{'createFunction'}->{'isDeferred'} );
    # Create the function.
    my $function =
    {
	type        => "void",
	name        => $code::class->{'name'}.ucfirst($code::member->{'name'})."CreateFunctionSet",
	description => "Set the create function for the {\\normalfont \\ttfamily ".$code::member->{'name'}."} implementation of the {\\normalfont \\ttfamily ".$code::class->{'name'}."} component class.",
	variables   =>
	    [
	     {
		 intrinsic  => "external",
		 variables  => [ "createFunction" ]
	     }
	    ]
    };
    $function->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
{$class->{'name'}.ucfirst($member->{'name'})}CreateFunction => createFunction
CODE
    # Insert a type-binding for this function.
    push(
	@{$build->{'types'}->{"nodeComponent".ucfirst($code::class->{'name'}).ucfirst($code::member->{'name'})}->{'boundFunctions'}},
	{
	    type        => "procedure",
	    descriptor  => $function,
	    name        => "createFunctionSet",
	    pass        => "nopass"
	}
	);
    # Create the associated function pointer.
    push(
	@{$build->{'variables'}},
	{
	    intrinsic  => "procedure",
	    type       => "",
	    attributes => [ "pointer" ],
	    variables  => [ $code::class->{'name'}.ucfirst($code::member->{'name'})."CreateFunction" ]
	}
	);
}

1;
