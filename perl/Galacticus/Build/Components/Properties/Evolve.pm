# Contains a Perl module which handles evolution of component properties during build.

package Galacticus::Build::Components::Properties::Evolve;
use strict;
use warnings;
use utf8;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use Data::Dumper;
use List::ExtraUtils;
use Galacticus::Build::Components::Utils qw(&isIntrinsic &offsetName);
use Text::Template 'fill_in_string';

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     propertiesEvolve => 
     {
	 propertyIteratedFunctions =>
	     [
	      \&Build_Count_Functions           ,
	      \&Build_Rate_Functions            ,
	      \&Build_Generic_Rate_Functions    ,
	      \&Build_Auto_Create_Rate_Functions,
	      \&Build_Scale_Functions
	     ]
     }
    );

sub Build_Count_Functions {
    # Build count functions for non-virtual, evolvable properties.
    my $build       = shift();
    my $class       = shift();
    my $member      = shift();
    $code::property = shift();
    # Skip this property if it is not evolvable, or is virtual.
    return
	if
	(
	    $code::property->{'attributes' }->{'isVirtual'  }
	 ||
	  ! $code::property->{'attributes' }->{'isEvolvable'}
	);
    # Build the function.
    my $implementationTypeName = "nodeComponent".ucfirst($class->{'name'}).ucfirst($member->{'name'});
    my $function =
    {
	type        => "integer",
	name        => $class->{'name'}.ucfirst($member->{'name'}).ucfirst($code::property->{'name'})."Count",
	description => "Return a count of the number of scalar properties in the {\\normalfont \\ttfamily ".$code::property->{'name'}."} property of an {\\normalfont \\ttfamily ".$member->{'name'}."} implementation of the {\\normalfont \\ttfamily ".$class->{'name'}."} component class.",
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => $implementationTypeName,
		 attributes => [ "intent(in   )" ],
		 variables  => [ "self" ]
	     }
	    ]
    };
    # Build the function.
    $code::typePrefix = $class->{'name'}.ucfirst($member->{'name'});
    if ( $code::property->{'data'}->{'rank'} == 0 ) {
	$function->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
!GCC$ attributes unused :: self
{$typePrefix.ucfirst($property->{'name'})}Count=1
CODE
    } elsif ( $code::property->{'data'}->{'rank'} == 1 ) {
	$function->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
if (allocated(self%{$property->{'name'}}Data)) then
   {$typePrefix.ucfirst($property->{'name'})}Count=size(self%{$property->{'name'}}Data)
else
   {$typePrefix.ucfirst($property->{'name'})}Count=0
end if
CODE
    }
    # Insert a type-binding for this function into the relevant type.
    push(
	@{$build->{'types'}->{$implementationTypeName}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    descriptor  => $function,
	    name        => $code::property->{'name'}."Count"
	}
	);  
}

sub Build_Rate_Functions {
    # Build rate setting functions for non-virtual, evolvable properties.
    my $build       = shift();
    my $class       = shift();
    my $member      = shift();
    $code::property = shift();
    # Skip this property if it is not evolvable, or is virtual.
    return
	if
	(
	    $code::property->{'attributes' }->{'isVirtual'  }
	 ||
	  ! $code::property->{'attributes' }->{'isEvolvable'}
	);
    # If rate function is deferred, then create an intrinsic version.
    my $intrinsicRate = grep {$_ eq "rate"} split(":",$code::property->{'attributes'}->{'isDeferred'});
    my $suffix        = $intrinsicRate ? "Intrinsic" : "";
    # Get the data descriptor for this propery.
    (my $propertyTypeDescriptor) = &Galacticus::Build::Components::DataTypes::dataObjectDefinition($code::property->{'data'},matchOnly => 1);
    push(@{$propertyTypeDescriptor->{'variables' }},"setValue"     );
    push(@{$propertyTypeDescriptor->{'attributes'}},"intent(in   )");
    # Build the function.
    my $implementationTypeName = "nodeComponent".ucfirst($class->{'name'}).ucfirst($member->{'name'});
    my $function =
    {
	type        => "void",
	name        => $class->{'name'}.ucfirst($member->{'name'}).ucfirst($code::property->{'name'})."Rate".$suffix,
	description => "Accumulate".($intrinsicRate ? " directly (i.e. circumventing any deferred function binding)" : "")." to the rate of change of the {\\normalfont \\ttfamily ".$code::property->{'name'}."} property of an {\\normalfont \\ttfamily ".$member->{'name'}."} implementation of the {\\normalfont \\ttfamily ".$class->{'name'}."} component class.",
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => $implementationTypeName,
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     },
	     $propertyTypeDescriptor,
	     {
		 intrinsic  => "logical",
		 attributes => [ "optional", "intent(inout)" ],
		 variables  => [ "interrupt" ]
	     },
	     {
		 intrinsic  => "procedure",
		 type       => "interruptTask",
		 attributes => [ "optional", "intent(inout)", "pointer" ],
		 variables  => [ "interruptProcedure" ]
	     }
	    ]
    };
    # For non-intrinsic types add a workspace variable.
    unless ( &isIntrinsic($code::property->{'data'}->{'type'}) ) {
	(my $currentTypeDescriptor) = &Galacticus::Build::Components::DataTypes::dataObjectDefinition($code::property->{'data'},matchOnly => 1);
	push(@{$currentTypeDescriptor->{'variables' }},"current");
	push(@{$function->{'variables'}},$currentTypeDescriptor)
    }
    # Add a count variable if necessary.
    push(
	@{$function->{'variables'}},
	{
	    intrinsic  => "integer",
	    variables  => [ "count" ]
	}				
	)
	if (
	    ( 
	      &isIntrinsic($code::property->{'data'}->{'type'})
	      &&
	                   $code::property->{'data'}->{'rank'} > 0
	    )
	    ||
	    ! &isIntrinsic($code::property->{'data'}->{'type'})   
	);
    # Determine which arguments are unused.
    @code::argumentsUnused = ( "interrupt", "interruptProcedure" );
    push(@code::argumentsUnused,"self")
	if ( &isIntrinsic($code::property->{'data'}->{'type'}) );
    # Build the function.
    $function->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
!GCC$ attributes unused :: {join(",",@argumentsUnused)}
CODE
    $code::offsetName = &offsetName($class->{'name'}.ucfirst($member->{'name'}),$code::property->{'name'});
    if ( &isIntrinsic($code::property->{'data'}->{'type'}) ) {
	if ( $code::property->{'data'}->{'rank'} == 0 ) {
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
nodeRates({$offsetName})=nodeRates({$offsetName})+setValue
CODE
	} else {
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
count=size(setValue)
nodeRates({$offsetName}:{$offsetName}+count-1)=nodeRates({$offsetName}:{$offsetName}+count-1)+setValue
CODE
	}
    } else {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
count=self%{$property->{'name'}}Data%serializeCount()
if (count > 0) then
   current=self%{$property->{'name'}}Data
   call current%deserialize(nodeRates({$offsetName}:{$offsetName}+count-1))
   call current%increment(setValue)
   call current%serialize(nodeRates({$offsetName}:{$offsetName}+count-1))
end if
CODE
    }
    # Insert a type-binding for this function into the relevant type.
    push(
	@{$build->{'types'}->{$implementationTypeName}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    descriptor  => $function,
	    name        => $code::property->{'name'}."Rate".$suffix
	}
	);  
}

sub Build_Generic_Rate_Functions {
    # Build rate setting functions for non-virtual, evolvable properties which need to be accessible via a generic nodeComponent class.
    my $build       = shift();
    $code::class    = shift();
    $code::member   = shift();
    $code::property = shift();
    # Skip this property if it is not evolvable, or is virtual.
    return
	if
	(
	   $code::property->{'attributes'}->{'isVirtual'  }
	 ||
	 ! $code::property->{'attributes'}->{'isEvolvable'}
	 ||
	 ! $code::property->{'attributes'}->{'makeGeneric'}
	);
    # Get the data descriptor for this propery.
    (my $propertyTypeDescriptor) = &Galacticus::Build::Components::DataTypes::dataObjectDefinition($code::property->{'data'},matchOnly => 1);
    push(@{$propertyTypeDescriptor->{'variables' }},"setValue"     );
    push(@{$propertyTypeDescriptor->{'attributes'}},"intent(in   )");
    # Build the function.
    my $function =
    {
	type        => "void",
	name        => $code::class->{'name'}.ucfirst($code::member->{'name'}).ucfirst($code::property->{'name'})."RateGeneric",
	description => "Set the rate of the {\\normalfont \\ttfamily ".$code::property->{'name'}."} property of the {\\normalfont \\ttfamily ".$code::property->{'name'}."} property of an {\\normalfont \\ttfamily ".$code::member->{'name'}."} implementation of the {\\normalfont \\ttfamily ".$code::class->{'name'}."} component class via a generic {\\normalfont \\ttfamily nodeComponent}.",
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent",
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     },
	     $propertyTypeDescriptor,
	     {
		 intrinsic  => "logical",
		 attributes => [ "optional", "intent(inout)" ],
		 variables  => [ "interrupt" ]
	     },
	     {
		 intrinsic  => "procedure",
		 type       => "interruptTask",
		 attributes => [ "optional", "intent(inout)", "pointer" ],
		 variables  => [ "interruptProcedure" ]
	     },
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($code::class->{'name'}),
		 attributes => [ "pointer" ],
		 variables  => [ $code::class->{'name'} ]
	     },
	     {
		 intrinsic  => "type",
		 type       => "treeNode",
		 attributes => [ "pointer" ],
		 variables  => [ "selfNode" ]
	     }
	    ]
    };
    # Add error reporting module if required.
    push(@{$function->{'modules'}},"Galacticus_Error")
	if ( $code::property->{'attributes'}->{'createIfNeeded'} );
    # Build the function.
    $function->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
selfNode => self%host()
{$class->{'name'}} => selfNode%{$class->{'name'}}()
CODE
    if ( $code::property->{'attributes'}->{'createIfNeeded'} ) {
    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
select type ({$class->{'name'}})
type is (nodeComponent{ucfirst($class->{'name'})})
  ! No specific component exists, we must interrupt and create one.
CODE
	if ( $code::property->{'data'}->{'rank'} == 0 ) {
	    if ( &isIntrinsic($code::property->{'data'}->{'type'}) ) {
    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
  if (setValue == 0.0d0) return
CODE
	    } else {
    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
  if (setValue%isZero()) return
CODE
	    }
	} else {
	    if ( &isIntrinsic($code::property->{'data'}->{'type'}) ) {
    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
  if (all(setValue == 0.0d0)) return
CODE
	    }
	}
    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
  if (.not.(present(interrupt).and.present(interruptProcedure))) call Galacticus_Error_Report('{$class->{'name'}.ucfirst($member->{'name'}).ucfirst($property->{'name'})}RateGeneric','interrupt required, but optional arguments missing')
  interrupt          =  .true.
  interruptProcedure => {$class->{'name'}}CreateByInterrupt
  return
end select
CODE
    }
    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
call {$class->{'name'}}%{$property->{'name'}}Rate(setValue,interrupt,interruptProcedure)
CODE
    # Insert the function into the function list. There is no type-binding for this function - the expectation is that it will be
    # bound to some deferred rate function at run-time.
    push(
	@{$build->{'functions'}},
	$function
	);
}

sub Build_Auto_Create_Rate_Functions {
    # Build rate setting functions for non-virtual, evolvable properties which are bound to the component class and so can be used to auto-create components.
    my $build       = shift();
    $code::class    = shift();
    my $member      = shift();
    $code::property = shift();
    # Skip this property if it is not evolvable, or is virtual.
    return
	if
	(
	   $code::property->{'attributes'}->{'isVirtual'     }
	 ||
	 ! $code::property->{'attributes'}->{'isEvolvable'   }
	 ||
	 ! $code::property->{'attributes'}->{'createIfNeeded'}
	);

    # Construct method name, skip this property if the method already exists.
    my $methodName = $code::property->{'name'}."Rate";
    return
	if ( grep {$_->{'name'} eq $methodName} @{$build->{'types'}->{"nodeComponent".ucfirst($code::class->{'name'})}->{'boundFunctions'}} );
    # Get the data descriptor for this propery.
    (my $propertyTypeDescriptor) = &Galacticus::Build::Components::DataTypes::dataObjectDefinition($code::property->{'data'},matchOnly => 1);
    push(@{$propertyTypeDescriptor->{'variables' }},"setValue"     );
    push(@{$propertyTypeDescriptor->{'attributes'}},"intent(in   )");
    # Build the function.
    my $function =
    {
	type        => "void",
	name        => $code::class->{'name'}.ucfirst($code::property->{'name'})."Rate",
	description => "Accept a rate set for the {\\normalfont \\ttfamily ".$code::property->{'name'}."} property of the {\\normalfont \\ttfamily ".$code::class->{'name'}."} component class. Trigger an interrupt to create the component.",
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
	     },
	     $propertyTypeDescriptor,
	     {
		 intrinsic  => "logical",
		 attributes => [ "optional", "intent(inout)" ],
		 variables  => [ "interrupt" ]
	     },
	     {
		 intrinsic  => "procedure",
		 type       => "interruptTask",
		 attributes => [ "optional", "intent(inout)", "pointer" ],
		 variables  => [ "interruptProcedure" ]
	     }
	    ]
    };
    # Build the function.
    $function->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
!GCC$ attributes unused :: self
! No specific component exists, so we must interrupt and create one unless the rate is zero.
CODE
    if ( $code::property->{'data'}->{'rank'} == 0 ) {
	if ( &isIntrinsic($code::property->{'data'}->{'type'}) ) {
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if (setValue == 0.0d0) return
CODE
	} else {
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if (setValue%isZero()) return
CODE
	}
    } else {
	if ( &isIntrinsic($code::property->{'data'}->{'type'}) ) {
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if (all(setValue == 0.0d0)) return
CODE
	}
    }
    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if (.not.(present(interrupt).and.present(interruptProcedure))) call Galacticus_Error_Report('{$class->{'name'}.ucfirst($property->{'name'})}Rate','interrupt required, but optional arguments missing')
interrupt=.true.
interruptProcedure => {$class->{'name'}}CreateByInterrupt
CODE
    # Insert a type-binding for this function into the relevant type.
    push(
	@{$build->{'types'}->{"nodeComponent".ucfirst($code::class->{'name'})}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    descriptor  => $function,
	    name        => $methodName
	}
	);  
}

sub Build_Scale_Functions {
    # Build scale setting functions for non-virtual, evolvable properties.
    my $build       = shift();
    my $class       = shift();
    my $member      = shift();
    $code::property = shift();
    # Skip this property if it is not evolvable, or is virtual.
    return
	if
	(
	    $code::property->{'attributes' }->{'isVirtual'  }
	 ||
	  ! $code::property->{'attributes' }->{'isEvolvable'}
	);
    # Get the data descriptor for this propery.
    (my $propertyTypeDescriptor) = &Galacticus::Build::Components::DataTypes::dataObjectDefinition($code::property->{'data'},matchOnly => 1);
    push(@{$propertyTypeDescriptor->{'variables' }},"setValue"     );
    push(@{$propertyTypeDescriptor->{'attributes'}},"intent(in   )");
    # Build the function.
    my $implementationTypeName = "nodeComponent".ucfirst($class->{'name'}).ucfirst($member->{'name'});
    my $function =
    {
	type        => "void",
	name        => $class->{'name'}.ucfirst($member->{'name'}).ucfirst($code::property->{'name'})."Scale",
	description => "Set the absolute scale of the {\\normalfont \\ttfamily ".$code::property->{'name'}."} property of an {\\normalfont \\ttfamily ".$member->{'name'}."} implementation of the {\\normalfont \\ttfamily ".$class->{'name'}."} component class.",
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => $implementationTypeName,
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     },
	     $propertyTypeDescriptor
	    ]
    };
    # Add a count variable if necessary.
    push(
	@{$function->{'variables'}},
	{
	    intrinsic  => "integer",
	    variables  => [ "count" ]
	}				
	)
	if ( ! &isIntrinsic($code::property->{'data'}->{'type'}) );
    # Build the function.
    $function->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
!GCC$ attributes unused :: self
CODE
    $code::offsetName = &offsetName($class->{'name'}.ucfirst($member->{'name'}),$code::property->{'name'});
    if ( &isIntrinsic($code::property->{'data'}->{'type'}) ) {
	if ( $code::property->{'data'}->{'rank'} == 0 ) {
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
nodeScales({$offsetName})=setValue
CODE
	} else {
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
nodeScales({$offsetName}:{$offsetName}+size(setValue))=setValue
CODE
	}
    } else {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
count=setValue%serializeCount()
if (count > 0) call setValue%serialize(nodeScales({$offsetName}:{$offsetName}+count-1))
CODE
    }
    # Insert a type-binding for this function into the relevant type.
    push(
	@{$build->{'types'}->{$implementationTypeName}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    descriptor  => $function,
	    name        => $code::property->{'name'}."Scale"
	}
	);  
}

1;
