# Contains a Perl module which handles evolution of component properties during build.

package Galacticus::Build::Components::Properties::Evolve;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
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
	      \&Build_Rate_Get_Functions        ,
	      \&Build_Rate_Functions            ,
	      \&Build_Generic_Rate_Functions    ,
	      \&Build_Auto_Create_Rate_Functions,
	      \&Build_Scale_Functions           ,
	      \&Build_Inactive_Functions        ,
	      \&Build_Analytic_Functions
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
!$GLC attributes unused :: self
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

sub Build_Rate_Get_Functions {
    # Build rate getting functions for non-virtual, evolvable properties.
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
    # Determine the return type of the function.
    (my $functionTypeDescriptor) = &Galacticus::Build::Components::DataTypes::dataObjectDefinition($code::property->{'data'});
    my $functionType = 
	                                                                    $functionTypeDescriptor->{'intrinsic' }            .
	(exists($functionTypeDescriptor->{'type'      }) ? "(" .            $functionTypeDescriptor->{'type'      }  .")" : "").
	(exists($functionTypeDescriptor->{'attributes'}) ? ", ".join(", ",@{$functionTypeDescriptor->{'attributes'}})     : "");
    # Build the function.
    my $implementationTypeName = "nodeComponent".ucfirst($class->{'name'}).ucfirst($member->{'name'});
    my $function =
    {
	type        => $functionType." => propertyRate",
	name        => $class->{'name'}.ucfirst($member->{'name'}).ucfirst($code::property->{'name'})."RateGet",
	description => "Get the rate of change of the {\\normalfont \\ttfamily ".$code::property->{'name'}."} property of an {\\normalfont \\ttfamily ".$member->{'name'}."} implementation of the {\\normalfont \\ttfamily ".$class->{'name'}."} component class.",
	modules     =>
	    [
	     "Error"
	    ],
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => $implementationTypeName,
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "integer",
		 variables  => [ "offset" ]
	     }
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
    @code::argumentsUnused = ( );
    push(@code::argumentsUnused,"self")
	if ( &isIntrinsic($code::property->{'data'}->{'type'}) );
    # Build the function.
    if ( scalar(@code::argumentsUnused) > 0 ) {
	$function->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
!$GLC attributes unused :: {join(",",@argumentsUnused)}
CODE
    }
    $code::offsetNameAll      = &offsetName('all'     ,$class->{'name'}.ucfirst($member->{'name'}),$code::property->{'name'});
    $code::offsetNameActive   = &offsetName('active'  ,$class->{'name'}.ucfirst($member->{'name'}),$code::property->{'name'});
    $code::offsetNameInactive = &offsetName('inactive',$class->{'name'}.ucfirst($member->{'name'}),$code::property->{'name'});
    $code::offset = fill_in_string(<<'CODE', PACKAGE => 'code');
if (nodeAnalytics({$offsetNameAll})) call Error_Report('rates are gettable only for numerically-solved variables'//\{introspection:location\})
if (rateComputeState == propertyTypeInactive) then
 if (nodeInactives({$offsetNameAll})) call Error_Report('rates are gettable only for active variables'//\{introspection:location\})
 offset={$offsetNameActive}
else
 offset=0
 call Error_Report('rates are gettable only during inactive variable integration'//\{introspection:location\})
end if
CODE
    if ( &isIntrinsic($code::property->{'data'}->{'type'}) ) {
	if ( $code::property->{'data'}->{'rank'} == 0 ) {
	    $function->{'content'} .= $code::offset;
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
propertyRate=nodeRatesActives(offset)
CODE
	} else {
	    $function->{'content'} .= $code::offset;
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
count=size(self%{$property->{'name'}}Data)
propertyRate=nodeRatesActives(offset:offset+count-1)
CODE
	}
    } else {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
count=self%{$property->{'name'}}Data%serializeCount()
if (count > 0) then
{$offset}
   propertyRate=self%{$property->{'name'}}Data
   call propertyRate%deserialize(nodeRatesActives(offset:offset+count-1))
end if
CODE
    }
    # Insert a type-binding for this function into the relevant type.
    push(
	@{$build->{'types'}->{$implementationTypeName}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    descriptor  => $function,
	    name        => $code::property->{'name'}."RateGet"
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
	modules     =>
	    [
	     "Error"
	    ],
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
	     },
	     {
		 intrinsic  => "integer",
		 variables  => [ "offset" ]
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
!$GLC attributes unused :: {join(",",@argumentsUnused)}
CODE
    $code::offsetNameAll      = &offsetName('all'     ,$class->{'name'}.ucfirst($member->{'name'}),$code::property->{'name'});
    $code::offsetNameActive   = &offsetName('active'  ,$class->{'name'}.ucfirst($member->{'name'}),$code::property->{'name'});
    $code::offsetNameInactive = &offsetName('inactive',$class->{'name'}.ucfirst($member->{'name'}),$code::property->{'name'});
    $code::offset = fill_in_string(<<'CODE', PACKAGE => 'code');
if (nodeAnalytics({$offsetNameAll})) call Error_Report('rates are only settable for numerically-solved variables'//\{introspection:location\})
if (rateComputeState == propertyTypeAll          ) then
 offset={$offsetNameAll}
else if (rateComputeState == propertyTypeActive  ) then
 if (     nodeInactives({$offsetNameAll})) return
 offset={$offsetNameActive}
else if (rateComputeState == propertyTypeInactive) then
 if (.not.nodeInactives({$offsetNameAll})) return
 offset={$offsetNameInactive}
else if (rateComputeState == propertyTypeNumerics) then
 offset={$offsetNameActive}
else
 return
end if
CODE
    if ( &isIntrinsic($code::property->{'data'}->{'type'}) ) {
	if ( $code::property->{'data'}->{'rank'} == 0 ) {
	    $function->{'content'} .= $code::offset;
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
nodeRates(offset)=nodeRates(offset)+setValue
CODE
	} else {
	    $function->{'content'} .= $code::offset;
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
count=size(setValue)
nodeRates(offset:offset+count-1)=nodeRates(offset:offset+count-1)+setValue
CODE
	}
    } else {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
count=self%{$property->{'name'}}Data%serializeCount()
if (count > 0) then
{$offset}
   current=self%{$property->{'name'}}Data
   call current%deserialize(nodeRates(offset:offset+count-1))
   call current%increment(setValue)
   call current%serialize(nodeRates(offset:offset+count-1))
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
	modules     =>
	    [
	     "Error"
	    ],
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
    push(@{$function->{'modules'}},"Error")
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
  if (.not.(present(interrupt).and.present(interruptProcedure))) call Error_Report('interrupt required, but optional arguments missing'//\{introspection:location\})
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
!$GLC attributes unused :: self
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
if (.not.(present(interrupt).and.present(interruptProcedure))) call Error_Report('interrupt required, but optional arguments missing'//\{introspection:location\})
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
	     {
		 intrinsic  => "integer",
		 variables  => [ "offset" ]
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
!$GLC attributes unused :: self
CODE
    $code::offsetNameAll      = &offsetName('all'     ,$class->{'name'}.ucfirst($member->{'name'}),$code::property->{'name'});
    $code::offsetNameActive   = &offsetName('active'  ,$class->{'name'}.ucfirst($member->{'name'}),$code::property->{'name'});
    $code::offsetNameInactive = &offsetName('inactive',$class->{'name'}.ucfirst($member->{'name'}),$code::property->{'name'});
    $code::offset = fill_in_string(<<'CODE', PACKAGE => 'code');
if (nodeAnalytics({$offsetNameAll})) return
if (rateComputeState == propertyTypeAll          ) then
 offset={$offsetNameAll}
else if (rateComputeState == propertyTypeActive  ) then
 if (     nodeInactives({$offsetNameAll})) return
 offset={$offsetNameActive}
else if (rateComputeState == propertyTypeInactive) then
 if (.not.nodeInactives({$offsetNameAll})) return
 offset={$offsetNameInactive}
else if (rateComputeState == propertyTypeNumerics) then
 offset={$offsetNameActive}
else
 return
end if
CODE
    if ( &isIntrinsic($code::property->{'data'}->{'type'}) ) {
	if ( $code::property->{'data'}->{'rank'} == 0 ) {
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
{$offset}
nodeScales(offset)=setValue
CODE
	} else {
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
{$offset}
nodeScales(offset:offset+size(setValue))=setValue
CODE
	}
    } else {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
count=setValue%serializeCount()
if (count > 0) then
{$offset}
   call setValue%serialize(nodeScales(offset:offset+count-1))
end if
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

sub Build_Inactive_Functions {
    # Build functions to indicate variables which are inactive (i.e. do not appear on the right-hand side of any differential equation being solved) for non-virtual, evolvable properties.
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
	type        => "void",
	name        => $class->{'name'}.ucfirst($member->{'name'}).ucfirst($code::property->{'name'})."JcbnZr",
	description => "Indicate that the {\\normalfont \\ttfamily ".$code::property->{'name'}."} property of an {\\normalfont \\ttfamily ".$member->{'name'}."} implementation of the {\\normalfont \\ttfamily ".$class->{'name'}."} component class is inactive for differential equation solving.",
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => $implementationTypeName,
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     }
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
!$GLC attributes unused :: self
CODE
    $code::offsetName = &offsetName('all',$class->{'name'}.ucfirst($member->{'name'}),$code::property->{'name'});
    if ( &isIntrinsic($code::property->{'data'}->{'type'}) ) {
	if ( $code::property->{'data'}->{'rank'} == 0 ) {
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
nodeInactives({$offsetName})=.true.
CODE
	} else {
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
nodeInactives({$offsetName}:{$offsetName}+size(self%{$property->{'name'}}Data)-1)=.true.
CODE
	}
    } else {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
count=self%{$property->{'name'}}Data%serializeCount()
if (count > 0) nodeInactives({$offsetName}:{$offsetName}+count-1)=.true.
CODE
    }
    # Insert a type-binding for this function into the relevant type.
    push(
	@{$build->{'types'}->{$implementationTypeName}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    descriptor  => $function,
	    name        => $code::property->{'name'}."Inactive"
	}
	);  
}

sub Build_Analytic_Functions {
    # Build functions to indicate variables which are to be solved analytically for non-virtual, evolvable properties.
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
	type        => "void",
	name        => $class->{'name'}.ucfirst($member->{'name'}).ucfirst($code::property->{'name'})."Alytc",
	description => "Indicate that the {\\normalfont \\ttfamily ".$code::property->{'name'}."} property of an {\\normalfont \\ttfamily ".$member->{'name'}."} implementation of the {\\normalfont \\ttfamily ".$class->{'name'}."} component class is to be solved analytically during differential evolution.",
	modules     =>
	    [
	     "Error"
	    ],
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => $implementationTypeName,
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     }
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
!$GLC attributes unused :: self
CODE
    $code::offsetName = &offsetName('all',$class->{'name'}.ucfirst($member->{'name'}),$code::property->{'name'});
    if ( &isIntrinsic($code::property->{'data'}->{'type'}) ) {
	if ( $code::property->{'data'}->{'rank'} == 0 ) {
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if (nodeAnalytics({$offsetName})) call Error_Report('property is already marked analytically-solvable'//\{introspection:location\})
nodeAnalytics({$offsetName})=.true.
CODE
	} else {
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if (any(nodeAnalytics({$offsetName}:{$offsetName}+size(self%{$property->{'name'}}Data)-1))) call Error_Report('property is already marked analytically-solvable'//\{introspection:location\})
nodeAnalytics({$offsetName}:{$offsetName}+size(self%{$property->{'name'}}Data)-1)=.true.
CODE
	}
    } else {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
count=self%{$property->{'name'}}Data%serializeCount()
if (count > 0) then
 if (any(nodeAnalytics({$offsetName}:{$offsetName}+count-1))) call Error_Report('property is already marked analytically-solvable'//\{introspection:location\})
 nodeAnalytics({$offsetName}:{$offsetName}+count-1)=.true.
end if
CODE
    }
    # Insert a type-binding for this function into the relevant type.
    push(
	@{$build->{'types'}->{$implementationTypeName}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    descriptor  => $function,
	    name        => $code::property->{'name'}."Analytic"
	}
	);  
}

1;
