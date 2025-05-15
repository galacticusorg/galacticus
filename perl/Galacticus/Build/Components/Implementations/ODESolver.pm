# Contains a Perl module which provides various ODE solver-related functions for component implementations.

package Galacticus::Build::Components::Implementations::ODESolver;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Text::Template 'fill_in_string';
use Data::Dumper;
use List::ExtraUtils;
use Galacticus::Build::Components::Utils qw(offsetName);
use Galacticus::Build::Components::DataTypes;
use Galacticus::Build::Components::Implementations::Utils;

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     implementationODESolver =>
     {
	 implementationIteratedFunctions =>
	     [
	      \&Implementation_ODE_Serialize_Count      ,
	      \&Implementation_ODE_Serialize_Values     ,
	      \&Implementation_ODE_Deserialize_Values   ,
	      \&Implementation_ODE_Serialize_NonNegative,
	      \&Implementation_ODE_Name_From_Index      ,
	      \&Implementation_ODE_Offsets              ,
	      \&Implementation_ODE_Offset_Variables
	     ],
	 functions =>
	     [
	      \&Implementation_ODE_Rate_Variables
	     ]
     }
    );

sub Implementation_ODE_Name_From_Index {
    # Generate a function to return the name of a property given the index of that property in the serialization of a component
    # implementation.
    my $build     = shift();
    $code::class  = shift();
    $code::member = shift();
    my $implementationTypeName = "nodeComponent".ucfirst($code::class->{'name'}).ucfirst($code::member->{'name'});
    my $function =
    {
	type        => "type(varying_string) => name",
	name        => $implementationTypeName."NameFromIndex",
	description => "Return the name of the property of given index for a {\\normalfont \\ttfamily ".$code::member->{'name'}."} implementation of the {\\normalfont \\ttfamily ".$code::class->{'name'}."} component class.",
	modules     =>
	    [
	     "ISO_Varying_String"
	    ],
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => $implementationTypeName,
		 attributes => [ "intent(in   )" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => [ "intent(inout)" ],
		 variables  => [ "count" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "propertyType" ]
	     }
	    ]
    };
    # Determine if "self" will be used. It is used iff the implementation extends another implementation, or if any
    # property is evolveable and is not a rank-0 double.
    undef(@code::unused);
    push(@code::unused,"self")
	unless (
	    exists($code::member->{'extends'})
	    ||
	    &Galacticus::Build::Components::Implementations::Utils::hasRealNonTrivialEvolvers($code::member)
	);
    # Determine if "count" will be used. It is used iff the implementation extends another implementation, or if any
    # property is evolveable, or if this is a non-extended type.
    push(@code::unused,"count","propertyType")
	unless (
	    exists($code::member->{'extends'})
	    ||
	    &Galacticus::Build::Components::Implementations::Utils::hasRealEvolvers          ($code::member)
	    ||
	    ! exists($code::member->{'extends'})
	);
    # Add an iterator variable for non-extended types.
    push(
	@{$function->{'variables'}},
	{
	    intrinsic  => "integer",
	    variables  => [ "i" ]
	}
	)
	if ( ! exists($code::member->{'extends'}) );
    # Build the function.
    $function->{'content'}  = "";
    if ( scalar(@code::unused) > 0 ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
!$GLC attributes unused :: {join(",",@unused)}
CODE
    }
    # If this component is an extension, first call on the extended type.
    if ( exists($code::member->{'extends'}) ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
name=self%nodeComponent{ucfirst($member->{'extends'}->{'class'}).ucfirst($member->{'extends'}->{'name'})}%nameFromIndex(count,propertyType)
if (count <= 0) return
CODE
    } elsif ( grep {$code::class->{'name'} eq $_} @{$build->{'componentClassListActive'}} ) {
	# Include meta-properties here.
	$code::offsetName = &offsetName('all',$code::class->{'name'},'floatRank0MetaProperties');
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if (allocated({$class->{'name'}}FloatRank0MetaPropertyNames)) then
 do i=1,size({$class->{'name'}}FloatRank0MetaPropertyNames)
   if (                                                                                    &
    &  {$class->{'name'}}FloatRank0MetaPropertyEvolvable(i)                                &
    & .and.                                                                                &
    &                                                 .not.nodeAnalytics({$offsetName}(i)) &
    & .and.                                                                                &
    &  (                                                                                   &
    &     propertyType == propertyTypeAll                                                  &
    &   .or.                                                                               &
    &    (propertyType == propertyTypeActive   .and. .not.nodeInactives({$offsetName}(i))) &
    &   .or.                                                                               &
    &    (propertyType == propertyTypeInactive .and.      nodeInactives({$offsetName}(i))) &
    &   .or.                                                                               &
    &     propertyType == propertyTypeNumerics                                             &
    &  )                                                                                   &
    & ) count=count-1
  if (count <= 0) then
   name={$class->{'name'}}FloatRank0MetaPropertyNames(i)
   return
  end if
 end do
end if
CODE
    }
    # Iterate over non-virtual, evolvable properties.
    foreach $code::property ( &Galacticus::Build::Components::Implementations::Utils::listRealEvolvers($code::member) ) {
	# Find condition for count update, including active/inactive status. For allocatable properties, condition is that the
	# object be allocated. For non-allocatable properties, always update the count.
	if ( $code::property->{'data'}->{'rank'} > 0 ) {
	    $code::conditionOpen  = "if (allocated(self%".$code::property->{'name'}."Data)) then";
	    $code::conditionClose = "end if";
	} else {
	    $code::conditionOpen  = "";
	    $code::conditionClose = "";
	}
	# Find the size of the object. Rank-0 double properties always have a count of 1. For other rank-0 types, call
	# their serialization count method. Rank>0 must be double, so simply use the array length.
	$code::count     = 
	    $code::property->{'data'}->{'rank'} > 0 
	    ?
	    "size(self%".$code::property->{'name'}."Data)"
	    :
	    (
	     $code::property->{'data'}->{'type'} eq "double"
	     ?
	     "1" 
	     : 
	     "self%".$code::property->{'name'}."Data%serializeCount()"
	    );
	$code::offsetName = &offsetName('all',$code::class,$code::member,$code::property);
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
{$conditionOpen}
 if (                                                                                                               &
  &                                                all(.not.nodeAnalytics({$offsetName}:{$offsetName}+{$count}-1))  &
  &  .and.                                                                                                          &
  &  (                                                                                                              &
  &     propertyType == propertyTypeAll                                                                             &
  &   .or.                                                                                                          &
  &    (propertyType == propertyTypeActive   .and. all(.not.nodeInactives({$offsetName}:{$offsetName}+{$count}-1))) &
  &   .or.                                                                                                          &
  &    (propertyType == propertyTypeInactive .and. all(     nodeInactives({$offsetName}:{$offsetName}+{$count}-1))) &
  &   .or.                                                                                                          &
  &     propertyType == propertyTypeNumerics                                                                        &
  &  )                                                                                                              &
  & ) count=count-{$count}
{$conditionClose}
CODE
	        $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if (count <= 0) then
  name='{$class->{'name'}}:{$member->{'name'}}:{$property->{'name'}}'
  return
end if
CODE
    }
    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
name='?'
CODE
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{$implementationTypeName}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    descriptor  => $function,
	    name        => "nameFromIndex"
	}
	);	    
}

sub Implementation_ODE_Serialize_Count {
    # Generate a function to return a count of serializable, evolvable properties of each component implementation for the ODE solver.
    my $build     = shift();
    $code::class  = shift();
    $code::member = shift();
    $code::implementationTypeName = "nodeComponent".ucfirst($code::class->{'name'}).ucfirst($code::member->{'name'});
    my $function =
    {
	type        => "integer",
	name        => $code::implementationTypeName."SerializeCount",
	description => "Return a count of the serialization of a ".$code::member->{'name'}." implementation of the ".$code::class->{'name'}." component.",
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => $code::implementationTypeName,
		 attributes => [ "intent(in   )" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "propertyType" ]
	     }
	    ]
    };
    push(
	@{$function->{'variables'}},
	{
	    intrinsic  => "integer",
	    variables  => [ "count" ]
	}
	)
	if ( &Galacticus::Build::Components::Implementations::Utils::hasRealEvolvers($code::member) || ! exists($code::member->{'extends'}) );
    push(
	@{$function->{'variables'}},
	{
	    intrinsic  => "integer",
	    variables  => [ "i" ]
	}
	)
	if (                                                                                           ! exists($code::member->{'extends'}) );
    # Determine if "self" will be used. It is used iff the implementation extends another implementation, or if any
    # property is evolveable and is not a rank-0 double.
    undef(@code::unused);
    push(@code::unused,"self")
	unless (
	    exists($code::member->{'extends'})
	    ||
	    &Galacticus::Build::Components::Implementations::Utils::hasRealNonTrivialEvolvers($code::member)
	);
    push(@code::unused,"propertyType")
	unless ( &Galacticus::Build::Components::Implementations::Utils::hasRealEvolvers($code::member) );
    # Build the function.
    $function->{'content'}  = "";
    if ( scalar(@code::unused) > 0 ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
!$GLC attributes unused :: {join(",",@unused)}
CODE
    }
    # If this component is an extension, first call on the extended type.
    if ( exists($code::member->{'extends'}) ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
{$implementationTypeName}SerializeCount=self%nodeComponent{ucfirst($member->{'extends'}->{'class'}).ucfirst($member->{'extends'}->{'name'})}%serializeCount(propertyType)
CODE
    } elsif ( grep {$code::class->{'name'} eq $_} @{$build->{'componentClassListActive'}} ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
{$implementationTypeName}SerializeCount=0
if (allocated({$class->{'name'}}FloatRank0MetaPropertyNames)) then
 count=0
 do i=1,size({$class->{'name'}}FloatRank0MetaPropertyNames)
  if (.not.{$class->{'name'}}FloatRank0MetaPropertyEvolvable(i)) cycle
  if      (propertyType == propertyTypeAll     ) then
                                                                                                                                                                  {$implementationTypeName}SerializeCount={$implementationTypeName}SerializeCount+1
  else if (propertyType == propertyTypeInactive) then
   if (     nodeInactives(offsetAll{$class->{'name'}}FloatRank0MetaProperties(i)).and..not.nodeAnalytics(offsetAll{$class->{'name'}}FloatRank0MetaProperties(i))) {$implementationTypeName}SerializeCount={$implementationTypeName}SerializeCount+1
  else if (propertyType == propertyTypeActive  ) then
   if (.not.nodeInactives(offsetAll{$class->{'name'}}FloatRank0MetaProperties(i)).and..not.nodeAnalytics(offsetAll{$class->{'name'}}FloatRank0MetaProperties(i))) {$implementationTypeName}SerializeCount={$implementationTypeName}SerializeCount+1
  else if (propertyType == propertyTypeNumerics) then
   if (                                                                               .not.nodeAnalytics(offsetAll{$class->{'name'}}FloatRank0MetaProperties(i))) {$implementationTypeName}SerializeCount={$implementationTypeName}SerializeCount+1
  end if
  count=count+1
 end do
end if
CODE
    }
    # Iterate over non-virtual, evolvable properties.
    foreach $code::property ( &Galacticus::Build::Components::Implementations::Utils::listRealEvolvers($code::member) ) {
	# Find the size of the object. Rank-0 double properties always have a count of 1. For other rank-0 types, call
	# their serialization count method. Rank>0 must be double, so simply use the array length.
	$code::count     = 
	    $code::property->{'data'}->{'rank'} > 0 
	    ?
	    "size(self%".$code::property->{'name'}."Data)"
	    :
	    (
	     $code::property->{'data'}->{'type'} eq "double"
	     ?
	     "1" 
	     : 
	     "self%".$code::property->{'name'}."Data%serializeCount()"
	    );
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
count={$count}
CODE
	# Find condition for count update. For allocatable properties, condition is that the object be allocated. For
	# non-allocatable properties, always update the count.
	$code::offsetName = &offsetName('all',$code::class,$code::member,$code::property);
	if ( $code::property->{'data'}->{'rank'} > 0 ) {
	    $code::conditionOpen  = "if (allocated(self%".$code::property->{'name'}."Data)) then";
	    $code::conditionClose = "end if";
	} else {
	    $code::conditionOpen  = "";
	    $code::conditionClose = "";
	}
	# Increment the count.
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
{$conditionOpen}
if (propertyType == propertyTypeAll) then
 {$implementationTypeName}SerializeCount={$implementationTypeName}SerializeCount+count
else if (((propertyType == propertyTypeActive .and. all(.not.nodeInactives({$offsetName}:{$offsetName}+count-1))) .or. (propertyType == propertyTypeInactive .and. all(nodeInactives({$offsetName}:{$offsetName}+count-1))) .or. propertyType == propertyTypeNumerics) .and. all(.not.nodeAnalytics({$offsetName}:{$offsetName}+count-1))) then
 {$implementationTypeName}SerializeCount={$implementationTypeName}SerializeCount+count
end if
{$conditionClose}
CODE
    }
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{$code::implementationTypeName}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    descriptor  => $function,
	    name        => "serializeCount"
	}
	);	    
}

sub Implementation_ODE_Serialize_Values {
    # Generate a function to serialize values of evolvable properties of each component implementation to array for the ODE solver.
    my $build     = shift();
    $code::class  = shift();
    $code::member = shift();
    my $implementationTypeName = "nodeComponent".ucfirst($code::class->{'name'}).ucfirst($code::member->{'name'});
    my $function =
    {
	type        => "void",
	name        => $implementationTypeName."SerializeValues",
	description => "Serialize evolvable properties of a {\\normalfont \\ttfamily ".$code::member->{'name'}."} implementation of the {\\normalfont \\ttfamily ".$code::class->{'name'}."} component to array.",
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => $implementationTypeName,
		 attributes => [ "intent(in   )" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "double precision",
		 attributes => [ "intent(  out)", "dimension(:)" ],
		 variables  => [ "array" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "propertyType" ]
	     }
	    ]
    };
    # Conditionally add "offset", "count" and indexing variables if they will be needed.
    my @requiredVariables = ( "offset" );
    push(@requiredVariables,"count")
	if
	(
	 exists($code::member->{'extends'})
	 ||
	 &Galacticus::Build::Components::Implementations::Utils::hasRealNonTrivialEvolvers($code::member)
	);
    push(@requiredVariables,"i")
	if ( ! exists($code::member->{'extends'}) );
    push(@{$function->{'variables'}},
	 {
	     intrinsic  => "integer",
	     variables  => \@requiredVariables
	 }
	)
	if ( scalar(@requiredVariables) > 0 );
    # Determine if the function arguments are unused.
    @code::unused = 
	(
	 exists($code::member->{'extends'})
	 ||
	 &Galacticus::Build::Components::Implementations::Utils::hasRealEvolvers          ($code::member)
	)
	?
	()
	:
	("array","propertyType");
    # Build the function.
    $function->{'content'} = "";
    if ( scalar(@code::unused) > 0 ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
!$GLC attributes unused :: {join(",",@unused)}
CODE
    }
    # Initialize offset if required.
    if ( grep {$_ eq "offset"} @requiredVariables ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');	
offset=1
CODE
    }    
    # If this component is an extension, call serialization on the extended type.
    if ( exists($code::member->{'extends'}) ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');	
count=self%nodeComponent{ucfirst($code::member->{'extends'}->{'class'}).ucfirst($code::member->{'extends'}->{'name'})}%serializeCount (      propertyType)
if (count > 0) then
 call self%nodeComponent{ucfirst($code::member->{'extends'}->{'class'}).ucfirst($code::member->{'extends'}->{'name'})}%serializeValues(array,propertyType)
 offset=offset+count
end if
CODE
    } elsif ( grep {$code::class->{'name'} eq $_} @{$build->{'componentClassListActive'}} ) {
	# For non-extended types serialize meta-properties.
	$code::offsetName = &offsetName('all',$code::class->{'name'},'floatRank0MetaProperties');
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');	
if (allocated({$class->{'name'}}FloatRank0MetaPropertyNames)) then
 do i=1,size({$class->{'name'}}FloatRank0MetaPropertyNames)
  if ({$class->{'name'}}FloatRank0MetaPropertyEvolvable(i).and..not.nodeAnalytics({$offsetName}(i)).and.(propertyType == propertyTypeAll .or. (propertyType == propertyTypeActive .and. .not.nodeInactives({$offsetName}(i))) .or. (propertyType == propertyTypeInactive .and. nodeInactives({$offsetName}(i))) .or. propertyType == propertyTypeNumerics)) then
   array(offset)=self%floatRank0MetaProperties(i)
   offset=offset+1
  end if
 end do
end if
CODE
    }
    # Iterate over non-virtual, evolvable properties.
    foreach $code::property ( &Galacticus::Build::Components::Implementations::Utils::listRealEvolvers($code::member) ) {
	$code::offsetName = &offsetName('all',$code::class,$code::member,$code::property);
	if ( $code::property->{'data'}->{'rank'} == 0 ) {
	    if ( $code::property->{'data'}->{'type'} eq "double" ) {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');	
if (.not.nodeAnalytics({$offsetName}) .and. (propertyType == propertyTypeAll .or. (propertyType == propertyTypeActive .and. .not.nodeInactives({$offsetName})) .or. (propertyType == propertyTypeInactive .and. nodeInactives({$offsetName})) .or. propertyType == propertyTypeNumerics)) then
 array(offset)=self%{$property->{'name'}}Data
 offset=offset+1
end if
CODE
	    } else {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');	
count=self%{$property->{'name'}}Data%serializeCount()
if (all(.not.nodeAnalytics({$offsetName}:{$offsetName}+count-1)) .and. (propertyType == propertyTypeAll .or. (propertyType == propertyTypeActive .and. all(.not.nodeInactives({$offsetName}:{$offsetName}+count-1))) .or. (propertyType == propertyTypeInactive .and. all(nodeInactives({$offsetName}:{$offsetName}+count-1))) .or. propertyType == propertyTypeNumerics)) then
 if (count > 0) call self%{$property->{'name'}}Data%serialize(array(offset:offset+count-1))
 offset=offset+count
end if
CODE
	    }
	} else {
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');	
if (allocated(self%{$property->{'name'}}Data)) then
   count=size(self%{$property->{'name'}}Data)
   if (all(.not.nodeAnalytics({$offsetName}:{$offsetName}+count-1)) .and. (propertyType == propertyTypeAll .or. (propertyType == propertyTypeActive .and. all(.not.nodeInactives({$offsetName}:{$offsetName}+count-1))) .or. (propertyType == propertyTypeInactive .and. all(nodeInactives({$offsetName}:{$offsetName}+count-1))) .or. propertyType == propertyTypeNumerics)) then
    array(offset:offset+count-1)=reshape(self%{$property->{'name'}}Data,[count])
    offset=offset+count
   end if
end if
CODE
	}
    }
    # Insert a type-binding for this function into the implementation type.
    push(
	@{$build->{'types'}->{$implementationTypeName}->{'boundFunctions'}},
	{
	    type        => "procedure",
	    descriptor  => $function,
	    name        => "serializeValues"
	}
	);	    
}

sub Implementation_ODE_Deserialize_Values {
    # Generate a function to deserialize values of evolvable properties of each component implementation from array for the ODE solver.
    my $build     = shift();
    $code::class  = shift();
    $code::member = shift();
    my $implementationTypeName = "nodeComponent".ucfirst($code::class->{'name'}).ucfirst($code::member->{'name'});
    my $function =
    {
	type        => "void",
	name        => $implementationTypeName."DeserializeValues",
	description => "Deserialize evolvable properties of a ".$code::member->{'name'}." implementation of the ".$code::class->{'name'}." component from array.",
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => $implementationTypeName,
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "double precision",
		 attributes => [ "intent(in   )", "dimension(:)" ],
		 variables  => [ "array" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "propertyType" ]
	     }
	    ]
    };
    # Conditionally add "offset", "count" and indexing variables if they will be needed.
    my @requiredVariables = ( "offset" );
    push(@requiredVariables,"count")
	if
	(
	 exists($code::member->{'extends'})
	 ||
	 &Galacticus::Build::Components::Implementations::Utils::hasRealNonTrivialEvolvers($code::member)
	);
    push(@{$function->{'variables'}},
	 {
	     intrinsic  => "integer",
	     variables  => \@requiredVariables
	 }
	)
	if ( scalar(@requiredVariables) > 0 );
    push(@{$function->{'variables'}},
	 {
	     intrinsic  => "integer",
	     variables  => [ "i" ]
	 }
	)
	if ( ! exists($code::member->{'extends'}) );
    # Determine if the function arguments are unused.
    @code::unused = 
	(
	 exists($code::member->{'extends'})
	 ||
	 &Galacticus::Build::Components::Implementations::Utils::hasRealEvolvers          ($code::member)
	)
	?
	()
	:
	("array","propertyType");
    # Build the function.
    $function->{'content'} = "";
    if ( scalar(@code::unused) > 0 ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
!$GLC attributes unused :: {join(",",@unused)}
CODE
    }
    # Initialize offset if required.
    if ( grep {$_ eq "offset"} @requiredVariables ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');	
offset=1
CODE
    }    
    # If this component is an extension, call deserialization on the extended type.
    if ( exists($code::member->{'extends'}) ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');	
count=self%nodeComponent{ucfirst($code::member->{'extends'}->{'class'}).ucfirst($code::member->{'extends'}->{'name'})}%serializeCount   (      propertyType)
if (count > 0) then
 call self%nodeComponent{ucfirst($code::member->{'extends'}->{'class'}).ucfirst($code::member->{'extends'}->{'name'})}%deserializeValues(array,propertyType)
 offset=offset+count
end if
CODE
    } elsif ( grep {$code::class->{'name'} eq $_} @{$build->{'componentClassListActive'}} ) {
	# For non-extended types serialize meta-properties.
	$code::offsetName = &offsetName('all',$code::class->{'name'},'floatRank0MetaProperties');
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');	
if (allocated({$class->{'name'}}FloatRank0MetaPropertyNames)) then
 do i=1,size({$class->{'name'}}FloatRank0MetaPropertyNames)
  if ({$class->{'name'}}FloatRank0MetaPropertyEvolvable(i).and..not.nodeAnalytics({$offsetName}(i)).and.(propertyType == propertyTypeAll .or. (propertyType == propertyTypeActive .and. .not.nodeInactives({$offsetName}(i))) .or. (propertyType == propertyTypeInactive .and. nodeInactives({$offsetName}(i))) .or. propertyType == propertyTypeNumerics)) then
   self%floatRank0MetaProperties(i)=array(offset)
   offset=offset+1
  end if
 end do
end if
CODE
    }
    # Iterate over non-virtual, evolvable properties.
    foreach $code::property ( &Galacticus::Build::Components::Implementations::Utils::listRealEvolvers($code::member) ) {
	$code::offsetName = &offsetName('all',$code::class,$code::member,$code::property);
	if ( $code::property->{'data'}->{'rank'} == 0 ) {
	    if ( $code::property->{'data'}->{'type'} eq "double" ) {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');	
if (.not.nodeAnalytics({$offsetName}) .and. (propertyType == propertyTypeAll .or. (propertyType == propertyTypeActive .and. .not.nodeInactives({$offsetName})) .or. (propertyType == propertyTypeInactive .and. nodeInactives({$offsetName})) .or. propertyType == propertyTypeNumerics)) then
 self%{$property->{'name'}}Data=array(offset)
 offset=offset+1
end if
CODE
	    } else {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');	
count=self%{$property->{'name'}}Data%serializeCount()
if (all(.not.nodeAnalytics({$offsetName}:{$offsetName}+count-1)) .and. (propertyType == propertyTypeAll .or. (propertyType == propertyTypeActive .and. all(.not.nodeInactives({$offsetName}:{$offsetName}+count-1))) .or. (propertyType == propertyTypeInactive .and. all(nodeInactives({$offsetName}:{$offsetName}+count-1))) .or. propertyType == propertyTypeNumerics)) then
 if (count > 0) call self%{$property->{'name'}}Data%deserialize(array(offset:offset+count-1))
 offset=offset+count
end if
CODE
	    }
	} else {
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');	
if (allocated(self%{$property->{'name'}}Data)) then
   count=size(self%{$property->{'name'}}Data)
   if (all(.not.nodeAnalytics({$offsetName}:{$offsetName}+count-1)) .and. (propertyType == propertyTypeAll .or. (propertyType == propertyTypeActive .and. all(.not.nodeInactives({$offsetName}:{$offsetName}+count-1))) .or. (propertyType == propertyTypeInactive .and. all(nodeInactives({$offsetName}:{$offsetName}+count-1))) .or. propertyType == propertyTypeNumerics)) then
    self%{$property->{'name'}}Data=reshape(array(offset:offset+count-1),shape(self%{$property->{'name'}}Data))
    offset=offset+count
   end if
end if
CODE
	}
    }
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{$implementationTypeName}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    descriptor  => $function,
	    name        => "deserializeValues"
	}
	);	    
}

sub Implementation_ODE_Serialize_NonNegative {
    # Generate a function to serialize non-negative status of evolvable properties of each component implementation to array.
    my $build     = shift();
    $code::class  = shift();
    $code::member = shift();
    my $implementationTypeName = "nodeComponent".ucfirst($code::class->{'name'}).ucfirst($code::member->{'name'});
    my $function =
    {
	type        => "void",
	name        => $implementationTypeName."SerializeNonNegative",
	description => "Serialize non-negative status of evolvable properties of a {\\normalfont \\ttfamily ".$code::member->{'name'}."} implementation of the {\\normalfont \\ttfamily ".$code::class->{'name'}."} component to array.",
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => $implementationTypeName,
		 attributes => [ "intent(in   )" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "logical",
		 attributes => [ "intent(  out)", "dimension(:)" ],
		 variables  => [ "array" ]
	     }
	    ]
    };
    # Conditionally add "offset", "count" and indexing variables if they will be needed.
    my @requiredVariables = ( "offset" );
    push(@requiredVariables,"count")
	if
	(
	 exists($code::member->{'extends'})
	 ||
	 &Galacticus::Build::Components::Implementations::Utils::hasRealNonTrivialEvolvers($code::member)
	);
    push(@requiredVariables,"i")
	if ( ! exists($code::member->{'extends'}) );
    push(@{$function->{'variables'}},
	 {
	     intrinsic  => "integer",
	     variables  => \@requiredVariables
	 }
	)
	if ( scalar(@requiredVariables) > 0 );
    # Determine if the function arguments are unused.
    @code::unused = 
	(
	 exists($code::member->{'extends'})
	 ||
	 &Galacticus::Build::Components::Implementations::Utils::hasRealEvolvers          ($code::member)
	)
	?
	()
	:
	("array");
    # Build the function.
    $function->{'content'} = "";
    if ( scalar(@code::unused) > 0 ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
!$GLC attributes unused :: {join(",",@unused)}
CODE
    }
    # Initialize offset if required.
    if ( grep {$_ eq "offset"} @requiredVariables ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');	
offset=1
CODE
    }    
    # If this component is an extension, call serialization on the extended type.
    if ( exists($code::member->{'extends'}) ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');	
count=self%nodeComponent{ucfirst($code::member->{'extends'}->{'class'}).ucfirst($code::member->{'extends'}->{'name'})}%serializeCount(propertyTypeActive)
if (count > 0) then
 call self%nodeComponent{ucfirst($code::member->{'extends'}->{'class'}).ucfirst($code::member->{'extends'}->{'name'})}%serializeNonNegative(array)
 offset=offset+count
end if
CODE
    } elsif ( grep {$code::class->{'name'} eq $_} @{$build->{'componentClassListActive'}} ) {
	# For non-extended types serialize meta-properties.
	$code::offsetName = &offsetName('all',$code::class->{'name'},'floatRank0MetaProperties');
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');	
if (allocated({$class->{'name'}}FloatRank0MetaPropertyNames)) then
 do i=1,size({$class->{'name'}}FloatRank0MetaPropertyNames)
  if ({$class->{'name'}}FloatRank0MetaPropertyEvolvable(i).and..not.nodeAnalytics({$offsetName}(i)) .and. .not.nodeInactives({$offsetName}(i))) then
   array(offset)=.false. ! Currently no support for non-negative meta-properties.
   offset=offset+1
  end if
 end do
end if
CODE
    }
    # Iterate over non-virtual, evolvable properties.
    foreach $code::property ( &Galacticus::Build::Components::Implementations::Utils::listRealEvolvers($code::member) ) {
	$code::offsetName    = &offsetName('all',$code::class,$code::member,$code::property);
	$code::isNonNegative = $code::property->{'attributes'}->{'isNonNegative'} ? ".true." : ".false.";
	if ( $code::property->{'data'}->{'rank'} == 0 ) {
	    if ( $code::property->{'data'}->{'type'} eq "double" ) {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');	
if (.not.nodeAnalytics({$offsetName}) .and. .not.nodeInactives({$offsetName})) then
 array(offset)={$isNonNegative}
 offset=offset+1
end if
CODE
	    } else {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');	
count=self%{$property->{'name'}}Data%serializeCount()
if (all(.not.nodeAnalytics({$offsetName}:{$offsetName}+count-1)) .and. all(.not.nodeInactives({$offsetName}:{$offsetName}+count-1))) then
 if (count > 0) array(offset:offset+count-1)={$isNonNegative}
 offset=offset+count
end if
CODE
	    }
	} else {
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');	
if (allocated(self%{$property->{'name'}}Data)) then
   count=size(self%{$property->{'name'}}Data)
   if (all(.not.nodeAnalytics({$offsetName}:{$offsetName}+count-1)) .and. all(.not.nodeInactives({$offsetName}:{$offsetName}+count-1))) then
    array(offset:offset+count-1)={$isNonNegative}
    offset=offset+count
   end if
end if
CODE
	}
    }
    # Insert a type-binding for this function into the implementation type.
    push(
	@{$build->{'types'}->{$implementationTypeName}->{'boundFunctions'}},
	{
	    type        => "procedure",
	    descriptor  => $function,
	    name        => "serializeNonNegative"
	}
	);	    
}

sub Implementation_ODE_Offsets {
    # Generate function to compute offsets into serialization arrays for component implementations.
    my $build     = shift();
    $code::class  = shift();
    $code::member = shift();
    my $implementationTypeName = "nodeComponent".ucfirst($code::class->{'name'}).ucfirst($code::member->{'name'});
    my $function =
    {
	type        => "void",
	name        => $implementationTypeName."SerializeOffsets",
	description => "Compute offsets into serialization arrays for a {\\normalfont \\ttfamily ".$code::member->{'name'}."} implementation of the {\\normalfont \\ttfamily ".$code::class->{'name'}."} component.",
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => $implementationTypeName,
		 attributes => [ "intent(in   )" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => [ "intent(inout)" ],
		 variables  => [ "count", "countSubset" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "propertyType" ]
	     }
	    ]
    };
    # Add an indexing variable for non-extended types needed for meta-property serialization.
    push(@{$function->{'variables'}},
	 {
	     intrinsic  => "integer",
	     variables  => [ "i" ]
	 }
	)
	if ( ! exists($code::member->{'extends'}) );
    # Determine if the function arguments are unused.
    @code::unused = ();
    push(@code::unused,"self")
	unless 
	(
	 exists($code::member->{'extends'})
	 ||
	 &Galacticus::Build::Components::Implementations::Utils::hasRealNonTrivialEvolvers($code::member)
	);
    push(@code::unused,"count","countSubset","propertyType")
	unless 
	(
	 exists($code::member->{'extends'})
	 ||
	 &Galacticus::Build::Components::Implementations::Utils::hasRealEvolvers           ($code::member)
	);
    push(@{$function->{'variables'}},{intrinsic => "integer", variables => [ "propertySize"] })
	if 
	(
	 &Galacticus::Build::Components::Implementations::Utils::hasRealEvolvers           ($code::member)
	);
    # Build the function.
    $function->{'content'} = "";
    if ( scalar(@code::unused) > 0 ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
!$GLC attributes unused :: {join(",",@unused)}
CODE
    }
    # If this component is an extension, compute offsets of the extended type.
    if ( exists($code::member->{'extends'}) ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
call self%nodeComponent{ucfirst($code::member->{'extends'}->{'class'}).ucfirst($code::member->{'extends'}->{'name'})}%serializationOffsets(count,countSubset,propertyType)
CODE
    } elsif ( grep {$code::class->{'name'} eq $_} @{$build->{'componentClassListActive'}} ) {
	# For non-extended types compute offsets for meta-properties.
	# Set the offset for this property to the current count plus 1 (since we haven't yet updated the count. 
	$code::offsetNameAll      = &offsetName('all'     ,$code::class->{'name'},'floatRank0MetaProperties');
	$code::offsetNameActive   = &offsetName('active'  ,$code::class->{'name'},'floatRank0MetaProperties');
	$code::offsetNameInactive = &offsetName('inactive',$code::class->{'name'},'floatRank0MetaProperties');
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');	
if (allocated({$class->{'name'}}FloatRank0MetaPropertyNames)) then
CODE
	foreach my $status ( "all", "active", "inactive" ) {
	    $code::offsetName = &offsetName($status,$code::class->{'name'},'floatRank0MetaProperties');
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');	
 if (.not.allocated({$offsetName})) then
  allocate({$offsetName}({$class->{'name'}}FloatRank0MetaPropertyCount))
 else if (size({$offsetName}) /= {$class->{'name'}}FloatRank0MetaPropertyCount) then
  deallocate({$offsetName})
  allocate({$offsetName}({$class->{'name'}}FloatRank0MetaPropertyCount))
 end if
CODE
	}
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');	
 do i=1,size({$class->{'name'}}FloatRank0MetaPropertyNames)
  if (.not.{$class->{'name'}}FloatRank0MetaPropertyEvolvable(i)) cycle
  if      (propertyType == propertyTypeAll     ) then
                                    {$offsetNameAll}     (i)=count      +1
  else if (propertyType == propertyTypeInactive) then
   if (     nodeInactives(count+1).and..not.nodeAnalytics(count+1)) {$offsetNameInactive}(i)=countSubset+1
  else if (propertyType == propertyTypeActive  ) then
   if (.not.nodeInactives(count+1).and..not.nodeAnalytics(count+1)) {$offsetNameActive}  (i)=countSubset+1
  else if (propertyType == propertyTypeNumerics) then
   if (                                .not.nodeAnalytics(count+1)) {$offsetNameActive}  (i)=countSubset+1
  end if
  if (propertyType /= propertyTypeAll) then
   if (((nodeInactives(count+1) .and. propertyType == propertyTypeInactive) .or. (.not.nodeInactives(count+1) .and. propertyType == propertyTypeActive) .or. propertyType == propertyTypeNumerics) .and. .not.nodeAnalytics(count+1)) countSubset=countSubset+1
  end if
  count=count+1
 end do
end if
CODE
    }
    # Iterate over non-virtual, evolvable properties.
    foreach $code::property ( &Galacticus::Build::Components::Implementations::Utils::listRealEvolvers($code::member) ) {
	# Set the offset for this property to the current count plus 1 (since we haven't yet updated the count. 
	$code::offsetNameAll      = &offsetName('all'     ,$code::class,$code::member,$code::property);
	$code::offsetNameActive   = &offsetName('active'  ,$code::class,$code::member,$code::property);
	$code::offsetNameInactive = &offsetName('inactive',$code::class,$code::member,$code::property);
	# Get the size of this property.
	if ( $code::property->{'data'}->{'rank'} == 0 ) {
	    if ( $code::property->{'data'}->{'type'} eq "double" ) {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');	
propertySize=1
CODE
	    } else {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');	
propertySize=self%{$property->{'name'}}Data%serializeCount()
CODE
	    }
        } else {
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');	
propertySize=size(self%{$property->{'name'}}Data)
CODE
        }
	# Set the offset into this property, but only if the property has non-zero size.
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');	
if (propertySize > 0) then
 if      (propertyType == propertyTypeAll     ) then
                                   {$offsetNameAll}     =count      +1
 else if (propertyType == propertyTypeInactive) then
  if (     nodeInactives(count+1).and..not.nodeAnalytics(count+1)) {$offsetNameInactive}=countSubset+1
 else if (propertyType == propertyTypeActive  ) then
  if (.not.nodeInactives(count+1).and..not.nodeAnalytics(count+1)) {$offsetNameActive}  =countSubset+1
 else if (propertyType == propertyTypeNumerics) then
  if (                                .not.nodeAnalytics(count+1)) {$offsetNameActive}  =countSubset+1
 end if
CODE
	# Update the counts by the size of this property.
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');	
 if (propertyType /= propertyTypeAll) then
  if (((nodeInactives(count+1) .and. propertyType == propertyTypeInactive) .or. (.not.nodeInactives(count+1) .and. propertyType == propertyTypeActive) .or. propertyType == propertyTypeNumerics) .and. .not.nodeAnalytics(count+1)) countSubset=countSubset+propertySize
 end if
 count=count+propertySize
end if
CODE
    }
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{$implementationTypeName}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    descriptor  => $function,
	    name        => "serializationOffsets"
	}
	);
}

sub Implementation_ODE_Offset_Variables {
    # Generate variables which store offsets into the ODE solver arrays.
    my $build  = shift();
    my $class  = shift();
    my $member = shift();
    # Include meta-properties for just the null class (since we need only one copy of these per-class).
    if ( grep {$class->{'name'} eq $_} @{$build->{'componentClassListActive'}} ) {
	if ( $member->{'name'} eq "null" ) {
	    foreach my $status ( "all", "active", "inactive" ) {
		my $offsetName = &offsetName($status,$class->{'name'},'floatRank0MetaProperties');
		push(
		    @{$build->{'variables'}},
		    {
			intrinsic  => "integer",
			ompPrivate => 1,
			attributes => [ "allocatable", "dimension(:)" ],
			variables  => [ $offsetName ]
		    }
		    );
	    }
	}
    }
    # Iterate over non-virtual, evolving properties.
    foreach my $property ( &Galacticus::Build::Components::Implementations::Utils::listRealEvolvers($member) ) {
	foreach my $status ( "all", "active", "inactive" ) {
	    my $offsetName = &offsetName($status,$class->{'name'}.$member->{'name'},$property->{'name'});
	    push(
		@{$build->{'variables'}},
		{
		    intrinsic  => "integer",
		    ompPrivate => 1,
		    variables  => [ $offsetName ]
		}
		);
	}
    }
}

sub Implementation_ODE_Rate_Variables {
    # Generate variables which store ODE solver variable rates and scales.
    my $build = shift();
    push(
	@{$build->{'variables'}},
	{
	    intrinsic  => "integer",
	    ompPrivate => 1,
	    variables  => [ "nodeSerializationCount", "nodeSerializationCountActive", "nodeSerializationCountInactive" ]
	},
	{
	    intrinsic  => "double precision",
	    attributes => [ "allocatable", "dimension(:)" ],
	    ompPrivate => 1,
	    variables  => [ "nodeScales", "nodeRates", "nodeRatesActives" ]
	},
	{
	    intrinsic  => "logical",
	    attributes => [ "allocatable", "dimension(:)" ],
	    ompPrivate => 1,
	    variables  => [ "nodeInactives", "nodeAnalytics" ]
	}
	);
}

1;
