# Contains a Perl module which provides various ODE solver-related functions for component implementations.

package Galacticus::Build::Components::Implementations::ODESolver;
use strict;
use warnings;
use utf8;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
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
	      \&Implementation_ODE_Serialize_Count   ,
	      \&Implementation_ODE_Serialize_Values  ,
	      \&Implementation_ODE_Deserialize_Values,
	      \&Implementation_ODE_Name_From_Index   ,
	      \&Implementation_ODE_Offsets           ,
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
    # property is evolveable.
    push(@code::unused,"count","propertyType")
	unless (
	    exists($code::member->{'extends'})
	    ||
	    &Galacticus::Build::Components::Implementations::Utils::hasRealEvolvers          ($code::member)
	);
    # Build the function.
    $function->{'content'}  = "";
    if ( scalar(@code::unused) > 0 ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
!GCC$ attributes unused :: {join(",",@unused)}
CODE
    }
    # If this component is an extension, first call on the extended type.
    if ( exists($code::member->{'extends'}) ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
name=self%nodeComponent{ucfirst($member->{'extends'}->{'class'}).ucfirst($member->{'extends'}->{'name'})}%nameFromIndex(count,propertyType)
if (count <= 0) return
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
 if (                                                                                             &
  &    propertyType == propertyTypeAll                                                            &
  &  .or.                                                                                         &
  &   (propertyType == propertyTypeActive   .and. all(.not.nodeInactives({$offsetName}:{$offsetName}+{$count}-1))) &
  &  .or.                                                                                         &
  &   (propertyType == propertyTypeInactive .and. all(     nodeInactives({$offsetName}:{$offsetName}+{$count}-1))) &
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
	if ( &Galacticus::Build::Components::Implementations::Utils::hasRealEvolvers($code::member) );
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
!GCC$ attributes unused :: {join(",",@unused)}
CODE
    }
    # If this component is an extension, first call on the extended type.
    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
{$implementationTypeName}SerializeCount={exists($member->{'extends'}) ? "self%nodeComponent".ucfirst($member->{'extends'}->{'class'}).ucfirst($member->{'extends'}->{'name'})."%serializeCount(propertyType)" : "0"}
CODE
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
else if ((propertyType == propertyTypeActive .and. all(.not.nodeInactives({$offsetName}:{$offsetName}+count-1))) .or. (propertyType == propertyTypeInactive .and. all(nodeInactives({$offsetName}:{$offsetName}+count-1)))) then
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
    # Conditionally add "offset" and "count" variables if they will be needed.
    my @requiredVariables;
    push(@requiredVariables,"count")
	if
	(
	 exists($code::member->{'extends'})
	 ||
	 &Galacticus::Build::Components::Implementations::Utils::hasRealNonTrivialEvolvers($code::member)
	);
    push(@requiredVariables,"offset")
	if
	(
	 exists($code::member->{'extends'})
	 ||
	 &Galacticus::Build::Components::Implementations::Utils::hasRealEvolvers          ($code::member)
	);	   
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
	("self","array","propertyType");
    # Build the function.
    $function->{'content'} = "";
    if ( scalar(@code::unused) > 0 ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
!GCC$ attributes unused :: {join(",",@unused)}
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
    }
    # Iterate over non-virtual, evolvable properties.
    foreach $code::property ( &Galacticus::Build::Components::Implementations::Utils::listRealEvolvers($code::member) ) {
	$code::offsetName = &offsetName('all',$code::class,$code::member,$code::property);
	if ( $code::property->{'data'}->{'rank'} == 0 ) {
	    if ( $code::property->{'data'}->{'type'} eq "double" ) {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');	
if (propertyType == propertyTypeAll .or. (propertyType == propertyTypeActive .and. .not.nodeInactives({$offsetName})) .or. (propertyType == propertyTypeInactive .and. nodeInactives({$offsetName}))) then
 array(offset)=self%{$property->{'name'}}Data
 offset=offset+1
end if
CODE
	    } else {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');	
count=self%{$property->{'name'}}Data%serializeCount()
if (propertyType == propertyTypeAll .or. (propertyType == propertyTypeActive .and. all(.not.nodeInactives({$offsetName}:{$offsetName}+count-1))) .or. (propertyType == propertyTypeInactive .and. all(nodeInactives({$offsetName}:{$offsetName}+count-1)))) then
 if (count > 0) call self%{$property->{'name'}}Data%serialize(array(offset:offset+count-1))
 offset=offset+count
end if
CODE
	    }
	} else {
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');	
if (allocated(self%{$property->{'name'}}Data)) then
   count=size(self%{$property->{'name'}}Data)
   if (propertyType == propertyTypeAll .or. (propertyType == propertyTypeActive .and. all(.not.nodeInactives({$offsetName}:{$offsetName}+count-1))) .or. (propertyType == propertyTypeInactive .and. all(nodeInactives({$offsetName}:{$offsetName}+count-1)))) then
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
    # Conditionally add "offset" and "count" variables if they will be needed.
    my @requiredVariables;
    push(@requiredVariables,"count")
	if
	(
	 exists($code::member->{'extends'})
	 ||
	 &Galacticus::Build::Components::Implementations::Utils::hasRealNonTrivialEvolvers($code::member)
	);
    push(@requiredVariables,"offset")
	if
	(
	 exists($code::member->{'extends'})
	 ||
	 &Galacticus::Build::Components::Implementations::Utils::hasRealEvolvers          ($code::member)
	);	   
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
	("self","array","propertyType");
    # Build the function.
    $function->{'content'} = "";
    if ( scalar(@code::unused) > 0 ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
!GCC$ attributes unused :: {join(",",@unused)}
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
    }
    # Iterate over non-virtual, evolvable properties.
    foreach $code::property ( &Galacticus::Build::Components::Implementations::Utils::listRealEvolvers($code::member) ) {
	$code::offsetName = &offsetName('all',$code::class,$code::member,$code::property);
	if ( $code::property->{'data'}->{'rank'} == 0 ) {
	    if ( $code::property->{'data'}->{'type'} eq "double" ) {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');	
if (propertyType == propertyTypeAll .or. (propertyType == propertyTypeActive .and. .not.nodeInactives({$offsetName})) .or. (propertyType == propertyTypeInactive .and. nodeInactives({$offsetName}))) then
 self%{$property->{'name'}}Data=array(offset)
 offset=offset+1
end if
CODE
	    } else {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');	
count=self%{$property->{'name'}}Data%serializeCount()
if (propertyType == propertyTypeAll .or. (propertyType == propertyTypeActive .and. all(.not.nodeInactives({$offsetName}:{$offsetName}+count-1))) .or. (propertyType == propertyTypeInactive .and. all(nodeInactives({$offsetName}:{$offsetName}+count-1)))) then
 if (count > 0) call self%{$property->{'name'}}Data%deserialize(array(offset:offset+count-1))
 offset=offset+count
end if
CODE
	    }
	} else {
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');	
if (allocated(self%{$property->{'name'}}Data)) then
   count=size(self%{$property->{'name'}}Data)
   if (propertyType == propertyTypeAll .or. (propertyType == propertyTypeActive .and. all(.not.nodeInactives({$offsetName}:{$offsetName}+count-1))) .or. (propertyType == propertyTypeInactive .and. all(nodeInactives({$offsetName}:{$offsetName}+count-1)))) then
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
	     },
	     {
		 intrinsic  => "integer",
		 variables  => [ "propertySize" ]
	     }
	    ]
    };
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
    push(@code::unused,"propertySize")
	unless 
	(
	 &Galacticus::Build::Components::Implementations::Utils::hasRealEvolvers           ($code::member)
	);
    # Build the function.
    $function->{'content'} = "";
    if ( scalar(@code::unused) > 0 ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
!GCC$ attributes unused :: {join(",",@unused)}
CODE
    }
    # If this component is an extension, compute offsets of the extended type.
    if ( exists($code::member->{'extends'}) ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
call self%nodeComponent{ucfirst($code::member->{'extends'}->{'class'}).ucfirst($code::member->{'extends'}->{'name'})}%serializationOffsets(count,countSubset,propertyType)
CODE
    }
    # Iterate over non-virtual, evolvable properties.
    foreach $code::property ( &Galacticus::Build::Components::Implementations::Utils::listRealEvolvers($code::member) ) {
	# Set the offset for this property to the current count plus 1 (since we haven't yet updated the count. 
	$code::offsetNameAll      = &offsetName('all'     ,$code::class,$code::member,$code::property);
	$code::offsetNameActive   = &offsetName('active'  ,$code::class,$code::member,$code::property);
	$code::offsetNameInactive = &offsetName('inactive',$code::class,$code::member,$code::property);
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');	
if      (propertyType == propertyTypeAll     ) then
                                  {$offsetNameAll}     =count      +1
else if (propertyType == propertyTypeInactive) then
 if (     nodeInactives(count+1)) {$offsetNameInactive}=countSubset+1
else if (propertyType == propertyTypeActive  ) then
 if (.not.nodeInactives(count+1)) {$offsetNameActive}  =countSubset+1
end if
CODE
	# Update the counts by the size of this property.
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
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');	
if (propertyType /= propertyTypeAll) then
 if ((nodeInactives(count+1) .and. propertyType == propertyTypeInactive) .or. (.not.nodeInactives(count+1) .and. propertyType == propertyTypeActive)) countSubset=countSubset+propertySize
end if
count=count+propertySize
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
	    variables  => [ "nodeScales", "nodeRates" ]
	},
	{
	    intrinsic  => "logical",
	    attributes => [ "allocatable", "dimension(:)" ],
	    ompPrivate => 1,
	    variables  => [ "nodeInactives" ]
	}
	);
}

1;
