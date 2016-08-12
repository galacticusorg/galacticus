# Contains a Perl module which provides various ODE solver-related functions for component implementations.

package ODESolver;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V094"}) ) {
    $galacticusPath = $ENV{"GALACTICUS_ROOT_V094"};
    $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
    $galacticusPath = "./";
}
unshift(@INC, $galacticusPath."perl"); 
use strict;
use warnings;
use utf8;
use Text::Template 'fill_in_string';
use Data::Dumper;
require List::ExtraUtils;
require Galacticus::Build::Components::Utils;
require Galacticus::Build::Components::DataTypes;
require Galacticus::Build::Components::Implementations::Utils;

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     implementationODESolver =>
     {
	 functions =>
	     [
	      \&Implementation_ODE_Serialize_Count   ,
	      \&Implementation_ODE_Serialize_Values  ,
	      \&Implementation_ODE_Deserialize_Values,
	      \&Implementation_ODE_Name_From_Index
	     ]
     }
    );

sub Implementation_ODE_Name_From_Index {
    # Generate a function to return the name of a property given the index of that property in the serialization of a component
    # implementation.
    my $build = shift();
    # Iterate over component classes.
    foreach $code::class ( &ExtraUtils::hashList($build->{'componentClasses'}) ) {
	# Iterate over class member implementations.
	foreach $code::member ( @{$code::class->{'members'}} ) {
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
	    push(@code::unused,"count")
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
name=self%nodeComponent{ucfirst($member->{'extends'}->{'class'}).ucfirst($member->{'extends'}->{'name'})}%nameFromIndex(count)
if (count <= 0) return
CODE
	    }
	    # Iterate over properties.
	    foreach $code::property ( &ExtraUtils::hashList($code::member->{'properties'}->{'property'}) ) {
		# Only evolvable, non-virtual properties are included in the ODE solver.
		next
		    unless
		    (
		     ! $code::property->{'attributes'}->{'isVirtual'  }
		     &&
		       $code::property->{'data'      }->{'isEvolvable'}
		    );
		# Find condition for count update. For allocatable properties, condition is that the object be allocated. For
		# non-allocatable properties, always update the count.
		$code::condition = 
		    $code::property->{'data'}->{'rank'} > 0 
		    ? 
		    "if (allocated(self%".$code::property->{'name'}."Data)) "
		    :
		    "";
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
{$condition}count=count-{$count}
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
		    name        => "nameFromIndex", 
		    returnType  => "\\textcolor{red}{\\textless varying\\_string\\textgreater}", 
		    arguments   => "\\intzero\\ index\\argin"
		}
		);	    
	}
    }
}

sub Implementation_ODE_Serialize_Count {
    # Generate a function to return a count of serializable, evolvable properties of each component implementation for the ODE solver.
    my $build = shift();
    # Iterate over component classes.
    foreach $code::class ( &ExtraUtils::hashList($build->{'componentClasses'}) ) {
	# Iterate over class member implementations.
	foreach $code::member ( @{$code::class->{'members'}} ) {
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
		    grep
		    {	
			! $_->{'attributes'}->{'isVirtual'  }
			&&	    
			  $_->{'data'      }->{'isEvolvable'}
			&&
			    (
			     $_->{'data'}->{'rank'       } >  0 
			     ||
			     $_->{'data'}->{'type'       } ne "double"
			    )			
		    }
		    &ExtraUtils::hashList($code::member->{'properties'}->{'property'})
		);
	    # Build the function.
	    $function->{'content'}  = "";
	    if ( scalar(@code::unused) > 0 ) {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
!GCC$ attributes unused :: {join(",",@unused)}
CODE
	    }
	    # If this component is an extension, first call on the extended type.
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
{$implementationTypeName}SerializeCount={exists($code::member->{'extends'}) ? "self%nodeComponent".ucfirst($extends->{'class'}).ucfirst($extends->{'name'})."%serializeCount()" : "0"}
CODE
	    # Initialize count of fixed, scalar properties to zero.
	    $code::scalarPropertyCount = 0;
	    # Iterate over properties.
	    foreach $code::property ( &ExtraUtils::hashList($code::member->{'properties'}->{'property'}) ) {
		# Only evolvable, non-virtual properties are included in the ODE solver.
		next
		    unless
		    (
		     ! $code::property->{'attributes'}->{'isVirtual'  }
		     &&
		       $code::property->{'data'      }->{'isEvolvable'}
		    );
		if ( $code::property->{'data'}->{'rank'} == 0 && $code::property->{'data'}->{'type'} eq "double" ) {
		    ++$code::scalarPropertyCount;
		} else {		
		    # Find condition for count update. For allocatable properties, condition is that the object be allocated. For
		    # non-allocatable properties, always update the count.
		    $code::condition = 
			$code::property->{'data'}->{'rank'} > 0 
			? 
			"if (allocated(self%".$code::property->{'name'}."Data)) "
			:
			"";
		    # Find the size of the object. Rank-0 double properties always have a count of 1. For other rank-0 types, call
		    # their serialization count method. Rank>0 must be double, so simply use the array length.
		    $code::count     = 
			$code::property->{'data'}->{'rank'} > 0 
			?
			"size(self%".$code::property->{'name'}."Data)"
			:
			"self%".$code::property->{'name'}."Data%serializeCount()";
		    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
{$condition}{$implementationTypeName}SerializeCount={$implementationTypeName}SerializeCount+{$count}
CODE
		}
	    }
	    # Add count of scalar properties.
	    if ( $code::scalarPropertyCount > 0 ) {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
{$implementationTypeName}SerializeCount={$implementationTypeName}SerializeCount+{$scalarPropertyCount}
CODE
	    }
	    # Insert a type-binding for this function into the treeNode type.
	    push(
		@{$build->{'types'}->{$code::implementationTypeName}->{'boundFunctions'}},
		{
		    type        => "procedure", 
		    descriptor  => $function,
		    name        => "serializeCount", 
		    returnType  => "\\intzero", 
		    arguments   => ""
		}
		);	    
	}
    }
}

sub Implementation_ODE_Serialize_Values {
    # Generate a function to serialize values of evolvable properties of each component implementation to array for the ODE solver.
    my $build = shift();
    # Iterate over component classes.
    foreach $code::class ( &ExtraUtils::hashList($build->{'componentClasses'}) ) {
	# Iterate over class member implementations.
	foreach $code::member ( @{$code::class->{'members'}} ) {
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
		("self","array");
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
count=self%nodeComponent{ucfirst($code::member->{'extends'}->{'class'}).ucfirst($code::member->{'extends'}->{'name'})}%serializeCount (     )
if (count > 0) then
 call self%nodeComponent{ucfirst($code::member->{'extends'}->{'class'}).ucfirst($code::member->{'extends'}->{'name'})}%serializeValues(array)
 offset=offset+count
end if
CODE
	    }
	    # Iterate over properties.
	    foreach $code::property ( &ExtraUtils::hashList($code::member->{'properties'}->{'property'}) ) {
		# Only evolvable, non-virtual properties are included in the ODE solver.
		next
		    unless
		    (
		     ! $code::property->{'attributes'}->{'isVirtual'  }
		     &&
		       $code::property->{'data'      }->{'isEvolvable'}
		    );
		if ( $code::property->{'data'}->{'rank'} == 0 ) {
		    if ( $code::property->{'data'}->{'type'} eq "double" ) {
			$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');	
array(offset)=self%{$property->{'name'}}Data
offset=offset+1
CODE
		    } else {
			$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');	
count=self%{$property->{'name'}}Data%serializeCount()
if (count > 0) call  self%{$property->{'name'}}Data%serialize(array(offset:offset+count-1))
offset=offset+count
CODE
		    }
		} else {
		    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');	
if (allocated(self%{$property->{'name'}}Data)) then
   count=size(self%{$property->{'name'}}Data)
   array(offset:offset+count-1)=reshape(self%{$property->{'name'}}Data,[count])
   offset=offset+count
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
		    name        => "serializeValues", 
		    returnType  => "\\void", 
		    arguments   => "\\doubleone array(:)\\argout"
		}
		);	    
	}
    }
}

sub Implementation_ODE_Deserialize_Values {
    # Generate a function to deserialize values of evolvable properties of each component implementation from array for the ODE solver.
    my $build = shift();
    # Iterate over component classes.
    foreach $code::class ( &ExtraUtils::hashList($build->{'componentClasses'}) ) {
	# Iterate over class member implementations.
	foreach $code::member ( @{$code::class->{'members'}} ) {
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
		("self","array");
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
count=self%nodeComponent{ucfirst($code::member->{'extends'}->{'class'}).ucfirst($code::member->{'extends'}->{'name'})}%serializeCount   (     )
if (count > 0) then
 call self%nodeComponent{ucfirst($code::member->{'extends'}->{'class'}).ucfirst($code::member->{'extends'}->{'name'})}%deserializeValues(array)
 offset=offset+count
end if
CODE
	    }
	    # Iterate over properties.
	    foreach $code::property ( &ExtraUtils::hashList($code::member->{'properties'}->{'property'}) ) {
		# Only evolvable, non-virtual properties are included in the ODE solver.
		next
		    unless
		    (
		     ! $code::property->{'attributes'}->{'isVirtual'  }
		     &&
		       $code::property->{'data'      }->{'isEvolvable'}
		    );
		if ( $code::property->{'data'}->{'rank'} == 0 ) {
		    if ( $code::property->{'data'}->{'type'} eq "double" ) {
			$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');	
self%{$property->{'name'}}Data=array(offset)
offset=offset+1
CODE
		    } else {
			$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');	
count=self%{$property->{'name'}}Data%serializeCount()
if (count > 0) call self%{$property->{'name'}}Data%deserialize(array(offset:offset+count-1))
offset=offset+count
CODE
		    }
		} else {
		    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');	
if (allocated(self%{$property->{'name'}}Data)) then
   count=size(self%{$property->{'name'}}Data)
   self%{$property->{'name'}}Data=reshape(array(offset:offset+count-1),shape(self%{$property->{'name'}}Data))
   offset=offset+count
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
		    name        => "deserializeValues", 
		    returnType  => "\\void", 
		    arguments   => "\\doubleone array(:)\\argin"
		}
		);	    
	}
    }
}

1;
