# Contains a Perl module which provides various ODE solver-related functions for tree nodes.

package Galacticus::Build::Components::TreeNodes::ODESolver;
use strict;
use warnings;
use utf8;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use Text::Template 'fill_in_string';
use List::ExtraUtils;
use Galacticus::Build::Components::Utils;
use Galacticus::Build::Components::DataTypes;

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     treeNodeODESolver =>
     {
	 functions =>
	     [
	      \&Tree_Node_ODE_Step_Initialize   ,
	      \&Tree_Node_ODE_Serialize_Count   ,
	      \&Tree_Node_ODE_Serialize_Values  ,
	      \&Tree_Node_ODE_Deserialize_Values,
	      \&Tree_Node_ODE_Serialize_Rates   ,
	      \&Tree_Node_ODE_Serialize_Scales  ,
	      \&Tree_Node_ODE_Serialize_Inactive,
	      \&Tree_Node_ODE_Name_From_Index   ,
	      \&Tree_Node_ODE_Offsets
	     ]
     }
    );

sub Tree_Node_ODE_Step_Initialize {
    # Generate a function to initialize ODE solver variables.
    my $build = shift();
    my @quantities =
	(
	 {
	     name  => "rate",
	     value => "0.0d0"
	 },
	 {
	     name  => "scale",
	     value => "1.0d0"
	 },
	 {
	     name  => "inactive",
	     value => ".false."
	 }
	);
    foreach $code::quantity ( @quantities ) {
	my $function =
	{
	    type        => "void",
	    name        => "treeNodeODEStep".ucfirst($code::quantity->{'name'})."sInitialize",
	    description => "Initialize the ".$code::quantity->{'name'}."s in components of tree node {\\normalfont \\ttfamily self} in preparation for an ODE solver step.",
	    variables   =>
		[
		 {
		     intrinsic  => "class",
		     type       => "treeNode",
		     attributes => [ "intent(in   )" ],
		     variables  => [ "self" ]
		 }
		]
	};    
	$function->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
!GCC$ attributes unused :: self
node{ucfirst($quantity->{'name'})}s={$quantity->{'value'}}
CODE
	# Insert a type-binding for this function into the treeNode type.
	push(
	    @{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	    {
		type        => "procedure", 
		descriptor  => $function,
		name        => "odeStep".ucfirst($code::quantity->{'name'})."sInitialize"
	    }
	    );
    }
}

sub Tree_Node_ODE_Serialize_Count {
    # Generate a function to return a count of serializable, evolvable properties of a node for the ODE solver.
    my $build = shift();
    my $function =
    {
	type        => "integer",
	name        => "treeNodeSerializeCount",
	description => "Return a count of the size of the node when serialized to an array.",
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => "treeNode",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "propertyType" ]
	     },
	     {
		 intrinsic  => "integer",
		 variables  => [ "i" ],
	     }
	    ]
    };    
    $function->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
treeNodeSerializeCount=0
CODE
    # Iterate over all component classes
    foreach $code::class ( &List::ExtraUtils::hashList($build->{'componentClasses'}) ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if (allocated(self%component{ucfirst($class->{'name'})})) then
  do i=1,size(self%component{ucfirst($class->{'name'})})
    treeNodeSerializeCount=treeNodeSerializeCount+self%component{ucfirst($class->{'name'})}(i)%serializeCount(propertyType)
  end do
end if
CODE
    }
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    descriptor  => $function,
	    name        => "serializeCount"
	}
	);
}

sub Tree_Node_ODE_Serialize_Values {
    # Generate a function to serialize evolvable properties of a node into an array for the ODE solver.
    my $build = shift();
    my $function =
    {
	type        => "void",
	name        => "treeNodeSerializeValuesToArray",
	description => "Serialize evolvable properties of a node into an array.",
	modules     =>
	    [
	     "Memory_Management"
	    ],
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => "treeNode",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "double precision",
		 attributes => [ "dimension(:)", "intent(  out)" ],
		 variables  => [ "array" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "propertyType" ],
	     },
	     {
		 intrinsic  => "integer",
		 variables  => [ "count", "offset", "i" ],
	     }
	    ]
    };    
    $function->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
offset=1
CODE
    # Iterate over all component classes
    foreach $code::class ( &List::ExtraUtils::hashList($build->{'componentClasses'}) ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if (allocated(self%component{ucfirst($class->{'name'})})) then
  do i=1,size(self%component{ucfirst($class->{'name'})})
    count=self%component{ucfirst($class->{'name'})}(i)%serializeCount(propertyType)
    if (count > 0) call self%component{ucfirst($class->{'name'})}(i)%serializeValues(array(offset:),propertyType)
    offset=offset+count
  end do
end if
CODE
    }
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    descriptor  => $function,
	    name        => "serializeValues"
	}
	);
}

sub Tree_Node_ODE_Deserialize_Values {
    # Generate a function to deserialize evolvable properties of a node from an array for the ODE solver.
    my $build = shift();
    my $function =
    {
	type        => "void",
	name        => "treeNodeDeserializeValuesFromArray",
	description => "Deserialize evolvable properties of a node from an array.",
	modules     =>
	    [
	     "Memory_Management"
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
		 intrinsic  => "double precision",
		 attributes => [ "dimension(:)", "intent(in   )" ],
		 variables  => [ "array" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "propertyType" ],
	     },
	     {
		 intrinsic  => "integer",
		 variables  => [ "count", "offset", "i" ],
	     }
	    ]
    };    
    $function->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
offset=1
CODE
    # Iterate over all component classes
    foreach $code::class ( &List::ExtraUtils::hashList($build->{'componentClasses'}) ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
select type (component => self%component{ucfirst($class->{'name'})})
type is (nodeComponent{ucfirst($class->{'name'})})
class default
  do i=1,size(component)
    count=component(i)%serializeCount(propertyType)
    if (count > 0) then
       call component(i)%deserializeValues(array(offset:),propertyType)
       offset=offset+count
    end if
  end do
end select
CODE
    }
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    descriptor  => $function,
	    name        => "deserializeValues"
	}
	);
}

sub Tree_Node_ODE_Serialize_Scales {
    # Generate a function to serialize scales of evolvable properties of a node into an array for the ODE solver.
    my $build = shift();
    # Define the function.
    my $function =
    {
	type        => "void",
	name        => "treeNodeSerializeScaleToArray",
	description => "Serialize scales of evolvable properties of a node into an array.",
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => "treeNode",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "double precision",
		 attributes => [ "dimension(:)", "intent(  out)" ],
		 variables  => [ "array" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "propertyType" ],
	     }
	    ]
    };    
    $function->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
!GCC$ attributes unused :: self
select case (propertyType)
case (propertyTypeAll     )
 array(1:nodeSerializationCount        )=     nodeScales(1:nodeSerializationCount)
case (propertyTypeActive  )
 array(1:nodeSerializationCountActive  )=pack(nodeScales(1:nodeSerializationCount),.not.nodeInactives(1:nodeSerializationCount))
case (propertyTypeInactive)
 array(1:nodeSerializationCountInactive)=pack(nodeScales(1:nodeSerializationCount),     nodeInactives(1:nodeSerializationCount))
end select
CODE
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    descriptor  => $function,
	    name        => "serializeScales"
	}
	);
}

sub Tree_Node_ODE_Serialize_Rates {
    # Generate a function to serialize rates of evolvable properties of a node into an array for the ODE solver.
    my $build = shift();
    # Define the function.
    my $function =
    {
	type        => "void",
	name        => "treeNodeSerializeRatesToArray",
	description => "Serialize rates of evolvable properties of a node into an array.",
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => "treeNode",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "double precision",
		 attributes => [ "dimension(:)", "intent(  out)" ],
		 variables  => [ "array" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "propertyType" ],
	     }
	    ]
    };    
    $function->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
!GCC$ attributes unused :: self
select case (propertyType)
case (propertyTypeAll     )
 array(1:nodeSerializationCount        )=nodeRates(1:nodeSerializationCount        )
case (propertyTypeActive  )
 array(1:nodeSerializationCountActive  )=nodeRates(1:nodeSerializationCountActive  )
case (propertyTypeInactive)
 array(1:nodeSerializationCountInactive)=nodeRates(1:nodeSerializationCountInactive)
end select
CODE
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    descriptor  => $function,
	    name        => "serializeRates"
	}
	);
}

sub Tree_Node_ODE_Serialize_Inactive {
    # Generate a function to serialize inactive status of evolvable properties of a node into an array for the ODE solver.
    my $build = shift();
    my $function =
    {
	type        => "void",
	name        => "treeNodeSerializeInactiveToArray",
	description => "Serialize inactive statuses of evolvable properties of a node into an array.",
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => "treeNode",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "logical",
		 attributes => [ "dimension(:)", "intent(  out)" ],
		 variables  => [ "array" ]
	     }
	    ]
    };    
    $function->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
!GCC$ attributes unused :: self
array(1:nodeSerializationCount)=nodeInactives(1:nodeSerializationCount)
CODE
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    descriptor  => $function,
	    name        => "serializeInactives"
	}
	);
}

sub Tree_Node_ODE_Name_From_Index {
    # Generate a function to return the name of a property given the index of that property in the ODE solver array.
    my $build = shift();
    my $function =
    {
	type        => "type(varying_string) => name",
	name        => "treeNodePropertyNameFromIndex",
	description => "Return the name of a property given its index within an array.",
	modules     =>
	    [
	     "ISO_Varying_String"
	    ],
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => "treeNode",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "index", "propertyType" ]
	     },
	     {
		 intrinsic  => "integer",
		 variables  => [ "count", "i" ]
	     }
	    ]
    };    
    $function->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
name='unknown'
count =index
CODE
    # Iterate over all component classes
    foreach $code::class ( &List::ExtraUtils::hashList($build->{'componentClasses'}) ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if (allocated(self%component{ucfirst($class->{'name'})})) then
  do i=1,size(self%component{ucfirst($class->{'name'})})
    name=self%component{ucfirst($class->{'name'})}(i)%nameFromIndex(count,propertyType)
    if (count <= 0) return
  end do
end if
CODE
    }
    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if (name == 'unknown') call Galacticus_Error_Report('property index out of range'//\{introspection:location\})
CODE
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    descriptor  => $function,
	    name        => "nameFromIndex"
	}
	);
}

sub Tree_Node_ODE_Offsets {
    # Generate a function to compute offsets into serialization arrays for components of a node.
    my $build = shift();
    my $function =
    {
	type        => "void",
	name        => "treeNodeSerializeOffsets",
	description => "Compute offsets into serialization arrays for {\\normalfont \\ttfamily treeNode} object.",
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => "treeNode",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "propertyType" ]
	     },
	     {
		 intrinsic  => "integer",
		 variables  => [ "i", "count", "countSubset" ]
	     }
	    ]
    };    
    $function->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
count      =0
countSubset=0
CODE
    # Iterate over all component classes
    foreach $code::class ( &List::ExtraUtils::hashList($build->{'componentClasses'}) ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if (allocated(self%component{ucfirst($class->{'name'})})) then
  do i=1,size(self%component{ucfirst($class->{'name'})})
    call self%component{ucfirst($class->{'name'})}(i)%serializationOffsets(count,countSubset,propertyType)
  end do
end if
CODE
    }
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if (.not.allocated(nodeScales)) then
   allocate  (nodeScales   (count))
   allocate  (nodeRates    (count))
   allocate  (nodeInactives(count))
else if (size(nodeScales) < count) then
   deallocate(nodeScales          )
   deallocate(nodeRates           )
   deallocate(nodeInactives       )
   allocate  (nodeScales   (count))
   allocate  (nodeRates    (count))
   allocate  (nodeInactives(count))
end if
nodeSerializationCount         =count
select case (propertyType)
case (propertyTypeInactive)
 nodeSerializationCountInactive=countSubset
case (propertyTypeActive  )
 nodeSerializationCountActive  =countSubset
end select
CODE
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    descriptor  => $function,
	    name        => "serializationOffsets"
	}
	);
}

1;
