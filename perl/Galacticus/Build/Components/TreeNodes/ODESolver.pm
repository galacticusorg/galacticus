# Contains a Perl module which provides various ODE solver-related functions for tree nodes.

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
require List::ExtraUtils;
require Galacticus::Build::Components::Utils;
require Galacticus::Build::Components::DataTypes;

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     treeNodeODESolver =>
     {
	 functions =>
	     [
	      \&Tree_Node_ODE_Serialize_Count      ,
	      \&Tree_Node_ODE_Serialize_Values     ,
	      \&Tree_Node_ODE_Deserialize_Values   ,
	      \&Tree_Node_ODE_Serialize_RatesScales,
	      \&Tree_Node_ODE_Name_From_Index      ,
	      \&Tree_Node_ODE_Offsets
	     ]
     }
    );

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
		 variables  => [ "i" ],
	     }
	    ]
    };    
    $function->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
treeNodeSerializeCount=0
CODE
    # Iterate over all component classes
    foreach $code::class ( &ExtraUtils::hashList($build->{'componentClasses'}) ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if (allocated(self%component{ucfirst($class->{'name'})})) then
  do i=1,size(self%component{ucfirst($class->{'name'})})
    treeNodeSerializeCount=treeNodeSerializeCount+self%component{ucfirst($class->{'name'})}(i)%serializeCount()
  end do
end if
CODE
    }
    # Add the function to the functions list.
    push(
	@{$build->{'functions'}},
	$function
	);
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    name        => "serializeCount", 
	    function    => "treeNodeSerializeCount", 
	    description => "Return a count of the number of evolvable properties of the serialized object.",
	    returnType  => "\\intzero", 
	    arguments   => ""
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
		 intrinsic  => "integer",
		 variables  => [ "count", "offset", "i" ],
	     },
	     {
		 intrinsic  => "double precision",
		 attributes => [ "dimension(:)", "intent(  out)" ],
		 variables  => [ "array" ]
	     }
	    ]
    };    
    $function->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
offset=1
CODE
    # Iterate over all component classes
    foreach $code::class ( &ExtraUtils::hashList($build->{'componentClasses'}) ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if (allocated(self%component{ucfirst($class->{'name'})})) then
  do i=1,size(self%component{ucfirst($class->{'name'})})
    count=self%component{ucfirst($class->{'name'})}(i)%serializeCount()
    if (count > 0) call self%component{ucfirst($class->{'name'})}(i)%serializeValues(array(offset:))
    offset=offset+count
  end do
end if
CODE
    }
    # Add the function to the functions list.
    push(
	@{$build->{'functions'}},
	$function
	);
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    name        => "serializeValues", 
	    function    => "treeNodeSerializeValuesToArray", 
	    description => "Serialize values to {\\normalfont \\ttfamily array}.",
	    returnType  => "\\void", 
	    arguments   => "\\doubleone\\ array\\argout"
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
		 intrinsic  => "integer",
		 variables  => [ "count", "offset", "i" ],
	     },
	     {
		 intrinsic  => "double precision",
		 attributes => [ "dimension(:)", "intent(in   )" ],
		 variables  => [ "array" ]
	     }
	    ]
    };    
    $function->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
offset=1
CODE
    # Iterate over all component classes
    foreach $code::class ( &ExtraUtils::hashList($build->{'componentClasses'}) ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if (allocated(self%component{ucfirst($class->{'name'})})) then
  do i=1,size(self%component{ucfirst($class->{'name'})})
    count=self%component{ucfirst($class->{'name'})}(i)%serializeCount()
    if (count > 0) call self%component{ucfirst($class->{'name'})}(i)%deserializeValues(array(offset:))
    offset=offset+count
  end do
end if
CODE
    }
    # Add the function to the functions list.
    push(
	@{$build->{'functions'}},
	$function
	);
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    name        => "deserializeValues", 
	    function    => "treeNodeDeserializeValuesFromArray", 
	    description => "Deserialize values from {\\normalfont \\ttfamily array}.",
	    returnType  => "\\void", 
	    arguments   => "\\doubleone\\ array\\argin"
	}
	);
}

sub Tree_Node_ODE_Serialize_RatesScales {
    # Generate a function to serialize rates and scales of evolvable properties of a node into an array for the ODE solver.
    my $build = shift();
    # Iterate over rates and scales.
    foreach $code::content ( "scale", "rate" ) {
	my $function =
	{
	    type        => "void",
	    name        => "treeNodeSerialize".ucfirst($code::content)."ToArray",
	    description => "Serialize ".$code::content."s of evolvable properties of a node into an array.",
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
		 }
		]
	};    
	$function->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
!GCC$ attributes unused :: self
array(1:nodeSerializationCount)=node{ucfirst($content)}s(1:nodeSerializationCount)
CODE
	# Add the function to the functions list.
	push(
	    @{$build->{'functions'}},
	    $function
	    );
	# Insert a type-binding for this function into the treeNode type.
	push(
	    @{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	    {
		type        => "procedure", 
		name        => "serialize".ucfirst($code::content)."s", 
		function    => "treeNodeSerialize".ucfirst($code::content)."ToArray",
		description => "Serialize ".$code::content."s to {\\normalfont \\ttfamily array}.",
		returnType  => "\\void", 
		arguments   => "\\doubleone\\ array\\argout"
	    }
	    );
    }
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
		 variables  => [ "index" ]
	     },
	     {
		 intrinsic  => "integer",
		 variables  => [ "count", "i" ]
	     }
	    ]
    };    
    $function->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
name='unknown'
count=index
CODE
    # Iterate over all component classes
    foreach $code::class ( &ExtraUtils::hashList($build->{'componentClasses'}) ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if (allocated(self%component{ucfirst($class->{'name'})})) then
  do i=1,size(self%component{ucfirst($class->{'name'})})
    call self%component{ucfirst($class->{'name'})}(i)%nameFromIndex(count,name)
    if (count <= 0) return
  end do
end if
CODE
    }
    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if (name == 'unknown') call Galacticus_Error_Report('treeNodePropertyNameFromIndex','property index out of range')
CODE
    # Add the function to the functions list.
    push(
	@{$build->{'functions'}},
	$function
	);
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    name        => "nameFromIndex", 
	    function    => "treeNodePropertyNameFromIndex", 
	    description => "Return the name of a property given its index in a node.",
	    returnType  => "\\textcolor{red}{\\textless varying\\_string\\textgreater}", 
	    arguments   => "\\intzero\\ index\\argin"
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
		 variables  => [ "i", "count" ]
	     }
	    ]
    };    
    $function->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
count=0
CODE
    # Iterate over all component classes
    foreach $code::class ( &ExtraUtils::hashList($build->{'componentClasses'}) ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if (allocated(self%component{ucfirst($class->{'name'})})) then
  do i=1,size(self%component{ucfirst($class->{'name'})})
    call self%component{ucfirst($class->{'name'})}(i)%serializationOffsets(count)
  end do
end if
CODE
    }
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if (.not.allocated(nodeScales)) then
   allocate  (nodeScales        (count))
   allocate  (nodeRates         (count))
   allocate  (nodeRatesIncrement(count))
else if (size(nodeScales) < count) then
   deallocate(nodeScales               )
   deallocate(nodeRates                )
   deallocate(nodeRatesIncrement       )
   allocate  (nodeScales        (count))
   allocate  (nodeRates         (count))
   allocate  (nodeRatesIncrement(count))
end if
nodeSerializationCount=count
CODE
    # Add the function to the functions list.
    push(
	@{$build->{'functions'}},
	$function
	);
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    name        => "serializationOffsets", 
	    function    => "treeNodeSerializeOffsets", 
	    description => "Compute offsets into serialization arrays for all properties.",
	    returnType  => "\\void", 
	    arguments   => ""
	}
	);
}

1;
