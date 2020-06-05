# Contains a Perl module which provides functions to map other functions over components of a node.

package Galacticus::Build::Components::TreeNodes::Map;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Text::Template 'fill_in_string';
use List::ExtraUtils;
use Data::Dumper;
use Galacticus::Build::Components::Utils;
use Galacticus::Build::Components::DataTypes;

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     treeNodeMap =>
     {
	 functions =>
	     [
	      \&Tree_Node_Map_Optimizations,
	      \&Tree_Node_Map_Void         ,
	      \&Tree_Node_Map_Double0      ,
	      \&Tree_Node_Map_Double1      ,
	      \&Tree_Node_Map_TensorR2D3
	     ]
     }
    );

sub Tree_Node_Map_Void {
    # Generate a function to map a void function over components of a node.
    my $build = shift();
    my $function =
    {
	type        => "void",
	name        => "treeNodeMapVoid",
	description => "Map a void function over components of the node.",
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => "treeNode",
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "procedure",
		 type       => "Node_Component_Null_Void0_InOut",
		 attributes => [ "pointer" ],
		 variables  => [ "mapFunction" ]
	     },
	     {
		 intrinsic  => "integer",
		 variables  => [ "i" ]
	     }
	    ]
    };    
    # Iterate over all component classes
    foreach $code::class ( &List::ExtraUtils::hashList($build->{'componentClasses'}) ) {
	next
	    unless ( grep {$code::class->{'name'} eq $_} @{$build->{'componentClassListActive'}} );
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
do i=1,size(self%component{ucfirst($class->{'name'})})
  call mapFunction(self%component{ucfirst($class->{'name'})}(i))
end do
CODE
    }
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    descriptor  => $function,
	    name        => "mapVoid"
	}
	);
}

sub Tree_Node_Map_Double0 {
    # Generate a function to map a rank-0, double function over components of a node.
    my $build = shift();
    my $function =
    {
	type        => "double precision",
	name        => "treeNodeMapDouble0",
	description => "Map a rank-0, double function over components of the node.",
	modules     =>
	    [
	     "Galacticus_Error"
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
		 intrinsic  => "procedure",
		 type       => "Node_Component_Null_Double0_InOut",
		 attributes => [ "pointer" ],
		 variables  => [ "mapFunction" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "reduction" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => [ "intent(in   )", "optional" ],
		 variables  => [ "optimizeFor" ]
	     },
	     {
		 intrinsic  => "double precision",
		 variables  => [ "componentValue" ]
	     },
	     {
		 intrinsic  => "integer",
		 variables  => [ "i" ]
	     }
	    ]
    };
    # Determine if we will build optimized versions of the map operator for specific bound functions.
    my $hasOptimizations = grep {exists($_->{'mappable'})} @{$build->{'types'}->{'nodeComponent'}->{'boundFunctions'}};
    # Scan through available node component methods and find ones which are mappable. Create optimized versions of this function
    # for them.
    $code::firstOptimization = 1;
    $code::reductionIdentity =
    {
	summation => "0.0d0",
	product   => "1.0d0"
    };
    $code::reductionOperator =
    {
	summation => "+",
	product   => "*"
    };
    foreach $code::boundFunction ( @{$build->{'types'}->{'nodeComponent'}->{'boundFunctions'}} ) {
	next
	    unless ( exists($code::boundFunction->{'mappable'}) );
	foreach $code::reduction ( split(/:/,$code::boundFunction->{'mappable'}) ) {
	    # Insert test for optimized case.
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
{$firstOptimization ? "" : "else "}if (present(optimizeFor).and.optimizeFor == optimizeFor{ucfirst($boundFunction->{'name'}).ucfirst($reduction)}) then
    if (reduction /= reduction{ucfirst($reduction)}) call Galacticus_Error_Report('reduction mismatch'//\{introspection:location\})
    treeNodeMapDouble0={$reductionIdentity->{$reduction}}
CODE
            # Iterate over classes.
            foreach $code::class ( &List::ExtraUtils::hashList($build->{'componentClasses'}) ) {
		next
		    unless ( grep {$code::class->{'name'} eq $_} @{$build->{'componentClassListActive'}} );
		# If any members of this class override the base class for this method then evaluate them.
		if ( 
		    grep {$_->{'method'} eq $code::boundFunction->{'name'}}
		    map {@{$_->{'bindings'}->{'binding'}}}
		    @{$code::class->{'members'}} 
		    ) {
		    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
do i=1,size(self%component{ucfirst($class->{'name'})})
  treeNodeMapDouble0=treeNodeMapDouble0{$reductionOperator->{$reduction}}mapFunction(self%component{ucfirst($class->{'name'})}(i))
end do
CODE
		}
            }
            # First optimization is completed.
            $code::firstOptimization = 0;
	}
    }    
    # Generate the generic, unoptimized function.
    if ( $hasOptimizations ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
else
CODE
    }
    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
select case (reduction)
case (reductionSummation)
  treeNodeMapDouble0=0.0d0
case (reductionProduct  )
  treeNodeMapDouble0=1.0d0
case default
  treeNodeMapDouble0=1.0d0
  call Galacticus_Error_Report('unknown reduction'//\{introspection:location\})
end select
CODE
    foreach $code::class ( &List::ExtraUtils::hashList($build->{'componentClasses'}) ) {
	next
	    unless ( grep {$code::class->{'name'} eq $_} @{$build->{'componentClassListActive'}} );
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
do i=1,size(self%component{ucfirst($class->{'name'})})
  componentValue=mapFunction(self%component{ucfirst($class->{'name'})}(i))
  select case (reduction)
  case (reductionSummation)
    treeNodeMapDouble0=treeNodeMapDouble0+componentValue
  case (reductionProduct  )
    treeNodeMapDouble0=treeNodeMapDouble0*componentValue
  end select
end do
CODE
    }
    if ( $hasOptimizations ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
end if
CODE
    }
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    descriptor  => $function,
	    name        => "mapDouble0"
	}
	);
}

sub Tree_Node_Map_Double1 {
    # Generate a function to map a rank-1, double function over components of a node.
    my $build = shift();
    my $function =
    {
	type        => "double precision, dimension(resultSize) => double1Result",
	name        => "treeNodeMapDouble1",
	description => "Map a rank-0, double function over components of the node.",
	modules     =>
	    [
	     "Galacticus_Error"
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
		 intrinsic  => "procedure",
		 type       => "Node_Component_Null_Double1_InOut",
		 attributes => [ "pointer" ],
		 variables  => [ "mapFunction" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "resultSize", "reduction" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => [ "intent(in   )", "optional" ],
		 variables  => [ "optimizeFor" ]
	     },
	     {
		 intrinsic  => "double precision",
		 attributes => [ "dimension(resultSize)" ],
		 variables  => [ "componentValue" ]
	     },
	     {
		 intrinsic  => "integer",
		 variables  => [ "i" ]
	     }
	    ]
    };
    # Determine if we will build optimized versions of the map operator for specific bound functions.
    my $hasOptimizations = grep {exists($_->{'mappable'})} @{$build->{'types'}->{'nodeComponent'}->{'boundFunctions'}};
    # Scan through available node component methods and find ones which are mappable. Create optimized versions of this function
    # for them.
    $code::firstOptimization = 1;
    $code::reductionIdentity =
    {
	summation => "0.0d0",
	product   => "1.0d0"
    };
    $code::reductionOperator =
    {
	summation => "+",
	product   => "*"
    };
    foreach $code::boundFunction ( @{$build->{'types'}->{'nodeComponent'}->{'boundFunctions'}} ) {
	next
	    unless ( exists($code::boundFunction->{'mappable'}) );
	foreach $code::reduction ( split(/:/,$code::boundFunction->{'mappable'}) ) {
	    # Insert test for optimized case.
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
{$firstOptimization ? "" : "else "}if (present(optimizeFor).and.optimizeFor == optimizeFor{ucfirst($boundFunction->{'name'}).ucfirst($reduction)}) then
    if (reduction /= reduction{ucfirst($reduction)}) call Galacticus_Error_Report('reduction mismatch'//\{introspection:location\})
    double1Result={$reductionIdentity->{$reduction}}
CODE
            # Iterate over classes.
            foreach $code::class ( &List::ExtraUtils::hashList($build->{'componentClasses'}) ) {
		next
		    unless ( grep {$code::class->{'name'} eq $_} @{$build->{'componentClassListActive'}} );
		# If any members of this class override the base class for this method then evaluate them.
		if ( 
		    grep {$_->{'method'} eq $code::boundFunction->{'name'}}
		    map {@{$_->{'bindings'}->{'binding'}}}
		    @{$code::class->{'members'}} 
		    ) {	
		    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
do i=1,size(self%component{ucfirst($class->{'name'})})
  double1Result=double1Result{$reductionOperator->{$reduction}}mapFunction(self%component{ucfirst($class->{'name'})}(i),resultSize)
end do
CODE
		}
            }
            # First optimization is completed.
            $code::firstOptimization = 0;
	}
    }    
    # Generate the generic, unoptimized function.
    if ( $hasOptimizations ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
else
CODE
    }
    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
select case (reduction)
case (reductionSummation)
  double1Result=0.0d0
case (reductionProduct  )
  double1Result=1.0d0
case default
  double1Result=1.0d0
  call Galacticus_Error_Report('unknown reduction'//\{introspection:location\})
end select
CODE
    foreach $code::class ( &List::ExtraUtils::hashList($build->{'componentClasses'}) ) {
	next
	    unless ( grep {$code::class->{'name'} eq $_} @{$build->{'componentClassListActive'}} );
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
do i=1,size(self%component{ucfirst($class->{'name'})})
  componentValue=mapFunction(self%component{ucfirst($class->{'name'})}(i),resultSize)
  select case (reduction)
  case (reductionSummation)
    double1Result=double1Result+componentValue
  case (reductionProduct  )
    double1Result=double1Result*componentValue
  end select
end do
CODE
    }
    if ( $hasOptimizations ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
end if
CODE
    }
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    descriptor  => $function,
	    name        => "mapDouble1"
	}
	);
}

sub Tree_Node_Map_Optimizations {
    # Generate an enumeration of optimizations for map functions.
    my $build = shift();   
    # Determine if we will build optimized versions of the map operator for specific bound functions.
    my $hasOptimizations = grep {exists($_->{'mappable'})} @{$build->{'types'}->{'nodeComponent'}->{'boundFunctions'}};
    # Create an enumeration for optimization types if necessary.
    if ( $hasOptimizations ) {
	my $optimizationEnumerationValue = -1;
	foreach my $boundFunction ( @{$build->{'types'}->{'nodeComponent'}->{'boundFunctions'}} ) {
	    next
		unless ( exists($boundFunction->{'mappable'}) );
	    foreach my $reduction ( split(/:/,$boundFunction->{'mappable'}) ) {
		# Insert enumeration entry for this optimization.
		++$optimizationEnumerationValue;
		push
		    (
		     @{$build->{'variables'}},
		     {
			 intrinsic  => "integer",
			 attributes => [ "public", "parameter" ],
			 variables  => [ "optimizeFor".ucfirst($boundFunction->{'name'}).ucfirst($reduction)."=".$optimizationEnumerationValue ]
		     }
		    );
	    }
	}
    }
}

sub Tree_Node_Map_TensorR2D3 {
    # Generate a function to map a rank-2, dimensio-3 tensor function over components of a node.
    my $build = shift();
    my $function =
    {
	type        => "type(tensorRank2Dimension3Symmetric) => tensorResult",
	name        => "treeNodeMapTensorR2D3",
	description => "Map a rank-2, dimension-3 tensor function over components of the node.",
	modules     =>
	    [
	     "Galacticus_Error"
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
		 intrinsic  => "procedure",
		 type       => "Node_Component_Null_TensorR2D3_InOut",
		 attributes => [ "pointer" ],
		 variables  => [ "mapFunction" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "reduction" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => [ "intent(in   )", "optional" ],
		 variables  => [ "optimizeFor" ]
	     },
	     {
		 intrinsic  => "type(tensorRank2Dimension3Symmetric)",
		 variables  => [ "componentValue" ]
	     },
	     {
		 intrinsic  => "integer",
		 variables  => [ "i" ]
	     }
	    ]
    };
    # Determine if we will build optimized versions of the map operator for specific bound functions.
    my $hasOptimizations = grep {exists($_->{'mappable'})} @{$build->{'types'}->{'nodeComponent'}->{'boundFunctions'}};
    # Scan through available node component methods and find ones which are mappable. Create optimized versions of this function
    # for them.
    $code::firstOptimization = 1;
    $code::reductionIdentity =
    {
	summation => "tensorNullR2D3Sym"
    };
    $code::reductionOperator =
    {
	summation => "+"
    };
    foreach $code::boundFunction ( @{$build->{'types'}->{'nodeComponent'}->{'boundFunctions'}} ) {
	next
	    unless ( exists($code::boundFunction->{'mappable'}) );
	foreach $code::reduction ( split(/:/,$code::boundFunction->{'mappable'}) ) {
	    # Insert test for optimized case.
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
{$firstOptimization ? "" : "else "}if (present(optimizeFor).and.optimizeFor == optimizeFor{ucfirst($boundFunction->{'name'}).ucfirst($reduction)}) then
    if (reduction /= reduction{ucfirst($reduction)}) call Galacticus_Error_Report('reduction mismatch'//\{introspection:location\})
    tensorResult={$reductionIdentity->{$reduction}}
CODE
            # Iterate over classes.
            foreach $code::class ( &List::ExtraUtils::hashList($build->{'componentClasses'}) ) {
		next
		    unless ( grep {$code::class->{'name'} eq $_} @{$build->{'componentClassListActive'}} );
		# If any members of this class override the base class for this method then evaluate them.
		if ( 
		    grep {$_->{'method'} eq $code::boundFunction->{'name'}}
		    map {@{$_->{'bindings'}->{'binding'}}}
		    @{$code::class->{'members'}} 
		    ) {	
		    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
do i=1,size(self%component{ucfirst($class->{'name'})})
  tensorResult=tensorResult{$reductionOperator->{$reduction}}mapFunction(self%component{ucfirst($class->{'name'})}(i))
end do
CODE
		}
            }
            # First optimization is completed.
            $code::firstOptimization = 0;
	}
    }    
    # Generate the generic, unoptimized function.
    if ( $hasOptimizations ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
else
CODE
    }
    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
select case (reduction)
case (reductionSummation)
  tensorResult=tensorNullR2D3Sym
case default
  tensorResult=tensorNullR2D3Sym
  call Galacticus_Error_Report('unknown reduction'//\{introspection:location\})
end select
CODE
    foreach $code::class ( &List::ExtraUtils::hashList($build->{'componentClasses'}) ) {
	next
	    unless ( grep {$code::class->{'name'} eq $_} @{$build->{'componentClassListActive'}} );
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
do i=1,size(self%component{ucfirst($class->{'name'})})
  componentValue=mapFunction(self%component{ucfirst($class->{'name'})}(i))
  select case (reduction)
  case (reductionSummation)
    tensorResult=tensorResult+componentValue
  end select
end do
CODE
    }
    if ( $hasOptimizations ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
end if
CODE
    }
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{'treeNode'}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    descriptor  => $function,
	    name        => "mapTensorR2D3"
	}
	);
}

1;
