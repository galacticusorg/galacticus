# Contains a Perl module which handles output of component classes.

package Galacticus::Build::Components::Classes::Output;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use NestedMap;
use Text::Template 'fill_in_string';
use List::ExtraUtils;
use Galacticus::Build::Components::Utils qw($fullyQualifiedNameLengthMax &isOutputIntrinsic);
use Galacticus::Build::Components::DataTypes;

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     classesOutput =>
     {
	 classIteratedFunctions =>
	     [
	      \&Class_Dump_ASCII  ,
	      \&Class_Output_Count,
	      \&Class_Output_Names,
	      \&Class_Post_Output ,
	      \&Class_Output
	     ]
     }
    );

sub Class_Dump_ASCII {
    # Generate a function to dump component classes.
    my $build    = shift();
    $code::class = shift();
    return
	unless ( grep {$code::class->{'name'} eq $_} @{$build->{'componentClassListActive'}} );
    my $function =
    {
	type        => "void",
	name        => "nodeComponent".ucfirst($code::class->{'name'})."DumpASCII",
	description => "Dump the content of a {\\normalfont \\ttfamily ".$code::class->{'name'}."} component.",
	modules     =>
	    [
	     "Display",
	     "ISO_Varying_String"
	    ],
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($code::class->{'name'}),
		 attributes => [ "intent(in   )" ],
		 variables  => [ "self" ]
	     }
	    ]
    };
    $code::padding = " " x ($fullyQualifiedNameLengthMax-length($code::class->{'name'}));
    $function->{'content'}  = fill_in_string(<<'CODE', PACKAGE => 'code');
!$GLC attributes unused :: self
call displayIndent('{$class->{'name'}}: {$padding}generic')
call displayUnindent('done')
CODE
    # Insert a type-binding for this function.
    push(
	@{$build->{'types'}->{"nodeComponent".ucfirst($code::class->{'name'})}->{'boundFunctions'}},
	{
	    type        => "procedure",
	    descriptor  => $function,
	    name        => "dumpASCII", 
	}
	);
}

sub Class_Output_Count {
    # Generate a function to return a count of the number of properties to be output from a generic node component.
    my $build    = shift();   
    $code::class = shift();
    return
	unless ( grep {$code::class->{'name'} eq $_} @{$build->{'componentClassListActive'}} );
    my $function =
    {
	type        => "void",
	name        => $code::class->{'name'}."OutputCount",
	description => "Increment the count of properties to output for a generic {\\normalfont \\ttfamily ".$code::class->{'name'}."} component.",
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
		 attributes => [ "intent(inout)" ],
		 variables  => [ "integerPropertyCount", "doublePropertyCount" ]
	     },
	     {
		 intrinsic  => "double precision",
		 attributes => [ "intent(in   )" ],
		 variables  => [ "time" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => ["intent(in   )" ],
		 variables  => [ "instance" ]
	     },
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($code::class->{'name'}),
		 attributes => [ "allocatable" ],
		 variables  => [ "selfDefault" ]
	     }
	    ]
    };
    $function->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
allocate(selfDefault,source=default{ucfirst($class->{'name'})}Component)
selfDefault%hostNode => self%hostNode
call selfDefault%outputCount(integerPropertyCount,doublePropertyCount,time,instance)
CODE
    # Insert a type-binding for this function into the node component class type.
    push(
	@{$build->{'types'}->{'nodeComponent'.ucfirst($code::class->{'name'})}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    descriptor  => $function,
	    name        => "outputCount"
	}
	);
}

sub Class_Output_Names {
    # Generate a function to return names of properties to be output from a generic node component.
    my $build    = shift();   
    $code::class = shift();
    return
	unless ( grep {$code::class->{'name'} eq $_} @{$build->{'componentClassListActive'}} );
    my $function =
    {
	type        => "void",
	name        => $code::class->{'name'}."OutputNames",
	description => "Establish the names of properties to output for a generic {\\normalfont \\ttfamily ".$code::class->{'name'}."} component.",
	modules     =>
	    [
	     "Merger_Tree_Outputter_Buffer_Types"
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
		 intrinsic  => "integer", 
		 attributes => [ "intent(inout)" ], 
		 variables  => [ "integerProperty" ]
	     },
	     {
		 intrinsic  => "type",
		 type       => "outputPropertyInteger",
		 attributes => [ "intent(inout)", "dimension(:)" ], 
		 variables  => [ "integerProperties" ]
	     },
	     {
		 intrinsic  => "integer", 
		 attributes => [ "intent(inout)" ], 
		 variables  => [ "doubleProperty" ]
	     },
	     {
		 intrinsic  => "type",
		 type       => "outputPropertyDouble",
		 attributes => [ "intent(inout)", "dimension(:)" ], 
		 variables  => [ "doubleProperties" ]
	     },
	     {
		 intrinsic  => "double precision",
		 attributes => [ "intent(in   )" ], 
		 variables  => [ "time" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => ["intent(in   )" ],
		 variables  => [ "instance" ]
	     },
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($code::class->{'name'}),
		 attributes => [ "allocatable" ],
		 variables  => [ "selfDefault" ]
	     }
	    ]
    };    
    $function->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
allocate(selfDefault,source=default{ucfirst($class->{'name'})}Component)
selfDefault%hostNode => self%hostNode
call selfDefault%outputNames(integerProperty,integerProperties,doubleProperty,doubleProperties,time,instance)
CODE
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{'nodeComponent'.ucfirst($code::class->{'name'})}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    descriptor  => $function,
	    name        => "outputNames"
	}
	);
}

sub Class_Output {
    # Generate a function to populate output buffers with data from a generic node component.
    my $build    = shift();   
    $code::class = shift();
    return
	unless ( grep {$code::class->{'name'} eq $_} @{$build->{'componentClassListActive'}} );
    my $function =
    {
	type        => "void",
	name        => $code::class->{'name'}."Output",
	description => "Populate output buffers with properties to output for a {\\normalfont \\ttfamily ".$code::class->{'name'}."} component.",
	content     => "",
	modules     =>
	    [
	     "Multi_Counters"                    ,
	     "Merger_Tree_Outputter_Buffer_Types"
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
		 intrinsic  => "integer", 
		 attributes => [ "intent(inout)" ], 
		 variables  => [ "integerProperty", "integerBufferCount" ]
	     },
	     {
		 intrinsic  => "type",
		 type       => "outputPropertyInteger",
		 attributes => [ "intent(inout)", "dimension(:)" ], 
		 variables  => [ "integerProperties" ]
	     },
	     {
		 intrinsic  => "integer", 
		 attributes => [ "intent(inout)" ], 
		 variables  => [ "doubleProperty", "doubleBufferCount" ]
	     },
	     {
		 intrinsic  => "type",
		 type       => "outputPropertyDouble",
		 attributes => [ "intent(inout)", "dimension(:)" ], 
		 variables  => [ "doubleProperties" ]
	     },
	     {
		 intrinsic  => "double precision",
		 attributes => [ "intent(in   )" ], 
		 variables  => [ "time" ]
	     },
	     {
		 intrinsic  => "type",
		 type       => "multiCounter",
		 attributes => [ "intent(in   )" ], 
		 variables  => [ "outputInstance" ]
	     },
	     {
		 intrinsic  => "integer",
		 attributes => ["intent(in   )" ],
		 variables  => [ "instance" ]
	     }
	    ]
    };
    # Find all derived types to be output, and any modules required.
    my %intrinsicTypeMap =
	(
	 double      => {label => "double" , intrinsic => "double precision"                     },
	 integer     => {label => "integer", intrinsic => "integer"         , type => "kind_int8"},
	 longInteger => {label => "integer", intrinsic => "integer"         , type => "kind_int8"}
	);
    my %argumentUsage =
	(
	 "self"     => 0,
	 "time"     => 0,
	 "instance" => 0,
	 "integer"  => 0,
	 "double"   => 0
	);
    my %tmpsAdded =
	(
	 "double"  => 0,
	 "integer" => 0,
	 "index"   => 0
	);
    my %outputDerivedTypes;
    my %modulesRequired;
    # Iterate over class member implementations.
    foreach my $member ( @{$code::class->{'members'}} ) {
	# Check if the instances argument is to be used.
	$argumentUsage{'instance'} = 1
	    if
	    (
	     exists($member->{'output'}               )           &&
	     exists($member->{'output'}->{'instances'})           &&
	            $member->{'output'}->{'instances'} eq "first"
	    );
	# Iterate over all properties belonging to this member.	
	foreach my $property ( &List::ExtraUtils::hashList($member->{'properties'}->{'property'},keyAs => 'name') ) {
	    next
		unless ( exists($property->{'output'}) );
	    $argumentUsage                   {'self'                                                     } = 1;
	    $argumentUsage                   {'time'                                                     } = 1
		unless ( &isOutputIntrinsic($property->{'data'}->{'type'}) );
	    $argumentUsage                   {$intrinsicTypeMap{$property->{'data'}->{'type'}}->{'label'}} = 1
		if     ( &isOutputIntrinsic($property->{'data'}->{'type'}) );
	    unless ( &isOutputIntrinsic($property->{'data'}->{'type'}) ) {
		$argumentUsage               {'integer'                                                  } = 1;
		$argumentUsage               {'double'                                                   } = 1;
	    }
	    $outputDerivedTypes              {                  $property->{'data'}->{'type'}}             = 1
		unless ( &isOutputIntrinsic($property->{'data'}->{'type'}) );
	    map {$modulesRequired{$_} = 1} split(/,/,$property->{'output'}->{'modules'})
		if ( exists($property->{'output'}->{'modules'}) );
	}
    }
    # Add output variables for all derived types to be output.
    push(
	@{$function->{'variables'}},
	map
	{
	    {
		intrinsic => "type",
		type      => $_,
		variables => [ "output".ucfirst($_) ]
	    }
	}
	&List::ExtraUtils::sortedKeys(\%outputDerivedTypes)
	);
    # Add all required modules.
    push(@{$function->{'modules'}},&List::ExtraUtils::sortedKeys(\%modulesRequired));
    # Determine unused arguments.
    @code::argumentsUnused = 
	nestedmap
        {
	 $argumentUsage{$NestedMap::stack[0]}
	 ?
	 ()
	 :
	 (
	  ($NestedMap::stack[0] eq "integer" || $NestedMap::stack[0] eq "double")
	  ?
	  nestedmap {$NestedMap::stack[1].$NestedMap::stack[0]} ( "Property", "BufferCount", "Properties" )
	  :
	  $NestedMap::stack[0]
	 )
	} &List::ExtraUtils::sortedKeys(\%argumentUsage);
    push(@code::argumentsUnused,"outputInstance");
    if ( scalar(@code::argumentsUnused) > 0 ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
!$GLC attributes unused :: {join(", ",@argumentsUnused)}
CODE
    }
    # Iterate over class member implementations.
    foreach $code::member ( @{$code::class->{'members'}} ) {
	# Skip members with no outputs.
	next
	    unless ( grep {exists($_->{'output'})} &List::ExtraUtils::hashList($code::member->{'properties'}->{'property'}) );
	# Build code for this member.
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if (default{ucfirst($class->{'name'})}Component%{$member->{'name'}}IsActive() {(exists($member->{'output'}->{'instances'}) && $member->{'output'}->{'instances'} eq "first") ? " .and. instance == 1" : ""}) then
CODE
	# Iterate over all properties belonging to this member which are to be output.
	foreach $code::property ( grep {exists($_->{'output'}) && ! $_->{'definedInParent'} } &List::ExtraUtils::hashList($code::member->{'properties'}->{'property'}, keyAs => 'name' ) ) {
	    # Determine output count.
	    if      ( $code::property->{'data'}->{'rank'} == 0 ) {
		$code::count = 1;
	    } elsif ( $code::property->{'data'}->{'rank'} == 1 ) {
		$code::count = (my @matches = $code::property->{'output'}->{'labels'} =~ m/^\[(.*)\]$/) ? ($1 =~ tr/,//)+1 : $code::property->{'output'}->{'count'};
	    }
	    # Determine which output buffer type to use.
	    $code::bufferType = $intrinsicTypeMap{$code::property->{'data'}->{'type'}}->{'label'};
	    # Increment the counters.
	    if ( &isOutputIntrinsic($code::property->{'data'}->{'type'}) ) {
		if ( $code::property->{'data'}->{'rank'} == 0 ) {
		    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
{$bufferType}Property={$bufferType}Property+1
{$bufferType}Properties({$bufferType}Property)%scalar({$bufferType}BufferCount)=self%{$property->{'name'}}()
CODE
		} else {
		    if ( $code::bufferType eq "integer" && ! $tmpsAdded{'integer'} ) {
			$tmpsAdded{'integer'} = 1;
			push(
			    @{$function->{'variables'}},
			    {
				intrinsic  => "integer",
				type       => "kind_int8",
				attributes => [ "allocatable", "dimension(:)" ],
				variables  => [ "integerOutputTmp" ]
			    }
			    );
		    }
		    if ( $code::bufferType eq "double" && ! $tmpsAdded{'double'} ) {
			$tmpsAdded{'double'} = 1;
			push(
			    @{$function->{'variables'}},
			    {
				intrinsic  => "double precision",
				attributes => [ "allocatable", "dimension(:)" ],
				variables  => [ "doubleOutputTmp" ]
			    }
			    );
		    }
		    if ( ! $tmpsAdded{'index'} ) {
			$tmpsAdded{'index'} = 1;
			push(
			    @{$function->{'variables'}},
			    {
				intrinsic => "integer",
				variables => [ "i" ]
			    }
			    );
		    }
		    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
{$bufferType}OutputTmp=reshape(self%{$property->{'name'}}(),[{$count}])
do i=1,{$count}
  {$bufferType}Properties({$bufferType}Property+i-1)%scalar({$bufferType}BufferCount)={$bufferType}OutputTmp(i)
end do
deallocate({$bufferType}OutputTmp)
{$bufferType}Property={$bufferType}Property+{$count}
CODE
		}
	    } else {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
output{ucfirst($property->{'data'}->{'type'})}=self%{$property->{'name'}}()
call output{ucfirst($property->{'data'}->{'type'})}%output(integerProperty,integerBufferCount,integerProperties,doubleProperty,doubleBufferCount,doubleProperties,time,outputInstance)
if (.not.same_type_as(self,{$class->{'name'}}Class)) call self%{$property->{'name'}}Set(output{ucfirst($property->{'data'}->{'type'})})
CODE
	    }
	}
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
end if
CODE
    }
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{'nodeComponent'.ucfirst($code::class->{'name'})}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    descriptor  => $function,
	    name        => "output"
	}
	);
}

sub Class_Post_Output {
    # Generate a function to perform post-output processing of a generic node component.
    my $build    = shift();   
    $code::class = shift();
    my $function =
    {
	type        => "void",
	name        => $code::class->{'name'}."PostOutput",
	description => "Perform post-output processing of a {\\normalfont \\ttfamily ".$code::class->{'name'}."} component.",
	content     => "!\$GLC attributes unused :: self, time\n",
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($code::class->{'name'}),
		 attributes => [ "intent(inout)" ],
		 variables  => [ "self" ]
	     },
	     {
		 intrinsic  => "double precision",
		 attributes => [ "intent(in   )" ], 
		 variables  => [ "time" ]
	     }
	    ]
    };
    # Insert a type-binding for this function into the treeNode type.
    push(
	@{$build->{'types'}->{'nodeComponent'.ucfirst($code::class->{'name'})}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    descriptor  => $function,
	    name        => "postOutput"
	}
	);
}

1;
