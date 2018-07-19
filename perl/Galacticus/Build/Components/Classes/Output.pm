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
	      \&Class_Dump_ASCII,
	      \&Class_Output_Count,
	      \&Class_Output_Names,
	      \&Class_Output
	     ]
     }
    );

sub Class_Dump_ASCII {
    # Generate a function to dump component classes.
    my $build    = shift();
    $code::class = shift();
    my $function =
    {
	type        => "void",
	name        => "nodeComponent".ucfirst($code::class->{'name'})."DumpASCII",
	description => "Dump the content of a {\\normalfont \\ttfamily ".$code::class->{'name'}."} component.",
	modules     =>
	    [
	     "Galacticus_Display",
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
!GCC$ attributes unused :: self
call Galacticus_Display_Indent('{$class->{'name'}}: {$padding}generic')
call Galacticus_Display_Unindent('done')
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
    my $function =
    {
	type        => "void",
	name        => $code::class->{'name'}."OutputNames",
	description => "Establish the names of properties to output for a generic {\\normalfont \\ttfamily ".$code::class->{'name'}."} component.",
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
		 intrinsic  => "character",
		 type       => "len=*",
		 attributes => [ "intent(inout)", "dimension(:)" ], 
		 variables  => [ "integerPropertyNames", "integerPropertyComments" ]
	     },
	     {
		 intrinsic  => "double precision",
		 attributes => [ "intent(inout)", "dimension(:)" ],
		 variables  => [ "integerPropertyUnitsSI" ]
	     },
	     {
		 intrinsic  => "integer", 
		 attributes => [ "intent(inout)" ], 
		 variables  => [ "doubleProperty" ]
	     },
	     {
		 intrinsic  => "character",
		 type       => "len=*",
		 attributes => [ "intent(inout)", "dimension(:)" ], 
		 variables  => [ "doublePropertyNames", "doublePropertyComments" ]
	     },
	     {
		 intrinsic  => "double precision",
		 attributes => [ "intent(inout)", "dimension(:)" ],
		 variables  => [ "doublePropertyUnitsSI" ]
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
call selfDefault%outputNames(integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI,doubleProperty,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time,instance)
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
    my $function =
    {
	type        => "void",
	name        => $code::class->{'name'}."Output",
	description => "Populate output buffers with properties to output for a {\\normalfont \\ttfamily ".$code::class->{'name'}."} component.",
	content     => "",
	modules     =>
	    [
	     "Multi_Counters"
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
		 intrinsic  => "integer",
		 type       => "kind=kind_int8",
		 attributes => [ "intent(inout)", "dimension(:,:)" ],
		 variables  => [ "integerBuffer" ]
	     },
	     {
		 intrinsic  => "integer", 
		 attributes => [ "intent(inout)" ], 
		 variables  => [ "doubleProperty", "doubleBufferCount" ]
	     },
	     {
		 intrinsic  => "double precision",
		 attributes => [ "intent(inout)", "dimension(:,:)" ],
		 variables  => [ "doubleBuffer" ]
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
    # Find all derived types to be output, whether we have any rank-1, conditional, intrinsic outputs, and any modules required.
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
    my %outputDerivedTypes;
    my %rank1ConditionalIntrinsicOutputs;
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
		unless ( &isOutputIntrinsic($property->{'data'}->{'type'})                                                                                       );
	    $argumentUsage                   {$intrinsicTypeMap{$property->{'data'}->{'type'}}->{'label'}} = 1
		if     ( &isOutputIntrinsic($property->{'data'}->{'type'})                                                                                       );
	    unless ( &isOutputIntrinsic($property->{'data'}->{'type'}) ) {
		$argumentUsage               {'integer'                                                  } = 1;
		$argumentUsage               {'double'                                                   } = 1;
	    }
	    $outputDerivedTypes              {                  $property->{'data'}->{'type'}}             = 1
		unless ( &isOutputIntrinsic($property->{'data'}->{'type'})                                                                                       );
	    $rank1ConditionalIntrinsicOutputs{                  $property->{'data'}->{'type'}}             = 1
		if ( &isOutputIntrinsic($property->{'data'}->{'type'}) && $property->{'data'}->{'rank'} == 1 && exists($property->{'output'}->{'condition'}) );
	    map {$modulesRequired{$_} = 1} split(/,/,$property->{'output'}->{'modules'})
		if ( exists($property->{'output'}->{'modules'}) );
	}
    }
    # Add a counter variable for rank-1, conditional, intrinsic outputs.
    push(
	@{$function->{'variables'}},
	{
	    intrinsic  => "integer",
	    variables  => [ "i" ]
	}
	)
	if ( %rank1ConditionalIntrinsicOutputs );
    # Add output arrays for rank-1, conditional, intrinsic outputs.
    push(
	@{$function->{'variables'}},
	map
	{
	    my %descriptor = %{$intrinsicTypeMap{$_}};
	    $descriptor{'attributes'} = [ "allocatable", "dimension(:)" ];
	    $descriptor{'variables' } = [ "outputRank1".ucfirst($intrinsicTypeMap{$_}->{'label'}) ];
	    \%descriptor
	}
        &List::ExtraUtils::sortedKeys(\%rank1ConditionalIntrinsicOutputs)
	);
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
	  nestedmap {$NestedMap::stack[1].$NestedMap::stack[0]} ( "Property", "BufferCount", "Buffer" )
	  :
	  $NestedMap::stack[0]
	 )
	} &List::ExtraUtils::sortedKeys(\%argumentUsage);
    push(@code::argumentsUnused,"outputInstance");
    if ( scalar(@code::argumentsUnused) > 0 ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
!GCC$ attributes unused :: {join(", ",@argumentsUnused)}
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
	    # Construct the output condition string.
	    $code::condition = exists($code::property->{'output'}->{'condition'}) ? $code::property->{'output'}->{'condition'} : "";
	    $code::condition =~ s/\[\[([^\]]+)\]\]/$1/g;
	    # Increment the counters.
	    if ( &isOutputIntrinsic($code::property->{'data'}->{'type'}) ) {
		if ( $code::property->{'data'}->{'rank'} == 0 ) {
		    if ( exists($code::property->{'output'}->{'condition'}) ) {
			$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if ({$condition}) then
CODE
		    }
		    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
{$bufferType}Property={$bufferType}Property+1
{$bufferType}Buffer({$bufferType}BufferCount,{$bufferType}Property)=self%{$property->{'name'}}()
CODE
		    if ( exists($code::property->{'output'}->{'condition'}) ) {
			$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
end if
CODE
		    }
		} else {
		    if ( exists($code::property->{'output'}->{'condition'}) ) {
			$code::condition =~ s/\{i\}/i/g;
			$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
outputRank1{ucfirst({$bufferType})}=self%{$property->{'name'}}()
do i=1,{$property->{'output'}->{'count'}}
  if ({$condition}) then
    {$bufferType}Property={$bufferType}Property+1
    {$bufferType}Buffer({$bufferType}BufferCount,{$bufferType}Property)=outputRank1{ucfirst({$bufferType})}(i)
  end if
end do
deallocate(outputRank1{ucfirst({$bufferType})})
CODE
		    } else {
			$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
{$bufferType}Buffer({$bufferType}BufferCount,{$bufferType}Property+1:{$bufferType}Property+{$count})=reshape(self%{$property->{'name'}}(),[{$count}])
{$bufferType}Property={$bufferType}Property+{$count}
CODE
		    }
		}
	    } else {
		if ( exists($code::property->{'output'}->{'condition'}) ) {
		    my $condition = $code::property->{'output'}->{'condition'};
		    $condition =~ s/\[\[([^\]]+)\]\]/$1/g;
		    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if ({$condition}) then
CODE
		}
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
output{ucfirst($property->{'data'}->{'type'})}=self%{$property->{'name'}}()
call output{ucfirst($property->{'data'}->{'type'})}%output(integerProperty,integerBufferCount,integerBuffer,doubleProperty,doubleBufferCount,doubleBuffer,time,outputInstance)
if (.not.same_type_as(self,{$class->{'name'}}Class)) call self%{$property->{'name'}}Set(output{ucfirst($property->{'data'}->{'type'})})
CODE
		if ( exists($code::property->{'output'}->{'condition'}) ) {
		    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
end if
CODE
		}
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

1;
