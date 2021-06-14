# Contains a Perl module which handles output of component implementations.

package Galacticus::Build::Components::Implementations::Output;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use NestedMap;
use Text::Template 'fill_in_string';
use List::ExtraUtils;
use List::Uniq ':all';
use Data::Dumper;
use Galacticus::Build::Components::Utils qw($fullyQualifiedNameLengthMax &isOutputIntrinsic %outputTypeMap);
use Galacticus::Build::Components::DataTypes;

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     implementationsOutput =>
     {
	 implementationIteratedFunctions =>
	     [
	      \&Implementation_Output_Count,
	      \&Implementation_Output_Names,
	      \&Implementation_Post_Output
	     ]
     }
    );

sub Implementation_Output_Count {
    # Generate a function to return a count of the number of properties to be output from a node component implementation.
    my $build     = shift();   
    $code::class  = shift();
    $code::member = shift();
    # Build the function.
    my $function =
    {
	type        => "void",
	name        => $code::class->{'name'}.ucfirst($code::member->{'name'})."OutputCount",
	description => "Increment the count of properties to output for a {\\normalfont \\ttfamily ".$code::member->{'name'}."} implementation of the {\\normalfont \\ttfamily ".$code::class->{'name'}."} component.",
	content     => "",
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($code::class->{'name'}).ucfirst($code::member->{'name'}),
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
	     }
	    ]
    };
    # Perform common output tasks (add modules, variables required, and identify unused arguments).
    &Implementation_Output_Common_Tasks($build,$code::class,$code::member,$function,[ "PropertyCount" ]);
    # Get the count of the parent type if necessary.
    if ( exists($code::member->{'extends'}) ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
call self%nodeComponent{ucfirst($member->{'extends'}->{'class'}).ucfirst($member->{'extends'}->{'name'})}%outputCount(integerPropertyCount,doublePropertyCount,time,instance)
CODE
    }

    # If only the first instance is to be output, check instance number and return if not first.
    if ( exists($code::member->{'output'}->{'instances'}) && $code::member->{'output'}->{'instances'} eq "first" ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if (instance > 1) return
CODE
    }
    
    # Initialize counts of fixed-size properties.
    %code::fixedSizeCount = ( integer => 0, double => 0 );
    # Iterate over properties, generating code to increment output counts appropriately.
    foreach $code::property ( &List::ExtraUtils::hashList($code::member->{'properties'}->{'property'}) ) {
	# Skip properties with no output or which are defined in their parent implementation.
	next
	    unless ( exists($code::property->{'output'}) && ! $code::property->{'definedInParent'} );
	# Detect whether intrinsic or not.
	if ( &isOutputIntrinsic($code::property->{'data'}->{'type'}) ) {
	    $code::outputType = $outputTypeMap{$code::property->{'data'}->{'type'}};
	    if ( $code::property->{'data'}->{'rank'} == 0 ) {
		++$code::fixedSizeCount{$code::outputType};
	    } elsif ( $code::property->{'data'}->{'rank'} == 1 ) {
		$code::count = (exists($code::property->{'output'}->{'labels'}) && $code::property->{'output'}->{'labels'} =~ m/^\[(.*)\]$/) ? ($1 =~ tr/,//)+1 : $code::property->{'output'}->{'count'};
		if ( $code::count =~ m/^\d/ ) {
		    $code::fixedSizeCount{$code::outputType} += $code::count;
		} else {
		    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
{$outputType}PropertyCount={$outputType}PropertyCount+{$count}
CODE
		}  
	    }
	} else {
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
output{ucfirst($property->{'data'}->{'type'})}=self%{$property->{'name'}}()
call output{ucfirst($property->{'data'}->{'type'})}%outputCount(integerPropertyCount,doublePropertyCount,time)
CODE
	}
    }	
    # Insert fixed-size-counts.
    foreach $code::type ( keys(%code::fixedSizeCount) ) {
	if ( $code::fixedSizeCount{$code::type} > 0 ) {
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
{$type}PropertyCount={$type}PropertyCount+{$fixedSizeCount{$type}}
CODE
	}
    }
    # Insert a type-binding for this function into the node component class type.
    push(
	@{$build->{'types'}->{'nodeComponent'.ucfirst($code::class->{'name'}).ucfirst($code::member->{'name'})}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    descriptor  => $function,
	    name        => "outputCount"
	}
	);
}

sub Implementation_Output_Names {
    # Generate a function to return the names of properties to be output from a node component implementation.
    my $build     = shift();   
    $code::class  = shift();
    $code::member = shift();
    # Build the function.
    my $function =
    {
	type        => "void",
	name        => $code::class->{'name'}.ucfirst($code::member->{'name'})."OutputNames",
	description => "Return the names of properties to output for a {\\normalfont \\ttfamily ".$code::member->{'name'}."} implementation of the {\\normalfont \\ttfamily ".$code::class->{'name'}."} component.",
	content     => "",
	modules     =>
	    [
	     "Merger_Tree_Outputter_Buffer_Types"
	    ],
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($code::class->{'name'}).ucfirst($code::member->{'name'}),
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
	     }
	    ]
    };
    # Perform common output tasks (add modules, variables required, and identify unused arguments).
    &Implementation_Output_Common_Tasks($build,$code::class,$code::member,$function,[ "Property", "Properties" ]);
    # Get the count of the parent type if necessary.
    if ( exists($code::member->{'extends'}) ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
call self%nodeComponent{ucfirst($member->{'extends'}->{'class'}).ucfirst($member->{'extends'}->{'name'})}%outputNames(integerProperty,integerProperties,doubleProperty,doubleProperties,time,instance)
CODE
    }
    # If only the first instance is to be output, check instance number and return if not first.
    if ( exists($code::member->{'output'}->{'instances'}) && $code::member->{'output'}->{'instances'} eq "first" ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
if (instance > 1) return
CODE
    }
    # Iterate over properties, generating code to set names appropriately.
    foreach $code::property ( &List::ExtraUtils::hashList($code::member->{'properties'}->{'property'}) ) {
	# Skip properties with no output or which are defined in their parent implementation.
	next
	    unless ( exists($code::property->{'output'}) && ! $code::property->{'definedInParent'} );
	# Detect whether intrinsic or not.
	if ( &isOutputIntrinsic($code::property->{'data'}->{'type'}) ) {
	    $code::outputType = $outputTypeMap{$code::property->{'data'}->{'type'}};
	    if ( $code::property->{'data'}->{'rank'} == 0 ) {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
{$outputType}Property                                   = {$outputType}Property+1
{$outputType}Properties({$outputType}Property)%name     ='{$class->{'name'}.ucfirst($property->{'name'})}'
{$outputType}Properties({$outputType}Property)%comment  ='{$property->{'output'}->{'comment'  }}'
{$outputType}Properties({$outputType}Property)%unitsInSI= {$property->{'output'}->{'unitsInSI'}}
CODE
	    } elsif ( $code::property->{'data'}->{'rank'} == 1 ) {
		if ( $code::property->{'output'}->{'labels'} =~ m/^\[(.*)\]$/ ) {
		    (my $labels = $1) =~ s/\s//g;
	            foreach $code::label ( split(",",$labels) ) {
			$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
{$outputType}Property                                   ={$outputType}Property+1
{$outputType}Properties({$outputType}Property)%name     ='{$class->{'name'}.ucfirst($property->{'name'}).$label}'
{$outputType}Properties({$outputType}Property)%comment  ='{$property->{'output'}->{'comment'  }} [{$label}]'
{$outputType}Properties({$outputType}Property)%unitsInSI= {$property->{'output'}->{'unitsInSI'}}
CODE
		    }
		} elsif ( exists($code::property->{'output'}->{'count'}) ) {
		    (my $label = $code::property->{'output'}->{'labels'}) =~ s/\{i\}/i/g;
		    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
do i=1,{$property->{'output'}->{'count'}}
   {$outputType}Property                                   ={$outputType}Property+1
   {$outputType}Properties({$outputType}Property)%name     ='{$class->{'name'}.ucfirst($propertyName)}'//{$label}
   {$outputType}Properties({$outputType}Property)%comment  ='{$property->{'output'}->{'comment'  }} [' //{$label}//']'
   {$outputType}Properties({$outputType}Property)%unitsInSI= {$property->{'output'}->{'unitsInSI'}}
end do
CODE
                }
	    }
	} else {
	    $code::unitsInSI = exists($code::property->{'output'}->{'unitsInSI'}) ? $code::property->{'output'}->{'unitsInSI'} : "0.0d0";
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
output{ucfirst($property->{'data'}->{'type'})}=self%{$property->{'name'}}()			   
call output{ucfirst($property->{'data'}->{'type'})}%outputNames(integerProperty,integerProperties,doubleProperty,doubleProperties,time,'{$class->{'name'}.ucfirst($property->{'name'})}','{$property->{'output'}->{'comment'}}',{$unitsInSI})
CODE
	}
    }
    # Insert a type-binding for this function into the node component class type.
    push(
	@{$build->{'types'}->{'nodeComponent'.ucfirst($code::class->{'name'}).ucfirst($code::member->{'name'})}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    descriptor  => $function,
	    name        => "outputNames"
	}
	);
}

sub Implementation_Post_Output {
    # Generate a function to perform post-output processing of a node component implementation.
    my $build     = shift();   
    $code::class  = shift();
    $code::member = shift();
    # Build the function.
    my $function =
    {
	type        => "void",
	name        => $code::class->{'name'}.ucfirst($code::member->{'name'})."PostOutput",
	description => "Perform post-output processing for a {\\normalfont \\ttfamily ".$code::member->{'name'}."} implementation of the {\\normalfont \\ttfamily ".$code::class->{'name'}."} component.",
	content     => "",
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($code::class->{'name'}).ucfirst($code::member->{'name'}),
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
    # Perform post-output processing of any parent type.
    if ( exists($code::member->{'extends'}) ) {
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
call self%nodeComponent{ucfirst($member->{'extends'}->{'class'}).ucfirst($member->{'extends'}->{'name'})}%postOutput(time)
CODE
    }
    # Iterate over properties, performing post-processing.
    foreach $code::property ( &List::ExtraUtils::hashList($code::member->{'properties'}->{'property'}, keyAs => 'name') ) {
	# Skip properties with no output or which are defined in their parent implementation.
	next
	    unless ( exists($code::property->{'output'}) && ! $code::property->{'definedInParent'} );
	# Detect whether intrinsic or not.
	unless ( &isOutputIntrinsic($code::property->{'data'}->{'type'}) ) {
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
call self%{$property->{'name'}}Data%postOutput(time)
CODE
	}
    }
    # Insert a type-binding for this function into the node component class type if it has any content.
    push(
	@{$build->{'types'}->{'nodeComponent'.ucfirst($code::class->{'name'}).ucfirst($code::member->{'name'})}->{'boundFunctions'}},
	{
	    type        => "procedure", 
	    descriptor  => $function,
	    name        => "postOutput"
	}
	)
	unless ( $function->{'content'} eq "" );
}

sub Implementation_Output_Common_Tasks {
    # Perform some common tasks required for all implementation output functions.
    my $build        =   shift() ;   
    my $class        =   shift() ;
    my $member       =   shift() ;
    my $function     =   shift() ;
    my @typeSuffixes = @{shift()};
    # Determine used arguments.
    my %argumentsUsed;
    $argumentsUsed{$_} = exists($member->{'extends'}) ? 1 : 0
	foreach ( 'self', 'time', 'instance', 'integer', 'double' );
    $argumentsUsed{'instance'} = 1
	if ( exists($member->{'output'}->{'instances'}) && $member->{'output'}->{'instances'} eq "first" );
    # Examine all properties. Determine modules and variables required.
    my @modulesRequired;
    my @typesRequired;
    foreach my $property ( &List::ExtraUtils::hashList($member->{'properties'}->{'property'}) ) {
	# Skip properties with no output or which are defined in their parent implementation.
	next
	    unless ( exists($property->{'output'}) && ! $code::property->{'definedInParent'} );
	# Extract modules required.
	if ( exists($property->{'output'}->{'modules'}) ) {	    
	    (my $moduleList = $property->{'output'}->{'modules'}) =~ s/\s//g;
	    push(@modulesRequired,split(/,/,$moduleList));
	}
	# For non-output-intrinsic types, record that a variable will be required.
	push(@typesRequired,$property->{'data'}->{'type'})
	    unless ( &isOutputIntrinsic($property->{'data'}->{'type'}) );
	# Record argument usage.
	unless (  &isOutputIntrinsic($property->{'data'}->{'type'}) ) {
	    $argumentsUsed{$_} = 1
		foreach ( 'self', 'time' );
	}
	if ( &isOutputIntrinsic($property->{'data'}->{'type'}) ) {
	    $argumentsUsed{$outputTypeMap{$property->{'data'}->{'type'}}} = 1;
	} else {
	    $argumentsUsed{'integer'} = 1;
	    $argumentsUsed{'double' } = 1;
	}
    }
    # Store all required modules in the function descriptor.
    push(@{$function->{'modules'}},uniq(sort(@modulesRequired)));
    # Append all required variables in the function descriptor.
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
	uniq(sort(@typesRequired))
	); 
    # Add record of unused arguments.
    $function->{'content'} .= "!\$GLC attributes unused :: ".join(",",nestedmap {($NestedMap::stack[0] eq "integer" || $NestedMap::stack[0] eq "double") ? nestedmap {$NestedMap::stack[1].$NestedMap::stack[0]} @typeSuffixes : $NestedMap::stack[0]} grep {! $argumentsUsed{$_}} keys(%argumentsUsed))."\n"
	if ( grep {! $argumentsUsed{$_}} keys(%argumentsUsed) );
}

1;
