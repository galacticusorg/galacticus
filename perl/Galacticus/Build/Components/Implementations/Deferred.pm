# Contains a Perl module which handles deferred methods of component implementations.

package Galacticus::Build::Components::Implementations::Deferred;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Data::Dumper;
use Text::Template 'fill_in_string';
use List::ExtraUtils;
use Galacticus::Build::Components::Utils qw(&isIntrinsic %intrinsicNulls %intrinsicTypes);
use Galacticus::Build::Components::NullFunctions qw(createNullFunction);

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     implementationsDeferred => 
     {
	 implementationIteratedFunctions =>
	     [
	      \&Implementation_Deferred_Binding_Pointers     ,
	      \&Implementation_Deferred_Binding_Attachers    ,
	      \&Implementation_Deferred_Binding_Attach_Status,
	      \&Implementation_Deferred_Binding_Wrappers
	     ]
     }
    );

sub Implementation_Deferred_Binding_Pointers {
    # Generate procedure pointers for deferred methods of component implementations.
    my $build  = shift();
    my $class  = shift();
    my $member = shift();
    # Iterate over bindings.
    foreach my $binding ( @{$member->{'bindings'}->{'binding'}} ) {
	# Skip non-deferred bindings.
	next
	    unless ( $binding->{'isDeferred'} );
	# Determine the name of the procedure template.
	my $baseMember = $member;
	while ( exists($baseMember->{'extends'}) && grep {$_->{'method'} eq $binding->{'method'}} @{$baseMember->{'extends'}->{'implementation'}->{'bindings'}->{'binding'}} ) {
	    $baseMember = $baseMember->{'extends'};
	}
	# Create a pointer at the implementation level.
	my $componentFunctionName = $class->{'name'}.ucfirst($member->{'name'}).ucfirst($binding->{'method'});
	push(
	    @{$build->{'variables'}},
	    {
		intrinsic  => "procedure",
		type       => $class->{'name'}.ucfirst($baseMember->{'name'}).ucfirst($binding->{'method'}),
		attributes => [ "pointer" ],
		variables  => [ $componentFunctionName."Deferred" ]
	    },
	    {
		intrinsic  => "logical",
		variables  => [ $componentFunctionName."IsSetValue=.false." ]
	    }
	    );
    }
}

sub Implementation_Deferred_Binding_Attachers {
    # Generate functions to attach functions to deferred methods of component implementations.
    my $build  = shift();
    my $class  = shift();
    my $member = shift();
    # Iterate over deferred bindings which bind at the class level.
    foreach my $binding ( grep {$_->{'isDeferred'}} @{$member->{'bindings'}->{'binding'}} ) {
	# Create the name of the function.
	$code::memberFunctionName = $class->{'name'}.ucfirst($member->{'name'}).ucfirst($binding->{'method'});
	# Determine the name of the procedure template.
	my $baseMember = $member;
	while ( exists($baseMember->{'extends'}) && grep {$_->{'method'} eq $binding->{'method'}} @{$baseMember->{'extends'}->{'implementation'}->{'bindings'}->{'binding'}} ) {
	    $baseMember = $baseMember->{'extends'};
	}
	# Create the function.
	my $function =
	{
	    type        => "void",
	    name        => $code::memberFunctionName."DfrrdFnctnSet",
	    description => "Set the function to be used for the {\\normalfont \\ttfamily ".$binding->{'method'}."} method of the {\\normalfont \\ttfamily ".$member->{'name'}."} implementation of the {\\normalfont \\ttfamily ".$class->{'name'}."} component class.",
	    variables   =>
		[
		 {
		     intrinsic  => "procedure",
		     type       => $class->{'name'}.ucfirst($baseMember->{'name'}).ucfirst($binding->{'method'}),
		     isArgument => 1,
		     variables  => [ "deferredFunction" ]
		 }
		]
	};
	$function->{'content'}  = fill_in_string(<<'CODE', PACKAGE => 'code');
{$memberFunctionName}Deferred   => deferredFunction
{$memberFunctionName}IsSetValue =  .true.
CODE
	# Insert a type-binding for this function into the node component class type.
	push(
	    @{$build->{'types'}->{'nodeComponent'.ucfirst($class->{'name'}).ucfirst($member->{'name'})}->{'boundFunctions'}},
	    {
		type        => "procedure", 
		pass        => "nopass",
		descriptor  => $function,
		name        => $binding->{'method'}."Function"
	    }
	);
    }
}

sub Implementation_Deferred_Binding_Attach_Status {
    # Generate functions to return the attachment status of deferred methods of component implementations.
    my $build  = shift();
    my $class  = shift();
    my $member = shift();
    # Iterate over deferred bindings which bind at the class level.
    foreach my $binding ( grep {$_->{'isDeferred'}} @{$member->{'bindings'}->{'binding'}} ) {
	# Create the name of the function.
	$code::memberFunctionName = $class->{'name'}.ucfirst($member->{'name'}).ucfirst($binding->{'method'});
	# Create the function.
	my $function =
	{
	    type        => "logical",
	    name        => $code::memberFunctionName."DfrrdFnctnIsSet",
	    description => "Return true if the deferred function for the {\\normalfont \\ttfamily ".$binding->{'method'}."} method of the {\\normalfont \\ttfamily ".$member->{'name'}."} implementation of the {\\normalfont \\ttfamily ".$class->{'name'}."} component class has been set."
	};
	$function->{'content'}  = fill_in_string(<<'CODE', PACKAGE => 'code');
{$memberFunctionName}DfrrdFnctnIsSet={$memberFunctionName}IsSetValue
CODE
	# Insert a type-binding for this function into the node component class type.
	push(
	    @{$build->{'types'}->{'nodeComponent'.ucfirst($class->{'name'}).ucfirst($member->{'name'})}->{'boundFunctions'}},
	    {
		type        => "procedure", 
		pass        => "nopass",
		descriptor  => $function,
		name        => $binding->{'method'}."FunctionIsSet"
	    }
	);
    }
}

sub Implementation_Deferred_Binding_Wrappers {
    # Generate wrapper functions to call deferred methods of component implementations.
    my $build     = shift();
    $code::class  = shift();
    $code::member = shift();
    # Iterate over deferred bindings which bind at the class level.
    foreach $code::binding ( grep {$_->{'isDeferred'}} @{$code::member->{'bindings'}->{'binding'}} ) {
	# Create the name of the function.
	$code::memberFunctionName  = $code::class->{'name'}.ucfirst($code::member->{'name'}).ucfirst($code::binding->{'method'});
	$code::returnName          = $code::memberFunctionName;
	my $specificType           = $code::binding->{'interface'}->{'type'} ne "void" && ! exists($intrinsicTypes{$code::binding->{'interface'}->{'type'}});
	my $isPointer              = $specificType && $code::binding->{'interface'}->{'type'} =~ m/,\s*pointer/;
	$code::returnName         .="_"
	    if ( $specificType );
	# Create the function.
	my $function =
	{
	    type        => $code::binding->{'interface'}->{'type'} eq "void" ? "void" : ($specificType ? $code::binding->{'interface'}->{'type'}." => ".$code::returnName : $intrinsicTypes{$code::binding->{'interface'}->{'type'}}),
	    name        => $code::memberFunctionName,
	    description => "Call the deferred function for the {\\normalfont \\ttfamily ".$code::binding->{'method'}."} method of the {\\normalfont \\ttfamily ".$code::class->{'name'}."} component class if it has been set.",
	    modules     =>
		[
		 "Error"
		]
	};
	push(@{$function->{'modules'}},&List::ExtraUtils::as_array($code::binding->{'interface'}->{'module'}))
	    if ( exists($code::binding->{'interface'}->{'module'}) );
	# Handle any rank/shape.
	if      ( exists($code::binding->{'interface'}->{'rank'}) && $code::binding->{'interface'}->{'rank'} > 0 ) {
	    die('can not specify both "rank" and "shape"')
		if ( exists($code::binding->{'interface'}->{'shape'}) );
	    $code::returnName .="_";
	    $function->{'type'} .= ", allocatable, dimension(".join(",",map {":"} 1..$code::binding->{'interface'}->{'rank'}).") => ".$code::returnName;
	} elsif ( exists($code::binding->{'interface'}->{'shape'}) ) {
	    die('can not specify both "rank" and "shape"')
		if ( exists($code::binding->{'interface'}->{'rank' }) );
	    $code::returnName .="_";
	    $function->{'type'} .= ", dimension(".$code::binding->{'interface'}->{'shape'}.") => ".$code::returnName;
	}
	# Add a "self" argument if required.
	push
	    (
	     @{$function->{'variables'}},
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($code::class->{'name'}).ucfirst($code::member->{'name'}),
		 attributes => [ "intent(".$code::binding->{'interface'}->{'self'}->{'intent'}.")" ],
		 variables  => [ "self" ]
	     }
	    ) if ( $code::binding->{'interface'}->{'self'}->{'pass'} );
	# Add any other arguments.
	push
	    (
	     @{$function->{'variables'}},
	     map {&Fortran::Utils::Unformat_Variables($_)} ( &List::ExtraUtils::as_array($code::binding->{'interface'}->{'argument'}) )
	    );
	# Extract a list of all argument names.
	@code::arguments = &Galacticus::Build::Components::Utils::argumentList(@{$function->{'variables'}});
	# Generate function code.	    
	$function->{'content'}  = fill_in_string(<<'CODE', PACKAGE => 'code');
select type (self)
class is (nodeComponent{ucfirst($class->{'name'}).ucfirst($member->{'name'})})
   if (self%{$binding->{'method'}}FunctionIsSet()) then
CODE
	if ( $code::binding->{'interface'}->{'type'} eq "void" ) {
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
      call {$memberFunctionName}Deferred({join(",",@arguments)})
CODE
	} else {
	    $code::assigner = $isPointer ? " => " : "=";
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
      {$returnName}{$assigner}{$memberFunctionName}Deferred({join(",",@arguments)})
CODE
	}
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
   else
CODE
	$code::parentType = undef();
	$code::parentType = "nodeComponent".ucfirst($code::member->{'extends'}->{'class'}).ucfirst($code::member->{'extends'}->{'name'})
	    if ( exists($code::member->{'extends'}) && grep {$_->{'method'} eq $code::binding->{'method'}} @{$code::member->{'extends'}->{'implementation'}->{'bindings'}->{'binding'}} );
	if ( defined($code::parentType) ) {
	    if ( $code::binding->{'interface'}->{'type'} eq "void" ) {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
      call self%{$parentType}%{$binding->{'method'}}({join(",",grep {$_ ne "self"} @arguments)})
CODE
	    } else {
		$code::assigner = $isPointer ? " => " : "=";
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
      {$returnName}{$assigner}self%{$parentType}%{$binding->{'method'}}({join(",",grep {$_ ne "self"} @arguments)})
CODE
	    }
	} else {
	    if ( &isIntrinsic($code::binding->{'interface'}->{'type'}) && $code::binding->{'interface'}->{'type'} ne "void" ) {
		$code::nullValue = $intrinsicNulls{$code::binding->{'interface'}->{'type'}};
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
      {$returnName}={$nullValue}
CODE
	    }
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
      call Error_Report('deferred function has not been assigned'//\{introspection:location\})
CODE
	}
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
   end if
class default
CODE
	if ( &isIntrinsic($code::binding->{'interface'}->{'type'}) && $code::binding->{'interface'}->{'type'} ne "void" ) {
	    $code::nullValue = $intrinsicNulls{$code::binding->{'interface'}->{'type'}};
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
      {$returnName}={$nullValue}
CODE
	}
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
   call Error_Report('incorrect class - this should not happen'//\{introspection:location\})
end select
CODE
	# Insert a type-binding for this function into the node component class type.
	push(
	    @{$build->{'types'}->{'nodeComponent'.ucfirst($code::class->{'name'}).ucfirst($code::member->{'name'})}->{'boundFunctions'}},
	    {
		type        => "procedure", 
		descriptor  => $function,
		name        => $code::binding->{'method'}
	    }
	    );
    }
}

1;
