# Contains a Perl module which handles deferred methods of component implementations.

package Galacticus::Build::Components::Implementations::Deferred;
use strict;
use warnings;
use utf8;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
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
	# Create a pointer at the implementation level.
	my $componentFunctionName = $class->{'name'}.ucfirst($member->{'name'}).ucfirst($binding->{'method'});
	push(
	    @{$build->{'variables'}},
	    {
		intrinsic  => "procedure",
		type       => $class->{'name'}.ucfirst($binding->{'method'})."Interface",
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
	# Create the function.
	my $function =
	{
	    type        => "void",
	    name        => $code::memberFunctionName."DeferredFunctionSet",
	    description => "Set the function to be used for the {\\normalfont \\ttfamily ".$binding->{'method'}."} method of the {\\normalfont \\ttfamily ".$member->{'name'}."} implementation of the {\\normalfont \\ttfamily ".$class->{'name'}."} component class.",
	    variables   =>
		[
		 {
		     intrinsic  => "procedure",
		     type       => $class->{'name'}.ucfirst($binding->{'method'})."Interface",
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
	$code::memberFunctionName = $code::class->{'name'}.ucfirst($code::member->{'name'}).ucfirst($code::binding->{'method'});
	# Create the function.
	my $function =
	{
	    type        => $code::binding->{'interface'}->{'type'} eq "void" ? "void" : $intrinsicTypes{$code::binding->{'interface'}->{'type'}},
	    name        => $code::memberFunctionName,
	    description => "Call the deferred function for the {\\normalfont \\ttfamily ".$code::binding->{'method'}."} method of the {\\normalfont \\ttfamily ".$code::class->{'name'}."} component class has been set.",
	    modules     =>
		[
		 "Galacticus_Error"
		]
	};
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
	     map {&Fortran::Utils::Unformat_Variables($_)} ( @{$code::binding->{'interface'}->{'argument'}} )
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
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
      {$memberFunctionName}={$memberFunctionName}Deferred({join(",",@arguments)})
CODE
	}
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
   else
CODE
	$code::parentType = undef();
	if ( exists($code::member->{'extends'}) ) {
	    $code::parentType = "nodeComponent".ucfirst($code::member->{'extends'}->{'class'}).ucfirst($code::member->{'extends'}->{'name'});
	} elsif ( $code::binding->{'bindsTo'} eq "componentClass" ) {
	    $code::parentType = "nodeComponent".ucfirst($code::class->{'name'});
	}
	if ( defined($code::parentType) ) {
	    if ( $code::binding->{'interface'}->{'type'} eq "void" ) {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
      call self%{$parentType}%{$binding->{'method'}}({join(",",grep {$_ ne "self"} @arguments)})
CODE
	    } else {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
      {$memberFunctionName}=self%{$parentType}%{$binding->{'method'}}({join(",",grep {$_ ne "self"} @arguments)})
CODE
	    }
	} else {
	    if ( &isIntrinsic($code::binding->{'interface'}->{'type'}) && $code::binding->{'interface'}->{'type'} ne "void" ) {
		$code::nullValue = $intrinsicNulls{$code::binding->{'interface'}->{'type'}};
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
      {$memberFunctionName}={$nullValue}
CODE
	    }
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
      call Galacticus_Error_Report('{$memberFunctionName}','deferred function has not been assigned')
CODE
	}
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
   end if
class default
CODE
	if ( &isIntrinsic($code::binding->{'interface'}->{'type'}) && $code::binding->{'interface'}->{'type'} ne "void" ) {
	    $code::nullValue = $intrinsicNulls{$code::binding->{'interface'}->{'type'}};
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
      {$memberFunctionName}={$nullValue}
CODE
	}
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
   call Galacticus_Error_Report('{$memberFunctionName}','incorrect class - this should not happen')
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
