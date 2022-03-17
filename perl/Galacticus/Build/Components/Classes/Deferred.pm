# Contains a Perl module which handles deferred methods of component classes.

package Galacticus::Build::Components::Classes::Deferred;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Data::Dumper;
use Text::Template 'fill_in_string';
use List::ExtraUtils;
use Fortran::Utils;
use Galacticus::Build::Components::Utils qw(&isIntrinsic %intrinsicNulls %intrinsicTypes);
use Galacticus::Build::Components::DataTypes;
use Galacticus::Build::Components::NullFunctions qw(createNullFunction);

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     classesDeferred => 
     {
	 implementationIteratedFunctions =>
	     [
	      \&Class_Deferred_Binding_Pointers     ,
	      \&Class_Deferred_Binding_Interfaces   ,
	      \&Class_Deferred_Binding_Attachers    ,
	      \&Class_Deferred_Binding_Attach_Status,
	      \&Class_Deferred_Binding_Wrappers
	     ]
     }
    );

# Record of pointers and interfaces made.
my %classPointers;
my %classInterfaces;
my $classFunctions;

sub Class_Deferred_Binding_Pointers {
    # Generate procedure pointers for deferred methods of component implementations.
    my $build  = shift();
    my $class  = shift();
    my $member = shift();
    # Iterate over bindings.
    foreach my $binding ( @{$member->{'bindings'}->{'binding'}} ) {
	# Skip non-deferred bindings.
	next
	    unless ( $binding->{'isDeferred'} );
	# Create a pointer for the component class level if needed.
	my $classFunctionName = $class->{'name'}.ucfirst($binding->{'method'});
	$classPointers{$classFunctionName} = 1;
	next
	    unless ( $binding->{'bindsTo'} eq 'componentClass' && ! exists($classPointers{$classFunctionName}) );
	# Create the pointer.
	push(
	    @{$build->{'variables'}},
	    {
		intrinsic  => "procedure",
		type       => $class->{'name'}.ucfirst($binding->{'method'})."Interface",
		attributes => [ "pointer" ],
		variables  => [ $classFunctionName."Deferred" ]
	    },
	    {
		intrinsic  => "logical",
		variables  => [ $classFunctionName."IsSetValue=.false." ]
	    }
	    );
    }
}

sub Class_Deferred_Binding_Interfaces {
    # Generate interfaces for deferred methods of component implementations.
    my $build  = shift();
    my $class  = shift();
    my $member = shift();
    # Iterate over bindings.
    foreach my $binding ( @{$member->{'bindings'}->{'binding'}} ) {
	# Skip non-deferred bindings.
	next
	    unless ( $binding->{'isDeferred'} );
	# Create an abstract interface for the deferred function.
	my $interfaceName = $class->{'name'}.ucfirst($binding->{'method'});
	next
	    if ( exists($classInterfaces{$interfaceName}) );
	$classInterfaces{$interfaceName} = 1;
	# Find the highest level in the node component class hierarchy at which this method should be bound.
	my $hierarchyBindPoint;
	if ( $binding->{'bindsTo'} eq "componentClass" ) {
	    $hierarchyBindPoint = $class->{'name'};
	} else {
	    my $parentImplementation = $member;
	    while ( exists($parentImplementation->{'extends'}) ) {
		$parentImplementation = $parentImplementation->{'extends'}->{'implementation'};
	    }
	    $hierarchyBindPoint = $parentImplementation->{'class'}.ucfirst($parentImplementation->{'name'});
	}
	# Determine the type of the interface.
	(my $type) = $binding->{'interface'}->{'type'} eq "void" ? ("void") : &Galacticus::Build::Components::DataTypes::dataObjectPrimitiveName($binding->{'interface'});
	# Construct the interface.
	my $interface =
	{
	    name      => $interfaceName."Interface",
	    intrinsic => $type,
	    comment   => "Interface for deferred function for {\\normalfont \\ttfamily ".$binding->{'method'}."} method of the {\\normalfont \\ttfamily ".$class->{'name'}."} class.",
	    data      => [ map {&Fortran::Utils::Unformat_Variables($_)} &List::ExtraUtils::as_array($binding->{'interface'}->{'argument'}) ]
	};
	# Add a self argument if required.
       	unshift
	    (
	     @{$interface->{'data'}},
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".$hierarchyBindPoint,
		 attributes => [ "intent(".$binding->{'interface'}->{'self'}->{'intent'}.")" ],
		 variables  => [ "self" ]
	     }
	    ) if ( $binding->{'interface'}->{'self'}->{'pass'} );
    	$build->{'interfaces'}->{$interfaceName} = $interface;
    }
}

sub Class_Deferred_Binding_Attachers {
    # Generate functions to attach functions to deferred methods of component implementations.
    my $build  = shift();
    my $class  = shift();
    my $member = shift();
    # Iterate over deferred bindings which bind at the class level.
    foreach my $binding ( grep {$_->{'isDeferred'} && $_->{'bindsTo'} eq 'componentClass'} @{$member->{'bindings'}->{'binding'}} ) {
	# Create the name of the function.
	$code::classFunctionName = $class->{'name'}.ucfirst($binding->{'method'});
	# Create the function if necesasary.
	next
	    if ( exists($classFunctions->{$code::classFunctionName}->{'attacher'}) );
	$classFunctions->{$code::classFunctionName}->{'attacher'} = 1;
	# Create the function.
	my $function =
	{
	    type        => "void",
	    name        => $code::classFunctionName."DeferredFunctionSet",
	    description => "Set the function to be used for the {\\normalfont \\ttfamily ".$binding->{'method'}."} method of the {\\normalfont \\ttfamily ".$class->{'name'}."} component class.",
	    variables   =>
		[
		 {
		     intrinsic  => "procedure",
		     type       => $code::classFunctionName."Interface",
		     isArgument => 1,
		     variables  => [ "deferredFunction" ]
		 }
		]
	};
	$function->{'content'}  = fill_in_string(<<'CODE', PACKAGE => 'code');
{$classFunctionName}Deferred   => deferredFunction
{$classFunctionName}IsSetValue =  .true.
CODE
	# Insert a type-binding for this function into the node component class type.
	push(
	    @{$build->{'types'}->{'nodeComponent'.ucfirst($class->{'name'})}->{'boundFunctions'}},
	    {
		type        => "procedure", 
		pass        => "nopass",
		descriptor  => $function,
		name        => $binding->{'method'}."Function"
	    }
	);
    }
}

sub Class_Deferred_Binding_Attach_Status {
    # Generate functions to return the attachment status of deferred methods of component implementations.
    my $build  = shift();
    my $class  = shift();
    my $member = shift();
    # Iterate over deferred bindings which bind at the class level.
    foreach my $binding ( grep {$_->{'isDeferred'} && $_->{'bindsTo'} eq 'componentClass'} @{$member->{'bindings'}->{'binding'}} ) {
	# Create the name of the function.
	$code::classFunctionName = $class->{'name'}.ucfirst($binding->{'method'});
	# Create the function if necesasary.
	next 
	    if ( exists($classFunctions->{$code::classFunctionName}->{'attachStatus'}) );
	$classFunctions->{$code::classFunctionName}->{'attachStatus'} = 1;
	# Create the function.
	my $function =
	{
	    type        => "logical",
	    name        => $code::classFunctionName."DfrrdFnctnIsSet",
	    description => "Return true if the deferred function for the {\\normalfont \\ttfamily ".$binding->{'method'}."} method of the {\\normalfont \\ttfamily ".$class->{'name'}."} component class has been set."
	};
	$function->{'content'}  = fill_in_string(<<'CODE', PACKAGE => 'code');
{$classFunctionName}DfrrdFnctnIsSet={$classFunctionName}IsSetValue
CODE
	# Insert a type-binding for this function into the node component class type.
	push(
	    @{$build->{'types'}->{'nodeComponent'.ucfirst($class->{'name'})}->{'boundFunctions'}},
	    {
		type        => "procedure", 
		pass        => "nopass",
		descriptor  => $function,
		name        => $binding->{'method'}."FunctionIsSet"
	    }
	    );
    }
}

sub Class_Deferred_Binding_Wrappers {
    # Generate wrapper functions to call deferred methods of component implementations.
    my $build  = shift();
    my $class  = shift();
    my $member = shift();
    # Iterate over deferred bindings which bind at the class level.
    foreach $code::binding ( grep {$_->{'isDeferred'} && $_->{'bindsTo'} eq 'componentClass'} @{$member->{'bindings'}->{'binding'}} ) {
	# Create the name of the function.
	$code::classFunctionName = $class->{'name'}.ucfirst($code::binding->{'method'});
	# Create the function if necesasary.
	next
	    if ( exists($classFunctions->{$code::classFunctionName}->{'wrapper'}) );
	$classFunctions->{$code::classFunctionName}->{'wrapper'} = 1;
	# Create the function.
	my $function =
	{
	    type        => $code::binding->{'interface'}->{'type'} eq "void" ? "void" : $intrinsicTypes{$code::binding->{'interface'}->{'type'}},
	    name        => $code::classFunctionName,
	    description => "Call the deferred function for the {\\normalfont \\ttfamily ".$code::binding->{'method'}."} method of the {\\normalfont \\ttfamily ".$class->{'name'}."} component class has been set.",
	    modules     =>
		[
		 "Error"
		]
	};
	# Add a "self" argument if required.
	push
	    (
	     @{$function->{'variables'}},
	     {
		 intrinsic  => "class",
		 type       => "nodeComponent".ucfirst($class->{'name'}),
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
if (self%{$binding->{'method'}}FunctionIsSet()) then
CODE
	if ( $code::binding->{'interface'}->{'type'} eq "void" ) {
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
   call {$classFunctionName}Deferred({join(",",@arguments)})
CODE
	} else {
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
   {$classFunctionName}={$classFunctionName}Deferred({join(",",@arguments)})
CODE
	}
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
else
CODE
	if ( &isIntrinsic($code::binding->{'interface'}->{'type'}) && $code::binding->{'interface'}->{'type'} ne "void" ) {
	    $code::nullValue = $intrinsicNulls{$code::binding->{'interface'}->{'type'}};
	    $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
   {$classFunctionName}={$nullValue}
CODE
	}
	$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
   call Error_Report('deferred function has not been assigned'//\{introspection:location\})
end if
CODE
	# Insert a type-binding for this function into the node component class type.
	push(
	    @{$build->{'types'}->{'nodeComponent'.ucfirst($class->{'name'})}->{'boundFunctions'}},
	    {
		type        => "procedure", 
		descriptor  => $function,
		name        => $code::binding->{'method'}
	    }
	    );
    }
}

1;
