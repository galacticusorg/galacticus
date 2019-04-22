# Contains a Perl module which provides state variables for component classes.

package Galacticus::Build::Components::Implementations::State;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Text::Template 'fill_in_string';
use List::ExtraUtils;
use Galacticus::Build::Components::Utils qw(&isIntrinsic);
use Galacticus::Build::Components::DataTypes;

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     implementationsState =>
     {
	 implementationIteratedFunctions =>
	     [
	      \&Implementation_State  ,
	      \&Implementation_Size_Of
	     ]
     }
    );

sub Implementation_State {
    # Generate variables which store active status for component implementations.
    my $build  = shift();
    my $class  = shift();
    my $member = shift();
    push
	(
	 @{$build->{'variables'}},
	 {
	     intrinsic  => "logical"                                                                                       ,
	     variables  => [ "nodeComponent".ucfirst($class->{'name'}).ucfirst($member->{'name'})."IsActiveValue=.false." ]
	 }
	);
}

sub Implementation_Size_Of {
    # Generate a function to return the size of the component implementation in bytes.
    my $build     = shift();
    $code::class  = shift();
    $code::member = shift();
    $code::implementationTypeName = "nodeComponent".ucfirst($code::class->{'name'}).ucfirst($code::member->{'name'});
    my $function =
    {
	type        => "integer(c_size_t)",
	name        => $code::implementationTypeName."SizeOf",
	description => "Return the size in bytes of a ".$code::member->{'name'}." implementation of the ".$code::class->{'name'}." component.",
	variables   =>
	    [
	     {
		 intrinsic  => "class",
		 type       => $code::implementationTypeName,
		 attributes => [ "intent(in   )" ],
		 variables  => [ "self" ]
	     }
	    ],
	content     => ""
    };
    # Initialize size to size of self.
    $function->{'content'} = $code::implementationTypeName."SizeOf=sizeof(self)\n";
    # Iterate over properties.
    my $loopIteratorRequired = 0;
    foreach $code::property ( grep {! $_->{'attributes'}->{'isVirtual'}} &List::ExtraUtils::hashList($code::member->{'properties'}->{'property'}) ) {
	if ( &isIntrinsic($code::property->{'data'}->{'type'}) ) {
	    # Add size of allocatable properties.
	    if ( $code::property->{'data'}->{'rank'} > 0 ) {
		$function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
{$implementationTypeName}SizeOf={$implementationTypeName}SizeOf+sizeof(self%{$property->{'name'}}Data)
CODE
	    }
	} else {
	    # Add non-static sizes of derived-type properties.
	    if ( $code::property->{'data'}->{'rank'} == 0 ) { 
	       $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
{$implementationTypeName}SizeOf={$implementationTypeName}SizeOf+self%{$property->{'name'}}Data%nonStaticSizeOf()
CODE
            } else {
               $loopIteratorRequired = 1;
	       $function->{'content'} .= fill_in_string(<<'CODE', PACKAGE => 'code');
do i=1,size(self%{$property->{'name'}}Data)
   {$implementationTypeName}SizeOf={$implementationTypeName}SizeOf+self%{$property->{'name'}}Data(i)%nonStaticSizeOf()
end do
CODE
           }
	}
    }
    push
    (
     @{$function->{'variables'}},
     {
	 intrinsic  => "integer",
	 type       => "c_size_t",
	 attributes => [ ],
         variables  => [ "i" ]
     }
    )
        if ( $loopIteratorRequired );
    # Insert a type-binding for this function.
    push(
	@{$build->{'types'}->{$code::implementationTypeName}->{'boundFunctions'}},
	{
	    type        => "procedure",
	    descriptor  => $function,
	    name        => "sizeOf", 
	}
	);
}

1;
