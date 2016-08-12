# Contains a Perl module which provides naming functions for component implementations.

package Names;
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
     implementationNames =>
     {
	 functions =>
	     [
	      \&Implementation_Type
	     ]
     }
    );

sub Implementation_Type {
    # Generate functions to provide type names for component implementations.
    my $build = shift();
    # Iterate over component classes.
    foreach $code::class ( &ExtraUtils::hashList($build->{'componentClasses'}) ) {
	# Iterate over class member implementations.
	foreach $code::member ( @{$code::class->{'members'}} ) {
	    my $implementationTypeName = "nodeComponent".ucfirst($code::class->{'name'}).ucfirst($code::member->{'name'});
	    my $function =
	    {
		type        => "type(varying_string) => name",
		name        => $implementationTypeName."Type",
		description => "Returns the type name for the ".$code::member->{'name'}." implementation of the ".$code::class->{'name'}." component class.",
		variables   =>
		    [
		     {
			 intrinsic  => "class",
			 type       => $implementationTypeName,
			 attributes => [ "intent(in   )" ],
			 variables  => [ "self" ]
		     }
		    ]
	    };    
	    $function->{'content'}  = fill_in_string(<<'CODE', PACKAGE => 'code');
!GCC$ attributes unused :: self
name='nodeComponent:{$class->{'name'}}:{$member->{'name'}}'
CODE
	    # Insert a type-binding for this function into the treeNode type.
	    push(
		@{$build->{'types'}->{$implementationTypeName}->{'boundFunctions'}},
		{
		    type        => "procedure", 
		    descriptor  => $function,
		    name        => "type", 
		    returnType  => "\\textcolor{red}{\\textless varying\\_string\\textgreater}", 
		    arguments   => ""
		}
		);
	}
    } 
}

1;
