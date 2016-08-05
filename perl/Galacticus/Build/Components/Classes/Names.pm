# Contains a Perl module which provides naming functions for component classes.

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
     classNames =>
     {
	 functions =>
	     [
	      \&Class_Type
	     ]
     }
    );

sub Class_Type {
    # Generate functions to provide type names for component classes.
    my $build = shift();
    # Iterate over component classes.
    foreach $code::class ( &ExtraUtils::hashList($build->{'componentClasses'}) ) {
	my $classTypeName = "nodeComponent".ucfirst($code::class->{'name'});
	my $function =
	{
	    type        => "type(varying_string) => name",
	    name        => $classTypeName."Type",
	    description => "Returns the type name for the ".$code::class->{'name'}." component class.",
	    variables   =>
		[
		 {
		     intrinsic  => "class",
		     type       => $classTypeName,
		     attributes => [ "intent(in   )" ],
		     variables  => [ "self" ]
		 }
		]
	};    
	$function->{'content'}  = fill_in_string(<<'CODE', PACKAGE => 'code');
!GCC$ attributes unused :: self
name='nodeComponent:{$class->{'name'}}'
CODE
	# Add the function to the functions list.
	push(
	    @{$build->{'functions'}},
	    $function
	    );
	# Insert a type-binding for this function into the treeNode type.
	push(
	    @{$build->{'types'}->{$classTypeName}->{'boundFunctions'}},
	    {
		type        => "procedure", 
		name        => "type", 
		function    => $classTypeName."Type", 
		description => "Returns the type name for the ".$code::class->{'name'}." component class.",
		returnType  => "\\textcolor{red}{\\textless varying\\_string\\textgreater}", 
		arguments   => ""
	    }
	    );
    }
}

1;
