# Contains a Perl module which handles serialization of component classes.

package Galacticus::Build::Components::Classes::Serialization;
use strict;
use warnings;
use utf8;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use Text::Template 'fill_in_string';
use List::ExtraUtils;
use Galacticus::Build::Components::Utils qw($fullyQualifiedNameLengthMax);
use Galacticus::Build::Components::DataTypes;

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     classesSerialization =>
     {
	 classIteratedFunctions =>
	     [
	      \&Class_Serialize_ASCII
	     ]
     }
    );

sub Class_Serialize_ASCII {
    # Generate a function to serialize component classes.
    my $build    = shift();
    $code::class = shift();
    my $function =
    {
	type        => "void",
	name        => "nodeComponent".ucfirst($code::class->{'name'})."SerializeASCII",
	description => "Serialize the content of a {\\normalfont \\ttfamily ".$code::class->{'name'}."} component to ASCII.",
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
	    name        => "serializeASCII", 
	}
	);
}

1;
