# Contains a Perl module which provides state variables for component classes.

package Galacticus::Build::Components::Classes::State;
use strict;
use warnings;
use utf8;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use Text::Template 'fill_in_string';
use List::ExtraUtils;
use Galacticus::Build::Components::Utils;
use Galacticus::Build::Components::DataTypes;

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     classState =>
     {
	 classIteratedFunctions =>
	     [
	      \&Class_State
	     ]
     }
    );

sub Class_State {
    # Generate variables which store state for component classes.
    my $build    = shift();
    $code::class = shift();
    my $classTypeName = "nodeComponent".ucfirst($code::class->{'name'});
    push
	(
	 @{$build->{'variables'}},
	 # Object that will record which implementation of the class is to be used by default.
	 {
	     intrinsic  => "class"                                                  ,
	     type       => $classTypeName                                           ,
	     attributes => [ "public", "allocatable" ]                              ,
	     variables  => [ "default".ucfirst($code::class->{'name'})."Component" ]
	 },
	 # Object used to allocate component.
	 {
	     intrinsic  => "type"                                                   ,
	     type       => $classTypeName                                           ,
	     variables  => [ $code::class->{'name'}."Class" ]
	 }
	);
}

1;
