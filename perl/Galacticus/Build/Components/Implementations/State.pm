# Contains a Perl module which provides state variables for component classes.

package Galacticus::Build::Components::Implementations::State;
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
     implementationsState =>
     {
	 implementationIteratedFunctions =>
	     [
	      \&Implementation_State
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

1;
