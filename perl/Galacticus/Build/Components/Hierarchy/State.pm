# Contains a Perl module which provides state for the hierarchy.

package Galacticus::Build::Components::Hierarchy::State;
use strict;
use warnings;
use utf8;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'}."/perl";
use Text::Template 'fill_in_string';
use List::ExtraUtils;
use Galacticus::Build::Components::Utils;
use Galacticus::Build::Components::DataTypes;

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     hierarchyState =>
     {
	 functions =>
	     [
	      \&Hierarchy_State
	     ]
     }
    );

sub Hierarchy_State {
    # Insert a variable which stores the initialization state of the class hierarchy.
    my $build = shift();
    push
	(
	 @{$build->{'variables'}},
	 {
	     intrinsic  => "integer",
	     variables  => [ "hierarchyInitialized=0" ]
	 }
	);
}

1;
