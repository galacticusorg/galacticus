# Contains a Perl module which provides utility functions related to component properties.

package Galacticus::Build::Components::Properties::Utils;
use strict;
use warnings;
use utf8;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use Sub::Identify ':all';
use Text::Template 'fill_in_string';
use Data::Dumper;
use List::ExtraUtils;
use Galacticus::Build::Components::Utils;
use Galacticus::Build::Components::DataTypes;

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     propertyUtils =>
     {
	 functions =>
	     [
	      \&Property_Function_Iterator
	     ]
     }
    );

# Adjectives associated with property attributes.
our %attributeAdjective =
    (
     get  => "isGettable" ,
     set  => "isSettable" ,
     rate => "isEvolvable"
    );

sub Property_Function_Iterator {
    # Iterates over component properties and calls registered functions for each property.
    my $build = shift();
    my @hooks = &List::ExtraUtils::hashList(\%Galacticus::Build::Component::Utils::componentUtils, keyAs => 'name');
    foreach my $hook ( @hooks ) {
	if ( exists($hook->{'propertyIteratedFunctions'}) ) {
	    my @functions = &List::ExtraUtils::as_array($hook->{'propertyIteratedFunctions'});
	    foreach my $function ( @functions ) {
		print "         --> ".$hook->{'name'}.(scalar(@functions) > 1 ? " {".sub_name($function)."}" : "")."\n";
		# Iterate over classes.
		foreach my $class ( &List::ExtraUtils::hashList($build->{'componentClasses'}) ) {
		    # Iterate over class member implementations.
		    foreach my $member ( @{$class->{'members'}} ) {
			# Iterate over all properties belonging to this member.	
			foreach my $property ( &List::ExtraUtils::hashList($member->{'properties'}->{'property'}, keyAs => 'name' ) ) {
			    &{$function}($build,$class,$member,$property);
			}
		    }
		}
	    }
	}
    }
}

1;
