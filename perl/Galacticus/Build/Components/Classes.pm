# Contains a Perl module which handles component classes during build.

package Classes;
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
use Data::Dumper;
require List::ExtraUtils;
require Galacticus::Build::Components::Utils;

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     classes => {
	 gather =>
	     [
	      \&Gather_Classes
	     ]
     }
    );

sub Gather_Classes {
    # Determine class membership of each component implementation, and build a list of component class names.
    my $build = shift();
    # Iterate over all component instances.
    foreach my $implementation ( &ExtraUtils::hashList($build->{'components'}) ) {
	# Get the name of the component class.
	my $className          = $implementation->{'class'};
	# Get the name of the implementation.
	my $implementationName = $implementation->{'name' };
	# Append the component name to the list of members for its class.
	push(
	    @{$build->{'componentClasses'}->{$className}->{'members'}},
	    $implementationName
	    );
    }
    # Construct a list of component class names.
    @{$build->{'componentClassList'}} = &ExtraUtils::sortedKeys($build->{'componentClasses'});
    # Output report if sufficiently verbose.
    if ( $Utils::verbosityLevel > 0 ) {
	print "         --> Found the following component classes and implementations:\n";
	foreach my $className ( @{$build->{'componentClassList'}} ) {
	    print "            --> ".$className."\n";
	    foreach my $implementationName ( @{$build->{'componentClasses'}->{$className}->{'members'}} ) {
		print "               --> ".$implementationName."\n";
	    }
	}
    }
}

1;
