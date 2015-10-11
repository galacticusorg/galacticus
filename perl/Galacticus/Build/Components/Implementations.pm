# Contains a Perl module which handles component implementations during build.

package Implementations;
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
use Sort::Topological qw(toposort);
require List::ExtraUtils;
require Galacticus::Build::Components::Utils;

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     implementations => 
     {
	 default =>
	     [
	      \&Null_Implementations
	     ],
	 gather =>
	     [
	      \&Implementation_Dependencies
	     ]
     }
    );

sub Null_Implementations {
    # Create a null implementation for each class.
    my $build = shift();
    # Iterate over components to determine which classes need a null case building.
    my %classes;
    foreach my $implementation ( &ExtraUtils::hashList($build->{'components'}) ) {
	# Initialize this class if it hasn't been seen before.
	unless ( exists($classes{$implementation->{'class'}}) ) {
	    $classes{$implementation->{'class'}}->{'hasNull'   } = 0;
	    $classes{$implementation->{'class'}}->{'hasDefault'} = 0;
	}
	# Record if a null component already exists.
	$classes{$implementation->{'class'}}->{'hasNull'   } = 1
	    if ( $implementation->{'name'     } eq "null" );
	# Record if a default is already specified.
	$classes{$implementation->{'class'}}->{'hasDefault'} = 1
	    if ( $implementation->{'isDefault'} eq "yes"  );
    }
    # Iterate over classes, creating null components as necessary.
    foreach my $class ( &ExtraUtils::sortedKeys(\%classes) ) {       
	# Test for pre-existing null component.
	if ( $classes{$class}->{'hasNull'} == 0 ) {
	    # No pre-existing null component is present, so simply insert one into the build data.
	    my $implementationName = ucfirst($class)."Null";
	    my $isDefault   = $classes{$class}->{'hasDefault'} ? "no" : "yes";
	    $build->{'components'}->{$implementationName}->{'class'    } = $class;
	    $build->{'components'}->{$implementationName}->{'name'     } = "null";
	    $build->{'components'}->{$implementationName}->{'isDefault'} = $isDefault;
	    # Append this new component ID to the component ID list.
	    push(@{$build->{'componentIdList'}},$implementationName);
	    # Display a message.
	    if ( $Utils::verbosityLevel >= 1 ) {
		print "         --> Adding null implementation ";
		print "as default "
		    unless ( $classes{$class}->{'hasDefault'} );
		print "for ".$class." class\n";
	    }
	} elsif ( $Utils::verbosityLevel >= 1 ) {
	    # Advise that null components don't need to be explicitly specified.
	    print "         --> INFO: a pre-existing null component was found for the '".$class."' class,\n";
	    print "                   but would be built automatically.\n";
	}
    }
}

sub Implementation_Dependencies {
    # Order class members such that parent classes come before child classes.
    my $build = shift();
    # Iterate over classes
    print "         --> Sorting implentations into parent->child order:\n"
	if ( $Utils::verbosityLevel >= 1 );
    foreach my $className ( @{$build->{'componentClassList'}} ) {
	my %dependencies;
	foreach my $implementationName ( @{$build->{'componentClasses'}->{$className}->{'members'}} ) {
	    my $implementationID = ucfirst($className).ucfirst($implementationName);
 	    my $implementation   = $build->{'components'}->{$implementationID};
	    push(@{$dependencies{$implementation->{'extends'}->{'name'}}},$implementationName)
		if ( exists($implementation->{'extends'}) );
	}
	@{$build->{'componentClasses'}->{$className}->{'members'}} =
	    toposort
	    (
	     sub { @{$dependencies{$_[0]} || []}; },
	     \@{$build->{'componentClasses'}->{$className}->{'members'}}
	    );
	if ( $Utils::verbosityLevel >= 1 ) {
	    print "            --> ".$className.":\n";
	    print "               --> ".$_."\n"
		foreach ( @{$build->{'componentClasses'}->{$className}->{'members'}} );
	}
    }
}

1;
