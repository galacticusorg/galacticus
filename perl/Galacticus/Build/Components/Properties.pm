# Contains a Perl module which handles component properties during build.

package Properties;
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
use List::Uniq ':all';
require List::ExtraUtils;
require Galacticus::Build::Components::Utils;

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     properties => 
     {
	 preValidate =>
	     [
	      \&Class_Defaults_Validate
	     ],
	 gather      =>
	     [
	      \&Class_Defaults_Gather
	     ],
	 scatter     =>
	     [
	      \&Class_Defaults_Scatter
	     ]
     }
    );

sub Class_Defaults_Validate {
    # Validate class default settings over all properties.
    my $build = shift();
    my $classDefaults;
    # Iterate over components.
    foreach my $component ( &ExtraUtils::hashList($build->{'components'}) ) {
	# Iterate over all properties belonging to this component.	
	foreach my $property ( &ExtraUtils::hashList($component->{'properties'}->{'property'}, keyAs => 'name' ) ) {
	    # Check for class defaults.
	    if ( exists($property->{'classDefault'}) ) {
		# Check that default code is the same for all implementations.
		my $code;
		if ( ref($property->{'classDefault'}) && exists($property->{'classDefault'}->{'content'}) ) {
		    $code = $property->{'classDefault'}->{'content'};
		} else {
		    $code = $property->{'classDefault'};
		}
		if ( exists($classDefaults->{$component->{'class'}}->{$property->{'name'}}->{'code'}) ) {
		    die
			(
			 "Class_Defaults: inconsistent class defaults for property '".
			 $property ->{'name' }                                       .
			 "' of class '"                                              .
			 $component->{'class'}                                       .
			 "'\n"
			)
			unless ( $code eq $classDefaults->{$component->{'class'}}->{$property->{'name'}}->{'code'} );
		} else {
		    $classDefaults->{$component->{'class'}}->{$property->{'name'}}->{'code'} = $code;
		}
		# Check that size of class defaults match.
		my $defaultCount;
		if ( ref($property->{'classDefault'}) && exists($property->{'classDefault'}->{'count'}) ) {
		    $defaultCount = $property->{'classDefault'}->{'count'};
		} elsif ( $property->{'classDefault'} =~ m/^\[.*\]$/ ) {
		    my @elements = split(/,/,$property->{'classDefault'});
		    $defaultCount = scalar(@elements);
		}
		if ( $defaultCount ) {
		    if ( exists($classDefaults->{$component->{'class'}}->{$property->{'name'}}->{'count'}) ) {
			die
			    (
			     "Class_Defaults: inconsistent class defaults count for property '".
			     $property ->{'name' }                                             .
			     "' of class '"                                                    .
			     $component->{'class'}                                             .
			     "'\n"
			    )
			    unless ( $defaultCount eq $classDefaults->{$component->{'class'}}->{$property->{'name'}}->{'count'} );
		    } else {
			$classDefaults->{$component->{'class'}}->{$property->{'name'}}->{'count'} = $defaultCount;
		    }
		}
	    }
	}
    }
}

sub Class_Defaults_Gather {
    # Gather class default settings.
    my $build = shift();
    # Iterate over components.
    foreach my $component ( &ExtraUtils::hashList($build->{'components'}) ) {
	# Iterate over all properties belonging to this component.	
	foreach my $property ( &ExtraUtils::hashList($component->{'properties'}->{'property'}, keyAs => 'name' ) ) {
	    # Check for class defaults.
	    if ( exists($property->{'classDefault'}) ) {
		# Get the default code.
		$build->{'classDefaults'}->{$component->{'class'}}->{$property->{'name'}}->{'code'} =
		    (
		     ref   ($property->{'classDefault'}             )
		     &&
		     exists($property->{'classDefault'}->{'content'}) 
		    ) 
		    ?
		    $property->{'classDefault'}->{'content'}
		    :
		    $property->{'classDefault'};
		# Get the default size.
		if ( ref($property->{'classDefault'}) && exists($property->{'classDefault'}->{'count'}) ) {
		    $build->{'classDefaults'}->{$component->{'class'}}->{$property->{'name'}}->{'count'} = $property->{'classDefault'}->{'count'};
		} elsif ( $property->{'classDefault'} =~ m/^\[.*\]$/ ) {
		    my @elements = split(/,/,$property->{'classDefault'});
		    $build->{'classDefaults'}->{$component->{'class'}}->{$property->{'name'}}->{'count'}  = scalar(@elements);
		}
		# Get any required modules.		
		@{$build->{'classDefaults'}->{$component->{'class'}}->{$property->{'name'}}->{'modules'}} = 
		    uniq(
			{sort => 1},
			(
			 @{$build->{'classDefaults'}->{$component->{'class'}}->{$property->{'name'}}->{'modules'}},
			 split(/\s*,\s*/,$property->{'classDefault'}->{'modules'})
			)
		    )
		    if ( ref($property->{'classDefault'}) && exists($property->{'classDefault'}->{'modules'}) );
	    }
	}
    }
}

sub Class_Defaults_Scatter {
    # Scatter class default settings to all components.
    my $build = shift();
    # Iterate over components.
    foreach my $component ( &ExtraUtils::hashList($build->{'components'}) ) {
	# Iterate over all properties belonging to this component.	
	foreach my $property ( &ExtraUtils::hashList($component->{'properties'}->{'property'}, keyAs => 'name' ) ) {
	    if ( exists($build->{'classDefaults'}->{$component->{'class'}}->{$property->{'name'}}) ) {
		print 
		    "         --> WARN: property '".
		    $property ->{'name' }                                                                        .
		    "' of component '"                                                                           .
		    $component->{'name' }                                                                        .
		    "' of class '"                                                                               .
		    $component->{'class'}                                                                        .
		    "' is being assigned a class default even though it does not have one explicitly declared\n"
		unless ( exists($property->{'classDefault'}) );
		$property->{'classDefault'} = $build->{'classDefaults'}->{$component->{'class'}}->{$property->{'name'}};
	    }
	}
    }
}

1;
