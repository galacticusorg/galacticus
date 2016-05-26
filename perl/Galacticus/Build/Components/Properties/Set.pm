# Contains a Perl module which handles setting of component properties during build.

package Set;
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
     propertiesSet => 
     {
	 functions =>
	     [
	      \&Build_Setters
	     ]
     }
    );

sub Build_Setters {
    # Validate that data type can be determined for each property.
    my $build = shift();
    # Iterate over components.
    foreach my $component ( &ExtraUtils::hashList($build->{'components'}) ) {
	# Iterate over all properties belonging to this component.	
	foreach my $property ( &ExtraUtils::hashList($component->{'properties'}->{'property'}, keyAs => 'name' ) ) {
	    # Insert an "isSettable" function into the base class.
	    my $functionName = $property->{'name'}."IsSettable";
	    unless ( grep {$_->{'name'} eq $functionName} @{$build->{'types'}->{'nodeComponent'.ucfirst($component->{'class'})}->{'boundFunctions'}} ) {
	    	push(
	    	    @{$build->{'types'}->{'nodeComponent'.ucfirst($component->{'class'})}->{'boundFunctions'}},
	    	    {
	    		type        => "procedure"    ,
	    		pass        => "nopass"       ,
	    		name        => $functionName  , 
	    		function    => "Boolean_False",
	    		returnType  => "\\logicalzero",
	    		arguments   => ""             ,
	    		description => "Specify whether the {\\normalfont \\ttfamily ".$property->{'name'}."} property of the {\\normalfont \\ttfamily ".$component->{'class'}."} component is settable."
	    	    }
	    	    );
	    }
	}
    }
}

1;
