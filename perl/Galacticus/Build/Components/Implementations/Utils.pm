# Contains a Perl module which provides utility functions related to component implementations.

package Galacticus::Build::Components::Implementations::Utils;
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
use Data::Dumper;
require List::ExtraUtils;
require Galacticus::Build::Components::Utils;
require Galacticus::Build::Components::DataTypes;

# Insert hooks for our functions.
%Galacticus::Build::Component::Utils::componentUtils = 
    (
     %Galacticus::Build::Component::Utils::componentUtils,
     implementationUtils =>
     {
	 functions =>
	     [
	      \&Implementation_Is_Active
	     ]
     }
    );

sub Implementation_Is_Active {
    # Generate functions which return true if a component implementation is active.
    my $build = shift();
    # Iterate over component classes.
    foreach $code::class ( &ExtraUtils::hashList($build->{'componentClasses'}) ) {
	# Iterate over class member implementations.
	foreach $code::member ( @{$code::class->{'members'}} ) {
	    $code::implementationTypeName = "nodeComponent".ucfirst($code::class->{'name'}).ucfirst($code::member->{'name'});
	    my $function =
	    {
		type        => "logical",
		name        => $code::implementationTypeName."IsActive",
		description => "Return true if the ".$code::member->{'name'}." implementation of the ".$code::class->{'name'}." component is the active choice."
	    };
	    # Build the function.
	    $function->{'content'} = fill_in_string(<<'CODE', PACKAGE => 'code');
{$implementationTypeName}IsActive={$implementationTypeName}IsActiveValue
CODE
	    # Insert a type-binding for this function into the treeNode type.
	    push(
		@{$build->{'types'}->{"nodeComponent".ucfirst($code::class->{'name'})}->{'boundFunctions'}},
		{
		    type        => "procedure", 
		    descriptor  => $function,
		    pass        => "nopass",
		    name        => $code::member->{'name'}."IsActive"
		}
		);	    
	}
    }
}

sub hasRealEvolvers {
    # Returns true if an implementation has real (i.e. non-virtual), evolvable properties.
    my $member = shift();
    return
	grep
         {	
	   ! $_->{'attributes'}->{'isVirtual'  }
	  &&	    
	     $_->{'data'      }->{'isEvolvable'}
         }
         &ExtraUtils::hashList($member->{'properties'}->{'property'});	
}

sub hasRealNonTrivialEvolvers {
    # Returns true if an implementation has real (i.e. non-virtual), evolvable, non-trivial (i.e. not rank-0, doubles) properties.
    my $member = shift();
    return
	grep
         {	
	   !   $_->{'attributes'}->{'isVirtual'  }
	  &&	    
	       $_->{'data'      }->{'isEvolvable'}
	  &&
             (
	       $_->{'data'      }->{'rank'       } >  0 
	      ||
	       $_->{'data'      }->{'type'       } ne "double"
	     )			
         }
         &ExtraUtils::hashList($member->{'properties'}->{'property'});	
}

1;
