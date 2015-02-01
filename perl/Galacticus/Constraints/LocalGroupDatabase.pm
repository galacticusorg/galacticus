# Contains a Perl module which provides an interface to the Local Group database file.

package LocalGroup;
use strict;
use warnings;
use XML::Simple;
use PDL;
use PDL::NiceSlice;
use Data::Dumper;

sub Select {
    # Get arguments.
    my @properties = @{shift()};
    my %options    = @_;
    # Specify defaults.
    $options{'badValue'      } = pdl -99
	unless ( exists($options{'badValue'      }) );
    $options{'excludeMissing'} = ""
	unless ( exists($options{'excludeMissing'}) );
    $options{'excludeCentral'} = 0
	unless ( exists($options{'excludeCentral'}) );
    # Encode memberships.
    $options{'membership'} = "G"
	if ( $options{'membership'} eq "Milky Way" || $options{'membership'} eq "MilkyWay" );
    $options{'membership'} = "A"
	if ( $options{'membership'} eq "M31"       );
    # Parse the database.
    my $xml        = new XML::Simple();
    my $localGroup = $xml->XMLin("data/observations/localGroup/localGroupSatellites.xml", KeyAttr => []);
    # Initialize result vectors.
    my @results;
    # Initialize exclusion vector.
    my $include = pdl ones(scalar(@{$localGroup->{'galaxies'}->{'galaxy'}}));
    # Define property attributes.
    my @attributes = ( "value", "error", "errorLow", "errorHigh" );
    # Iterate over properties.
    foreach my $property ( @properties ) {
	my $result;
	$result->{'property'} = $property;
	if ( $property ne "name" ) {
	    $result->{$_} = pdl []
		foreach ( @attributes );
	}
	# Iterate over galaxies.
	my $galaxyCounter = -1;
	foreach my $galaxy ( @{$localGroup->{'galaxies'}->{'galaxy'}} ) {
	    ++$galaxyCounter;
	    # Check membership.
	    if ( exists($options{'membership'}) ) {
		if ( exists($galaxy->{'membership'}->{'value'}) ) {
		    $include->(($galaxyCounter)) .= 0
			unless ( $galaxy->{'membership'}->{'value'} =~ m/$options{'membership'}/ );
		} else {
		    $include->(($galaxyCounter)) .= 0;
		}
	    }
	    # Check for central.
	    if ( $options{'excludeCentral'} == 1 ) {
		$include->(($galaxyCounter)) .= 0
		    if ( 
			$galaxy->{'name'} eq "The Galaxy"
			||
			$galaxy->{'name'} eq "Andromeda"
		    );
	    }
	    # Handle name as a special case.
	    if ( $property eq "name" ) {
		push(@{$result->{'value'}},$galaxy->{'name'});
	    } else {
		# Extract property for this galaxy.
		my $galaxyData;
		foreach ( @attributes ) {
		    if ( exists($galaxy->{$property}->{$_}) && $galaxy->{$property}->{$_} !~ m/\?/ ) {
			$galaxyData->{$_} = $galaxy->{$property}->{$_}
		    } else {
			$galaxyData->{$_} = undef;
		    }
		}
		# Set symmetric error.
		$galaxyData->{'error'} = 0.5*($galaxyData->{'errorLow'}+$galaxyData->{'errorHigh'})
		    if ( ! defined($galaxyData->{'error'}) && defined($galaxyData->{'errorLow'}) && defined($galaxyData->{'errorHigh'}) );
		# Check for missing property.
		if ( $options{'excludeMissing'} eq $property ) {
		    $include->(($galaxyCounter)) .= 0
			unless ( defined($galaxyData->{'value'}) && defined($galaxyData->{'error'}) );
		}
		# Accumulate data.
		foreach ( @attributes ) {
		    if ( defined($galaxyData->{$_}) ) {
			$result->{$_} = $result->{$_}->append($galaxyData->{$_}   );
		    } else {
			$result->{$_} = $result->{$_}->append($options{'badValue'});
		    }
		}
	    }
	}
	push(@results,$result);
    }
    # Exclude galaxies.
    my $includeWhich = which($include == 1);
    foreach my $result ( @results ) {
	if ( $result->{'property'} eq "name" ) {
	    my @includedNames;
	    for(my $i=0;$i<nelem($includeWhich);++$i) {
		push(@includedNames,$result->{'value'}->[$includeWhich->(($i))]);		
	    }
	    @{$result->{'value'}} = @includedNames;
	} else {
	    foreach ( @attributes ) {
		my $includedData = $result->{$_}->($includeWhich);
		$result->{$_} = $includedData;
	    }
	}
    }
    # Add property definitions.
    foreach my $result ( @results ) {
	foreach my $property ( @{$localGroup->{'properties'}} ) {
	    $result->{'meta'} = $property
		if ( $property->{'name'} eq $result->{'property'} );
	}    
    }
    # Return the results.
    return @results;
}

1;
