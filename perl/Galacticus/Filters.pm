# Contains a Perl module which implements loading of filters for Galacticus analysis packages.

# Contributions to this file from: Andrew Benson; Christoph Behrens.

package Filters;
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V094"}) ) {
    $galacticusPath = $ENV{"GALACTICUS_ROOT_V094"};
    $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
    $galacticusPath = "./";
}
use utf8;
use PDL;
use XML::Simple;

sub Load {
    # Get the name of the requested filter.
    my $filterName = shift;
    # Load a filter response from file.
    my $filterFile        = $galacticusPath."data/filters/".$filterName.".xml";
    my $xml               = new XML::Simple;
    my $filter            = $xml->XMLin($filterFile);
    my $filterWavelengths = pdl [];
    my $filterResponse    = pdl [];
    foreach my $datum ( @{$filter->{'response'}->{'datum'}} ) {
	$datum =~ s/^\s*//;
	$datum =~ s/\s*$//;
	my @columns = split(/\s+/,$datum);
	$filterWavelengths = append($filterWavelengths,$columns[0]);
	$filterResponse    = append($filterResponse   ,$columns[1]);
    }
    # Return the filter response.
    return ($filterWavelengths,$filterResponse);
}

1;
