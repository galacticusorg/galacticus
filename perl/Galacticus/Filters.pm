# Contains a Perl module which implements loading of filters for Galacticus analysis packages.

# Contributions to this file from: Andrew Benson; Christoph Behrens.

package Galacticus::Filters;
use strict;
use warnings;
use utf8;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use PDL;
use XML::Simple;
use Galacticus::Path;

sub Load {
    # Get the name of the requested filter.
    my $filterName        = shift();
    my $filterWavelengths          ;
    my $filterResponse             ; 
    # Determine filter type.
    if ( $filterName =~ m/^topHat_(\d+)_(\d+)$/ ) {
	# Handle top hat filters.
	my $wavelengthCentral = $1;
	my $resolution        = $2;
	my $cutOffResolution  = 1.0e4;
	$filterWavelengths    = pdl
	    [
	     $wavelengthCentral*(sqrt(4.0*$resolution**2+1.0)-1.0)/2.0/$resolution/(1.0+1.0/$cutOffResolution),
	     $wavelengthCentral*(sqrt(4.0*$resolution**2+1.0)-1.0)/2.0/$resolution                            ,
	     $wavelengthCentral*(sqrt(4.0*$resolution**2+1.0)+1.0)/2.0/$resolution                            ,
	     $wavelengthCentral*(sqrt(4.0*$resolution**2+1.0)+1.0)/2.0/$resolution*(1.0+1.0/$cutOffResolution) 
	    ];
	$filterResponse       = pdl [ 0.0, 1.0, 1.0, 0.0 ];
    } else {
	# Load a filter response from file.
	my $filterFile     = &galacticusPath()."data/filters/".$filterName.".xml";
	my $xml            = new XML::Simple;
	my $filter         = $xml->XMLin($filterFile);
	$filterWavelengths = pdl [];
	$filterResponse    = pdl [];
	foreach my $datum ( @{$filter->{'response'}->{'datum'}} ) {
	    $datum =~ s/^\s*//;
	    $datum =~ s/\s*$//;
	    my @columns = split(/\s+/,$datum);
	    $filterWavelengths = append($filterWavelengths,$columns[0]);
	    $filterResponse    = append($filterResponse   ,$columns[1]);
	}
    }
    # Return the filter response.
    return ($filterWavelengths,$filterResponse);
}

1;
