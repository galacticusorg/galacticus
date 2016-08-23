# Contains a Perl module which implements calculations of binned percentiles in weighted data.

package Stats::Percentiles;
use strict;
use warnings;
use PDL;
use PDL::NiceSlice;
use PDL::Ufunc;

sub BinnedPercentiles {
    # Distribute input data into specified bins, find the cumulative distribution of (weighted) values and determine the specified
    # percentiles of that distribution.

    # Get the arguments.
    my $binCenters  = shift;
    my $xValues     = shift;
    my $yValues     = shift;
    my $weights     = shift;
    my $percentiles = shift;
    my %options;
    if ( $#_ >= 1 ) {(%options) = @_};

    # Minimum number of points in a bin.
    my $binCountMinimum = exists($options{'binCountMinimum'}) ? $options{'binCountMinimum'} : 2;
    
    # Compute bin size.
    my $binWidth = ($binCenters->index(nelem($binCenters)-1)-$binCenters->index(0))/(nelem($binCenters)-1);

    # Compute bin ranges.
    my $binMinimum = $binCenters-0.5*$binWidth;
    my $binMaximum = $binCenters+0.5*$binWidth;

    # Create a PDL for results.
    my $results = pdl zeroes(nelem($binCenters),nelem($percentiles));

    # Loop through bins.
    for(my $iBin=0;$iBin<nelem($binCenters);++$iBin) {
	# Select properties in this bin.
	my $yValuesSelected = where($yValues,($xValues >= $binMinimum->index($iBin)) & ($xValues < $binMaximum->index($iBin)) );
	my $weightsSelected = where($weights,($xValues >= $binMinimum->index($iBin)) & ($xValues < $binMaximum->index($iBin)) );

        # Only compute results for cases where we have more than one entry.
	if ( nelem($yValuesSelected) >= $binCountMinimum ) {	

	    # Sort the selected values.
	    my $sortIndex = qsorti $yValuesSelected;
	    
	    # Get the cumulative weight and normalize to 100%.
	    my $cumulativeWeightsSelected  = cumusumover($weightsSelected->index($sortIndex));
	    $cumulativeWeightsSelected *= 100.0/$cumulativeWeightsSelected->index(nelem($cumulativeWeightsSelected)-1);
	    
	    # Interpolate to the desired percentiles to get the corresponding y values.
	    $results(($iBin),:) .= interpol($percentiles,$cumulativeWeightsSelected,$yValuesSelected->index($sortIndex));

	}

    }

    # Return the results.
    return $results;
}

1;
