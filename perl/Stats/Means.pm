# Contains a Perl module which implements calculations of binned means (and dispersions) of weighted data.

package Stats::Means;
use strict;
use warnings;
use PDL;
use PDL::NiceSlice;

sub AdaptiveBinnedMean {
    # Compute the mean of data binned into bins containing equal numbers of points.

    # Get the arguments.
    my $x           = shift();
    my $y           = shift();
    my $countPerBin = shift();
    my %options     = @_;

    # Construct output arrays.
    my $binCount    = int(nelem($x)/$countPerBin);
    ++$binCount
	unless ( nelem($x) % $countPerBin == 0 );
    my $xMean       = pdl [];
    my $yMean       = pdl [];
    my $yScatter    = pdl [];

    # Order the x values.
    my $xIndex      = $x->qsorti();

    # Iterate through bins, computing the means.
    my $jMinimum = -1;
    my $jMaximum = -1;
    while ($jMaximum < nelem($x)-1) {
	# Find minimum and maximum indices to include in this bin.
	$jMinimum = $jMaximum+1;
	$jMaximum = $jMinimum+$countPerBin-1;
	$jMaximum = nelem($x)-1
	    if ( $jMaximum > nelem($x)-1 );
	if ( exists($options{'binWidthMinimum'}) && $jMinimum > 0 ) {
	    while ( $jMaximum < nelem($x)-1 && $x->($xIndex)->(($jMaximum))-$x->($xIndex)->(($jMinimum)) < $options{'binWidthMinimum'} ) {
		++$jMaximum;	 
	    }
	}
	# Compute the means.
	$xMean    = $xMean   ->append($x->($xIndex)->($jMinimum:$jMaximum)->average());
	$yMean    = $yMean   ->append($y->($xIndex)->($jMinimum:$jMaximum)->average());
	# Compute the root-variance.
	$yScatter = $yScatter->append(sqrt(average(($y->($xIndex)->($jMinimum:$jMaximum)-$yMean->((-1)))**2)));
    }

    # Return the results.
    return ($xMean, $yMean, $yScatter);
}

sub BinnedMean {
    # Distribute input data into specified bins, find the total weight and the error.

    # Get the arguments.
    my $binCenters  = shift;
    my $xValues     = shift;
    my $yValues     = shift;
    my $weights     = shift;

    # Validate input.
    die("Stats::Means::BinnedMean: xValues must be a 1D array")
	unless ( $xValues->ndims() == 1 );
    die("Stats::Means::BinnedMean: yValues must be a 1D array")
	unless ( $yValues->ndims() == 1 );
    die("Stats::Means::BinnedMean: weights must be a 1D array")
	unless ( $weights->ndims() == 1 );
    die("Stats::Means::BinnedMean: yValues must have same number of elements as xValues")
	unless ( $yValues->dim(0) == $xValues->dim(0) );
    die("Stats::Means::BinnedMean: weights must have same number of elements as xValues")
	unless ( $weights->dim(0) == $xValues->dim(0) );
    
    # Compute bin size.
    my $binWidth = ($binCenters->index(nelem($binCenters)-1)-$binCenters->index(0))/(nelem($binCenters)-1);

    # Compute bin ranges.
    my $binMinimum = $binCenters-0.5*$binWidth;
    my $binMaximum = $binCenters+0.5*$binWidth;

    # Create a PDL for mean and dispersion.
    my $mean            = pdl zeroes(nelem($binCenters));
    my $meanError       = pdl zeroes(nelem($binCenters));
    my $dispersion      = pdl zeroes(nelem($binCenters));
    my $dispersionError = pdl zeroes(nelem($binCenters));

    # Loop through bins.
    for(my $iBin=0;$iBin<nelem($binCenters);++$iBin) {
	# Select properties in this bin.
	my $yValuesSelected = where($yValues,($xValues >= $binMinimum->index($iBin)) & ($xValues < $binMaximum->index($iBin)) );
	my $weightsSelected = where($weights,($xValues >= $binMinimum->index($iBin)) & ($xValues < $binMaximum->index($iBin)) );
        # Only compute results for cases where we have more than zero entries.
	if ( nelem($weightsSelected) == 1 ) {
	    $mean           ->index($iBin) .= sum($yValuesSelected);
	    $meanError      ->index($iBin) .= sum($yValuesSelected);
	    $dispersion     ->index($iBin) .= 0.0;
	    $dispersionError->index($iBin) .= 0.0;
	} elsif ( nelem($weightsSelected) > 1 ) {	    
	    # Compute weighted sums.
	    my $sumYWeight  = sum($yValuesSelected   *$weightsSelected   );
	    my $sumY2Weight = sum($yValuesSelected**2*$weightsSelected   );
	    my $sumWeight   = sum(                    $weightsSelected   );
	    my $sumWeight2  = sum(                    $weightsSelected**2);
	    # Compute mean and dispersion.
	    $mean           ->index($iBin) .= $sumYWeight/$sumWeight;
	    # Source: http://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Statistical_properties
	    $meanError      ->index($iBin) .= ($sumY2Weight/$sumWeight-($sumYWeight/$sumWeight)**2)/($sumWeight**2/$sumWeight2);
	    if ( $meanError->index($iBin) > 0.0 ) {
		$meanError->index($iBin) .= sqrt($meanError->index($iBin));
	    } else {
		$meanError->index($iBin) .= 0.0;
	    }
	    # Source: http://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Weighted_sample_variance
	    $dispersion     ->index($iBin) .= ($sumY2Weight+$mean->index($iBin)**2*$sumWeight-2.0*$mean->index($iBin)*$sumYWeight)*$sumWeight/($sumWeight**2-$sumWeight2);
	    if ( $dispersion->index($iBin) > 0.0 ) {
		$dispersion->index($iBin) .= sqrt($dispersion->index($iBin));
	    } else {
		$dispersion->index($iBin) .= 0.0;
	    }
	    # Source: http://mcs.une.edu.au/~stat354/notes/node63.html
	    $dispersionError->index($iBin) .= $dispersion->index($iBin)
		*0.5
		*sqrt(
		    2.0
		    *$sumWeight**2/(1.0*$sumWeight**2-$sumWeight2)
		    *$sumWeight**2/(2.0*$sumWeight**2-$sumWeight2)
		);

	} else {

	    # No values in this bin - return all zeroes.
	    $mean           ->index($iBin) .= 0.0;
	    $meanError      ->index($iBin) .= 0.0;
	    $dispersion     ->index($iBin) .= 0.0;
	    $dispersionError->index($iBin) .= 0.0;

	}

    }

    # Return the results.
    return ($mean,$meanError,$dispersion,$dispersionError);
}

1;
