# Contains a Perl module which implements calculations of binned means (and dispersions) of weighted data.

package Means;
use PDL;

my $status = 1;
$status;

sub BinnedMean {
    # Distribute input data into specified bins, find the total weight and the error.

    # Get the arguments.
    $binCenters  = shift;
    $xValues     = shift;
    $yValues     = shift;
    $weights     = shift;

    # Compute bin size.
    $binWidth = ($binCenters->index(nelem($binCenters)-1)-$binCenters->index(0))/(nelem($binCenters)-1);

    # Compute bin ranges.
    $binMinimum = $binCenters-0.5*$binWidth;
    $binMaximum = $binCenters+0.5*$binWidth;

    # Create a PDL for mean and dispersion.
    $mean            = pdl zeroes(nelem($binCenters));
    $meanError       = pdl zeroes(nelem($binCenters));
    $dispersion      = pdl zeroes(nelem($binCenters));
    $dispersionError = pdl zeroes(nelem($binCenters));

    # Loop through bins.
    for($iBin=0;$iBin<nelem($binCenters);++$iBin) {
	# Select properties in this bin.
	$yValuesSelected = where($yValues,$xValues >= $binMinimum->index($iBin) & $xValues < $binMaximum->index($iBin) );
	$weightsSelected = where($weights,$xValues >= $binMinimum->index($iBin) & $xValues < $binMaximum->index($iBin) );
        # Only compute results for cases where we have more than zero entries.
	if ( nelem($weightsSelected) == 1 ) {
	    $mean           ->index($iBin) .= sum($yValuesSelected);
	    $meanError      ->index($iBin) .= sum($yValuesSelected);
	    $dispersion     ->index($iBin) .= 0.0;
	    $dispersionError->index($iBin) .= 0.0;
	} elsif ( nelem($weightsSelected) > 1 ) {	    
	    # Compute weighted sums.
	    $sumYWeight  = sum($yValuesSelected   *$weightsSelected   );
	    $sumY2Weight = sum($yValuesSelected**2*$weightsSelected   );
	    $sumWeight   = sum(                    $weightsSelected   );
	    $sumWeight2  = sum(                    $weightsSelected**2);
	    # Compute mean and dispersion.
	    $mean           ->index($iBin) .= $sumYWeight/$sumWeight;
	    $meanError      ->index($iBin) .= ($sumY2Weight/$sumWeight-($sumYWeight/$sumWeight)**2)/($sumWeight**2/$sumWeight2);
	    if ( $meanError->index($iBin) > 0.0 ) {
		$meanError->index($iBin) .= sqrt($meanError->index($iBin));
	    } else {
		$meanError->index($iBin) .= 0.0;
	    }
	    $dispersion     ->index($iBin) .= $sumY2Weight/$sumWeight-($sumYWeight/$sumWeight)**2;
	    if ( $dispersion->index($iBin) > 0.0 ) {
		$dispersion->index($iBin) .= sqrt($dispersion->index($iBin));
	    } else {
		$dispersion->index($iBin) .= 0.0;
	    }
	    $dispersionError->index($iBin) .= $dispersion->index($iBin)*sqrt(2.0*$sumWeight2/$sumWeight**2);

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
