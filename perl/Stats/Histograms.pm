# Contains a Perl module which implements calculations of binned histograms of weighted data.

package Histograms;
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V092"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V092"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
}
unshift(@INC,$galacticusPath."perl"); 
use PDL;
use PDL::NiceSlice;
use Data::Dumper;

sub Histogram {
    # Distribute input data into specified bins, find the total weight and the error.

    # Get the arguments.
    my $binCenters  = shift;
    my $xValues     = shift;
    my $weights     = shift;
    my %options;
    if ( $#_ >= 1 ) {(%options) = @_};

    my $binMinimum = pdl zeroes(nelem($binCenters));
    my $binMaximum = pdl zeroes(nelem($binCenters));
    for(my $iBin=0;$iBin<nelem($binCenters);++$iBin) {
	if ( $iBin == 0 ) {
	    $binMinimum->(($iBin)) .= $binCenters->(($iBin))-0.5*($binCenters->(($iBin+1))-$binCenters->(($iBin)));
	} else {
	    $binMinimum->(($iBin)) .= 0.5*($binCenters->(($iBin))+$binCenters->(($iBin-1)));
	}
	if ( $iBin == nelem($binCenters)-1 ) {
	    $binMaximum->(($iBin)) .= $binCenters->(($iBin))+0.5*($binCenters->(($iBin))-$binCenters->(($iBin-1)));
	} else {
	    $binMaximum->(($iBin)) .= 0.5*($binCenters->(($iBin))+$binCenters->(($iBin+1)));
	}
    }
    my $binWidth = $binMaximum-$binMinimum;

    # Create a PDL for histogram and errors.
    my $histogram = pdl zeroes(nelem($binCenters));
    my $errors    = pdl zeroes(nelem($binCenters));
    my$covariance = pdl zeroes(nelem($binCenters),nelem($binCenters));

    # Method for constructing histogram depends on whether points are being smoothed.
    my $smooth       = 0;
    my $smoothRandom = 0;
    $smooth = 1
	if ( exists($options{'gaussianSmooth'}) );
    $smoothRandom = 1
	if ( exists($options{'gaussianSmoothRandom'}) );
    if ( $smooth == 1 && $smoothRandom == 0 ) {
	# Loop over points.
	my $sigma = $options{'gaussianSmooth'};
	for(my $i=0;$i<nelem($weights);++$i) {
	    my $fraction =
		(
		 +erf(($binMaximum-$xValues(($i)))/$sigma(($i))/sqrt(2.0))
		 -erf(($binMinimum-$xValues(($i)))/$sigma(($i))/sqrt(2.0))
		)/2.0;
	    $histogram  +=  $weights(($i))*$fraction    ;
	    $errors     += ($weights(($i))*$fraction)**2;
	    $covariance +=  $weights(($i))**2*(transpose($fraction) x $fraction);
	}
	$errors .= sqrt($errors);
    } else {
	# Use direct binning.
	# Add random offsets to x values if necessary.
	my $xValuesUsed = $xValues->copy();
	$xValuesUsed += grandom(nelem($xValuesUsed))*$options{'gaussianSmooth'}
	    if ( $smooth == 1 && $smoothRandom == 1 );
	# Loop through bins.
	for(my $iBin=0;$iBin<nelem($binCenters);++$iBin) {
	    # Select properties in this bin.
	    my $weightsSelected = where($weights,($xValuesUsed >= $binMinimum->index($iBin)) & ($xValuesUsed < $binMaximum->index($iBin)) );
	    # Only compute results for cases where we have at least one entry.
	    if ( nelem($weightsSelected) >= 1 ) {	
		
		# Sum up the weights in the bin.
		$histogram ->index($iBin) .= sum($weightsSelected);
		$errors    ->index($iBin) .= sqrt(sum($weightsSelected**2));
		$covariance->index($iBin) .= sum($weightsSelected**2);
	    } else {
		
		# No values in this bin - return all zeroes.
		$histogram ->index($iBin) .= 0.0;
		$errors    ->index($iBin) .= 0.0;
		$covariance->index($iBin) .= 0.0;
		
	    }
	    
	}

    }

    # Process the histogram according to any options specified.
    if ( exists($options{'normalized'}) && $options{'normalized'} == 1 ) {
	# Find the total weight.
	my $total = sum($weights);
	# Normalize the curve to unit area.
	$errors     /= $total;
	$histogram  /= $total;
	$covariance /= $total**2;
    }
    if ( exists($options{'differential'}) && $options{'differential'} == 1 ) {
	# Divide by the bin width to get a differential distribution.
	$errors     /= $binWidth;
	$histogram  /= $binWidth;
	$covariance /= $binWidth**2;
    }

    # Return the results.
    return ($histogram,$errors,$covariance);
}

1;
