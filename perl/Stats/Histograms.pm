# Contains a Perl module which implements calculations of binned histograms of weighted data.

package Histograms;
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V093"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V093"};
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
	if ( exists($options{'binWidths'}) ) {
	    # Bin widths were specified - use them.
	    $binMinimum->(($iBin)) .= $binCenters->(($iBin))-0.5*$options{'binWidths'}->(($iBin));
	    $binMaximum->(($iBin)) .= $binCenters->(($iBin))+0.5*$options{'binWidths'}->(($iBin));
	} else {
	    # Bin widths were not specified - assume bin boundaries are midway between adjacent bin centers and extrapolate the
	    # bounaries of the edge bins.
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
    }
    my $binWidth = $binMaximum-$binMinimum;

    # Create a PDL for histogram and errors.
    my $histogram  = pdl zeroes(nelem($binCenters));
    my $errors     = pdl zeroes(nelem($binCenters));
    my $covariance = pdl zeroes(nelem($binCenters),nelem($binCenters));

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
	    my $weightsSelected = where($weights,($xValuesUsed >= $binMinimum->(($iBin))) & ($xValuesUsed < $binMaximum->(($iBin))) );
	    # Only compute results for cases where we have at least one entry.
	    if ( nelem($weightsSelected) >= 1 ) {	
		# Sum up the weights in the bin.
		$histogram ->(($iBin)        ) .= sum($weightsSelected);
		$errors    ->(($iBin)        ) .= sqrt(sum($weightsSelected**2));
		$covariance->(($iBin),($iBin)) .= sum($weightsSelected**2);
	    } else {
		# No values in this bin - return all zeroes.
		$histogram ->(($iBin)        ) .= 0.0;
		$errors    ->(($iBin)        ) .= 0.0;
		$covariance->(($iBin),($iBin)) .= 0.0;	
	    }   
	}	
    }
    
    # Process the histogram according to any options specified.
    if ( exists($options{'normalized'}) && $options{'normalized'} == 1 ) {
	# Find the total weight.
	my $total;
	if ( exists($options{'normalizeBy'}) && $options{'normalizeBy'} eq "histogram" ) {
	    $total = sum($histogram);   
	} else {
	    $total = sum($weights  );
	}
	# Normalize the curve to unit area.
	$errors     /= $total;
	$histogram  /= $total;
	# Dividing through by the sum of the histogram correlates the bins, leading to a modification of the covariance
	# matrix. Account for this by constructing the Jacobian of the transformation from the unnormalized to normalized
	# histogram and using this to propagate the original covariance matrix.
	my $jacobian = pdl zeroes(nelem($histogram),nelem($histogram));
	for(my $i=0;$i<nelem($histogram);++$i) {
	    $jacobian->(($i),  : ) .= -$histogram;
	    $jacobian->(($i),($i)) += 1.0;
	}
	$jacobian /= $total;
	my $covarianceCorrelated = $jacobian x $covariance x transpose($jacobian);
	$covariance .= $covarianceCorrelated;
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

sub Histogram2D {
    # Distribute input data into specified bins in a 2D grid, find the total weight and the error.
    # Get the arguments.
    my $xBinCentersIn = shift;
    my $yBinCenters   = shift;
    my $xValues       = shift;
    my $yValues       = shift;
    my $weights       = shift;
    my %options;
    if ( $#_ >= 1 ) {(%options) = @_};
    # If x-bins are indepenent of y, construct a matrix of x-bin centers.
    my $xBinCenters;
    if ( $xBinCentersIn->ndims() == 2 ) {
	$xBinCenters = $xBinCentersIn;
    } else {
	$xBinCenters = pdl zeroes($xBinCentersIn->dim(0),$yBinCenters->dim(0));
	for(my $i=0;$i<$yBinCenters->dim(0);++$i) {
	    $xBinCenters->(:,($i)) .= $xBinCentersIn;
	}
    }
    # Construct x-bins.
    my $xBinMinimum = pdl zeroes($xBinCenters->dim(0),$yBinCenters->dim(0));
    my $xBinMaximum = pdl zeroes($xBinCenters->dim(0),$yBinCenters->dim(0));
    for(my $j=0;$j<$yBinCenters->dim(0);++$j) {
	for(my $i=0;$i<$xBinCenters->dim(0);++$i) {
	    if ( $i == 0 ) {
		$xBinMinimum->(($i),($j)) .= $xBinCenters->(($i),($j))-0.5*($xBinCenters->(($i+1),($j))-$xBinCenters->(($i),($j)));
	    } else {
		$xBinMinimum->(($i),($j)) .= 0.5*($xBinCenters->(($i),($j))+$xBinCenters->(($i-1),($j)));
	    }
	    if ( $i == $xBinCenters->dim(0)-1 ) {
		$xBinMaximum->(($i),($j)) .= $xBinCenters->(($i),($j))+0.5*($xBinCenters->(($i),($j))-$xBinCenters->(($i-1),($j)));
	    } else {
		$xBinMaximum->(($i),($j)) .= 0.5*($xBinCenters->(($i),($j))+$xBinCenters->(($i+1),($j)));
	    }
	}
    }
    my $xBinWidth = $xBinMaximum-$xBinMinimum; 
    my $xBinCount = $xBinCenters->dim(0);
    # Construct y-bins.
    my $yBinMinimum = pdl zeroes($yBinCenters->dim(0));
    my $yBinMaximum = pdl zeroes($yBinCenters->dim(0));
	for(my $j=0;$j<$yBinCenters->dim(0);++$j) {
	    if ( $j == 0 ) {
	    $yBinMinimum->(($j)) .= $yBinCenters->(($j))-0.5*($yBinCenters->(($j+1))-$yBinCenters->(($j)));
	} else {
	    $yBinMinimum->(($j)) .= 0.5*($yBinCenters->(($j))+$yBinCenters->(($j-1)));
	}
	if ( $j == $yBinCenters->dim(0)-1 ) {
	    $yBinMaximum->(($j)) .= $yBinCenters->(($j))+0.5*($yBinCenters->(($j))-$yBinCenters->(($j-1)));
	} else {
	    $yBinMaximum->(($j)) .= 0.5*($yBinCenters->(($j))+$yBinCenters->(($j+1)));
	}
    }
	my $yBinWidth = $yBinMaximum-$yBinMinimum;
    my $yBinCount = $yBinCenters->dim(0); 
    # Check that x bin array size is consistent with y bin array size.
    die('Histograms::Histogram2D: size of x-bin centers is mismatched with size of y-bin centers')
	unless ( $xBinCenters->dim(1) == $yBinCenters->dim(0) );
    # Set covariance model and check all required information is available.
    my $covarianceModel;
    $covarianceModel->{'type'} = "poisson";
    $covarianceModel->{'type'} = $options{'covarianceModel'}
       if ( exists($options{'covarianceModel'}) );
    if ( $covarianceModel->{'type'} eq "binomial" ) {
	foreach ( "mainBranchStatus", "nodeMass", "haloMassBinsPerDecade", "haloMassBinsMinimum", "haloMassBinsMaximum" ) {
	    die("Histograms::Histogram2D: '".$_."' option must be specified for 'binomial' covariance model")
		unless ( exists($options{$_}) );
	}
	# Construct a set of halo mass bins.
	$covarianceModel->{'haloMassMinimumLogarithmic'        } =      log10($options{'haloMassBinsMinimum'});
	$covarianceModel->{'haloMassBinsCount'                 } = long(log10($options{'haloMassBinsMaximum'}/$options{'haloMassBinsMinimum'})*$options{'haloMassBinsPerDecade'}+0.5);
	$covarianceModel->{'haloMassIntervalLogarithmicInverse'} = $covarianceModel->{'haloMassBinsCount'}/log10($options{'haloMassBinsMaximum'}/$options{'haloMassBinsMinimum'});
	$covarianceModel->{'mainBranchGalaxyWeights'           } = pdl zeroes($xBinCount*$yBinCount,$covarianceModel->{'haloMassBinsCount'});
	$covarianceModel->{'mainBranchGalaxyWeightsSquared'    } = pdl zeroes($xBinCount*$yBinCount,$covarianceModel->{'haloMassBinsCount'});
    }
    # Create a PDL for histogram and errors.
    my $histogram  = pdl zeroes($xBinCount           ,           $yBinCount);
    my $errors     = pdl zeroes($xBinCount           ,           $yBinCount);
    my $covariance = pdl zeroes($xBinCount*$yBinCount,$xBinCount*$yBinCount);
    # Create PDLs for total weight applied to each row and each column.
    my $rowWeights    = pdl zeroes($xBinCount,$yBinCount);
    my $columnWeights = pdl zeroes(           $yBinCount);
    # Method for constructing histogram depends on whether points are being smoothed.
    my $smooth       = 0;
    my $smoothRandom = 0;
    $smooth = 1
	if ( exists($options{'gaussianSmoothX'}) );
    die("Histogram2D: gaussianSmoothY must be specified if gaussianSmoothX is given")
	if ( $smooth == 1 && ! exists($options{'gaussianSmoothY'}) );
    $smoothRandom = 1
	if ( exists($options{'gaussianSmoothRandom'}) );
    if ( $smooth == 1 && $smoothRandom == 0 ) {
	# Loop over points.
	my $sigmaX = $options{'gaussianSmoothX'};
	my $sigmaY = $options{'gaussianSmoothY'};
	for(my $i=0;$i<nelem($weights);++$i) {
	    my $fractionX =
		(
		 +erf(($xBinMaximum-$xValues(($i)))/$sigmaX(($i))/sqrt(2.0))
		 -erf(($xBinMinimum-$xValues(($i)))/$sigmaX(($i))/sqrt(2.0))
		)/2.0;
	    my $nonPositiveX = which($fractionX < 0.0);
	    $fractionX->flat()->($nonPositiveX) .= 0.0;
	    my $fractionY =
		(
		 +erf(($yBinMaximum-$yValues(($i)))/$sigmaY(($i))/sqrt(2.0))
		 -erf(($yBinMinimum-$yValues(($i)))/$sigmaY(($i))/sqrt(2.0))
		)/2.0;
	    my $nonPositiveY = which($fractionY < 0.0);
	    $fractionY->($nonPositiveY) .= 0.0;
	    my $fraction    = $fractionX->copy();
	    for(my $j=0;$j<$yBinCenters->dim(0);++$j) {
		$fraction->(:,($j)) *= $fractionY->(($j));
	    }
	    $rowWeights    +=  $weights(($i))*$fractionX;
	    $columnWeights +=  $weights(($i))*$fractionY;
	    $histogram     +=  $weights(($i))*$fraction    ;
	    $errors        += ($weights(($i))*$fraction)**2;
	    my $thisCovarianceType = "poisson";
	    $thisCovarianceType = "binomial"
		if ( $covarianceModel->{'type'} eq "binomial" && $options{'mainBranchStatus'}->(($i)) eq 1 );
	    if ( $thisCovarianceType eq "poisson" ) {
		$covariance    +=  $weights(($i))**2*(transpose($fraction->flat()) x $fraction->flat());
	    } elsif ( $thisCovarianceType eq "binomial" ) {
		my $haloMassBin = floor((log10($options{'nodeMass'}->(($i)))-$covarianceModel->{'haloMassMinimumLogarithmic'})*$covarianceModel->{'haloMassIntervalLogarithmicInverse'});
		if ( $haloMassBin >= 1 && $haloMassBin <= $covarianceModel->{'haloMassBinsCount'}) {
		    $covarianceModel->{'mainBranchGalaxyWeights'       }->(:,($haloMassBin)) +=  $weights(($i))*$fraction->flat()    ;
		    $covarianceModel->{'mainBranchGalaxyWeightsSquared'}->(:,($haloMassBin)) += ($weights(($i))*$fraction->flat())**2;
		}
	    }
	}
	$errors .= sqrt($errors);
    } else {
	# Use direct binning.
	# Add random offsets to x and y values if necessary.
	my $xValuesUsed = $xValues->copy();
	my $yValuesUsed = $yValues->copy();
	if ( $smooth == 1 && $smoothRandom == 1 ) {
	    $xValuesUsed += grandom(nelem($xValuesUsed))*$options{'gaussianSmoothX'};
	    $yValuesUsed += grandom(nelem($yValuesUsed))*$options{'gaussianSmoothY'};
	}
	# Loop through bins.
	for(my $i=0;$i<$xBinCenters->dim(0);++$i) {
	    for(my $j=0;$j<$yBinCenters->dim(0);++$j) {
		my $k = $i+$j*$xBinCount;
		# Select properties in this bin.
		my $weightsSelected = 
		    where
		    (
		     $weights,
		     ($xValuesUsed >= $xBinMinimum->(($i),($j))) & ($xValuesUsed < $xBinMaximum->(($i),($j)))
		     &
		     ($yValuesUsed >= $yBinMinimum->(     ($j))) & ($yValuesUsed < $yBinMaximum->(     ($j)))
		    );
		# Only compute results for cases where we have at least one entry.
		if ( nelem($weightsSelected) >= 1 ) {	
		    # Sum up the weights in the bin.
		    $histogram    ->(($i),($j)) .=      sum($weightsSelected   ) ;
		    $errors       ->(($i),($j)) .= sqrt(sum($weightsSelected**2));
		    $covariance   ->(($k),($k)) .=      sum($weightsSelected**2) ;
		} else {		
		    # No values in this bin - return all zeroes.
		    $histogram ->(($i),($j)) .= 0.0;
		    $errors    ->(($i),($j)) .= 0.0;
		    $covariance->(($k),($k)) .= 0.0;
		}
	    }
	}
	for(my $j=0;$j<$yBinCenters->dim(0);++$j) {
	    my $weightsSelected = 
		where
		(
		 $weights,
		 ($yValuesUsed >= $yBinMinimum->(($j))) & ($yValuesUsed < $yBinMaximum->(($j)))
		);
	    $columnWeights->(($j)) += sum($weightsSelected) ;
	    for(my $i=0;$i<$xBinCenters->dim(0);++$i) {
		my $weightsSelected = 
		    where
		    (
		     $weights,
		     ($xValuesUsed >= $xBinMinimum->(($i),($j))) & ($xValuesUsed < $xBinMaximum->(($i),($j)))
		    );
		$rowWeights->(($i),($j)) += sum($weightsSelected) ;
	    }
	}
    }
    # Complete the covariance model.
    if ( $covarianceModel->{'type'} eq "binomial" ) {
	# Add the contribution from main branch galaxies to the covariance matrix.
	for(my $m=0;$m<$covarianceModel->{'haloMassBinsCount'};++$m) {
	    my $haloWeightBinTotal = sum($covarianceModel->{'mainBranchGalaxyWeights'}->(:,($m)));
	    if ( $haloWeightBinTotal > 0.0 ) {
                for(my $i=0;$i<nelem($histogram);++$i) {
		    $covariance->(($i),($i)) += 
			+(1.0-$covarianceModel->{'mainBranchGalaxyWeights'       }->(($i),($m))/$haloWeightBinTotal)
			*     $covarianceModel->{'mainBranchGalaxyWeightsSquared'}->(($i),($m));
		    for(my $j=0;$j<nelem($histogram);++$j) {
			next
			    if ( $i == $j );
			$covariance->(($i),($j)) += 
			    -($covarianceModel->{'mainBranchGalaxyWeights'       }->(($j),($m))/$haloWeightBinTotal) 
			    * $covarianceModel->{'mainBranchGalaxyWeightsSquared'}->(($i),($m));
		    }
		}
	    }
	}
    }
    # Process the histogram according to any options specified.
    if ( exists($options{'normalized'}) ) {
	if ( $options{'normalized'} eq "xy" ) {
	    # Find the total weight.
	    my $total;
	    if ( exists($options{'normalizeBy'}) && $options{'normalizeBy'} eq "histogram" ) {
		$total = sum($histogram);   
	    } else {
		$total = sum($weights  );
	    }
	    $total = 1.0
		if ( $total == 0.0 );
	    # Normalize the curve to unit area.
	    $errors     /= $total;
	    $histogram  /= $total;
	    # Dividing through by the sum of the histogram correlates the bins, leading to a modification of the covariance
	    # matrix. Account for this by constructing the Jacobian of the transformation from the unnormalized to normalized
	    # histogram and using this to propagate the original covariance matrix.
	    my $jacobian = pdl zeroes(nelem($histogram),nelem($histogram));
	    for(my $k=0;$k<nelem($histogram);++$k) {
		$jacobian->(($k),  : ) .= -$histogram->flat();
		$jacobian->(($k),($k)) += 1.0/$total;
	    }
	    my $covarianceCorrelated  = $jacobian x $covariance x transpose($jacobian);
	    $covariance              .= $covarianceCorrelated;
	} elsif ( $options{'normalized'} eq "x" ) {
	    # Normalize each row.
	    my $jacobian = pdl zeroes(nelem($histogram),nelem($histogram));
	    for(my $j=0;$j<$histogram->dim(1);++$j) {
		# Find the total weight.
		my $total;
		if ( exists($options{'normalizeBy'}) && $options{'normalizeBy'} eq "histogram" ) {
		    $total = sum($histogram->(:,($j)));   
		} else {
		    $total = $columnWeights->(($j));
		}
		$total = 1.0
		    if ( $total == 0.0 );
		# Normalize the curve to unit area.
		$errors   ->(:,($j)) /= $total;
		$histogram->(:,($j)) /= $total;
		# Dividing through by the sum of the histogram correlates the bins, leading to a modification of the covariance
		# matrix. Account for this by constructing the Jacobian of the transformation from the unnormalized to normalized
		# histogram and using this to propagate the original covariance matrix.
		for(my $i=0;$i<$histogram->dim(0);++$i) {
		    my $k = $i+$j*$xBinCount;
		    for(my $ii=0;$ii<$xBinCount;++$ii) {
			$jacobian->(($k),($ii+$j*$xBinCount)) += -$histogram->(($ii),($j))/$total;
		    }
		    $jacobian->(($k),($k)) += 1.0/$total;
		}
	    }
	    my $covarianceCorrelated  = $jacobian x $covariance x transpose($jacobian);
	    $covariance              .= $covarianceCorrelated;
	} elsif ( $options{'normalized'} eq "y" ) {
	    # Normalize each column.
	    my $jacobian = pdl zeroes(nelem($histogram),nelem($histogram));
	    for(my $i=0;$i<$histogram->dim(0);++$i) {
		# Find the total weight.
		my $total;
		if ( exists($options{'normalizeBy'}) && $options{'normalizeBy'} eq "histogram" ) {
		    $total = sum($histogram->(($i),:));   
		} else {
		    $total = $rowWeights->(($i));
		}
		$total = 1.0
		    if ( $total == 0.0 );
		# Normalize the curve to unit area.
		$errors   ->(($i),:) /= $total;
		$histogram->(($i),:) /= $total;
		# Dividing through by the sum of the histogram correlates the bins, leading to a modification of the covariance
		# matrix. Account for this by constructing the Jacobian of the transformation from the unnormalized to normalized
		# histogram and using this to propagate the original covariance matrix.
		for(my $j=0;$j<$histogram->dim(1);++$j) {
		    my $k = $i+$j*$xBinCount;
		    for(my $jj=0;$jj<$yBinCount;++$jj) {
			$jacobian->(($k),($i+$jj*$xBinCount)) .= -$histogram->(($i),($jj))/$total;
		    }
		    $jacobian->(($k),($k)) += 1.0/$total;
		}
	    }
	    my $covarianceCorrelated  = $jacobian x $covariance x transpose($jacobian);
	    $covariance              .= $covarianceCorrelated;
	} else {
	    die("Histogram2D: unrecognized 'normalized' option");
	}
    }
    if ( exists($options{'differential'}) ) {
	# Divide by the bin width to get a differential distribution.
	my $widths = pdl ones($xBinCount,$yBinCount);
	$widths .= $xBinWidth
	    if ( $options{'differential'} =~ m/x/ );
	if ( $options{'differential'} =~ m/y/ ) {
	    for(my $j=0;$j<$yBinCenters->dim(0);++$j) {
		$widths->(:,($j)) *= $yBinWidth->(($j));
	    }
	}
	$errors     /= $widths;
	$histogram  /= $widths;
	$covariance /= transpose($widths->flat()) x $widths->flat();
    }
    # Return the results.
    return ($histogram,$errors,$covariance);
}

1;
