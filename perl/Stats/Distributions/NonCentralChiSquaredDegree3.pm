package Stats::Distributions::NonCentralChiSquaredDegree3;
use strict;
use warnings;
use Cwd;
use lib exists($ENV{'GALACTICUS_ROOT_V094'}) ? $ENV{'GALACTICUS_ROOT_V094'}.'/perl' : cwd().'/perl';
use PDL;
use PDL::NiceSlice;
use PDL::Constants qw(PI);

# Functions related to the non-central chi-squared distribution for 3 degrees of freedom.
# Andrew Benson (22-March-2017)

sub variance {
    # Returns the variance of a k=3 degrees of freedom non-central χ² distribution with given non-centrality parameter.
    my $nonCentrality = shift();
    my $variance      = 3.0+$nonCentrality-&mean($nonCentrality)**2;
    return $variance;
}

sub mean {
    # Returns the mean of √x with a k=3 degrees of freedom non-central χ² distribution with given non-centrality parameter. Note
    # that we want the mean of √x here as we are interested in the mean of λ while it is λ² that is distribution as a non-central
    # χ² distribution. Result was computed using Mathematica. For large non-centrality parameter the answer converges to √x, while
    # it converges to 2√(2/π)≅1.596 as x → 0.
    my $nonCentrality = shift();
    my $mean          = sqrt(2.0/PI)*exp(-0.5*$nonCentrality)+(1.0+$nonCentrality)*erf(sqrt(0.5*$nonCentrality))/sqrt($nonCentrality);
    return $mean;
}

sub median {
    # Compute the median value of a k=3 degrees of freedom non-central χ² distribution. This must be done numerically.
    my $nonCentrality = shift();
    my $mean          = &mean($nonCentrality);
    my $medians       = pdl zeroes(nelem($nonCentrality));
    for(my $i=0;$i<nelem($nonCentrality);++$i) {
	print "Finding median ".$i." of ".nelem($nonCentrality)."\n";
	# For large non-centrality parameters, return the limiting behavior.
	if ( $nonCentrality->(($i)) > 100.0 ) {
	    $medians->(($i)) .= sqrt($nonCentrality->(($i)));
	} else {
	    # For non-large non-centrality parameters, find the solution numerically.
	    # First guess for the median is the mean.
	    my $lambda = $nonCentrality->(($i));
	    my $median = 3.0*$lambda;
	    # Iteratively seek an improved solution. We use a simple approach - evaluate the CMF, and adjust the median by a small
	    # factor up or down until we get close to CMF=1/2.
	    my $factor = pdl 1.0;
	    my $cmf    = pdl 0.0;
	    while ( abs($cmf-0.5) > 0.001 ) {
		$median *= $factor;	    
		$cmf    .= exp(-0.5*(sqrt($median)+sqrt($lambda))**2)*($median*$lambda)**0.75*(2.0-2.0*exp(2.0*sqrt($median*$lambda))+exp(0.5*(sqrt($median)+sqrt($lambda))**2)*sqrt(2.0*PI*$lambda)*(erf((sqrt($median)-sqrt($lambda))/sqrt(2.0))+erf((sqrt($median)+sqrt($lambda))/sqrt(2.0))))/sqrt(2.0*PI)/$lambda**2/($median/$lambda)**0.75/2.0;
		if ( $cmf > 0.5 ) {
		    $factor .= 0.99;
		} else {
		    $factor .= 1.01;
		}
	    }
	    # Convert from x to √x.
	    $medians->(($i)) .= sqrt($median);
	}
    }
    return $medians;
}

1;
