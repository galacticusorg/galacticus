# Contains a Perl module which provides various systematics models for
# use in model discrepancy modeling.

package DiscrepancySystematics;
use strict;
use warnings;
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V093"}) ) {
    $galacticusPath  = $ENV{"GALACTICUS_ROOT_V093"};
    $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
    $galacticusPath  = "./";
}
unshift(@INC,$galacticusPath."perl"); 
use PDL;
use PDL::NiceSlice;
use Data::Dumper;

our %models = 
    (
     MassShift => \&MassFunctionMassShift
    );

sub MassFunctionMassShift {
    # Implements a model discrepancy systematic suitable for mass
    # functions in which a simple model is applied to the masses
    # utilizing an abundance matching approach.
    my %arguments   = %{shift()};
    my $constraint  = shift;
    my $X           = shift;
    my $discrepantY = shift;
    my $discrepantC = shift;
    my $trueY       = shift;
    my $trueC       = shift;
    # Extract model parameters from arguments.
    my $zeroPoint;
    $zeroPoint = $arguments{'systematicMassShiftZeroPoint'};
    die("MassFunctionMassShift: systematicMassShiftZeroPoint argument must be specified")
	unless ( defined($zeroPoint) );
    if ( $zeroPoint eq "constraint" ) {
	die("MassFunctionMassShift: constraint must specify a mass systematic zero point")
	    unless ( exists($constraint->{'massSystematicZeroPoint'}) );
	$zeroPoint = $constraint->{'massSystematicZeroPoint'};
    }
    my $systematicMassShiftOrder;
    $systematicMassShiftOrder = $arguments{'systematicMassShiftOrder'};
    die("MassFunctionMassShift: systematicMassShiftOrder argument must be specified")
	unless ( defined($systematicMassShiftOrder) );
    if ( $systematicMassShiftOrder eq "constraint" ) {
	die("MassFunctionMassShift: constraint must specify a mass systematic order")
	    unless ( exists($constraint->{'massSystematicOrder'}) );
	$zeroPoint = $constraint->{'massSystematicOrder'};
    }
    # First find cumulative mass functions.
    my $discrepantCumulative         = $discrepantY               ->(-1:0)->cumusumover()->(-1:0);
    my $discrepantCumulativeVariance = $discrepantC->diagonal(0,1)->(-1:0)->cumusumover()->(-1:0);
    my $trueCumulative               = $trueY                     ->(-1:0)->cumusumover()->(-1:0);
    my $trueCumulativeVariance       = $trueC      ->diagonal(0,1)->(-1:0)->cumusumover()->(-1:0);
    # Find the errors in the x-values.
    my $deltaX                       = $X->((1))-$X->((0));
    my $discrepantGradient           = ($discrepantY/$deltaX);
    my $trueGradient                 = ($trueY      /$deltaX);
    my $discrepantXVariance          = $discrepantCumulativeVariance/$discrepantGradient**2;
    my $trueXVariance                = $trueCumulativeVariance      /$trueGradient      **2;
    # Search for the best-fit parameters of the systematic model.
    my $testStatisticMinimum         = pdl 1.0e30;
    my $coefficientBest              = pdl zeroes($systematicMassShiftOrder+1);
    my $coefficientRange             = pdl zeroes($systematicMassShiftOrder+1);
    my $coefficientStep              = pdl zeroes($systematicMassShiftOrder+1);
    for(my $i=0;$i<=$systematicMassShiftOrder;++$i) {
	die("MassFunctionMassShift: systematicMassShiftCoefficientRange".$i." argument must be present")
	    unless ( exists($arguments{'systematicMassShiftCoefficientRange'.$i}) );
	die("MassFunctionMassShift: systematicMassShiftCoefficientStep".$i." argument must be present")
	    unless ( exists($arguments{'systematicMassShiftCoefficientStep'.$i}) );
	$coefficientRange->(($i)) .= $arguments{'systematicMassShiftCoefficientRange'.$i};
	$coefficientStep ->(($i)) .= $arguments{'systematicMassShiftCoefficientStep' .$i};
    }
    my $coefficient = -$coefficientRange->copy();
    my $done = 0;
    while ( $done == 0 ) {
	# Shift the x-axis in the discrepant model using systematic model.
	my $discrepantXShifted                          = $X->copy();
	for(my $i=0;$i<=$systematicMassShiftOrder;++$i) {
	    $discrepantXShifted += $coefficient->(($i))*($X-$zeroPoint)**$i;
	}
	# Interpolate these shifted x-values by matching abundances between discrepant and true models.
	(my $discrepantXShiftedInterpolated, my $error) = interpolate($trueCumulative,$discrepantCumulative,$discrepantXShifted);
	# Compute a test statistic.
	my $offset                                      = $discrepantXShiftedInterpolated-$X;
	my $variance                                    = $discrepantXVariance+$trueXVariance;
	my $testStatistic                               = sum($offset**2/$variance);
	# Minimize the test statistic.
	if ( $testStatistic < $testStatisticMinimum ) {
	    $testStatisticMinimum .= $testStatistic;
	    $coefficientBest      .= $coefficient->copy();
	}
	# Shift to next set of coefficients.
	my $i = 0;
	while ( $i >= 0 ) {
	    $coefficient->(($i)) += $coefficientStep->(($i));
	    if ( $coefficient->(($i)) > $coefficientRange->(($i))+0.5*$coefficientStep->(($i)) ) {
		$coefficient->(($i)) .= -$coefficientRange->(($i));
		++$i;
		if ( $i > $systematicMassShiftOrder ) {
		    $i    = -1;
		    $done =  1;
		}
	    } else {
		$i = -1;
	    }
	}
    }
    # Report on best-fit systematic model.
    print "Best fit systematic model coefficients: ".$coefficientBest."\n";
    # Shift the discrepant model x-values by the best fit systematic model.
    my $discrepantXShiftedBest  = $X->copy();
    for(my $i=0;$i<=$systematicMassShiftOrder;++$i) {
	$discrepantXShiftedBest += $coefficientBest->(($i))*($X-$zeroPoint)**$i;
    }
    # Interpolate the y-values from the shifted discrepant model onto the original x-values.
    my $discrepantYShiftedBest  = interpol($X,$discrepantXShiftedBest,$discrepantY);
    # Compute multiplicative correction.
    my $discrepantYMultiplier   = $discrepantYShiftedBest/$discrepantY;
    # Replace the original y-values of the discrepant model.
    $discrepantY               .= $discrepantYShiftedBest;
    # Scale the covariance appropriately.
    $discrepantC               *= outer($discrepantYMultiplier,$discrepantYMultiplier);
    # Return model parameters.
    my %results;
    for(my $i=0;$i<=$systematicMassShiftOrder;++$i) {
	my $label = "systematicMassShiftCoefficient".$i;
	$results{$label} = $coefficientBest->(($i))->sclr();
    }
    return %results;
}

1;
