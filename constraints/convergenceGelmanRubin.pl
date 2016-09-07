#!/usr/bin/perl
use strict;
use warnings;
use lib './perl';
use PDL;
use PDL::NiceSlice;
use PDL::IO::Misc;
use PDL::Stats;
use Data::Dumper;
require Galacticus::Options;
require Galacticus::Constraints::Parameters;

# Compute the Gelman & Rubin Rhat convergence statistic using a statelog
# from the Bayesian Inference Engine. Based on Brooks & Gelman (1998;
# Journal of Computational and Graphical Statistics; 7; 4; 434-455) and
# Gelman & Rubin (1992; Statistical Science; 7; 4; 457-472).

# Get command line arguments.
die("Usage: convergenceGelmanRubin.pl <configFile> [options...]")
    unless ( scalar(@ARGV) >= 1 );
my $configFile = $ARGV[0];
# Create a hash of named arguments.
my %arguments =
    (
     outliersMaximum  => 10,
     useUnconverged   => "yes",
     parametersMapped => "yes",
     randomSample     => "no" ,
     includePrevious  => "yes"
    );
&Options::Parse_Options(\@ARGV,\%arguments);

# Parse the constraint config file.
my $config = &Parameters::Parse_Config($configFile);

# Validate the config file.
die("convergenceGelmanRubin.pl: workDirectory must be specified in config file")
    unless ( exists($config->{'likelihood'}->{'workDirectory' }) );

# Compute the number of parallel chains.
my $chainCount = &Parameters::Chains_Count($config,\%arguments);
print "Found ".$chainCount." chains\n";

# Construct mask of outlier chains.
my $outlierMask = pdl zeroes($chainCount);
if ( exists($arguments{'outliers'}) ) {
    foreach ( split(/\s*,\s*/,$arguments{'outliers'}) ) {
	$outlierMask->(($_)) .= 1;
    }
}

# Determine the number of parameters.
my $parameterCount = &Parameters::Parameters_Count($config,\%arguments);
print "Found ".$parameterCount." parameters\n";

# Get the chain state matrix.
my @chainData;
for(my $i=0;$i<$chainCount;++$i) {
    $arguments{'selectChain'} = $i;
    push(@chainData,&Parameters::Sample_Matrix($config,\%arguments));
}

# Detect outlier chains.
my @currentState;
for(my $j=0;$j<$parameterCount;++$j) {
    $currentState[$j] = pdl [];
    for(my $i=0;$i<$chainCount;++$i) {
	$currentState[$j] = $currentState[$j]->append($chainData[$i]->(($j),(-1)));
    }
}
while ( $outlierMask->sum() < $arguments{'outliersMaximum'} ) {
    my $activeChains     = which($outlierMask == 0);
    my $activeChainCount = nelem($activeChains);
    my @currentMean;
    my @currentMeanSquared;
    my @currentVariance;
    my @deviation;
    for(my $j=0;$j<$parameterCount;++$j) {
	$currentMean       [$j] = average($currentState[$j]->($activeChains)   );
	$currentMeanSquared[$j] = average($currentState[$j]->($activeChains)**2);
	$currentVariance   [$j] = ($currentMeanSquared[$j]-$currentMean[$j]**2)*$activeChainCount/($activeChainCount-1);
	$deviation         [$j] = abs($currentState[$j]->($activeChains)-$currentMean[$j])/sqrt($currentVariance[$j]);
    }
    my $outlierSignificance     = pdl 0.05;
    my $tStatisticCriticalValue = &gsl_cdf_tdist_Pinv($outlierSignificance/2.0/$activeChainCount,$activeChainCount-2);
    my $grubbsCriticalValue     = ($activeChainCount-1)/sqrt($activeChainCount)
	*sqrt($tStatisticCriticalValue**2/($activeChainCount-2.0+$tStatisticCriticalValue**2));
    my $deviationMaximum        = 0.0;
    my $deviationMaximumChain;
    my $deviationMaximumParameter;
    for(my $i=0;$i<$activeChainCount;++$i) {
	for(my $j=0;$j<$parameterCount;++$j) {
	    if ($deviation[$j]->(($i)) > $grubbsCriticalValue && $deviation[$j]->(($i)) > $deviationMaximum) {
		$deviationMaximum          = $deviation[$j]->(($i));
		$deviationMaximumChain     = $activeChains->(($i));
		$deviationMaximumParameter = $j;
	    }
	}
    }
    if ( defined($deviationMaximumChain) && defined($deviationMaximumParameter) ) {
	print "Outlier:\n";
	print "         G = ".$grubbsCriticalValue."\n";
	print "         D = ".$deviationMaximum."\n";
	print "     chain = ".$deviationMaximumChain."\n";
	print " parameter = ".$deviationMaximumParameter."\n";
	$outlierMask->(($deviationMaximumChain)) .= 1;
	my $outlierCount = $outlierMask->sum();
	print "  outliers = ".$outlierCount."\n";
    } else {
	last;
    }
}
my $outliers = which($outlierMask == 1);
print "convergenceGelmanRubin.pl: maximum number of outliers exceeded\n"
    if ( nelem($outliers) > $arguments{'outliersMaximum'} );

# Discard outlier chains.
my @chains;
my $ii = -1;
for(my $i=0;$i<$chainCount;++$i) {
    if ( $outlierMask->(($i)) == 0 ) {
	++$ii;
	for(my $j=0;$j<$parameterCount;++$j) {
	    $chains[$ii][$j] = $chainData[$i]->(($j),:);
	}
    }
}
my $stepCount = nelem($chains[0][0]); # Determine number of available steps.
$chainCount = $ii+1;                  # Reset number of chains to total minus outliers.
print "Found ".$stepCount." states in ".$chainCount." non-outlier chains\n";

# Find the mean within each chain.
my @mean;
for(my $i=0;$i<$chainCount;++$i) {
    for(my $j=0;$j<$parameterCount;++$j) {
	$mean[$i][$j] = average($chains[$i][$j]);
    }
}

# Find the variance within each chain.
my @variance;
for(my $i=0;$i<$chainCount;++$i) {
    for(my $j=0;$j<$parameterCount;++$j) {
	$variance[$i][$j] = average(($chains[$i][$j]-$mean[$i][$j])**2)*$stepCount/($stepCount-1);
    }
}

# Iterate over parameters.
my @convergence;
my $iParameter = -1;
for(my $j=0;$j<$parameterCount;++$j) {
    # Find the next active parameter in our list.
    ++$iParameter;
    while ( ! exists($config->{'parameters'}->{'parameter'}->[$iParameter]->{'prior'}) ) {
	++$iParameter;
    }
    # Find the mean over all chains for this parameter.
    my $interChainMean  = pdl 0.0;
    my $interChainMean2 = pdl 0.0;
    for(my $i=0;$i<$chainCount;++$i) {
	$interChainMean  += $mean[$i][$j];
	$interChainMean2 += $mean[$i][$j]**2;
    }
    $interChainMean  /= $chainCount;
    $interChainMean2 /= $chainCount;
    # Find the variance between chain means for this parameter.
    my $interChainVariance = pdl 0.0;
    for(my $i=0;$i<$chainCount;++$i) {
	$interChainVariance += ($mean[$i][$j]-$interChainMean)**2;
    }
    $interChainVariance /= $chainCount-1;
    # Compute B from Brooks & Gelman (section 1.2).
    my $B = $stepCount*$interChainVariance;
    # Compute W from Brooks & Gelman (section 1.2).
    my $W = pdl 0.0;
    for(my $i=0;$i<$chainCount;++$i) {
	$W += $variance[$i][$j];
    }
    $W /= $chainCount;  
    # Compute variance of chain variances.
    my $varW = pdl 0.0;
    for(my $i=0;$i<$chainCount;++$i) {
	$varW += ($variance[$i][$j]-$W)**2;
    }
    $varW /= $chainCount-1;
    # Find the covariance of chain variances and means.
    my $covWx  = pdl 0.0;
    my $covWx2 = pdl 0.0;
    for(my $i=0;$i<$chainCount;++$i) {
	$covWx  += ($variance[$i][$j]-$W)*($mean[$i][$j]   -$interChainMean );
	$covWx2 += ($variance[$i][$j]-$W)*($mean[$i][$j]**2-$interChainMean2);
    }
    $covWx  /= $chainCount-2; # We've used up two degrees of freedom in estimating the mean and variance.
    $covWx2 /= $chainCount-2;
    # Estimate variance in Vhat from Gelman & Rubin (eqn. 4).
    my $varVhat = 
	(($stepCount-1)/$stepCount)**2*$varW/$chainCount
	+(($chainCount+1)/$chainCount/$stepCount)**2*(2/($chainCount-1))*$B**2
	+2*($chainCount+1)*($stepCount-1)/$chainCount/$stepCount**2
	*($stepCount/$chainCount)*($covWx2-2*$interChainMean*$covWx);
    my $Vhat = ($stepCount-1)/$stepCount*$W+$B/$stepCount+$B/$stepCount/$chainCount;
    # Estimate degrees of freedom (e.g. section 1.3 of Brooks & Gelman).
    my $d = 2*$Vhat**2/$varVhat;
    # Compute R statistic (Brooks & Gelman, eqn. 1.1).
    my $Rhat = sqrt((($chainCount+1)/$chainCount)*((($stepCount-1)/$stepCount)*$W+$B/$stepCount)/$W-($stepCount-1)/$stepCount/$chainCount);
    # Compute corrected statistic (Brooks & Gelman, eqn. at end of section 1.3).
    my $Rhatc = ($d+3)*$Rhat/($d+1);
    # For comparison we also compute the interval statistic proposed by Brooks & Gelman (1998, Section 4.3,
    # http://www.statslab.cam.ac.uk/Reports/1996/1996-4.ps.gz) since it makes no assumptions of normality in the parameter
    # distributions.
    my $alpha           = 0.15;
    my $x               = pdl [ $alpha/2.0, 1.0-$alpha/2.0 ];
    my $intervalLengths = pdl zeroes($chainCount);
    my $mixedChains     = pdl [];
    for(my $i=0;$i<$chainCount;++$i) {
	# Accumulate chain to sample of all chains.
	$mixedChains                       = $mixedChains->append($chains[$i][$j]);
	# Get index into an ordered list.
	my $rank                           = $chains[$i][$j]->qsorti();
	# Identify interval.
	my $uniform                        = pdl (sequence(nelem($rank))+1)/nelem($rank);
	(my $interval, my $intervalError)  = interpolate($x,$uniform,$chains[$i][$j]->($rank));
	$intervalLengths->(($i))          .= $interval->((1))-$interval->((0));
    }
    # Get index into an ordered list.
    my $rank                          = $mixedChains->qsorti();
    # Identify interval.
    my $uniform                       = pdl (sequence(nelem($rank))+1)/nelem($rank);
    (my $interval, my $intervalError) = interpolate($x,$uniform,$mixedChains->($rank));
    my $mixedIntervalLength           = $interval->((1))-$interval->((0));
    my $Rinterval                     = $mixedIntervalLength/$intervalLengths->average();
    # Store results.
    push(
	@convergence,
	{
	    parameter        => $j                 ,
	    Rhatc            => $Rhatc             ,
	    Rinterval        => $Rinterval         ,
	    mean             => $interChainMean    ,
	    between          => $interChainVariance,
	    within           => $W                 ,
	    degreesOfFreedom => $d                 ,
	    label            => $config->{'parameters'}->{'parameter'}->[$iParameter]->{'name'}
	}
	);
}
my @sortedConvergence = sort {$a->{'Rhatc'    }->sclr() <=> $b->{'Rhatc'    }->sclr()} @convergence;
my @sortedInterval    = sort {$a->{'Rinterval'}->sclr() <=> $b->{'Rinterval'}->sclr()} @convergence;
foreach ( @sortedConvergence ) {
    print $_->{'parameter'}."\t".$_->{'mean'}."\t".$_->{'between'}."\t".$_->{'within'}."\t".$_->{'Rhatc'}."\t".$_->{'Rinterval'}."\t".$_->{'degreesOfFreedom'}."\t".$_->{'label'}."\n";
}
print "\n";
print "Minimum/Maximum Rhat = ".$sortedConvergence[0]->{'Rhatc'    }."\t".$sortedConvergence[-1]->{'Rhatc'    }."\n";
print "Minimum/Maximum Rint = ".$sortedInterval   [0]->{'Rinterval'}."\t".$sortedInterval   [-1]->{'Rinterval'}."\n\n";
print "Outliers: ".join(",",$outliers->list())."\n";
print "Chain length: ".$stepCount."\n";

exit;
