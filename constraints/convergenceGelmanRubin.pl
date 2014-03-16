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
     outliersMaximum => 10
    );
&Options::Parse_Options(\@ARGV,\%arguments);

# Parse the constraint config file.
my $config = &Parameters::Parse_Config($configFile);

# Validate the config file.
die("convergenceGelmanRubin.pl: workDirectory must be specified in config file" ) unless ( exists($config->{'workDirectory' }) );
die("convergenceGelmanRubin.pl: nodes must be specified in config file"         ) unless ( exists($config->{'nodes'         }) );
die("convergenceGelmanRubin.pl: threadsPerNode must be specified in config file") unless ( exists($config->{'threadsPerNode'}) );

# Compute the number of parallel chains.
my $chainCount = $config->{'nodes'}*$config->{'threadsPerNode'};

# Construct mask of outlier chains.
my $outlierMask = pdl zeroes($chainCount);
if ( exists($arguments{'outliers'}) ) {
    foreach ( split(/\s*,\s*/,$arguments{'outliers'}) ) {
	$outlierMask->(($_)) .= 1;
    }
}

# Read the header line.
open(iHndl,$config->{'workDirectory'}."/mcmc/galacticusBIE.statelog");
my $header = <iHndl>;
close(iHndl);
chomp($header);
$header =~ s/^"//;
$header =~ s/"$//;
my @labels = split(/"\s*"/,$header);

# Determine the number of parameters.
my $parameterCount = scalar(@labels)-5;

# Read the chains.
my @columns        = 5..($parameterCount+4);
my @chainData      = rcols($config->{'workDirectory'}."/mcmc/galacticusBIE.statelog",@columns,{IGNORE => '/"/'});
my @likelihoodData = rcols($config->{'workDirectory'}."/mcmc/galacticusBIE.statelog",2       ,{IGNORE => '/"/'});

# Determine how many steps to burn.
my $burnCount   = nelem($chainData[0])*0;
my $chainLength = nelem($chainData[0])-$burnCount;

# Detect outlier chains.
my @currentState;
for(my $j=0;$j<$parameterCount;++$j) {
    $currentState[$j] = $chainData[$j]->(-$chainCount:-1;|);
}
while ( $outlierMask->sum() <= $arguments{'outliersMaximum'} ) {
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
	    if ($j == 10 && $currentState[$j]->($activeChains)->(($i)) > 0.0) {
		$deviationMaximum          = 10.0;
		$deviationMaximumChain     = $activeChains->(($i));
		$deviationMaximumParameter = $j;
	    }
	}
    }
    if ( defined($deviationMaximumChain) && defined($deviationMaximumParameter) ) {
	my $offset = $likelihoodData[0]->(($deviationMaximumChain))-max($likelihoodData[0]);
	print "Outlier:\n";
	print "         G = ".$grubbsCriticalValue."\n";
	print "         D = ".$deviationMaximum."\n";
	print "     chain = ".$deviationMaximumChain."\n";
	print " parameter = ".$deviationMaximumParameter."\n";
	print "    offset = ".$offset."\n";
	$outlierMask->(($deviationMaximumChain)) .= 1;
	my $outlierCount = $outlierMask->sum();
	print "  outliers = ".$outlierCount."\n";
    } else {
	last;
    }
}
my $outliers = which($outlierMask == 1);
die("convergenceGelmanRubin.pl: maximum number of outliers exceeded")
    if ( nelem($outliers) > $arguments{'outliersMaximum'} );

# Extract the unburned section of each chain and discard outlier chains.
my @chains;
my $ii = -1;
for(my $i=0;$i<$chainCount;++$i) {
    if ( $outlierMask->(($i)) == 0 ) {
	++$ii;
	for(my $j=0;$j<$parameterCount;++$j) {
	    $chains[$ii][$j] = $chainData[$j]->($burnCount+$i:-1:$chainCount;|);
	}
    }
}
my $stepCount = nelem($chains[0][0]); # Determine number of available steps.
$chainCount = $ii+1;                  # Reset number of chains to total minus outliers.

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
for(my $j=0;$j<$parameterCount;++$j) {
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
    # Store results.
    push(
	@convergence,
	{
	    parameter        => $j                 ,
	    Rhatc            => $Rhatc             ,
	    mean             => $interChainMean    ,
	    between          => $interChainVariance,
	    within           => $W                 ,
	    degreesOfFreedom => $d                 ,
	    label            => $labels[$j+5]
	}
	);
}
my @sortedConvergence = sort {$a->{'Rhatc'}->sclr() <=> $b->{'Rhatc'}->sclr()} @convergence;
foreach ( @sortedConvergence ) {
    print $_->{'parameter'}."\t".$_->{'mean'}."\t".$_->{'between'}."\t".$_->{'within'}."\t".$_->{'Rhatc'}."\t".$_->{'degreesOfFreedom'}."\t".$_->{'label'}."\n";
}
print "\n";
print "Minimum/Maximum Rhat = ".$sortedConvergence[0]->{'Rhatc'}."\t".$sortedConvergence[-1]->{'Rhatc'}."\n\n";
print "Outliers: ".join(",",$outliers->list())."\n";
print "Chain length: ".$stepCount."\n";

exit;
