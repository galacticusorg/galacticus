#!/usr/bin/env perl
use strict;
use warnings;
use Cwd;
use lib $ENV{'GALACTICUS_EXEC_PATH'         }."/perl";
use lib $ENV{'GALACTICUS_ANALYSIS_PERL_PATH'}."/perl";
use PDL;
use PDL::NiceSlice;
use Galacticus::Constraints::Parameters;
use Galacticus::Options;
use GnuPlot::PrettyPlots;
use GnuPlot::LaTeX;

# Compute acceptance rates from an MCMC simulation.
# Andrew Benson (04-November-2019)

# Get command line arguments.
die("Usage: acceptanceRate.pl <parameterFile> [options...]")
    unless ( scalar(@ARGV) >= 1 );
my $parameterFile = $ARGV[0];
# Create a hash of named options.
my %options =
    (
     stateSampleCount => 30
    );
&Galacticus::Options::Parse_Options(\@ARGV,\%options);

# Parse the constraint config file.
my $config = &Galacticus::Constraints::Parameters::parseConfig($parameterFile);

# Compute acceptance rates.
my $chainCount      = &Galacticus::Constraints::Parameters::chainCount($config,\%options);
my $stepCount       = &Galacticus::Constraints::Parameters::stepCount ($config,\%options);
die('acceptanceRate.pl: too few steps to compute acceptance rates')
    if ( $stepCount < $options{'stateSampleCount'} );
my $acceptanceRates = pdl zeros($chainCount,$stepCount-$options{'stateSampleCount'});
for(my $chain=0;$chain<$chainCount;++$chain) {
    my $parameterStates = &Galacticus::Constraints::Parameters::parameterMatrix($config,$chain,\%options);
    for(my $i=0;$i<$parameterStates->dim(1)-$options{'stateSampleCount'};++$i) {
	my $accepted = 0;
	for(my $j=0;$j<=$options{'stateSampleCount'};++$j) {
	    my $k = $i+$j;
	    $accepted++
		if ( $k > 0 &&  any($parameterStates->(:,($k-1)) != $parameterStates->(:,($k))) );
	}
	$acceptanceRates->(($chain),($i)) .= $accepted/$options{'stateSampleCount'};
    }	
}

# Compute mean acceptance rate across chains.
my $acceptanceRateMean = $acceptanceRates->average();

# Compute mean acceptance rates within chains.
my $acceptanceRateChains = $acceptanceRates     ->xchg(0,1)->average    ();
my $chainBest            = $acceptanceRateChains           ->maximum_ind();
my $chainWorst           = $acceptanceRateChains           ->minimum_ind();

# Make a plot of acceptance rates.
my $plot;
my $gnuPlot;
my $plotFileTeX = "acceptanceRates.tex";
open($gnuPlot,"|gnuplot");
print $gnuPlot "set terminal cairolatex pdf standalone color lw 2 size 4in,4in\n";
print $gnuPlot "set output '".$plotFileTeX."'\n";
print $gnuPlot "set xlabel 'Step\n";
print $gnuPlot "set ylabel 'Acceptance rate'\n";
print $gnuPlot "set lmargin screen 0.15\n";
print $gnuPlot "set rmargin screen 0.95\n";
print $gnuPlot "set bmargin screen 0.15\n";
print $gnuPlot "set tmargin screen 0.95\n";
print $gnuPlot "set xrange [".($options{'stateSampleCount'}-1).":".$stepCount."]\n";
print $gnuPlot "set yrange [".(-$acceptanceRates->flat()->maximum()*0.01).":".($acceptanceRates->flat()->maximum()*1.01)."]\n";
print $gnuPlot "set pointsize 1.0\n";
my $step                     = pdl sequence($stepCount-$options{'stateSampleCount'})+$options{'stateSampleCount'};
my $rankValue                = pdl zeros($chainCount);
$rankValue->(($chainBest )) .= 1;
$rankValue->(($chainWorst)) .= 1;
my $rank                     = $rankValue->qsorti();
for(my $chain=0;$chain<$chainCount;++$chain) {
    my $color = "peachPuff";
    $color = "mediumSeaGreen"
	if ( $rank->(($chain)) == $chainBest  );
    $color = "indianRed"
	if ( $rank->(($chain)) == $chainWorst );
    &GnuPlot::PrettyPlots::Prepare_Dataset(
	 \$plot                                                   ,
	 $step                                                    ,
	 $acceptanceRates->(($rank->(($chain))),:)                ,
	 style        => "line"                                   ,
	 weight       => [1,1]                                    ,
	 color        => $GnuPlot::PrettyPlots::colorPairs{$color}
	);
}
&GnuPlot::PrettyPlots::Prepare_Dataset(
    \$plot                                                        ,
    $step                                                         ,
    $acceptanceRateMean                                           ,
    style        => "line"                                        ,
    weight       => [1,1]                                         ,
    color        => $GnuPlot::PrettyPlots::colorPairs{'redYellow'}
    );
&GnuPlot::PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
close($gnuPlot);
&GnuPlot::LaTeX::GnuPlot2PDF($plotFileTeX);

# Report.
print "Chain with highest acceptance rate is #".$chainBest ." with rate of ".$acceptanceRateChains->(($chainBest ))."\n";
print "Chain with lowest  acceptance rate is #".$chainWorst." with rate of ".$acceptanceRateChains->(($chainWorst))."\n";

exit 0;
