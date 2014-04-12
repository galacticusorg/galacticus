#!/usr/bin/env perl
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
use PDL::IO::Misc;
use PDL::Fit::Polynomial;
require Stats::Percentiles;
require GnuPlot::PrettyPlots;
require GnuPlot::LaTeX;

# Construct the median line-width vs. HI mass relation for the ALFALFA survey.
# Andrew Benson (9-July-2013)

# Define working directory.
my $workDirectory = "constraints/dataAnalysis/hiMassFunction_ALFALFA_z0.00/";

# Retrieve the ALFALFA 40% dataset if we don't already have it.
system("wget http://egg.astro.cornell.edu/alfalfa/data/a40files/a40.datafile1.csv -O ".$workDirectory."alfalfaSurveyData.csv")
    unless ( -e $workDirectory."alfalfaSurveyData.csv" );

# Load the data.
(my $W50, my $errorW50, my $logarithmicMassHI) = rcols($workDirectory."alfalfaSurveyData.csv",7,8,14,{ IGNORE => '/^AGCNr/', COLSEP => ','});

# Find the binned median line width vs. logarithmic mass.
my $usableGalaxies          = 
    which(
	($logarithmicMassHI > 0.0)
	&
	($errorW50          > 0.0)
    );
my $sortedMassIndex         = $logarithmicMassHI->($usableGalaxies)->qsorti();
my $logarithmicMassMinimum  = $logarithmicMassHI->($usableGalaxies)->($sortedMassIndex)->(( 0));
my $logarithmicMassMaximum  = $logarithmicMassHI->($usableGalaxies)->($sortedMassIndex)->((-1));
my $logarithmicMassCount    = nelem($usableGalaxies);
my $logarithmicMassBinCount = long($logarithmicMassCount/500);
my $logarithmicMassBins     = pdl sequence($logarithmicMassBinCount->sclr())*($logarithmicMassMaximum-$logarithmicMassMinimum)/($logarithmicMassBinCount-1)+$logarithmicMassMinimum;
my $percentiles             = pdl ( 50.0 );
my $weight                  = 1.0/$errorW50**2;
my $median                  = &Percentiles::BinnedPercentiles($logarithmicMassBins,$logarithmicMassHI->($usableGalaxies),log10($W50->($usableGalaxies)),$weight->($usableGalaxies),$percentiles);

# Fit a polynomial to the relation.
my $nonZeroBins = which($median->(:,(0)) > 0.0);
(my $fit, my $coeffs) = fitpoly1d($logarithmicMassBins->index($nonZeroBins),$median->(:,(0))->index($nonZeroBins),2);

# Generate a fit to the data.
my $fitW50 = pdl zeroes(nelem($logarithmicMassBins));
my $closing     = "";
for(my $i=0;$i<nelem($coeffs);++$i) {
    $fitW50 += $coeffs->(($i))*($logarithmicMassBins**$i);
    print $coeffs->(($i))."d0";
    if ( $i == nelem($coeffs)-1 ) {
	print $closing."\n";
    } else {
	print "+logarithmicMass*(";
	$closing .= ")";
    }
}

# Create a plot.
my $plot;
my $gnuPlot;
my $plotFile = $workDirectory."lineWidthMassRelation.pdf";
(my $plotFileEPS = $plotFile) =~ s/\.pdf$/.eps/;
open($gnuPlot,"|gnuplot 1>/dev/null 2>&1");
print $gnuPlot "set terminal epslatex color colortext lw 2 solid 7\n";
print $gnuPlot "set output '".$plotFileEPS."'\n";
print $gnuPlot "set title 'Line width vs. HI gas mass for ALFALFA 40\\\% (\$\\alpha.40\$) sample'\n";
print $gnuPlot "set xlabel 'HI gas mass; \$M_{\\rm HI} [M_\\odot]\$'\n";
print $gnuPlot "set ylabel 'Line width; \$W_{50}\$ [km/s]'\n";
print $gnuPlot "set lmargin screen 0.15\n";
print $gnuPlot "set rmargin screen 0.95\n";
print $gnuPlot "set bmargin screen 0.15\n";
print $gnuPlot "set tmargin screen 0.95\n";
print $gnuPlot "set key spacing 1.2\n";
print $gnuPlot "set key at screen 0.275,0.76\n";
print $gnuPlot "set key left\n";
print $gnuPlot "set key bottom\n";
print $gnuPlot "set logscale xy\n";
print $gnuPlot "set mxtics 10\n";
print $gnuPlot "set mytics 10\n";
print $gnuPlot "set format x '\$10^{\%L}\$'\n";
print $gnuPlot "set format y '\$10^{\%L}\$'\n";
print $gnuPlot "set xrange [1.0e8:1.0e11]\n";
print $gnuPlot "set yrange [10.0:600.0]\n";
print $gnuPlot "set pointsize 2.0\n";
my $r             = pdl random(nelem($usableGalaxies));
my $plotSelection = which($r < 0.01);
&PrettyPlots::Prepare_Dataset(
    \$plot,
    10.0**$logarithmicMassHI->($usableGalaxies)->($plotSelection),
    $W50                    ->($usableGalaxies)->($plotSelection),
    errorUp    => $errorW50 ->($usableGalaxies)->($plotSelection),
    errorDown  => $errorW50 ->($usableGalaxies)->($plotSelection),
    style      => "point",
    weight     => [2,1],
    symbol     => [6,7],
    pointSize  => 0.1,
    color      => $PrettyPlots::colorPairs{'redYellow'},
    title      => '\$\\\\alpha.40\$'
    );
&PrettyPlots::Prepare_Dataset(
    \$plot,
    10.0**$logarithmicMassBins,
    10.0**$median->(:,(0)),
    style      => "point",
    weight     => [5,3],
    symbol     => [6,7],
    pointSize  => 2.0,
    color      => $PrettyPlots::colorPairs{'peachPuff'},
    title      => 'Median'
    );
&PrettyPlots::Prepare_Dataset(
    \$plot,
    10.0**$logarithmicMassBins,
    10.0**$fitW50,
    style      => "line",
    weight     => [5,3],
    color      => $PrettyPlots::colorPairs{'mediumSeaGreen'},
    title      => 'Fit'
    );
&PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
close($gnuPlot);
&LaTeX::GnuPlot2PDF($plotFileEPS);

exit;
