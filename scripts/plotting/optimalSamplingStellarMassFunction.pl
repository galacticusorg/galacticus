#!/usr/bin/env perl
use strict;
use warnings;
use lib "./perl";
use PDL;
use XML::Simple;
use GnuPlot::PrettyPlots;
use GnuPlot::LaTeX;
use Data::Dumper;
use PDL::IO::HDF5;

# Make a plot of the tree sampling function.
# Andrew Benson (20-November-2010)

# Read arguments.
die("Usage: optimalSamplingStellarMassFunction.pl <parameterFile> <plotFile>")
    unless ( scalar(@ARGV) == 2 );
my $parameterFile = $ARGV[0];
my $plotFile      = $ARGV[1];

# Compile the sampling code.
system("make optimal_sampling.stellar_mass_function.exe");

# Run the sampling code.
system("optimal_sampling.stellar_mass_function.exe ".$parameterFile);

# Parse the parameter file.
my $xml = new XML::Simple;
my $parameters = $xml->XMLin($parameterFile);
my $outputFile = $parameters->{'parameter'}->{'optimalSamplingDensityOutputFileName'}->{'value'};

# Read data from the tree weighting file.
my $samplingData = new PDL::IO::HDF5($outputFile);
my $haloMass     = $samplingData->dataset('haloMass'        )->get();
my $gamma        = $samplingData->dataset('samplingDensity' )->get();
my $haloMF       = $samplingData->dataset('haloMassFunction')->get();
my $timeTree     = $samplingData->dataset('treeTiming'      )->get();

# Normalize to unit time.
my $deltaLogM  = log($haloMass->index(1)/$haloMass->index(0));
my $gammaTime  = sum($gamma *$timeTree*$deltaLogM);
my $haloMFTime = sum($haloMF*$timeTree*$deltaLogM);
$gamma        /= $gammaTime;
$haloMF       /= $haloMFTime;

# Declare variables for GnuPlot;
my ($gnuPlot, $plot);

# Open a pipe to GnuPlot.
(my $plotFileEPS = $plotFile) =~ s/\.pdf$/\.eps/;
open($gnuPlot,"|gnuplot");
print $gnuPlot "set terminal epslatex color colortext lw 2 solid 7\n";
print $gnuPlot "set output '".$plotFileEPS."'\n";
print $gnuPlot "set lmargin screen 0.15\n";
print $gnuPlot "set rmargin screen 0.95\n";
print $gnuPlot "set bmargin screen 0.15\n";
print $gnuPlot "set tmargin screen 0.95\n";
print $gnuPlot "set key spacing 1.2\n";
print $gnuPlot "set key at screen 0.4,0.2\n";
print $gnuPlot "set key left\n";
print $gnuPlot "set key bottom\n";
print $gnuPlot "set logscale xy\n";
print $gnuPlot "set mxtics 10\n";
print $gnuPlot "set mytics 10\n";
print $gnuPlot "set format x '\$10^{\%L}\$'\n";
print $gnuPlot "set format y '\$10^{\%L}\$'\n";
print $gnuPlot "set xrange [1.0e10:1.0e15]\n";
print $gnuPlot "set yrange [1.0e-7:3.0]\n";
print $gnuPlot "set title 'Optimal sampling density'\n";
print $gnuPlot "set xlabel '\$M_{\\rm halo} [M_\\odot]\$'\n";
print $gnuPlot "set ylabel '\$\\gamma(M_{\\rm halo})\$'\n";
print $gnuPlot "set pointsize 2.0\n";
&PrettyPlots::Prepare_Dataset(\$plot,
			      $haloMass,$haloMF,
			      style     => "line",
			      weight    => [5,3],
			      color     => $PrettyPlots::colorPairs{'cornflowerBlue'},
			      title     => 'Halo mass function weight'
    );
&PrettyPlots::Prepare_Dataset(\$plot,
			      $haloMass,$gamma,
			      style     => "line",
			      weight    => [5,3],
			      color     => $PrettyPlots::colorPairs{'redYellow'},
			      title     => 'Optimal weight'
    );
&PrettyPlots::Plot_Datasets($gnuPlot,\$plot);
close($gnuPlot);
&LaTeX::GnuPlot2PDF($plotFileEPS);

exit;
